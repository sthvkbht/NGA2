!> Various definitions and tools for running an NGA2 simulation
module simulation
  use string,             only: str_medium
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use lpt_class,         only: lpt
   use lowmach_class,     only: lowmach
   use sgsmodel_class,    only: sgsmodel
   use vdscalar_class,    only: vdscalar
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use datafile_class,    only: datafile
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single low Mach flow solver, scalar solver, particle solver, and corresponding time tracker
   type(hypre_str),   public :: ps
   type(ddadi),       public :: vs,ss
   type(lowmach),     public :: fs
   type(vdscalar),    public :: sc
   type(lpt),         public :: lp
   type(sgsmodel),    public :: sgs
   type(timetracker), public :: time

   !> Provide a datafile and an event tracker for saving restarts
   type(event)    :: save_evt
   type(datafile) :: df
   logical :: restarted
   
   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,lptfile,consfile,pfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resSC
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,rho0
   real(WP), dimension(:,:,:), allocatable :: srcUlp,srcVlp,srcWlp
   real(WP), dimension(:,:,:), allocatable :: tmp1,tmp2,tmp3
   logical, dimension(:,:,:), allocatable :: flag

   !> Max timestep size for LPT
   integer ::  lp_iter
   real(WP) :: lp_dt,lp_dt_max
   
   !> Equation of state
   real(WP) :: T,P,Schmidt
   real(WP) :: minS,maxS
   
   !> Inlet parameters
   real(WP) :: Sjet,Djet,Vjet,Hjet
   
   !> Integral of pressure residual
   real(WP) :: int_RP=0.0_WP
   
contains
   
   
   !> Obtain density from equation of state based on salt water
   subroutine get_rho()
      implicit none
      integer :: i,j,k
      do k=sc%cfg%kmino_,sc%cfg%kmaxo_
         do j=sc%cfg%jmino_,sc%cfg%jmaxo_
            do i=sc%cfg%imino_,sc%cfg%imaxo_
               sc%rho(i,j,k)=saltwater_eos(sc%SC(i,j,k),T,P)
            end do
         end do
      end do
   end subroutine get_rho

   
   !> Salt water EOS
   function saltwater_eos(S_in,T_in,P_in) result(rho)
     real(WP), intent(in) :: S_in ! Salinity (g·kg-1)
     real(WP), intent(in) :: T_in ! Temperature (K)
     real(WP), intent(in) :: P_in ! Pressure (Pa)
     real(WP) :: Tc   ! Temperature (°C)
     real(WP) :: s    ! Salinity (kg·kg-1)
     real(WP) :: Pc   ! Pressure (MPa)
     real(WP) :: rho  ! Density (kg·m-3)
     real(WP), parameter :: a1 =  9.9992293295E+2_WP
     real(WP), parameter :: a2 =  2.0341179217E-2_WP
     real(WP), parameter :: a3 = -6.1624591598E-3_WP
     real(WP), parameter :: a4 =  2.2614664708E-5_WP
     real(WP), parameter :: a5 = -4.6570659168E-8_WP
     real(WP), parameter :: b1 = 8.0200240891E+2_WP
     real(WP), parameter :: b2 =-2.0005183488E+0_WP
     real(WP), parameter :: b3 = 1.6771024982E-2_WP
     real(WP), parameter :: b4 =-3.0600536746E-5_WP
     real(WP), parameter :: b5 =-1.6132224742E-5_WP
     real(WP), parameter :: c1 =  5.0792E-04_WP
     real(WP), parameter :: c2 = -3.4168E-06_WP
     real(WP), parameter :: c3 =  5.6931E-08_WP
     real(WP), parameter :: c4 = -3.7263E-10_WP
     real(WP), parameter :: c5 =  1.4465E-12_WP
     real(WP), parameter :: c6 = -1.7058E-15_WP
     real(WP), parameter :: c7 = -1.3389E-06_WP
     real(WP), parameter :: c8 =  4.8603E-09_WP
     real(WP), parameter :: c9 = -6.8039E-13_WP
     real(WP), parameter :: d1 = -1.1077E-06_WP
     real(WP), parameter :: d2 =  5.5584E-09_WP
     real(WP), parameter :: d3 = -4.2539E-11_WP
     real(WP), parameter :: d4 =  8.3702e-09_WP
     real(WP) rho0, Fp0, Fp1, Fp
     Tc = T_in - 273.15_WP
     s  = 1.0e-3_WP*max(S_in,0.0_WP)
     Pc = 1.0e-6_WP*P_in
     rho0 = a1 + Tc*(a2 + Tc*(a3 + Tc*(a4 + Tc*a5))) + &
          &   s*(b1 + Tc*(b2 + Tc*(b3 + Tc*b4 + s*b5)))
     Fp0 = c1 + Tc*(c2 + Tc*(c3 + Tc*(c4 + Tc*(c5 + Tc*c6)))) &
          + S_in*(d1 + Tc*(d2 + Tc*d3))
     Fp1 = c7 + Tc*(c8 + Tc*Tc*c9) + d4*S_in
     Fp = exp( Pc*Fp0 + 0.5*Pc*Pc*Fp1 )
     rho = rho0 * Fp
   end function saltwater_eos


   !> Obtain viscosity of salt water
   subroutine get_visc()
      implicit none
      integer :: i,j,k
      do k=sc%cfg%kmino_,sc%cfg%kmaxo_
         do j=sc%cfg%jmino_,sc%cfg%jmaxo_
            do i=sc%cfg%imino_,sc%cfg%imaxo_
               fs%visc(i,j,k)=saltwater_visc(sc%SC(i,j,k),T)
               sc%diff(i,j,k)=fs%visc(i,j,k)/Schmidt
            end do
         end do
      end do
   end subroutine get_visc
   
   !> Computer viscosity based on salinity
   function saltwater_visc(S_in,T_in) result(visc)
     real(WP), intent(in) :: T_in            !< Temperature [K]
     real(WP), intent(in) :: S_in            !< Salinity    [g/kg]
     real(WP), parameter :: A0 = 2.414E-5_WP !< [Pa.s]
     real(WP), parameter :: B0 = 247.8_WP    !< [K]
     real(WP), parameter :: C0 = 140.0_WP    !< [K]
     real(WP), parameter :: T0 = 273.15_WP   !< [K]
     real(WP) :: visc,visc_w,A,B,Tc,Sc
     Tc=T_in-273.15_WP
     Sc=max(0.0_WP,S_in)
     visc_w=A0*10.0_WP**(B0/(Tc+T0-C0))
     A=1.541E-3_WP+1.998E-5_WP*Tc-9.52E-8_WP*Tc**2
     B=7.974E-6_WP-7.561E-8_WP*Tc+4.724E-10_WP*Tc**2
     visc=visc_w*(1.0_WP+A*Sc+B*Sc**2)
   end function saltwater_visc


   !> Initialization of problem solver
   subroutine simulation_init
     use param, only: param_read
     implicit none


     ! Initialize time tracker with 2 subiterations
     initialize_timetracker: block
       time=timetracker(amRoot=cfg%amRoot,name="plume")
       call param_read('Max timestep size',time%dtmax)
       call param_read('Max cfl number',time%cflmax)
       time%dt=time%dtmax
       time%itmax=2
     end block initialize_timetracker

     
     ! Handle restart/saves here
     restart_and_save: block
       character(len=str_medium) :: timestamp
       ! Create event for saving restart files
       save_evt=event(time,'Restart output')
       call param_read('Restart output period',save_evt%tper)
       ! Check if we are restarting
       call param_read(tag='Restart from',val=timestamp,short='r',default='')
       restarted=.false.; if (len_trim(timestamp).gt.0) restarted=.true.
       if (restarted) then
          ! If we are, read the name of the directory
          call param_read('Restart from',timestamp,'r')
          ! Read the datafile
          df=datafile(pg=cfg,fdata='restart/data_'//trim(adjustl(timestamp)))
       else
          ! Prepare a new directory for storing files for restart
          call execute_command_line('mkdir -p restart')
          ! If we are not restarting, we will still need a datafile for saving restart files
          df=datafile(pg=cfg,filename=trim(cfg%name),nval=2,nvar=5)
          df%valname(1)='t'
          df%valname(2)='dt'
          df%varname(1)='U'
          df%varname(2)='V'
          df%varname(3)='W'
          df%varname(4)='P'
          df%varname(5)='S'
       end if
     end block restart_and_save


     ! Revisit timetracker to adjust time and time step values if this is a restart
     update_timetracker: block
       if (restarted) then
          call df%pullval(name='t' ,val=time%t )
          call df%pullval(name='dt',val=time%dt)
          time%told=time%t-time%dt
       end if
     end block update_timetracker
    
     
     ! Create a low-Mach flow solver with bconds
     create_velocity_solver: block
       use hypre_str_class, only: pcg_pfmg
       use lowmach_class,   only: dirichlet,clipped_neumann
       ! Create flow solver
       fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
       ! Assign acceleration of gravity
       call param_read('Gravity',fs%gravity)
       ! Read in EOS parameters
       call param_read('Temperature',T,default=293.15_WP)
       call param_read('Pressure',P,default=101325.0_WP)
       ! Read in inlet parameters
       call param_read('S jet',Sjet)
       call param_read('D jet',Djet)
       call param_read('V jet',Vjet)
       call param_read('H jet',Hjet)
       ! Define boundary conditions
       call fs%add_bcond(name='bottom',type=dirichlet      ,face='y',dir=-1,canCorrect=.false.,locator=ym_locator)
       call fs%add_bcond(name='top'   ,type=clipped_neumann,face='y',dir=+1,canCorrect=.true. ,locator=yp_locator)
       call fs%add_bcond(name='left'  ,type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
       call fs%add_bcond(name='right ',type=dirichlet,face='x',dir=+1,canCorrect=.false.,locator=xp_locator)
       ! Prepare and configure pressure solver
       ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
       call param_read('Pressure iteration',ps%maxit)
       call param_read('Pressure tolerance',ps%rcvg)
       ! Configure implicit velocity solver
       vs=ddadi(cfg=cfg,name='Velocity',nst=7)
       ! Setup the solver
       call fs%setup(pressure_solver=ps,implicit_solver=vs)
     end block create_velocity_solver

     
     ! Create a scalar solver
     create_scalar: block
       use vdscalar_class, only: dirichlet,neumann,bquick
       ! Create scalar solver
       sc=vdscalar(cfg=cfg,scheme=bquick,name='Salinity')
       ! Define boundary conditions
       call sc%add_bcond(name='bottom',type=dirichlet,locator=ym_locator_sc)
       call sc%add_bcond(name='top'   ,type=neumann  ,locator=yp_locator   ,dir='+y')
       call sc%add_bcond(name='left'  ,type=neumann  ,locator=xm_locator_sc,dir='-x')
       call sc%add_bcond(name='right' ,type=neumann  ,locator=xp_locator   ,dir='+x')
       ! Read in Schmidt number
       call param_read('Schmidt number',Schmidt)
       ! Configure implicit scalar solver
       ss=ddadi(cfg=cfg,name='Scalar',nst=13)
       ! Setup the solver
       call sc%setup(implicit_solver=ss)
     end block create_scalar

     
     ! Allocate work arrays
     allocate_work_arrays: block
       ! Flow solver
       allocate(gradU(1:3,1:3,fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
       allocate(resU(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
       allocate(resV(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
       allocate(resW(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
       allocate(Ui  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
       allocate(Vi  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
       allocate(Wi  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
       ! Scalar solver
       allocate(resSC(sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_))
       ! Particle solver
       allocate(srcUlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
       allocate(srcVlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
       allocate(srcWlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
       allocate(rho0    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
       allocate(tmp1    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
       allocate(tmp2    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
       allocate(tmp3    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
       allocate(flag    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
     end block allocate_work_arrays

     
     ! Initialize our LPT
     initialize_lpt: block
       use random, only: random_uniform
       character(len=str_medium) :: timestamp
       ! Create solver
       lp=lpt(cfg=cfg,name='LPT')
       ! Get particle density from the input
       call param_read('Particle density',lp%rho)
       ! Set gravity
       call param_read('Gravity',lp%gravity)
       ! Set filter scale to 3.5*dx
       lp%filter_width=3.5_WP*cfg%min_meshsize
       ! Initialize particles
       if (restarted) then
          call param_read('Restart from',timestamp,'r')
          ! Read the part file
          call lp%read(filename='restart/part_'//trim(adjustl(timestamp)))
       else
          ! Start with zero particles
          call lp%resize(0)
       end if
       ! Get initial particle volume fraction
       call lp%update_VF()
       ! Maximum timestep size used for particles
       call param_read('Particle timestep size',lp_dt_max,default=huge(1.0_WP))
       lp_dt=lp_dt_max
       ! Set collision timescale
       call param_read('Collision timescale',lp%tau_col,default=15.0_WP*lp_dt_max)
       ! Set coefficient of restitution
       call param_read('Coefficient of restitution',lp%e_n)
       call param_read('Wall restitution',lp%e_w)
       call param_read('Friction coefficient',lp%mu_f)
       ! Injection parameters
       call param_read('Particle mass flow rate',lp%mfr)
       call param_read('Particle velocity',lp%inj_vel)
       call param_read('Particle mean diameter',lp%inj_dmean)
       call param_read('Particle standard deviation',lp%inj_dsd,default=0.0_WP)
       call param_read('Particle min diameter',lp%inj_dmin,default=tiny(1.0_WP))
       call param_read('Particle max diameter',lp%inj_dmax,default=huge(1.0_WP))
       call param_read('Particle diameter shift',lp%inj_dshift,default=0.0_WP)
       if (lp%inj_dsd.le.epsilon(1.0_WP)) then
          lp%inj_dmin=lp%inj_dmean
          lp%inj_dmax=lp%inj_dmean
       end if
       lp%inj_d=Djet-lp%inj_dmax
       lp%inj_pos(2)=lp%cfg%y(lp%cfg%jmin)+lp%inj_dmax
       lp%inj_pos(1)=0.0_WP; lp%inj_pos(3)=0.0_WP
       ! Update Gib to be consistent with LPT collisions
       cfg%Gib=-cfg%Gib
       cfg%Nib=-cfg%Nib
     end block initialize_lpt


     ! Create partmesh object for Lagrangian particle output
     create_pmesh: block
       integer :: i
       pmesh=partmesh(nvar=1,nvec=1,name='lpt')
       pmesh%varname(1)='diameter'
       pmesh%vecname(1)='velocity'
       call lp%update_partmesh(pmesh)
       do i=1,lp%np_
          pmesh%var(1,i)=lp%p(i)%d
          pmesh%vec(:,1,i)=lp%p(i)%vel
       end do
     end block create_pmesh


     ! Initialize our mixture fraction field
     initialize_scalar: block
       use vdscalar_class, only: bcond
       integer :: n,i,j,k
       type(bcond), pointer :: mybc
       real(WP) :: y
       if (restarted) then
          call df%pullvar(name='S',var=sc%SC)
       else
          ! Initialize stratified salinity profile
          do k=sc%cfg%kmino_,sc%cfg%kmaxo_
             do j=sc%cfg%jmino_,sc%cfg%jmaxo_
                do i=sc%cfg%imino_,sc%cfg%imaxo_
                   y=sc%cfg%ym(j)
                   sc%SC(i,j,k)=max(0.003153_WP*y**3-167.9_WP*y**2+4.758_WP*y+43.29_WP,0.0_WP)
                end do
             end do
          end do
       end if
       y=sc%cfg%ym(sc%cfg%jmin); minS=max(0.003153_WP*y**3-167.9_WP*y**2+4.758_WP*y+43.29_WP,0.0_WP)
       y=sc%cfg%ym(sc%cfg%jmax); maxS=max(0.003153_WP*y**3-167.9_WP*y**2+4.758_WP*y+43.29_WP,0.0_WP)
       ! Apply BCs
       call sc%get_bcond('bottom',mybc)
       do n=1,mybc%itr%no_
          i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
          sc%SC(i,j,k)=Sjet
       end do
       ! Compute density
       call get_rho()
       ! Compute viscosity and diffusivity
       call get_visc()
     end block initialize_scalar

     
     ! Initialize our velocity field
     initialize_velocity: block
       use lowmach_class, only: bcond
       integer :: n,i,j,k
       type(bcond), pointer :: mybc
       real(WP) :: radius,myS
       ! Zero initial field
       if (restarted) then
          call df%pullvar(name='U',var=fs%U)
          call df%pullvar(name='V',var=fs%V)
          call df%pullvar(name='W',var=fs%W)
          call df%pullvar(name='P',var=fs%P)
       else
          fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP; fs%P=0.0_WP
       end if
       ! Set density from scalar
       rho0=sc%rho
       sc%rho=sc%rho*(1.0_WP-lp%VF); fs%rho=sc%rho

       ! Form momentum
       call fs%rho_multiply
       ! Apply BCs
       call fs%apply_bcond(time%t,time%dt)
       call fs%get_bcond('bottom',mybc)
       do n=1,mybc%itr%no_
          i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
          radius=norm2([fs%cfg%xm(i),fs%cfg%zm(k)]-[0.0_WP,0.0_WP])
          myS           =Sjet
          fs%V(i,j,k)   =Vjet
          fs%rhoV(i,j,k)=fs%V(i,j,k)*saltwater_eos(myS,T,P)
       end do
       ! Get cell-centered velocities and continuity residual
       call fs%interp_vel(Ui,Vi,Wi)
       resSC=0.0_WP; call fs%get_div(drhodt=resSC)
       ! Compute MFR through all boundary conditions
       call fs%get_mfr()
     end block initialize_velocity

     ! Create an LES model
     create_sgs: block
       sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
     end block create_sgs

     ! Add Ensight output
     create_ensight: block
       ! Create Ensight output from cfg
       ens_out=ensight(cfg=cfg,name='plume')
       ! Create event for Ensight output
       ens_evt=event(time=time,name='Ensight output')
       call param_read('Ensight output period',ens_evt%tper)
       ! Add variables to output
       call ens_out%add_particle('particles',pmesh)
       call ens_out%add_scalar('levelset',cfg%Gib)
       call ens_out%add_scalar('pressure',fs%P)
       call ens_out%add_vector('velocity',Ui,Vi,Wi)
       call ens_out%add_scalar('density',rho0)
       call ens_out%add_scalar('viscosity',fs%visc)
       call ens_out%add_scalar('salinity',sc%SC)
       call ens_out%add_scalar('epsp',lp%VF)
       call ens_out%add_scalar('visc_sgs',sgs%visc)
       ! Output to ensight
       if (ens_evt%occurs()) call ens_out%write_data(time%t)
     end block create_ensight

     ! Create a monitor file
     create_monitor: block
       ! Prepare some info about fields
       call fs%get_cfl(time%dt,time%cfl)
       call fs%get_max()
       call sc%get_max()
       call sc%get_int()
       ! Create simulation monitor
       mfile=monitor(fs%cfg%amRoot,'simulation')
       call mfile%add_column(time%n,'Timestep number')
       call mfile%add_column(time%t,'Time')
       call mfile%add_column(time%dt,'Timestep size')
       call mfile%add_column(time%cfl,'Maximum CFL')
       call mfile%add_column(fs%Umax,'Umax')
       call mfile%add_column(fs%Vmax,'Vmax')
       call mfile%add_column(fs%Wmax,'Wmax')
       call mfile%add_column(sc%SCmax,'Smax')
       call mfile%add_column(sc%SCmin,'Smin')
       call mfile%add_column(sc%rhomax,'RHOmax')
       call mfile%add_column(sc%rhomin,'RHOmin')
       call mfile%write()
       ! Create pressure monitor
       pfile=monitor(fs%cfg%amRoot,'pressure')
       call pfile%add_column(time%n,'Timestep number')
       call pfile%add_column(time%t,'Time')
       call pfile%add_column(fs%Pmax,'Pmax')
       call pfile%add_column(int_RP,'Int(RP)')
       call pfile%add_column(fs%divmax,'Maximum divergence')
       call pfile%add_column(fs%psolv%it,'Pressure iteration')
       call pfile%add_column(fs%psolv%rerr,'Pressure error')
       call pfile%write()
       ! Create CFL monitor
       cflfile=monitor(fs%cfg%amRoot,'cfl')
       call cflfile%add_column(time%n,'Timestep number')
       call cflfile%add_column(time%t,'Time')
       call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
       call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
       call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
       call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
       call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
       call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
       call cflfile%add_column(lp%CFL_col,'Collision CFL')
       call cflfile%write()
       ! Create conservation monitor
       consfile=monitor(fs%cfg%amRoot,'conservation')
       call consfile%add_column(time%n,'Timestep number')
       call consfile%add_column(time%t,'Time')
       call consfile%add_column(sc%SCint,'SC integral')
       call consfile%add_column(sc%rhoint,'RHO integral')
       call consfile%add_column(sc%rhoSCint,'rhoSC integral')
       call consfile%write()
       ! Create LPT monitor
       lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
       call lptfile%add_column(time%n,'Timestep number')
       call lptfile%add_column(time%t,'Time')
       call lptfile%add_column(lp_dt,'Particle dt')
       call lptfile%add_column(lp_iter,'Particle iter')
       call lptfile%add_column(lp%np,'Particle number')
       call lptfile%add_column(lp%VFmean,'VFp mean')
       call lptfile%add_column(lp%VFmax,'VFp max')
       call lptfile%add_column(lp%Umin,'Particle Umin')
       call lptfile%add_column(lp%Umax,'Particle Umax')
       call lptfile%add_column(lp%Vmin,'Particle Vmin')
       call lptfile%add_column(lp%Vmax,'Particle Vmax')
       call lptfile%add_column(lp%Wmin,'Particle Wmin')
       call lptfile%add_column(lp%Wmax,'Particle Wmax')
       call lptfile%add_column(lp%dmin,'Particle dmin')
       call lptfile%add_column(lp%dmax,'Particle dmax')
       call lptfile%write()
     end block create_monitor

   contains

      !> Function that localizes the x- boundary
      function xm_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (pg%ym(j).gt.Hjet.and.i.eq.pg%imin) isIn=.true.
      end function xm_locator

      !> Function that localizes the x- boundary for SC
      function xm_locator_sc(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (pg%ym(j).gt.Hjet.and.i.eq.pg%imin-1) isIn=.true.
      end function xm_locator_sc
       
      !> Function that localizes the x+ boundary
      function xp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (pg%ym(j).gt.Hjet.and.i.eq.pg%imax+1) isIn=.true.
      end function xp_locator
       
      !> Function that localizes y- boundary
      function ym_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         real(WP) :: radius
         logical :: isIn
         isIn=.false.
         radius=norm2([pg%xm(i),pg%zm(k)]-[0.0_WP,0.0_WP])
         if (j.eq.pg%jmin.and.radius.le.Djet) isIn=.true.
      end function ym_locator

      !> Function that localizes the y- boundary for SC
      function ym_locator_sc(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         real(WP) :: radius
         logical :: isIn
         isIn=.false.
         radius=norm2([pg%xm(i),pg%zm(k)]-[0.0_WP,0.0_WP])
         if (j.eq.pg%jmin-1.and.radius.le.Djet) isIn=.true.
      end function ym_locator_sc
       
      !> Function that localizes y+ boundary
      function yp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmax+1) isIn=.true.
      end function yp_locator
      
      !> Function that localizes z- boundary
      function zm_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (k.eq.pg%kmin) isIn=.true.
      end function zm_locator
      
      !> Function that localizes z+ boundary
      function zp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (k.eq.pg%kmax+1) isIn=.true.
      end function zp_locator
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old scalar
         sc%rhoold=sc%rho
         sc%SCold =sc%SC
         
         ! Remember old velocity and momentum
         fs%rhoold=fs%rho
         fs%Uold=fs%U; fs%rhoUold=fs%rhoU
         fs%Vold=fs%V; fs%rhoVold=fs%rhoV
         fs%Wold=fs%W; fs%rhoWold=fs%rhoW

         ! Particle update
         lpt: block
           real(WP) :: dt_done,mydt,cfl
           integer :: i
           ! Inject particles
           call lp%inject(dt=time%dt,face='y',avoid_overlap=.true.)
           ! Get fluid stress
           call fs%get_div_stress(resU,resV,resW)
           ! Zero-out LPT source terms
           srcUlp=0.0_WP; srcVlp=0.0_WP; srcWlp=0.0_WP
           ! Sub-iterate
           call lp%get_cfl(lp_dt,cflc=cfl,cfl=cfl)
           if (cfl.gt.0.0_WP) lp_dt=min(lp_dt*time%cflmax/cfl,lp_dt_max)
           dt_done=0.0_WP; lp_iter=0
           do while (dt_done.lt.time%dtmid)
              ! Decide the timestep size
              mydt=min(lp_dt,time%dtmid-dt_done)
              ! Collide and advance particles
              call lp%collide(dt=mydt,Gib=cfg%Gib,Nxib=cfg%Nib(1,:,:,:),Nyib=cfg%Nib(2,:,:,:),Nzib=cfg%Nib(3,:,:,:))
              call lp%advance(dt=mydt,U=fs%U,V=fs%V,W=fs%W,rho=rho0,visc=fs%visc,stress_x=resU,stress_y=resV,stress_z=resW,&
                   srcU=tmp1,srcV=tmp2,srcW=tmp3)
              srcUlp=srcUlp+tmp1
              srcVlp=srcVlp+tmp2
              srcWlp=srcWlp+tmp3
              ! Increment
              dt_done=dt_done+mydt
              lp_iter=lp_iter+1
           end do
         end block lpt

         ! Turbulence modeling
         sgs_modeling: block
           use sgsmodel_class, only: vreman,WALE
           call fs%get_gradu(gradU)
           call sgs%get_visc(type=WALE,dt=time%dtold,rho=rho0,gradu=gradU)
         end block sgs_modeling

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! ============= SCALAR SOLVER =======================

            ! Reset interpolation metrics to QUICK scheme
            call sc%metric_reset()
            
            ! Build mid-time scalar
            sc%SC=0.5_WP*(sc%SC+sc%SCold)
            
            ! Explicit calculation of drhoSC/dt from scalar equation
            call sc%get_drhoSCdt(resSC,fs%rhoU,fs%rhoV,fs%rhoW)

            ! Perform bquick procedure
            bquick: block
               integer :: i,j,k
               ! Assemble explicit residual
               resSC=time%dt*resSC-(2.0_WP*sc%rho*sc%SC-(sc%rho+sc%rhoold)*sc%SCold)
               ! Apply it to get explicit scalar prediction
               tmp1=2.0_WP*sc%SC-sc%SCold+resSC/sc%rho
               ! Check cells that require bquick
               do k=sc%cfg%kmino_,sc%cfg%kmaxo_
                  do j=sc%cfg%jmino_,sc%cfg%jmaxo_
                     do i=sc%cfg%imino_,sc%cfg%imaxo_
                        if (tmp1(i,j,k).le.minS.or.tmp1(i,j,k).ge.maxS) then
                           flag(i,j,k)=.true.
                        else
                           flag(i,j,k)=.false.
                        end if
                     end do
                  end do
               end do
               ! Adjust metrics
               call sc%metric_adjust(flag)
               ! Recompute drhoSC/dt
               call sc%get_drhoSCdt(resSC,fs%rhoU,fs%rhoV,fs%rhoW)
            end block bquick
            
            ! Assemble explicit residual
            resSC=time%dt*resSC-(2.0_WP*sc%rho*sc%SC-(sc%rho+sc%rhoold)*sc%SCold)
            
            ! Form implicit residual
            call sc%solve_implicit(time%dt,resSC,fs%rhoU,fs%rhoV,fs%rhoW)
            
            ! Advance scalar field
            sc%SC=2.0_WP*sc%SC-sc%SCold+resSC
            
            ! Apply all other boundary conditions on the resulting field
            call sc%apply_bcond(time%t,time%dt)
            dirichlet_scalar: block
               use vdscalar_class, only: bcond
               integer :: n,i,j,k
               type(bcond), pointer :: mybc
               call sc%get_bcond('bottom',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  sc%SC(i,j,k)=Sjet
               end do
             end block dirichlet_scalar

             ! Apply IBM to enforce Neumann BC at the walls
             ibm_sc: block
              integer :: i,j,k,ii,jj,kk
              real(WP) :: sum_VF,sum_VFS
              do k=fs%cfg%kmin_,fs%cfg%kmax_
                 do j=fs%cfg%jmin_,fs%cfg%jmax_
                    do i=fs%cfg%imin_,fs%cfg%imax_
                       sum_VF=0.0_WP; sum_VFS=0.0_WP
                       do kk=-1,1; do jj=-1,1; do ii=-1,1
                          if (ii.eq.0.and.jj.eq.0.and.kk.eq.0) cycle
                          sum_VF =sum_VF +cfg%VF(i+ii,j+jj,k+kk)
                          sum_VFS=sum_VFS+cfg%VF(i+ii,j+jj,k+kk)*sc%SC(i+ii,j+jj,k+kk)
                       end do; end do; end do
                       if (sum_VF.gt.0.0_WP) then
                          sc%SC(i,j,k)=cfg%VF(i,j,k)*sc%SC(i,j,k)+(1.0_WP-cfg%VF(i,j,k))*sum_VFS/sum_VF
                       end if
                    end do
                 end do
              end do
              call sc%cfg%sync(sc%SC)
            end block ibm_sc
            ! ===================================================
            
            ! ============ UPDATE PROPERTIES ====================
            ! Update density
            call get_rho()
            rho0=sc%rho
            sc%rho=sc%rho*(1.0_WP-lp%VF)
            
            ! Update the transport variables
            call get_visc()
            fs%visc=fs%visc+sgs%visc
            sc%diff=sc%diff+sgs%visc/Schmidt
            ! ===================================================
            
            ! ============ VELOCITY SOLVER ======================
            
            ! Build n+1 density
            fs%rho=0.5_WP*(sc%rho+sc%rhoold)
            
            ! Build mid-time velocity and momentum
            fs%U=0.5_WP*(fs%U+fs%Uold); fs%rhoU=0.5_WP*(fs%rhoU+fs%rhoUold)
            fs%V=0.5_WP*(fs%V+fs%Vold); fs%rhoV=0.5_WP*(fs%rhoV+fs%rhoVold)
            fs%W=0.5_WP*(fs%W+fs%Wold); fs%rhoW=0.5_WP*(fs%rhoW+fs%rhoWold)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)

            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=time%dtmid*resU-(2.0_WP*fs%rhoU-2.0_WP*fs%rhoUold)
            resV=time%dtmid*resV-(2.0_WP*fs%rhoV-2.0_WP*fs%rhoVold)
            resW=time%dtmid*resW-(2.0_WP*fs%rhoW-2.0_WP*fs%rhoWold)

            ! Add momentum source term from lpt
            add_lpt_src: block
              integer :: i,j,k
              do k=fs%cfg%kmin_,fs%cfg%kmax_
                 do j=fs%cfg%jmin_,fs%cfg%jmax_
                    do i=fs%cfg%imin_,fs%cfg%imax_
                       if (fs%umask(i,j,k).eq.0) resU(i,j,k)=resU(i,j,k)+sum(fs%itpr_x(:,i,j,k)*srcUlp(i-1:i,j,k))
                       if (fs%vmask(i,j,k).eq.0) resV(i,j,k)=resV(i,j,k)+sum(fs%itpr_y(:,i,j,k)*srcVlp(i,j-1:j,k))
                       if (fs%wmask(i,j,k).eq.0) resW(i,j,k)=resW(i,j,k)+sum(fs%itpr_z(:,i,j,k)*srcWlp(i,j,k-1:k))
                    end do
                 end do
              end do
            end block add_lpt_src
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dtmid,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Update momentum
            call fs%rho_multiply()

            ! Apply boundary conditions
            call fs%apply_bcond(time%tmid,time%dtmid)
            dirichlet_velocity: block
               use lowmach_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k
               real(WP) :: myS
               call fs%get_bcond('bottom',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  myS           =Sjet
                  fs%V(i,j,k)   =Vjet
                  fs%rhoV(i,j,k)=fs%V(i,j,k)*saltwater_eos(myS,T,P)
               end do
            end block dirichlet_velocity

            ! Apply IBM to enforce no-slip BC at the walls
            ibm_vel: block
              integer :: i,j,k
              do k=fs%cfg%kmin_,fs%cfg%kmax_
                 do j=fs%cfg%jmin_,fs%cfg%jmax_
                    do i=fs%cfg%imin_,fs%cfg%imax_
                       fs%U(i,j,k)=sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))*fs%U(i,j,k)
                       fs%V(i,j,k)=sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))*fs%V(i,j,k)
                       fs%W(i,j,k)=sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))*fs%W(i,j,k)
                    end do
                 end do
              end do
              call fs%cfg%sync(fs%U)
              call fs%cfg%sync(fs%V)
              call fs%cfg%sync(fs%W)
              call fs%rho_multiply()
            end block ibm_vel
            
            ! Solve Poisson equation
            call sc%get_drhodt(dt=time%dt,drhodt=resSC)
            call fs%correct_mfr(drhodt=resSC)
            call fs%get_div(drhodt=resSC)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dtmid
            call cfg%integrate(A=fs%psolv%rhs,integral=int_RP)
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct momentum and rebuild velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%rhoU=fs%rhoU-time%dtmid*resU
            fs%rhoV=fs%rhoV-time%dtmid*resV
            fs%rhoW=fs%rhoW-time%dtmid*resW
            call fs%rho_divide
            ! ===================================================
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call sc%get_drhodt(dt=time%dt,drhodt=resSC)
         call fs%get_div(drhodt=resSC)

         ! Output to ensight
         if (ens_evt%occurs()) then
            update_pmesh: block
              integer :: i
              call lp%update_partmesh(pmesh)
              do i=1,lp%np_
                 pmesh%var(1,i)=lp%p(i)%d
                 pmesh%vec(:,1,i)=lp%p(i)%vel
              end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if

         ! Perform and output monitoring
         call fs%get_max()
         call sc%get_max()
         call sc%get_int()
         call lp%get_max()
         call mfile%write()
         call pfile%write()
         call cflfile%write()
         call consfile%write()
         call lptfile%write()

         ! Finally, see if it's time to save restart files
         if (save_evt%occurs()) then
            save_restart: block
              character(len=str_medium) :: timestamp
              ! Prefix for files
              write(timestamp,'(es12.5)') time%t
              ! Populate df and write it
              call df%pushval(name='t' ,val=time%t )
              call df%pushval(name='dt',val=time%dt)
              call df%pushvar(name='U' ,var=fs%U   )
              call df%pushvar(name='V' ,var=fs%V   )
              call df%pushvar(name='W' ,var=fs%W   )
              call df%pushvar(name='P' ,var=fs%P   )
              call df%pushvar(name='S' ,var=sc%SC  )
              call df%write(fdata='restart/data_'//trim(adjustl(timestamp)))
              ! Write particle file
              call lp%write(filename='restart/part_'//trim(adjustl(timestamp)))
            end block save_restart
         end if
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Deallocate work arrays
      deallocate(resSC,resU,resV,resW,Ui,Vi,Wi,srcUlp,srcVlp,srcWlp,rho0,gradU,tmp1,tmp2,tmp3,flag)
   end subroutine simulation_final
   

end module simulation
