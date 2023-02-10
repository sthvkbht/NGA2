!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   !use hypre_str_class,   only: hypre_str
   use fourier3d_class,   only: fourier3d
   use incomp_class,      only: incomp
   use lpt_class,         only: lpt
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use sgsmodel_class,    only: sgsmodel
   use datafile_class,    only: datafile
   use string,            only: str_medium
   implicit none
   private
   
   !> Single-phase incompressible flow solver, pressure and implicit solvers, and a time tracker
   !type(hypre_str),   public :: ps
   type(fourier3d),   public :: ps
   type(incomp),      public :: fs
   type(timetracker), public :: time
   type(lpt),         public :: lp
   type(sgsmodel),    public :: sgs
   
   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt

   !> Provide a datafile and an event tracker for saving restarts
   type(event)    :: save_evt
   type(datafile) :: df
   logical :: restarted
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,hitfile,lptfile,sgsfile,tfile,ssfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:), allocatable :: SR
   real(WP), dimension(:,:,:,:,:), allocatable :: gradu
   
   !> Fluid, forcing, and particle parameters
   real(WP) :: visc,meanU,meanV,meanW
   real(WP) :: Urms0,TKE0,EPS0,Re_max
   real(WP) :: TKE,URMS,sgsTKE,EPSp
   real(WP) :: stk,tau_eta
   real(WP) :: meanvisc,Lx,N
   real(WP) :: tauinf,G,Gdtau,Gdtaui,dx
   logical  :: linforce,use_sgs,maxRe,use_crw
   integer  :: sgs_type

   !> For monitoring
   real(WP) :: EPS,sgsEPSp
   real(WP) :: Re_L,Re_lambda
   real(WP) :: eta,ell
   real(WP) :: dx_eta,ell_Lx,Re_ratio,eps_ratio,tke_ratio,nondtime

   !> Wallclock time for monitoring
   type :: timer
      real(WP) :: time_in
      real(WP) :: time
      real(WP) :: percent
   end type timer
   type(timer) :: wt_total,wt_vel,wt_pres,wt_lpt,wt_rest,wt_sgs,wt_stat,wt_force
   
contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(SR  (1:6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(gradu(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max iter',time%nmax)
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
           ! If we are not restarting, we will still need a datafile for saving restart files
           df=datafile(pg=cfg,filename=trim(cfg%name),nval=2,nvar=4)
           df%valname(1)='t'
           df%valname(2)='dt'
           df%varname(1)='U'
           df%varname(2)='V'
           df%varname(3)='W'
           df%varname(4)='P'
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

      ! Initialize timers
      initialize_timers: block
         wt_total%time=0.0_WP; wt_total%percent=0.0_WP
         wt_vel%time=0.0_WP;   wt_vel%percent=0.0_WP
         wt_pres%time=0.0_WP;  wt_pres%percent=0.0_WP
         wt_lpt%time=0.0_WP;   wt_lpt%percent=0.0_WP
         wt_rest%time=0.0_WP;  wt_rest%percent=0.0_WP
         wt_sgs%time=0.0_WP;   wt_sgs%percent=0.0_WP
         wt_stat%time=0.0_WP;  wt_stat%percent=0.0_WP
         wt_force%time=0.0_WP; wt_force%percent=0.0_WP
      end block initialize_timers
      
      ! Create a single-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg
         ! Create flow solver
         fs=incomp(cfg=cfg,name='NS solver')
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign constant density
         call param_read('Density',fs%rho)
         ! Prepare and configure pressure solver
         !ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         !ps%maxlevel=10
         !call param_read('Pressure iteration',ps%maxit)
         !call param_read('Pressure tolerance',ps%rcvg)
         ps=fourier3d(cfg=cfg,name='Pressure',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps)
      end block create_and_initialize_flow_solver

      ! Create an LES model
      create_sgs: block
         call param_read(tag='Use SGS model',val=use_sgs,default=.false.)
         if (use_sgs) call param_read('SGS model type',sgs_type)
         sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
         sgs%Cs_ref=0.1_WP
      end block create_sgs

      ! Prepare initial velocity field
      initialize_velocity: block
         use random,   only: random_normal
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
         use parallel, only: MPI_REAL_WP
         use mathtools,only: Pi
         integer :: i,j,k
         ! Read in forcing, grid, and initial velocity field parameters
         call param_read('Linear forcing',linforce)
         call param_read('Force to maximum Re_lambda',maxRe)
         call param_read('Forcing constant (G)', G)
         call param_read('Lx', Lx)
         call param_read('nx', N)
         dx=Lx/N
         if (linforce) then
            if (maxRE) then
                EPS0 = (visc/fs%rho)**3*(Pi*N/(1.5_WP*Lx))**4
                TKE0 = 1.5_WP*(0.2_WP*Lx*EPS0)**(0.6667_WP)
            else
                call param_read('Steady-state TKE',TKE0)
                EPS0 = 5.0_WP*(0.6667_WP*TKE0)**1.5_WP / Lx
            end if
            Re_max = sqrt(15.0_WP*sqrt(0.6667_WP*TKE0)*0.2_WP*Lx/visc)
            tauinf = 2.0_WP*TKE0/(3.0_WP*EPS0)
            Gdtau = G/tauinf
            Gdtaui= 1/Gdtau
         else
            TKE0 = 0.0_WP
            tauinf = 99999.0_WP
         end if

         if (restarted) then
            call df%pullvar(name='U',var=fs%U)
            call df%pullvar(name='V',var=fs%V)
            call df%pullvar(name='W',var=fs%W)
            call df%pullvar(name='P',var=fs%P)
         else
            Urms0 = sqrt(0.6667_WP*TKE0)
            ! Gaussian initial field
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do j=fs%cfg%jmin_,fs%cfg%jmax_
                  do i=fs%cfg%imin_,fs%cfg%imax_
                     fs%U(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
                     fs%V(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
                     fs%W(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
                  end do
               end do
            end do
            call fs%cfg%sync(fs%U)
            call fs%cfg%sync(fs%V)
            call fs%cfg%sync(fs%W)
         end if
         
         ! Compute mean and remove it from the velocity field to obtain <U>=0
         call fs%cfg%integrate(A=fs%U,integral=meanU); meanU=meanU/fs%cfg%vol_total
         call fs%cfg%integrate(A=fs%V,integral=meanV); meanV=meanV/fs%cfg%vol_total
         call fs%cfg%integrate(A=fs%W,integral=meanW); meanW=meanW/fs%cfg%vol_total
 
         fs%U = fs%U - meanU
         fs%V = fs%V - meanV
         fs%W = fs%W - meanW
 
         ! Project to ensure divergence-free
         call fs%get_div()
         fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
         fs%psolv%sol=0.0_WP
         call fs%psolv%solve()
         call fs%shift_p(fs%psolv%sol)
         call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
         fs%P=fs%P+fs%psolv%sol
         fs%U=fs%U-time%dt*resU/fs%rho
         fs%V=fs%V-time%dt*resV/fs%rho
         fs%W=fs%W-time%dt*resW/fs%rho
 
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block initialize_velocity
      
      ! Compute turbulence statistics for monitor at time=0
      compute_stats: block 
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
         use parallel, only: MPI_REAL_WP
         real(WP) :: mysgstke,myvisc,myTKE,myEPS,gu2buf
         real(WP) :: mysgsEPSp,myEPSp
         integer :: i,j,k,ierr

         myvisc=0.0_WP; myTKE=0.0_WP; myEPS=0.0_WP
         sgsTKE=0.0_WP; mysgstke=0.0_WP; 
         mysgsEPSp=0.0_WP; myEPSp=0.0_WP

         call fs%get_strainrate(SR=SR)
         call fs%get_gradu(gradu=gradu)

         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  ! LES stats
                  mysgsTKE = mysgsTKE + fs%cfg%vol(i,j,k)*(sgs%visc(i,j,k)/0.067_WP/sgs%delta(i,j,k)/fs%rho)**2
                  myvisc = myvisc+fs%visc(i,j,k)*fs%cfg%vol(i,j,k)

                  ! Resolved TKE
                  myTKE = myTKE+0.5_WP*((Ui(i,j,k)-meanU)**2+(Vi(i,j,k)-meanV)**2+(Wi(i,j,k)-meanW)**2)*fs%cfg%vol(i,j,k)

                  ! Dissipation
                  myEPS = myEPS + 2.0_WP*fs%visc(i,j,k)*fs%cfg%vol(i,j,k)*(SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)**2+2.0_WP*(SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2))/fs%rho

                  ! Pseudo-dissipation (sgs, total)
                  gu2buf = gradu(1,1,i,j,k)**2+gradu(1,2,i,j,k)**2+gradu(1,3,i,j,k)**2 + &
                           gradu(2,1,i,j,k)**2+gradu(2,2,i,j,k)**2+gradu(2,3,i,j,k)**2 + &
                           gradu(3,1,i,j,k)**2+gradu(3,2,i,j,k)**2+gradu(3,3,i,j,k)**2

                  mysgsEPSp=mysgsEPSp+fs%cfg%vol(i,j,k)*sgs%visc(i,j,k)*gu2buf
                  myEPSp=myEPSp+fs%cfg%vol(i,j,k)*fs%visc(i,j,k)*gu2buf
               end do
            end do
         end do

         call MPI_ALLREDUCE(myvisc,meanvisc,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanvisc=meanvisc/fs%cfg%vol_total

         call MPI_ALLREDUCE(mysgstke,sgsTKE,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); sgsTKE=sgsTKE/fs%cfg%vol_total
         call MPI_ALLREDUCE(myTKE,TKE,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr      ); TKE=TKE/fs%cfg%vol_total

         call MPI_ALLREDUCE(myEPS,EPS,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr        ); EPS=EPS/fs%cfg%vol_total
         call MPI_ALLREDUCE(myEPSp,EPSp,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr      ); EPSp=EPSp/fs%cfg%vol_total
         call MPI_ALLREDUCE(mysgsEPSp,sgsEPSp,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); sgsEPSp=sgsEPSp/fs%cfg%vol_total

         URMS = sqrt(2.0_WP/3.0_WP*(TKE+sgsTKE))
         Re_L = (TKE+sgsTKE)**2.0_WP/EPS/meanvisc 
         Re_lambda = sqrt(20.0_WP*Re_L/3.0_WP)
         eta = (meanvisc**3.0_WP/EPS)**0.25_WP
         ell = (0.6667_WP*TKE)**1.5_WP / EPS
         tau_eta = sqrt(visc/EPS0)
         
         nondtime  = time%t/tauinf
         dx_eta    = dx/eta
         eps_ratio = EPS/EPS0
         tke_ratio = TKE/TKE0
         ell_Lx    = ell/Lx
         Re_ratio  = Re_lambda/Re_max
      end block compute_stats
      
      ! Initialize LPT solver
      initialize_lpt: block
         use random, only: random_uniform, random_normal
         real(WP) :: dp
         integer :: i,np
         character(len=str_medium) :: timestamp
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')
         ! Get drag model from the inpit
         call param_read('Drag model',lp%drag_model,default='Schiller-Naumann')
         ! Get particle density from the input
         call param_read('Particle density',lp%rho)
         ! Get particle diameter from the input
         call param_read('Particle Stokes number', stk)
         ! Get number of particles
         call param_read('Number of particles',np)
         ! Check if a stochastic SGS model is used
         call param_read('Use CRW', use_crw)
         if (.false.) then
            call param_read('Restart from',timestamp,'r')
            ! Read the part file
            call lp%read(filename='restart/part_'//trim(adjustl(timestamp)))
         else
         ! Root process initializes np particles randomly
            if (lp%cfg%amRoot) then
               dp = sqrt(18.0_WP*visc*stk*tau_eta/lp%rho)
               call lp%resize(np)
               do i=1,np
                  ! Give id
                  lp%p(i)%id=int(i,8)
                  ! Set the diameter
                  lp%p(i)%d=dp
                  ! Assign random position in the domain
                  lp%p(i)%pos=[random_uniform(lp%cfg%x(lp%cfg%imin),lp%cfg%x(lp%cfg%imax+1)),&
                  &            random_uniform(lp%cfg%y(lp%cfg%jmin),lp%cfg%y(lp%cfg%jmax+1)),&
                  &            random_uniform(lp%cfg%z(lp%cfg%kmin),lp%cfg%z(lp%cfg%kmax+1))]
                  ! Give zero velocity
                  lp%p(i)%vel=0.0_WP
                  ! OU
                  lp%p(i)%uf=0.0_WP
                  if (use_crw) then
                     lp%p(i)%uf= [random_normal(m=0.0_WP,sd=0.1_WP),& 
                                  random_normal(m=0.0_WP,sd=0.1_WP),&    
                                  random_normal(m=0.0_WP,sd=0.1_WP)] 
                  end if
                  ! Give zero dt
                  lp%p(i)%dt=0.0_WP
                  lp%p(i)%a_crw=0.0_WP
                  lp%p(i)%b_crw=0.0_WP
                  ! Locate the particle on the mesh
                  lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
                  ! Activate the particle
                  lp%p(i)%flag=0
               end do
            end if
            ! Distribute particles
            call lp%sync()
         end if
      end block initialize_lpt
      

      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: np,i
         call param_read('Number of particles',np)
         pmesh=partmesh(nvar=1,nvec=3,name='lpt')
         pmesh%varname(1)="id"
         pmesh%vecname(1)="vel"
         pmesh%vecname(2)="fld"
         pmesh%vecname(3)="uf"
         call lp%resize(np)
         call lp%update_partmesh(pmesh)
         do i = 1,lp%np_
            pmesh%var(1,i) = lp%p(i)%id
            pmesh%vec(:,1,i) = lp%p(i)%vel
            pmesh%vec(:,2,i) = lp%cfg%get_velocity(pos=lp%p(i)%pos,i0=lp%p(i)%ind(1),j0=lp%p(i)%ind(2),k0=lp%p(i)%ind(3),U=Ui,V=Vi,W=Wi)
            pmesh%vec(:,3,i) = lp%p(i)%uf 
         end do
      end block create_pmesh
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='HIT')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('visc',sgs%visc)
         call ens_out%add_particle('particles',pmesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(meanU,'Umean')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(meanV,'Vmean')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(meanW,'Wmean')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(TKE,'Kinetic energy')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
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
         call cflfile%write()
         ! Create hit monitor
         hitfile=monitor(fs%cfg%amRoot,'hit')
         call hitfile%add_column(time%n,'Timestep number')
         call hitfile%add_column(time%t,'Time')
         call hitfile%add_column(Re_L,'Re_L')
         call hitfile%add_column(Re_lambda,'Re_lambda')
         call hitfile%add_column(meanvisc,'mean visc') 
         call hitfile%add_column(eta,'eta')
         call hitfile%add_column(sgsTKE,'sgsTKE')
         call hitfile%add_column(TKE,'TKE')
         call hitfile%add_column(URMS,'URMS')
         call hitfile%add_column(EPS,'EPS')
         call hitfile%add_column(EPSp,'EPSp')
         call hitfile%add_column(ell,'Integral lengthscale')
         call hitfile%write()
         ! Create hit convergence monitor
         ssfile=monitor(fs%cfg%amRoot,'convergence')
         call ssfile%add_column(time%n,'Timestep number')
         call ssfile%add_column(time%t,'Time')
         call ssfile%add_column(nondtime,'Time/t_int')
         call ssfile%add_column(Re_ratio,'Re_ratio')
         call ssfile%add_column(eps_ratio,'EPS_ratio')
         call ssfile%add_column(tke_ratio,'TKE_ratio')
         call ssfile%add_column(dx_eta,'dx/eta')
         call ssfile%add_column(ell_Lx,'ell/Lx')
         call ssfile%write()
         ! Create LPT monitor
         call lp%get_max()
         lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
         call lptfile%add_column(time%n,'Timestep number')
         call lptfile%add_column(time%t,'Time')
         call lptfile%add_column(lp%np,'Particle number')
         call lptfile%add_column(lp%Umin,'Particle Umin')
         call lptfile%add_column(lp%Umax,'Particle Umax')
         call lptfile%add_column(lp%Vmin,'Particle Vmin')
         call lptfile%add_column(lp%Vmax,'Particle Vmax')
         call lptfile%add_column(lp%Wmin,'Particle Wmin')
         call lptfile%add_column(lp%Wmax,'Particle Wmax')
         call lptfile%add_column(lp%dmin,'Particle dmin')
         call lptfile%add_column(lp%dmax,'Particle dmax')
         call lptfile%write()
         ! Create SGS monitor
         sgsfile=monitor(fs%cfg%amroot,'sgs')
         call sgsfile%add_column(time%n,'Timestep number')
         call sgsfile%add_column(time%t,'Time')
         call sgsfile%add_column(sgs%min_visc,'Min eddy visc')
         call sgsfile%add_column(sgs%max_visc,'Max eddy visc')
         call sgsfile%add_column(sgsEPSp,'SGSEPS')
         call sgsfile%add_column(sgsTKE,'sgsTKE')
         call sgsfile%write()
         ! Create timing monitor
         tfile=monitor(amroot=fs%cfg%amRoot,name='timing')
         call tfile%add_column(time%n,'Timestep number')
         call tfile%add_column(time%t,'Time')
         call tfile%add_column(wt_total%time,'Total [s]')
         call tfile%add_column(wt_vel%time,'Velocity [s]')
         call tfile%add_column(wt_vel%percent,'Velocity [%]')
         call tfile%add_column(wt_pres%time,'Pressure [s]')
         call tfile%add_column(wt_pres%percent,'Pressure [%]')
         call tfile%add_column(wt_lpt%time,'LPT [s]')
         call tfile%add_column(wt_lpt%percent,'LPT [%]')
         call tfile%add_column(wt_sgs%time,'SGS [s]')
         call tfile%add_column(wt_sgs%percent,'SGS [%]')
         call tfile%add_column(wt_stat%time,'Stats [s]')
         call tfile%add_column(wt_stat%percent,'Stats [%]')
         call tfile%add_column(wt_force%time,'Forcing [s]')
         call tfile%add_column(wt_force%percent,'Forcing [%]')
         call tfile%add_column(wt_rest%time,'Rest [s]')
         call tfile%add_column(wt_rest%percent,'Rest [%]')
         call tfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Time integrate our problem
   subroutine simulation_run
      use parallel,       only: parallel_time
      implicit none
      integer :: ii
      
      ! Perform time integration
      do while (.not.time%done())

         ! init wallclock
         wt_total%time_in=parallel_time()
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()


         wt_lpt%time_in=parallel_time()
         ! Advance particles by dt
         resU=fs%rho; resV=fs%visc
         call lp%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=resU,visc=resV,use_crw=use_crw,sgs=sgs)
         wt_lpt%time=wt_lpt%time+parallel_time()-wt_lpt%time_in
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W

         ! Turbulence modeling
         call fs%interp_vel(Ui,Vi,Wi)
         wt_sgs%time_in=parallel_time()
         if (use_sgs) then
            fs%visc=visc
            call fs%get_strainrate(SR=SR)
            call fs%get_gradu(gradu=gradu)
            resU=fs%rho
            call sgs%get_visc(type=sgs_type,rho=resU,gradu=gradu,dt=time%dt,SR=SR,Ui=Ui,Vi=Vi,Wi=Wi)
            where (sgs%visc.lt.-fs%visc)
               sgs%visc=-fs%visc
            end where
            fs%visc=fs%visc+sgs%visc
         end if
         wt_sgs%time=wt_sgs%time+parallel_time()-wt_sgs%time_in
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            wt_vel%time_in=parallel_time()
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dt*resU
            resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dt*resV
            resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dt*resW
            
            wt_vel%time=wt_vel%time+parallel_time()-wt_vel%time_in
            wt_force%time_in=parallel_time()
            if (linforce) then
               ! Add linear forcing term
               ! See Bassenne et al. (2016)
               linear_forcing: block
                  use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
                  use parallel, only: MPI_REAL_WP
                  use messager, only: die
                  real(WP) :: myTKE,A,myvisc,myEPSp
                  integer :: i,j,k,ierr
   
                  ! Calculate mean velocity
                  call fs%cfg%integrate(A=fs%U,integral=meanU); meanU=meanU/fs%cfg%vol_total
                  call fs%cfg%integrate(A=fs%V,integral=meanV); meanV=meanV/fs%cfg%vol_total
                  call fs%cfg%integrate(A=fs%W,integral=meanW); meanW=meanW/fs%cfg%vol_total
   
                  ! Calculate TKE and EPS
                  call fs%interp_vel(Ui,Vi,Wi)
                  call fs%get_gradu(gradu=gradu)
                  myvisc=0.0_WP; myTKE=0.0_WP; myEPSp=0.0_WP
                  do k=fs%cfg%kmin_,fs%cfg%kmax_
                     do j=fs%cfg%jmin_,fs%cfg%jmax_
                        do i=fs%cfg%imin_,fs%cfg%imax_  
                           myTKE=myTKE+0.5_WP*((Ui(i,j,k)-meanU)**2+(Vi(i,j,k)-meanV)**2+(Wi(i,j,k)-meanW)**2)*fs%cfg%vol(i,j,k)

                           ! Pseudo-dissipation 
                           myEPSp=myEPSp+fs%cfg%vol(i,j,k)*fs%visc(i,j,k)*(                       &
                                    gradu(1,1,i,j,k)**2+gradu(1,2,i,j,k)**2+gradu(1,3,i,j,k)**2 + &
                                    gradu(2,1,i,j,k)**2+gradu(2,2,i,j,k)**2+gradu(2,3,i,j,k)**2 + &
                                    gradu(3,1,i,j,k)**2+gradu(3,2,i,j,k)**2+gradu(3,3,i,j,k)**2)
                        end do
                     end do
                  end do
                  call MPI_ALLREDUCE(myTKE,TKE,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr      ); TKE = TKE/fs%cfg%vol_total
                  call MPI_ALLREDUCE(myEPSp,EPSp,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr    ); EPSp = EPSp/fs%cfg%vol_total/fs%rho
   
                  if (Gdtaui.lt.time%dt) call die("[linear_forcing] Controller time constant less than timestep")
                  A   = (EPSp - Gdtau*(TKE-TKE0))/(2.0_WP*TKE)*fs%rho ! - Eq. (7) (forcing constant TKE)
   
                  resU=resU+time%dt*(fs%U-meanU)*A
                  resV=resV+time%dt*(fs%V-meanV)*A
                  resW=resW+time%dt*(fs%W-meanW)*A
               end block linear_forcing
            end if
            wt_force%time=wt_force%time+parallel_time()-wt_force%time_in
            wt_vel%time_in=parallel_time()

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU/fs%rho
            fs%V=2.0_WP*fs%V-fs%Vold+resV/fs%rho
            fs%W=2.0_WP*fs%W-fs%Wold+resW/fs%rho
            
            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)
            wt_vel%time=wt_vel%time+parallel_time()-wt_vel%time_in
            
            ! Solve Poisson equation
            wt_pres%time_in=parallel_time()
            call fs%correct_mfr()
            call fs%get_div()
            fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho
            fs%V=fs%V-time%dt*resV/fs%rho
            fs%W=fs%W-time%dt*resW/fs%rho
            wt_pres%time=wt_pres%time+parallel_time()-wt_pres%time_in
            
            ! Increment sub-iteration counter
            time%it=time%it+1
   
         end do
         
         wt_vel%time_in=parallel_time()
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         wt_vel%time=wt_vel%time+parallel_time()-wt_vel%time_in

         wt_stat%time_in=parallel_time()
         ! Compute turbulence statistics for monitor
         compute_stats: block 
            use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
            use parallel, only: MPI_REAL_WP
            real(WP) :: mysgstke,myvisc,myTKE,myEPS
            real(WP) :: mysgsEPSp,myEPSp
            integer :: i,j,k,ierr

            myvisc=0.0_WP; myTKE=0.0_WP; myEPS=0.0_WP
            sgsTKE=0.0_WP; mysgstke=0.0_WP; 
            mysgsEPSp=0.0_WP

            call fs%get_strainrate(SR=SR)
            call fs%get_gradu(gradu=gradu)

            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do j=fs%cfg%jmin_,fs%cfg%jmax_
                  do i=fs%cfg%imin_,fs%cfg%imax_
                     ! LES stats
                     mysgsTKE = mysgsTKE + fs%cfg%vol(i,j,k)*(sgs%visc(i,j,k)/0.067_WP/sgs%delta(i,j,k)/fs%rho)**2
                     myvisc = myvisc+fs%visc(i,j,k)*fs%cfg%vol(i,j,k)

                     ! Resolved TKE
                     myTKE = myTKE+0.5_WP*((Ui(i,j,k)-meanU)**2+(Vi(i,j,k)-meanV)**2+(Wi(i,j,k)-meanW)**2)*fs%cfg%vol(i,j,k)

                     ! Dissipation
                     myEPS = myEPS + 2.0_WP*fs%visc(i,j,k)*fs%cfg%vol(i,j,k)*(SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)**2+2.0_WP*(SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2))/fs%rho

                     ! Pseudo-dissipation (sgs, resolved, total)
                     mysgsEPSp=mysgsEPSp+fs%cfg%vol(i,j,k)*sgs%visc(i,j,k)*(                   &
                                 gradu(1,1,i,j,k)**2+gradu(1,2,i,j,k)**2+gradu(1,3,i,j,k)**2 + &
                                 gradu(2,1,i,j,k)**2+gradu(2,2,i,j,k)**2+gradu(2,3,i,j,k)**2 + &
                                 gradu(3,1,i,j,k)**2+gradu(3,2,i,j,k)**2+gradu(3,3,i,j,k)**2)
                  end do
               end do
            end do

            call MPI_ALLREDUCE(myvisc,meanvisc,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanvisc=meanvisc/fs%cfg%vol_total

            call MPI_ALLREDUCE(mysgstke,sgsTKE,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); sgsTKE=sgsTKE/fs%cfg%vol_total
            call MPI_ALLREDUCE(myTKE,TKE,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr      ); TKE=TKE/fs%cfg%vol_total

            call MPI_ALLREDUCE(myEPS,EPS,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr        ); EPS=EPS/fs%cfg%vol_total
            call MPI_ALLREDUCE(mysgsEPSp,sgsEPSp,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); sgsEPSp=sgsEPSp/fs%cfg%vol_total

            URMS = sqrt(2.0_WP/3.0_WP*(TKE+sgsTKE))
            Re_L = (TKE+sgsTKE)**2.0_WP/EPS/meanvisc 
            Re_lambda = sqrt(20.0_WP*Re_L/3.0_WP)
            eta = (meanvisc**3.0_WP/EPS)**0.25_WP
            ell = (0.6667_WP*TKE)**1.5_WP / EPS
            
            nondtime  = time%t/tauinf
            dx_eta    = dx/eta
            eps_ratio = EPS/EPS0
            tke_ratio = TKE/TKE0
            ell_Lx    = ell/Lx
            Re_ratio  = Re_lambda/Re_max
         end block compute_stats
         wt_stat%time=wt_stat%time+parallel_time()-wt_stat%time_in
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            call lp%update_partmesh(pmesh)
            do ii = 1,lp%np_
               pmesh%var(1,ii) = lp%p(ii)%id
               pmesh%vec(:,1,ii) = lp%p(ii)%vel
               pmesh%vec(:,2,ii) = lp%cfg%get_velocity(pos=lp%p(ii)%pos,i0=lp%p(ii)%ind(1),j0=lp%p(ii)%ind(2),k0=lp%p(ii)%ind(3),U=Ui,V=Vi,W=Wi)
               pmesh%vec(:,3,ii) = lp%p(ii)%uf
            end do
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call mfile%write()
         call cflfile%write()
         call hitfile%write()
         call lp%get_max()
         call lptfile%write()
         call sgsfile%write()
         call ssfile%write()

         ! Monitor timing
         wt_total%time=parallel_time()-wt_total%time_in
         wt_vel%percent=wt_vel%time/wt_total%time*100.0_WP
         wt_pres%percent=wt_pres%time/wt_total%time*100.0_WP
         wt_lpt%percent=wt_lpt%time/wt_total%time*100.0_WP
         wt_sgs%percent=wt_sgs%time/wt_total%time*100.0_WP
         wt_stat%percent=wt_stat%time/wt_total%time*100.0_WP
         wt_force%percent=wt_force%time/wt_total%time*100.0_WP
         wt_rest%time=wt_total%time-wt_vel%time-wt_pres%time-wt_lpt%time-wt_sgs%time-wt_stat%time-wt_force%time
         wt_rest%percent=wt_rest%time/wt_total%time*100.0_WP
         call tfile%write()
         wt_total%time=0.0_WP; wt_total%percent=0.0_WP
         wt_vel%time=0.0_WP;   wt_vel%percent=0.0_WP
         wt_pres%time=0.0_WP;  wt_pres%percent=0.0_WP
         wt_lpt%time=0.0_WP;   wt_lpt%percent=0.0_WP
         wt_sgs%time=0.0_WP;   wt_sgs%percent=0.0_WP
         wt_stat%time=0.0_WP;  wt_stat%percent=0.0_WP
         wt_force%time=0.0_WP; wt_force%percent=0.0_WP
         wt_rest%time=0.0_WP;  wt_rest%percent=0.0_WP

         ! Finally, see if it's time to save restart files
         if (save_evt%occurs()) then
            save_restart: block
              character(len=str_medium) :: timestamp
              ! Prefix for files
              write(timestamp,'(es12.5)') time%t
              ! Prepare a new directory
              if (fs%cfg%amRoot) call execute_command_line('mkdir -p restart')
              ! Populate df and write it
              call df%pushval(name='t' ,val=time%t     )
              call df%pushval(name='dt',val=time%dt    )
              call df%pushvar(name='U' ,var=fs%U       )
              call df%pushvar(name='V' ,var=fs%V       )
              call df%pushvar(name='W' ,var=fs%W       )
              call df%pushvar(name='P' ,var=fs%P       )
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
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi,SR,gradu)
      
   end subroutine simulation_final
   
   
   
   
   
end module simulation
