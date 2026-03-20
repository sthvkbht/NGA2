!> Various definitions and tools for running an NGA2 simulation
module simulation
   use string,            only: str_medium
   use precision,         only: WP
   use geometry,          only: cfg
   use spcomp_class,      only: spcomp
   use lpt_class,         only: lpt
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use datafile_class,    only: datafile
   implicit none
   private
   
   !> Get a couple linear solvers, an incompressible flow solver and corresponding time tracker
   type(spcomp),      public :: fs
   type(lpt),         public :: lp
   type(timetracker), public :: time

   !> Provide a datafile and an event tracker for saving restarts
   type(event)    :: save_evt
   type(datafile) :: df
   logical :: restarted
   
   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,consfile,lptfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: dQdt
   real(WP), dimension(:,:,:)    , allocatable :: Ui,Vi,Wi,Ma,beta,visc,visc_t,div
   real(WP), dimension(:,:,:)    , allocatable :: srcUlp,srcVlp,srcWlp,srcIlp
   real(WP), dimension(:,:,:)    , allocatable :: stressx,stressy,stressz,stressI

   !> Equations of state
   real(WP) :: Pinf,Gamma,Cv,Prandtl,Tib

   !> Flow parameters
   real(WP) :: U0,T0,rhoU0,meanRhoU,meanT
   
 contains


   !> P=EOS(RHO,I)
   pure real(WP) function get_P(RHO,I)
     implicit none
     real(WP), intent(in) :: RHO,I
     get_P=RHO*I*(Gamma-1.0_WP)-Gamma*Pinf
   end function get_P
   !> T=f(RHO,P)
   pure real(WP) function get_T(RHO,P)
     implicit none
     real(WP), intent(in) :: RHO,P
     get_T=(P+Pinf)/(Cv*RHO*(Gamma-1.0_WP))
   end function get_T
   !> RHO = f(P,T)
   pure real(WP) function get_RHO(P,T)
     implicit none
     real(WP), intent(in) :: P,T
     get_RHO=(P+Pinf)/(Cv*T*(Gamma-1.0_WP))
   end function get_RHO
   !> I=EOS(RHO,P)
   pure real(WP) function get_I(RHO,P)
     implicit none
     real(WP), intent(in) :: RHO,P
     get_I=(P+Gamma*Pinf)/(RHO*(Gamma-1.0_WP))
   end function get_I
   !> C=f(RHO,P)
   pure real(WP) function get_C(RHO,P)
     implicit none
     real(WP), intent(in) :: RHO,P
     get_C=sqrt(Gamma*(P+Pinf)/RHO)
   end function get_C
   !> S=f(RHO,P)
   pure real(WP) function get_S(RHO,P)
     implicit none
     real(WP), intent(in) :: RHO,P
     get_S=Cv*log((P+Pinf)/RHO**Gamma)
   end function get_S


   !> Calculate viscosities
   subroutine prepare_viscosities()
     implicit none
     integer :: i,j,k
     real(WP), parameter :: T0=350.0_WP
     real(WP), parameter :: S=1064.0_WP
     ! Get viscosity from Sutherland's law for water vapor
     do k=fs%cfg%kmino_,fs%cfg%kmaxo_
        do j=fs%cfg%jmino_,fs%cfg%jmaxo_
           do i=fs%cfg%imino_,fs%cfg%imaxo_
              visc(i,j,k)=1.12e-5_WP*(fs%T(i,j,k)/T0)**1.5_WP*(T0+S)/(fs%T(i,j,k)+S)
           end do
        end do
     end do
     ! Get LAD
     call fs%get_viscartif(dt=time%dt,beta=beta); fs%BETA=fs%Q(:,:,:,1)*beta
     ! Get eddy viscosity
     call fs%get_vreman   (dt=time%dt,visc=visc_t); fs%VISC=fs%Q(:,:,:,1)*visc_t+visc
     ! Recompute thermal conductivity
     fs%diff=Gamma*Cv*fs%visc/Prandtl
   end subroutine prepare_viscosities


   !> Calculate velocity divergence
   subroutine get_div()
     implicit none
     integer :: i,j,k
     do k=fs%cfg%kmino_,fs%cfg%kmaxo_-1; do j=fs%cfg%jmino_,fs%cfg%jmaxo_-1; do i=fs%cfg%imino_,fs%cfg%imaxo_-1
        div(i,j,k)=fs%dxi*(fs%U(i+1,j,k)-fs%U(i,j,k))+fs%dyi*(fs%V(i,j+1,k)-fs%V(i,j,k))+fs%dzi*(fs%W(i,j,k+1)-fs%W(i,j,k))
     end do; end do; end do
     call fs%cfg%sync(div)
     if (.not.fs%cfg%xper.and.fs%cfg%iproc.eq.fs%cfg%npx) div(fs%cfg%imaxo,:,:)=div(fs%cfg%imaxo-1,:,:)
     if (.not.fs%cfg%yper.and.fs%cfg%jproc.eq.fs%cfg%npy) div(:,fs%cfg%jmaxo,:)=div(:,fs%cfg%jmaxo-1,:)
     if (.not.fs%cfg%zper.and.fs%cfg%kproc.eq.fs%cfg%npz) div(:,:,fs%cfg%kmaxo)=div(:,:,fs%cfg%kmaxo-1)
   end subroutine get_div


   !> Compute mean momentum and temperature
   subroutine get_bodyforce()
     use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE,MPI_IN_PLACE
     use parallel, only: MPI_REAL_WP
     integer :: i,j,k,ierr
     real(WP) :: vol,Uvol
     Uvol=0.0_WP; meanRhoU=0.0_WP; meanT=0.0_WP
     do k=fs%cfg%kmin_,fs%cfg%kmax_
        do j=fs%cfg%jmin_,fs%cfg%jmax_
           do i=fs%cfg%imin_,fs%cfg%imax_
              vol=fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*0.5_WP*sum(cfg%VF(i-1:i,j,k))
              Uvol=Uvol+vol
              meanRhoU=meanRhoU+vol*fs%Q(i,j,k,3)
              meanT=meanT+vol*get_T(fs%Q(i,j,k,1),fs%P(i,j,k))
           end do
        end do
     end do
     call MPI_ALLREDUCE(MPI_IN_PLACE,Uvol,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,meanRhoU,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanRhoU=meanRhoU/Uvol
     call MPI_ALLREDUCE(MPI_IN_PLACE,meanT,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanT=meanT/Uvol
   end subroutine get_bodyforce


   !> Overwrite cosnerved variables using volume-of-solid IBM
   subroutine apply_ibm()
     implicit none
     integer :: i,j,k,ii,jj,kk
     real(WP) :: sum_VF,sum_VFP
     call fs%get_primitive(1.0_WP-lp%VF)
     do k=cfg%kmin_,cfg%kmax_
        do j=cfg%jmin_,cfg%jmax_
           do i=cfg%imin_,cfg%imax_
              ! No-slip
              fs%Q(i,j,k,3)=0.5_WP*(cfg%VF(i-1,j,k)+cfg%VF(i,j,k))*fs%Q(i,j,k,3)
              fs%Q(i,j,k,4)=0.5_WP*(cfg%VF(i,j-1,k)+cfg%VF(i,j,k))*fs%Q(i,j,k,4)
              fs%Q(i,j,k,5)=0.5_WP*(cfg%VF(i,j,k-1)+cfg%VF(i,j,k))*fs%Q(i,j,k,5)
              if (cfg%VF(i,j,k).gt.1.0_WP-epsilon(1.0_WP)) cycle
              ! Isothermal
              fs%T(i,j,k)=cfg%VF(i,j,k)*fs%T(i,j,k)+(1.0_WP-cfg%VF(i,j,k))*Tib
              ! Neumann: VF-weighted neighbor average for pressure
              sum_VF=0.0_WP; sum_VFP=0.0_WP
              do kk=-1,1; do jj=-1,1; do ii=-1,1
                 if (ii.eq.0.and.jj.eq.0.and.kk.eq.0) cycle
                 sum_VF =sum_VF +cfg%VF(i+ii,j+jj,k+kk)
                 sum_VFP=sum_VFP+cfg%VF(i+ii,j+jj,k+kk)*fs%P(i+ii,j+jj,k+kk)
              end do; end do; end do
              if (sum_VF.gt.0.0_WP) then
                 fs%P(i,j,k)=cfg%VF(i,j,k)*fs%P(i,j,k)+(1.0_WP-cfg%VF(i,j,k))*sum_VFP/sum_VF
              end if
              ! Reconstruct other variables
              fs%Q(i,j,k,1)=get_RHO(fs%P(i,j,k),fs%T(i,j,k))
              fs%I(i,j,k)=get_I(fs%Q(i,j,k,1),fs%P(i,j,k))
              fs%Q(i,j,k,1)=fs%Q(i,j,k,1)*(1.0_WP-lp%VF(i,j,k))
              fs%Q(i,j,k,2)=fs%Q(i,j,k,1)*fs%I(i,j,k)
           end do
        end do
     end do
     ! Communicate
     call fs%cfg%sync(fs%Q(:,:,:,1))
     call fs%cfg%sync(fs%Q(:,:,:,2))
     call fs%cfg%sync(fs%Q(:,:,:,3))
     call fs%cfg%sync(fs%Q(:,:,:,4))
     call fs%cfg%sync(fs%Q(:,:,:,5))
     ! Recompute primitive variables
     call fs%get_primitive(1.0_WP-lp%VF)
   end subroutine apply_ibm
   

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      
      ! Create compressible flow solver
      create_flow_solver: block
        call fs%initialize(cfg=cfg,name='Compressible NS')
      end block create_flow_solver


      ! Allocate work arrays
      allocate_work_arrays: block
        allocate(dQdt   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:fs%nQ,1:4))
        allocate(Ui     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Vi     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Wi     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Ma     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(beta   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(visc   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(visc_t(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(div    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(stressx(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(stressy(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(stressz(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(stressI(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(srcUlp (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(srcVlp (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(srcWlp (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(srcIlp (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays


      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
        time=timetracker(amRoot=cfg%amRoot)
        call param_read('Max timestep size',time%dtmax)
        call param_read('Max cfl number',time%cflmax)
        call param_read('Max time',time%tmax)
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
           df=datafile(pg=cfg,filename=trim(cfg%name),nval=4,nvar=5)
           df%valname(1)='t'
           df%valname(2)='dt'
           df%valname(3)='meanRhoU'
           df%valname(4)='meanT'
           df%varname(1)='Q1'
           df%varname(2)='Q2'
           df%varname(3)='Q3'
           df%varname(4)='Q4'
           df%varname(5)='Q5'
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


      ! Initialize our LPT solver
      initialize_lpt: block
        use random, only: random_lognormal,random_uniform
        use mathtools, only: Pi,twoPi
        use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE,MPI_INTEGER
        use parallel, only: MPI_REAL_WP
        character(len=str_medium) :: timestamp
        real(WP) :: VFavg,Vol_,sumVolp,dp,Tp,buf
        integer :: i,j,k,ii,jj,kk,nn,ip,jp,kp,np,offset,ierr
        integer, dimension(:,:,:), allocatable :: npic      !< Number of particle in cell
        integer, dimension(:,:,:,:), allocatable :: ipic    !< Index of particle in cell
        logical :: overlap,outside
        ! Create solver
        lp=lpt(cfg=cfg,name='LPT')
        ! Get mean volume fraction from input
        call param_read('Particle volume fraction',VFavg)
        ! Set gravity
        call param_read('Gravity',lp%gravity)
        ! Read in bulk velocity
        call param_read('Bulk velocity',U0)
        ! Get particle density from input
        call param_read('Particle density',lp%rho)
        ! Get particle specific heat from input
        call param_read('Particle heat capacity',lp%Cp)
        ! Get particle diameter from input
        call param_read('Particle diameter',dp)
        ! Get particle temperature from input
        call param_read('Particle temperature',Tp)
        ! Set collision timescale
        call param_read('Collision timescale',lp%tau_col,default=15.0_WP*time%dt)
        ! Set coefficient of restitution
        call param_read('Coefficient of restitution',lp%e_n)
        call param_read('Friction coefficient',lp%mu_f,default=0.0_WP)
        ! Set filter scale to 3.5*dx
        lp%filter_width=3.5_WP*cfg%min_meshsize
        ! Initialize particles
        if (restarted) then
           call param_read('Restart from',timestamp,'r')
           ! Read the part file
           call lp%read(filename='restart/part_'//trim(adjustl(timestamp)))
        else
           ! Get volume of domain belonging to this proc
           Vol_=0.0_WP
           do k=lp%cfg%kmin_,lp%cfg%kmax_
              do j=lp%cfg%jmin_,lp%cfg%jmax_
                 do i=lp%cfg%imin_,fs%cfg%imax_
                    Vol_=Vol_+lp%cfg%vol(i,j,k)*lp%cfg%VF(i,j,k)
                 end do
              end do
           end do
           ! Get particle diameters
           np=5*ceiling(VFavg*Vol_/(Pi*dp**3/6.0_WP))
           sumVolp=0.0_WP; np=0
           do while(sumVolp.lt.VFavg*Vol_)
              np=np+1
              sumVolp=sumVolp+Pi/6.0_WP*dp**3
           end do
           call lp%resize(np)
           ! Allocate particle in cell arrays
           allocate(npic(     lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); npic=0
           allocate(ipic(1:40,lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); ipic=0
           ! Distribute particles
           sumVolp=0.0_WP
           do i=1,np
              ! Set the diameter
              lp%p(i)%d=dp
              ! Give position (avoid overlap)
              overlap=.true.
              do while (overlap)
                 outside=.true.
                 do while (outside)
                    lp%p(i)%pos=[random_uniform(lp%cfg%x(lp%cfg%imin_),lp%cfg%x(lp%cfg%imax_+1)-dp),&
                         &       random_uniform(lp%cfg%y(lp%cfg%jmin_),lp%cfg%y(lp%cfg%jmax_+1)-dp),&
                         &       random_uniform(lp%cfg%z(lp%cfg%kmin_),lp%cfg%z(lp%cfg%kmax_+1)-dp)]
                    if (lp%cfg%nz.eq.1) lp%p(i)%pos(3)=0.0_WP
                    lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
                    buf=cfg%get_scalar(pos=lp%p(i)%pos,i0=lp%p(i)%ind(1),j0=lp%p(i)%ind(2),k0=lp%p(i)%ind(3),S=cfg%Gib,bc='n')
                    if (buf.lt.-3.0_WP*cfg%min_meshsize) outside=.false.
                 end do
                 overlap=.false.
                 do kk=lp%p(i)%ind(3)-1,lp%p(i)%ind(3)+1
                    do jj=lp%p(i)%ind(2)-1,lp%p(i)%ind(2)+1
                       do ii=lp%p(i)%ind(1)-1,lp%p(i)%ind(1)+1
                          do nn=1,npic(ii,jj,kk)
                             j=ipic(nn,ii,jj,kk)
                             if (sqrt(sum((lp%p(i)%pos-lp%p(j)%pos)**2)).lt.0.5_WP*(lp%p(i)%d+lp%p(j)%d)) overlap=.true.
                          end do
                       end do
                    end do
                 end do
              end do

              ! Activate the particle
              lp%p(i)%flag=0
              ip=lp%p(i)%ind(1); jp=lp%p(i)%ind(2); kp=lp%p(i)%ind(3)
              npic(ip,jp,kp)=npic(ip,jp,kp)+1
              ipic(npic(ip,jp,kp),ip,jp,kp)=i
              ! Give temperature
              lp%p(i)%T=Tp
              ! Give zero velocity
              buf=0.2_WP
              lp%p(i)%vel=[random_uniform(-buf*U0,buf*U0),&
                   &       random_uniform(-buf*U0,buf*U0),&
                   &       random_uniform(-buf*U0,buf*U0)]
              !lp%p(i)%vel=0.0_WP
              lp%p(i)%vel(1)=lp%p(i)%vel(1)+U0
              ! Give zero collision force
              lp%p(i)%Acol=0.0_WP
              lp%p(i)%Tcol=0.0_WP
              ! Sum up volume
              sumVolp=sumVolp+Pi/6.0_WP*lp%p(i)%d**3
           end do
           deallocate(npic,ipic)
           call lp%sync()
           ! Set ID
           offset=0
           do i=1,lp%cfg%rank
              offset=offset+lp%np_proc(i)
           end do
           do i=1,lp%np_
              lp%p(i)%id=int(i+offset,8)
           end do
           ! Get mean diameter and volume fraction
           call MPI_ALLREDUCE(sumVolp,VFavg,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); VFavg=VFavg/lp%cfg%fluid_vol
           if (lp%cfg%amRoot) then
              print*,"===== Particle Setup Description ====="
              print*,'Number of particles', lp%np
              print*,'Mean volume fraction',VFavg
           end if
        end if
        ! Get initial particle volume fraction
        call lp%update_VF()
      end block initialize_lpt


      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
        integer :: i
        pmesh=partmesh(nvar=2,nvec=1,name='lpt')
        pmesh%varname(1)='diameter'
        pmesh%varname(2)='temperature'
        pmesh%vecname(1)='velocity'
        call lp%update_partmesh(pmesh)
        do i=1,lp%np_
           pmesh%var(1,i)=lp%p(i)%d
           pmesh%var(2,i)=lp%p(i)%T
           pmesh%vec(:,1,i)=lp%p(i)%vel
        end do
      end block create_pmesh


      ! Initialize variables
      initialize_variables: block
        integer :: i,j,k
        real(WP) :: rho0,p0,Rgas
        ! Provide thermodynamic model
        fs%getP=>get_P; fs%getC=>get_C; fs%getS=>get_S; fs%getT=>get_T
        ! Set Pinf to zero
        Pinf=0.0_WP
        ! Assign acceleration of gravity
        call param_read('Gravity',fs%gravity)
        ! Assign gas constant
        call param_read('Gas constant',Rgas)
        ! Read in Gamma
        call param_read('Gamma',Gamma)
        ! Read in Prandtl number
        call param_read('Prandtl number',Prandtl)
        ! Assign temperature
        call param_read('Temperature',Tib)
        ! Assign pressure
        call param_read('Pressure',p0)
        ! Set heat capacities
        Cv=Rgas/(Gamma-1.0_WP)
        ! Get density
        rho0=get_RHO(p0,Tib)
        if (restarted) then
           call df%pullvar(name='Q1',var=fs%Q(:,:,:,1))
           call df%pullvar(name='Q2',var=fs%Q(:,:,:,2))
           call df%pullvar(name='Q3',var=fs%Q(:,:,:,3))
           call df%pullvar(name='Q4',var=fs%Q(:,:,:,4))
           call df%pullvar(name='Q5',var=fs%Q(:,:,:,5))
           call df%pullval(name='meanRhoU',val=meanRhoU)
           call df%pullval(name='meanT',val=meanT)
           rhoU0=meanRhoU; T0=meanT
        else
           ! Initialize primary variables
           do k=cfg%kmino_,cfg%kmaxo_
              do j=cfg%jmino_,cfg%jmaxo_
                 do i=cfg%imino_,cfg%imaxo_
                    fs%U(i,j,k)  =U0
                    fs%V(i,j,k)  =0.0_WP
                    fs%W(i,j,k)  =0.0_WP
                    fs%Q(i,j,k,1)=rho0
                    fs%P(i,j,k)  =p0
                    fs%I(i,j,k)  =get_I(fs%Q(i,j,k,1),fs%P(i,j,k))
                 end do
              end do
           end do
           ! Initialize conserved variables
           fs%Q(:,:,:,2)=fs%Q(:,:,:,1)*fs%I
           ! Multiply by volume fraction
           do i=1,fs%nQ
              fs%Q(:,:,:,i)=fs%Q(:,:,:,i)*(1.0_WP-lp%VF)
           end do
           call fs%get_momentum()
           ! Get target momentum and energy
           call get_bodyforce()
           rhoU0=meanRhoU; T0=meanT
           if (fs%cfg%amRoot) then
              print*,"===== Fluid Setup Description ====="
              print*,'Mach number', U0/maxval(fs%C)
              print*,'Density',rho0
              print*,'Pressure',p0
              print*,'Temperature',Tib
              print*, 'Cp',Gamma*Cv
              print*, 'Cv',Cv
              print*,"==================================="
           end if
        end if
        ! Rebuild primitive variables
        call fs%get_primitive(1.0_WP-lp%VF)
        ! Interpolate velocity
        call fs%interp_vel(Ui,Vi,Wi)
        ! Compute local Mach number
        Ma=sqrt(Ui**2+Vi**2+Wi**2)/fs%C
        ! Compute viscosities
        call prepare_viscosities()
        ! Compute dilatation
        call get_div()
        !> Perform and output monitoring
        call fs%get_info()
        call mfile%write()
        call cflfile%write()
        call consfile%write()
      end block initialize_variables


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='conduit')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_particle('particles',pmesh)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('P',fs%P)
         call ens_out%add_scalar('T',fs%T)
         call ens_out%add_scalar('Mach',Ma)
         call ens_out%add_scalar('beta',beta)
         call ens_out%add_scalar('visc',visc)
         call ens_out%add_scalar('visc_t',visc_t)
         call ens_out%add_scalar('div',div) 
         call ens_out%add_scalar('epsp',lp%VF)
         call ens_out%add_scalar('levelset',cfg%Gib)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create monitor files
      create_monitor: block
        real(WP) :: cfl
        call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
        call fs%get_cfl(time%dt,cfl); time%cfl=max(time%cfl,cfl)
        call fs%get_info()
        call lp%get_max()
        ! Create simulation monitor
        mfile=monitor(fs%cfg%amRoot,'simulation')
        call mfile%add_column(time%n,'Timestep number')
        call mfile%add_column(time%t,'Time')
        call mfile%add_column(time%dt,'Timestep size')
        call mfile%add_column(time%cfl,'Maximum CFL')
        call mfile%add_column(meanRhoU,'MFR')
        call mfile%add_column(meanT,'<T>')
        call mfile%add_column(fs%Umax,'Umax')
        call mfile%add_column(fs%Vmax,'Vmax')
        call mfile%add_column(fs%Wmax,'Wmax')
        call mfile%add_column(fs%RHOmax,'max(RHO)')
        call mfile%add_column(fs%RHOmin,'min(RHO)')
        call mfile%add_column(fs%Pmax  ,'max(P)'  )
        call mfile%add_column(fs%Pmin  ,'min(P)'  )
        call mfile%add_column(fs%Tmax  ,'max(T)'  )
        call mfile%add_column(fs%Tmin  ,'min(T)'  )
        call mfile%write()
        ! Create CFL monitor
        cflfile=monitor(fs%cfg%amRoot,'cfl')
        call cflfile%add_column(time%n,'Timestep number')
        call cflfile%add_column(time%t,'Time')
        call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
        call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
        call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
        call cflfile%add_column(fs%CFLa_x,'Acoustic xCFL')
        call cflfile%add_column(fs%CFLa_y,'Acoustic yCFL')
        call cflfile%add_column(fs%CFLa_z,'Acoustic zCFL')
        call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
        call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
        call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
        call cflfile%add_column(lp%CFLp_x,'Particle xCFL')
        call cflfile%add_column(lp%CFLp_y,'Particle yCFL')
        call cflfile%add_column(lp%CFLp_z,'Particle zCFL')
        call cflfile%add_column(lp%CFL_col,'Particle cCFL')
        call cflfile%write()
        ! Create conservation monitor
        consfile=monitor(fs%cfg%amRoot,'conservation')
        call consfile%add_column(time%n,'Timestep number')
        call consfile%add_column(time%t,'Time')
        call consfile%add_column(fs%Qint(1),'Mass')
        call consfile%add_column(fs%Qint(2),'Energy')
        call consfile%add_column(fs%Qint(3),'U Momentum')
        call consfile%add_column(fs%Qint(4),'V Momentum')
        call consfile%add_column(fs%Qint(5),'W Momentum')
        call consfile%add_column(fs%RHOKint,'Kinetic Energy')
        call consfile%add_column(fs%RHOSint,'Entropy')
        call consfile%write()
        ! Create LPT monitor
        lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
        call lptfile%add_column(time%n,'Timestep number')
        call lptfile%add_column(time%t,'Time')
        call lptfile%add_column(lp%np,'Particle number')
        call lptfile%add_column(lp%VFmean,'mean(VFp)')
        call lptfile%add_column(lp%VFmax,'max(VFp)')
        call lptfile%add_column(lp%VFvar,'var(VFp)')
        call lptfile%add_column(lp%Umin,'min(U)')
        call lptfile%add_column(lp%Umax,'max(U)')
        call lptfile%add_column(lp%Vmin,'min(V)')
        call lptfile%add_column(lp%Vmax,'max(V)')
        call lptfile%add_column(lp%Wmin,'min(W)')
        call lptfile%add_column(lp%Wmax,'max(W)')
        call lptfile%add_column(lp%Remax,'max(Re)')
        call lptfile%add_column(lp%Mamax,'max(Ma)')
        call lptfile%add_column(lp%Knmax,'max(Kn)')
        call lptfile%write()
      end block create_monitor

    end subroutine simulation_init


    !> Perform an NGA2 simulation
    subroutine simulation_run
      implicit none
      real(WP) :: cfl

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
         call fs%get_cfl(time%dt,cfl); time%cfl=max(time%cfl,cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember conserved variables
         fs%Qold=fs%Q

         ! Get mass flow rate
         call get_bodyforce()

         ! Prepare SGS viscosity models
         call prepare_viscosities()

         ! First RK step ====================================================================================
         ! Particle increment
         call fs%get_div_stress(stressx,stressy,stressz,stressI)
         call lp%collide(dt=time%dt,Gib=cfg%Gib,Nxib=cfg%Nib(1,:,:,:),Nyib=cfg%Nib(2,:,:,:),Nzib=cfg%Nib(3,:,:,:))
         !call lp%collide(dt=time%dt)
         call lp%substep_rk4(stage=1,dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=fs%rho,visc=fs%visc,T=fs%T,C=fs%C,&
              stress_x=stressx,stress_y=stressy,stress_z=stressz,heat_flux=stressI,srcU=srcUlp,srcV=srcVlp,srcW=srcWlp,srcI=srcIlp)
         ! Get non-SL RHS and increment
         call fs%rhs(VF=lp%VF,VFU=lp%VFU,VFV=lp%VFV,VFW=lp%VFW,dQdt=dQdt(:,:,:,:,1))
         ! LPT source
         dQdt(:,:,:,2,1)=dQdt(:,:,:,2,1)+srcIlp
         dQdt(:,:,:,3,1)=dQdt(:,:,:,3,1)+srcUlp
         dQdt(:,:,:,4,1)=dQdt(:,:,:,4,1)+srcVlp
         dQdt(:,:,:,5,1)=dQdt(:,:,:,5,1)+srcWlp
         ! Add body forcing to momentum and energy
         dQdt(:,:,:,2,1)=dQdt(:,:,:,2,1)+Cv*fs%Q(:,:,:,1)*(T0-meanT)/time%dt
         dQdt(:,:,:,3,1)=dQdt(:,:,:,3,1)+(rhoU0-meanRhoU)/time%dt
         ! Advance
         fs%Q=fs%Qold+0.5_WP*time%dt*dQdt(:,:,:,:,1)
         ! Apply IBM
         call apply_ibm()

         ! Second RK step ===================================================================================
         ! Particle increment
         call fs%get_div_stress(stressx,stressy,stressz,stressI)
         call lp%substep_rk4(stage=2,dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=fs%rho,visc=fs%visc,T=fs%T,C=fs%C,&
              stress_x=stressx,stress_y=stressy,stress_z=stressz,heat_flux=stressI,srcU=srcUlp,srcV=srcVlp,srcW=srcWlp,srcI=srcIlp)
         ! Get non-SL RHS and increment
         call fs%rhs(VF=lp%VF,VFU=lp%VFU,VFV=lp%VFV,VFW=lp%VFW,dQdt=dQdt(:,:,:,:,2))
         ! LPT source
         dQdt(:,:,:,2,2)=dQdt(:,:,:,2,2)+srcIlp
         dQdt(:,:,:,3,2)=dQdt(:,:,:,3,2)+srcUlp
         dQdt(:,:,:,4,2)=dQdt(:,:,:,4,2)+srcVlp
         dQdt(:,:,:,5,2)=dQdt(:,:,:,5,2)+srcWlp
         ! Add body forcing to momentum and energy
         dQdt(:,:,:,2,2)=dQdt(:,:,:,2,2)+Cv*fs%Q(:,:,:,1)*(T0-meanT)/time%dt
         dQdt(:,:,:,3,2)=dQdt(:,:,:,3,2)+(rhoU0-meanRhoU)/time%dt
         ! Advance
         fs%Q=fs%Qold+0.5_WP*time%dt*dQdt(:,:,:,:,2)
         ! Apply IBM
         call apply_ibm()

         ! Third RK step ====================================================================================
         ! Particle increment
         call fs%get_div_stress(stressx,stressy,stressz,stressI)
         call lp%substep_rk4(stage=3,dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=fs%rho,visc=fs%visc,T=fs%T,C=fs%C,&
              stress_x=stressx,stress_y=stressy,stress_z=stressz,heat_flux=stressI,srcU=srcUlp,srcV=srcVlp,srcW=srcWlp,srcI=srcIlp)
         ! Get non-SL RHS and increment
         call fs%rhs(VF=lp%VF,VFU=lp%VFU,VFV=lp%VFV,VFW=lp%VFW,dQdt=dQdt(:,:,:,:,3))
         ! LPT source
         dQdt(:,:,:,2,3)=dQdt(:,:,:,2,3)+srcIlp
         dQdt(:,:,:,3,3)=dQdt(:,:,:,3,3)+srcUlp
         dQdt(:,:,:,4,3)=dQdt(:,:,:,4,3)+srcVlp
         dQdt(:,:,:,5,3)=dQdt(:,:,:,5,3)+srcWlp
         ! Add body forcing to momentum and energy
         dQdt(:,:,:,2,3)=dQdt(:,:,:,2,3)+Cv*fs%Q(:,:,:,1)*(T0-meanT)/time%dt
         dQdt(:,:,:,3,3)=dQdt(:,:,:,3,3)+(rhoU0-meanRhoU)/time%dt
         ! Advance
         fs%Q=fs%Qold+1.0_WP*time%dt*dQdt(:,:,:,:,3)
         ! Apply IBM
         call apply_ibm()

         ! Fourth RK step ===================================================================================
         ! Particle increment
         call fs%get_div_stress(stressx,stressy,stressz,stressI)
         call lp%substep_rk4(stage=4,dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=fs%rho,visc=fs%visc,T=fs%T,C=fs%C,&
              stress_x=stressx,stress_y=stressy,stress_z=stressz,heat_flux=stressI,srcU=srcUlp,srcV=srcVlp,srcW=srcWlp,srcI=srcIlp)
         ! Get non-SL RHS and increment
         call fs%rhs(VF=lp%VF,VFU=lp%VFU,VFV=lp%VFV,VFW=lp%VFW,dQdt=dQdt(:,:,:,:,4))
         ! LPT source
         dQdt(:,:,:,2,4)=dQdt(:,:,:,2,4)+srcIlp
         dQdt(:,:,:,3,4)=dQdt(:,:,:,3,4)+srcUlp
         dQdt(:,:,:,4,4)=dQdt(:,:,:,4,4)+srcVlp
         dQdt(:,:,:,5,4)=dQdt(:,:,:,5,4)+srcWlp
         ! Add body forcing to momentum and energy
         dQdt(:,:,:,2,4)=dQdt(:,:,:,2,4)+Cv*fs%Q(:,:,:,1)*(T0-meanT)/time%dt
         dQdt(:,:,:,3,4)=dQdt(:,:,:,3,4)+(rhoU0-meanRhoU)/time%dt
         ! Advance
         fs%Q=fs%Qold+time%dt/6.0_WP*(dQdt(:,:,:,:,1)+2.0_WP*dQdt(:,:,:,:,2)+2.0_WP*dQdt(:,:,:,:,3)+dQdt(:,:,:,:,4))
         ! Apply IBM
         call apply_ibm()

         ! Interpolate velocity
         call fs%interp_vel(Ui,Vi,Wi)

         ! Compute local Mach number
         Ma=sqrt(Ui**2+Vi**2+Wi**2)/fs%C

         ! Compute dilatation
         call get_div()

         !> Perform and output monitoring
         call get_bodyforce()
         call fs%get_info()
         call lp%get_max()
         call mfile%write()
         call cflfile%write()
         call consfile%write()
         call lptfile%write()

         ! Output to ensight
         if (ens_evt%occurs()) then
            update_pmesh: block
              integer :: i
              call lp%update_partmesh(pmesh)
              do i=1,lp%np_
                 pmesh%var(1,i)=lp%p(i)%d
                 pmesh%var(2,i)=lp%p(i)%T
                 pmesh%vec(:,1,i)=lp%p(i)%vel
              end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if

         ! Finally, see if it's time to save restart files
         if (save_evt%occurs()) then
            save_restart: block
              character(len=str_medium) :: timestamp
              ! Prefix for files
              write(timestamp,'(es12.5)') time%t
              ! Populate df and write it
              call df%pushval(name='t' ,val=time%t        )
              call df%pushval(name='dt',val=time%dt       )
              call df%pushval(name='meanRhoU',val=rhoU0   )
              call df%pushval(name='meanT',val=T0         )
              call df%pushvar(name='Q1' ,var=fs%Q(:,:,:,1))
              call df%pushvar(name='Q2' ,var=fs%Q(:,:,:,2))
              call df%pushvar(name='Q3' ,var=fs%Q(:,:,:,3))
              call df%pushvar(name='Q4' ,var=fs%Q(:,:,:,4))
              call df%pushvar(name='Q5' ,var=fs%Q(:,:,:,5))
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
      deallocate(dQdt,Ui,Vi,Wi,Ma,beta,visc,visc_t,div,srcUlp,srcVlp,srcWlp,srcIlp,stressx,stressy,stressz,stressI)
      
   end subroutine simulation_final
   
   
end module simulation
