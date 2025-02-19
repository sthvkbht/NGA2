!> Various definitions and tools for running an NGA2 simulation
module simulation
  use string,            only: str_medium
  use precision,         only: WP
  use geometry,          only: cfg,D,get_VF
  use lowmach_class,     only: lowmach
  use fourier3d_class,   only: fourier3d
  use hypre_str_class,   only: hypre_str
  use lpt_class,         only: lpt
  use timetracker_class, only: timetracker
  use ensight_class,     only: ensight
  use partmesh_class,    only: partmesh
  use event_class,       only: event
  use datafile_class,    only: datafile
  use monitor_class,     only: monitor
  implicit none
  private

  !> Get an LPT solver, a lowmach solver, and corresponding time tracker
  type(lowmach),     public :: fs
  type(fourier3d),   public :: ps
  type(hypre_str),   public :: vs
  type(lpt),         public :: lp
  type(timetracker), public :: time

  !> Provide a datafile and an event tracker for saving restarts
  type(event)    :: save_evt
  type(datafile) :: df
  logical :: restarted

  !> Ensight postprocessing
  type(ensight)  :: ens_out
  type(partmesh) :: pmesh
  type(event)    :: ens_evt

  !> Event for post-processing
  type(event) :: ppevt

  !> Simulation monitor file
  type(monitor) :: mfile,cflfile,lptfile,tfile

  public :: simulation_init,simulation_run,simulation_final

  !> Work arrays and fluid properties
  real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
  real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,rho0,dRHOdt
  real(WP), dimension(:,:,:), allocatable :: srcUlp,srcVlp,srcWlp
  real(WP), dimension(:,:,:), allocatable :: G
  real(WP), dimension(:,:,:,:), allocatable :: Gnorm
  real(WP) :: visc,rho,mfr,mfr_target,bforce

  !> Wallclock time for monitoring
  type :: timer
     real(WP) :: time_in
     real(WP) :: time
     real(WP) :: percent
  end type timer
  type(timer) :: wt_total,wt_vel,wt_pres,wt_lpt,wt_rest

contains


    !> Specialized subroutine that outputs mean statistics
   subroutine postproc()
      use string,    only: str_medium
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      implicit none
      integer :: iunit,ierr,i,j,k
      real(WP), dimension(:), allocatable :: Uavg,Uavg_,vol,vol_
      character(len=str_medium) :: filename,timestamp
!!$      ! Allocate radial line storage
!!$      allocate(Uavg (fs%cfg%jmin:fs%cfg%jmax)); Uavg =0.0_WP
!!$      allocate(Uavg_(fs%cfg%jmin:fs%cfg%jmax)); Uavg_=0.0_WP
!!$      allocate(vol_ (fs%cfg%jmin:fs%cfg%jmax)); vol_ =0.0_WP
!!$      allocate(vol  (fs%cfg%jmin:fs%cfg%jmax)); vol  =0.0_WP
!!$      ! Integrate all data over x and z
!!$      do k=fs%cfg%kmin_,fs%cfg%kmax_
!!$         do j=fs%cfg%jmin_,fs%cfg%jmax_
!!$            do i=fs%cfg%imin_,fs%cfg%imax_
!!$               vol_(j) = vol_(j)+fs%cfg%vol(i,j,k)
!!$               Uavg_(j)=Uavg_(j)+fs%cfg%vol(i,j,k)*fs%U(i,j,k)
!!$            end do
!!$         end do
!!$      end do
!!$      ! All-reduce the data
!!$      call MPI_ALLREDUCE( vol_, vol,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
!!$      call MPI_ALLREDUCE(Uavg_,Uavg,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
!!$      do j=fs%cfg%jmin,fs%cfg%jmax
!!$         if (vol(j).gt.0.0_WP) then
!!$            Uavg(j)=Uavg(j)/vol(j)
!!$         else
!!$            Uavg(j)=0.0_WP
!!$         end if
!!$      end do
!!$      ! If root, print it out
!!$      if (fs%cfg%amRoot) then
!!$         filename='Uavg_'
!!$         write(timestamp,'(es12.5)') time%t
!!$         open(newunit=iunit,file=trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
!!$         write(iunit,'(a12,3x,a12)') 'Height','Uavg'
!!$         do j=fs%cfg%jmin,fs%cfg%jmax
!!$            write(iunit,'(es12.5,3x,es12.5)') fs%cfg%ym(j),Uavg(j)
!!$         end do
!!$         close(iunit)
!!$      end if
!!$      ! Deallocate work arrays
!!$      deallocate(Uavg,Uavg_,vol,vol_)
    end subroutine postproc
    

  !> Compute massflow rate
   function get_bodyforce_mfr(srcU) result(mfr)
      use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      real(WP), dimension(fs%cfg%imino_:,fs%cfg%jmino_:,fs%cfg%kmino_:), intent(in), optional :: srcU !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,ierr
      real(WP) :: vol,myRhoU,myUvol,Uvol,mfr
      myRhoU=0.0_WP; myUvol=0.0_WP
      if (present(srcU)) then
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  vol=fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*get_VF(i,j,k,'U')
                  myUvol=myUvol+vol
                  myRhoU=myRhoU+vol*(fs%rhoU(i,j,k)+srcU(i,j,k))
               end do
            end do
         end do
      else
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  vol=fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*get_VF(i,j,k,'U')
                  myUvol=myUvol+vol
                  myRhoU=myRhoU+vol*fs%rhoU(i,j,k)
               end do
            end do
         end do
      end if
      call MPI_ALLREDUCE(myUvol,Uvol,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(myRhoU,mfr ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); mfr=mfr/Uvol
   end function get_bodyforce_mfr


  !> Initialization of problem solver
  subroutine simulation_init
    use param, only: param_read
    implicit none


    ! Initialize time tracker with 1 subiterations
    initialize_timetracker: block
      time=timetracker(amRoot=cfg%amRoot)
      call param_read('Max timestep size',time%dtmax)
      call param_read('Max time',time%tmax)
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
         ! If we are not restarting, we will still need a datafile for saving restart files
         df=datafile(pg=cfg,filename=trim(cfg%name),nval=4,nvar=4)
         df%valname(1)='t'
         df%valname(2)='dt'
         df%valname(3)='mfr'
         df%valname(4)='bforce'
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

    
    ! Initialize wallclock timers
    initialize_timers: block
      wt_total%time=0.0_WP; wt_total%percent=0.0_WP
      wt_vel%time=0.0_WP;   wt_vel%percent=0.0_WP
      wt_pres%time=0.0_WP;  wt_pres%percent=0.0_WP
      wt_lpt%time=0.0_WP;   wt_lpt%percent=0.0_WP
      wt_rest%time=0.0_WP;  wt_rest%percent=0.0_WP
    end block initialize_timers


    ! Create a low Mach flow solver with bconds
    create_flow_solver: block
      use hypre_str_class, only: pcg_pfmg
      ! Create flow solver
      fs=lowmach(cfg=cfg,name='Variable density low Mach NS')

      ! Assign constant density
      call param_read('Density',rho); fs%rho=rho
      ! Assign constant viscosity
      call param_read('Dynamic viscosity',visc); fs%visc=visc
      ! Assign acceleration of gravity
      call param_read('Gravity',fs%gravity)
      ! Configure pressure solver
      ps=fourier3d(cfg=cfg,name='Pressure',nst=7)
      ! Configure implicit velocity solver
      vs=hypre_str(cfg=cfg,name='Velocity',method=pcg_pfmg,nst=7)
      call param_read('Implicit iteration',vs%maxit)
      call param_read('Implicit tolerance',vs%rcvg)
      ! Setup the solver
      call fs%setup(pressure_solver=ps,implicit_solver=vs)
    end block create_flow_solver


    ! Allocate work arrays
    allocate_work_arrays: block
      allocate(dRHOdt  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resU    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resV    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resW    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(srcUlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(srcVlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(srcWlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Ui      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Vi      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Wi      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(rho0    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(G       (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Gnorm   (1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
    end block allocate_work_arrays


    ! Initialize our LPT solver
    initialize_lpt: block
      use random, only: random_lognormal,random_uniform
      use mathtools, only: Pi,twoPi
      use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE,MPI_INTEGER
      use parallel, only: MPI_REAL_WP
      real(WP) :: VFavg,Volp,Vol_,sumVolp,Tp,meand,dist,r,buf
      real(WP) :: dmean,dsd,dmin,dmax,dshift
      real(WP), dimension(:), allocatable :: dp
      integer :: i,j,k,ii,jj,kk,nn,ip,jp,kp,np,offset,ierr
      integer, dimension(:,:,:), allocatable :: npic      !< Number of particle in cell
      integer, dimension(:,:,:,:), allocatable :: ipic    !< Index of particle in cell
      character(len=str_medium) :: timestamp
      logical :: overlap
      ! Create solver
      lp=lpt(cfg=cfg,name='LPT')
      ! Get pipe diameter from input
      ! Get mean volume fraction from input
      call param_read('Particle volume fraction',VFavg)
      ! Get drag model from input
      call param_read('Drag model',lp%drag_model,default='Tenneti')
      ! Get particle density from input
      call param_read('Particle density',lp%rho)
      ! Get particle diameter from input
      call param_read('Particle mean diameter',dmean)
      call param_read('Particle standard deviation',dsd,default=0.0_WP)
      call param_read('Particle min diameter',dmin,default=tiny(1.0_WP))
      call param_read('Particle max diameter',dmax,default=huge(1.0_WP))
      call param_read('Particle diameter shift',dshift,default=0.0_WP)
      if (dsd.le.epsilon(1.0_WP)) then
         dmin=dmean
         dmax=dmean
      end if
      ! Get particle temperature from input
      call param_read('Particle temperature',Tp,default=298.15_WP)
      ! Set collision timescale
      call param_read('Collision timescale',lp%tau_col,default=15.0_WP*time%dt)
      ! Set coefficient of restitution
      call param_read('Coefficient of restitution',lp%e_n)
      call param_read('Wall restitution',lp%e_w,default=lp%e_n)
      call param_read('Friction coefficient',lp%mu_f,default=0.0_WP)
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
         ! Get volume of domain belonging to this proc
         Vol_=0.0_WP
         do k=lp%cfg%kmin_,lp%cfg%kmax_
            do j=lp%cfg%jmin_,lp%cfg%jmax_
               do i=lp%cfg%imin_,fs%cfg%imax_
                  Vol_=Vol_+lp%cfg%dx(i)*lp%cfg%dy(j)*lp%cfg%dz(k)
               end do
            end do
         end do
         ! Get particle diameters
         np=5*VFavg*Vol_/(pi*dmean**3/6.0_WP)
         allocate(dp(np))
         sumVolp=0.0_WP; np=0
         do while(sumVolp.lt.VFavg*Vol_)
            np=np+1
            dp(np)=random_lognormal(m=dmean-dshift,sd=dsd)+dshift
            do while (dp(np).gt.dmax+epsilon(1.0_WP).or.dp(np).lt.dmin-epsilon(1.0_WP))
               dp(np)=random_lognormal(m=dmean-dshift,sd=dsd)+dshift
            end do
            sumVolp=sumVolp+Pi/6.0_WP*dp(np)**3
         end do
         call lp%resize(np)
         ! Allocate particle in cell arrays
         allocate(npic(     lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); npic=0
         allocate(ipic(1:40,lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); ipic=0
         ! Distribute particles
         sumVolp=0.0_WP; meand=0.0_WP
         do i=1,np
            ! Set the diameter
            lp%p(i)%d=dp(i)
            ! Give position (avoid overlap)
            overlap=.true.
            do while (overlap)
               lp%p(i)%pos=[random_uniform(lp%cfg%x(lp%cfg%imin_),lp%cfg%x(lp%cfg%imax_+1)-dmax),&
                    &       random_uniform(lp%cfg%y(lp%cfg%jmin_),lp%cfg%y(lp%cfg%jmax_+1)-dmax),&
                    &       random_uniform(lp%cfg%z(lp%cfg%kmin_),lp%cfg%z(lp%cfg%kmax_+1)-dmax)]
               lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
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
            ! Check if particle is within the pipe
            r=sqrt(lp%p(i)%pos(2)**2+lp%p(i)%pos(3)**2)
            if (r.le.0.5_WP*(D-lp%p(i)%d)) then
               ! Activate the particle
               lp%p(i)%flag=0
               ip=lp%p(i)%ind(1); jp=lp%p(i)%ind(2); kp=lp%p(i)%ind(3)
               npic(ip,jp,kp)=npic(ip,jp,kp)+1
               ipic(npic(ip,jp,kp),ip,jp,kp)=i
               ! Set the temperature
               lp%p(i)%T=Tp
               ! Give zero velocity
               lp%p(i)%vel=0.0_WP
               ! Give zero collision force
               lp%p(i)%Acol=0.0_WP
               lp%p(i)%Tcol=0.0_WP
               ! Give zero dt
               lp%p(i)%dt=0.0_WP
               ! Sum up volume
               sumVolp=sumVolp+Pi/6.0_WP*lp%p(i)%d**3
               meand=meand+lp%p(i)%d
            else
               lp%p(i)%flag=1
            end if
         end do
         deallocate(dp,npic,ipic)
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
         call MPI_ALLREDUCE(meand,buf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meand=buf/real(lp%np,WP)
         call MPI_ALLREDUCE(sumVolp,VFavg,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); VFavg=VFavg/(Pi*D**2/4.0_WP*lp%cfg%xL)
         if (lp%cfg%amRoot) then
            print*,"===== Particle Setup Description ====="
            print*,'Number of particles', lp%np
            print*,'Mean diameter', meand
            print*,'Mean volume fraction',VFavg
         end if
      end if
      ! Get initial particle volume fraction
      call lp%update_VF()
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


    ! Initialize our velocity field
    initialize_velocity: block
      real(WP) :: Ubulk
      if (restarted) then
         call df%pullvar(name='U',var=fs%U)
         call df%pullvar(name='V',var=fs%V)
         call df%pullvar(name='W',var=fs%W)
         call df%pullvar(name='P',var=fs%P)
      else
         ! Initial velocity
         call param_read('Bulk velocity',Ubulk)
         fs%U=Ubulk; fs%V=0.0_WP; fs%W=0.0_WP; fs%P=0.0_WP
      end if
      ! Set density from particle volume fraction and store initial density
      fs%rho=rho*(1.0_WP-lp%VF)
      rho0=rho
      ! Form momentum
      call fs%rho_multiply
      fs%rhoUold=fs%rhoU
      ! Apply all other boundary conditions
      call fs%interp_vel(Ui,Vi,Wi)
      call fs%get_div(drhodt=dRHOdt)
      ! Compute MFR through all boundary conditions
      call fs%get_mfr()
      ! Set initial MFR and body force
      if (restarted) then
         call df%pullval(name='mfr',val=mfr)
         call df%pullval(name='bforce',val=bforce)
      else
         mfr=get_bodyforce_mfr()
         bforce=0.0_WP
      end if
      mfr_target=mfr
    end block initialize_velocity


    ! Initialize levelset and its normal
    initialize_G: block
      integer :: i,j,k
      real(WP) :: buf
      real(WP), dimension(3) :: n12
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               G(i,j,k)=0.5_WP*D-sqrt(fs%cfg%ym(j)**2+fs%cfg%zm(k)**2)
               n12(1)=0.0_WP
               n12(2)=-fs%cfg%ym(j)
               n12(3)=-fs%cfg%zm(k)
               buf = sqrt(sum(n12*n12))+epsilon(1.0_WP)
               Gnorm(:,i,j,k)=n12/buf
            end do
         end do
      end do
      call fs%cfg%sync(G)
      call fs%cfg%sync(Gnorm)
    end block initialize_G

    ! Add Ensight output
    create_ensight: block
      ! Create Ensight output from cfg
      ens_out=ensight(cfg=cfg,name='riser')
      ! Create event for Ensight output
      ens_evt=event(time=time,name='Ensight output')
      call param_read('Ensight output period',ens_evt%tper)
      ! Add variables to output
      call ens_out%add_particle('particles',pmesh)
      call ens_out%add_vector('velocity',Ui,Vi,Wi)
      call ens_out%add_scalar('levelset',G)
      call ens_out%add_scalar('ibm_vf',fs%cfg%VF)
      call ens_out%add_scalar('pressure',fs%P)
      call ens_out%add_scalar('epsp',lp%VF)
      ! Output to ensight
      if (ens_evt%occurs()) call ens_out%write_data(time%t)
    end block create_ensight

    ! Create monitor filea
    create_monitor: block
      real(WP) :: cfl
      ! Prepare some info about fields
      call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
      call fs%get_cfl(time%dt,cfl); time%cfl=max(time%cfl,cfl)
      call fs%get_max()
      call lp%get_max()
      ! Create simulation monitor
      mfile=monitor(fs%cfg%amRoot,'simulation')
      call mfile%add_column(time%n,'Timestep number')
      call mfile%add_column(time%t,'Time')
      call mfile%add_column(time%dt,'Timestep size')
      call mfile%add_column(time%cfl,'Maximum CFL')
      call mfile%add_column(mfr,'MFR')
      call mfile%add_column(bforce,'Body force')
      call mfile%add_column(fs%Umax,'Umax')
      call mfile%add_column(fs%Vmax,'Vmax')
      call mfile%add_column(fs%Wmax,'Wmax')
      call mfile%add_column(fs%Pmax,'Pmax')
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
      call cflfile%add_column(lp%CFL_col,'Collision CFL')
      call cflfile%write()
      ! Create LPT monitor
      lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
      call lptfile%add_column(time%n,'Timestep number')
      call lptfile%add_column(time%t,'Time')
      call lptfile%add_column(lp%VFmean,'VFp mean')
      call lptfile%add_column(lp%VFmax,'VFp max')
      call lptfile%add_column(lp%np,'Particle number')
      call lptfile%add_column(lp%ncol,'Collision number')
      call lptfile%add_column(lp%Umin,'Particle Umin')
      call lptfile%add_column(lp%Umax,'Particle Umax')
      call lptfile%add_column(lp%Vmin,'Particle Vmin')
      call lptfile%add_column(lp%Vmax,'Particle Vmax')
      call lptfile%add_column(lp%Wmin,'Particle Wmin')
      call lptfile%add_column(lp%Wmax,'Particle Wmax')
      call lptfile%add_column(lp%dmin,'Particle dmin')
      call lptfile%add_column(lp%dmax,'Particle dmax')
      call lptfile%write()
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
      call tfile%add_column(wt_rest%time,'Rest [s]')
      call tfile%add_column(wt_rest%percent,'Rest [%]')
      call tfile%write()
    end block create_monitor


    ! Create a specialized post-processing file
    create_postproc: block
      ! Create event for data postprocessing
      ppevt=event(time=time,name='Postproc output')
      call param_read('Postproc output period',ppevt%tper)
      ! Perform the output
      if (ppevt%occurs()) call postproc()
    end block create_postproc

  end subroutine simulation_init


  !> Perform an NGA2 simulation
  subroutine simulation_run
    use parallel, only: parallel_time
    implicit none
    real(WP) :: cfl

    ! Perform time integration
    do while (.not.time%done())

       ! Initial wallclock time
       wt_total%time_in=parallel_time()

       ! Increment time
       call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
       call fs%get_cfl(time%dt,cfl); time%cfl=max(time%cfl,cfl)
       call time%adjust_dt()
       call time%increment()

       ! Remember old density, velocity, and momentum
       fs%rhoold=fs%rho
       fs%Uold=fs%U; fs%rhoUold=fs%rhoU
       fs%Vold=fs%V; fs%rhoVold=fs%rhoV
       fs%Wold=fs%W; fs%rhoWold=fs%rhoW

       ! Get fluid stress (include mean body force from mfr forcing)
       wt_lpt%time_in=parallel_time()
       call fs%get_div_stress(resU,resV,resW)
       resU=resU+bforce

       ! Collide and advance particles
       call lp%collide(dt=time%dtmid,Gib=G,Nxib=Gnorm(1,:,:,:),Nyib=Gnorm(2,:,:,:),Nzib=Gnorm(3,:,:,:))
       call lp%advance(dt=time%dtmid,U=fs%U,V=fs%V,W=fs%W,rho=rho0,visc=fs%visc,stress_x=resU,stress_y=resV,stress_z=resW,&
            srcU=srcUlp,srcV=srcVlp,srcW=srcWlp)

       ! Update density based on particle volume fraction
       fs%rho=rho*(1.0_WP-lp%VF)
       dRHOdt=(fs%RHO-fs%RHOold)/time%dtmid
       wt_lpt%time=wt_lpt%time+parallel_time()-wt_lpt%time_in

       ! Perform sub-iterations
       do while (time%it.le.time%itmax)

          wt_vel%time_in=parallel_time()

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
                     resU(i,j,k)=resU(i,j,k)+sum(fs%itpr_x(:,i,j,k)*srcUlp(i-1:i,j,k))
                     resV(i,j,k)=resV(i,j,k)+sum(fs%itpr_y(:,i,j,k)*srcVlp(i,j-1:j,k))
                     resW(i,j,k)=resW(i,j,k)+sum(fs%itpr_z(:,i,j,k)*srcWlp(i,j,k-1:k))
                  end do
               end do
            end do
            call fs%cfg%sync(resU)
            call fs%cfg%sync(resV)
            call fs%cfg%sync(resW)
          end block add_lpt_src

          ! Apply direct forcing to enforce BC at the pipe walls
          ibm_correction: block
            integer :: i,j,k
            real(WP) :: VFx,VFy,VFz,RHOx,RHOy,RHOz
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do j=fs%cfg%jmin_,fs%cfg%jmax_
                  do i=fs%cfg%imin_,fs%cfg%imax_
                     VFx=get_VF(i,j,k,'U')
                     VFy=get_VF(i,j,k,'V')
                     VFz=get_VF(i,j,k,'W')
                     RHOx=sum(fs%itpr_x(:,i,j,k)*fs%rho(i-1:i,j,k))
                     RHOy=sum(fs%itpr_y(:,i,j,k)*fs%rho(i,j-1:j,k))
                     RHOz=sum(fs%itpr_z(:,i,j,k)*fs%rho(i,j,k-1:k))
                     resU(i,j,k)=resU(i,j,k)-(1.0_WP-VFx)*RHOx*fs%U(i,j,k)
                     resV(i,j,k)=resV(i,j,k)-(1.0_WP-VFy)*RHOy*fs%V(i,j,k)
                     resW(i,j,k)=resW(i,j,k)-(1.0_WP-VFz)*RHOz*fs%W(i,j,k)
                  end do
               end do
            end do
          end block ibm_correction

          ! Add body forcing
          bodyforcing: block
            mfr=get_bodyforce_mfr(resU)
            bforce=(mfr_target-mfr)/time%dtmid
            resU=resU+time%dtmid*bforce
          end block bodyforcing

          ! Form implicit residuals
          call fs%solve_implicit(time%dtmid,resU,resV,resW)

          ! Apply these residuals
          fs%U=2.0_WP*fs%U-fs%Uold+resU
          fs%V=2.0_WP*fs%V-fs%Vold+resV
          fs%W=2.0_WP*fs%W-fs%Wold+resW

!!$          ! Apply direct forcing to enforce BC at the pipe walls
!!$          ibm_correction: block
!!$            integer :: i,j,k
!!$            real(WP) :: VFx,VFy,VFz
!!$            do k=fs%cfg%kmin_,fs%cfg%kmax_
!!$               do j=fs%cfg%jmin_,fs%cfg%jmax_
!!$                  do i=fs%cfg%imin_,fs%cfg%imax_
!!$                     VFx=get_VF(i,j,k,'U')
!!$                     VFy=get_VF(i,j,k,'V')
!!$                     VFz=get_VF(i,j,k,'W')
!!$                     fs%U(i,j,k)=fs%U(i,j,k)*VFx
!!$                     fs%V(i,j,k)=fs%V(i,j,k)*VFy
!!$                     fs%W(i,j,k)=fs%W(i,j,k)*VFz
!!$                  end do
!!$               end do
!!$            end do
!!$            call fs%cfg%sync(fs%U)
!!$            call fs%cfg%sync(fs%V)
!!$            call fs%cfg%sync(fs%W)
!!$          end block ibm_correction

          ! Apply other boundary conditions and update momentum
          call fs%apply_bcond(time%tmid,time%dtmid)
          call fs%rho_multiply()
          call fs%apply_bcond(time%tmid,time%dtmid)

          wt_vel%time=wt_vel%time+parallel_time()-wt_vel%time_in

          ! Solve Poisson equation
          wt_pres%time_in=parallel_time()
          call fs%correct_mfr(drhodt=dRHOdt)
          call fs%get_div(drhodt=dRHOdt)
          fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dtmid
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
          wt_pres%time=wt_pres%time+parallel_time()-wt_pres%time_in

          ! Increment sub-iteration counter
          time%it=time%it+1

       end do

       ! Recompute interpolated velocity and divergence
       wt_vel%time_in=parallel_time()
       call fs%interp_vel(Ui,Vi,Wi)
       call fs%get_div(drhodt=dRHOdt)
       wt_vel%time=wt_vel%time+parallel_time()-wt_vel%time_in

       ! Recompute massflow rate
       mfr=get_bodyforce_mfr()

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

       ! Specialized post-processing
       if (ppevt%occurs()) call postproc()

       ! Perform and output monitoring
       call fs%get_max()
       call lp%get_max()
       call mfile%write()
       call cflfile%write()
       call lptfile%write()

       ! Monitor timing
       wt_total%time=parallel_time()-wt_total%time_in
       wt_vel%percent=wt_vel%time/wt_total%time*100.0_WP
       wt_pres%percent=wt_pres%time/wt_total%time*100.0_WP
       wt_lpt%percent=wt_lpt%time/wt_total%time*100.0_WP
       wt_rest%time=wt_total%time-wt_vel%time-wt_pres%time-wt_lpt%time
       wt_rest%percent=wt_rest%time/wt_total%time*100.0_WP
       call tfile%write()
       wt_total%time=0.0_WP; wt_total%percent=0.0_WP
       wt_vel%time=0.0_WP;   wt_vel%percent=0.0_WP
       wt_pres%time=0.0_WP;  wt_pres%percent=0.0_WP
       wt_lpt%time=0.0_WP;   wt_lpt%percent=0.0_WP
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
            call df%pushval(name='mfr',val=mfr_target)
            call df%pushval(name='bforce',val=bforce )
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
    ! timetracker

    ! Deallocate work arrays
    deallocate(resU,resV,resW,srcUlp,srcVlp,srcWlp,Ui,Vi,Wi,dRHOdt,G,Gnorm)

  end subroutine simulation_final

end module simulation
