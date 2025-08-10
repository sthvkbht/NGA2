!> Various definitions and tools for running an NGA2 simulation
module simulation
  use string,            only: str_medium
  use precision,         only: WP
  use geometry,          only: cfg,R0
  use multisphere_class, only: multisphere,get_positions,max_part
  use timetracker_class, only: timetracker
  use ensight_class,     only: ensight
  use partmesh_class,    only: partmesh
  use event_class,       only: event
  use monitor_class,     only: monitor
  implicit none
  private

  !> Get a multisphere solver and corresponding time tracker, plus a couple of linear solvers
  type(multisphere),  public :: ms
  type(timetracker),  public :: time

  !> Ensight postprocessing
  type(ensight)  :: ens_out
  type(partmesh) :: pmesh
  type(event)    :: ens_evt

  !> Simulation monitor file
  type(monitor) :: mfile,pfile

  public :: simulation_init,simulation_run,simulation_final

  !> Work arrays and fluid properties
  real(WP), dimension(:,:,:), allocatable :: dpdx,dpdy,dpdz
  real(WP), dimension(:,:,:), allocatable :: U,V,W,P
  real(WP) :: dp,visc,rho,inlet_velocity,Wmax,Pmax

contains


  !> Subroutine that computes fluid velocity, pressure, and pressure gradient
  subroutine get_pressure()
    use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_IN_PLACE
    use parallel, only: MPI_REAL_WP
    implicit none
    integer :: i,j,k,ierr
    real(WP) :: eps,A,B,dpdz0,W0
    real(WP), parameter :: eps_ref=0.999_WP
    ! Use constant pressure gradient for desired inlet velocity
    A=150.0_WP*visc/dp**2*(1.0_WP-eps_ref)**2/eps_ref**3
    B=1.75_WP*rho/dp*(1.0_WP-eps_ref)/eps_ref**3
    dpdz0=-(A*inlet_velocity+B*inlet_velocity**2)
    dpdz=dpdz0
    ! Compute pressure based on Egun's equaton
    do k=cfg%kmin_,cfg%kmax_
       do j=cfg%jmin_,cfg%jmax_
          do i=cfg%imin_,cfg%imax_
             ! Voidage
             eps=1.0_WP-ms%VF(i,j,k)
             if (eps.lt.eps_ref) then
                ! Get velocity from Ergun equation
                A=150.0_WP*visc/dp**2*(1.0_WP-eps)**2/(eps**3+epsilon(1.0_WP))
                B=1.75_WP*rho/dp*(1.0_WP-eps)/(eps**3+epsilon(1.0_WP))
                W(i,j,k)=(-A+sqrt(A**2-4.0_WP*B*dpdz0))/(2.0_WP*B)
             else
                W(i,j,k)=inlet_velocity
             end if
             ! Get pressure
             P(i,j,k)=0.5_WP*rho*(inlet_velocity**2-W(i,j,k)**2)
!!$               ! Velocity
!!$               W(i,j,k)=inlet_velocity/(eps+epsilon(1.0_WP))
!!$               ! Pressure
!!$               P(i,j,k)=150.0_WP*visc*ms%Hp/dp**2*(1.0_WP-eps)**2/(eps**3+epsilon(1.0_WP))*inlet_velocity+&
!!$                    1.75_WP*ms%Hp*rho/dp*(1.0_WP-eps)/(eps**3+epsilon(1.0_WP))*inlet_velocity**2
!!$               ! Pressure drop
!!$               dPdz(i,j,k)=P(i,j,k)/ms%HP
          end do
       end do
    end do
    ! Ensure conservation
    call cfg%integrate(A=W,integral=W0); W0=W0/cfg%vol_total
    W=W*inlet_velocity/W0
    ! Synchronize it
    call cfg%sync(W)
    call cfg%sync(P)
    ! Compute the gradient
    do k=cfg%kmin_,cfg%kmax_
       do j=cfg%jmin_,cfg%jmax_
          do i=cfg%imin_,cfg%imax_
             dPdx(i,j,k)=sum(ms%grd_x(:,i,j,k)*P(i-1:i,j,k))
             dPdy(i,j,k)=sum(ms%grd_y(:,i,j,k)*P(i,j-1:j,k))
          end do
       end do
    end do
    ! Synchronize it
    call cfg%sync(dPdx)
    call cfg%sync(dPdy)
    ! Get max values
    Wmax=maxval(W); Pmax=maxval(P)
    call MPI_ALLREDUCE(MPI_IN_PLACE,Wmax,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,Pmax,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
  end subroutine get_pressure


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


    ! Allocate work arrays
    allocate_work_arrays: block
      allocate(dpdx(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(dpdy(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(dpdz(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(P   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(U   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(V   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(W   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
    end block allocate_work_arrays


    ! Initialize our multisphere solver
    initialize_multisphere: block
      use random, only: random_lognormal,random_uniform
      use mathtools, only: Pi,twoPi
      use messager, only: die
      real(WP) :: L,dx,dy,Ub
      real(WP), dimension(:,:), allocatable :: pos
      integer :: i,j,nb,npb,nx,ny,ix,iy
      character(len=str_medium) :: shape
      ! Create solver
      ms=multisphere(cfg=cfg,name='multi-sphere')
      ! Get shape of rigid body
      call param_read('Body shape',shape)
      ! Get number of rigid bodies
      call param_read('Number of bodies',nb)
      ! Get number of particles per body
      call param_read('Particles per body',npb)
      if (npb.gt.max_part) call die('Number of particles per body is too large')
      ! Side length of rigid body
      call param_read('Body side length',L)
      ! Set initial body velocity
      call param_read('Initial body velocity',Ub,default=0.0_WP)
      ! Get particle density from input
      call param_read('Particle density',ms%rho)
      ! Get particle diameter from input
      call param_read('Particle diameter',dp)
      ! Get particle height from input
      call param_read('Particle height',ms%Hp)
      ! Set collision timescale
      call param_read('Collision timescale',ms%tau_col,default=15.0_WP*time%dt)
      ! Set coefficient of restitution
      call param_read('Coefficient of restitution',ms%e_n)
      call param_read('Friction coefficient',ms%mu_f,default=0.0_WP)
      ! Set filter scale
      call param_read('Filter width',ms%filter_width,default=3.5_WP*cfg%min_meshsize)
      ! Initialize rigid bodies
      if (ms%cfg%amRoot) then
         nx=ceiling(sqrt(real(nb,WP))); ny=ceiling(real(nb,WP)/real(nx,WP))
         dx=cfg%xL/real(nx,WP); dy=cfg%yL/real(ny,WP)
         ! Generate the rigid bodies
         call ms%resize(nb)
         allocate(pos(2,npb))
         do j=1,nb
            ! Give ID
            ms%b(j)%id=j
            ! Give random velocity
            ms%b(j)%vel(1:2)=[random_uniform(-Ub,Ub),random_uniform(-Ub,Ub)]
            ms%b(j)%vel(3)=0.0_WP
            ! Give zero angular velocity
            ms%b(j)%angVel=0.0_WP
            ! Give zero force and torque
            ms%b(j)%force=0.0_WP
            ms%b(j)%torque=0.0_WP
            ! Set number of particles within this body
            ms%b(j)%np=npb
            ! Give center-of-mass position
            ix=mod(j-1,nx)+1
            iy=(j-1)/nx+1
            ms%b(j)%pos(1)=ms%cfg%x(ms%cfg%imin)+(real(ix,WP)-0.5_WP)*dx
            ms%b(j)%pos(2)=ms%cfg%y(ms%cfg%jmin)+(real(iy,WP)-0.5_WP)*dy
            ms%b(j)%pos(3)=ms%cfg%z(ms%cfg%kmin)+0.5_WP*ms%Hp
            ! Locate the rigid body on the mesh
            ms%b(j)%ind=ms%cfg%get_ijk_global(ms%b(j)%pos,[ms%cfg%imin,ms%cfg%jmin,ms%cfg%kmin])
            ! Give random angle
            ms%b(j)%theta=random_uniform(0.0_WP,twoPi)
            ! Get positions of each particle
            pos=get_positions(npb,ms%b(j)%pos(1),ms%b(j)%pos(2),ms%b(j)%theta,L,trim(adjustl(shape)))
            do i=1,npb
               ! Give diameter
               ms%b(j)%p_d(i)=dp
               ! Give zero collision force
               ms%b(j)%p_col(:,i)=0.0_WP
               ! Give position
               ms%b(j)%p_pos(1:2,i)=pos(:,i)
               ms%b(j)%p_pos(3,i)=ms%b(j)%pos(3)
               ! Locate the particle on the mesh
               ms%b(j)%p_ind(:,i)=ms%cfg%get_ijk_global(ms%b(j)%p_pos(:,i),[ms%cfg%imin,ms%cfg%jmin,ms%cfg%kmin])
            end do
            ! Activate the rigid body
            ms%b(j)%flag=0
         end do
         deallocate(pos)
      end if
      call ms%sync()
      call ms%update_rigid_body()
      call ms%update_VF()
      if (ms%cfg%amRoot) then
         print*,"===== Body Setup Description ====="
         print*,'Number of rigid bodies', ms%nb
         print*,'Number of particles', ms%np
      end if
    end block initialize_multisphere


    ! Initialize fluid
    initialize_velocity: block
      integer :: n,i,j,k
      ! Zero initial field
      P=0.0_WP; U=0.0_WP; V=0.0_WP; W=0.0_WP
      ! Set inflow velocity/momentum
      call param_read('Inlet velocity',inlet_velocity)
      call param_read('Density',rho)
      call param_read('Dynamic viscosity',visc)
      call get_pressure()
    end block initialize_velocity


    ! Create partmesh object for Lagrangian particle output
    create_pmesh: block
      integer :: i,j,k
      pmesh=partmesh(nvar=2,nvec=2,name='multi-sphere')
      pmesh%varname(1)='id'
      pmesh%varname(2)='diameter'
      pmesh%vecname(1)='velocity'
      pmesh%vecname(2)='Fcol'
      call ms%update_partmesh(pmesh)
      k=0
      do j=1,ms%nb_
         do i=1,ms%b(j)%np
            k=k+1
            pmesh%var(1,k)=real(ms%b(j)%id,WP)
            pmesh%var(2,k)=ms%b(j)%p_d(i)
            pmesh%vec(:,1,k)=ms%b(j)%p_vel(:,i)
            pmesh%vec(:,2,k)=ms%b(j)%p_col(:,i)
         end do
      end do
    end block create_pmesh


    ! Add Ensight output
    create_ensight: block
      ! Create Ensight output from cfg
      ens_out=ensight(cfg=cfg,name='fluidized_bed')
      ! Create event for Ensight output
      ens_evt=event(time=time,name='Ensight output')
      call param_read('Ensight output period',ens_evt%tper)
      ! Add variables to output
      call ens_out%add_particle('particles',pmesh)
      call ens_out%add_vector('velocity',U,V,W)
      call ens_out%add_scalar('epsp',ms%VF)
      call ens_out%add_scalar('pressure',P)
      ! Output to ensight
      if (ens_evt%occurs()) call ens_out%write_data(time%t)
    end block create_ensight
    

    ! Create monitor filea
    create_monitor: block
      ! Prepare some info about fields
      real(WP) :: cfl
      call ms%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
      call ms%get_max()
      ! Create simulation monitor
      mfile=monitor(ms%cfg%amRoot,'simulation')
      call mfile%add_column(time%n,'Timestep number')
      call mfile%add_column(time%t,'Time')
      call mfile%add_column(time%dt,'Timestep size')
      call mfile%add_column(time%cfl,'Maximum CFL')
      call mfile%add_column(Wmax,'Wmax')
      call mfile%add_column(Pmax,'Pmax')
      call mfile%write()
      ! Create particle monitor
      pfile=monitor(amroot=ms%cfg%amRoot,name='particle')
      call pfile%add_column(time%n,'Timestep number')
      call pfile%add_column(time%t,'Time')
      call pfile%add_column(ms%VFmean,'VFp mean')
      call pfile%add_column(ms%VFmax,'VFp max')
      call pfile%add_column(ms%nb,'N bodies')
      call pfile%add_column(ms%np,'N particles')
      call pfile%add_column(ms%ncol,'N collisions')
      call pfile%add_column(ms%Umin,'Particle Umin')
      call pfile%add_column(ms%Umax,'Particle Umax')
      call pfile%add_column(ms%Vmin,'Particle Vmin')
      call pfile%add_column(ms%Vmax,'Particle Vmax')
      call pfile%write()
    end block create_monitor

  end subroutine simulation_init


  !> Perform an NGA2 simulation
  subroutine simulation_run
    use parallel, only: parallel_time
    implicit none

    ! Perform time integration
    do while (.not.time%done())

       ! Increment time
       call ms%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
       call time%adjust_dt()
       call time%increment()

       ! Get fluid fields
       call get_pressure()

       ! Update particles
       call ms%collide(dt=time%dtmid)
       call ms%advance(dt=time%dtmid,rho=rho,visc=visc,U=U,V=V,W=W,dpdx=dpdx,dpdy=dpdy,dpdz=dpdz)

       ! Output to ensight
       if (ens_evt%occurs()) then
          update_pmesh: block
            integer :: i,j,k
            call ms%update_partmesh(pmesh)
            k=0
            do j=1,ms%nb_
               do i=1,ms%b(j)%np
                  k=k+1
                  pmesh%var(1,k)=real(ms%b(j)%id,WP)
                  pmesh%var(2,k)=ms%b(j)%p_d(i)
                  pmesh%vec(:,1,k)=ms%b(j)%p_vel(:,i)
                  pmesh%vec(:,2,k)=ms%b(j)%p_col(:,i)
               end do
            end do
          end block update_pmesh
          call ens_out%write_data(time%t)
       end if

       ! Perform and output monitoring
       call ms%get_max()
       call mfile%write()
       call pfile%write()

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
    deallocate(U,V,W,P,dpdx,dpdy,dpdz)

  end subroutine simulation_final

end module simulation
