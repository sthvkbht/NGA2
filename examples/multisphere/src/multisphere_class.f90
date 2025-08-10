!> Basic multi-sphere solver class:
!> Provides support for Lagrangian-transported objects
module multisphere_class
  use precision,      only: WP
  use string,         only: str_medium
  use config_class,   only: config
  use ddadi_class,    only: ddadi
  use mpi_f08,        only: MPI_Datatype,MPI_INTEGER8,MPI_INTEGER,MPI_DOUBLE_PRECISION
  implicit none
  private

  ! Expose type/constructor/methods
  public :: multisphere

  ! Function to generate rigid body shapes
  public :: get_positions

  !> Memory adaptation parameter
  real(WP), parameter :: coeff_up=1.3_WP      !< Particle array size increase factor
  real(WP), parameter :: coeff_dn=0.7_WP      !< Particle array size decrease factor

  !> Maximum number of particles per body
  integer, parameter, public :: max_part=40

  !> I/O chunk size to read at a time
  integer, parameter :: part_chunk_size=1000  !< Read 1000 particles at a time before redistributing

  
  !> Basic body object definition
  type :: body
     !> MPI_DOUBLE_PRECISION data
     real(WP) :: m                            !< Total mass
     real(WP) :: theta                        !< Angle of rotation
     real(WP) :: I                            !< Moment of inertia (2D)
     real(WP), dimension(3) :: pos            !< Center-of-mass position
     real(WP), dimension(3) :: vel            !< Linear velocity
     real(WP), dimension(3) :: angVel         !< Angular velocity (body frame)
     real(WP), dimension(3) :: force          !< Net force
     real(WP), dimension(3) :: torque         !< Net torque
     real(WP), dimension(max_part) :: p_d     !< Particle diameter
     real(WP), dimension(3,max_part) :: p_pos !< Particle position
     real(WP), dimension(3,max_part) :: p_vel !< Particle velocity
     real(WP), dimension(3,max_part) :: p_col !< Collision force
     real(WP), dimension(3,max_part) :: p_r0  !< Relative position
     !> MPI_INTEGER data
     integer :: id                            !< Body ID
     integer :: np                            !< Number of particles in this body (<= max_part)
     integer :: flag                          !< Control parameter (0=normal, 1=done->will be removed)
     integer, dimension(3) :: ind             !< Index of cell containing body center of mass
     integer, dimension(3,max_part) :: p_ind  !< Index of cell containing particle center
  end type body

  !> Number of blocks, block lengths, and block types in a body
  integer, parameter                             :: body_nblock = 2
  integer           , dimension(body_nblock)     :: body_lblock = [18+13*max_part,6+3*max_part]
  type(MPI_Datatype), dimension(body_nblock)     :: body_tblock = [MPI_DOUBLE_PRECISION, MPI_INTEGER]
  !> MPI_BODY derived datatype and size
  type(MPI_Datatype) :: MPI_BODY
  integer            :: MPI_BODY_SIZE
  

  !> Multi-sphere solver object definition
  type :: multisphere

     ! This is our underlying config
     class(config), pointer :: cfg                       !< This is the config the solver is build for

     type(ddadi) :: implicit                             !< Implicit solver for filtering

     ! This is the name of the solver
     character(len=str_medium) :: name='UNNAMED_MULTISPHERE' !< Solver name (default=UNNAMED_MULTISPHERE)

     ! Particle data
     integer :: np                                       !< Global number of particles
     integer :: np_                                      !< Local number of particles

     ! Body data
     integer :: nb                                       !< Global number of bodies
     integer :: nb_                                      !< Local number of bodies
     integer, dimension(:), allocatable :: nb_proc       !< Number of rigid bodies on each processor
     integer, dimension(:), allocatable :: np_proc       !< Number of particles on each processor
     type(body), dimension(:), allocatable :: b          !< Array of rigid multi-sphere bodies

     ! Overlap body (i.e., ghost) data
     integer :: ng_                                      !< Local number of ghosts
     type(body), dimension(:), allocatable :: g          !< Array of ghosts of type body

     ! CFL numbers
     real(WP) :: CFLp_x,CFLp_y,CFLp_z,CFL_col            !< CFL numbers

     ! Maximum bounding radius
     real(WP) :: max_R0

     ! Particle density
     real(WP) :: rho                                     !< Density of particle

     ! Particle height
     real(WP) :: Hp                                      !< Height of particle

     ! Gravitational acceleration
     real(WP), dimension(3) :: gravity=0.0_WP            !< Acceleration of gravity

     ! Solver parameters
     real(WP) :: nstep=1                                 !< Number of substeps (default=1)
     
     ! Collisional parameters
     real(WP) :: tau_col                                 !< Characteristic collision time scale
     real(WP) :: e_n=1.0_WP                              !< Normal restitution coefficient
     real(WP) :: e_w=1.0_WP                              !< Wall restitution coefficient
     real(WP) :: mu_f=0.0_WP                             !< Friction coefficient
     real(WP) :: clip_col=0.2_WP                         !< Maximum allowable overlap
     real(WP), dimension(:,:,:), allocatable :: xwall    !< Location in x of closest wall in x
     real(WP), dimension(:,:,:), allocatable :: ywall    !< Location in y of closest wall in y
     real(WP), dimension(:,:,:), allocatable :: zwall    !< Location in z of closest wall in z

     ! Monitoring info
     real(WP) :: VFmin,VFmax,VFmean,VFvar                !< Volume fraction info
     real(WP) :: Umin,Umax,Umean,Uvar                    !< U velocity info
     real(WP) :: Vmin,Vmax,Vmean,Vvar                    !< V velocity info
     real(WP) :: Wmin,Wmax,Wmean,Wvar                    !< W velocity info
     integer  :: ncol=0                                  !< Number of collisions

     ! Particle volume fraction
     real(WP), dimension(:,:,:), allocatable :: VF       !< Particle volume fraction, cell-centered

     ! Filtering operation
     logical :: implicit_filter                          !< Solve implicitly
     real(WP) :: filter_width                            !< Characteristic filter width
     real(WP), dimension(:,:,:,:), allocatable :: div_x,div_y,div_z !< Divergence operator
     real(WP), dimension(:,:,:,:), allocatable :: grd_x,grd_y,grd_z !< Gradient operator

     ! Pseudo turbulence
     real(WP), dimension(:,:,:,:), allocatable :: itpr_x,itpr_y,itpr_z !< Interpolation from cell face to cell center

   contains
     procedure :: update_partmesh                        !< Update a partmesh object using current particles
     procedure :: collide                                !< Evaluate interparticle collision force
     procedure :: advance                                !< Step forward the particle ODEs
     procedure :: get_force                              !< Compute force and torque acting on each rigid body
     procedure :: resize                                 !< Resize particle array to given size
     procedure :: resize_ghost                           !< Resize ghost array to given size
     procedure :: recycle                                !< Recycle particle array by removing flagged particles
     procedure :: sync                                   !< Synchronize particles across interprocessor boundaries
     procedure :: share                                  !< Share particles across interprocessor boundaries
     procedure :: read                                   !< Parallel read particles from file
     procedure :: write                                  !< Parallel write particles to file
     procedure :: get_max                                !< Extract various monitoring data
     procedure :: get_cfl                                !< Calculate maximum CFL
     procedure :: update_VF                              !< Compute particle volume fraction
     procedure :: update_rigid_body                      !< Update various properties based on rigid body mechanics
     procedure :: filter                                 !< Apply volume filtering to field
  end type multisphere


  !> Declare multisphere solver constructor
  interface multisphere
     procedure constructor
  end interface multisphere

contains


  !> Default constructor for multisphere solver
  function constructor(cfg,name) result(self)
    implicit none
    type(multisphere) :: self
    class(config), target, intent(in) :: cfg
    character(len=*), optional :: name
    integer :: i,j,k,l

    ! Set the name for the solver
    if (present(name)) self%name=trim(adjustl(name))

    ! Point to pgrid object
    self%cfg=>cfg

    ! Allocate variables
    allocate(self%nb_proc(1:self%cfg%nproc)); self%nb_proc=0
    allocate(self%np_proc(1:self%cfg%nproc)); self%np_proc=0
    self%nb_=0; self%nb=0
    call self%resize(0)
    self%np_=0; self%np=0

    ! Initialize MPI derived datatype for a particle
    call prepare_mpi_body()

    ! Allocate VF on cfg mesh
    allocate(self%VF(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VF=0.0_WP

    ! Allocate finite volume divergence operators
    allocate(self%div_x(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
    allocate(self%div_y(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
    allocate(self%div_z(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
    ! Create divergence operator to cell center [xm,ym,zm]
    do k=self%cfg%kmin_,self%cfg%kmax_
       do j=self%cfg%jmin_,self%cfg%jmax_
          do i=self%cfg%imin_,self%cfg%imax_
             self%div_x(:,i,j,k)=self%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< Divergence from [x ,ym,zm]
             self%div_y(:,i,j,k)=self%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< Divergence from [xm,y ,zm]
             self%div_z(:,i,j,k)=self%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< Divergence from [xm,ym,z ]
          end do
       end do
    end do

    ! Allocate finite difference velocity gradient operators
    allocate(self%grd_x(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< X-face-centered
    allocate(self%grd_y(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Y-face-centered
    allocate(self%grd_z(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Z-face-centered
    ! Create gradient coefficients to cell faces
    do k=self%cfg%kmin_,self%cfg%kmax_+1
       do j=self%cfg%jmin_,self%cfg%jmax_+1
          do i=self%cfg%imin_,self%cfg%imax_+1
             self%grd_x(:,i,j,k)=self%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< Gradient in x from [xm,ym,zm] to [x,ym,zm]
             self%grd_y(:,i,j,k)=self%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< Gradient in y from [xm,ym,zm] to [xm,y,zm]
             self%grd_z(:,i,j,k)=self%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< Gradient in z from [xm,ym,zm] to [xm,ym,z]
          end do
       end do
    end do

    ! Allocate finite difference density interpolation coefficients
    allocate(self%itpr_x(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< X-face-centered
    allocate(self%itpr_y(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Y-face-centered
    allocate(self%itpr_z(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Z-face-centered
    ! Create density interpolation coefficients to cell face
    do k=self%cfg%kmin_,self%cfg%kmax_+1
       do j=self%cfg%jmin_,self%cfg%jmax_+1
          do i=self%cfg%imin_,self%cfg%imax_+1
             self%itpr_x(:,i,j,k)=self%cfg%dxmi(i)*[self%cfg%xm(i)-self%cfg%x(i),self%cfg%x(i)-self%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
             self%itpr_y(:,i,j,k)=self%cfg%dymi(j)*[self%cfg%ym(j)-self%cfg%y(j),self%cfg%y(j)-self%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
             self%itpr_z(:,i,j,k)=self%cfg%dzmi(k)*[self%cfg%zm(k)-self%cfg%z(k),self%cfg%z(k)-self%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
          end do
       end do
    end do
    
    ! Loop over the domain and zero divergence in walls
    do k=self%cfg%kmin_,self%cfg%kmax_
       do j=self%cfg%jmin_,self%cfg%jmax_
          do i=self%cfg%imin_,self%cfg%imax_
             if (self%cfg%VF(i,j,k).eq.0.0_WP) then
                self%div_x(:,i,j,k)=0.0_WP
                self%div_y(:,i,j,k)=0.0_WP
                self%div_z(:,i,j,k)=0.0_WP
             end if
          end do
       end do
    end do

    ! Zero out gradient to wall faces
    do k=self%cfg%kmin_,self%cfg%kmax_+1
       do j=self%cfg%jmin_,self%cfg%jmax_+1
          do i=self%cfg%imin_,self%cfg%imax_+1
             if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i-1,j,k).eq.0.0_WP) self%grd_x(:,i,j,k)=0.0_WP
             if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i,j-1,k).eq.0.0_WP) self%grd_y(:,i,j,k)=0.0_WP
             if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i,j,k-1).eq.0.0_WP) self%grd_z(:,i,j,k)=0.0_WP
          end do
       end do
    end do

    ! Adjust metrics to account for lower dimensionality
    if (self%cfg%nx.eq.1) then
       self%div_x=0.0_WP
       self%grd_x=0.0_WP
    end if
    if (self%cfg%ny.eq.1) then
       self%div_y=0.0_WP
       self%grd_y=0.0_WP
    end if
    if (self%cfg%nz.eq.1) then
       self%div_z=0.0_WP
       self%grd_z=0.0_WP
    end if

    ! Create implicit solver object for filtering
    self%implicit=ddadi(cfg=self%cfg,name='Filter',nst=7)
    self%implicit%stc(1,:)=[ 0, 0, 0]
    self%implicit%stc(2,:)=[+1, 0, 0]
    self%implicit%stc(3,:)=[-1, 0, 0]
    self%implicit%stc(4,:)=[ 0,+1, 0]
    self%implicit%stc(5,:)=[ 0,-1, 0]
    self%implicit%stc(6,:)=[ 0, 0,+1]
    self%implicit%stc(7,:)=[ 0, 0,-1]
    call self%implicit%init()

    ! Set filter width to zero by default
    self%filter_width=0.0_WP

    ! Solve implicitly by default
    self%implicit_filter=.true.

    ! Generate a wall data for collisions
    allocate(self%xwall(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%xwall=huge(1.0_WP)
    allocate(self%ywall(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%ywall=huge(1.0_WP)
    allocate(self%zwall(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%zwall=huge(1.0_WP)
    ! Get closest wall location in x
    do k=self%cfg%kmino_,self%cfg%kmaxo_
       do j=self%cfg%jmino_,self%cfg%jmaxo_
          do i=self%cfg%imin_,self%cfg%imax_
             do l=-1,+2
                if ((self%cfg%VF(i+l-1,j,k).eq.0.0_WP.and.self%cfg%VF(i+l,j,k).eq.1.0_WP).or.&
                     &   (self%cfg%VF(i+l-1,j,k).eq.1.0_WP.and.self%cfg%VF(i+l,j,k).eq.0.0_WP)) then
                   self%xwall(i,j,k)=min(self%xwall(i,j,k),self%cfg%x(i+l))
                end if
             end do
          end do
       end do
    end do
    call self%cfg%sync(self%xwall)
    ! Get closest wall location in y
    do k=self%cfg%kmino_,self%cfg%kmaxo_
       do j=self%cfg%jmin_,self%cfg%jmax_
          do i=self%cfg%imino_,self%cfg%imaxo_
             do l=-1,+2
                if ((self%cfg%VF(i,j+l-1,k).eq.0.0_WP.and.self%cfg%VF(i,j+l,k).eq.1.0_WP).or.&
                     &   (self%cfg%VF(i,j+l-1,k).eq.1.0_WP.and.self%cfg%VF(i,j+l,k).eq.0.0_WP)) then
                   self%ywall(i,j,k)=min(self%ywall(i,j,k),self%cfg%y(j+l))
                end if
             end do
          end do
       end do
    end do
    call self%cfg%sync(self%ywall)
    ! Get closest wall location in z
    do k=self%cfg%kmin_,self%cfg%kmax_
       do j=self%cfg%jmino_,self%cfg%jmaxo_
          do i=self%cfg%imino_,self%cfg%imaxo_
             do l=-1,+2
                if ((self%cfg%VF(i,j,k+l-1).eq.0.0_WP.and.self%cfg%VF(i,j,k+l).eq.1.0_WP).or.&
                     &   (self%cfg%VF(i,j,k+l-1).eq.1.0_WP.and.self%cfg%VF(i,j,k+l).eq.0.0_WP)) then
                   self%zwall(i,j,k)=min(self%zwall(i,j,k),self%cfg%z(k+l))
                end if
             end do
          end do
       end do
    end do
    call self%cfg%sync(self%zwall)

    ! Log/screen output
    logging: block
      use, intrinsic :: iso_fortran_env, only: output_unit
      use param,    only: verbose
      use messager, only: log
      use string,   only: str_long
      character(len=str_long) :: message
      if (self%cfg%amRoot) then
         write(message,'("Particle solver [",a,"] on partitioned grid [",a,"]")') trim(self%name),trim(self%cfg%name)
         if (verbose.gt.1) write(output_unit,'(a)') trim(message)
         if (verbose.gt.0) call log(message)
      end if
    end block logging

  end function constructor
  

  !> Resolve collisional interaction between particles, walls, and an optional IB level set
  !> Requires tau_col, e_n, e_w and mu_f to be set beforehand
  subroutine collide(this,dt,Gib,Nxib,Nyib,Nzib)
   implicit none
   class(multisphere), intent(inout) :: this
   real(WP), intent(inout) :: dt  !< Timestep size over which to advance
   real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: Gib  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
   real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: Nxib !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
   real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: Nyib !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
   real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: Nzib !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
   integer, dimension(:,:,:), allocatable :: npic      !< Number of particle in cell
   integer, dimension(:,:,:,:), allocatable :: ipic    !< Index of particle in cell
   integer, dimension(:,:,:,:), allocatable :: jpic    !< Index of body in cell

   ! Check if all IB parameters are present
   check_G: block
     use messager, only: die
     if (present(Gib).and.(.not.present(Nxib).or..not.present(Nyib).or..not.present(Nzib))) &
          call die('[multisphere collide] IB collisions need Gib, Nxib, Nyib, AND Nzib')
   end block check_G

   ! Start by zeroing out the collision force
   zero_force: block
     integer :: j
     do j=1,this%nb_
        this%b(j)%p_col=0.0_WP
     end do
   end block zero_force
   
   ! Then share bodies across overlap
   call this%share(ceiling(2.0_WP*this%max_R0/this%cfg%min_meshsize))

   ! We can now assemble particle-in-cell information
   pic_prep: block
     use mpi_f08
     integer :: i,j,ip,jp,kp,ierr
     integer :: mymax_npic,max_npic

     ! Allocate number of particle in cell
     allocate(npic(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); npic=0

     ! Count particles and ghosts per cell
     do j=1,this%nb_
        do i=1,this%b(j)%np
           ip=this%b(j)%p_ind(1,i); jp=this%b(j)%p_ind(2,i); kp=this%b(j)%p_ind(3,i)
           npic(ip,jp,kp)=npic(ip,jp,kp)+1
        end do
     end do
     do j=1,this%ng_
        do i=1,this%g(j)%np
           ip=this%g(j)%p_ind(1,i); jp=this%g(j)%p_ind(2,i); kp=this%g(j)%p_ind(3,i)
           if ( ip.ge.this%cfg%imino_.and.ip.le.this%cfg%imaxo_.and.&
                jp.ge.this%cfg%jmino_.and.jp.le.this%cfg%jmaxo_.and.&
                kp.ge.this%cfg%kmino_.and.kp.le.this%cfg%kmaxo_) npic(ip,jp,kp)=npic(ip,jp,kp)+1
        end do
     end do

     ! Get maximum number of particle in cell
     mymax_npic=maxval(npic); call MPI_ALLREDUCE(mymax_npic,max_npic,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)

     ! Allocate pic map
     allocate(ipic(1:max_npic,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); ipic=0
     allocate(jpic(1:max_npic,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); jpic=0

     ! Assemble pic map
     npic=0
     do j=1,this%nb_
        do i=1,this%b(j)%np
           ip=this%b(j)%p_ind(1,i); jp=this%b(j)%p_ind(2,i); kp=this%b(j)%p_ind(3,i)
           npic(ip,jp,kp)=npic(ip,jp,kp)+1
           ipic(npic(ip,jp,kp),ip,jp,kp)=i
           jpic(npic(ip,jp,kp),ip,jp,kp)=j
        end do
     end do
     do j=1,this%ng_
        do i=1,this%g(j)%np
           ip=this%g(j)%p_ind(1,i); jp=this%g(j)%p_ind(2,i); kp=this%g(j)%p_ind(3,i)
           if ( ip.ge.this%cfg%imino_.and.ip.le.this%cfg%imaxo_.and.&
                jp.ge.this%cfg%jmino_.and.jp.le.this%cfg%jmaxo_.and.&
                kp.ge.this%cfg%kmino_.and.kp.le.this%cfg%kmaxo_) then
              npic(ip,jp,kp)=npic(ip,jp,kp)+1
              ipic(npic(ip,jp,kp),ip,jp,kp)=-i
              jpic(npic(ip,jp,kp),ip,jp,kp)=-j
           end if
        end do
     end do

   end block pic_prep

   ! Calculate collision force
   collision_force: block
     use mpi_f08
     use mathtools, only: Pi,normalize,cross_product
     integer :: i1,i2,j1,j2,ii,jj,kk,nn,ierr
     real(WP) :: d1,m1,d2,m2,d12,m12,buf
     real(WP), dimension(3) :: r1,v1,w1,r2,v2,w2,v12,n12,f_n,t12,f_t
     real(WP) :: k_n,eta_n,k_coeff,eta_coeff,k_coeff_w,eta_coeff_w,rnv,r_influ,delta_n,rtv
     real(WP), parameter :: aclipnorm=1.0e-6_WP
     real(WP), parameter :: acliptan=1.0e-9_WP
     real(WP), parameter :: rcliptan=0.05_WP

     ! Reset collision counter
     this%ncol=0

     ! Precompute coefficients for k and eta
     k_coeff=(Pi**2+log(this%e_n)**2)/this%tau_col**2
     eta_coeff=-2.0_WP*log(this%e_n)/this%tau_col
     k_coeff_w=(Pi**2+log(this%e_w)**2)/this%tau_col**2
     eta_coeff_w=-2.0_WP*log(this%e_w)/this%tau_col

     ! Loop over all local bodies and their particles
     collision: do j1=1,this%nb_
        do i1=1,this%b(j1)%np

           ! Store particle data
           r1=this%b(j1)%p_pos(:,i1)
           v1=this%b(j1)%p_vel(:,i1)
           w1=this%b(j1)%angVel
           d1=this%b(j1)%p_d(i1)
           m1=this%rho*0.25_WP*Pi*d1**2*this%Hp

           ! Collide with walls in x
           d12=abs(this%xwall(this%b(j1)%p_ind(1,i1),this%b(j1)%p_ind(2,i1),this%b(j1)%p_ind(3,i1))-this%b(j1)%p_pos(1,i1))
           n12=[sign(1.0_WP,this%xwall(this%b(j1)%p_ind(1,i1),this%b(j1)%p_ind(2,i1),this%b(j1)%p_ind(3,i1))-this%b(j1)%p_pos(1,i1)),0.0_WP,0.0_WP]
           rnv=dot_product(v1,n12)
           r_influ=min(2.0_WP*abs(rnv)*dt,0.2_WP*d1)
           delta_n=min(0.5_WP*d1+r_influ-d12,this%clip_col*0.5_WP*d1)
           ! Assess if there is collision
           if (delta_n.gt.0.0_WP) then
              ! Normal collision
              k_n=m1*k_coeff_w
              eta_n=m1*eta_coeff_w
              f_n=-k_n*delta_n*n12-eta_n*rnv*n12
              ! Tangential collision
              f_t=0.0_WP
              if (this%mu_f.gt.0.0_WP) then
                 t12 = v1-rnv*n12+cross_product(0.5_WP*d1*w1,n12)
                 rtv = sqrt(sum(t12*t12))
                 if (rnv*dt/d1.gt.aclipnorm) then
                    if (rtv/rnv.lt.rcliptan) rtv=0.0_WP
                 else
                    if (rtv*dt/d1.lt.acliptan) rtv=0.0_WP
                 end if
                 if (rtv.gt.0.0_WP) f_t=-this%mu_f*sqrt(sum(f_n*f_n))*t12/rtv
              end if
              ! Update collision force
              this%b(j1)%p_col(:,i1)=this%b(j1)%p_col(:,i1)+f_n+f_t
           end if

           ! Collide with walls in y
           d12=abs(this%ywall(this%b(j1)%p_ind(1,i1),this%b(j1)%p_ind(2,i1),this%b(j1)%p_ind(3,i1))-this%b(j1)%p_pos(2,i1))
           n12=[0.0_WP,sign(1.0_WP,this%ywall(this%b(j1)%p_ind(1,i1),this%b(j1)%p_ind(2,i1),this%b(j1)%p_ind(3,i1))-this%b(j1)%p_pos(2,i1)),0.0_WP]
           rnv=dot_product(v1,n12)
           r_influ=min(2.0_WP*abs(rnv)*dt,0.2_WP*d1)
           delta_n=min(0.5_WP*d1+r_influ-d12,this%clip_col*0.5_WP*d1)
           ! Assess if there is collision
           if (delta_n.gt.0.0_WP) then
              ! Normal collision
              k_n=m1*k_coeff_w
              eta_n=m1*eta_coeff_w
              f_n=-k_n*delta_n*n12-eta_n*rnv*n12
              ! Tangential collision
              f_t=0.0_WP
              if (this%mu_f.gt.0.0_WP) then
                 t12 = v1-rnv*n12+cross_product(0.5_WP*d1*w1,n12)
                 rtv = sqrt(sum(t12*t12))
                 if (rnv*dt/d1.gt.aclipnorm) then
                    if (rtv/rnv.lt.rcliptan) rtv=0.0_WP
                 else
                    if (rtv*dt/d1.lt.acliptan) rtv=0.0_WP
                 end if
                 if (rtv.gt.0.0_WP) f_t=-this%mu_f*sqrt(sum(f_n*f_n))*t12/rtv
              end if
              ! Update collision force
              this%b(j1)%p_col(:,i1)=this%b(j1)%p_col(:,i1)+f_n+f_t
           end if

           ! Collide with walls in z
           d12=abs(this%zwall(this%b(j1)%p_ind(1,i1),this%b(j1)%p_ind(2,i1),this%b(j1)%p_ind(3,i1))-this%b(j1)%p_pos(3,i1))
           n12=[0.0_WP,0.0_WP,sign(1.0_WP,this%zwall(this%b(j1)%p_ind(1,i1),this%b(j1)%p_ind(2,i1),this%b(j1)%p_ind(3,i1))-this%b(j1)%p_pos(3,i1))]
           rnv=dot_product(v1,n12)
           r_influ=min(2.0_WP*abs(rnv)*dt,0.2_WP*d1)
           delta_n=min(0.5_WP*d1+r_influ-d12,this%clip_col*0.5_WP*d1)
           ! Assess if there is collision
           if (delta_n.gt.0.0_WP) then
              ! Normal collision
              k_n=m1*k_coeff_w
              eta_n=m1*eta_coeff_w
              f_n=-k_n*delta_n*n12-eta_n*rnv*n12
              ! Tangential collision
              f_t=0.0_WP
              if (this%mu_f.gt.0.0_WP) then
                 t12 = v1-rnv*n12+cross_product(0.5_WP*d1*w1,n12)
                 rtv = sqrt(sum(t12*t12))
                 if (rnv*dt/d1.gt.aclipnorm) then
                    if (rtv/rnv.lt.rcliptan) rtv=0.0_WP
                 else
                    if (rtv*dt/d1.lt.acliptan) rtv=0.0_WP
                 end if
                 if (rtv.gt.0.0_WP) f_t=-this%mu_f*sqrt(sum(f_n*f_n))*t12/rtv
              end if
              ! Update collision force
              this%b(j1)%p_col(:,i1)=this%b(j1)%p_col(:,i1)+f_n+f_t
           end if

           ! Collide with IB
           if (present(Gib)) then
              d12=this%cfg%get_scalar(pos=this%b(j1)%p_pos(:,i1),i0=this%b(j1)%p_ind(1,i1),j0=this%b(j1)%p_ind(2,i1),k0=this%b(j1)%p_ind(3,i1),S=Gib,bc='n')
              n12(1)=this%cfg%get_scalar(pos=this%b(j1)%p_pos(:,i1),i0=this%b(j1)%p_ind(1,i1),j0=this%b(j1)%p_ind(2,i1),k0=this%b(j1)%p_ind(3,i1),S=Nxib,bc='n')
              n12(2)=this%cfg%get_scalar(pos=this%b(j1)%p_pos(:,i1),i0=this%b(j1)%p_ind(1,i1),j0=this%b(j1)%p_ind(2,i1),k0=this%b(j1)%p_ind(3,i1),S=Nyib,bc='n')
              n12(3)=this%cfg%get_scalar(pos=this%b(j1)%p_pos(:,i1),i0=this%b(j1)%p_ind(1,i1),j0=this%b(j1)%p_ind(2,i1),k0=this%b(j1)%p_ind(3,i1),S=Nzib,bc='n')
              buf = sqrt(sum(n12*n12))+epsilon(1.0_WP)
              n12 = -n12/buf
              rnv=dot_product(v1,n12)
              r_influ=min(2.0_WP*abs(rnv)*dt,0.2_WP*d1)
              delta_n=min(0.5_WP*d1+r_influ-d12,this%clip_col*0.5_WP*d1)

              ! Assess if there is collision
              if (delta_n.gt.0.0_WP) then
                 ! Normal collision
                 k_n=m1*k_coeff_w
                 eta_n=m1*eta_coeff_w
                 f_n=-k_n*delta_n*n12-eta_n*rnv*n12
                 ! Tangential collision
                 f_t=0.0_WP
                 if (this%mu_f.gt.0.0_WP) then
                    t12 = v1-rnv*n12+cross_product(0.5_WP*d1*w1,n12)
                    rtv = sqrt(sum(t12*t12))
                    if (rnv*dt/d1.gt.aclipnorm) then
                       if (rtv/rnv.lt.rcliptan) rtv=0.0_WP
                    else
                       if (rtv*dt/d1.lt.acliptan) rtv=0.0_WP
                    end if
                    if (rtv.gt.0.0_WP) f_t=-this%mu_f*sqrt(sum(f_n*f_n))*t12/rtv
                 end if
                 ! Update collision force
                 this%b(j1)%p_col(:,i1)=this%b(j1)%p_col(:,i1)+f_n+f_t
              end if
           end if

           ! Loop over nearest cells
           do kk=this%b(j1)%p_ind(3,i1)-1,this%b(j1)%p_ind(3,i1)+1
              do jj=this%b(j1)%p_ind(2,i1)-1,this%b(j1)%p_ind(2,i1)+1
                 do ii=this%b(j1)%p_ind(1,i1)-1,this%b(j1)%p_ind(1,i1)+1

                    ! Loop over particles in that cell
                    do nn=1,npic(ii,jj,kk)

                       ! Get index of neighbor body and particle
                       j2=jpic(nn,ii,jj,kk)
                       i2=ipic(nn,ii,jj,kk)

                       ! Get relevant data from correct storage
                       if (i2.gt.0) then
                          r2=this%b(j2)%p_pos(:,i2)
                          v2=this%b(j2)%p_vel(:,i2)
                          w2=this%b(j2)%angVel
                          d2=this%b(j2)%p_d(i2)
                          m2=this%rho*0.25_WP*Pi*d2**2*this%Hp
                       else if (i2.lt.0) then
                          i2=-i2; j2=-j2
                          r2=this%g(j2)%p_pos(:,i2)
                          v2=this%g(j2)%p_vel(:,i2)
                          w2=this%g(j2)%angVel
                          d2=this%g(j2)%p_d(i2)
                          m2=this%rho*0.25_WP*Pi*d2**2*this%Hp
                       end if

                       ! Avoid collision with particles from same body
                       if (j1.ne.j2) then

                          ! Compute relative information
                          d12=norm2(r1-r2)
                          if (d12.lt.10.0_WP*epsilon(d12)) cycle !< this should skip auto-collision
                          n12=(r2-r1)/d12
                          v12=v1-v2
                          rnv=dot_product(v12,n12)
                          r_influ=min(abs(rnv)*dt,0.1_WP*(d1+d2))
                          delta_n=min(0.5_WP*(d1+d2)+r_influ-d12,this%clip_col*0.5_WP*(d1+d2))

                          ! Assess if there is collision (avoid collision with self)
                          if (delta_n.gt.0.0_WP) then
                             ! Normal collision
                             m12=m1*m2/(m1+m2)
                             k_n=m12*k_coeff
                             eta_n=m12*eta_coeff
                             f_n=-k_n*delta_n*n12-eta_n*rnv*n12
                             ! Tangential collision
                             f_t=0.0_WP
                             if (this%mu_f.gt.0.0_WP) then
                                t12 = v12-rnv*n12+cross_product(0.5_WP*(d1*w1+d2*w2),n12)
                                rtv = sqrt(sum(t12*t12))
                                if (rnv*dt*2.0_WP/(d1+d2).gt.aclipnorm) then
                                   if (rtv/rnv.lt.rcliptan) rtv=0.0_WP
                                else
                                   if (rtv*dt*2.0_WP/(d1+d2).lt.acliptan) rtv=0.0_WP
                                end if
                                if (rtv.gt.0.0_WP) f_t=-this%mu_f*sqrt(sum(f_n*f_n))*t12/rtv
                             end if
                             ! Update the collision force
                             this%b(j1)%p_col(:,i1)=this%b(j1)%p_col(:,i1)+f_n+f_t
                             ! Add up the collisions
                             this%ncol=this%ncol+1
                          end if
                       end if

                    end do

                 end do
              end do
           end do

           ! Deal with dimensionality
           if (this%cfg%nx.eq.1) this%b(j1)%p_col(1,i1)=0.0_WP
           if (this%cfg%ny.eq.1) this%b(j1)%p_col(2,i1)=0.0_WP
           if (this%cfg%nz.eq.1) this%b(j1)%p_col(3,i1)=0.0_WP

        end do
     end do collision

     ! Determine total number of collisions
     call MPI_ALLREDUCE(this%ncol,nn,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr); this%ncol=nn/2

   end block collision_force

   ! Clean up
   if (allocated(npic)) deallocate(npic)
   if (allocated(ipic)) deallocate(ipic)
   if (allocated(jpic)) deallocate(jpic)

 end subroutine collide
  

 !> Advance the particle equations by a specified time step dt
  subroutine advance(this,dt,rho,visc,U,V,W,dpdx,dpdy,dpdz)
    use mathtools, only: Pi,twoPi
    use mpi_f08, only : MPI_SUM,MPI_INTEGER
    implicit none
    class(multisphere), intent(inout) :: this
    real(WP), intent(inout) :: dt,rho,visc
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: dpdx  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: dpdy  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: dpdz  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    integer :: i,j
    real(WP), dimension(3) :: r
    real(WP), dimension(2,2) :: R2D
    type(body) :: myb,bold
    
    do j=1,this%nb_
       ! Create local copy of body
       myb=this%b(j)
       ! Remember the particle
       bold=myb
       ! Advance with Euler prediction
       call this%get_force(dt=dt,rho=rho,visc=visc,U=U,V=V,W=W,dpdx=dpdx,dpdy=dpdy,dpdz=dpdz,b=myb)
       myb%pos=bold%pos+0.5_WP*dt*myb%vel
       myb%vel=bold%vel+0.5_WP*dt*myb%force/myb%m
       myb%theta=bold%theta+0.5_WP*dt*myb%angVel(3)
       myb%theta=mod(myb%theta,twoPi)
       myb%angVel(1:2)=0.0_WP
       myb%angVel(3)=bold%angVel(3)+0.5_WP*dt*myb%torque(3)/myb%I
       ! Build 2D rotation matrix from theta
       R2D(1,1)=cos(myb%theta); R2D(1,2)=-sin(myb%theta)
       R2D(2,1)=sin(myb%theta); R2D(2,2)= cos(myb%theta)
       ! Update particle positions and velocities using rotated r0
       do i=1,myb%np
          ! Rotate initial relative position vector
          r(1:2) = matmul(R2D,myb%p_r0(1:2,i))
          ! Update position based on rotated relative position
          myb%p_pos(1:2,i)=myb%pos(1:2)+r(1:2)
          ! Update particle velocity from rigid body motion
          myb%p_vel(1,i)=myb%vel(1)-myb%angVel(3)*r(2)
          myb%p_vel(2,i)=myb%vel(2)+myb%angVel(3)*r(1)
          myb%p_vel(3,i)=0.0_WP
       end do
       ! Correct with midpoint rule
       call this%get_force(dt=dt,rho=rho,visc=visc,U=U,V=V,W=W,dpdx=dpdx,dpdy=dpdy,dpdz=dpdz,b=myb)
       myb%pos=bold%pos+dt*myb%vel
       myb%vel=bold%vel+dt*myb%force/myb%m
       myb%theta=bold%theta+dt*myb%angVel(3)
       myb%theta=mod(myb%theta,twoPi)
       myb%angVel(1:2)=0.0_WP
       myb%angVel(3)=bold%angVel(3)+dt*myb%torque(3)/myb%I
       ! Build 2D rotation matrix from theta
       R2D(1,1)=cos(myb%theta); R2D(1,2)=-sin(myb%theta)
       R2D(2,1)=sin(myb%theta); R2D(2,2)= cos(myb%theta)
       ! Update particle positions and velocities using rotated r0
       do i=1,myb%np
          ! Rotate initial relative position vector
          r(1:2) = matmul(R2D,myb%p_r0(1:2,i))
          ! Update position based on rotated relative position
          myb%p_pos(1:2,i)=myb%pos(1:2)+r(1:2)
          ! Update particle velocity from rigid body motion
          myb%p_vel(1,i)=myb%vel(1)-myb%angVel(3)*r(2)
          myb%p_vel(2,i)=myb%vel(2)+myb%angVel(3)*r(1)
          myb%p_vel(3,i)=0.0_WP
       end do
       ! Relocalize
       myb%ind=this%cfg%get_ijk_global(myb%pos,myb%ind)
       do i=1,myb%np
          myb%p_ind(:,i)=this%cfg%get_ijk_global(myb%p_pos(:,i),myb%p_ind(:,i))
       end do
       ! Correct the position to take into account periodicity
       if (this%cfg%xper) myb%pos(1)=this%cfg%x(this%cfg%imin)+modulo(myb%pos(1)-this%cfg%x(this%cfg%imin),this%cfg%xL)
       if (this%cfg%yper) myb%pos(2)=this%cfg%y(this%cfg%jmin)+modulo(myb%pos(2)-this%cfg%y(this%cfg%jmin),this%cfg%yL)
       if (this%cfg%zper) myb%pos(3)=this%cfg%z(this%cfg%kmin)+modulo(myb%pos(3)-this%cfg%z(this%cfg%kmin),this%cfg%zL)
       ! Handle particles that have left the domain
       if (myb%pos(1).lt.this%cfg%x(this%cfg%imin).or.myb%pos(1).gt.this%cfg%x(this%cfg%imax+1)) myb%flag=1
       if (myb%pos(2).lt.this%cfg%y(this%cfg%jmin).or.myb%pos(2).gt.this%cfg%y(this%cfg%jmax+1)) myb%flag=1
       if (myb%pos(3).lt.this%cfg%z(this%cfg%kmin).or.myb%pos(3).gt.this%cfg%z(this%cfg%kmax+1)) myb%flag=1
       ! Relocalize
       myb%ind=this%cfg%get_ijk_global(myb%pos,myb%ind)
       do i=1,myb%np
          myb%p_ind(:,i)=this%cfg%get_ijk_global(myb%p_pos(:,i),myb%p_ind(:,i))
       end do
       ! Copy back to body
       this%b(j)=myb
    end do
    ! Communicate particles
    call this%sync()

    ! Recompute volume fraction
    call this%update_VF()

    ! Log/screen output
    logging: block
      use, intrinsic :: iso_fortran_env, only: output_unit
      use param,    only: verbose
      use messager, only: log
      use string,   only: str_long
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Particle solver [",a,"] on partitioned grid [",a,"]: ",i0," particles were advanced")') trim(this%name),trim(this%cfg%name),this%np
         if (verbose.gt.1) write(output_unit,'(a)') trim(message)
         if (verbose.gt.0) call log(message)
      end if
    end block logging

  end subroutine advance


  !> Compute net force and torque on the rigid bodies
  subroutine get_force(this,dt,rho,visc,U,V,W,dpdx,dpdy,dpdz,b)
    use mathtools, only: Pi,cross_product
    implicit none
    class(multisphere), intent(inout) :: this
    real(WP), intent(inout) :: dt,rho,visc
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: dpdx  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: dpdy  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: dpdz  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    type(body), intent(inout) :: b
    integer :: i
    real(WP) :: m,pVF,fVF
    real(WP), dimension(3) :: acc,F,gradP
    b%force=0.0_WP; b%torque=0.0_WP
    do i=1,b%np
       ! Interpolate fluid quantities to particle location
       gradP=this%cfg%get_velocity(pos=b%p_pos(:,i),i0=b%p_ind(1,i),j0=b%p_ind(2,i),k0=b%p_ind(3,i),U=dpdx,V=dpdy,W=dpdz)
       pVF=this%cfg%get_scalar(pos=b%p_pos(:,i),i0=b%p_ind(1,i),j0=b%p_ind(2,i),k0=b%p_ind(3,i),S=this%VF,bc='n')
       fVF=1.0_WP-pVF
       ! Compute acceleration due to pressure
       acc=-gradP/this%rho
       ! Update rigid body force and torque (only in 2D)
       m=this%rho*0.25_WP*Pi*b%p_d(i)**2*this%Hp
       F=m*acc+b%p_col(:,i)
       b%force(1:2)=b%force(1:2)+F(1:2)
       b%torque=b%torque+cross_product(b%p_pos(:,i)-b%pos,F)
    end do
  end subroutine get_force


  !> Update particle volume fraction using our current particles
  subroutine update_VF(this)
    use mathtools, only: Pi
    implicit none
    class(multisphere), intent(inout) :: this
    integer :: i,j
    real(WP) :: Vp
    ! Reset volume fraction
    this%VF=0.0_WP
    ! Transfer particle volume
    do j=1,this%nb_
       ! Skip inactive body
       if (this%b(j)%flag.eq.1) cycle
       do i=1,this%b(j)%np
          ! Transfer particle volume
          Vp=0.25_WP*Pi*this%b(j)%p_d(i)**2*this%Hp
          call this%cfg%set_scalar(Sp=Vp,pos=this%b(j)%p_pos(:,i),i0=this%b(j)%p_ind(1,i),j0=this%b(j)%p_ind(2,i),k0=this%b(j)%p_ind(3,i),S=this%VF,bc='n')
       end do
    end do
    this%VF=this%VF/this%cfg%vol
    ! Sum at boundaries
    call this%cfg%syncsum(this%VF)
    ! Apply volume filter
    call this%filter(this%VF)
  end subroutine update_VF


  !> Update properties based on rigid body mechanics (should be called only once)
  subroutine update_rigid_body(this)
    use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
    use parallel, only: MPI_REAL_WP
    use mathtools, only: Pi
    implicit none
    class(multisphere), intent(inout) :: this
    integer :: i,j,ierr
    real(WP) :: m
    do j=1,this%nb_
       ! Compute center-of-mass
       this%b(j)%m=0.0_WP; this%b(j)%pos=0.0_WP
       do i=1,this%b(j)%np
          m=this%rho*0.25_WP*Pi*this%b(j)%p_d(i)**2*this%Hp
          this%b(j)%m=this%b(j)%m+m
          this%b(j)%pos=this%b(j)%pos+m*this%b(j)%p_pos(:,i)
       end do
       this%b(j)%pos=this%b(j)%pos/this%b(j)%m
       ! Initialize relative position vector and moment of inertia
       this%b(j)%I=0.0_WP
       this%max_R0=0.0_WP
       do i=1,this%b(j)%np
          m=this%rho*0.25_WP*Pi*this%b(j)%p_d(i)**2*this%Hp
          ! Compute relative position r0
          this%b(j)%p_r0(:,i)=this%b(j)%p_pos(:,i)-this%b(j)%pos
          ! Update moment of inertia
          this%b(j)%I=this%b(j)%I+m*(this%b(j)%p_r0(1,i)**2+this%b(j)%p_r0(2,i)**2)
          ! Update bounding radius
          this%max_R0=max(this%max_R0,sqrt(sum(this%b(j)%p_r0(1:2,i)**2)))
          ! Update particle velocity based on rigid body motion
          this%b(j)%p_vel(1,i)=this%b(j)%vel(1)-this%b(j)%angVel(3)*this%b(j)%p_r0(2,i)
          this%b(j)%p_vel(2,i)=this%b(j)%vel(2)+this%b(j)%angVel(3)*this%b(j)%p_r0(1,i)
          this%b(j)%p_vel(3,i)=0.0_WP
       end do
    end do    
    call MPI_ALLREDUCE(this%max_R0,m,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%max_R0=m
  end subroutine update_rigid_body


  !> Laplacian filtering operation
  subroutine filter(this,A)
    implicit none
    class(multisphere), intent(inout) :: this
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: A     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP) :: filter_coeff
    integer :: i,j,k,n,nstep
    real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ

    ! Return without filtering if filter width is zero
    if (this%filter_width.le.0.0_WP) return

    ! Recompute filter coefficient
    filter_coeff=max(this%filter_width**2-this%cfg%min_meshsize**2,0.0_WP)/(16.0_WP*log(2.0_WP))
    if (filter_coeff.le.0.0_WP) return

    ! Allocate flux arrays
    allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
    allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
    allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))

    if (this%implicit_filter) then  !< Apply filter implicitly
       if (.not.this%implicit%setup_done) then
          ! Prepare diffusive operator (only need to do this once)
          do k=this%cfg%kmin_,this%cfg%kmax_
             do j=this%cfg%jmin_,this%cfg%jmax_
                do i=this%cfg%imin_,this%cfg%imax_
                   this%implicit%opr(1,i,j,k)=1.0_WP-(this%div_x(+1,i,j,k)*filter_coeff*this%grd_x(-1,i+1,j,k)+&
                   &                                  this%div_x( 0,i,j,k)*filter_coeff*this%grd_x( 0,i  ,j,k)+&
                   &                                  this%div_y(+1,i,j,k)*filter_coeff*this%grd_y(-1,i,j+1,k)+&
                   &                                  this%div_y( 0,i,j,k)*filter_coeff*this%grd_y( 0,i,j  ,k)+&
                   &                                  this%div_z(+1,i,j,k)*filter_coeff*this%grd_z(-1,i,j,k+1)+&
                   &                                  this%div_z( 0,i,j,k)*filter_coeff*this%grd_z( 0,i,j,k  ))
                   this%implicit%opr(2,i,j,k)=      -(this%div_x(+1,i,j,k)*filter_coeff*this%grd_x( 0,i+1,j,k))
                   this%implicit%opr(3,i,j,k)=      -(this%div_x( 0,i,j,k)*filter_coeff*this%grd_x(-1,i  ,j,k))
                   this%implicit%opr(4,i,j,k)=      -(this%div_y(+1,i,j,k)*filter_coeff*this%grd_y( 0,i,j+1,k))
                   this%implicit%opr(5,i,j,k)=      -(this%div_y( 0,i,j,k)*filter_coeff*this%grd_y(-1,i,j  ,k))
                   this%implicit%opr(6,i,j,k)=      -(this%div_z(+1,i,j,k)*filter_coeff*this%grd_z( 0,i,j,k+1))
                   this%implicit%opr(7,i,j,k)=      -(this%div_z( 0,i,j,k)*filter_coeff*this%grd_z(-1,i,j,k  ))
                end do
             end do
          end do
       end if
       ! Explicit step
       do k=this%cfg%kmin_,this%cfg%kmax_+1
          do j=this%cfg%jmin_,this%cfg%jmax_+1
             do i=this%cfg%imin_,this%cfg%imax_+1
                FX(i,j,k)=filter_coeff*sum(this%grd_x(:,i,j,k)*A(i-1:i,j,k))
                FY(i,j,k)=filter_coeff*sum(this%grd_y(:,i,j,k)*A(i,j-1:j,k))
                FZ(i,j,k)=filter_coeff*sum(this%grd_z(:,i,j,k)*A(i,j,k-1:k))
             end do
          end do
       end do
       do k=this%cfg%kmin_,this%cfg%kmax_
          do j=this%cfg%jmin_,this%cfg%jmax_
             do i=this%cfg%imin_,this%cfg%imax_
                this%implicit%rhs(i,j,k)=sum(this%div_x(:,i,j,k)*FX(i:i+1,j,k))+sum(this%div_y(:,i,j,k)*FY(i,j:j+1,k))+sum(this%div_z(:,i,j,k)*FZ(i,j,k:k+1))
             end do
          end do
       end do
       ! Implicit step
       call this%implicit%setup()
       this%implicit%sol=0.0_WP
       call this%implicit%solve()
       A=A+this%implicit%sol
       call this%cfg%sync(A)
       
    else  !< Apply filter explicitly
       nstep=ceiling(6.0_WP*filter_coeff/this%cfg%min_meshsize**2)
       filter_coeff=filter_coeff/real(nstep,WP)
       do n=1,nstep
          ! Diffusive flux of A
          do k=this%cfg%kmin_,this%cfg%kmax_+1
             do j=this%cfg%jmin_,this%cfg%jmax_+1
                do i=this%cfg%imin_,this%cfg%imax_+1
                   FX(i,j,k)=filter_coeff*sum(this%grd_x(:,i,j,k)*A(i-1:i,j,k))
                   FY(i,j,k)=filter_coeff*sum(this%grd_y(:,i,j,k)*A(i,j-1:j,k))
                   FZ(i,j,k)=filter_coeff*sum(this%grd_z(:,i,j,k)*A(i,j,k-1:k))
                end do
             end do
          end do
          ! Divergence of fluxes
          do k=this%cfg%kmin_,this%cfg%kmax_
             do j=this%cfg%jmin_,this%cfg%jmax_
                do i=this%cfg%imin_,this%cfg%imax_
                   A(i,j,k)=A(i,j,k)+sum(this%div_x(:,i,j,k)*FX(i:i+1,j,k))+sum(this%div_y(:,i,j,k)*FY(i,j:j+1,k))+sum(this%div_z(:,i,j,k)*FZ(i,j,k:k+1))
                end do
             end do
          end do
          ! Sync A
          call this%cfg%sync(A)
       end do
    end if

    ! Deallocate flux arrays
    deallocate(FX,FY,FZ)

  end subroutine filter


  !> Read in number of particles, center of mass (x0,y0), angle (theta), diameter (dp)
  !> and length (L), return individual particle positions
  function get_positions(n,x0,y0,theta,L,shape) result(pos)
    use mathtools, only: Pi
    use messager, only: die
    implicit none
    integer, intent(in) :: n
    real(WP), intent(in) :: x0,y0,theta,L
    character(len=*), intent(in) :: shape
    real(WP), dimension(2,n) :: pos
    integer :: i,j,k,total_side_particles
    real(WP) :: dx,x_local,x_rot,y_rot
    real(WP) :: angle,R,frac
    real(WP) :: x_start,y_start,x_end,y_end,x_seg,y_seg
    real(WP), dimension(5,2) :: verts
    select case (trim(adjustl(shape)))
    case ('rod')
       ! ========== ROD ==========
       if (n.le.1) then
          pos(1,1) = x0
          pos(2,1) = y0
          return
       end if
       dx=L/real(n-1,WP)
       do i=1,n
          x_local=-L/2.0_WP+dx*real(i-1,WP)
          x_rot=cos(theta)*x_local
          y_rot=sin(theta)*x_local
          pos(1,i) = y0 + x_rot
          pos(2,i) = x0 + y_rot
       end do
    case ('pentagon')
       ! ========== PENTAGON ==========
       ! Radius of circumscribed circle
       R=L/(2.0_WP*sin(Pi/5.0_WP))
       ! Compute 5 vertices of pentagon
       do k=1,5
          angle=theta+2.0_WP*pi*real(k-1,WP)/5.0_WP
          verts(k,1)=R*cos(angle)
          verts(k,2)=R*sin(angle)
       end do
       ! Distribute particles evenly along 5 sides
       i=1
       do k=1,5
          x_start=verts(k,1)
          y_start=verts(k,2)
          if (k.lt.5) then
             x_end=verts(k+1,1)
             y_end=verts(k+1,2)
          else
             x_end=verts(1,1)
             y_end=verts(1,2)
          end if
          ! Total number of particles on this side (even split)
          total_side_particles=n/5
          if (k.le.mod(n,5)) total_side_particles = total_side_particles + 1
          do j=0,total_side_particles-1
             if (i.gt.n) exit
             frac=real(j,WP)/real(max(1,total_side_particles-1),WP)
             x_seg=(1.0_WP-frac)*x_start+frac*x_end
             y_seg=(1.0_WP-frac)*y_start+frac*y_end
             pos(1,i)=x0+x_seg
             pos(2,i)=y0+y_seg
             i=i+1
          end do
       end do
       ! Pad remainder at center if n < 5
       do while (i.le.n)
          pos(1,i) = x0
          pos(2,i) = y0
          i=i+1
       end do
    case default
       call die('[get_positions] unknown shape: '//trim(shape))
    end select
  end function get_positions

  
  !> Calculate the CFL
  subroutine get_cfl(this,dt,cflc,cfl)
    use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
    use parallel, only: MPI_REAL_WP
    implicit none
    class(multisphere), intent(inout) :: this
    real(WP), intent(in)  :: dt
    real(WP), intent(out) :: cflc
    real(WP), optional :: cfl
    integer :: i,j,ierr
    real(WP) :: my_CFLp_x,my_CFLp_y,my_CFLp_z,my_CFL_col

    ! Set the CFLs to zero
    my_CFLp_x=0.0_WP; my_CFLp_y=0.0_WP; my_CFLp_z=0.0_WP; my_CFL_col=0.0_WP
    do j=1,this%nb_
       do i=1,this%b(j)%np
          my_CFLp_x=max(my_CFLp_x,abs(this%b(j)%p_vel(1,i))*this%cfg%dxi(this%b(j)%p_ind(1,i)))
          my_CFLp_y=max(my_CFLp_y,abs(this%b(j)%p_vel(2,i))*this%cfg%dyi(this%b(j)%p_ind(2,i)))
          my_CFLp_z=max(my_CFLp_z,abs(this%b(j)%p_vel(3,i))*this%cfg%dzi(this%b(j)%p_ind(3,i)))
          my_CFL_col=max(my_CFL_col,sqrt(sum(this%b(j)%p_vel(:,i)**2))/this%b(j)%p_d(i))
       end do
    end do
    my_CFLp_x=my_CFLp_x*dt; my_CFLp_y=my_CFLp_y*dt; my_CFLp_z=my_CFLp_z*dt

    ! Get the parallel max
    call MPI_ALLREDUCE(my_CFLp_x,this%CFLp_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
    call MPI_ALLREDUCE(my_CFLp_y,this%CFLp_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
    call MPI_ALLREDUCE(my_CFLp_z,this%CFLp_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)

    ! Return the maximum convective CFL
    cflc=max(this%CFLp_x,this%CFLp_y,this%CFLp_z)

    ! Compute collision CFL
    my_CFL_col=10.0_WP*my_CFL_col*dt
    call MPI_ALLREDUCE(my_CFL_col,this%CFL_col,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)

    ! If asked for, also return the maximum overall CFL
    if (present(CFL)) cfl=max(cflc,this%CFL_col)
    
  end subroutine get_cfl


  !> Extract various monitoring data from particle field
  subroutine get_max(this)
    use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN,MPI_SUM
    use parallel, only: MPI_REAL_WP
    implicit none
    class(multisphere), intent(inout) :: this
    real(WP) :: buf,safe_nb
    integer :: i,j,k,ierr

    ! Create safe nb
    safe_nb=real(max(this%nb,1),WP)

    ! Diameter and velocity min/max/mean
    this%Umin=huge(1.0_WP); this%Umax=-huge(1.0_WP); this%Umean=0.0_WP
    this%Vmin=huge(1.0_WP); this%Vmax=-huge(1.0_WP); this%Vmean=0.0_WP
    this%Wmin=huge(1.0_WP); this%Wmax=-huge(1.0_WP); this%Wmean=0.0_WP
    do j=1,this%nb_
       this%Umin=min(this%Umin,this%b(j)%vel(1)); this%Umax=max(this%Umax,this%b(j)%vel(1)); this%Umean=this%Umean+this%b(j)%vel(1)
       this%Vmin=min(this%Vmin,this%b(j)%vel(2)); this%Vmax=max(this%Vmax,this%b(j)%vel(2)); this%Vmean=this%Vmean+this%b(j)%vel(2)
       this%Wmin=min(this%Wmin,this%b(j)%vel(3)); this%Wmax=max(this%Wmax,this%b(j)%vel(3)); this%Wmean=this%Wmean+this%b(j)%vel(3)
    end do
    call MPI_ALLREDUCE(this%Umin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Umin =buf
    call MPI_ALLREDUCE(this%Umax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Umax =buf
    call MPI_ALLREDUCE(this%Umean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Umean=buf/safe_nb
    call MPI_ALLREDUCE(this%Vmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Vmin =buf
    call MPI_ALLREDUCE(this%Vmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Vmax =buf
    call MPI_ALLREDUCE(this%Vmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Vmean=buf/safe_nb
    call MPI_ALLREDUCE(this%Wmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Wmin =buf
    call MPI_ALLREDUCE(this%Wmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Wmax =buf
    call MPI_ALLREDUCE(this%Wmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Wmean=buf/safe_nb

    ! Diameter and velocity variance
    this%Uvar=0.0_WP
    this%Vvar=0.0_WP
    this%Wvar=0.0_WP
    do j=1,this%nb_
       this%Uvar=this%Uvar+(this%b(j)%vel(1)-this%Umean)**2.0_WP
       this%Vvar=this%Vvar+(this%b(j)%vel(2)-this%Vmean)**2.0_WP
       this%Wvar=this%Wvar+(this%b(j)%vel(3)-this%Wmean)**2.0_WP
    end do
    call MPI_ALLREDUCE(this%Uvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Uvar=buf/safe_nb
    call MPI_ALLREDUCE(this%Vvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Vvar=buf/safe_nb
    call MPI_ALLREDUCE(this%Wvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Wvar=buf/safe_nb

    ! Get mean, max, and min volume fraction
    this%VFmean=0.0_WP
    this%VFmax =-huge(1.0_WP)
    this%VFmin =+huge(1.0_WP)
    do k=this%cfg%kmin_,this%cfg%kmax_
       do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
             this%VFmean=this%VFmean+this%cfg%VF(i,j,k)*this%cfg%vol(i,j,k)*this%VF(i,j,k)
             this%VFmax =max(this%VFmax,this%VF(i,j,k))
             this%VFmin =min(this%VFmin,this%VF(i,j,k))
          end do
       end do
    end do
    call MPI_ALLREDUCE(this%VFmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%VFmean=buf/this%cfg%fluid_vol
    call MPI_ALLREDUCE(this%VFmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%VFmax =buf
    call MPI_ALLREDUCE(this%VFmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%VFmin =buf

    ! Get volume fraction variance
    this%VFvar=0.0_WP
    do k=this%cfg%kmin_,this%cfg%kmax_
       do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
             this%VFvar=this%VFvar+this%cfg%VF(i,j,k)*this%cfg%vol(i,j,k)*(this%VF(i,j,k)-this%VFmean)**2
          end do
       end do
    end do
    call MPI_ALLREDUCE(this%VFvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%VFvar=buf/this%cfg%fluid_vol

  end subroutine get_max


  !> Update particle mesh using our current particles
  subroutine update_partmesh(this,pmesh)
    use partmesh_class, only: partmesh
    implicit none
    class(multisphere), intent(inout) :: this
    class(partmesh), intent(inout) :: pmesh
    integer :: i,j,k
    ! Reset particle mesh storage
    call pmesh%reset()
    if (this%nb_.gt.0) then
       ! Copy particle info
       call pmesh%set_size(this%np_)
       k=0
       do j=1,this%nb_
          do i=1,this%b(j)%np
             k=k+1
             pmesh%pos(:,k)=this%b(j)%p_pos(:,i)
          end do
       end do
    end if
    ! Root adds a particle if there are none
    if (this%nb.eq.0.and.this%cfg%amRoot) then
       call pmesh%set_size(1)
       pmesh%pos(1,1)=this%cfg%x(this%cfg%imin)
       pmesh%pos(2,1)=this%cfg%y(this%cfg%jmin)
       pmesh%pos(3,1)=this%cfg%z(this%cfg%kmin)
    end if
  end subroutine update_partmesh


  !> Creation of the MPI datatype for rigid body
  subroutine prepare_mpi_body()
    use mpi_f08
    use messager, only: die
    implicit none
    integer(MPI_ADDRESS_KIND), dimension(body_nblock) :: disp
    integer(MPI_ADDRESS_KIND) :: lb,extent
    type(MPI_Datatype) :: MPI_BODY_TMP
    integer :: i,mysize,ierr
    ! Prepare the displacement array
    disp(1)=0
    do i=2,body_nblock
       call MPI_Type_size(body_tblock(i-1),mysize,ierr)
       disp(i)=disp(i-1)+int(mysize,MPI_ADDRESS_KIND)*int(body_lblock(i-1),MPI_ADDRESS_KIND)
    end do
    ! Create and commit the new type
    call MPI_Type_create_struct(body_nblock,body_lblock,disp,body_tblock,MPI_BODY_TMP,ierr)
    call MPI_Type_get_extent(MPI_BODY_TMP,lb,extent,ierr)
    call MPI_Type_create_resized(MPI_BODY_TMP,lb,extent,MPI_BODY,ierr)
    call MPI_Type_commit(MPI_BODY,ierr)
    ! If a problem was encountered, say it
    if (ierr.ne.0) call die('[multisphere prepare_mpi_body] MPI Body type creation failed')
    ! Get the size of this type
    call MPI_type_size(MPI_BODY,MPI_BODY_SIZE,ierr)
  end subroutine prepare_mpi_body


  !> Synchronize rigid body arrays across processors
  subroutine sync(this)
    use mpi_f08
    implicit none
    class(multisphere), intent(inout) :: this
    integer, dimension(0:this%cfg%nproc-1) :: nsend_proc,nrecv_proc
    integer, dimension(0:this%cfg%nproc-1) :: nsend_disp,nrecv_disp
    integer :: i,j,n,prank,ierr
    type(body), dimension(:), allocatable :: buf_send
    ! Recycle first to minimize communication load
    call this%recycle()
    ! Prepare information about what to send
    nsend_proc=0
    do n=1,this%nb_
       prank=this%cfg%get_rank(this%b(n)%ind)
       nsend_proc(prank)=nsend_proc(prank)+1
    end do
    nsend_proc(this%cfg%rank)=0
    ! Inform processors of what they will receive
    call MPI_ALLtoALL(nsend_proc,1,MPI_INTEGER,nrecv_proc,1,MPI_INTEGER,this%cfg%comm,ierr)
    ! Prepare displacements for all-to-all
    nsend_disp(0)=0
    nrecv_disp(0)=this%nb_   !< Directly add bodies at the end of main array
    do n=1,this%cfg%nproc-1
       nsend_disp(n)=nsend_disp(n-1)+nsend_proc(n-1)
       nrecv_disp(n)=nrecv_disp(n-1)+nrecv_proc(n-1)
    end do
    ! Allocate buffer to send rigid bodies
    allocate(buf_send(sum(nsend_proc)))
    ! Pack the rigid bodies in the send buffer
    nsend_proc=0
    do n=1,this%nb_
       ! Get the rank
       prank=this%cfg%get_rank(this%b(n)%ind)
       ! Skip bodies still inside
       if (prank.eq.this%cfg%rank) cycle
       ! Pack up for sending
       nsend_proc(prank)=nsend_proc(prank)+1
       buf_send(nsend_disp(prank)+nsend_proc(prank))=this%b(n)
       ! Flag body for removal
       this%b(n)%flag=1
    end do
    ! Allocate buffer for receiving bodies
    call this%resize(this%nb_+sum(nrecv_proc))
    ! Perform communication
    call MPI_ALLtoALLv(buf_send,nsend_proc,nsend_disp,MPI_BODY,this%b,nrecv_proc,nrecv_disp,MPI_BODY,this%cfg%comm,ierr)
    ! Deallocate buffer
    deallocate(buf_send)
    ! Recycle to remove duplicate bodies
    call this%recycle()
    ! Integer values in part object may have been corrupted
    do j=1,this%nb_
       do i=1,this%b(j)%np
          this%b(j)%p_ind(:,i)=this%cfg%get_ijk_global(this%b(j)%p_pos(:,i),[this%cfg%imin,this%cfg%jmin,this%cfg%kmin])
       end do
    end do
  end subroutine sync


   !> Share rigid bodies across processor boundaries
   subroutine share(this,nover)
      use mpi_f08
      use messager, only: warn,die
      implicit none
      class(multisphere), intent(inout) :: this
      integer, optional :: nover
      type(body), dimension(:), allocatable :: tosend
      type(body), dimension(:), allocatable :: torecv
      integer :: i,j,n,no,nsend,nrecv
      type(MPI_Status) :: status
      integer :: isrc,idst,ierr
      
      ! Check overlap size
      if (present(nover)) then
         no=nover
         if (no.gt.this%cfg%no) then
            call warn('[multisphere_class share] Specified overlap is larger than that of cfg - reducing no')
            no=this%cfg%no
         else if (no.le.0) then
            call die('[multisphere_class share] Specified overlap cannot be less or equal to zero')
         end if
      else
         no=1
      end if
      
      ! Clean up ghost array
      call this%resize_ghost(n=0); this%ng_=0
      
      ! Share ghost bodies in -x (no ghosts are sent here)
      nsend=0
      do n=1,this%nb_
         if (this%b(n)%ind(1).lt.this%cfg%imin_+no) nsend=nsend+1
      end do
      allocate(tosend(nsend))
      nsend=0
      do n=1,this%nb_
         if (this%b(n)%ind(1).lt.this%cfg%imin_+no) then
            nsend=nsend+1
            tosend(nsend)=this%b(n)
            if (this%cfg%xper.and.tosend(nsend)%ind(1).lt.this%cfg%imin+no) then
               tosend(nsend)%pos(1)=tosend(nsend)%pos(1)+this%cfg%xL
               tosend(nsend)%ind(1)=tosend(nsend)%ind(1)+this%cfg%nx
            end if
         end if
      end do
      nrecv=0
      call MPI_CART_SHIFT(this%cfg%comm,0,-1,isrc,idst,ierr)
      call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
      allocate(torecv(nrecv))
      call MPI_SENDRECV(tosend,nsend,MPI_BODY,idst,0,torecv,nrecv,MPI_BODY,isrc,0,this%cfg%comm,status,ierr)
      call this%resize_ghost(this%ng_+nrecv)
      this%g(this%ng_+1:this%ng_+nrecv)=torecv
      this%ng_=this%ng_+nrecv
      if (allocated(tosend)) deallocate(tosend)
      if (allocated(torecv)) deallocate(torecv)
      
      ! Share ghost bodies in +x (no ghosts are sent here)
      nsend=0
      do n=1,this%nb_
         if (this%b(n)%ind(1).gt.this%cfg%imax_-no) nsend=nsend+1
      end do
      allocate(tosend(nsend))
      nsend=0
      do n=1,this%nb_
         if (this%b(n)%ind(1).gt.this%cfg%imax_-no) then
            nsend=nsend+1
            tosend(nsend)=this%b(n)
            if (this%cfg%xper.and.tosend(nsend)%ind(1).gt.this%cfg%imax-no) then
               tosend(nsend)%pos(1)=tosend(nsend)%pos(1)-this%cfg%xL
               tosend(nsend)%ind(1)=tosend(nsend)%ind(1)-this%cfg%nx
            end if
         end if
      end do
      nrecv=0
      call MPI_CART_SHIFT(this%cfg%comm,0,+1,isrc,idst,ierr)
      call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
      allocate(torecv(nrecv))
      call MPI_SENDRECV(tosend,nsend,MPI_BODY,idst,0,torecv,nrecv,MPI_BODY,isrc,0,this%cfg%comm,status,ierr)
      call this%resize_ghost(this%ng_+nrecv)
      this%g(this%ng_+1:this%ng_+nrecv)=torecv
      this%ng_=this%ng_+nrecv
      if (allocated(tosend)) deallocate(tosend)
      if (allocated(torecv)) deallocate(torecv)
      
      ! Share ghost bodies in -y (ghosts need to be sent now)
      nsend=0
      do n=1,this%nb_
         if (this%b(n)%ind(2).lt.this%cfg%jmin_+no) nsend=nsend+1
      end do
      do n=1,this%ng_
         if (this%g(n)%ind(2).lt.this%cfg%jmin_+no) nsend=nsend+1
      end do
      allocate(tosend(nsend))
      nsend=0
      do n=1,this%nb_
         if (this%b(n)%ind(2).lt.this%cfg%jmin_+no) then
            nsend=nsend+1
            tosend(nsend)=this%b(n)
            if (this%cfg%yper.and.tosend(nsend)%ind(2).lt.this%cfg%jmin+no) then
               tosend(nsend)%pos(2)=tosend(nsend)%pos(2)+this%cfg%yL
               tosend(nsend)%ind(2)=tosend(nsend)%ind(2)+this%cfg%ny
            end if
         end if
      end do
      do n=1,this%ng_
         if (this%g(n)%ind(2).lt.this%cfg%jmin_+no) then
            nsend=nsend+1
            tosend(nsend)=this%g(n)
            if (this%cfg%yper.and.tosend(nsend)%ind(2).lt.this%cfg%jmin+no) then
               tosend(nsend)%pos(2)=tosend(nsend)%pos(2)+this%cfg%yL
               tosend(nsend)%ind(2)=tosend(nsend)%ind(2)+this%cfg%ny
            end if
         end if
      end do
      nrecv=0
      call MPI_CART_SHIFT(this%cfg%comm,1,-1,isrc,idst,ierr)
      call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
      allocate(torecv(nrecv))
      call MPI_SENDRECV(tosend,nsend,MPI_BODY,idst,0,torecv,nrecv,MPI_BODY,isrc,0,this%cfg%comm,status,ierr)
      call this%resize_ghost(this%ng_+nrecv)
      this%g(this%ng_+1:this%ng_+nrecv)=torecv
      this%ng_=this%ng_+nrecv
      if (allocated(tosend)) deallocate(tosend)
      if (allocated(torecv)) deallocate(torecv)
      
      ! Share ghost bodies in +y (ghosts need to be sent now - but not newly received ghosts!)
      nsend=0
      do n=1,this%nb_
         if (this%b(n)%ind(2).gt.this%cfg%jmax_-no) nsend=nsend+1
      end do
      do n=1,this%ng_-nrecv
         if (this%g(n)%ind(2).gt.this%cfg%jmax_-no) nsend=nsend+1
      end do
      allocate(tosend(nsend))
      nsend=0
      do n=1,this%nb_
         if (this%b(n)%ind(2).gt.this%cfg%jmax_-no) then
            nsend=nsend+1
            tosend(nsend)=this%b(n)
            if (this%cfg%yper.and.tosend(nsend)%ind(2).gt.this%cfg%jmax-no) then
               tosend(nsend)%pos(2)=tosend(nsend)%pos(2)-this%cfg%yL
               tosend(nsend)%ind(2)=tosend(nsend)%ind(2)-this%cfg%ny
            end if
         end if
      end do
      do n=1,this%ng_-nrecv
         if (this%g(n)%ind(2).gt.this%cfg%jmax_-no) then
            nsend=nsend+1
            tosend(nsend)=this%g(n)
            if (this%cfg%yper.and.tosend(nsend)%ind(2).gt.this%cfg%jmax-no) then
               tosend(nsend)%pos(2)=tosend(nsend)%pos(2)-this%cfg%yL
               tosend(nsend)%ind(2)=tosend(nsend)%ind(2)-this%cfg%ny
            end if
         end if
      end do
      nrecv=0
      call MPI_CART_SHIFT(this%cfg%comm,1,+1,isrc,idst,ierr)
      call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
      allocate(torecv(nrecv))
      call MPI_SENDRECV(tosend,nsend,MPI_BODY,idst,0,torecv,nrecv,MPI_BODY,isrc,0,this%cfg%comm,status,ierr)
      call this%resize_ghost(this%ng_+nrecv)
      this%g(this%ng_+1:this%ng_+nrecv)=torecv
      this%ng_=this%ng_+nrecv
      if (allocated(tosend)) deallocate(tosend)
      if (allocated(torecv)) deallocate(torecv)
      
      ! Share ghost bodies in -z (ghosts need to be sent now)
      if (this%cfg%nz.gt.1) then
         nsend=0
         do n=1,this%nb_
            if (this%b(n)%ind(3).lt.this%cfg%kmin_+no) nsend=nsend+1
         end do
         do n=1,this%ng_
            if (this%g(n)%ind(3).lt.this%cfg%kmin_+no) nsend=nsend+1
         end do
         allocate(tosend(nsend))
         nsend=0
         do n=1,this%nb_
            if (this%b(n)%ind(3).lt.this%cfg%kmin_+no) then
               nsend=nsend+1
               tosend(nsend)=this%b(n)
               if (this%cfg%zper.and.tosend(nsend)%ind(3).lt.this%cfg%kmin+no) then
                  tosend(nsend)%pos(3)=tosend(nsend)%pos(3)+this%cfg%zL
                  tosend(nsend)%ind(3)=tosend(nsend)%ind(3)+this%cfg%nz
               end if
            end if
         end do
         do n=1,this%ng_
            if (this%g(n)%ind(3).lt.this%cfg%kmin_+no) then
               nsend=nsend+1
               tosend(nsend)=this%g(n)
               if (this%cfg%zper.and.tosend(nsend)%ind(3).lt.this%cfg%kmin+no) then
                  tosend(nsend)%pos(3)=tosend(nsend)%pos(3)+this%cfg%zL
                  tosend(nsend)%ind(3)=tosend(nsend)%ind(3)+this%cfg%nz
               end if
            end if
         end do
         nrecv=0
         call MPI_CART_SHIFT(this%cfg%comm,2,-1,isrc,idst,ierr)
         call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
         allocate(torecv(nrecv))
         call MPI_SENDRECV(tosend,nsend,MPI_BODY,idst,0,torecv,nrecv,MPI_BODY,isrc,0,this%cfg%comm,status,ierr)
         call this%resize_ghost(this%ng_+nrecv)
         this%g(this%ng_+1:this%ng_+nrecv)=torecv
         this%ng_=this%ng_+nrecv
         if (allocated(tosend)) deallocate(tosend)
         if (allocated(torecv)) deallocate(torecv)

         ! Share ghost bodies in +z (ghosts need to be sent now - but not newly received ghosts!)
         nsend=0
         do n=1,this%nb_
            if (this%b(n)%ind(3).gt.this%cfg%kmax_-no) nsend=nsend+1
         end do
         do n=1,this%ng_-nrecv
            if (this%g(n)%ind(3).gt.this%cfg%kmax_-no) nsend=nsend+1
         end do
         allocate(tosend(nsend))
         nsend=0
         do n=1,this%nb_
            if (this%b(n)%ind(3).gt.this%cfg%kmax_-no) then
               nsend=nsend+1
               tosend(nsend)=this%b(n)
               if (this%cfg%zper.and.tosend(nsend)%ind(3).gt.this%cfg%kmax-no) then
                  tosend(nsend)%pos(3)=tosend(nsend)%pos(3)-this%cfg%zL
                  tosend(nsend)%ind(3)=tosend(nsend)%ind(3)-this%cfg%nz
               end if
            end if
         end do
         do n=1,this%ng_-nrecv
            if (this%g(n)%ind(3).gt.this%cfg%kmax_-no) then
               nsend=nsend+1
               tosend(nsend)=this%g(n)
               if (this%cfg%zper.and.tosend(nsend)%ind(3).gt.this%cfg%kmax-no) then
                  tosend(nsend)%pos(3)=tosend(nsend)%pos(3)-this%cfg%zL
                  tosend(nsend)%ind(3)=tosend(nsend)%ind(3)-this%cfg%nz
               end if
            end if
         end do
         nrecv=0
         call MPI_CART_SHIFT(this%cfg%comm,2,+1,isrc,idst,ierr)
         call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
         allocate(torecv(nrecv))
         call MPI_SENDRECV(tosend,nsend,MPI_BODY,idst,0,torecv,nrecv,MPI_BODY,isrc,0,this%cfg%comm,status,ierr)
         call this%resize_ghost(this%ng_+nrecv)
         this%g(this%ng_+1:this%ng_+nrecv)=torecv
         this%ng_=this%ng_+nrecv
         if (allocated(tosend)) deallocate(tosend)
         if (allocated(torecv)) deallocate(torecv)
      end if

      ! Integer values in part object may have been corrupted
      do j=1,this%nb_
         do i=1,this%b(j)%np
            this%b(j)%p_ind(:,i)=this%cfg%get_ijk_global(this%b(j)%p_pos(:,i),[this%cfg%imin,this%cfg%jmin,this%cfg%kmin])
         end do
      end do
      do j=1,this%ng_
         do i=1,this%g(j)%np
            this%g(j)%p_ind(:,i)=this%cfg%get_ijk_global(this%g(j)%p_pos(:,i),[this%cfg%imin,this%cfg%jmin,this%cfg%kmin])
         end do
      end do
      
   end subroutine share
   
   
   !> Adaptation of body array size
   subroutine resize(this,n)
      implicit none
      class(multisphere), intent(inout) :: this
      integer, intent(in) :: n
      type(body), dimension(:), allocatable :: tmp
      integer :: size_now,size_new
      ! Resize body array to size n
      if (.not.allocated(this%b)) then
         ! Allocate directly to size n
         allocate(this%b(n))
         this%b(1:n)%flag=1
      else
         ! Update from a non-zero size to another non-zero size
         size_now=size(this%b,dim=1)
         if (n.gt.size_now) then
            size_new=max(n,int(real(size_now,WP)*coeff_up))
            allocate(tmp(size_new))
            tmp(1:size_now)=this%b
            tmp(size_now+1:)%flag=1
            call move_alloc(tmp,this%b)
         else if (n.lt.int(real(size_now,WP)*coeff_dn)) then
            allocate(tmp(n))
            tmp(1:n)=this%b(1:n)
            call move_alloc(tmp,this%b)
         end if
      end if
   end subroutine resize


  !> Adaptation of ghost array size
  subroutine resize_ghost(this,n)
    implicit none
    class(multisphere), intent(inout) :: this
    integer, intent(in) :: n
    type(body), dimension(:), allocatable :: tmp
    integer :: size_now,size_new
    ! Resize ghost array to size n
    if (.not.allocated(this%g)) then
       ! Allocate directly to size n
       allocate(this%g(n))
       this%g(1:n)%flag=1
    else
       ! Update from a non-zero size to another non-zero size
       size_now=size(this%g,dim=1)
       if (n.gt.size_now) then
          size_new=max(n,int(real(size_now,WP)*coeff_up))
          allocate(tmp(size_new))
          tmp(1:size_now)=this%g
          tmp(size_now+1:)%flag=1
          call move_alloc(tmp,this%g)
       else if (n.lt.int(real(size_now,WP)*coeff_dn)) then
          allocate(tmp(n))
          tmp(1:n)=this%g(1:n)
          call move_alloc(tmp,this%g)
       end if
    end if
  end subroutine resize_ghost


  !> Clean-up of particle array by removing flag=1 particles
  subroutine recycle(this)
    use mpi_f08, only : MPI_INTEGER
    implicit none
    class(multisphere), intent(inout) :: this
    integer :: new_size,i,ierr
    ! Compact all active particles at the beginning of the array
    new_size=0
    if (allocated(this%b)) then
       do i=1,size(this%b,dim=1)
          if (this%b(i)%flag.ne.1) then
             new_size=new_size+1
             if (i.ne.new_size) then
                this%b(new_size)=this%b(i)
                this%b(i)%flag=1
             end if
          end if
       end do
    end if
    ! Resize to new size
    call this%resize(new_size)
    ! Update number of rigid bodies
    this%nb_=new_size
    call MPI_ALLGATHER(this%nb_,1,MPI_INTEGER,this%nb_proc,1,MPI_INTEGER,this%cfg%comm,ierr)
    this%nb=sum(this%nb_proc)
    ! Update number of particles
    this%np_=0
    do i=1,this%nb_
       this%np_=this%np_+this%b(i)%np
    end do
    call MPI_ALLGATHER(this%np_,1,MPI_INTEGER,this%np_proc,1,MPI_INTEGER,this%cfg%comm,ierr)
    this%np=sum(this%np_proc)
  end subroutine recycle


  !> Parallel write rigid bodies to file
  subroutine write(this,filename)
    use mpi_f08
    use messager, only: die
    use parallel, only: info_mpiio
    implicit none
    class(multisphere), intent(inout) :: this
    character(len=*), intent(in) :: filename
    type(MPI_File) :: ifile
    type(MPI_Status):: status
    integer(kind=MPI_OFFSET_KIND) :: offset
    integer :: i,ierr,iunit

    ! Root serial-writes the file header
    if (this%cfg%amRoot) then
       ! Open the file
       open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
       if (ierr.ne.0) call die('[multisphere write] Problem encountered while serial-opening data file: '//trim(filename))
       ! Number of bodies and body object size
       write(iunit) this%nb,MPI_BODY_SIZE
       ! Done with the header
       close(iunit)
    end if

    ! The rest is done in parallel
    call MPI_FILE_OPEN(this%cfg%comm,trim(filename),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
    if (ierr.ne.0) call die('[multisphere write] Problem encountered while parallel-opening data file: '//trim(filename))

    ! Get current position
    call MPI_FILE_GET_POSITION(ifile,offset,ierr)

    ! Compute the offset and write
    do i=1,this%cfg%rank
       offset=offset+int(this%nb_proc(i),MPI_OFFSET_KIND)*int(MPI_BODY_SIZE,MPI_OFFSET_KIND)
    end do
    if (this%nb_.gt.0) call MPI_FILE_WRITE_AT(ifile,offset,this%b,this%nb_,MPI_BODY,status,ierr)

    ! Close the file
    call MPI_FILE_CLOSE(ifile,ierr)

    ! Log/screen output
    logging: block
      use, intrinsic :: iso_fortran_env, only: output_unit
      use param,    only: verbose
      use messager, only: log
      use string,   only: str_long
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Wrote ",i0," bodies to file [",a,"] on partitioned grid [",a,"]")') this%nb,trim(filename),trim(this%cfg%name)
         if (verbose.gt.2) write(output_unit,'(a)') trim(message)
         if (verbose.gt.1) call log(message)
      end if
    end block logging

  end subroutine write


  !> Parallel read rigid bodies to file
  subroutine read(this,filename)
    use mpi_f08
    use messager, only: die
    use parallel, only: info_mpiio
    implicit none
    class(multisphere), intent(inout) :: this
    character(len=*), intent(in) :: filename
    type(MPI_File) :: ifile
    type(MPI_Status):: status
    integer(kind=MPI_OFFSET_KIND) :: offset,header_offset
    integer :: i,j,k,ierr,nbadd,psize,nchunk,cnt
    integer, dimension(:,:), allocatable :: ppp

    ! First open the file in parallel
    call MPI_FILE_OPEN(this%cfg%comm,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
    if (ierr.ne.0) call die('[multisphere read] Problem encountered while reading data file: '//trim(filename))

    ! Read file header first
    call MPI_FILE_READ_ALL(ifile,nbadd,1,MPI_INTEGER,status,ierr)
    call MPI_FILE_READ_ALL(ifile,psize,1,MPI_INTEGER,status,ierr)

    ! Remember current position
    call MPI_FILE_GET_POSITION(ifile,header_offset,ierr)

    ! Check compatibility of body type
    if (psize.ne.MPI_BODY_SIZE) call die('[multisphere read] Body type unreadable')

    ! Naively share reading task among all processors
    nchunk=int(nbadd/(this%cfg%nproc*part_chunk_size))+1
    allocate(ppp(this%cfg%nproc,nchunk))
    ppp=int(nbadd/(this%cfg%nproc*nchunk))
    cnt=0
    out:do j=1,nchunk
       do i=1,this%cfg%nproc
          cnt=cnt+1
          if (cnt.gt.mod(nbadd,this%cfg%nproc*nchunk)) exit out
          ppp(i,j)=ppp(i,j)+1
       end do
    end do out

    ! Read by chunk
    do j=1,nchunk
       ! Find offset
       offset=header_offset+int(MPI_BODY_SIZE,MPI_OFFSET_KIND)*int(sum(ppp(1:this%cfg%rank,:))+sum(ppp(this%cfg%rank+1,1:j-1)),MPI_OFFSET_KIND)
       ! Resize particle array
       call this%resize(this%nb_+ppp(this%cfg%rank+1,j))
       ! Read this file
       call MPI_FILE_READ_AT(ifile,offset,this%b(this%nb_+1:this%nb_+ppp(this%cfg%rank+1,j)),ppp(this%cfg%rank+1,j),MPI_BODY,status,ierr)
       ! Most general case: relocate every droplet
       do i=this%nb_+1,this%nb_+ppp(this%cfg%rank+1,j)
          this%b(i)%ind=this%cfg%get_ijk_global(this%b(i)%pos,this%b(i)%ind)
          do k=1,this%b(i)%np
             this%b(i)%p_ind(:,k)=this%cfg%get_ijk_global(this%b(i)%p_pos(:,k),this%b(i)%p_ind(:,k))
          end do
       end do
       ! Exchange all that
       call this%sync()
    end do

    ! Close the file
    call MPI_FILE_CLOSE(ifile,ierr)

    ! Log/screen output
    logging: block
      use, intrinsic :: iso_fortran_env, only: output_unit
      use param,    only: verbose
      use messager, only: log
      use string,   only: str_long
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Read ",i0," bodies from file [",a,"] on partitioned grid [",a,"]")') nbadd,trim(filename),trim(this%cfg%name)
         if (verbose.gt.2) write(output_unit,'(a)') trim(message)
         if (verbose.gt.1) call log(message)
      end if
    end block logging

  end subroutine read


end module multisphere_class
