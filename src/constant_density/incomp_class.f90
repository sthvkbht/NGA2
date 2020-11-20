!> Incompressible flow solver class:
!> Provides support for various BC, RHS calculation,
!> implicit solver, and pressure solution
!> Assumes constant viscosity and density.
module incomp_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use ils_class,      only: ils
   use iterator_class, only: iterator
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: incomp,bcond
   
   ! List of known available bcond for this solver
   integer, parameter, public :: wall=1              !< Dirichlet at zero condition
   integer, parameter, public :: dirichlet=2         !< Dirichlet condition
   integer, parameter, public :: neumann=3           !< Zero normal gradient
   integer, parameter, public :: convective=4        !< Convective outflow condition
   integer, parameter, public :: clipped_neumann=5   !< Clipped Neumann condition (outflow only)
   
   !> Boundary conditions for the incompressible solver
   type :: bcond
      type(bcond), pointer :: next                        !< Linked list of bconds
      character(len=str_medium) :: name='UNNAMED_BCOND'   !< Bcond name (default=UNNAMED_BCOND)
      integer :: type                                     !< Bcond type
      integer :: dir                                      !< Bcond direction (1 to 6)
      integer :: si,sj,sk                                 !< Index shift in the outward normal direction
      logical :: canCorrect                               !< Can this bcond be corrected for global conservation?
      logical :: inBalance                                !< Is this bcond contributing to global conservation?
      type(iterator) :: itr                               !< This is the iterator for the bcond
   end type bcond
   
   !> Bcond shift value
   integer, dimension(3,6), parameter :: shift=reshape([+1,0,0,-1,0,0,0,+1,0,0,-1,0,0,0,+1,0,0,-1],shape(shift))
   
   !> Incompressible solver object definition
   type :: incomp
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_INCOMP'  !< Solver name (default=UNNAMED_INCOMP)
      
      ! Constant property fluid
      real(WP) :: rho                                     !< This is our constant fluid density
      real(WP) :: visc                                    !< These is our constant fluid dynamic viscosity
      
      ! Boundary condition list
      integer :: nbc                                      !< Number of bcond for our solver
      real(WP), dimension(:), allocatable :: mfr          !< MFR through each bcond
      real(WP), dimension(:), allocatable :: area         !< Area for each bcond
      real(WP) :: correctable_area                        !< Area of bcond that can be corrected
      type(bcond), pointer :: first_bc                    !< List of bcond for our solver
      
      ! Flow variables
      real(WP), dimension(:,:,:), allocatable :: U        !< U velocity array
      real(WP), dimension(:,:,:), allocatable :: V        !< V velocity array
      real(WP), dimension(:,:,:), allocatable :: W        !< W velocity array
      real(WP), dimension(:,:,:), allocatable :: P        !< Pressure array
      
      ! Old flow variables
      real(WP), dimension(:,:,:), allocatable :: Uold     !< Uold velocity array
      real(WP), dimension(:,:,:), allocatable :: Vold     !< Vold velocity array
      real(WP), dimension(:,:,:), allocatable :: Wold     !< Wold velocity array
      
      ! Flow divergence
      real(WP), dimension(:,:,:), allocatable :: div      !< Divergence array
      
      ! Pressure solver
      type(ils) :: psolv                                  !< Iterative linear solver object for the pressure Poisson equation
      
      ! Implicit velocity solver
      type(ils) :: implicit                               !< Iterative linear solver object for an implicit prediction of the NS residual
      
      ! Metrics
      real(WP), dimension(:,:,:,:), allocatable :: itpu_x,itpu_y,itpu_z   !< Interpolation for U
      real(WP), dimension(:,:,:,:), allocatable :: itpv_x,itpv_y,itpv_z   !< Interpolation for V
      real(WP), dimension(:,:,:,:), allocatable :: itpw_x,itpw_y,itpw_z   !< Interpolation for W
      real(WP), dimension(:,:,:,:), allocatable :: divp_x,divp_y,divp_z   !< Divergence for P-cell
      real(WP), dimension(:,:,:,:), allocatable :: divu_x,divu_y,divu_z   !< Divergence for U-cell
      real(WP), dimension(:,:,:,:), allocatable :: divv_x,divv_y,divv_z   !< Divergence for V-cell
      real(WP), dimension(:,:,:,:), allocatable :: divw_x,divw_y,divw_z   !< Divergence for W-cell
      real(WP), dimension(:,:,:,:), allocatable :: grdu_x,grdu_y,grdu_z   !< Velocity gradient for U
      real(WP), dimension(:,:,:,:), allocatable :: grdv_x,grdv_y,grdv_z   !< Velocity gradient for V
      real(WP), dimension(:,:,:,:), allocatable :: grdw_x,grdw_y,grdw_z   !< Velocity gradient for W
      
      ! Masking info for metric modification
      integer, dimension(:,:,:), allocatable :: umask                     !< Integer array used for modifying U metrics
      integer, dimension(:,:,:), allocatable :: vmask                     !< Integer array used for modifying V metrics
      integer, dimension(:,:,:), allocatable :: wmask                     !< Integer array used for modifying W metrics
      
      ! CFL numbers
      real(WP) :: CFLc_x,CFLc_y,CFLc_z                                    !< Convective CFL numbers
      real(WP) :: CFLv_x,CFLv_y,CFLv_z                                    !< Viscous CFL numbers
      
      ! Monitoring quantities
      real(WP) :: Umax,Vmax,Wmax,Pmax,divmax                              !< Maximum velocity, pressure, divergence
      
   contains
      procedure :: print=>incomp_print                    !< Output solver to the screen
      procedure :: setup                                  !< Finish configuring the flow solver
      procedure :: add_bcond                              !< Add a boundary condition
      procedure :: get_bcond                              !< Get a boundary condition
      procedure :: apply_bcond                            !< Apply all boundary conditions
      procedure :: init_metrics                           !< Initialize metrics
      procedure :: adjust_metrics                         !< Adjust metrics
      procedure :: get_dmomdt                             !< Calculate dmom/dt
      procedure :: get_div                                !< Calculate velocity divergence
      procedure :: get_pgrad                              !< Calculate pressure gradient
      procedure :: get_cfl                                !< Calculate maximum CFL
      procedure :: get_max                                !< Calculate maximum field values
      procedure :: interp_vel                             !< Calculate interpolated velocity
      procedure :: get_mfr                                !< Calculate outgoing MFR through each bcond
      procedure :: correct_mfr                            !< Correct for mfr mismatch to ensure global conservation
      procedure :: shift_p                                !< Shift pressure to have zero average
      procedure :: solve_implicit                         !< Solve for the velocity residuals implicitly
   end type incomp
   
   
   !> Declare incompressible solver constructor
   interface incomp
      procedure constructor
   end interface incomp
   
contains
   
   
   !> Default constructor for incompressible flow solver
   function constructor(cfg,name) result(self)
      implicit none
      type(incomp) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Nullify bcond list
      self%nbc=0
      self%first_bc=>NULL()
      
      ! Allocate flow variables
      allocate(self%U(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%U=0.0_WP
      allocate(self%V(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%V=0.0_WP
      allocate(self%W(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%W=0.0_WP
      allocate(self%P(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%P=0.0_WP
      
      ! Allocate flow divergence
      allocate(self%div(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%div=0.0_WP
      
      ! Allocate old flow variables
      allocate(self%Uold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Uold=0.0_WP
      allocate(self%Vold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Vold=0.0_WP
      allocate(self%Wold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Wold=0.0_WP
      
      ! Create pressure solver object
      self%psolv=ils(cfg=self%cfg,name='Pressure Poisson Solver')
      
      ! Create implicit velocity solver object
      self%implicit=ils(cfg=self%cfg,name='Implicit NS residual')
      
      ! Prepare default metrics
      call self%init_metrics()
      
      ! Prepare face mask for U
      allocate(self%umask(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%umask=0
      do k=self%cfg%kmino_  ,self%cfg%kmaxo_
         do j=self%cfg%jmino_  ,self%cfg%jmaxo_
            do i=self%cfg%imino_+1,self%cfg%imaxo_
               if (minval(self%cfg%VF(i-1:i,j,k)).eq.0.0_WP) self%umask(i,j,k)=1
            end do
         end do
      end do
      call self%cfg%sync(self%umask)
      if (.not.self%cfg%xper.and.self%cfg%iproc.eq.1) self%umask(self%cfg%imino,:,:)=self%umask(self%cfg%imino+1,:,:)
      
      ! Prepare face mask for V
      allocate(self%vmask(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%vmask=0
      do k=self%cfg%kmino_  ,self%cfg%kmaxo_
         do j=self%cfg%jmino_+1,self%cfg%jmaxo_
            do i=self%cfg%imino_  ,self%cfg%imaxo_
               if (minval(self%cfg%VF(i,j-1:j,k)).eq.0.0_WP) self%vmask(i,j,k)=1
            end do
         end do
      end do
      call self%cfg%sync(self%vmask)
      if (.not.self%cfg%yper.and.self%cfg%jproc.eq.1) self%vmask(:,self%cfg%jmino,:)=self%vmask(:,self%cfg%jmino+1,:)
      
      ! Prepare face mask for W
      allocate(self%wmask(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%wmask=0
      do k=self%cfg%kmino_+1,self%cfg%kmaxo_
         do j=self%cfg%jmino_  ,self%cfg%jmaxo_
            do i=self%cfg%imino_  ,self%cfg%imaxo_
               if (minval(self%cfg%VF(i,j,k-1:k)).eq.0.0_WP) self%wmask(i,j,k)=1
            end do
         end do
      end do
      call self%cfg%sync(self%wmask)
      if (.not.self%cfg%zper.and.self%cfg%kproc.eq.1) self%wmask(:,:,self%cfg%kmino)=self%wmask(:,:,self%cfg%kmino+1)
      
      ! Adjust face masks to reflect domain boundaries
      if (.not.self%cfg%xper) then
         if (self%cfg%iproc.eq.           1) self%umask(self%cfg%imin  ,:,:)=2
         if (self%cfg%iproc.eq.self%cfg%npx) self%umask(self%cfg%imax+1,:,:)=2
      end if
      if (.not.self%cfg%yper) then
         if (self%cfg%jproc.eq.           1) self%vmask(:,self%cfg%jmin  ,:)=2
         if (self%cfg%jproc.eq.self%cfg%npy) self%vmask(:,self%cfg%jmax+1,:)=2
      end if
      if (.not.self%cfg%zper) then
         if (self%cfg%kproc.eq.           1) self%wmask(:,:,self%cfg%kmin  )=2
         if (self%cfg%kproc.eq.self%cfg%npz) self%wmask(:,:,self%cfg%kmax+1)=2
      end if
      
   end function constructor
      
   
   !> Metric initialization with no awareness of walls nor bcond
   subroutine init_metrics(this)
      implicit none
      class(incomp), intent(inout) :: this
      integer :: i,j,k
      
      ! Allocate finite difference velocity interpolation coefficients
      allocate(this%itpu_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%itpv_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%itpw_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%itpv_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (xy)
      allocate(this%itpw_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (zx)
      allocate(this%itpu_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (xy)
      allocate(this%itpw_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (yz)
      allocate(this%itpu_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Edge-centered (zx)
      allocate(this%itpv_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Edge-centered (yz)
      ! Create velocity interpolation coefficients to cell center [xm,ym,zm]
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               this%itpu_x(:,i,j,k)=[+0.5_WP,+0.5_WP] !< Linear interpolation in x of U from [x ,ym,zm]
               this%itpv_y(:,i,j,k)=[+0.5_WP,+0.5_WP] !< Linear interpolation in y of V from [xm,y ,zm]
               this%itpw_z(:,i,j,k)=[+0.5_WP,+0.5_WP] !< Linear interpolation in z of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create velocity interpolation coefficients to cell edge in x
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%itpv_x(:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x of V from [xm,y ,zm]
               this%itpw_x(:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create velocity interpolation coefficients to cell edge in y
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%itpu_y(:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y of U from [x ,ym,zm]
               this%itpw_y(:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create velocity interpolation coefficients to cell edge in z
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%itpu_z(:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z of U from [x ,ym,zm]
               this%itpv_z(:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z of V from [xm,y ,zm]
            end do
         end do
      end do
      
      ! Allocate finite volume divergence operators
      allocate(this%divp_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%divp_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%divp_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%divu_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (x)
      allocate(this%divu_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (x)
      allocate(this%divu_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (x)
      allocate(this%divv_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (y)
      allocate(this%divv_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (y)
      allocate(this%divv_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (y)
      allocate(this%divw_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (z)
      allocate(this%divw_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (z)
      allocate(this%divw_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Face-centered (z)
      ! Create divergence operator to cell center [xm,ym,zm] or tangent to cell face
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               this%divp_x(:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,zm]
               this%divp_y(:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,zm]
               this%divp_z(:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,z ]
               
               this%divu_y(:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,y ,zm]
               this%divu_z(:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,z ]
               
               this%divv_x(:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,y ,zm]
               this%divv_z(:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,z ]
               
               this%divw_x(:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,z ]
               this%divw_y(:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,z ]
            end do
         end do
      end do
      ! Create divergence operator perpendicular to cell face [x ,ym,zm]
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%divu_x(:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,zm]
            end do
         end do
      end do
      ! Create divergence operator perpendicular to cell face [xm,y ,zm]
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%divv_y(:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,zm]
            end do
         end do
      end do
      ! Create divergence operator perpendicular to cell face [xm,ym,z ]
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%divw_z(:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,zm]
            end do
         end do
      end do
      
      ! Allocate finite difference velocity gradient operators
      allocate(this%grdu_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%grdv_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%grdw_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%grdv_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (xy)
      allocate(this%grdw_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (zx)
      allocate(this%grdu_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (xy)
      allocate(this%grdw_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (yz)
      allocate(this%grdu_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Edge-centered (zx)
      allocate(this%grdv_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Edge-centered (yz)
      ! Create gradient coefficients to cell center [xm,ym,zm]
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               this%grdu_x(:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FD gradient in x of U from [x ,ym,zm]
               this%grdv_y(:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FD gradient in y of V from [xm,y ,zm]
               this%grdw_z(:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FD gradient in z of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create gradient coefficients to cell edge in x
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%grdv_x(:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient in x of V from [xm,y ,zm]
               this%grdw_x(:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient in x of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create gradient coefficients to cell edge in y
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%grdu_y(:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FD gradient in y of U from [x ,ym,zm]
               this%grdw_y(:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FD gradient in y of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create gradient coefficients to cell edge in z
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%grdu_z(:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient in z of U from [x ,ym,zm]
               this%grdv_z(:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient in z of V from [xm,y ,zm]
            end do
         end do
      end do
      
   end subroutine init_metrics
   
   
   !> Metric adjustment accounting for bconds and walls
   subroutine adjust_metrics(this)
      implicit none
      class(incomp), intent(inout) :: this
      integer :: i,j,k
      real(WP) :: delta
      
      ! Sync up u/v/wmasks
      call this%cfg%sync(this%umask)
      call this%cfg%sync(this%vmask)
      call this%cfg%sync(this%wmask)
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) this%umask(this%cfg%imino,:,:)=this%umask(this%cfg%imino+1,:,:)
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) this%vmask(:,this%cfg%jmino,:)=this%vmask(:,this%cfg%jmino+1,:)
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) this%wmask(:,:,this%cfg%kmino)=this%wmask(:,:,this%cfg%kmino+1)
      
      ! Loop over the domain and adjust divergence for P cell
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%cfg%VF(i,j,k).eq.0.0_WP) then
                  this%divp_x(:,i,j,k)=0.0_WP
                  this%divp_y(:,i,j,k)=0.0_WP
                  this%divp_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Loop over the domain and apply masked conditions to U metrics
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               if (this%umask(i,j,k).gt.0) then
                  this%divu_x(:,i,j,k)=0.0_WP
                  this%divu_y(:,i,j,k)=0.0_WP
                  this%divu_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Loop over the domain and apply masked conditions to V metrics
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               if (this%vmask(i,j,k).gt.0) then
                  this%divv_x(:,i,j,k)=0.0_WP
                  this%divv_y(:,i,j,k)=0.0_WP
                  this%divv_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Loop over the domain and apply masked conditions to W metrics
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               if (this%wmask(i,j,k).gt.0) then
                  this%divw_x(:,i,j,k)=0.0_WP
                  this%divw_y(:,i,j,k)=0.0_WP
                  this%divw_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Adjust gradient coefficients to cell edge in x
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               ! FD gradient in x of V from [xm,y ,zm]
               if (maxval(this%vmask(i-1:i,j,k)).gt.0) then
                  delta=0.0_WP
                  if (this%vmask(i  ,j,k).eq.0) delta=delta+(this%cfg%xm(i)-this%cfg%x (i  ))
                  if (this%vmask(i-1,j,k).eq.0) delta=delta+(this%cfg%x (i)-this%cfg%xm(i-1))
                  if (delta.gt.0.0_WP) then
                     this%grdv_x(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdv_x(:,i,j,k)=0.0_WP
                  end if
               end if
               ! FD gradient in x of W from [xm,ym,z ]
               if (maxval(this%wmask(i-1:i,j,k)).gt.0) then
                  delta=0.0_WP
                  if (this%wmask(i  ,j,k).eq.0) delta=delta+(this%cfg%xm(i)-this%cfg%x (i  ))
                  if (this%wmask(i-1,j,k).eq.0) delta=delta+(this%cfg%x (i)-this%cfg%xm(i-1))
                  if (delta.gt.0.0_WP) then
                     this%grdw_x(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdw_x(:,i,j,k)=0.0_WP
                  end if
               end if
            end do
         end do
      end do
      
      ! Adjust gradient coefficients to cell edge in y
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               ! FD gradient in y of U from [x ,ym,zm]
               if (maxval(this%umask(i,j-1:j,k)).gt.0) then
                  delta=0.0_WP
                  if (this%umask(i,j  ,k).eq.0) delta=delta+(this%cfg%ym(j)-this%cfg%y (j  ))
                  if (this%umask(i,j-1,k).eq.0) delta=delta+(this%cfg%y (j)-this%cfg%ym(j-1))
                  if (delta.gt.0.0_WP) then
                     this%grdu_y(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdu_y(:,i,j,k)=0.0_WP
                  end if
               end if
               ! FD gradient in y of W from [xm,ym,z ]
               if (maxval(this%wmask(i,j-1:j,k)).gt.0) then
                  delta=0.0_WP
                  if (this%wmask(i,j  ,k).eq.0) delta=delta+(this%cfg%ym(j)-this%cfg%y (j  ))
                  if (this%wmask(i,j-1,k).eq.0) delta=delta+(this%cfg%y (j)-this%cfg%ym(j-1))
                  if (delta.gt.0.0_WP) then
                     this%grdw_y(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdw_y(:,i,j,k)=0.0_WP
                  end if
               end if
            end do
         end do
      end do
      
      ! Adjust gradient coefficients to cell edge in z
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               ! FD gradient in z of U from [x ,ym,zm]
               if (maxval(this%umask(i,j,k-1:k)).gt.0) then
                  delta=0.0_WP
                  if (this%umask(i,j,k  ).eq.0) delta=delta+(this%cfg%zm(k)-this%cfg%z (k  ))
                  if (this%umask(i,j,k-1).eq.0) delta=delta+(this%cfg%z (k)-this%cfg%zm(k-1))
                  if (delta.gt.0.0_WP) then
                     this%grdu_z(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdu_z(:,i,j,k)=0.0_WP
                  end if
               end if
               ! FD gradient in z of V from [xm,y ,zm]
               if (maxval(this%vmask(i,j,k-1:k)).gt.0) then
                  delta=0.0_WP
                  if (this%vmask(i,j,k  ).eq.0) delta=delta+(this%cfg%zm(k)-this%cfg%z (k  ))
                  if (this%vmask(i,j,k-1).eq.0) delta=delta+(this%cfg%z (k)-this%cfg%zm(k-1))
                  if (delta.gt.0.0_WP) then
                     this%grdv_z(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdv_z(:,i,j,k)=0.0_WP
                  end if
               end if
            end do
         end do
      end do
      
      ! Adjust interpolation coefficients to cell edge in x
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               ! Linear interpolation in x of V from [xm,y ,zm]
               if (this%vmask(i,j,k).eq.0.and.this%vmask(i-1,j,k).gt.0) this%itpv_x(:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%vmask(i,j,k).gt.0.and.this%vmask(i-1,j,k).eq.0) this%itpv_x(:,i,j,k)=[0.0_WP,1.0_WP]
               ! Linear interpolation in x of W from [xm,ym,z ]
               if (this%wmask(i,j,k).eq.0.and.this%wmask(i-1,j,k).gt.0) this%itpw_x(:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%wmask(i,j,k).gt.0.and.this%wmask(i-1,j,k).eq.0) this%itpw_x(:,i,j,k)=[0.0_WP,1.0_WP]
            end do
         end do
      end do
      
      ! Adjust interpolation coefficients to cell edge in y
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               ! Linear interpolation in y of U from [x ,ym,zm]
               if (this%umask(i,j,k).eq.0.and.this%umask(i,j-1,k).gt.0) this%itpu_y(:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%umask(i,j,k).gt.0.and.this%umask(i,j-1,k).eq.0) this%itpu_y(:,i,j,k)=[0.0_WP,1.0_WP]
               ! Linear interpolation in y of W from [xm,ym,z ]
               if (this%wmask(i,j,k).eq.0.and.this%wmask(i,j-1,k).gt.0) this%itpw_y(:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%wmask(i,j,k).gt.0.and.this%wmask(i,j-1,k).eq.0) this%itpw_y(:,i,j,k)=[0.0_WP,1.0_WP]
            end do
         end do
      end do
      
      ! Adjust interpolation coefficients to cell edge in z
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               ! Linear interpolation in z of U from [x ,ym,zm]
               if (this%umask(i,j,k).eq.0.and.this%umask(i,j,k-1).gt.0) this%itpu_z(:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%umask(i,j,k).gt.0.and.this%umask(i,j,k-1).eq.0) this%itpu_z(:,i,j,k)=[0.0_WP,1.0_WP]
               !  Linear interpolation in z of V from [xm,y ,zm]
               if (this%vmask(i,j,k).eq.0.and.this%vmask(i,j,k-1).gt.0) this%itpv_z(:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%vmask(i,j,k).gt.0.and.this%vmask(i,j,k-1).eq.0) this%itpv_z(:,i,j,k)=[0.0_WP,1.0_WP]
            end do
         end do
      end do
      
   end subroutine adjust_metrics
   
   
   !> Finish setting up the flow solver now that bconds have been defined
   subroutine setup(this,pressure_ils,implicit_ils)
      implicit none
      class(incomp), intent(inout) :: this
      integer, intent(in) :: pressure_ils
      integer, intent(in) :: implicit_ils
      integer :: i,j,k
      
      ! Adjust metrics based on bcflag array
      call this%adjust_metrics()
      
      ! Set 7-pt stencil map for the pressure solver
      this%psolv%stc(1,:)=[ 0, 0, 0]
      this%psolv%stc(2,:)=[+1, 0, 0]
      this%psolv%stc(3,:)=[-1, 0, 0]
      this%psolv%stc(4,:)=[ 0,+1, 0]
      this%psolv%stc(5,:)=[ 0,-1, 0]
      this%psolv%stc(6,:)=[ 0, 0,+1]
      this%psolv%stc(7,:)=[ 0, 0,-1]
      
      ! Setup the scaled Laplacian operator from incomp metrics: lap(*)=-vol*div(grad(*))
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Set Laplacian
               this%psolv%opr(1,i,j,k)=this%divp_x(1,i,j,k)*this%divu_x(-1,i+1,j,k)+&
               &                       this%divp_x(0,i,j,k)*this%divu_x( 0,i  ,j,k)+&
               &                       this%divp_y(1,i,j,k)*this%divv_y(-1,i,j+1,k)+&
               &                       this%divp_y(0,i,j,k)*this%divv_y( 0,i,j  ,k)+&
               &                       this%divp_z(1,i,j,k)*this%divw_z(-1,i,j,k+1)+&
               &                       this%divp_z(0,i,j,k)*this%divw_z( 0,i,j,k  )
               this%psolv%opr(2,i,j,k)=this%divp_x(1,i,j,k)*this%divu_x( 0,i+1,j,k)
               this%psolv%opr(3,i,j,k)=this%divp_x(0,i,j,k)*this%divu_x(-1,i  ,j,k)
               this%psolv%opr(4,i,j,k)=this%divp_y(1,i,j,k)*this%divv_y( 0,i,j+1,k)
               this%psolv%opr(5,i,j,k)=this%divp_y(0,i,j,k)*this%divv_y(-1,i,j  ,k)
               this%psolv%opr(6,i,j,k)=this%divp_z(1,i,j,k)*this%divw_z( 0,i,j,k+1)
               this%psolv%opr(7,i,j,k)=this%divp_z(0,i,j,k)*this%divw_z(-1,i,j,k  )
               ! Scale it by the cell volume
               this%psolv%opr(:,i,j,k)=-this%psolv%opr(:,i,j,k)*this%cfg%vol(i,j,k)
            end do
         end do
      end do
      
      ! Initialize the pressure Poisson solver
      call this%psolv%init_solver(pressure_ils)
      call this%psolv%update_solver()
      
      ! Set 7-pt stencil map for the velocity solver
      this%implicit%stc(1,:)=[ 0, 0, 0]
      this%implicit%stc(2,:)=[+1, 0, 0]
      this%implicit%stc(3,:)=[-1, 0, 0]
      this%implicit%stc(4,:)=[ 0,+1, 0]
      this%implicit%stc(5,:)=[ 0,-1, 0]
      this%implicit%stc(6,:)=[ 0, 0,+1]
      this%implicit%stc(7,:)=[ 0, 0,-1]
      
      ! Set the diagonal to 1 to make sure all cells participate in solver
      this%implicit%opr(1,:,:,:)=1.0_WP
      
      ! Initialize the implicit velocity solver
      call this%implicit%init_solver(implicit_ils)
      
   end subroutine setup
   
   
   !> Add a boundary condition
   subroutine add_bcond(this,name,type,dir,canCorrect,locator,inBalance)
      use string,   only: lowercase
      use messager, only: die
      implicit none
      class(incomp), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer,  intent(in) :: type
      character(len=2), intent(in) :: dir
      logical,  intent(in) :: canCorrect
      interface
         logical function locator(pargrid,ind1,ind2,ind3)
            use pgrid_class, only: pgrid
            class(pgrid), intent(in) :: pargrid
            integer, intent(in) :: ind1,ind2,ind3
         end function locator
      end interface
      logical, optional :: inBalance
      type(bcond), pointer :: new_bc
      integer :: i,j,k,n
      
      ! Prepare new bcond
      allocate(new_bc)
      new_bc%name=trim(adjustl(name))
      new_bc%type=type
      new_bc%canCorrect=canCorrect
      select case (lowercase(dir))
      case ('+x','x+','xp','px'); new_bc%dir=1
      case ('-x','x-','xm','mx'); new_bc%dir=2
      case ('+y','y+','yp','py'); new_bc%dir=3
      case ('-y','y-','ym','my'); new_bc%dir=4
      case ('+z','z+','zp','pz'); new_bc%dir=5
      case ('-z','z-','zm','mz'); new_bc%dir=6
      case default; call die('[incomp add_bcond] Unknown bcond direction')
      end select
      new_bc%itr=iterator(this%cfg,new_bc%name,locator)
      if (present(inBalance)) then
         new_bc%inBalance=inBalance
      else
         new_bc%inBalance=.true.
      end if
      
      ! Insert it up front
      new_bc%next=>this%first_bc
      this%first_bc=>new_bc
      
      ! Increment bcond counter
      this%nbc=this%nbc+1
      
      ! Now adjust the metrics accordingly
      select case (new_bc%type)
      case (dirichlet)
         ! Implement based on bcond direction
         select case (new_bc%dir)
         case (1) ! +x
            do n=1,new_bc%itr%n_
               i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
               this%umask(i+1,j,k)=2
            end do
         case (2) ! -x
            do n=1,new_bc%itr%n_
               i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
               this%umask(i,j,k)=2
            end do
         case (3) ! +y
            do n=1,new_bc%itr%n_
               i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
               this%vmask(i,j+1,k)=2
            end do
         case (4) ! -y
            do n=1,new_bc%itr%n_
               i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
               this%vmask(i,j,k)=2
            end do
         case (5) ! +z
            do n=1,new_bc%itr%n_
               i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
               this%wmask(i,j,k+1)=2
            end do
         case (6) ! -z
            do n=1,new_bc%itr%n_
               i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
               this%wmask(i,j,k)=2
            end do
         end select
         
      case (neumann)
         ! Not yet implemented
      case (clipped_neumann)
         ! Not yet implemented
      case (convective)
         ! Not yet implemented
      case default
         call die('[incomp apply_bcond] Unknown bcond type')
      end select
   
   end subroutine add_bcond
   
   
   !> Get a boundary condition
   subroutine get_bcond(this,name,my_bc)
      implicit none
      class(incomp), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(bcond), pointer, intent(out) :: my_bc
      my_bc=>this%first_bc
      search: do while (associated(my_bc))
         if (trim(my_bc%name).eq.trim(name)) exit search
         my_bc=>my_bc%next
      end do search
   end subroutine get_bcond
   
   
   !> Enforce boundary condition
   subroutine apply_bcond(this,t,dt)
      use messager, only: die
      use mpi_f08,  only: MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), intent(in) :: t,dt
      integer :: i,j,k,n,ierr
      type(bcond), pointer :: my_bc
      real(WP) :: conv_vel,my_conv_vel
      
      ! ! First enfore zero velocity at walls
      ! do k=this%cfg%kmin_,this%cfg%kmax_
      !    do j=this%cfg%jmin_,this%cfg%jmax_
      !       do i=this%cfg%imin_,this%cfg%imax_
      !          if (minval(this%cfg%VF(i-1:i,j,k)).lt.10.0_WP*epsilon(1.0_WP)) this%U(i,j,k)=0.0_WP
      !          if (minval(this%cfg%VF(i,j-1:j,k)).lt.10.0_WP*epsilon(1.0_WP)) this%V(i,j,k)=0.0_WP
      !          if (minval(this%cfg%VF(i,j,k-1:k)).lt.10.0_WP*epsilon(1.0_WP)) this%W(i,j,k)=0.0_WP
      !       end do
      !    end do
      ! end do
      ! ! Sync fields
      ! call this%cfg%sync(this%U)
      ! call this%cfg%sync(this%V)
      ! call this%cfg%sync(this%W)
      
      ! Traverse bcond list
      my_bc=>this%first_bc
      do while (associated(my_bc))
         
         ! Only processes inside the bcond work here
         if (my_bc%itr%amIn) then
            
            ! Select appropriate action based on the bcond type
            select case (my_bc%type)
               
            case (dirichlet)           ! Apply Dirichlet conditions
               
               ! This is done by the user directly
               ! Unclear whether we want to do this within the solver...
               
            case (neumann)             ! Apply Neumann condition
               
               ! Implement based on bcond direction
               select case (my_bc%dir)
               case (1) ! Neumann in +x
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i+1,j,k)=this%U(i,j,k)
                     this%V(i+1,j:j+1,k)=this%V(i,j:j+1,k)
                     this%W(i+1,j,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (2) ! Neumann in -x
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i,j,k)=this%U(i+1,j,k)
                     this%V(i-1,j:j+1,k)=this%V(i,j:j+1,k)
                     this%W(i-1,j,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (3) ! Neumann in +y
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j+1,k)=this%U(i:i+1,j,k)
                     this%V(i,j+1,k)=this%V(i,j,k)
                     this%W(i,j+1,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (4) ! Neumann in -y
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j-1,k)=this%U(i:i+1,j,k)
                     this%V(i,j,k)=this%V(i,j+1,k)
                     this%W(i,j-1,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (5) ! Neumann in +z
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j,k+1)=this%U(i:i+1,j,k)
                     this%V(i,j:j+1,k+1)=this%V(i,j:j+1,k)
                     this%W(i,j,k+1)=this%W(i,j,k)
                  end do
               case (6) ! Neumann in -z
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j,k-1)=this%U(i:i+1,j,k)
                     this%V(i,j:j+1,k-1)=this%V(i,j:j+1,k)
                     this%W(i,j,k)=this%W(i,j,k+1)
                  end do
               end select
               
            case (clipped_neumann)     ! Apply clipped Neumann condition
               
               ! Implement based on bcond direction
               select case (my_bc%dir)
               case (1) ! Neumann in +x
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i+1,j,k)=max(this%U(i,j,k),0.0_WP)
                     this%V(i+1,j:j+1,k)=this%V(i,j:j+1,k)
                     this%W(i+1,j,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (2) ! Neumann in -x
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i,j,k)=min(this%U(i+1,j,k),0.0_WP)
                     this%V(i-1,j:j+1,k)=this%V(i,j:j+1,k)
                     this%W(i-1,j,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (3) ! Neumann in +y
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j+1,k)=this%U(i:i+1,j,k)
                     this%V(i,j+1,k)=max(this%V(i,j,k),0.0_WP)
                     this%W(i,j+1,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (4) ! Neumann in -y
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j-1,k)=this%U(i:i+1,j,k)
                     this%V(i,j,k)=min(this%V(i,j+1,k),0.0_WP)
                     this%W(i,j-1,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (5) ! Neumann in +z
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j,k+1)=this%U(i:i+1,j,k)
                     this%V(i,j:j+1,k+1)=this%V(i,j:j+1,k)
                     this%W(i,j,k+1)=max(this%W(i,j,k),0.0_WP)
                  end do
               case (6) ! Neumann in -z
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j,k-1)=this%U(i:i+1,j,k)
                     this%V(i,j:j+1,k-1)=this%V(i,j:j+1,k)
                     this%W(i,j,k)=min(this%W(i,j,k+1),0.0_WP)
                  end do
               end select
            
            case (convective)   ! Apply convective condition
               
               ! Implement based on bcond direction
               select case (my_bc%dir)
               case (1) ! Neumann in +x
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i+1,j,k)=this%U(i,j,k)
                     this%V(i+1,j:j+1,k)=this%V(i,j:j+1,k)
                     this%W(i+1,j,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (2) ! Neumann in -x
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i,j,k)=this%U(i+1,j,k)
                     this%V(i-1,j:j+1,k)=this%V(i,j:j+1,k)
                     this%W(i-1,j,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (3) ! Neumann in +y
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j+1,k)=this%U(i:i+1,j,k)
                     this%V(i,j+1,k)=this%V(i,j,k)
                     this%W(i,j+1,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (4) ! Neumann in -y
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j-1,k)=this%U(i:i+1,j,k)
                     this%V(i,j,k)=this%V(i,j+1,k)
                     this%W(i,j-1,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (5) ! Neumann in +z
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j,k+1)=this%U(i:i+1,j,k)
                     this%V(i,j:j+1,k+1)=this%V(i,j:j+1,k)
                     this%W(i,j,k+1)=this%W(i,j,k)
                  end do
               case (6) ! Neumann in -z
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j,k-1)=this%U(i:i+1,j,k)
                     this%V(i,j:j+1,k-1)=this%V(i,j:j+1,k)
                     this%W(i,j,k)=this%W(i,j,k+1)
                  end do
               end select
               
               ! ! Get convective velocity
               ! my_conv_vel=0.0_WP
               ! do n=1,my_bc%itr%no_
               !    ! Get the indices
               !    i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
               !    ! Find maximum velocity
               !    my_conv_vel(1)=max(my_conv_vel(1),this%U(i+si,j+sj,k+sk))
               !    my_conv_vel(2)=max(my_conv_vel(2),this%V(i+si,j+sj,k+sk))
               !    my_conv_vel(3)=max(my_conv_vel(3),this%W(i+si,j+sj,k+sk))
               ! end do
               ! call MPI_ALLREDUCE(my_conv_vel,conv_vel,3,MPI_REAL_WP,MPI_MAX,my_bc%itr%comm,ierr)
               
               
            case default
               call die('[incomp apply_bcond] Unknown bcond type')
            end select
            
         end if
         
         ! Sync full fields after each bcond - this should be optimized
         call this%cfg%sync(this%U)
         call this%cfg%sync(this%V)
         call this%cfg%sync(this%W)
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
   end subroutine apply_bcond
   
   
   !> Calculate the explicit momentum time derivative based on U/V/W/P
   subroutine get_dmomdt(this,drhoUdt,drhoVdt,drhoWdt)
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoUdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoVdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoWdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,ii,jj,kk
      real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ
      
      ! Allocate flux arrays
      allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! Flux of rhoU
      do kk=this%cfg%kmin_,this%cfg%kmax_+1
         do jj=this%cfg%jmin_,this%cfg%jmax_+1
            do ii=this%cfg%imin_,this%cfg%imax_+1
               ! Fluxes on x-face
               i=ii-1; j=jj; k=kk
               FX(i,j,k)=-this%rho * sum(this%itpu_x(:,i,j,k)*this%U(i:i+1,j,k))*sum(this%itpu_x(:,i,j,k)*this%U(i:i+1,j,k)) &
               &         +this%visc*(sum(this%grdu_x(:,i,j,k)*this%U(i:i+1,j,k))+sum(this%grdu_x(:,i,j,k)*this%U(i:i+1,j,k)))&
               &         -this%P(i,j,k)
               ! Fluxes on y-face
               i=ii  ; j=jj; k=kk
               FY(i,j,k)=-this%rho * sum(this%itpu_y(:,i,j,k)*this%U(i,j-1:j,k))*sum(this%itpv_x(:,i,j,k)*this%V(i-1:i,j,k)) &
               &         +this%visc*(sum(this%grdu_y(:,i,j,k)*this%U(i,j-1:j,k))+sum(this%grdv_x(:,i,j,k)*this%V(i-1:i,j,k)))
               ! Fluxes on z-face
               i=ii  ; j=jj; k=kk
               FZ(i,j,k)=-this%rho * sum(this%itpu_z(:,i,j,k)*this%U(i,j,k-1:k))*sum(this%itpw_x(:,i,j,k)*this%W(i-1:i,j,k)) &
               &         +this%visc*(sum(this%grdu_z(:,i,j,k)*this%U(i,j,k-1:k))+sum(this%grdw_x(:,i,j,k)*this%W(i-1:i,j,k)))
            end do
         end do
      end do
      ! Time derivative of rhoU
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               drhoUdt(i,j,k)=sum(this%divu_x(:,i,j,k)*FX(i-1:i,j,k))+&
               &              sum(this%divu_y(:,i,j,k)*FY(i,j:j+1,k))+&
               &              sum(this%divu_z(:,i,j,k)*FZ(i,j,k:k+1))
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(drhoUdt)
      
      ! Flux of rhoV
      do kk=this%cfg%kmin_,this%cfg%kmax_+1
         do jj=this%cfg%jmin_,this%cfg%jmax_+1
            do ii=this%cfg%imin_,this%cfg%imax_+1
               ! Fluxes on x-face
               i=ii; j=jj  ; k=kk
               FX(i,j,k)=-this%rho * sum(this%itpv_x(:,i,j,k)*this%V(i-1:i,j,k))*sum(this%itpu_y(:,i,j,k)*this%U(i,j-1:j,k)) &
               &         +this%visc*(sum(this%grdv_x(:,i,j,k)*this%V(i-1:i,j,k))+sum(this%grdu_y(:,i,j,k)*this%U(i,j-1:j,k)))
               ! Fluxes on y-face
               i=ii; j=jj-1; k=kk
               FY(i,j,k)=-this%rho * sum(this%itpv_y(:,i,j,k)*this%V(i,j:j+1,k))*sum(this%itpv_y(:,i,j,k)*this%V(i,j:j+1,k)) &
               &         +this%visc*(sum(this%grdv_y(:,i,j,k)*this%V(i,j:j+1,k))+sum(this%grdv_y(:,i,j,k)*this%V(i,j:j+1,k)))&
               &         -this%P(i,j,k)
               ! Fluxes on z-face
               i=ii; j=jj  ; k=kk
               FZ(i,j,k)=-this%rho * sum(this%itpv_z(:,i,j,k)*this%V(i,j,k-1:k))*sum(this%itpw_y(:,i,j,k)*this%W(i,j-1:j,k)) &
               &         +this%visc*(sum(this%grdv_z(:,i,j,k)*this%V(i,j,k-1:k))+sum(this%grdw_y(:,i,j,k)*this%W(i,j-1:j,k)))
            end do
         end do
      end do
      ! Time derivative of rhoV
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               drhoVdt(i,j,k)=sum(this%divv_x(:,i,j,k)*FX(i:i+1,j,k))+&
               &              sum(this%divv_y(:,i,j,k)*FY(i,j-1:j,k))+&
               &              sum(this%divv_z(:,i,j,k)*FZ(i,j,k:k+1))
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(drhoVdt)
      
      ! Flux of rhoW
      do kk=this%cfg%kmin_,this%cfg%kmax_+1
         do jj=this%cfg%jmin_,this%cfg%jmax_+1
            do ii=this%cfg%imin_,this%cfg%imax_+1
               ! Fluxes on x-face
               i=ii; j=jj; k=kk
               FX(i,j,k)=-this%rho * sum(this%itpw_x(:,i,j,k)*this%W(i-1:i,j,k))*sum(this%itpu_z(:,i,j,k)*this%U(i,j,k-1:k)) &
               &         +this%visc*(sum(this%grdw_x(:,i,j,k)*this%W(i-1:i,j,k))+sum(this%grdu_z(:,i,j,k)*this%U(i,j,k-1:k)))
               ! Fluxes on y-face
               i=ii; j=jj; k=kk
               FY(i,j,k)=-this%rho * sum(this%itpw_y(:,i,j,k)*this%W(i,j-1:j,k))*sum(this%itpv_z(:,i,j,k)*this%V(i,j,k-1:k)) &
               &         +this%visc*(sum(this%grdw_y(:,i,j,k)*this%W(i,j-1:j,k))+sum(this%grdv_z(:,i,j,k)*this%V(i,j,k-1:k)))
               ! Fluxes on z-face
               i=ii; j=jj; k=kk-1
               FZ(i,j,k)=-this%rho * sum(this%itpw_z(:,i,j,k)*this%W(i,j,k:k+1))*sum(this%itpw_z(:,i,j,k)*this%W(i,j,k:k+1)) &
               &         +this%visc*(sum(this%grdw_z(:,i,j,k)*this%W(i,j,k:k+1))+sum(this%grdw_z(:,i,j,k)*this%W(i,j,k:k+1)))&
               &         -this%P(i,j,k)
            end do
         end do
      end do
      ! Time derivative of rhoW
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               drhoWdt(i,j,k)=sum(this%divw_x(:,i,j,k)*FX(i:i+1,j,k))+&
               &              sum(this%divw_y(:,i,j,k)*FY(i,j:j+1,k))+&
               &              sum(this%divw_z(:,i,j,k)*FZ(i,j,k-1:k))
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(drhoWdt)
      
      ! Deallocate flux arrays
      deallocate(FX,FY,FZ)
      
   end subroutine get_dmomdt
   
   
   !> Calculate the velocity divergence based on U/V/W
   subroutine get_div(this)
      implicit none
      class(incomp), intent(inout) :: this
      integer :: i,j,k
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%div(i,j,k)=sum(this%divp_x(:,i,j,k)*this%U(i:i+1,j,k))+&
               &               sum(this%divp_y(:,i,j,k)*this%V(i,j:j+1,k))+&
               &               sum(this%divp_z(:,i,j,k)*this%W(i,j,k:k+1))
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(this%div)
   end subroutine get_div
   
   
   !> Calculate the pressure gradient based on P
   subroutine get_pgrad(this,P,Pgradx,Pgrady,Pgradz)
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: P      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Pgradx !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Pgrady !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Pgradz !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               Pgradx(i,j,k)=sum(this%divu_x(:,i,j,k)*P(i-1:i,j,k))
               Pgrady(i,j,k)=sum(this%divv_y(:,i,j,k)*P(i,j-1:j,k))
               Pgradz(i,j,k)=sum(this%divw_z(:,i,j,k)*P(i,j,k-1:k))
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(Pgradx)
      call this%cfg%sync(Pgrady)
      call this%cfg%sync(Pgradz)
   end subroutine get_pgrad
   
   
   !> Calculate the interpolated velocity
   subroutine interp_vel(this,Ui,Vi,Wi)
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Ui !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Vi !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Wi !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               Ui(i,j,k)=sum(this%itpu_x(:,i,j,k)*this%U(i:i+1,j,k))
               Vi(i,j,k)=sum(this%itpv_y(:,i,j,k)*this%V(i,j:j+1,k))
               Wi(i,j,k)=sum(this%itpw_z(:,i,j,k)*this%W(i,j,k:k+1))
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(Ui)
      call this%cfg%sync(Vi)
      call this%cfg%sync(Wi)
   end subroutine interp_vel
   
   
   !> Calculate the CFL
   subroutine get_cfl(this,dt,cflc,cfl)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), intent(out) :: cflc
      real(WP), optional :: cfl
      integer :: i,j,k,ierr
      real(WP) :: my_CFLc_x,my_CFLc_y,my_CFLc_z,my_CFLv_x,my_CFLv_y,my_CFLv_z
      
      ! Set the CFLs to zero
      my_CFLc_x=0.0_WP; my_CFLc_y=0.0_WP; my_CFLc_z=0.0_WP
      my_CFLv_x=0.0_WP; my_CFLv_y=0.0_WP; my_CFLv_z=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               my_CFLc_x=max(my_CFLc_x,abs(this%U(i,j,k))*this%cfg%dxmi(i))
               my_CFLc_y=max(my_CFLc_y,abs(this%V(i,j,k))*this%cfg%dymi(j))
               my_CFLc_z=max(my_CFLc_z,abs(this%W(i,j,k))*this%cfg%dzmi(k))
               my_CFLv_x=max(my_CFLv_x,4.0_WP*this%visc*this%cfg%dxi(i)**2/this%rho)
               my_CFLv_y=max(my_CFLv_y,4.0_WP*this%visc*this%cfg%dyi(j)**2/this%rho)
               my_CFLv_z=max(my_CFLv_z,4.0_WP*this%visc*this%cfg%dzi(k)**2/this%rho)
            end do
         end do
      end do
      my_CFLc_x=my_CFLc_x*dt; my_CFLc_y=my_CFLc_y*dt; my_CFLc_z=my_CFLc_z*dt
      my_CFLv_x=my_CFLv_x*dt; my_CFLv_y=my_CFLv_y*dt; my_CFLv_z=my_CFLv_z*dt
      
      ! Get the parallel max
      call MPI_ALLREDUCE(my_CFLc_x,this%CFLc_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLc_y,this%CFLc_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLc_z,this%CFLc_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLv_x,this%CFLv_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLv_y,this%CFLv_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLv_z,this%CFLv_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
      ! Return the maximum convective CFL
      cflc=max(this%CFLc_x,this%CFLc_y,this%CFLc_z)
      
      ! If asked for, also return the maximum overall CFL
      if (present(CFL)) cfl =max(this%CFLc_x,this%CFLc_y,this%CFLc_z,this%CFLv_x,this%CFLv_y,this%CFLv_z)
      
   end subroutine get_cfl
   
   
   !> Calculate the max of our fields
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(incomp), intent(inout) :: this
      integer :: i,j,k,ierr
      real(WP) :: my_Umax,my_Vmax,my_Wmax,my_Pmax,my_divmax
      
      ! Set all to zero
      my_Umax=0.0_WP; my_Vmax=0.0_WP; my_Wmax=0.0_WP; my_Pmax=0.0_WP; my_divmax=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               my_Umax  =max(my_Umax  ,abs(this%U(i,j,k)  ))
               my_Vmax  =max(my_Vmax  ,abs(this%V(i,j,k)  ))
               my_Wmax  =max(my_Wmax  ,abs(this%W(i,j,k)  ))
               my_Pmax  =max(my_Pmax  ,abs(this%P(i,j,k)  ))
               my_divmax=max(my_divmax,abs(this%div(i,j,k)))
            end do
         end do
      end do
      
      ! Get the parallel max
      call MPI_ALLREDUCE(my_Umax  ,this%Umax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Vmax  ,this%Vmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Wmax  ,this%Wmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Pmax  ,this%Pmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_divmax,this%divmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
   end subroutine get_max
   
   
   !> Compute MFR through all bcs
   subroutine get_mfr(this)
      use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(incomp), intent(inout) :: this
      integer :: i,j,k,n,ibc,ierr
      type(bcond), pointer :: my_bc
      real(WP), dimension(:), allocatable :: my_mfr,my_area
      real(WP), dimension(:), allocatable :: canCorrect
      
      ! Ensure this%mfr is of proper size
      if (size(this%mfr).ne.this%nbc.or..not.allocated(this%mfr)) then
         if (allocated(this%mfr)) deallocate(this%mfr)
         allocate(this%mfr(this%nbc))
      end if
      
      ! Ensure this%area is of proper size
      if (size(this%area).ne.this%nbc.or..not.allocated(this%area)) then
         if (allocated(this%area)) deallocate(this%area)
         allocate(this%area(this%nbc))
      end if
      
      ! Allocate temp array for communication
      allocate(my_mfr(this%nbc))
      allocate(my_area(this%nbc))
      allocate(canCorrect(this%nbc))
      
      ! Traverse bcond list and integrate local outgoing MFR
      my_bc=>this%first_bc; ibc=1
      do while (associated(my_bc))
         
         ! Set zero local MFR and area
         my_mfr(ibc)=0.0_WP
         my_area(ibc)=0.0_WP
         if (my_bc%canCorrect) then
            canCorrect(ibc)=1.0_WP
         else
            canCorrect(ibc)=0.0_WP
         end if
         
         ! Check that the bcond participates in volume balance
         if (my_bc%inBalance) then
            
            ! Only processes inside the bcond have a non-zero MFR
            if (my_bc%itr%amIn) then
               
               ! Implement based on bcond direction, loop over interior only
               select case (my_bc%dir)
               case (1) ! BC in +x
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     my_mfr(ibc)=my_mfr(ibc)+this%rho*this%U(i+1,j,k)*this%cfg%dy(j)*this%cfg%dz(k)
                     my_area(ibc)=my_area(ibc)+this%cfg%dy(j)*this%cfg%dz(k)
                  end do
               case (2) ! BC in -x
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     my_mfr(ibc)=my_mfr(ibc)-this%rho*this%U(i  ,j,k)*this%cfg%dy(j)*this%cfg%dz(k)
                     my_area(ibc)=my_area(ibc)+this%cfg%dy(j)*this%cfg%dz(k)
                  end do
               case (3) ! BC in +y
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     my_mfr(ibc)=my_mfr(ibc)+this%rho*this%V(i,j+1,k)*this%cfg%dx(i)*this%cfg%dz(k)
                     my_area(ibc)=my_area(ibc)+this%cfg%dx(i)*this%cfg%dz(k)
                  end do
               case (4) ! BC in -y
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     my_mfr(ibc)=my_mfr(ibc)-this%rho*this%V(i,j  ,k)*this%cfg%dx(i)*this%cfg%dz(k)
                     my_area(ibc)=my_area(ibc)+this%cfg%dx(i)*this%cfg%dz(k)
                  end do
               case (5) ! BC in +z
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     my_mfr(ibc)=my_mfr(ibc)+this%rho*this%W(i,j,k+1)*this%cfg%dy(j)*this%cfg%dx(i)
                     my_area(ibc)=my_area(ibc)+this%cfg%dy(j)*this%cfg%dx(i)
                  end do
               case (6) ! BC in -z
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     my_mfr(ibc)=my_mfr(ibc)-this%rho*this%W(i,j,k  )*this%cfg%dy(j)*this%cfg%dx(i)
                     my_area(ibc)=my_area(ibc)+this%cfg%dy(j)*this%cfg%dx(i)
                  end do
               end select
               
            end if
            
         end if
         
         ! Move on to the next bcond
         my_bc=>my_bc%next; ibc=ibc+1
         
      end do
      
      ! Sum up all values
      call MPI_ALLREDUCE(my_mfr ,this%mfr ,this%nbc,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_area,this%area,this%nbc,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      
      ! Compute the correctable area
      this%correctable_area=sum(this%area*canCorrect)
      
      ! Deallocate temp array
      deallocate(my_mfr,my_area,canCorrect)
      
   end subroutine get_mfr
   
   
   !> Correct MFR through correctable bconds
   subroutine correct_mfr(this)
      use mpi_f08,  only: MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      class(incomp), intent(inout) :: this
      real(WP) :: mfr_error,vel_correction
      integer :: i,j,k,n
      type(bcond), pointer :: my_bc
      
      ! Evaluate MFR mismatch and velocity correction
      call this%get_mfr()
      mfr_error=sum(this%mfr)
      if (abs(mfr_error).lt.10.0_WP*epsilon(1.0_WP).or.abs(this%correctable_area).lt.10.0_WP*epsilon(1.0_WP)) return
      vel_correction=-mfr_error/(this%rho*this%correctable_area)
      
      ! Traverse bcond list and correct bcond MFR
      my_bc=>this%first_bc
      do while (associated(my_bc))
         
         ! Only processes inside correctable bcond need to work
         if (my_bc%itr%amIn.and.my_bc%canCorrect) then
            
            ! Implement based on bcond direction, loop over all cell
            select case (my_bc%dir)
            case (1) ! BC in +x
               do n=1,my_bc%itr%no_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  this%U(i+1,j,k)=this%U(i+1,j,k)+vel_correction
               end do
            case (2) ! BC in -x
               do n=1,my_bc%itr%no_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  this%U(i,j,k)=this%U(i,j,k)-vel_correction
               end do
            case (3) ! BC in +y
               do n=1,my_bc%itr%no_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  this%V(i,j+1,k)=this%V(i,j+1,k)+vel_correction
               end do
            case (4) ! BC in -y
               do n=1,my_bc%itr%no_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  this%V(i,j,k)=this%V(i,j,k)-vel_correction
               end do
            case (5) ! BC in +z
               do n=1,my_bc%itr%no_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  this%W(i,j,k+1)=this%W(i,j,k+1)+vel_correction
               end do
            case (6) ! BC in -z
               do n=1,my_bc%itr%no_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  this%W(i,j,k)=this%W(i,j,k)-vel_correction
               end do
            end select
            
         end if
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
   end subroutine correct_mfr
   
   
   !> Shift pressure to ensure zero average
   subroutine shift_p(this,pressure)
      use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(incomp), intent(in) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: pressure !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP) :: vol_tot,pressure_tot,my_vol_tot,my_pressure_tot
      integer :: i,j,k,ierr
      
      ! Loop over domain and integrate volume and pressure
      my_vol_tot=0.0_WP
      my_pressure_tot=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               my_vol_tot     =my_vol_tot     +this%cfg%vol(i,j,k)*this%cfg%VF(i,j,k)
               my_pressure_tot=my_pressure_tot+this%cfg%vol(i,j,k)*this%cfg%VF(i,j,k)*pressure(i,j,k)
            end do
         end do
      end do
      call MPI_ALLREDUCE(my_vol_tot     ,vol_tot     ,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_pressure_tot,pressure_tot,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      pressure_tot=pressure_tot/vol_tot
      
      ! Shift the pressure
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%cfg%VF(i,j,k).gt.0.0_WP) pressure(i,j,k)=pressure(i,j,k)-pressure_tot
            end do
         end do
      end do
      call this%cfg%sync(pressure)
      
   end subroutine shift_p
   
   
   !> Solve for implicit velocity residual
   subroutine solve_implicit(this,dt,resU,resV,resW)
      use ils_class, only: amg
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resU !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resV !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resW !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP) :: rhoUp,rhoUm,rhoVp,rhoVm,rhoWp,rhoWm
      
      ! Solve implicit U problem
      this%implicit%opr(1,:,:,:)=this%rho; this%implicit%opr(2:,:,:,:)=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               rhoUp=this%rho*sum(this%itpu_x(:,i  ,j,k)*this%U(i  :i+1,j,k))*2.0_WP
               rhoUm=this%rho*sum(this%itpu_x(:,i-1,j,k)*this%U(i-1:i  ,j,k))*2.0_WP
               rhoVp=this%rho*sum(this%itpv_x(:,i,j+1,k)*this%V(i-1:i,j+1,k))
               rhoVm=this%rho*sum(this%itpv_x(:,i,j  ,k)*this%V(i-1:i,j  ,k))
               rhoWp=this%rho*sum(this%itpw_x(:,i,j,k+1)*this%W(i-1:i,j,k+1))
               rhoWm=this%rho*sum(this%itpw_x(:,i,j,k  )*this%W(i-1:i,j,k  ))
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+0.5_WP*dt*(this%divu_x( 0,i,j,k)*this%itpu_x( 0,i  ,j,k)*rhoUp+&
               &                                                                this%divu_x(-1,i,j,k)*this%itpu_x(+1,i-1,j,k)*rhoUm+&
               &                                                                this%divu_y(+1,i,j,k)*this%itpu_y(-1,i,j+1,k)*rhoVp+&
               &                                                                this%divu_y( 0,i,j,k)*this%itpu_y( 0,i,j  ,k)*rhoVm+&
               &                                                                this%divu_z(+1,i,j,k)*this%itpu_z(-1,i,j,k+1)*rhoWp+&
               &                                                                this%divu_z( 0,i,j,k)*this%itpu_z( 0,i,j,k  )*rhoWm)
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)+0.5_WP*dt*(this%divu_x( 0,i,j,k)*this%itpu_x(+1,i  ,j,k)*rhoUp)
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)+0.5_WP*dt*(this%divu_x(-1,i,j,k)*this%itpu_x( 0,i-1,j,k)*rhoUm)
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)+0.5_WP*dt*(this%divu_y(+1,i,j,k)*this%itpu_y( 0,i,j+1,k)*rhoVp)
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)+0.5_WP*dt*(this%divu_y( 0,i,j,k)*this%itpu_y(-1,i,j  ,k)*rhoVm)
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)+0.5_WP*dt*(this%divu_z(+1,i,j,k)*this%itpu_z( 0,i,j,k+1)*rhoWp)
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)+0.5_WP*dt*(this%divu_z( 0,i,j,k)*this%itpu_z(-1,i,j,k  )*rhoWm)
            end do
         end do
      end do
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%divu_x( 0,i,j,k)*this%visc*this%grdu_x( 0,i  ,j,k)+&
               &                                                                this%divu_x(-1,i,j,k)*this%visc*this%grdu_x(+1,i-1,j,k)+&
               &                                                                this%divu_y(+1,i,j,k)*this%visc*this%grdu_y(-1,i,j+1,k)+&
               &                                                                this%divu_y( 0,i,j,k)*this%visc*this%grdu_y( 0,i,j  ,k)+&
               &                                                                this%divu_z(+1,i,j,k)*this%visc*this%grdu_z(-1,i,j,k+1)+&
               &                                                                this%divu_z( 0,i,j,k)*this%visc*this%grdu_z( 0,i,j,k  ))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%divu_x( 0,i,j,k)*this%visc*this%grdu_x(+1,i  ,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%divu_x(-1,i,j,k)*this%visc*this%grdu_x( 0,i-1,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%divu_y(+1,i,j,k)*this%visc*this%grdu_y( 0,i,j+1,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%divu_y( 0,i,j,k)*this%visc*this%grdu_y(-1,i,j  ,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%divu_z(+1,i,j,k)*this%visc*this%grdu_z( 0,i,j,k+1))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%divu_z( 0,i,j,k)*this%visc*this%grdu_z(-1,i,j,k  ))
            end do
         end do
      end do
      call this%implicit%update_solver()
      this%implicit%rhs=resU
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resU=this%implicit%sol
      
      ! Solve implicit V problem
      this%implicit%opr(1,:,:,:)=this%rho; this%implicit%opr(2:,:,:,:)=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               rhoUp=this%rho*sum(this%itpu_y(:,i+1,j,k)*this%U(i+1,j-1:j,k))
               rhoUm=this%rho*sum(this%itpu_y(:,i  ,j,k)*this%U(i  ,j-1:j,k))
               rhoVp=this%rho*sum(this%itpv_y(:,i,j  ,k)*this%V(i,j  :j+1,k))*2.0_WP
               rhoVm=this%rho*sum(this%itpv_y(:,i,j-1,k)*this%V(i,j-1:j  ,k))*2.0_WP
               rhoWp=this%rho*sum(this%itpw_y(:,i,j,k+1)*this%W(i,j-1:j,k+1))
               rhoWm=this%rho*sum(this%itpw_y(:,i,j,k  )*this%W(i,j-1:j,k  ))
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+0.5_WP*dt*(this%divv_x(+1,i,j,k)*this%itpv_x(-1,i+1,j,k)*rhoUp+&
               &                                                                this%divv_x( 0,i,j,k)*this%itpv_x( 0,i  ,j,k)*rhoUm+&
               &                                                                this%divv_y( 0,i,j,k)*this%itpv_y( 0,i,j  ,k)*rhoVp+&
               &                                                                this%divv_y(-1,i,j,k)*this%itpv_y(+1,i,j-1,k)*rhoVm+&
               &                                                                this%divv_z(+1,i,j,k)*this%itpv_z(-1,i,j,k+1)*rhoWp+&
               &                                                                this%divv_z( 0,i,j,k)*this%itpv_z( 0,i,j,k  )*rhoWm)
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)+0.5_WP*dt*(this%divv_x(+1,i,j,k)*this%itpv_x( 0,i+1,j,k)*rhoUp)
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)+0.5_WP*dt*(this%divv_x( 0,i,j,k)*this%itpv_x(-1,i  ,j,k)*rhoUm)
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)+0.5_WP*dt*(this%divv_y( 0,i,j,k)*this%itpv_y(+1,i,j  ,k)*rhoVp)
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)+0.5_WP*dt*(this%divv_y(-1,i,j,k)*this%itpv_y( 0,i,j-1,k)*rhoVm)
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)+0.5_WP*dt*(this%divv_z(+1,i,j,k)*this%itpv_z( 0,i,j,k+1)*rhoWp)
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)+0.5_WP*dt*(this%divv_z( 0,i,j,k)*this%itpv_z(-1,i,j,k  )*rhoWm)
            end do
         end do
      end do
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%divv_x(+1,i,j,k)*this%visc*this%grdv_x(-1,i+1,j,k)+&
               &                                                                this%divv_x( 0,i,j,k)*this%visc*this%grdv_x( 0,i  ,j,k)+&
               &                                                                this%divv_y( 0,i,j,k)*this%visc*this%grdv_y( 0,i,j  ,k)+&
               &                                                                this%divv_y(-1,i,j,k)*this%visc*this%grdv_y(+1,i,j-1,k)+&
               &                                                                this%divv_z(+1,i,j,k)*this%visc*this%grdv_z(-1,i,j,k+1)+&
               &                                                                this%divv_z( 0,i,j,k)*this%visc*this%grdv_z( 0,i,j,k  ))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%divv_x(+1,i,j,k)*this%visc*this%grdv_x( 0,i+1,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%divv_x( 0,i,j,k)*this%visc*this%grdv_x(-1,i  ,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%divv_y( 0,i,j,k)*this%visc*this%grdv_y(+1,i,j  ,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%divv_y(-1,i,j,k)*this%visc*this%grdv_y( 0,i,j-1,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%divv_z(+1,i,j,k)*this%visc*this%grdv_z( 0,i,j,k+1))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%divv_z( 0,i,j,k)*this%visc*this%grdv_z(-1,i,j,k  ))
            end do
         end do
      end do
      call this%implicit%update_solver()
      this%implicit%rhs=resV
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resV=this%implicit%sol
      
      ! Solve implicit W problem
      this%implicit%opr(1,:,:,:)=this%rho; this%implicit%opr(2:,:,:,:)=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               rhoUp=this%rho*sum(this%itpu_z(:,i+1,j,k)*this%U(i+1,j,k-1:k))
               rhoUm=this%rho*sum(this%itpu_z(:,i  ,j,k)*this%U(i  ,j,k-1:k))
               rhoVp=this%rho*sum(this%itpv_z(:,i,j+1,k)*this%V(i,j+1,k-1:k))
               rhoVm=this%rho*sum(this%itpv_z(:,i,j  ,k)*this%V(i,j  ,k-1:k))
               rhoWp=this%rho*sum(this%itpw_z(:,i,j,k  )*this%W(i,j,k  :k+1))*2.0_WP
               rhoWm=this%rho*sum(this%itpw_z(:,i,j,k-1)*this%W(i,j,k-1:k  ))*2.0_WP
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+0.5_WP*dt*(this%divw_x(+1,i,j,k)*this%itpw_x(-1,i+1,j,k)*rhoUp+&
               &                                                                this%divw_x( 0,i,j,k)*this%itpw_x( 0,i  ,j,k)*rhoUm+&
               &                                                                this%divw_y(+1,i,j,k)*this%itpw_y(-1,i,j+1,k)*rhoVp+&
               &                                                                this%divw_y( 0,i,j,k)*this%itpw_y( 0,i,j  ,k)*rhoVm+&
               &                                                                this%divw_z( 0,i,j,k)*this%itpw_z( 0,i,j,k  )*rhoWp+&
               &                                                                this%divw_z(-1,i,j,k)*this%itpw_z(+1,i,j,k-1)*rhoWm)
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)+0.5_WP*dt*(this%divw_x(+1,i,j,k)*this%itpw_x( 0,i+1,j,k)*rhoUp)
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)+0.5_WP*dt*(this%divw_x( 0,i,j,k)*this%itpw_x(-1,i  ,j,k)*rhoUm)
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)+0.5_WP*dt*(this%divw_y(+1,i,j,k)*this%itpw_y( 0,i,j+1,k)*rhoVp)
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)+0.5_WP*dt*(this%divw_y( 0,i,j,k)*this%itpw_y(-1,i,j  ,k)*rhoVm)
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)+0.5_WP*dt*(this%divw_z( 0,i,j,k)*this%itpw_z(+1,i,j,k  )*rhoWp)
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)+0.5_WP*dt*(this%divw_z(-1,i,j,k)*this%itpw_z( 0,i,j,k-1)*rhoWm)
            end do
         end do
      end do
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%divw_x(+1,i,j,k)*this%visc*this%grdw_x(-1,i+1,j,k)+&
               &                                                                this%divw_x( 0,i,j,k)*this%visc*this%grdw_x( 0,i  ,j,k)+&
               &                                                                this%divw_y(+1,i,j,k)*this%visc*this%grdw_y(-1,i,j+1,k)+&
               &                                                                this%divw_y( 0,i,j,k)*this%visc*this%grdw_y( 0,i,j  ,k)+&
               &                                                                this%divw_z( 0,i,j,k)*this%visc*this%grdw_z( 0,i,j,k  )+&
               &                                                                this%divw_z(-1,i,j,k)*this%visc*this%grdw_z(+1,i,j,k-1))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%divw_x(+1,i,j,k)*this%visc*this%grdw_x( 0,i+1,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%divw_x( 0,i,j,k)*this%visc*this%grdw_x(-1,i  ,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%divw_y(+1,i,j,k)*this%visc*this%grdw_y( 0,i,j+1,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%divw_y( 0,i,j,k)*this%visc*this%grdw_y(-1,i,j  ,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%divw_z( 0,i,j,k)*this%visc*this%grdw_z(+1,i,j,k  ))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%divw_z(-1,i,j,k)*this%visc*this%grdw_z( 0,i,j,k-1))
            end do
         end do
      end do
      call this%implicit%update_solver()
      this%implicit%rhs=resW
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resW=this%implicit%sol
      
      ! Sync up all residuals
      call this%cfg%sync(resU)
      call this%cfg%sync(resV)
      call this%cfg%sync(resW)
      
   end subroutine solve_implicit
   
   
   !> Print out info for incompressible flow solver
   subroutine incomp_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(incomp), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Incompressible solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
         write(output_unit,'(" >   density = ",es12.5)') this%rho
         write(output_unit,'(" > viscosity = ",es12.5)') this%visc
      end if
      
   end subroutine incomp_print
   
   
end module incomp_class