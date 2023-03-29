!> 3D FFT pressure solver concept is defined by extension of the linsol class
!> This solver is specifically intended to be a FFT-based pressure Poisson solver
!> for 3D periodic uniform computational domains decomposed in at most 2 directions
!> uses fft3d for Fourier transform

module fftsolver3d_class
  use precision,    only: WP
  use config_class, only: config
  use fft3d_class,  only: fft3d
  use string,       only: str_short
  use linsol_class, only: linsol
  use, intrinsic :: iso_c_binding
  implicit none
  private


  ! Expose type/constructor/methods
  public :: fftsolver3d


  !> fftsolver object definition
  type, extends(linsol) :: fftsolver3d

    ! FFT object
    class(fft3d), allocatable :: dft

    !> Unstrided arrays
    complex(WP), dimension(:,:,:), allocatable :: factored_operator
    complex(WP), dimension(:,:,:), allocatable :: transformed_rhs

  contains

    procedure :: print_short=>fftsolver3d_print_short !< One-line printing of solver status
    procedure :: print=>fftsolver3d_print             !< Long-form printing of solver status
    procedure :: log=>fftsolver3d_log                 !< Long-form logging of solver status
    procedure :: init=>fftsolver3d_init               !< Grid and stencil initialization - done once for the grid and stencil
    procedure :: setup=>fftsolver3d_setup             !< Solver setup (every time the operator changes)
    procedure :: solve=>fftsolver3d_solve             !< Execute solver (assumes new RHS and initial guess at every call)
    procedure :: destroy=>fftsolver3d_destroy

  end type fftsolver3d


  !> Declare fftsolver3d constructor
  interface fftsolver3d; procedure fftsolver3d_from_args; end interface fftsolver3d;


contains


  !> Constructor for a fftsolver3d object
  function fftsolver3d_from_args(cfg,name,nst) result(self)
    use messager, only: die
    implicit none
    type(fftsolver3d) :: self
    class(config), target, intent(in) :: cfg
    character(len=*), intent(in) :: name
    integer, intent(in) :: nst

    ! Link the config and store the name
    self%cfg=>cfg
    self%name=trim(adjustl(name))

    ! Set solution method - not used
    self%method=0

    ! Set up stencil size and map
    self%nst=nst
    allocate(self%stc(1:self%nst,1:3))
    self%stc=0

    ! Allocate operator, rhs, and sol arrays
    allocate(self%opr(self%nst,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%opr=0.0_WP
    allocate(self%rhs(         self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rhs=0.0_WP
    allocate(self%sol(         self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%sol=0.0_WP

    ! Zero out some info
    self%it=1
    self%aerr=0.0_WP
    self%rerr=0.0_WP

    ! Allocate unstrided arrays
    allocate(self%factored_operator(self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_))
    allocate(self%transformed_rhs  (self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_))

    ! Setup is not done
    self%setup_done=.false.

  end function fftsolver3d_from_args


  !> Initialize grid and stencil - done at start-up only as long as the stencil or cfg does not change
  !> When calling, zero values of this%opr(1,:,:,:) indicate cells that do not participate in the solver
  !> Only the stencil needs to be defined at this point
  subroutine fftsolver3d_init(this)
    use messager, only: die
    implicit none
    class(fftsolver3d), intent(inout) :: this
    integer :: st,stx1,stx2,sty1,sty2,stz1,stz2
    include 'fftw3.f03'

    ! From the provided stencil, generate an inverse map
    stx1=minval(this%stc(:,1)); stx2=maxval(this%stc(:,1));
    sty1=minval(this%stc(:,2)); sty2=maxval(this%stc(:,2));
    stz1=minval(this%stc(:,3)); stz2=maxval(this%stc(:,3));
    allocate(this%stmap(stx1:stx2,sty1:sty2,stz1:stz2)); this%stmap=0
    do st=1,this%nst
      this%stmap(this%stc(st,1),this%stc(st,2),this%stc(st,3))=st
    end do

    ! Initialize DFT
    this%dft=fft3d(this%cfg)

  end subroutine fftsolver3d_init


  !> Setup solver - done everytime the operator changes
  subroutine fftsolver3d_setup(this)
    use messager, only: die
    use mpi_f08,  only: MPI_BCAST,MPI_ALLREDUCE,MPI_INTEGER,MPI_SUM
    use parallel, only: MPI_REAL_WP
    implicit none
    class(fftsolver3d), intent(inout) :: this
    real(WP), dimension(1:this%nst) :: ref_opr
    integer :: i,j,k,n,ierr

    ! Check operator is circulant
    checkcirc: block
      logical :: circulant
      real(WP) :: circtol
      if (this%cfg%amRoot) ref_opr=this%opr(:,this%cfg%imin,this%cfg%jmin,this%cfg%kmin)
      call MPI_BCAST(ref_opr,this%nst,MPI_REAL_WP,0,this%cfg%comm,ierr)
      circulant=.true.
      circtol = 6.0_WP*epsilon(1.0_WP)/this%cfg%min_meshsize**4
      do k=this%cfg%kmin_,this%cfg%kmax_
        do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
            if (any(abs(this%opr(:,i,j,k)-ref_opr).gt.circtol)) circulant=.false.
          end do
        end do
      end do
      if (.not.circulant) call die('[fftsolver3d setup] operator must be uniform in space')
    end block checkcirc

    ! Build the operator
    this%factored_operator=0.0_WP
    do n=1,this%nst
      i=modulo(this%stc(n,1)-this%cfg%imin+1,this%cfg%nx)+this%cfg%imin
      j=modulo(this%stc(n,2)-this%cfg%jmin+1,this%cfg%ny)+this%cfg%jmin
      k=modulo(this%stc(n,3)-this%cfg%kmin+1,this%cfg%nz)+this%cfg%kmin
      if (this%cfg%imin_.le.i.and.i.le.this%cfg%imax_.and.                    &
        & this%cfg%jmin_.le.j.and.j.le.this%cfg%jmax_.and.                    &
        & this%cfg%kmin_.le.k.and.k.le.this%cfg%kmax_)                        &
      this%factored_operator(i,j,k)=this%factored_operator(i,j,k)+ref_opr(n)
    end do

    ! Take transform of operator
    call this%dft%forward_transform(this%factored_operator)

    ! Make zero wavenumber not zero
    ! Setting this to one has the nice side effect of returning a solution
    ! with the same integral
    if (this%dft%oddball) this%factored_operator(this%cfg%imin_,this%cfg%jmin_,&
      this%cfg%kmin_)=(1.0_WP, 0.0_WP)

    ! Make sure other wavenumbers are not close to zero
    i=count(abs(this%factored_operator).lt.1000_WP*epsilon(1.0_WP))
    call MPI_ALLREDUCE(i,j,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr)
    if (j.gt.0) then
      write(*,*) j
      call die('[fftsolver3d setup] elements of transformed operator near zero')
    end if

    ! Divide now instead of later
    this%factored_operator=1.0_WP/this%factored_operator

    ! Check for division issues
    i=count(isnan(abs(this%factored_operator)))
    call MPI_ALLREDUCE(i,j,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr)
    if (j.gt.0) call die('[fftsolver3d setup] elements of transformed operator are NaN')

    ! Set setup-flag to true
    this%setup_done=.true.

  end subroutine fftsolver3d_setup


  !> Solve the linear system iteratively
  subroutine fftsolver3d_solve(this)
    use messager, only: die
    use param,    only: verbose
    implicit none
    class(fftsolver3d), intent(inout) :: this

    ! Check that setup was done
    if (.not.this%setup_done) call die('[fftsolver3d solve] Solver has not &
      &been setup.')

    ! Copy to unstrided array
    this%transformed_rhs=this%rhs(this%cfg%imin_:this%cfg%imax_,              &
      this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)

    ! Forward transform
    call this%dft%forward_transform(this%transformed_rhs)

    ! Divide
    this%transformed_rhs=this%transformed_rhs*this%factored_operator

    ! Backward transform
    call this%dft%backward_transform(this%transformed_rhs)

    ! Copy to strided output
    this%sol(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,     &
      this%cfg%kmin_:this%cfg%kmax_)=realpart(this%transformed_rhs)

    ! Sync the solution vector
    call this%cfg%sync(this%sol)

    ! If verbose run, log and or print info
    if (verbose.gt.0) call this%log
    if (verbose.gt.1) call this%print_short

  end subroutine fftsolver3d_solve


  !> Log fftsolver3d info
  subroutine fftsolver3d_log(this)
    use string,   only: str_long
    use messager, only: log
    implicit none
    class(fftsolver3d), intent(in) :: this
    character(len=str_long) :: message

    if (this%cfg%amRoot) then
      write(message,'("fftsolver3d solver [",a,"] for config [",a,"]")')    &
        trim(this%name), trim(this%cfg%name)
      call log(message)
    end if

  end subroutine fftsolver3d_log


  !> Print fftsolver3d info to the screen
  subroutine fftsolver3d_print(this)
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    class(fftsolver3d), intent(in) :: this

    if (this%cfg%amRoot) then
      write(output_unit,'("fftsolver3d solver [",a,"] for config [",a,"]")') &
        trim(this%name), trim(this%cfg%name)
    end if

  end subroutine fftsolver3d_print


  !> Short print of fftsolver3d info to the screen
  subroutine fftsolver3d_print_short(this)
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    class(fftsolver3d), intent(in) :: this

    if (this%cfg%amRoot) write(output_unit,'("fftsolver3d solver [",a16,"] &
      &for config [",a16,"]")') trim(this%name),trim(this%cfg%name)

  end subroutine fftsolver3d_print_short

  subroutine fftsolver3d_destroy(this)
    implicit none
    class(fftsolver3d), intent(inout) :: this

    deallocate(this%stc,this%opr,this%rhs,this%sol,this%factored_operator,    &
      this%transformed_rhs,this%stmap)

  end subroutine fftsolver3d_destroy

end module fftsolver3d_class
