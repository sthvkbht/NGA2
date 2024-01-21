!> Ghost point immersed boundary class
!> Provides support for Neumann boundary conditions with IBM
module gp_class
  use precision,      only: WP
  use ibconfig_class, only: ibconfig
  implicit none
  private

  ! Expose type/constructor/methods
  public :: gpibm

  ! List of known available bcond types for this solver
  integer, parameter, public :: dirichlet=1         !< Dirichlet condition
  integer, parameter, public :: neumann=2           !< Normal gradient

  !> Basic image point definition
  type :: image
     integer, dimension(3,2) :: st_sc                    !< Interpolation stencil extent for scalars
     integer , dimension(3) :: ind                       !< Index of cell containing image point
     real(WP), dimension(3) :: pos                       !< Coordinates of image point
     real(WP), dimension(2,2,2) :: itp_sc                !< Interpolation weights for scalars
  end type image

  !> Basic ghost point definition
  type :: ghost
     integer, dimension(3) :: ind                        !< Index of cell containing image point
     type(image) :: im                                   !< Associated image point
  end type ghost

  !> SGS model object definition
  type :: gpibm

     ! This is our config
     class(ibconfig), pointer :: cfg                     !< This is the config the model is build for

     ! These are the ghost points
     integer :: no                                       !< Number of overlapping ghost points
     integer :: ngp                                      !< Number of ghost points for our solver
     type(ghost), dimension(:), allocatable :: gp        !< Array of ghost points
     integer, dimension(:,:,:), allocatable :: label     !< Integer array used for labeling ghost/image points

   contains

     procedure :: update                                 !< Update ghost point information
     procedure :: apply_bcond                            !< Apply specified boundary conditions

  end type gpibm


  !> Declare model constructor
  interface gpibm
     procedure constructor
  end interface gpibm

contains


  !> Default constructor for model
  function constructor(cfg,no) result(self)
    implicit none
    type(gpibm) :: self
    class(ibconfig), target, intent(in) :: cfg
    integer, intent(in) :: no
    integer :: i,j,k,ii,jj,kk,i1,i2,j1,j2,k1,k2,n

    ! Point to config object
    self%cfg=>cfg

    ! Reset ghost point list
    self%no=no
    self%ngp=0
    if (allocated(self%gp)) deallocate(self%gp)

    ! Allocate label array (0=fluid cell, 1=ghost point, 2=image point)
    allocate(self%label(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%label=0

  end function constructor


  !> Updates list of ghost points and associated image point data
  subroutine update(this)
    use messager, only: die
    use param, only: verbose
    implicit none
    class(gpibm), intent(inout) :: this
    integer :: i,j,k,ii,jj,kk,i1,i2,j1,j2,k1,k2,n
    real(WP) :: dist,buf
    real(WP), dimension(3) :: pos,pos_im
    real(WP), dimension(2,2,2) :: alpha3D,delta3D,dist3D,eta3D
    real(WP), parameter :: eps=1.0e-10_WP
    logical :: success

    ! Reset ghost points
    this%ngp=0
    if (allocated(this%gp)) deallocate(this%gp)

    ! Identify ghost points
    this%ngp=0
    do k=this%cfg%kmin_,this%cfg%kmax_
       do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
             if (this%cfg%Gib(i,j,k).ge.0.0_WP) cycle
             i1=i-this%no; i2=i+this%no; i1=max(i1,this%cfg%imino); i2=min(i2,this%cfg%imaxo)
             j1=j-this%no; j2=j+this%no; j1=max(j1,this%cfg%jmino); j2=min(j2,this%cfg%jmaxo)
             k1=k-this%no; k2=k+this%no; k1=max(k1,this%cfg%kmino); k2=min(k2,this%cfg%kmaxo)
             success=.false.
             do kk = k1, k2
                do jj = j1, j2
                   do ii = i1, i2
                      if (this%cfg%Gib(ii,jj,kk).gt.0.0_WP.and.this%cfg%VF(ii,jj,kk).gt.0.0_WP) success=.true.
                   end do
                end do
             end do
             if (success) this%ngp=this%ngp+1
          end do
       end do
    end do
    if (this%ngp.gt.0) allocate(this%gp(this%ngp))
    ! Store ghost points and associated image points
    n=0
    do k=this%cfg%kmin_,this%cfg%kmax_
       do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
             if (this%cfg%Gib(i,j,k).ge.0.0_WP) cycle
             i1=i-this%no; i2=i+this%no; i1=max(i1,this%cfg%imino); i2=min(i2,this%cfg%imaxo)
             j1=j-this%no; j2=j+this%no; j1=max(j1,this%cfg%jmino); j2=min(j2,this%cfg%jmaxo)
             k1=k-this%no; k2=k+this%no; k1=max(k1,this%cfg%kmino); k2=min(k2,this%cfg%kmaxo)
             success=.false.
             do kk = k1, k2
                do jj = j1, j2
                   do ii = i1, i2
                      if (this%cfg%Gib(ii,jj,kk).gt.0.0_WP.and.this%cfg%VF(ii,jj,kk).gt.0.0_WP) success=.true.
                   end do
                end do
             end do
             if (success) then
                n=n+1
                ! Store index of ghost points
                this%gp(n)%ind(1)=i; this%gp(n)%ind(2)=j; this%gp(n)%ind(3)=k
                pos=(/this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)/)
                ! Get position and index of image point
                dist=abs(this%cfg%Gib(i,j,k))
                pos_im=pos+2.0_WP*dist*this%cfg%Nib(:,i,j,k)
                ! Find right i index
                i1=i
                do while (pos_im(1)-this%cfg%xm(i1  ).lt.0.0_WP.and.i1  .gt.this%cfg%imino_); i1=i1-1; end do
                do while (pos_im(1)-this%cfg%xm(i1+1).ge.0.0_WP.and.i1+1.lt.this%cfg%imaxo_); i1=i1+1; end do
                i2=i1+1     
                ! Find right j index
                j1=j
                do while (pos_im(2)-this%cfg%ym(j1  ).lt.0.0_WP.and.j1  .gt.this%cfg%jmino_); j1=j1-1; end do
                do while (pos_im(2)-this%cfg%ym(j1+1).ge.0.0_WP.and.j1+1.lt.this%cfg%jmaxo_); j1=j1+1; end do
                j2=j1+1      
                ! Find right k index
                k1=k
                do while (pos_im(3)-this%cfg%zm(k1  ).lt.0.0_WP.and.k1  .gt.this%cfg%kmino_); k1=k1-1; end do
                do while (pos_im(3)-this%cfg%zm(k1+1).ge.0.0_WP.and.k1+1.lt.this%cfg%kmaxo_); k1=k1+1; end do
                k2=k1+1      
                ! Check if the image point is inside the levelset
                alpha3D = 1.0_WP
                if (this%cfg%Gib(i1,j1,k1).le.0.0_WP) alpha3D(1,1,1) = 0.0_WP
                if (this%cfg%Gib(i2,j1,k1).le.0.0_WP) alpha3D(2,1,1) = 0.0_WP
                if (this%cfg%Gib(i1,j2,k1).le.0.0_WP) alpha3D(1,2,1) = 0.0_WP
                if (this%cfg%Gib(i2,j2,k1).le.0.0_WP) alpha3D(2,2,1) = 0.0_WP
                if (this%cfg%Gib(i1,j1,k2).le.0.0_WP) alpha3D(1,1,2) = 0.0_WP
                if (this%cfg%Gib(i2,j1,k2).le.0.0_WP) alpha3D(2,1,2) = 0.0_WP
                if (this%cfg%Gib(i1,j2,k2).le.0.0_WP) alpha3D(1,2,2) = 0.0_WP
                if (this%cfg%Gib(i2,j2,k2).le.0.0_WP) alpha3D(2,2,2) = 0.0_WP
                buf = sum(alpha3D)
                if (buf.gt.0.0_WP) then
                   ! Get interpolation weights at image point (Chaudhuri et al. 2011, JCP)
                   dist3D(1,1,1)=sqrt((this%cfg%xm(i1)-pos_im(1))**2+(this%cfg%ym(j1)-pos_im(2))**2+(this%cfg%zm(k1)-pos_im(3))**2)
                   dist3D(2,1,1)=sqrt((this%cfg%xm(i2)-pos_im(1))**2+(this%cfg%ym(j1)-pos_im(2))**2+(this%cfg%zm(k1)-pos_im(3))**2)
                   dist3D(1,2,1)=sqrt((this%cfg%xm(i1)-pos_im(1))**2+(this%cfg%ym(j2)-pos_im(2))**2+(this%cfg%zm(k1)-pos_im(3))**2)
                   dist3D(2,2,1)=sqrt((this%cfg%xm(i2)-pos_im(1))**2+(this%cfg%ym(j2)-pos_im(2))**2+(this%cfg%zm(k1)-pos_im(3))**2)
                   dist3D(1,1,2)=sqrt((this%cfg%xm(i1)-pos_im(1))**2+(this%cfg%ym(j1)-pos_im(2))**2+(this%cfg%zm(k2)-pos_im(3))**2)
                   dist3D(2,1,2)=sqrt((this%cfg%xm(i2)-pos_im(1))**2+(this%cfg%ym(j1)-pos_im(2))**2+(this%cfg%zm(k2)-pos_im(3))**2)
                   dist3D(1,2,2)=sqrt((this%cfg%xm(i1)-pos_im(1))**2+(this%cfg%ym(j2)-pos_im(2))**2+(this%cfg%zm(k2)-pos_im(3))**2)
                   dist3D(2,2,2)=sqrt((this%cfg%xm(i2)-pos_im(1))**2+(this%cfg%ym(j2)-pos_im(2))**2+(this%cfg%zm(k2)-pos_im(3))**2)
                   delta3D=0.0_WP
                   if (dist3D(1,1,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,1,1)=1.0_WP
                   else if (dist3D(2,1,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,1,1)=1.0_WP
                   else if (dist3D(1,2,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,2,1)=1.0_WP
                   else if (dist3D(2,2,1).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,2,1)=1.0_WP
                   else if (dist3D(1,1,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,1,2)=1.0_WP
                   else if (dist3D(2,1,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,1,2)=1.0_WP
                   else if (dist3D(1,2,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(1,2,2)=1.0_WP
                   else if (dist3D(2,2,2).le.eps*this%cfg%min_meshsize) then
                      delta3D(2,2,2)=1.0_WP
                   else
                      eta3D=1.0_WP/dist3D**2
                      buf=sum(eta3D*alpha3D)
                      delta3D=alpha3D*eta3D/buf
                   end if
                   ! Store image point data
                   this%gp(n)%im%pos=pos_im
                   this%gp(n)%im%ind(1)=i1; this%gp(n)%im%ind(2)=j1; this%gp(n)%im%ind(3)=k1
                   this%gp(n)%im%ind=this%cfg%get_ijk_global(this%gp(n)%im%pos,this%gp(n)%im%ind)
                   this%gp(n)%im%st_sc(1,1)=i1; this%gp(n)%im%st_sc(1,2)=i2
                   this%gp(n)%im%st_sc(2,1)=j1; this%gp(n)%im%st_sc(2,2)=j2
                   this%gp(n)%im%st_sc(3,1)=k1; this%gp(n)%im%st_sc(3,2)=k2
                   this%gp(n)%im%itp_sc=delta3D
                else
                   call die('[gpibm update] Problem computing interpolation weight for image point')
                end if
             end if
          end do
       end do
    end do

    ! Update label array
    this%label=0
    do n=1,this%ngp
       i=this%gp(n)%ind(1); j=this%gp(n)%ind(2); k=this%gp(n)%ind(3)
       this%label(i,j,k)=1 !< Ghost point
       i=this%gp(n)%im%ind(1); j=this%gp(n)%im%ind(2); k=this%gp(n)%im%ind(3)
       this%label(i,j,k)=2 !< Image points
    end do

    ! Log/screen output
    logging: block
      use, intrinsic :: iso_fortran_env, only: output_unit
      use param,    only: verbose
      use messager, only: log
      use string,   only: str_long
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Ghost point solver on partitioned grid [",a,"]: ",i0," ghost points found")') trim(this%cfg%name),this%ngp
         if (verbose.gt.1) write(output_unit,'(a)') trim(message)
         if (verbose.gt.0) call log(message)
      end if
    end block logging

  end subroutine update


  !> Enforce boundary condition at ghost points of A using boundary point BP
  subroutine apply_bcond(this,type,BP,A)
    use messager,       only: die
    implicit none
    class(gpibm), intent(inout) :: this
    integer, intent(in) :: type
    real(WP), intent(in) :: BP
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: A !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    integer :: i,j,k,n
    real(WP) :: IP,dist
    select case (type)
    case (dirichlet)
       do n=1,this%ngp
          i=this%gp(n)%ind(1); j=this%gp(n)%ind(2); k=this%gp(n)%ind(3)
          IP=sum(this%gp(n)%im%itp_sc*A(this%gp(n)%im%st_sc(1,1:2),this%gp(n)%im%st_sc(2,1:2),this%gp(n)%im%st_sc(3,1:2)))
          A(i,j,k)=2.0_WP*BP-IP
       end do
    case (neumann)
       do n=1,this%ngp
          i=this%gp(n)%ind(1); j=this%gp(n)%ind(2); k=this%gp(n)%ind(3)
          IP=sum(this%gp(n)%im%itp_sc*A(this%gp(n)%im%st_sc(1,1:2),this%gp(n)%im%st_sc(2,1:2),this%gp(n)%im%st_sc(3,1:2)))
          dist=abs(this%cfg%Gib(i,j,k))
          A(i,j,k)=IP-2.0_WP*dist*BP
       end do
    case default
       call die('[gpibm apply_bcond] Unknown bcond type')
    end select
    call this%cfg%sync(A)
  end subroutine apply_bcond


end module gp_class
