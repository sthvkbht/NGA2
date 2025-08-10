!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   use string,       only: str_medium
   implicit none
   private
   
   !> Single config
   type(config), public :: cfg

   !> Max bounding radius
   real(WP), public :: R0
   
   public :: geometry_init
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      use mathtools, only: Pi
      implicit none
      type(sgrid) :: grid
      
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz,no
         real(WP) :: Lx,Ly,Lz,L,dx
         real(WP), dimension(:), allocatable :: x,y,z
         character(len=str_medium) :: shape
         
         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))

         ! Read in rigid body definition
         call param_read('Body shape',shape)
         call param_read('Body side length',L)
         
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz
         end do

         ! Get max extend for communication
         select case(trim(adjustl(shape)))
         case ('rod')
            R0=0.5_WP*L
         case('pentagon')
            R0=0.5_WP*L/sin(Pi/5.0_WP)
         end select
         dx=min(Lx/real(nx,WP),Ly/real(ny,WP))
         no=ceiling(2.0_WP*R0/dx)
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=no,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.true.,name='Box')
         
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         
         ! Read in partition
         call param_read('Partition',partition,short='p')
         
         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)
         
      end block create_cfg
      
      
      ! Create masks for this config
      create_walls: block
        integer :: i,j,k
        cfg%VF=0.0_WP
        cfg%VF(cfg%imin_:cfg%imax_,cfg%jmin_:cfg%jmax_,cfg%kmin_:cfg%kmax_)=1.0_WP
        call cfg%sync(cfg%VF)
      end block create_walls
      
      
   end subroutine geometry_init
   
   
end module geometry
