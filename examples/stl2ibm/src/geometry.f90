!> Various definitions and tools for initializing NGA2 config
module geometry
   use ibconfig_class, only: ibconfig
   use precision,      only: WP,SP
   implicit none
   private
   
   !> Single config
   type(ibconfig), public :: cfg

   public :: geometry_init
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      integer :: nt
      real(WP) :: Lx,Ly,Lz
      real(WP) :: xmin,xmax,ymin,ymax,zmin,zmax
      type triangle_type
         real(WP), dimension(3) :: norm
         real(WP), dimension(3) :: v1
         real(WP), dimension(3) :: v2
         real(WP), dimension(3) :: v3
      end type triangle_type
      type(triangle_type), dimension(:), allocatable :: t

      ! Read the STL file and get domain extents and levelset
      read_stl: block
        use messager, only: die
        use mathtools, only: normalize
        integer :: i,iunit,ierr
        real(WP), parameter :: scaling=1.0e-3_WP
        real(SP) :: rbuf
        character(len=80) :: stlfile,cbuf
        character(len= 2) :: padding
        ! Open the file
        call param_read('STL file',stlfile)
        open(newunit=iunit,file=trim(stlfile),status="old",action="read",access="stream",iostat=ierr)
        if (ierr.ne.0) call die('[read_stl] Could not open file: '//trim(stlfile))
        ! Read the header
        read(iunit) cbuf
        read(iunit) nt
        ! Read the triangle data and scale to meters
        allocate(t(nt))
        do i=1,nt
           read(iunit) rbuf; t(i)%norm(1)=real(rbuf,WP)
           read(iunit) rbuf; t(i)%norm(2)=real(rbuf,WP)
           read(iunit) rbuf; t(i)%norm(3)=real(rbuf,WP)
           read(iunit) rbuf; t(i)%v1(1)=real(rbuf,WP)
           read(iunit) rbuf; t(i)%v1(2)=real(rbuf,WP)
           read(iunit) rbuf; t(i)%v1(3)=real(rbuf,WP)
           read(iunit) rbuf; t(i)%v2(1)=real(rbuf,WP)
           read(iunit) rbuf; t(i)%v2(2)=real(rbuf,WP)
           read(iunit) rbuf; t(i)%v2(3)=real(rbuf,WP)
           read(iunit) rbuf; t(i)%v3(1)=real(rbuf,WP)
           read(iunit) rbuf; t(i)%v3(2)=real(rbuf,WP)
           read(iunit) rbuf; t(i)%v3(3)=real(rbuf,WP)
           read(iunit) padding
        end do
        ! Close the file
        close(iunit)
        ! Scale
        do i=1,nt
           t(i)%v1=t(i)%v1*scaling
           t(i)%v2=t(i)%v2*scaling
           t(i)%v3=t(i)%v3*scaling
           t(i)%norm=normalize(t(i)%norm*scaling)
        end do
        ! Get extents
        xmin=minval(t(:)%v1(1)); xmin=min(xmin,minval(t(:)%v2(1))); xmin=min(xmin,minval(t(:)%v3(1)))
        xmax=maxval(t(:)%v1(1)); xmax=max(xmin,maxval(t(:)%v2(1))); xmax=max(xmax,maxval(t(:)%v3(1)))
        ymin=minval(t(:)%v1(2)); ymin=min(ymin,minval(t(:)%v2(2))); ymin=min(ymin,minval(t(:)%v3(2)))
        ymax=maxval(t(:)%v1(2)); ymax=max(ymin,maxval(t(:)%v2(2))); ymax=max(ymax,maxval(t(:)%v3(2)))
        zmin=minval(t(:)%v1(3)); zmin=min(zmin,minval(t(:)%v2(3))); zmin=min(zmin,minval(t(:)%v3(3)))
        zmax=maxval(t(:)%v1(3)); zmax=max(zmin,maxval(t(:)%v2(3))); zmax=max(zmax,maxval(t(:)%v3(3)))
        Lx=xmax-xmin; Ly=ymax-ymin; Lz=zmax-zmin
      end block read_stl
      
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz
         real(WP), dimension(:), allocatable :: x,y,z
         
         ! Read in grid definition
         call param_read('nx',nx); allocate(x(nx+1))
         call param_read('ny',ny); allocate(y(ny+1))
         call param_read('nz',nz); allocate(z(nz+1))
         
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx+xmin
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly+ymin
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz+zmin
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='vts')
         
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid
         cfg=ibconfig(grp=group,decomp=partition,grid=grid)
      end block create_cfg
      
      
      ! Create IB walls for this config
      create_walls: block
         use parallel,       only: amroot
         use mathtools,      only: triproj
         use ibconfig_class, only: bigot,sharp
         integer :: i,j,k,n,nod,count,myloc,togo,prog,iratio_new,iratio_old
         real(WP) ::mydist,newdist,tmp,ang
         real(WP), dimension(3) :: c,mynorm,myproj,newproj,v1,v2
         real(WP), dimension(3,1500) :: proj,norm
         real(WP), dimension(1500) :: proj_count
         real(WP), parameter :: eps=1.0e-9_WP
         ! Create IB field
         prog=0
         togo=cfg%nxo_*cfg%nyo_*cfg%nzo_
         iratio_old=-1
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Prepare projections
                  count=0
                  mydist=huge(1.0_WP)

                  ! Get centroid information
                  c=(/cfg%xm(i),cfg%ym(j),cfg%zm(k)/)

                  ! Loop over triangles and check distance
                  do n=1,nt

                     ! If all vertices of triangle > mydist+eps, no need to project
                     if (sum((c - t(n)%v1)**2) .gt. (mydist+eps)**2 .and.     &
                         sum((c - t(n)%v2)**2) .gt. (mydist+eps)**2 .and.     &
                         sum((c - t(n)%v3)**2) .gt. (mydist+eps)**2) cycle

                     ! Normalize triangle normal
                     tmp=sqrt(dot_product(t(n)%norm,t(n)%norm))
                     t(n)%norm=t(n)%norm/tmp

                     call triproj(c,t(n)%v1,t(n)%v2,t(n)%v3,newproj)
                     newdist=sqrt(dot_product(c-newproj,c-newproj))

                     ! Check point
                     if (newdist.lt.mydist-eps) then
                        ! new closest point
                        mydist=newdist
                        myproj=newproj
                        mynorm=t(n)%norm
                     else if (newdist.lt.mydist+eps) then
                        ! Choose better normal
                        if ( sqrt(abs(dot_product(t(n)%norm,(newproj-c)/sqrt(dot_product(newproj-c,newproj-c))))) .gt. &
                             sqrt(abs(dot_product(mynorm,   (myproj -c)/sqrt(dot_product(myproj -c,myproj -c))))) ) then
                           mynorm=t(n)%norm
                           myproj=newproj
                           mydist=newdist
                        end if
                     end if

                  end do

                  ! Postprocess distance
                  cfg%Gib(i,j,k)=mydist

                  ! Get sign based on normal
                  tmp=dot_product(myproj-c,mynorm)
                  if (tmp.gt.0.0_WP) cfg%Gib(i,j,k)=-cfg%Gib(i,j,k)
                  !cfg%Gib(i,j,k)=-cfg%Gib(i,j,k)

                  ! Add point to counter
                  if (amRoot) then
                     prog=prog+1
                     iratio_new=int(real(prog,WP)/real(togo,WP)*100.0_WP)
                     if (iratio_new.gt.iratio_old) then
                        iratio_old=iratio_new
                        write(*,'(i3,x,a1)') iratio_new,'%'
                     end if
                  end if
               end do
            end do
         end do
         ! Get normal vector
         call cfg%calculate_normal()
         ! Get VF field
         call cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)
      end block create_walls

      ! Clean up
      deallocate(t)
      
   end subroutine geometry_init
   
   
end module geometry
