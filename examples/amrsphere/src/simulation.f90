!> AMR Sphere - Flow over a sphere with Immersed Boundary
!> Inflow/outflow in X, periodic in Y/Z
module simulation
   use precision,         only: WP
   use amrviz_class,      only: amrviz
   use amrgrid_class,     only: amrgrid
   use amrincomp_class,   only: amrincomp
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use messager,          only: log
   implicit none
   private

   public :: simulation_init, simulation_run, simulation_final

   ! Grid
   type(amrgrid), target :: amr

   ! Solver data
   type(timetracker) :: time
   type(amrincomp), target :: fs
   type(amrdata) :: resU,resV,resW

   ! Fluid volume fraction for IB forcing
   type(amrdata) :: VF

   ! Visualization
   type(amrviz) :: viz
   type(event) :: viz_evt

   ! Regrid parameters
   type(event) :: regrid_evt
   real(WP) :: Re_tag=huge(1.0_WP)

   ! Monitoring
   type(monitor) :: mfile,cflfile,gridfile

contains

   !> Levelset function for sphere
   function sphere_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=sqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)-0.5_WP
   end function sphere_levelset

   !> Tagger for this case based on velocity gradient magnitude and distance to sphere surface
   subroutine my_tagger(solver,lvl,tags_ptr,time)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrincomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags_ptr
      real(WP), intent(in) :: time
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW
      real(WP) :: dx,dy,dz,dist,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,gradU_mag,Re_cell
      integer :: i,j,k
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      tags=tags_ptr
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         bx=mfi%tilebox()
         tagarr=>tags%dataPtr(mfi)
         pU=>solver%U%mf(lvl)%dataptr(mfi)
         pV=>solver%V%mf(lvl)%dataptr(mfi)
         pW=>solver%W%mf(lvl)%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Skip if close to outflow (within 5 units)
            if (solver%amr%xlo+(real(i,WP)+0.5_WP)*dx.gt.solver%amr%xhi-5.0_WP) cycle
            ! Diagonal gradients
            dudx=(pU(i+1,j,k,1)-pU(i,j,k,1))/dx
            dvdy=(pV(i,j+1,k,1)-pV(i,j,k,1))/dy
            dwdz=(pW(i,j,k+1,1)-pW(i,j,k,1))/dz
            ! Off-diagonal gradients
            dudy=0.25_WP*(pU(i,j+1,k,1)+pU(i+1,j+1,k,1)-pU(i,j-1,k,1)-pU(i+1,j-1,k,1))/(2.0_WP*dy)
            dudz=0.25_WP*(pU(i,j,k+1,1)+pU(i+1,j,k+1,1)-pU(i,j,k-1,1)-pU(i+1,j,k-1,1))/(2.0_WP*dz)
            dvdx=0.25_WP*(pV(i+1,j,k,1)+pV(i+1,j+1,k,1)-pV(i-1,j,k,1)-pV(i-1,j+1,k,1))/(2.0_WP*dx)
            dvdz=0.25_WP*(pV(i,j,k+1,1)+pV(i,j+1,k+1,1)-pV(i,j,k-1,1)-pV(i,j+1,k-1,1))/(2.0_WP*dz)
            dwdx=0.25_WP*(pW(i+1,j,k,1)+pW(i+1,j,k+1,1)-pW(i-1,j,k,1)-pW(i-1,j,k+1,1))/(2.0_WP*dx)
            dwdy=0.25_WP*(pW(i,j+1,k,1)+pW(i,j+1,k+1,1)-pW(i,j-1,k,1)-pW(i,j-1,k+1,1))/(2.0_WP*dy)
            ! |∇u| = sqrt(sum of all gradients squared)
            gradU_mag=sqrt(dudx**2+dudy**2+dudz**2+dvdx**2+dvdy**2+dvdz**2+dwdx**2+dwdy**2+dwdz**2)
            ! Normalize into a local Reynolds number
            Re_cell=solver%rho*gradU_mag*min(dx,dy,dz)**2/solver%visc
            ! Tagged based on cell Re value
            if (Re_cell.gt.Re_tag) tagarr(i,j,k,1)=SETtag
            ! Also tag based on closeness to sphere surface
            dist=sphere_levelset([solver%amr%xlo+(real(i,WP)+0.5_WP)*dx,solver%amr%ylo+(real(j,WP)+0.5_WP)*dy,solver%amr%zlo+(real(k,WP)+0.5_WP)*dz],time)
            if (dist.lt.5.0_WP*dx.and.dist.gt.-dx) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> Dirichlet BC: uniform inflow Uin at xlo/xhi for U
   subroutine dirichlet_velocity(solver,bx,p,comp,face,time,geom)
      use amrex_amr_module, only: amrex_box,amrex_geometry
      class(amrincomp), intent(in) :: solver
      type(amrex_box), intent(in) :: bx
      real(WP), dimension(:,:,:,:), pointer, intent(inout) :: p
      character(len=1), intent(in) :: comp
      integer, intent(in) :: face
      real(WP), intent(in) :: time
      type(amrex_geometry), intent(in) :: geom
      integer :: i,j,k
      ! face: 1=xlo, 2=xhi, 3=ylo, 4=yhi, 5=zlo, 6=zhi
      select case (comp)
       case ('U')
         select case (face)
          case (1)  ! xlo - Inflow Dirichlet condition
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=1.0_WP
            end do; end do; end do
         end select
      end select
   end subroutine dirichlet_velocity

   !> Initialize fluid volume fraction
   subroutine init_VF(data,lvl,time,ba,dm)
      use mms_geom, only: initialize_volume_moments
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
      class(amrdata), intent(inout) :: data
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      real(WP), dimension(3) :: BL,BG  ! Dummy barycenters
      real(WP) :: dx,dy,dz
      integer :: i,j,k
      real(WP), parameter :: VFlo=1.0e-12_WP
      dx=data%amr%dx(lvl); dy=data%amr%dy(lvl); dz=data%amr%dz(lvl)
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         bx=mfi%growntilebox(data%ng)
         pVF=>data%mf(lvl)%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            call initialize_volume_moments(lo=[data%amr%xlo+real(i  ,WP)*dx,data%amr%ylo+real(j  ,WP)*dy,data%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[data%amr%xlo+real(i+1,WP)*dx,data%amr%ylo+real(j+1,WP)*dy,data%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=sphere_levelset,time=time,level=4,VFlo=VFlo,VF=pVF(i,j,k,1),BL=BL,BG=BG)
            pVF(i,j,k,1)=max(pVF(i,j,k,1),VFlo)
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine init_VF

   !> Initialization hook
   subroutine simulation_init()
      use param, only: param_read
      implicit none
      
      ! Create amrgrid
      create_amrgrid: block
         amr%name='amrsphere'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=-5.0_WP; amr%xhi=+15.0_WP
         amr%ylo=-10.0_WP; amr%yhi=+10.0_WP
         amr%zlo=-10.0_WP; amr%zhi=+10.0_WP
         amr%xper=.false.; amr%yper=.true.; amr%zper=.true.
         call param_read('Max level',amr%maxlvl)
         call amr%initialize()
      end block create_amrgrid

      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      ! Create flow solver
      create_flow_solver: block
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap
         ! Create flow solver
         call fs%initialize(amr)
         ! Set viscosity
         call param_read('Reynolds number',fs%visc)
         fs%visc=1.0_WP/fs%visc
         ! Set pressure convergence
         fs%psolver%max_iter=20
         fs%psolver%tol_rel=1.0e-5_WP
         ! Set boundary conditions
         fs%U%lo_bc(1,1)=amrex_bc_ext_dir
         fs%V%lo_bc(1,1)=amrex_bc_foextrap
         fs%W%lo_bc(1,1)=amrex_bc_foextrap
         fs%U%hi_bc(1,1)=amrex_bc_foextrap
         fs%V%hi_bc(1,1)=amrex_bc_foextrap
         fs%W%hi_bc(1,1)=amrex_bc_foextrap
         fs%user_dirichlet=>dirichlet_velocity
      end block create_flow_solver

      ! Create workspace arrays
      create_workspace: block
         use amrdata_class, only: amrex_interp_none
         call resU%initialize(amr,name='resU',ncomp=1,ng=0,nodal=[.true. ,.false.,.false.],interp=amrex_interp_none); call resU%register()
         call resV%initialize(amr,name='resV',ncomp=1,ng=0,nodal=[.false.,.true. ,.false.],interp=amrex_interp_none); call resV%register()
         call resW%initialize(amr,name='resW',ncomp=1,ng=0,nodal=[.false.,.false.,.true. ],interp=amrex_interp_none); call resW%register()
      end block create_workspace

      ! Create VF for IB forcing
      create_VF: block
         use amrdata_class, only: amrex_interp_reinit
         call VF%initialize(amr,name='VF',ncomp=1,ng=1,interp=amrex_interp_reinit)
         VF%user_init=>init_VF
         call VF%register()
      end block create_VF

      ! Initialize regridding
      init_regridding: block
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_tagging=>my_tagger
         call param_read('Tagging Reynolds',Re_tag)
         ! Create initial grid
         call amr%init_from_scratch(time=time%t)
      end block init_regridding

      ! Initialize visualization
      create_visualization: block
         ! Create visualization object
         call viz%initialize(amr,'amrsphere')
         call viz%add_scalar(fs%U,1,'U')
         call viz%add_scalar(fs%V,1,'V')
         call viz%add_scalar(fs%W,1,'W')
         call viz%add_scalar(fs%P,1,'pressure')
         call viz%add_scalar(fs%div,1,'divergence')
         call viz%add_scalar(VF,1,'VF')
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) call viz%write(time=time%t)
      end block create_visualization

      ! Create monitor
      create_monitor: block
         ! Get solver info and cfl
         call fs%get_info()
         call fs%get_cfl(time%dt,time%cfl)
         ! Create simulation monitor
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'dt')
         call mfile%add_column(fs%CFL,'CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%divmax,'Divergence')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLc_x,'CFLc_x')
         call cflfile%add_column(fs%CFLc_y,'CFLc_y')
         call cflfile%add_column(fs%CFLc_z,'CFLc_z')
         call cflfile%add_column(fs%CFLv_x,'CFLv_x')
         call cflfile%add_column(fs%CFLv_y,'CFLv_y')
         call cflfile%add_column(fs%CFLv_z,'CFLv_z')
         call cflfile%write()
         ! Create grid monitor
         gridfile=monitor(amRoot=amr%amRoot,name='grid')
         call gridfile%add_column(time%n,'Timestep')
         call gridfile%add_column(time%t,'Time')
         call gridfile%add_column(amr%nlevels,'Nlvl')
         call gridfile%add_column(amr%nboxes,'Nbox')
         call gridfile%add_column(amr%ncells,'Ncell')
         call gridfile%add_column(amr%compression,'Compression')
         call gridfile%add_column(amr%maxRSS,'Maximum RSS')
         call gridfile%add_column(amr%minRSS,'Minimum RSS')
         call gridfile%add_column(amr%avgRSS,'Average RSS')
         call gridfile%write()
      end block create_monitor

   end subroutine simulation_init

   !> Run the simulation
   subroutine simulation_run()
      implicit none

      ! Time integration loop
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Store old velocity
         call fs%Uold%copy(src=fs%U)
         call fs%Vold%copy(src=fs%V)
         call fs%Wold%copy(src=fs%W)

         ! Sub-iterations
         do while (time%it.le.time%itmax)

            ! Build mid-time velocity: U^{mid} = 0.5*(U + Uold)
            call fs%U%lincomb(a=0.5_WP,src1=fs%U,b=0.5_WP,src2=fs%Uold)
            call fs%V%lincomb(a=0.5_WP,src1=fs%V,b=0.5_WP,src2=fs%Vold)
            call fs%W%lincomb(a=0.5_WP,src1=fs%W,b=0.5_WP,src2=fs%Wold)

            ! Compute advective momentum RHS
            call fs%get_dmomdt(U=fs%U,V=fs%V,W=fs%W,drhoUdt=resU,drhoVdt=resV,drhoWdt=resW,time=time%t)

            ! Increment velocity
            call fs%U%lincomb(a=1.0_WP,src1=fs%Uold,b=time%dt/fs%rho,src2=resU)
            call fs%V%lincomb(a=1.0_WP,src1=fs%Vold,b=time%dt/fs%rho,src2=resV)
            call fs%W%lincomb(a=1.0_WP,src1=fs%Wold,b=time%dt/fs%rho,src2=resW)

            ! Apply IB direct forcing
            ibforcing: block
               use amrex_amr_module, only: amrex_mfiter,amrex_box
               type(amrex_mfiter) :: mfi
               type(amrex_box) :: bx
               real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW,pVF
               integer :: i,j,k,lvl
               do lvl=0,amr%clvl()
                  call amr%mfiter_build(lvl,mfi)
                  do while (mfi%next())
                     ! Get data pointers
                     pU=>fs%U%mf(lvl)%dataptr(mfi)
                     pV=>fs%V%mf(lvl)%dataptr(mfi)
                     pW=>fs%W%mf(lvl)%dataptr(mfi)
                     pVF=>VF%mf(lvl)%dataptr(mfi)
                     ! U at x-faces
                     bx=mfi%nodaltilebox(1)
                     do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                        pU(i,j,k,1)=0.5_WP*sum(pVF(i-1:i,j,k,1))*pU(i,j,k,1)
                     end do; end do; end do
                     ! V at y-faces
                     bx=mfi%nodaltilebox(2)
                     do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                        pV(i,j,k,1)=0.5_WP*sum(pVF(i,j-1:j,k,1))*pV(i,j,k,1)
                     end do; end do; end do
                     ! W at z-faces
                     bx=mfi%nodaltilebox(3)
                     do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                        pW(i,j,k,1)=0.5_WP*sum(pVF(i,j,k-1:k,1))*pW(i,j,k,1)
                     end do; end do; end do
                  end do
                  call amr%mfiter_destroy(mfi)
               end do
            end block ibforcing

            ! Average down and fill ghosts
            call fs%average_down_velocity()
            call fs%fill_velocity(time%t)

            ! Correct outflow for mass conservation
            call fs%correct_outflow()

            ! Compute divergence
            call fs%get_div()

            ! Solve pressure Poisson
            call fs%div%mult(val=fs%rho/time%dt)
            call fs%psolver%solve(rhs=fs%div,phi=fs%P)

            ! Get gradients and correct velocity
            call fs%psolver%get_fluxes(phi=fs%P,flux_x=resU,flux_y=resV,flux_z=resW)
            call fs%U%saxpy(a=time%dt/fs%rho,src=resU)
            call fs%V%saxpy(a=time%dt/fs%rho,src=resV)
            call fs%W%saxpy(a=time%dt/fs%rho,src=resW)

            ! Average down and fill ghosts
            call fs%average_down_velocity()
            call fs%fill_velocity(time%t)

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Monitor output
         call fs%get_info()
         call mfile%write()
         call cflfile%write()

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time=time%t)
         
      end do

   end subroutine simulation_run

   !> Finalization hook
   subroutine simulation_final()
      implicit none
      ! Finalize time
      call time%finalize()
      ! Finalize grid
      call amr%finalize()
      call regrid_evt%finalize()
      ! Finalize solver
      call fs%finalize()
      call resU%finalize()
      call resV%finalize()
      call resW%finalize()
      call VF%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call gridfile%finalize()
   end subroutine simulation_final

end module simulation
