!> Test VOF advection with time loop
module mod_test_amrvof
   use precision,         only: WP
   use amrviz_class,      only: amrviz
   use amrgrid_class,     only: amrgrid
   use amrvof_class,      only: amrvof, VFlo
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use messager,          only: log
   use string,            only: itoa, rtoa
   use amrex_amr_module,  only: amrex_multifab, amrex_multifab_destroy
   implicit none
   private
   public :: test_amrvof

   ! Grid
   type(amrgrid), target :: amr

   ! VOF solver
   type(timetracker) :: time
   type(amrvof), target :: vof

   ! Velocity MultiFabs (finest level only, staggered)
   type(amrex_multifab) :: U, V, W
   integer :: vel_ng = 2  ! Ghost cells for velocity

   ! Sphere parameters
   real(WP) :: sphere_xc, sphere_yc, sphere_zc, sphere_radius

   ! Visualization
   type(amrviz) :: viz
   type(event) :: viz_evt

   ! Regrid parameters (disabled by default)
   type(event) :: regrid_evt

   ! Monitoring
   type(monitor) :: mfile, gridfile

contains

   !> Initialize VF field with sphere using levelset-based moments
   subroutine sphere_init(solver, lvl, time, ba, dm)
      use amrex_amr_module, only: amrex_mfiter, amrex_box, amrex_boxarray, amrex_distromap, &
      &                           amrex_mfiter_build, amrex_mfiter_destroy
      use mms_geom,         only: initialize_volume_moments
      class(amrvof), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas
      real(WP) :: dx, dy, dz
      real(WP), dimension(3) :: BL, BG
      integer :: i, j, k
      
      dx = solver%amr%dx(lvl)
      dy = solver%amr%dy(lvl)
      dz = solver%amr%dz(lvl)
      
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi, ba, dm, tiling=.false.)
      do while (mfi%next())
         bx = mfi%growntilebox(solver%VF%ng)  ! Include ghost cells
         pVF => solver%VF%mf(lvl)%dataptr(mfi)
         pCliq => solver%Cliq%mf(lvl)%dataptr(mfi)
         pCgas => solver%Cgas%mf(lvl)%dataptr(mfi)
         
         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)
                  ! Compute VF and barycenters from levelset with 4 levels of refinement
                  call initialize_volume_moments(lo=[solver%amr%xlo + real(i  ,WP)*dx, solver%amr%ylo + real(j  ,WP)*dy, solver%amr%zlo + real(k  ,WP)*dz], &
                  &                              hi=[solver%amr%xlo + real(i+1,WP)*dx, solver%amr%ylo + real(j+1,WP)*dy, solver%amr%zlo + real(k+1,WP)*dz], &
                  &                              levelset=sphere_levelset, time=time, level=4, VFlo=VFlo, VF=pVF(i,j,k,1), BL=BL, BG=BG)
                  ! Store barycenters
                  pCliq(i,j,k,1:3) = BL
                  pCgas(i,j,k,1:3) = BG
               end do
            end do
         end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine sphere_init
   
   !> Sphere levelset function with periodic distance
   function sphere_levelset(xyz, t) result(G)
      implicit none
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      real(WP), dimension(3) :: d, L
      L = [amr%xhi - amr%xlo, amr%yhi - amr%ylo, amr%zhi - amr%zlo]
      d = xyz - [sphere_xc, sphere_yc, sphere_zc]
      d = d - L * nint(d / L)  ! Nearest image
      G = sphere_radius - sqrt(sum(d**2))
   end function sphere_levelset


   !> Main test routine
   subroutine test_amrvof()
      implicit none

      call log("=== VOF Advection Test ===")

      ! Create amrgrid
      create_amrgrid: block
         use param, only: param_read
         amr%name = 'vof_advect'
         call param_read('Base nx', amr%nx)
         call param_read('Base ny', amr%ny)
         call param_read('Base nz', amr%nz)
         amr%xlo = 0.0_WP; amr%xhi = 1.0_WP
         amr%ylo = 0.0_WP; amr%yhi = 1.0_WP
         amr%zlo = 0.0_WP; amr%zhi = 1.0_WP
         amr%xper = .true.; amr%yper = .true.; amr%zper = .true.
         call param_read('Max level', amr%maxlvl)
         call amr%initialize()
      end block create_amrgrid

      ! Initialize time tracker
      initialize_timetracker: block
         use param, only: param_read
         time = timetracker(amRoot=amr%amRoot)
         call param_read('Max time', time%tmax)
         call param_read('Max dt', time%dtmax)
         call param_read('Max CFL', time%cflmax)
         time%dt = time%dtmax
         time%itmax = 1  ! No sub-iterations
      end block initialize_timetracker

      ! Setup sphere parameters
      setup_sphere: block
         sphere_radius = 0.15_WP
         sphere_xc = 0.35_WP
         sphere_yc = 0.35_WP
         sphere_zc = 0.35_WP
      end block setup_sphere

      ! Create VOF solver
      create_vof_solver: block
         vof%user_init => sphere_init
         call vof%initialize(amr, name='sphere_vof')
      end block create_vof_solver

      ! Initialize grid
      initialize: block
         call amr%init_from_scratch(time=time%t)
      end block initialize

      ! Build initial PLIC and reset moments for consistency
      build_plic_and_reset_moments: block
         call vof%build_plic(time%t)
         call vof%reset_moments()
         call log("  Initial PLIC constructed")
      end block build_plic_and_reset_moments

      ! Create visualization
      create_visualization: block
         use param, only: param_read
         call viz%initialize(amr, 'vof_advect')
         call viz%add_scalar(vof%VF, 1, 'VF')
         call viz%add_surfmesh(vof%smesh, 'plic')
         ! Create visualization output event
         viz_evt = event(time=time, name='Visualization output')
         call param_read('Output period', viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) call viz%write(time=time%t)
      end block create_visualization

      ! Regridding setup (disabled by default)
      regrid_setup: block
         use param, only: param_read
         regrid_evt = event(time=time, name='Regrid')
         call param_read('Regrid nsteps', regrid_evt%nper, default=0)
      end block regrid_setup

      ! Create monitor
      create_monitor: block
         ! Create VOF monitor
         call vof%get_info()
         mfile = monitor(amRoot=amr%amRoot, name='simulation')
         call mfile%add_column(time%n, 'Timestep')
         call mfile%add_column(time%t, 'Time')
         call mfile%add_column(time%dt, 'dt')
         call mfile%add_column(time%cfl, 'CFL')
         call mfile%add_column(vof%VFint, 'VFint')
         call mfile%add_column(vof%VFmin, 'VFmin')
         call mfile%add_column(vof%VFmax, 'VFmax')
         call mfile%write()
         ! Create grid monitor
         gridfile = monitor(amRoot=amr%amRoot, name='grid')
         call gridfile%add_column(time%n, 'Timestep')
         call gridfile%add_column(time%t, 'Time')
         call gridfile%add_column(amr%nlevels, 'Nlvl')
         call gridfile%add_column(amr%nboxes, 'Nbox')
         call gridfile%add_column(amr%ncells, 'Ncell')
         call gridfile%add_column(amr%compression, 'Compression')
         call gridfile%add_column(amr%maxRSS,'Maximum RSS')
         call gridfile%add_column(amr%minRSS,'Minimum RSS')
         call gridfile%add_column(amr%avgRSS,'Average RSS')
         call gridfile%write()
      end block create_monitor

      ! Time integration loop
      time_loop: do while (.not.time%done())

         ! Build velocity and set LeVeque vortex field
         set_velocity: block
            use amrex_amr_module, only: amrex_mfiter, amrex_box, amrex_multifab
            use mathtools, only: Pi, twoPi
            type(amrex_mfiter) :: mfi
            type(amrex_box) :: bx
            type(amrex_multifab) :: A  ! Cell-centered vector potential (3 components)
            real(WP), dimension(:,:,:,:), contiguous, pointer :: pU, pV, pW, pA
            real(WP) :: x, y, z, dx, dy, dz, dxi, dyi, dzi, T, tfac
            integer :: i, j, k
            T = 3.0_WP
            tfac = cos(Pi*time%t/T) / Pi
            dx = amr%dx(amr%clvl()); dxi=1.0_WP/dx
            dy = amr%dy(amr%clvl()); dyi=1.0_WP/dy
            dz = amr%dz(amr%clvl()); dzi=1.0_WP/dz
            ! Build velocity MultiFabs
            call amr%mfab_build(amr%clvl(), U, ncomp=1, nover=vel_ng, atface=[.true. , .false., .false.])
            call amr%mfab_build(amr%clvl(), V, ncomp=1, nover=vel_ng, atface=[.false., .true. , .false.])
            call amr%mfab_build(amr%clvl(), W, ncomp=1, nover=vel_ng, atface=[.false., .false., .true. ])
            ! Build cell-centered MultiFab for vector potential (Ax, Ay, Az)
            call amr%mfab_build(amr%clvl(), A, ncomp=3, nover=vel_ng+1)
            ! Fill vector potential at cell centers
            call amr%mfiter_build(amr%clvl(), mfi)
            do while (mfi%next())
               pA => A%dataptr(mfi)
               bx = mfi%growntilebox(vel_ng+1)
               do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
                  x = amr%xlo + (real(i,WP)+0.5_WP)*dx
                  y = amr%ylo + (real(j,WP)+0.5_WP)*dy
                  z = amr%zlo + (real(k,WP)+0.5_WP)*dz
                  pA(i,j,k,1) =  0.0_WP
                  pA(i,j,k,2) = -tfac * sin(Pi*x)**2 * sin(twoPi*y) * sin(Pi*z)**2
                  pA(i,j,k,3) = +tfac * sin(Pi*x)**2 * sin(Pi*y)**2 * sin(twoPi*z)
               end do; end do; end do
            end do
            call amr%mfiter_destroy(mfi)
            ! Compute velocity from curl of vector potential using 4-point averaging
            call amr%mfiter_build(amr%clvl(), mfi)
            do while (mfi%next())
               pU => U%dataptr(mfi); pV => V%dataptr(mfi); pW => W%dataptr(mfi); pA => A%dataptr(mfi)
               ! U = dAz/dy - dAy/dz at X-faces
               bx = mfi%grownnodaltilebox(1-1, vel_ng) ! bug amrex's fortran interface!
               do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
                  pU(i,j,k,1) = 0.25_WP*dyi*(pA(i-1,j+1,k,3)+pA(i,j+1,k,3)-pA(i-1,j-1,k,3)-pA(i,j-1,k,3)) &
                  &           - 0.25_WP*dzi*(pA(i-1,j,k+1,2)+pA(i,j,k+1,2)-pA(i-1,j,k-1,2)-pA(i,j,k-1,2))
               end do; end do; end do
               ! V = dAx/dz - dAz/dx at Y-faces
               bx = mfi%grownnodaltilebox(2-1, vel_ng) ! bug amrex's fortran interface!
               do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
                  pV(i,j,k,1) = 0.25_WP*dzi*(pA(i,j-1,k+1,1)+pA(i,j,k+1,1)-pA(i,j-1,k-1,1)-pA(i,j,k-1,1)) &
                  &           - 0.25_WP*dxi*(pA(i+1,j-1,k,3)+pA(i+1,j,k,3)-pA(i-1,j-1,k,3)-pA(i-1,j,k,3))
               end do; end do; end do
               ! W = dAy/dx - dAx/dy at Z-faces
               bx = mfi%grownnodaltilebox(3-1, vel_ng) ! bug amrex's fortran interface!
               do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
                  pW(i,j,k,1) = 0.25_WP*dxi*(pA(i+1,j,k-1,2)+pA(i+1,j,k,2)-pA(i-1,j,k-1,2)-pA(i-1,j,k,2)) &
                  &           - 0.25_WP*dyi*(pA(i,j+1,k-1,1)+pA(i,j+1,k,1)-pA(i,j-1,k-1,1)-pA(i,j-1,k,1))
               end do; end do; end do
            end do
            call amr%mfiter_destroy(mfi)
            call amrex_multifab_destroy(A)
         end block set_velocity

         ! Compute CFL and update dt based on CFL constraint
         call vof%get_cfl(U=U, V=V, W=W, dt=time%dt, cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Store old state
         call vof%VFold%copy(src=vof%VF)
         call vof%Cliqold%copy(src=vof%Cliq)
         call vof%Cgasold%copy(src=vof%Cgas)
         call vof%PLICold%copy(src=vof%PLIC)

         ! Advect VOF
         call vof%advance_vof(U=U, V=V, W=W, dt=time%dt, time=time%t)

         ! Destroy velocity
         call amrex_multifab_destroy(U)
         call amrex_multifab_destroy(V)
         call amrex_multifab_destroy(W)

         ! Rebuild PLIC and reset moments for consistency
         call vof%build_plic(time%t)
         call vof%reset_moments()

         ! Regrid if event triggers (disabled by default)
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0, time=time%t)
            call gridfile%write()
         end if

         ! Monitor output
         call vof%get_info()
         call mfile%write()

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time=time%t)
         
      end do time_loop

      ! Final summary
      final_summary: block
         call log("  Final VFint: "//trim(rtoa(vof%VFint)))
         call log("  Final VFmin: "//trim(rtoa(vof%VFmin)))
         call log("  Final VFmax: "//trim(rtoa(vof%VFmax)))
      end block final_summary

      ! Cleanup
      cleanup: block
         call viz%finalize()
         call vof%finalize()
         call amr%finalize()
         call time%finalize()
         call mfile%finalize()
         call viz_evt%finalize()
         call regrid_evt%finalize()
      end block cleanup

      call log("=== VOF Advection Test Complete ===")

   end subroutine test_amrvof

end module mod_test_amrvof
