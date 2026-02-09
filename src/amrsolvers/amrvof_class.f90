!> AMR VOF Solver class
!> Provides Volume-of-Fluid advection for two-phase flow with AMReX
!> IRL-free implementation using native cutting geometry
module amrvof_class
   use iso_c_binding,    only: c_ptr, c_null_ptr, c_loc, c_f_pointer, c_char
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata
   use amrsolver_class,  only: amrsolver
   use surfmesh_class,   only: surfmesh
   use amrex_amr_module, only: amrex_multifab, amrex_mfiter, amrex_box, &
   &                           amrex_boxarray, amrex_distromap, amrex_geometry
   implicit none
   private

   ! Expose type and dispatchers
   public :: amrvof
   public :: amrvof_on_init, amrvof_on_coarse, amrvof_on_remake
   public :: amrvof_on_clear, amrvof_tagging, amrvof_postregrid

   ! PLIC boundary condition types
   integer, parameter, public :: BC_LIQ     = 1  !< All liquid in ghost
   integer, parameter, public :: BC_GAS     = 2  !< All gas in ghost
   integer, parameter, public :: BC_REFLECT = 3  !< Symmetry (mirror across boundary)
   integer, parameter, public :: BC_USER    = 4  !< User-defined callback

   ! Default parameters for volume fraction
   real(WP), parameter, public :: VFlo = 1.0e-12_WP    !< Minimum VF value considered
   real(WP), parameter, public :: VFhi = 1.0_WP-VFlo   !< Maximum VF value considered
   real(WP), parameter, public :: vol_eps = 1.0e-8_WP  !< Volume epsilon for division by zero

   !> AMR VOF solver type
   type, extends(amrsolver) :: amrvof
      ! User-configurable callbacks
      procedure(vof_init_iface), pointer, nopass :: user_init => null()
      procedure(vof_tagging_iface), pointer, nopass :: user_tagging => null()
      procedure(vof_bc_iface), pointer, nopass :: user_vof_bc => null()

      ! PLIC boundary conditions (per face, only used if direction is non-periodic)
      integer :: vof_lo_bc(3) = BC_REFLECT
      integer :: vof_hi_bc(3) = BC_REFLECT

      ! VOF data (solver owns these - 4 MultiFabs as per plan)
      type(amrdata) :: VF           !< Volume fraction (cell-centered)
      type(amrdata) :: Cliq         !< Liquid barycenter (3 components)
      type(amrdata) :: Cgas         !< Gas barycenter (3 components)
      type(amrdata) :: PLIC         !< PLIC plane (4 components: nx, ny, nz, d)

      ! Old data for time stepping
      type(amrdata) :: VFold
      type(amrdata) :: Cliqold
      type(amrdata) :: Cgasold
      type(amrdata) :: PLICold

      ! Tagging parameter
      integer :: regrid_buffer = 10  !< Number of cells to buffer around interface for tagging

      ! Monitoring quantities
      real(WP) :: VFmin = 0.0_WP    !< Minimum VF
      real(WP) :: VFmax = 0.0_WP    !< Maximum VF
      real(WP) :: VFint = 0.0_WP    !< Integral of VF (liquid volume)
      
      ! Surface mesh for visualization (finest level polygons)
      type(surfmesh) :: smesh

   contains
      procedure :: initialize
      procedure :: finalize
      ! Override internal type-bound callbacks from amrsolver
      procedure :: on_init
      procedure :: on_coarse
      procedure :: on_remake
      procedure :: on_clear
      procedure :: post_regrid
      ! Deferred from amrsolver base class
      procedure :: get_info
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
      ! VOF-specific procedures
      procedure :: build_plic             !< Reconstruct PLIC from VF and barycenters
      procedure :: advance_vof            !< Advect VF using staggered or collocated velocity
      procedure :: fill_moments_lvl       !< Fill VF/Cliq/Cgas ghosts at level (sync + BC)
      procedure :: sync_moments_lvl       !< Sync VF/Cliq/Cgas ghosts at level + fix periodic barycenters
      procedure :: sync_moments           !< Sync VF/Cliq/Cgas ghosts all levels + fix periodic barycenters
      procedure :: sync_plic_lvl          !< Sync PLIC ghosts at level + fix periodic plane distance
      procedure :: sync_plic              !< Sync PLIC ghosts all levels + fix periodic plane distance
      procedure :: fill_plic_lvl          !< Fill PLIC ghosts at level (sync + BC)
      procedure :: average_down           !< Average down VF/Cliq/Cgas to coarse levels, and clean up PLIC at coarse levels
      procedure :: reset_moments          !< Recompute VF/barycenters from PLIC
      procedure :: get_cfl                !< Compute advective CFL at finest level
      procedure :: print => amrvof_print  !< Print solver info
   end type amrvof

   !> Abstract interface for user-overridable on_init callback
   abstract interface
      subroutine vof_init_iface(solver, lvl, time, ba, dm)
         import :: amrvof, WP, amrex_boxarray, amrex_distromap
         class(amrvof), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine vof_init_iface
   end interface

   !> Abstract interface for user-overridable tagging callback
   abstract interface
      subroutine vof_tagging_iface(solver, lvl, tags, time)
         import :: amrvof, c_ptr, WP
         class(amrvof), intent(inout) :: solver
         integer, intent(in) :: lvl
         type(c_ptr), intent(in) :: tags
         real(WP), intent(in) :: time
      end subroutine vof_tagging_iface
   end interface

   !> Abstract interface for user-defined VOF boundary condition
   !> User must set VF, Cliq, Cgas, and PLIC consistently in ghost cells
   !> what=1 -> fill PLIC, what=2 -> fill moments
   abstract interface
      subroutine vof_bc_iface(solver, bx, pVF, pCliq, pCgas, pPLIC, face, time, what)
         import :: amrvof, amrex_box, WP
         class(amrvof), intent(inout) :: solver
         type(amrex_box), intent(in) :: bx           !< Ghost region to fill
         real(WP), dimension(:,:,:,:), contiguous, intent(inout) :: pVF, pCliq, pCgas, pPLIC
         integer, intent(in) :: face                 !< 1=xlo,2=xhi,3=ylo,4=yhi,5=zlo,6=zhi
         real(WP), intent(in) :: time
         integer, intent(in) :: what                 !< 1=PLIC, 2=moments
      end subroutine vof_bc_iface
   end interface

contains

   ! ============================================================================
   ! DISPATCHERS (module-level) - recover concrete amrvof type
   ! ============================================================================

   !> Dispatch on_init: calls type-bound method then user callback
   subroutine amrvof_on_init(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrvof), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_init(lvl, time, ba, dm)
      if (associated(this%user_init)) call this%user_init(this, lvl, time, ba, dm)
   end subroutine amrvof_on_init

   !> Dispatch on_coarse: calls type-bound method
   subroutine amrvof_on_coarse(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrvof), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_coarse(lvl, time, ba, dm)
   end subroutine amrvof_on_coarse

   !> Dispatch on_remake: calls type-bound method
   subroutine amrvof_on_remake(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrvof), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_remake(lvl, time, ba, dm)
   end subroutine amrvof_on_remake

   !> Dispatch on_clear: calls type-bound method
   subroutine amrvof_on_clear(ctx, lvl)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrvof), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_clear(lvl)
   end subroutine amrvof_on_clear

   !> Dispatch tagging: tag cells near interface with regrid_buffer layer growth
   subroutine amrvof_tagging(ctx, lvl, tags, time)
      use amrex_amr_module, only: amrex_tagboxarray, amrex_mfiter, amrex_box, amrex_multifab
      use amrgrid_class, only: SETtag
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags
      real(WP), intent(in) :: time
      type(amrvof), pointer :: this
      type(amrex_tagboxarray) :: tba
      type(amrex_multifab) :: band
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pBand
      integer :: i, j, k, dir, n, layer
      integer, dimension(3) :: ind
      integer :: eff_buffer
   
      call c_f_pointer(ctx, this)
      tba = tags
   
      ! Build band MultiFab with 1 ghost cell
      call this%amr%mfab_build(lvl=lvl, mfab=band, ncomp=1, nover=1)
      call band%setval(0.0_WP)

      ! Effective buffer
      eff_buffer = max(1, this%regrid_buffer / (2**(this%amr%clvl() - lvl)))
   
      ! Pass 1: Mark interface cells (band=1)
      call this%amr%mfiter_build(lvl, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         pVF => this%VF%mf(lvl)%dataptr(mfi)
         pBand => band%dataptr(mfi)
         do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
            ! Flag all obvious mixture cells
            if (pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi) then
               pBand(i,j,k,1)=1.0_WP
               cycle
            end if
            ! We may have missed implicit interfaces, check those
            do dir=1,3; do n=-1,+1,2
               ind=[i,j,k]; ind(dir)=ind(dir)+n
               if (pVF(i,j,k,1).lt.VFlo.and.pVF(ind(1),ind(2),ind(3),1).gt.VFhi.or.&
               &   pVF(i,j,k,1).gt.VFhi.and.pVF(ind(1),ind(2),ind(3),1).lt.VFlo) then
                  pBand(i,j,k,1)=1.0_WP
                  cycle
               end if
            end do; end do
         end do; end do; end do
      end do
      call this%amr%mfiter_destroy(mfi)
      call band%fill_boundary(this%amr%geom(lvl))
   
      ! Pass 2: Grow band by effective buffer layers
      do layer = 2, eff_buffer
         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            pBand => band%dataptr(mfi)
            do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
               if (pBand(i,j,k,1).eq.0.0_WP.and.any(pBand(i-1:i+1,j-1:j+1,k-1:k+1,1).eq.real(layer-1,WP))) pBand(i,j,k,1)=real(layer,WP)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         call band%fill_boundary(this%amr%geom(lvl))
      end do
   
      ! Pass 3: Set tags from band
      call this%amr%mfiter_build(lvl, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         tagarr => tba%dataPtr(mfi)
         pBand => band%dataptr(mfi)
         do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
            if (pBand(i,j,k,1) .gt. 0.0_WP) tagarr(i,j,k,1) = SETtag
         end do; end do; end do
      end do
      call this%amr%mfiter_destroy(mfi)
   
      ! Cleanup
      call this%amr%mfab_destroy(band)
   
      ! Call user tagging if provided
      if (associated(this%user_tagging)) call this%user_tagging(this, lvl, tags, time)
   end subroutine amrvof_tagging

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrvof_postregrid(ctx, lbase, time)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrvof), pointer :: this
      call c_f_pointer(ctx, this)
      call this%post_regrid(lbase, time)
   end subroutine amrvof_postregrid

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize the VOF solver
   subroutine initialize(this, amr, name)
      class(amrvof), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name

      ! Set name
      if (present(name)) then
         this%name = trim(name)
      else
         this%name = 'UNNAMED_VOF'
      end if

      ! Store amrgrid pointer
      this%amr => amr

      ! Initialize VOF data (cell-centered)
      call this%VF%initialize(amr, name='VF', ncomp=1, ng=1)
      call this%Cliq%initialize(amr, name='Cliq', ncomp=3, ng=1)
      call this%Cgas%initialize(amr, name='Cgas', ncomp=3, ng=1)
      call this%PLIC%initialize(amr, name='PLIC', ncomp=4, ng=2)

      ! Initialize old data
      call this%VFold%initialize(amr, name='VFold', ncomp=1, ng=1)
      call this%Cliqold%initialize(amr, name='Cliqold', ncomp=3, ng=1)
      call this%Cgasold%initialize(amr, name='Cgasold', ncomp=3, ng=1)
      call this%PLICold%initialize(amr, name='PLICold', ncomp=4, ng=2)
      
      ! Initialize surface mesh for visualization
      this%smesh%name = trim(this%name)//'_plic'

      ! Set parent pointers for callback context access
      this%VF%parent => this
      this%Cliq%parent => this
      this%Cgas%parent => this
      this%PLIC%parent => this
      this%VFold%parent => this
      this%Cliqold%parent => this
      this%Cgasold%parent => this
      this%PLICold%parent => this

      ! Register all 6 callbacks with amrgrid using concrete dispatchers
      ! NOTE: Standard amrdata interpolation is fine IF the tagger ensures the
      !       interface never escapes the finest level between regrids. Tag with
      !       sufficient buffer: interface cells + CFL*dt*regrid_interval.
      select type (this)
       type is (amrvof)
         call this%amr%add_on_init   (amrvof_on_init,    c_loc(this))
         call this%amr%add_on_coarse (amrvof_on_coarse,  c_loc(this))
         call this%amr%add_on_remake (amrvof_on_remake,  c_loc(this))
         call this%amr%add_on_clear  (amrvof_on_clear,   c_loc(this))
         call this%amr%add_tagging   (amrvof_tagging,    c_loc(this))
         call this%amr%add_postregrid(amrvof_postregrid, c_loc(this))
      end select

      ! Print solver info
      call this%print()

   end subroutine initialize

   ! ============================================================================
   ! INTERNAL CALLBACK OVERRIDES
   ! ============================================================================

   !> Override on_init: reset levels and set to zero
   subroutine on_init(this, lvl, time, ba, dm)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Reset level layouts
      call this%VF%reset_level(lvl, ba, dm)
      call this%Cliq%reset_level(lvl, ba, dm)
      call this%Cgas%reset_level(lvl, ba, dm)
      call this%PLIC%reset_level(lvl, ba, dm)
      call this%VFold%reset_level(lvl, ba, dm)
      call this%Cliqold%reset_level(lvl, ba, dm)
      call this%Cgasold%reset_level(lvl, ba, dm)
      call this%PLICold%reset_level(lvl, ba, dm)
      ! Set to zero
      call this%VF%setval(val=0.0_WP, lvl=lvl)
      call this%Cliq%setval(val=0.0_WP, lvl=lvl)
      call this%Cgas%setval(val=0.0_WP, lvl=lvl)
      call this%PLIC%setval(val=0.0_WP, lvl=lvl)
      call this%VFold%setval(val=0.0_WP, lvl=lvl)
      call this%Cliqold%setval(val=0.0_WP, lvl=lvl)
      call this%Cgasold%setval(val=0.0_WP, lvl=lvl)
      call this%PLICold%setval(val=0.0_WP, lvl=lvl)
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse
   subroutine on_coarse(this, lvl, time, ba, dm)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Allocate and interpolate moments from coarse
      call this%VF%on_coarse(this%VF, lvl, time, ba, dm)
      call this%Cliq%on_coarse(this%Cliq, lvl, time, ba, dm)
      call this%Cgas%on_coarse(this%Cgas, lvl, time, ba, dm)
      ! PLIC: just allocate, then set to trivial planes (new cells are pure, away from interface)
      call this%PLIC%reset_level(lvl, ba, dm)
      call set_trivial_plic()
      ! Old data just needs geometry
      call this%VFold%reset_level(lvl, ba, dm)
      call this%Cliqold%reset_level(lvl, ba, dm)
      call this%Cgasold%reset_level(lvl, ba, dm)
      call this%PLICold%reset_level(lvl, ba, dm)
      ! Fill moment and PLIC ghosts
      call this%fill_moments_lvl(lvl, time)
      call this%fill_plic_lvl(lvl, time)
   contains
      !> Set PLIC to trivial planes based on VF
      subroutine set_trivial_plic()
         use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_box
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pPLIC
         integer :: i, j, k
         call amrex_mfiter_build(mfi, this%PLIC%mf(lvl), tiling=.false.)
         do while (mfi%next())
            bx = mfi%tilebox()
            pVF => this%VF%mf(lvl)%dataptr(mfi)
            pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     pPLIC(i,j,k,1:3) = 0.0_WP
                     pPLIC(i,j,k,4) = sign(1.0e10_WP, pVF(i,j,k,1) - 0.5_WP)
                  end do
               end do
            end do
         end do
         call amrex_mfiter_destroy(mfi)
      end subroutine set_trivial_plic
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid
   subroutine on_remake(this, lvl, time, ba, dm)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Remake all data (copies existing + fills from coarse for new areas)
      call this%VF%on_remake(this%VF, lvl, time, ba, dm)
      call this%Cliq%on_remake(this%Cliq, lvl, time, ba, dm)
      call this%Cgas%on_remake(this%Cgas, lvl, time, ba, dm)
      call this%PLIC%on_remake(this%PLIC, lvl, time, ba, dm)
      ! Fix PLIC for pure cells (new cells from coarse have interpolated PLIC which is wrong)
      call fix_pure_plic()
      ! Old data just needs geometry
      call this%VFold%reset_level(lvl, ba, dm)
      call this%Cliqold%reset_level(lvl, ba, dm)
      call this%Cgasold%reset_level(lvl, ba, dm)
      call this%PLICold%reset_level(lvl, ba, dm)
      ! Fill moment and PLIC ghosts
      call this%fill_moments_lvl(lvl, time)
      call this%fill_plic_lvl(lvl, time)
   contains
      !> Set PLIC to trivial for pure cells
      subroutine fix_pure_plic()
         use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_box
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pPLIC
         integer :: i, j, k
         real(WP) :: vf
         call amrex_mfiter_build(mfi, this%PLIC%mf(lvl), tiling=.false.)
         do while (mfi%next())
            bx = mfi%tilebox()
            pVF => this%VF%mf(lvl)%dataptr(mfi)
            pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     vf = pVF(i,j,k,1)
                     if (vf.lt.VFlo .or. vf.gt.VFhi) then
                        pPLIC(i,j,k,1:3) = 0.0_WP
                        pPLIC(i,j,k,4) = sign(1.0e10_WP, vf - 0.5_WP)
                     end if
                  end do
               end do
            end do
         end do
         call amrex_mfiter_destroy(mfi)
      end subroutine fix_pure_plic
   end subroutine on_remake

   !> Override on_clear: delete level
   subroutine on_clear(this, lvl)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%VF%clear_level(lvl)
      call this%Cliq%clear_level(lvl)
      call this%Cgas%clear_level(lvl)
      call this%PLIC%clear_level(lvl)
      call this%VFold%clear_level(lvl)
      call this%Cliqold%clear_level(lvl)
      call this%Cgasold%clear_level(lvl)
      call this%PLICold%clear_level(lvl)
   end subroutine on_clear

   !> Override post_regrid: average down for C/F consistency
   subroutine post_regrid(this, lbase, time)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      call this%average_down(lbase)
   end subroutine post_regrid

   !> Average down VF/Cliq/Cgas from finest to lbase, then sync ghost cells
   !> Clean up PLIC at coarse levels and sync ghost cells
   subroutine average_down(this, lbase)
      use amrex_interface, only: amrmfab_average_down_cell
      class(amrvof), intent(inout) :: this
      integer, intent(in), optional :: lbase
      integer :: lvl, lb
      lb = 0; if (present(lbase)) lb = lbase
      ! Average valid cells from fine to coarse
      do lvl = this%amr%clvl()-1, lb, -1
         call amrmfab_average_down_cell(fmf=this%VF%mf(lvl+1)  , cmf=this%VF%mf(lvl)  , rr=this%amr%rref(lvl), cgeom=this%amr%geom(lvl))
         call amrmfab_average_down_cell(fmf=this%Cliq%mf(lvl+1), cmf=this%Cliq%mf(lvl), rr=this%amr%rref(lvl), cgeom=this%amr%geom(lvl))
         call amrmfab_average_down_cell(fmf=this%Cgas%mf(lvl+1), cmf=this%Cgas%mf(lvl), rr=this%amr%rref(lvl), cgeom=this%amr%geom(lvl))
      end do
      ! Sync ghost cells on all levels + fix periodic barycenters
      call this%sync_moments()
      ! Clean up PLIC at coarse levels
      do lvl = this%amr%clvl()-1, lb, -1
         call set_trivial_plic()
         call this%sync_plic_lvl(lvl)
      end do
   contains
      !> Set PLIC to trivial planes based on VF
      subroutine set_trivial_plic()
         use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_box
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pPLIC
         integer :: i, j, k
         call amrex_mfiter_build(mfi, this%PLIC%mf(lvl), tiling=.false.)
         do while (mfi%next())
            bx = mfi%tilebox()
            pVF => this%VF%mf(lvl)%dataptr(mfi)
            pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     pPLIC(i,j,k,1:3) = 0.0_WP
                     pPLIC(i,j,k,4) = sign(1.0e10_WP, pVF(i,j,k,1) - 0.5_WP)
                  end do
               end do
            end do
         end do
         call amrex_mfiter_destroy(mfi)
      end subroutine set_trivial_plic
   end subroutine average_down

   !> Sync VF/Cliq/Cgas ghosts on all levels + fix periodic barycenters
   subroutine sync_moments(this)
      class(amrvof), intent(inout) :: this
      integer :: lvl
      do lvl = 0, this%amr%clvl()
         call this%sync_moments_lvl(lvl)
      end do
   end subroutine sync_moments

   !> Sync VF/Cliq/Cgas ghosts at level + fix periodic barycenters
   subroutine sync_moments_lvl(this, lvl)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_box, amrex_geometry
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter) :: mfi
      type(amrex_geometry) :: geom
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pCliq, pCgas
      integer :: ig, jg, kg, ilo, ihi, jlo, jhi, klo, khi
      integer :: dlo(3), dhi(3)
      real(WP) :: xL, yL, zL
      
      ! Sync ghosts (periodic + MPI exchange)
      call this%VF%sync_lvl(lvl)
      call this%Cliq%sync_lvl(lvl)
      call this%Cgas%sync_lvl(lvl)
      
      ! Fix barycenter positions in periodic ghost cells
      geom = this%amr%geom(lvl)
      dlo = geom%domain%lo
      dhi = geom%domain%hi
      xL = this%amr%xhi - this%amr%xlo
      yL = this%amr%yhi - this%amr%ylo
      zL = this%amr%zhi - this%amr%zlo
      
      call amrex_mfiter_build(mfi, this%Cliq%mf(lvl), tiling=.false.)
      do while (mfi%next())
         pCliq => this%Cliq%mf(lvl)%dataptr(mfi)
         pCgas => this%Cgas%mf(lvl)%dataptr(mfi)
         ilo = lbound(pCliq,1); ihi = ubound(pCliq,1)
         jlo = lbound(pCliq,2); jhi = ubound(pCliq,2)
         klo = lbound(pCliq,3); khi = ubound(pCliq,3)
         ! X-periodic
         if (this%amr%xper) then
            if (ilo .lt. dlo(1)) then
               do kg = klo, khi; do jg = jlo, jhi; do ig = ilo, dlo(1)-1
                  pCliq(ig,jg,kg,1) = pCliq(ig,jg,kg,1) - xL
                  pCgas(ig,jg,kg,1) = pCgas(ig,jg,kg,1) - xL
               end do; end do; end do
            end if
            if (ihi .gt. dhi(1)) then
               do kg = klo, khi; do jg = jlo, jhi; do ig = dhi(1)+1, ihi
                  pCliq(ig,jg,kg,1) = pCliq(ig,jg,kg,1) + xL
                  pCgas(ig,jg,kg,1) = pCgas(ig,jg,kg,1) + xL
               end do; end do; end do
            end if
         end if
         ! Y-periodic
         if (this%amr%yper) then
            if (jlo .lt. dlo(2)) then
               do kg = klo, khi; do jg = jlo, dlo(2)-1; do ig = ilo, ihi
                  pCliq(ig,jg,kg,2) = pCliq(ig,jg,kg,2) - yL
                  pCgas(ig,jg,kg,2) = pCgas(ig,jg,kg,2) - yL
               end do; end do; end do
            end if
            if (jhi .gt. dhi(2)) then
               do kg = klo, khi; do jg = dhi(2)+1, jhi; do ig = ilo, ihi
                  pCliq(ig,jg,kg,2) = pCliq(ig,jg,kg,2) + yL
                  pCgas(ig,jg,kg,2) = pCgas(ig,jg,kg,2) + yL
               end do; end do; end do
            end if
         end if
         ! Z-periodic
         if (this%amr%zper) then
            if (klo .lt. dlo(3)) then
               do kg = klo, dlo(3)-1; do jg = jlo, jhi; do ig = ilo, ihi
                  pCliq(ig,jg,kg,3) = pCliq(ig,jg,kg,3) - zL
                  pCgas(ig,jg,kg,3) = pCgas(ig,jg,kg,3) - zL
               end do; end do; end do
            end if
            if (khi .gt. dhi(3)) then
               do kg = dhi(3)+1, khi; do jg = jlo, jhi; do ig = ilo, ihi
                  pCliq(ig,jg,kg,3) = pCliq(ig,jg,kg,3) + zL
                  pCgas(ig,jg,kg,3) = pCgas(ig,jg,kg,3) + zL
               end do; end do; end do
            end if
         end if
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine sync_moments_lvl

   !> Sync PLIC ghosts at level + fix periodic plane distance
   !> The plane distance d must be corrected in periodic ghost cells:
   !> d ← d ± n·L where L is domain length and n is normal component
   subroutine sync_plic_lvl(this, lvl)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_geometry
      implicit none
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter) :: mfi
      type(amrex_geometry) :: geom
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pP
      integer :: ig, jg, kg, ilo, ihi, jlo, jhi, klo, khi
      integer :: dlo(3), dhi(3)
      real(WP) :: xL, yL, zL
      
      ! Sync ghosts (periodic + MPI exchange)
      call this%PLIC%sync_lvl(lvl)
      
      ! Get geometry and domain bounds
      geom = this%amr%geom(lvl)
      dlo = geom%domain%lo
      dhi = geom%domain%hi
      
      xL = this%amr%xhi - this%amr%xlo
      yL = this%amr%yhi - this%amr%ylo
      zL = this%amr%zhi - this%amr%zlo
      
      ! Fix periodic plane distance
      call amrex_mfiter_build(mfi, this%PLIC%mf(lvl), tiling=.false.)
      do while (mfi%next())
         pP => this%PLIC%mf(lvl)%dataptr(mfi)
         ilo = lbound(pP,1); ihi = ubound(pP,1)
         jlo = lbound(pP,2); jhi = ubound(pP,2)
         klo = lbound(pP,3); khi = ubound(pP,3)
         
         ! X-periodic
         if (this%amr%xper) then
            if (ilo .lt. dlo(1)) then
               do kg = klo, khi; do jg = jlo, jhi; do ig = ilo, dlo(1)-1
                  pP(ig,jg,kg,4) = pP(ig,jg,kg,4) - pP(ig,jg,kg,1)*xL
               end do; end do; end do
            end if
            if (ihi .gt. dhi(1)) then
               do kg = klo, khi; do jg = jlo, jhi; do ig = dhi(1)+1, ihi
                  pP(ig,jg,kg,4) = pP(ig,jg,kg,4) + pP(ig,jg,kg,1)*xL
               end do; end do; end do
            end if
         end if
         
         ! Y-periodic
         if (this%amr%yper) then
            if (jlo .lt. dlo(2)) then
               do kg = klo, khi; do jg = jlo, dlo(2)-1; do ig = ilo, ihi
                  pP(ig,jg,kg,4) = pP(ig,jg,kg,4) - pP(ig,jg,kg,2)*yL
               end do; end do; end do
            end if
            if (jhi .gt. dhi(2)) then
               do kg = klo, khi; do jg = dhi(2)+1, jhi; do ig = ilo, ihi
                  pP(ig,jg,kg,4) = pP(ig,jg,kg,4) + pP(ig,jg,kg,2)*yL
               end do; end do; end do
            end if
         end if
         
         ! Z-periodic
         if (this%amr%zper) then
            if (klo .lt. dlo(3)) then
               do kg = klo, dlo(3)-1; do jg = jlo, jhi; do ig = ilo, ihi
                  pP(ig,jg,kg,4) = pP(ig,jg,kg,4) - pP(ig,jg,kg,3)*zL
               end do; end do; end do
            end if
            if (khi .gt. dhi(3)) then
               do kg = dhi(3)+1, khi; do jg = jlo, jhi; do ig = ilo, ihi
                  pP(ig,jg,kg,4) = pP(ig,jg,kg,4) + pP(ig,jg,kg,3)*zL
               end do; end do; end do
            end if
         end if
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine sync_plic_lvl

   !> Sync PLIC ghosts on all levels + fix periodic plane distance
   subroutine sync_plic(this)
      class(amrvof), intent(inout) :: this
      integer :: lvl
      do lvl = 0, this%amr%clvl()
         call this%sync_plic_lvl(lvl)
      end do
   end subroutine sync_plic

   !> Fill PLIC ghosts at a level (fill + physical BC)
   !> Handles BC_LIQ (trivial d=+∞), BC_GAS (trivial d=-∞), 
   !> BC_REFLECT (mirror + flip normal), BC_USER (user callback)
   subroutine fill_plic_lvl(this, lvl, time)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_geometry, amrex_box
      implicit none
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_mfiter) :: mfi
      type(amrex_geometry) :: geom
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas, pPLIC
      integer :: ig, jg, kg, ilo, ihi, jlo, jhi, klo, khi
      integer :: dlo(3), dhi(3)
      real(WP) :: xL, yL, zL, dx, dy, dz
      
      ! Sync ghosts (periodic + MPI exchange + C/F)
      call this%PLIC%fill_lvl(lvl, time)
      
      ! Get geometry and domain bounds
      geom = this%amr%geom(lvl)
      dlo = geom%domain%lo
      dhi = geom%domain%hi
      xL = this%amr%xhi - this%amr%xlo
      yL = this%amr%yhi - this%amr%ylo
      zL = this%amr%zhi - this%amr%zlo
      dx = this%amr%dx(lvl)
      dy = this%amr%dy(lvl)
      dz = this%amr%dz(lvl)
      
      ! Fix periodic plane distance + apply physical BC
      call amrex_mfiter_build(mfi, this%PLIC%mf(lvl), tiling=.false.)
      do while (mfi%next())
         pVF   => this%VF%mf(lvl)%dataptr(mfi)
         pCliq => this%Cliq%mf(lvl)%dataptr(mfi)
         pCgas => this%Cgas%mf(lvl)%dataptr(mfi)
         pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
         ilo = lbound(pPLIC,1); ihi = ubound(pPLIC,1)
         jlo = lbound(pPLIC,2); jhi = ubound(pPLIC,2)
         klo = lbound(pPLIC,3); khi = ubound(pPLIC,3)
         
         ! X-periodic: fix plane distance
         if (this%amr%xper) then
            if (ilo .lt. dlo(1)) then
               do kg = klo, khi; do jg = jlo, jhi; do ig = ilo, dlo(1)-1
                  pPLIC(ig,jg,kg,4) = pPLIC(ig,jg,kg,4) - pPLIC(ig,jg,kg,1)*xL
               end do; end do; end do
            end if
            if (ihi .gt. dhi(1)) then
               do kg = klo, khi; do jg = jlo, jhi; do ig = dhi(1)+1, ihi
                  pPLIC(ig,jg,kg,4) = pPLIC(ig,jg,kg,4) + pPLIC(ig,jg,kg,1)*xL
               end do; end do; end do
            end if
         end if
         ! Y-periodic: fix plane distance
         if (this%amr%yper) then
            if (jlo .lt. dlo(2)) then
               do kg = klo, khi; do jg = jlo, dlo(2)-1; do ig = ilo, ihi
                  pPLIC(ig,jg,kg,4) = pPLIC(ig,jg,kg,4) - pPLIC(ig,jg,kg,2)*yL
               end do; end do; end do
            end if
            if (jhi .gt. dhi(2)) then
               do kg = klo, khi; do jg = dhi(2)+1, jhi; do ig = ilo, ihi
                  pPLIC(ig,jg,kg,4) = pPLIC(ig,jg,kg,4) + pPLIC(ig,jg,kg,2)*yL
               end do; end do; end do
            end if
         end if
         ! Z-periodic: fix plane distance
         if (this%amr%zper) then
            if (klo .lt. dlo(3)) then
               do kg = klo, dlo(3)-1; do jg = jlo, jhi; do ig = ilo, ihi
                  pPLIC(ig,jg,kg,4) = pPLIC(ig,jg,kg,4) - pPLIC(ig,jg,kg,3)*zL
               end do; end do; end do
            end if
            if (khi .gt. dhi(3)) then
               do kg = dhi(3)+1, khi; do jg = jlo, jhi; do ig = ilo, ihi
                  pPLIC(ig,jg,kg,4) = pPLIC(ig,jg,kg,4) + pPLIC(ig,jg,kg,3)*zL
               end do; end do; end do
            end if
         end if
         
         ! Apply physical BC for PLIC
         if (.not.this%amr%xper) then
            if (ilo.lt.dlo(1)) call apply_bc_face(1, -1, this%vof_lo_bc(1), ilo, dlo(1)-1, jlo, jhi, klo, khi, dlo(1), this%amr%xlo)
            if (ihi.gt.dhi(1)) call apply_bc_face(1, +1, this%vof_hi_bc(1), dhi(1)+1, ihi, jlo, jhi, klo, khi, dhi(1), this%amr%xhi)
         end if
         if (.not.this%amr%yper) then
            if (jlo.lt.dlo(2)) call apply_bc_face(2, -1, this%vof_lo_bc(2), ilo, ihi, jlo, dlo(2)-1, klo, khi, dlo(2), this%amr%ylo)
            if (jhi.gt.dhi(2)) call apply_bc_face(2, +1, this%vof_hi_bc(2), ilo, ihi, dhi(2)+1, jhi, klo, khi, dhi(2), this%amr%yhi)
         end if
         if (.not.this%amr%zper) then
            if (klo.lt.dlo(3)) call apply_bc_face(3, -1, this%vof_lo_bc(3), ilo, ihi, jlo, jhi, klo, dlo(3)-1, dlo(3), this%amr%zlo)
            if (khi.gt.dhi(3)) call apply_bc_face(3, +1, this%vof_hi_bc(3), ilo, ihi, jlo, jhi, dhi(3)+1, khi, dhi(3), this%amr%zhi)
         end if
      end do
      call amrex_mfiter_destroy(mfi)
      
   contains
      
      !> Apply BC to PLIC on a single face
      !> Uses host association for pVF, pCliq, pCgas, pPLIC
      subroutine apply_bc_face(dir, side, bc_type, i1, i2, j1, j2, k1, k2, bnd, x_bnd)
         integer, intent(in) :: dir, side, bc_type, i1, i2, j1, j2, k1, k2, bnd
         real(WP), intent(in) :: x_bnd
         integer :: ig, jg, kg, isrc, jsrc, ksrc, face
         type(amrex_box) :: bc_bx
         
         select case (bc_type)
         
          case (BC_LIQ)
            ! Trivial PLIC: full liquid
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               pPLIC(ig,jg,kg,1:3) = 0.0_WP
               pPLIC(ig,jg,kg,4) = 1.0e10_WP
            end do; end do; end do
            
          case (BC_GAS)
            ! Trivial PLIC: full gas
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               pPLIC(ig,jg,kg,1:3) = 0.0_WP
               pPLIC(ig,jg,kg,4) = -1.0e10_WP
            end do; end do; end do
            
          case (BC_REFLECT)
            ! Mirror PLIC from interior + flip normal component
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               isrc = ig; jsrc = jg; ksrc = kg
               if (dir.eq.1) isrc = 2*bnd - ig - side
               if (dir.eq.2) jsrc = 2*bnd - jg - side
               if (dir.eq.3) ksrc = 2*bnd - kg - side
               ! Copy plane
               pPLIC(ig,jg,kg,1:4) = pPLIC(isrc,jsrc,ksrc,1:4)
               ! Flip normal component
               pPLIC(ig,jg,kg,dir) = -pPLIC(ig,jg,kg,dir)
               ! Correct plane distance
               pPLIC(ig,jg,kg,4) = pPLIC(ig,jg,kg,4) - 2.0_WP*pPLIC(isrc,jsrc,ksrc,dir)*x_bnd
            end do; end do; end do
            
          case (BC_USER)
            ! User callback sets PLIC
            if (associated(this%user_vof_bc)) then
               bc_bx = amrex_box([i1, j1, k1], [i2, j2, k2])
               face = 2*dir - 1 + (1+side)/2  ! dir=1,side=-1 -> 1; dir=1,side=+1 -> 2; etc.
               call this%user_vof_bc(this, bc_bx, pVF, pCliq, pCgas, pPLIC, face, time, 1)
            end if
            
          case default
            ! Do nothing
            
         end select
         
      end subroutine apply_bc_face
      
   end subroutine fill_plic_lvl

   !> Fill moments ghosts at a level (fill + physical BC)
   !> Handles all BC types except mirroring/user PLIC
   subroutine fill_moments_lvl(this, lvl, time)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_box, amrex_geometry
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_mfiter) :: mfi
      type(amrex_geometry) :: geom
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas, pPLIC
      integer :: ig, jg, kg, ilo, ihi, jlo, jhi, klo, khi
      integer :: dlo(3), dhi(3)
      real(WP) :: xL, yL, zL, dx, dy, dz
      
      ! Sync ghosts (periodic + MPI exchange + C/F)
      call this%VF%fill_lvl(lvl, time)
      call this%Cliq%fill_lvl(lvl, time)
      call this%Cgas%fill_lvl(lvl, time)
      
      ! Get geometry info
      geom = this%amr%geom(lvl)
      dlo = geom%domain%lo
      dhi = geom%domain%hi
      xL = this%amr%xhi - this%amr%xlo
      yL = this%amr%yhi - this%amr%ylo
      zL = this%amr%zhi - this%amr%zlo
      dx = this%amr%dx(lvl)
      dy = this%amr%dy(lvl)
      dz = this%amr%dz(lvl)
      
      ! Fix barycenter positions in periodic ghost cells + apply physical BC
      call amrex_mfiter_build(mfi, this%VF%mf(lvl), tiling=.false.)
      do while (mfi%next())
         pVF   => this%VF%mf(lvl)%dataptr(mfi)
         pCliq => this%Cliq%mf(lvl)%dataptr(mfi)
         pCgas => this%Cgas%mf(lvl)%dataptr(mfi)
         pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
         ilo = lbound(pVF,1); ihi = ubound(pVF,1)
         jlo = lbound(pVF,2); jhi = ubound(pVF,2)
         klo = lbound(pVF,3); khi = ubound(pVF,3)
         
         ! X-periodic: shift barycenters
         if (this%amr%xper) then
            if (ilo .lt. dlo(1)) then
               do kg = klo, khi; do jg = jlo, jhi; do ig = ilo, dlo(1)-1
                  pCliq(ig,jg,kg,1) = pCliq(ig,jg,kg,1) - xL
                  pCgas(ig,jg,kg,1) = pCgas(ig,jg,kg,1) - xL
               end do; end do; end do
            end if
            if (ihi .gt. dhi(1)) then
               do kg = klo, khi; do jg = jlo, jhi; do ig = dhi(1)+1, ihi
                  pCliq(ig,jg,kg,1) = pCliq(ig,jg,kg,1) + xL
                  pCgas(ig,jg,kg,1) = pCgas(ig,jg,kg,1) + xL
               end do; end do; end do
            end if
         end if
         ! Y-periodic: shift barycenters
         if (this%amr%yper) then
            if (jlo .lt. dlo(2)) then
               do kg = klo, khi; do jg = jlo, dlo(2)-1; do ig = ilo, ihi
                  pCliq(ig,jg,kg,2) = pCliq(ig,jg,kg,2) - yL
                  pCgas(ig,jg,kg,2) = pCgas(ig,jg,kg,2) - yL
               end do; end do; end do
            end if
            if (jhi .gt. dhi(2)) then
               do kg = klo, khi; do jg = dhi(2)+1, jhi; do ig = ilo, ihi
                  pCliq(ig,jg,kg,2) = pCliq(ig,jg,kg,2) + yL
                  pCgas(ig,jg,kg,2) = pCgas(ig,jg,kg,2) + yL
               end do; end do; end do
            end if
         end if
         ! Z-periodic: shift barycenters
         if (this%amr%zper) then
            if (klo .lt. dlo(3)) then
               do kg = klo, dlo(3)-1; do jg = jlo, jhi; do ig = ilo, ihi
                  pCliq(ig,jg,kg,3) = pCliq(ig,jg,kg,3) - zL
                  pCgas(ig,jg,kg,3) = pCgas(ig,jg,kg,3) - zL
               end do; end do; end do
            end if
            if (khi .gt. dhi(3)) then
               do kg = dhi(3)+1, khi; do jg = jlo, jhi; do ig = ilo, ihi
                  pCliq(ig,jg,kg,3) = pCliq(ig,jg,kg,3) + zL
                  pCgas(ig,jg,kg,3) = pCgas(ig,jg,kg,3) + zL
               end do; end do; end do
            end if
         end if
         
         ! Apply physical BC for moments
         if (.not.this%amr%xper) then
            if (ilo.lt.dlo(1)) call apply_bc_face(1, -1, this%vof_lo_bc(1), ilo, dlo(1)-1, jlo, jhi, klo, khi, dlo(1), this%amr%xlo)
            if (ihi.gt.dhi(1)) call apply_bc_face(1, +1, this%vof_hi_bc(1), dhi(1)+1, ihi, jlo, jhi, klo, khi, dhi(1), this%amr%xhi)
         end if
         if (.not.this%amr%yper) then
            if (jlo.lt.dlo(2)) call apply_bc_face(2, -1, this%vof_lo_bc(2), ilo, ihi, jlo, dlo(2)-1, klo, khi, dlo(2), this%amr%ylo)
            if (jhi.gt.dhi(2)) call apply_bc_face(2, +1, this%vof_hi_bc(2), ilo, ihi, dhi(2)+1, jhi, klo, khi, dhi(2), this%amr%yhi)
         end if
         if (.not.this%amr%zper) then
            if (klo.lt.dlo(3)) call apply_bc_face(3, -1, this%vof_lo_bc(3), ilo, ihi, jlo, jhi, klo, dlo(3)-1, dlo(3), this%amr%zlo)
            if (khi.gt.dhi(3)) call apply_bc_face(3, +1, this%vof_hi_bc(3), ilo, ihi, jlo, jhi, dhi(3)+1, khi, dhi(3), this%amr%zhi)
         end if
      end do
      call amrex_mfiter_destroy(mfi)
      
   contains
      
      !> Apply BC to VF/Cliq/Cgas on a single face
      !> Uses host association for pVF, pCliq, pCgas, pPLIC, dx, dy, dz
      subroutine apply_bc_face(dir, side, bc_type, i1, i2, j1, j2, k1, k2, bnd, x_bnd)
         integer, intent(in) :: dir, side, bc_type, i1, i2, j1, j2, k1, k2, bnd
         real(WP), intent(in) :: x_bnd
         integer :: ig, jg, kg, isrc, jsrc, ksrc, face
         real(WP), dimension(3) :: center
         type(amrex_box) :: bc_bx
         
         select case (bc_type)
         
          case (BC_LIQ)
            ! Full liquid: VF=1, barycenters at cell center
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               center = [this%amr%xlo + (real(ig,WP)+0.5_WP)*dx, &
               &         this%amr%ylo + (real(jg,WP)+0.5_WP)*dy, &
               &         this%amr%zlo + (real(kg,WP)+0.5_WP)*dz]
               pVF(ig,jg,kg,1) = 1.0_WP
               pCliq(ig,jg,kg,1:3) = center
               pCgas(ig,jg,kg,1:3) = center
            end do; end do; end do
            
          case (BC_GAS)
            ! Full gas: VF=0, barycenters at cell center
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               center = [this%amr%xlo + (real(ig,WP)+0.5_WP)*dx, &
               &         this%amr%ylo + (real(jg,WP)+0.5_WP)*dy, &
               &         this%amr%zlo + (real(kg,WP)+0.5_WP)*dz]
               pVF(ig,jg,kg,1) = 0.0_WP
               pCliq(ig,jg,kg,1:3) = center
               pCgas(ig,jg,kg,1:3) = center
            end do; end do; end do
            
          case (BC_REFLECT)
            ! Mirror VF/Cliq/Cgas from interior + reflect barycenters
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               isrc = ig; jsrc = jg; ksrc = kg
               if (dir.eq.1) isrc = 2*bnd - ig - side
               if (dir.eq.2) jsrc = 2*bnd - jg - side
               if (dir.eq.3) ksrc = 2*bnd - kg - side
               ! Copy VF
               pVF(ig,jg,kg,1) = pVF(isrc,jsrc,ksrc,1)
               ! Copy and reflect barycenters
               pCliq(ig,jg,kg,1:3) = pCliq(isrc,jsrc,ksrc,1:3)
               pCgas(ig,jg,kg,1:3) = pCgas(isrc,jsrc,ksrc,1:3)
               pCliq(ig,jg,kg,dir) = 2.0_WP*x_bnd - pCliq(isrc,jsrc,ksrc,dir)
               pCgas(ig,jg,kg,dir) = 2.0_WP*x_bnd - pCgas(isrc,jsrc,ksrc,dir)
            end do; end do; end do
            
          case (BC_USER)
            ! User callback sets moments
            if (associated(this%user_vof_bc)) then
               bc_bx = amrex_box([i1, j1, k1], [i2, j2, k2])
               face = 2*dir - 1 + (1+side)/2  ! dir=1,side=-1 -> 1; dir=1,side=+1 -> 2; etc.
               call this%user_vof_bc(this, bc_bx, pVF, pCliq, pCgas, pPLIC, face, time, 2)
            end if
            
          case default
            ! Do nothing
            
         end select
         
      end subroutine apply_bc_face
      
   end subroutine fill_moments_lvl

   !> Finalize the VOF solver
   subroutine finalize(this)
      class(amrvof), intent(inout) :: this
      call this%VF%finalize()
      call this%Cliq%finalize()
      call this%Cgas%finalize()
      call this%PLIC%finalize()
      call this%VFold%finalize()
      call this%Cliqold%finalize()
      call this%Cgasold%finalize()
      call this%PLICold%finalize()
      call this%smesh%reset()
      nullify(this%amr)
      nullify(this%user_init)
      nullify(this%user_tagging)
   end subroutine finalize

   ! ============================================================================
   ! VOF-SPECIFIC METHODS (STUBS)
   ! ============================================================================

   !> Build PLIC reconstruction from VF and barycenters using PLICnet
   !> Also extracts PLIC polygons and accumulates them to smesh
   subroutine build_plic(this, time)
      use plicnet, only: get_normal, reflect_moments
      use mathtools, only: normalize
      use amrvof_geometry, only: get_plane_dist, cut_hex_polygon
      class(amrvof), intent(inout) :: this
      real(WP), intent(in) :: time
      integer :: lvl
      real(WP) :: dx, dy, dz, dxi, dyi, dzi

      ! Only build at finest level
      lvl = this%amr%clvl()

      ! Get cell size at this level
      dx = this%amr%dx(lvl); dxi = 1.0_WP / dx
      dy = this%amr%dy(lvl); dyi = 1.0_WP / dy
      dz = this%amr%dz(lvl); dzi = 1.0_WP / dz

      ! ========== Pass 1: Compute PLIC planes ==========
      plic_reconstruction: block
         integer :: i, j, k, ii, jj, kk, direction, direction2
         real(WP), dimension(0:188) :: moments
         real(WP), dimension(3) :: normal, center, lo, hi
         real(WP) :: m000, m100, m010, m001, temp, vf_cell
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas, pPLIC
         logical :: flip
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx

         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()

            ! Get pointers (with ghost cells for stencil access)
            pVF   => this%VF%mf(lvl)%dataptr(mfi)
            pCliq => this%Cliq%mf(lvl)%dataptr(mfi)
            pCgas => this%Cgas%mf(lvl)%dataptr(mfi)
            pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)

            ! Loop over cells in this box
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)

                     vf_cell = pVF(i,j,k,1)

                     ! Handle full cells: set trivial plane
                     if (vf_cell.lt.VFlo .or. vf_cell.gt.VFhi) then
                        pPLIC(i,j,k,1) = 0.0_WP  ! nx
                        pPLIC(i,j,k,2) = 0.0_WP  ! ny
                        pPLIC(i,j,k,3) = 0.0_WP  ! nz
                        pPLIC(i,j,k,4) = sign(1.0e10_WP, vf_cell - 0.5_WP)  ! d
                        cycle
                     end if

                     ! Liquid-gas symmetry
                     flip = .false.
                     if (vf_cell.ge.0.5_WP) flip = .true.

                     ! Initialize geometric moments
                     m000 = 0.0_WP; m100 = 0.0_WP; m010 = 0.0_WP; m001 = 0.0_WP

                     ! Construct neighborhood of volume moments (3x3x3 stencil)
                     if (flip) then
                        do kk = k-1, k+1
                           do jj = j-1, j+1
                              do ii = i-1, i+1
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+0) = 1.0_WP - pVF(ii,jj,kk,1)
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1) = (pCgas(ii,jj,kk,1) - (this%amr%xlo + (real(ii,WP)+0.5_WP)*dx)) * dxi
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2) = (pCgas(ii,jj,kk,2) - (this%amr%ylo + (real(jj,WP)+0.5_WP)*dy)) * dyi
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3) = (pCgas(ii,jj,kk,3) - (this%amr%zlo + (real(kk,WP)+0.5_WP)*dz)) * dzi
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4) = (pCliq(ii,jj,kk,1) - (this%amr%xlo + (real(ii,WP)+0.5_WP)*dx)) * dxi
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5) = (pCliq(ii,jj,kk,2) - (this%amr%ylo + (real(jj,WP)+0.5_WP)*dy)) * dyi
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6) = (pCliq(ii,jj,kk,3) - (this%amr%zlo + (real(kk,WP)+0.5_WP)*dz)) * dzi
                                 m000 = m000 +  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                                 m100 = m100 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                                 m010 = m010 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                                 m001 = m001 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                              end do
                           end do
                        end do
                     else
                        do kk = k-1, k+1
                           do jj = j-1, j+1
                              do ii = i-1, i+1
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+0) = pVF(ii,jj,kk,1)
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1) = (pCliq(ii,jj,kk,1) - (this%amr%xlo + (real(ii,WP)+0.5_WP)*dx)) * dxi
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2) = (pCliq(ii,jj,kk,2) - (this%amr%ylo + (real(jj,WP)+0.5_WP)*dy)) * dyi
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3) = (pCliq(ii,jj,kk,3) - (this%amr%zlo + (real(kk,WP)+0.5_WP)*dz)) * dzi
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4) = (pCgas(ii,jj,kk,1) - (this%amr%xlo + (real(ii,WP)+0.5_WP)*dx)) * dxi
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5) = (pCgas(ii,jj,kk,2) - (this%amr%ylo + (real(jj,WP)+0.5_WP)*dy)) * dyi
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6) = (pCgas(ii,jj,kk,3) - (this%amr%zlo + (real(kk,WP)+0.5_WP)*dz)) * dzi
                                 m000 = m000 +  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                                 m100 = m100 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                                 m010 = m010 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                                 m001 = m001 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                              end do
                           end do
                        end do
                     end if

                     ! Geometric center of neighborhood
                     if (m000.gt.tiny(1.0_WP)) then
                        center = [m100, m010, m001] / m000
                     else
                        center = 0.0_WP
                     end if

                     ! Apply symmetry (48 symmetries via reflect_moments)
                     call reflect_moments(moments, center, direction, direction2)

                     ! Get normal from neural network
                     call get_normal(moments, normal)
                     normal = normalize(normal)
                  
                     ! Undo direction2 rotation (axis permutation)
                     if (direction2.eq.1) then
                        temp = normal(1); normal(1) = normal(2); normal(2) = temp
                     else if (direction2.eq.2) then
                        temp = normal(2); normal(2) = normal(3); normal(3) = temp
                     else if (direction2.eq.3) then
                        temp = normal(1); normal(1) = normal(3); normal(3) = temp
                     else if (direction2.eq.4) then
                        temp = normal(2); normal(2) = normal(3); normal(3) = temp
                        temp = normal(1); normal(1) = normal(2); normal(2) = temp
                     else if (direction2.eq.5) then
                        temp = normal(1); normal(1) = normal(3); normal(3) = temp
                        temp = normal(1); normal(1) = normal(2); normal(2) = temp
                     end if
                  
                     ! Undo direction reflection (octant)
                     if (direction.eq.1) then
                        normal(1) = -normal(1)
                     else if (direction.eq.2) then
                        normal(2) = -normal(2)
                     else if (direction.eq.3) then
                        normal(3) = -normal(3)
                     else if (direction.eq.4) then
                        normal(1) = -normal(1); normal(2) = -normal(2)
                     else if (direction.eq.5) then
                        normal(1) = -normal(1); normal(3) = -normal(3)
                     else if (direction.eq.6) then
                        normal(2) = -normal(2); normal(3) = -normal(3)
                     else if (direction.eq.7) then
                        normal(1) = -normal(1); normal(2) = -normal(2); normal(3) = -normal(3)
                     end if
                  
                     ! Undo liquid-gas flip
                     if (.not.flip) normal = -normal
                  
                     ! Renormalize
                     normal = normalize(normal)
                  
                     ! Cell bounds
                     lo = [this%amr%xlo + real(i  ,WP)*dx, this%amr%ylo + real(j  ,WP)*dy, this%amr%zlo + real(k  ,WP)*dz]
                     hi = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
                  
                     ! Store PLIC plane: (nx, ny, nz, d)
                     pPLIC(i,j,k,1) = normal(1)
                     pPLIC(i,j,k,2) = normal(2)
                     pPLIC(i,j,k,3) = normal(3)
                     pPLIC(i,j,k,4) = get_plane_dist(normal, lo, hi, vf_cell)
                  
                  end do
               end do
            end do
         
         end do
         call this%amr%mfiter_destroy(mfi)
      end block plic_reconstruction
      
      ! Fill PLIC ghosts (sync + periodic correction + physical BC)
      call this%fill_plic_lvl(lvl, time)
      
      ! ========== Pass 2: Per-FAB polygon extraction ==========
      call this%smesh%reset()
      
      polygon_extraction: block
         integer :: i, j, k, ii, jj, kk
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pPLIC
         real(WP), dimension(3) :: lo, hi
         real(WP), dimension(4) :: plane
         real(WP), dimension(3,8) :: hex
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx, gbx
         
         ! Per-FAB polygon storage (allocatable, indexed by cell)
         real(WP), dimension(:,:,:,:,:), allocatable :: polygon_local  ! (3, 6, ilo:ihi, jlo:jhi, klo:khi)
         integer, dimension(:,:,:), allocatable :: poly_nv_local       ! (ilo:ihi, jlo:jhi, klo:khi)
         real(WP), dimension(3,6) :: poly_verts
         integer :: poly_nv
         integer :: glo(3), ghi(3)
         
         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            gbx = mfi%growntilebox(2)  ! Grown by 2 for 5x5x5 stencil
            pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
            
            ! Get bounds for per-FAB allocation
            glo = gbx%lo
            ghi = gbx%hi
            
            ! ----- Step A: Allocate per-FAB polygon storage -----
            allocate(polygon_local(1:3, 1:6, glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3)))
            allocate(poly_nv_local(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3)))
            polygon_local = 0.0_WP
            poly_nv_local = 0
            
            ! ----- Step B: Extract polygons (grown box including ghosts) -----
            do kk = glo(3), ghi(3)
               do jj = glo(2), ghi(2)
                  do ii = glo(1), ghi(1)
                     
                     ! Skip cells with no interface
                     if (abs(pPLIC(ii,jj,kk,4)) .gt. 1.0e+9_WP) cycle
                     
                     ! Build hex and plane
                     lo = [this%amr%xlo + real(ii  ,WP)*dx, this%amr%ylo + real(jj  ,WP)*dy, this%amr%zlo + real(kk  ,WP)*dz]
                     hi = [this%amr%xlo + real(ii+1,WP)*dx, this%amr%ylo + real(jj+1,WP)*dy, this%amr%zlo + real(kk+1,WP)*dz]
                     plane = [pPLIC(ii,jj,kk,1), pPLIC(ii,jj,kk,2), pPLIC(ii,jj,kk,3), pPLIC(ii,jj,kk,4)]
                     hex(:,1) = [hi(1), lo(2), lo(3)]
                     hex(:,2) = [hi(1), hi(2), lo(3)]
                     hex(:,3) = [hi(1), hi(2), hi(3)]
                     hex(:,4) = [hi(1), lo(2), hi(3)]
                     hex(:,5) = [lo(1), lo(2), lo(3)]
                     hex(:,6) = [lo(1), hi(2), lo(3)]
                     hex(:,7) = [lo(1), hi(2), hi(3)]
                     hex(:,8) = [lo(1), lo(2), hi(3)]
                     
                     call cut_hex_polygon(hex, plane, poly_nv, poly_verts)
                     
                     ! Store in per-FAB array
                     poly_nv_local(ii,jj,kk) = poly_nv
                     if (poly_nv.ge.3) then
                        polygon_local(:,1:poly_nv,ii,jj,kk) = poly_verts(:,1:poly_nv)
                     end if
                  end do
               end do
            end do
            
            ! ----- Step C: Compute curvature (valid cells, stencil access) -----
            ! TODO: curvature = f(polygon_local stencil around i,j,k)
            ! For now, skip curvature computation
            
            ! ----- Step D: Append to smesh (valid cells only) -----
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     poly_nv = poly_nv_local(i,j,k)
                     if (poly_nv.ge.3) then
                        call this%smesh%add_polygon(polygon_local(:,1:poly_nv,i,j,k), poly_nv)
                     end if
                  end do
               end do
            end do
            
            ! ----- Step E: Deallocate per-FAB storage -----
            deallocate(polygon_local, poly_nv_local)
            
         end do
         call this%amr%mfiter_destroy(mfi)
      end block polygon_extraction
      
   end subroutine build_plic
   
   !> Reset VF and barycenters from PLIC plane to ensure consistency
   !> Computes in valid + ghost cells from PLIC (which is already filled)
   !> Then averages down to coarse levels
   subroutine reset_moments(this)
      use amrvof_geometry, only: cut_hex_vol
      class(amrvof), intent(inout) :: this
      integer :: lvl,i,j,k
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCliq,pCgas,pPLIC
      real(WP), dimension(3,8) :: hex
      real(WP), dimension(4) :: plane
      real(WP) :: vol_liq,vol_gas,cell_vol,dx,dy,dz
      real(WP), dimension(3) :: bary_liq,bary_gas,cell_center
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      
      ! Only work at finest level
      lvl=this%amr%clvl()
      dx=this%amr%dx(lvl); dy=this%amr%dy(lvl); dz=this%amr%dz(lvl)
      cell_vol=dx*dy*dz
      
      call this%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         pVF  =>this%VF%mf(lvl)%dataptr(mfi)
         pCliq=>this%Cliq%mf(lvl)%dataptr(mfi)
         pCgas=>this%Cgas%mf(lvl)%dataptr(mfi)
         pPLIC=>this%PLIC%mf(lvl)%dataptr(mfi)
         ! Loop over tiles grown by 1 (matching VF and Cliq/Cgas)
         bx=mfi%growntilebox(1)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Cell center
            cell_center=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx, &
            &            this%amr%ylo+(real(j,WP)+0.5_WP)*dy, &
            &            this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
            ! Build hex cell (8 vertices)
            hex(:,1)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hex(:,2)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hex(:,3)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hex(:,4)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hex(:,5)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            hex(:,6)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            hex(:,7)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            hex(:,8)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            ! Get plane from PLIC
            plane(1:3)=pPLIC(i,j,k,1:3)
            plane(4)  =pPLIC(i,j,k,4)
            ! Skip cutting for full cells (trivial PLIC with large distance)
            if (abs(plane(4)).ge.1.0e9_WP) then
               if (plane(4).gt.0.0_WP) then
                  pVF(i,j,k,1)=1.0_WP
               else
                  pVF(i,j,k,1)=0.0_WP
               end if
               pCliq(i,j,k,1:3)=cell_center
               pCgas(i,j,k,1:3)=cell_center
               cycle
            end if
            ! Cut hex by plane
            call cut_hex_vol(hex,plane,vol_liq,vol_gas,bary_liq,bary_gas)
            ! Update VF and barycenters
            pVF(i,j,k,1)=vol_liq/cell_vol
            pCliq(i,j,k,1:3)=bary_liq
            pCgas(i,j,k,1:3)=bary_gas
            ! Clean up edge cases
            if (pVF(i,j,k,1).lt.VFlo) then
               pVF(i,j,k,1)=0.0_WP
               pCliq(i,j,k,1:3)=cell_center
               pCgas(i,j,k,1:3)=cell_center
            end if
            if (pVF(i,j,k,1).gt.VFhi) then
               pVF(i,j,k,1)=1.0_WP
               pCliq(i,j,k,1:3)=cell_center
               pCgas(i,j,k,1:3)=cell_center
            end if
         end do; end do; end do
      end do
      call this%amr%mfiter_destroy(mfi)
      
      ! Average down to coarse levels (uses ghost-capable average_down)
      call this%average_down()
      
   end subroutine reset_moments

   !> Advect VF using velocity field (staggered or collocated, auto-detected from nodality)
   !> User must provide MultiFabs at finest level with >= 2 ghost cells filled
   subroutine advance_vof(this,U,V,W,dt,time)
      use amrex_amr_module, only: amrex_multifab
      implicit none
      class(amrvof), intent(inout) :: this
      type(amrex_multifab), intent(in) :: U,V,W
      real(WP), intent(in) :: dt
      real(WP), intent(in) :: time
      type(amrex_multifab) :: band,Fx,Fy,Fz
      logical :: is_staggered
      integer :: lvl
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,vol

      ! Shared variables for internal functions
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW ! Velocity used for project
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pPLICold ! PLIC old used in tet2flux_plic

      ! Level at which we're working
      lvl=this%amr%clvl()

      ! Mesh info
      dx=this%amr%dx(lvl); dxi=1.0_WP/dx
      dy=this%amr%dy(lvl); dyi=1.0_WP/dy
      dz=this%amr%dz(lvl); dzi=1.0_WP/dz
      vol=dx*dy*dz

      ! Check velocity centering and ghost cell requirements
      check_velocity: block
         use messager, only: die
         logical, dimension(3) :: nodal_U,nodal_V,nodal_W
         nodal_U=U%nodal_type()
         nodal_V=V%nodal_type()
         nodal_W=W%nodal_type()
         if (all(nodal_U .eqv. [.true. ,.false.,.false.]) .and. & 
         &   all(nodal_V .eqv. [.false.,.true. ,.false.]) .and. & 
         &   all(nodal_W .eqv. [.false.,.false.,.true. ])) then
            is_staggered=.true.
         else if (.not.any(nodal_U).and..not.any(nodal_V).and..not.any(nodal_W)) then
            is_staggered=.false.
         else
            call die('[advance_vof] velocity must be either staggered (face-centered) or collocated (cell-centered)')
         end if
         if (U%nghost().lt.2.or.V%nghost().lt.2.or.W%nghost().lt.2) then
            call die('[advance_vof] velocity requires >= 2 ghost cells')
         end if
      end block check_velocity
      
      ! Build transport band to localize computation
      build_band: block
         integer :: dir,n,i,j,k
         integer, dimension(3) :: ind
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pBand,pVFold
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         ! Build MultiFab with 1 ghost cell
         call this%amr%mfab_build(lvl=lvl,mfab=band,ncomp=1,nover=1)
         ! Pass 1: Mark interface cells (band=1)
         call band%setval(0.0_WP)
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            bx=mfi%tilebox()
            pVFold=>this%VFold%mf(lvl)%dataptr(mfi)
            pBand =>band%dataptr(mfi)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Mixed cell
               if (pVFold(i,j,k,1).ge.VFlo.and.pVFold(i,j,k,1).le.VFhi) then
                  pBand(i,j,k,1)=1.0_WP
               ! Implicit interface: pure cell adjacent to opposite phase
               else
                  do dir=1,3; do n=-1,+1,2
                     ind=[i,j,k]; ind(dir)=ind(dir)+n
                     if (pVFold(i,j,k,1).lt.VFlo.and.pVFold(ind(1),ind(2),ind(3),1).gt.VFhi.or.&
                     &   pVFold(i,j,k,1).gt.VFhi.and.pVFold(ind(1),ind(2),ind(3),1).lt.VFlo) then
                        pBand(i,j,k,1)=1.0_WP
                        cycle
                     end if
                  end do; end do
               end if
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         ! Nullify VFold and pBand
         nullify(pVFold,pBand)
         ! Synchronize within level
         call band%fill_boundary(this%amr%geom(lvl))
         ! Pass 2: Extend by 1 layer (band=2)
         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx=mfi%tilebox()
            pBand=>band%dataptr(mfi)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               if (pBand(i,j,k,1).eq.0.0_WP.and.any(pBand(i-1:i+1,j-1:j+1,k-1:k+1,1).eq.1.0_WP)) pBand(i,j,k,1)=2.0_WP
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         ! Synchronize within level
         call band%fill_boundary(this%amr%geom(lvl))
      end block build_band

      ! Build face-centered flux MultiFabs (8 components: Lvol, Gvol, Lbar(3), Gbar(3))
      call this%amr%mfab_build(lvl=lvl,mfab=Fx,ncomp=8,nover=0,atface=[.true. ,.false.,.false.]); call Fx%setval(0.0_WP)
      call this%amr%mfab_build(lvl=lvl,mfab=Fy,ncomp=8,nover=0,atface=[.false.,.true. ,.false.]); call Fy%setval(0.0_WP)
      call this%amr%mfab_build(lvl=lvl,mfab=Fz,ncomp=8,nover=0,atface=[.false.,.false.,.true. ]); call Fz%setval(0.0_WP)
      
      ! Phase 1: Compute all fluxes
      compute_fluxes: block
         use amrvof_geometry, only: tet_sign,tet_map,correct_flux_poly
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: fbx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pBand,pFx,pFy,pFz
         integer :: i,j,k,n,nn
         real(WP), dimension(3,9) :: face
         real(WP), dimension(3,4) :: tet
         integer , dimension(3,4) :: ijk
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get data pointers: PLICold, band, velocity, fluxes
            pPLICold=>this%PLICold%mf(lvl)%dataptr(mfi)
            pBand=>band%dataptr(mfi)
            pU =>U%dataptr(mfi)
            pV =>V%dataptr(mfi)
            pW =>W%dataptr(mfi)
            pFx=>Fx%dataptr(mfi)
            pFy=>Fy%dataptr(mfi)
            pFz=>Fz%dataptr(mfi)
            ! X-fluxes: loop over nodaltilebox(1)
            fbx=mfi%nodaltilebox(1)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               ! Skip if outside band
               if (maxval(pBand(i-1:i,j,k,1)).eq.0.0_WP) cycle
               ! Face vertices: 1-4 current position, 5-8 projected back, 9 at backface barycenter
               face(:,1)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,5)=project(face(:,1),-dt)
               face(:,2)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,6)=project(face(:,2),-dt)
               face(:,3)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,7)=project(face(:,3),-dt)
               face(:,4)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,8)=project(face(:,4),-dt)
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               ! Set 9th vertex to inforce target volume
               call correct_flux_poly(poly=face,target_volume=dt*dy*dz*merge(pU(i,j,k,1),0.5_WP*(pU(i-1,j,k,1)+pU(i,j,k,1)),is_staggered))
               ! Compute sign of each tet and accumulate flux
               pFx(i,j,k,1:8)=0.0_WP
               do n=1,8
                  ! Get the vertices and indices
                  do nn=1,4
                     tet(:,nn)=face(:,tet_map(nn,n))
                     ijk(:,nn)=floor([(tet(1,nn)-this%amr%xlo)*dxi,(tet(2,nn)-this%amr%ylo)*dyi,(tet(3,nn)-this%amr%zlo)*dzi])
                  end do
                  ! Cut tet and accumulate
                  pFx(i,j,k,1:8)=pFx(i,j,k,1:8)+tet_sign(tet)*tet2flux(tet,ijk)
               end do
            end do; end do; end do
            ! Y-fluxes: loop over nodaltilebox(2)
            fbx=mfi%nodaltilebox(2)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               ! Skip if outside band
               if (maxval(pBand(i,j-1:j,k,1)).eq.0.0_WP) cycle
               ! Face vertices: 1-4 current position, 5-8 projected back, 9 at backface barycenter
               face(:,1)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,5)=project(face(:,1),-dt)
               face(:,2)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,6)=project(face(:,2),-dt)
               face(:,3)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,7)=project(face(:,3),-dt)
               face(:,4)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,8)=project(face(:,4),-dt)
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               ! Set 9th vertex to inforce target volume
               call correct_flux_poly(poly=face,target_volume=dt*dz*dx*merge(pV(i,j,k,1),0.5_WP*(pV(i,j-1,k,1)+pV(i,j,k,1)),is_staggered))
               ! Compute sign of each tet and accumulate flux
               pFy(i,j,k,1:8)=0.0_WP
               do n=1,8
                  ! Get the vertices and indices
                  do nn=1,4
                     tet(:,nn)=face(:,tet_map(nn,n))
                     ijk(:,nn)=floor([(tet(1,nn)-this%amr%xlo)*dxi,(tet(2,nn)-this%amr%ylo)*dyi,(tet(3,nn)-this%amr%zlo)*dzi])
                  end do
                  ! Cut tet and accumulate
                  pFy(i,j,k,1:8)=pFy(i,j,k,1:8)+tet_sign(tet)*tet2flux(tet,ijk)
               end do
            end do; end do; end do
            ! Z-fluxes: loop over nodaltilebox(3)
            fbx=mfi%nodaltilebox(3)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               ! Skip if outside band
               if (maxval(pBand(i,j,k-1:k,1)).eq.0.0_WP) cycle
               ! Face vertices: 1-4 current position, 5-8 projected back, 9 at backface barycenter
               face(:,1)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,5)=project(face(:,1),-dt)
               face(:,2)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,6)=project(face(:,2),-dt)
               face(:,3)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,7)=project(face(:,3),-dt)
               face(:,4)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,8)=project(face(:,4),-dt)
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               ! Set 9th vertex to inforce target volume
               call correct_flux_poly(poly=face,target_volume=dt*dx*dy*merge(pW(i,j,k,1),0.5_WP*(pW(i,j,k-1,1)+pW(i,j,k,1)),is_staggered))
               ! Compute sign of each tet and accumulate flux
               pFz(i,j,k,1:8)=0.0_WP
               do n=1,8
                  ! Get the vertices and indices
                  do nn=1,4
                     tet(:,nn)=face(:,tet_map(nn,n))
                     ijk(:,nn)=floor([(tet(1,nn)-this%amr%xlo)*dxi,(tet(2,nn)-this%amr%ylo)*dyi,(tet(3,nn)-this%amr%zlo)*dzi])
                  end do
                  ! Cut tet and accumulate
                  pFz(i,j,k,1:8)=pFz(i,j,k,1:8)+tet_sign(tet)*tet2flux(tet,ijk)
               end do
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         ! Nullify pointers
         nullify(pU,pV,pW,pPLICold)
      end block compute_fluxes

      ! Phase 2: Update VF from fluxes
      update_vf: block
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         integer :: i,j,k
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pBand,pFx,pFy,pFz
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCliq,pCgas,pVFold,pCliqold,pCgasold
         real(WP) :: Lvol_old,Lvol_new,Lvol_flux
         real(WP) :: Gvol_old,Gvol_new,Gvol_flux
         real(WP), dimension(3) :: Lbar_old,Lbar_new,Lbar_flux
         real(WP), dimension(3) :: Gbar_old,Gbar_new,Gbar_flux
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get data pointers
            pBand=>band%dataptr(mfi)
            pVF=>this%VF%mf(lvl)%dataptr(mfi)
            pCliq=>this%Cliq%mf(lvl)%dataptr(mfi)
            pCgas=>this%Cgas%mf(lvl)%dataptr(mfi)
            pVFold=>this%VFold%mf(lvl)%dataptr(mfi)
            pCliqold=>this%Cliqold%mf(lvl)%dataptr(mfi)
            pCgasold=>this%Cgasold%mf(lvl)%dataptr(mfi)
            pU=>U%dataptr(mfi)
            pV=>V%dataptr(mfi)
            pW=>W%dataptr(mfi)
            pFx=>Fx%dataptr(mfi)
            pFy=>Fy%dataptr(mfi)
            pFz=>Fz%dataptr(mfi)
            ! Loop over tilebox
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Skip if cell not in band
               if (pBand(i,j,k,1).eq.0.0_WP) cycle
               ! Old phasic moments
               Lvol_old=(       pVFold(i,j,k,1))*vol
               Gvol_old=(1.0_WP-pVFold(i,j,k,1))*vol
               Lbar_old=pCliqold(i,j,k,1:3)
               Gbar_old=pCgasold(i,j,k,1:3)
               ! Net flux (outflow positive)
               Lvol_flux=pFx(i+1,j,k,1)  -pFx(i,j,k,1)  +pFy(i,j+1,k,1)  -pFy(i,j,k,1)  +pFz(i,j,k+1,1)  -pFz(i,j,k,1)
               Gvol_flux=pFx(i+1,j,k,2)  -pFx(i,j,k,2)  +pFy(i,j+1,k,2)  -pFy(i,j,k,2)  +pFz(i,j,k+1,2)  -pFz(i,j,k,2)
               Lbar_flux=pFx(i+1,j,k,3:5)-pFx(i,j,k,3:5)+pFy(i,j+1,k,3:5)-pFy(i,j,k,3:5)+pFz(i,j,k+1,3:5)-pFz(i,j,k,3:5)
               Gbar_flux=pFx(i+1,j,k,6:8)-pFx(i,j,k,6:8)+pFy(i,j+1,k,6:8)-pFy(i,j,k,6:8)+pFz(i,j,k+1,6:8)-pFz(i,j,k,6:8)
               ! New phasic volumes
               Lvol_new=Lvol_old-Lvol_flux
               Gvol_new=Gvol_old-Gvol_flux
               ! New volume fraction and default barycenters
               pVF(i,j,k,1)=Lvol_new/(Lvol_new+Gvol_new)
               pCliq(i,j,k,1:3) = [this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
               pCgas(i,j,k,1:3) = [this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
               ! Clip and update barycenters
               if (pVF(i,j,k,1).lt.VFlo) then
                  pVF(i,j,k,1)=0.0_WP
               else if (pVF(i,j,k,1).gt.VFhi) then
                  pVF(i,j,k,1)=1.0_WP
               else
                  ! Update barycenters from moment conservation and project forward
                  if (Lvol_new/(Lvol_new+Gvol_new).gt.vol_eps) then; Lbar_new=(Lbar_old*Lvol_old-Lbar_flux)/Lvol_new; pCliq(i,j,k,1:3)=project(Lbar_new,dt); end if
                  if (Gvol_new/(Lvol_new+Gvol_new).gt.vol_eps) then; Gbar_new=(Gbar_old*Gvol_old-Gbar_flux)/Gvol_new; pCgas(i,j,k,1:3)=project(Gbar_new,dt); end if
               end if
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         ! Nullify pointers
         nullify(pU,pV,pW)
      end block update_vf

      ! Cleanup flux multifabs
      call this%amr%mfab_destroy(Fx)
      call this%amr%mfab_destroy(Fy)
      call this%amr%mfab_destroy(Fz)

      ! Clean up band multifab
      call this%amr%mfab_destroy(band)
      
      ! Sync and apply BC
      call this%fill_moments_lvl(lvl,time)

   contains
      
      !> Recursive function that cuts a tet by grid planes to compute fluxes
      recursive function tet2flux(mytet, myind) result(myflux)
         use amrvof_geometry, only: cut_side, cut_v1, cut_v2, cut_vtet, cut_ntets, cut_nvert
         real(WP), dimension(3,4), intent(in) :: mytet
         integer, dimension(3,4), intent(in) :: myind
         real(WP), dimension(8) :: myflux
         integer :: dir, cut_ind, icase, n1, n2, v1, v2
         real(WP), dimension(4) :: dd
         real(WP), dimension(3,8) :: vert
         integer, dimension(3,8,2) :: vert_ind
         real(WP) :: mu, my_vol
         real(WP), dimension(3,4) :: newtet
         integer, dimension(3,4) :: newind
         real(WP), dimension(3) :: a, b, c
         real(WP) :: xcut, ycut, zcut
         
         myflux = 0.0_WP
         
         ! Determine if tet spans multiple cells and needs cutting
         if (maxval(myind(1,:)) - minval(myind(1,:)) .gt. 0) then
            ! Cut by x planes
            dir = 1
            cut_ind = maxval(myind(1,:))
            xcut = this%amr%xlo + real(cut_ind,WP) * dx
            dd(:) = mytet(1,:) - xcut
         else if (maxval(myind(2,:)) - minval(myind(2,:)) .gt. 0) then
            ! Cut by y planes
            dir = 2
            cut_ind = maxval(myind(2,:))
            ycut = this%amr%ylo + real(cut_ind,WP) * dy
            dd(:) = mytet(2,:) - ycut
         else if (maxval(myind(3,:)) - minval(myind(3,:)) .gt. 0) then
            ! Cut by z planes
            dir = 3
            cut_ind = maxval(myind(3,:))
            zcut = this%amr%zlo + real(cut_ind,WP) * dz
            dd(:) = mytet(3,:) - zcut
         else
            ! All vertices in same cell - cut by PLIC and return
            myflux = tet2flux_plic(mytet, myind(1,1), myind(2,1), myind(3,1))
            return
         end if
         
         ! Find cut case (1-indexed: 1-16)
         icase = 1 + int(0.5_WP + sign(0.5_WP, dd(1))) &
               + 2 * int(0.5_WP + sign(0.5_WP, dd(2))) &
               + 4 * int(0.5_WP + sign(0.5_WP, dd(3))) &
               + 8 * int(0.5_WP + sign(0.5_WP, dd(4)))
         
         ! Copy vertices and indices
         do n1 = 1, 4
            vert(:, n1) = mytet(:, n1)
            vert_ind(:, n1, 1) = myind(:, n1)
            vert_ind(:, n1, 2) = myind(:, n1)
            ! Enforce boundedness at cut plane
            vert_ind(dir, n1, 1) = min(vert_ind(dir, n1, 1), cut_ind - 1)
            vert_ind(dir, n1, 2) = max(vert_ind(dir, n1, 1), cut_ind)
         end do
         
         ! Create interpolated vertices on cut plane
         do n1 = 1, cut_nvert(icase)
            v1 = cut_v1(n1, icase); v2 = cut_v2(n1, icase)
            mu = min(1.0_WP,max(0.0_WP,-dd(v1)/(sign(abs(dd(v2)-dd(v1))+epsilon(1.0_WP),dd(v2)-dd(v1)))))
            vert(:, 4 + n1) = (1.0_WP - mu) * vert(:, v1) + mu * vert(:, v2)
            ! Compute index for interpolated vertex
            vert_ind(1, 4+n1, 1) = floor((vert(1, 4+n1) - this%amr%xlo) * dxi)
            vert_ind(2, 4+n1, 1) = floor((vert(2, 4+n1) - this%amr%ylo) * dyi)
            vert_ind(3, 4+n1, 1) = floor((vert(3, 4+n1) - this%amr%zlo) * dzi)
            ! Enforce boundedness
            vert_ind(:, 4+n1, 1) = max(vert_ind(:, 4+n1, 1), min(vert_ind(:, v1, 1), vert_ind(:, v2, 1)))
            vert_ind(:, 4+n1, 1) = min(vert_ind(:, 4+n1, 1), max(vert_ind(:, v1, 1), vert_ind(:, v2, 1)))
            ! Set +/- indices in cut direction
            vert_ind(:, 4+n1, 2) = vert_ind(:, 4+n1, 1)
            vert_ind(dir, 4+n1, 1) = cut_ind - 1
            vert_ind(dir, 4+n1, 2) = cut_ind
         end do
         
         ! Create and process sub-tets
         do n1 = 1, cut_ntets(icase)
            do n2 = 1, 4
               newtet(:, n2) = vert(:, cut_vtet(n2, n1, icase))
               newind(:, n2) = vert_ind(:, cut_vtet(n2, n1, icase), cut_side(n1, icase))
            end do
            ! Check for zero-volume tet
            a = newtet(:,1) - newtet(:,4)
            b = newtet(:,2) - newtet(:,4)
            c = newtet(:,3) - newtet(:,4)
            my_vol = abs(a(1)*(b(2)*c(3)-c(2)*b(3)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
            if (my_vol .lt. 1.0e-15_WP * vol) cycle
            ! Recursively process sub-tet
            myflux = myflux + tet2flux(newtet, newind)
         end do
         
      end function tet2flux
      
      !> Cut tet by PLIC and compute flux (base case of recursion) - uses pPLICold
      function tet2flux_plic(mytet, i0, j0, k0) result(myflux)
         use amrvof_geometry, only: cut_v1, cut_v2, cut_vtet, cut_ntets, cut_nvert, cut_nntet, tet_vol
         use messager, only: die
         real(WP), dimension(3,4), intent(in) :: mytet
         integer, intent(in) :: i0, j0, k0
         real(WP), dimension(8) :: myflux
         integer :: icase, n1, v1, v2
         real(WP), dimension(4) :: dd
         real(WP), dimension(3,8) :: vert
         real(WP), dimension(3) :: a, b, c, bary, normal
         real(WP) :: mu, my_vol, dist
         
         myflux = 0.0_WP

         ! Check indices are within PLICold bounds
         if (i0 .lt. lbound(pPLICold,1) .or. i0 .gt. ubound(pPLICold,1) .or. &
             j0 .lt. lbound(pPLICold,2) .or. j0 .gt. ubound(pPLICold,2) .or. &
             k0 .lt. lbound(pPLICold,3) .or. k0 .gt. ubound(pPLICold,3)) then
            call die('[tet2flux_plic] Index out of bounds - check CFL or ghost cells')
         end if

         ! Pure cell shortcut - skip PLIC cutting
         if (pPLICold(i0,j0,k0,4).gt.+1.0e9_WP) then
            ! Pure liquid - all volume goes to liquid phase
            my_vol=abs(tet_vol(mytet))
            bary=0.25_WP*(mytet(:,1)+mytet(:,2)+mytet(:,3)+mytet(:,4))
            myflux( 1 )=my_vol
            myflux(3:5)=my_vol*bary
            return
         else if (pPLICold(i0,j0,k0,4).lt.-1.0e9_WP) then
            ! Pure gas - all volume goes to gas phase
            my_vol=abs(tet_vol(mytet))
            bary=0.25_WP*(mytet(:,1)+mytet(:,2)+mytet(:,3)+mytet(:,4))
            myflux( 2 )=my_vol
            myflux(6:8)=my_vol*bary
            return
         end if
         
         ! Get PLIC from this cell
         normal = pPLICold(i0, j0, k0, 1:3)
         dist = pPLICold(i0, j0, k0, 4)
         
         ! Compute signed distance to plane for each vertex
         dd(1) = normal(1)*mytet(1,1) + normal(2)*mytet(2,1) + normal(3)*mytet(3,1) - dist
         dd(2) = normal(1)*mytet(1,2) + normal(2)*mytet(2,2) + normal(3)*mytet(3,2) - dist
         dd(3) = normal(1)*mytet(1,3) + normal(2)*mytet(2,3) + normal(3)*mytet(3,3) - dist
         dd(4) = normal(1)*mytet(1,4) + normal(2)*mytet(2,4) + normal(3)*mytet(3,4) - dist
         
         ! Find cut case
         icase = 1 + int(0.5_WP + sign(0.5_WP, dd(1))) &
               + 2 * int(0.5_WP + sign(0.5_WP, dd(2))) &
               + 4 * int(0.5_WP + sign(0.5_WP, dd(3))) &
               + 8 * int(0.5_WP + sign(0.5_WP, dd(4)))
         
         ! Copy vertices
         vert(:, 1:4) = mytet(:, 1:4)
         
         ! Create interpolated vertices on cut plane
         do n1 = 1, cut_nvert(icase)
            v1 = cut_v1(n1, icase); v2 = cut_v2(n1, icase)
            mu = min(1.0_WP,max(0.0_WP,-dd(v1)/(sign(abs(dd(v2)-dd(v1))+epsilon(1.0_WP),dd(v2)-dd(v1)))))
            vert(:, 4 + n1) = (1.0_WP - mu) * vert(:, v1) + mu * vert(:, v2)
         end do
         
         ! Gas tets: from 1 to cut_nntet-1
         do n1 = 1, cut_nntet(icase) - 1
            a = vert(:, cut_vtet(1, n1, icase)) - vert(:, cut_vtet(4, n1, icase))
            b = vert(:, cut_vtet(2, n1, icase)) - vert(:, cut_vtet(4, n1, icase))
            c = vert(:, cut_vtet(3, n1, icase)) - vert(:, cut_vtet(4, n1, icase))
            my_vol = abs(a(1)*(b(2)*c(3)-c(2)*b(3)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
            bary = 0.25_WP * (vert(:, cut_vtet(1, n1, icase)) + vert(:, cut_vtet(2, n1, icase)) &
            &               + vert(:, cut_vtet(3, n1, icase)) + vert(:, cut_vtet(4, n1, icase)))
            myflux( 2 ) = myflux( 2 ) + my_vol
            myflux(6:8) = myflux(6:8) + my_vol * bary
         end do
         
         ! Liquid tets: from cut_ntets down to cut_nntet
         do n1 = cut_ntets(icase), cut_nntet(icase), -1
            a = vert(:, cut_vtet(1, n1, icase)) - vert(:, cut_vtet(4, n1, icase))
            b = vert(:, cut_vtet(2, n1, icase)) - vert(:, cut_vtet(4, n1, icase))
            c = vert(:, cut_vtet(3, n1, icase)) - vert(:, cut_vtet(4, n1, icase))
            my_vol = abs(a(1)*(b(2)*c(3)-c(2)*b(3)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
            bary = 0.25_WP * (vert(:, cut_vtet(1, n1, icase)) + vert(:, cut_vtet(2, n1, icase)) &
            &               + vert(:, cut_vtet(3, n1, icase)) + vert(:, cut_vtet(4, n1, icase)))
            myflux( 1 ) = myflux( 1 ) + my_vol
            myflux(3:5) = myflux(3:5) + my_vol * bary
         end do
         
      end function tet2flux_plic
      
      !> RK2 vertex projection back in time
      function project(p1,mydt) result(p2)
         implicit none
         real(WP), dimension(3), intent(in) :: p1
         real(WP), dimension(3)             :: p2
         real(WP),               intent(in) :: mydt
         p2=p1+mydt*interp_velocity(        p1    )
         p2=p1+mydt*interp_velocity(0.5_WP*(p1+p2))
      end function project
      
      !> Trilinear interpolation of velocity (handles staggered or collocated) - uses pU,pV,pW
      function interp_velocity(pos) result(vel)
         implicit none
         real(WP), dimension(3), intent(in) :: pos
         real(WP), dimension(3) :: vel
         integer  :: ipc, jpc, kpc   ! Cell-centered indices
         integer  :: ipu, jpv, kpw   ! Face-centered indices
         real(WP) :: wxc1, wyc1, wzc1, wxc2, wyc2, wzc2  ! Cell-centered weights
         real(WP) :: wxu1, wyv1, wzw1, wxu2, wyv2, wzw2  ! Face-centered weights
         if (is_staggered) then
            ! Compute raw indices
            ipc = floor((pos(1) - this%amr%xlo) * dxi - 0.5_WP)
            jpc = floor((pos(2) - this%amr%ylo) * dyi - 0.5_WP)
            kpc = floor((pos(3) - this%amr%zlo) * dzi - 0.5_WP)
            ipu = floor((pos(1) - this%amr%xlo) * dxi)
            jpv = floor((pos(2) - this%amr%ylo) * dyi)
            kpw = floor((pos(3) - this%amr%zlo) * dzi)
            ! Clamp to array bounds
            if (ipu<lbound(pU,1).or.ipu>ubound(pU,1)-1.or.jpc<lbound(pU,2).or.jpc>ubound(pU,2)-1.or. &
                kpc<lbound(pU,3).or.kpc>ubound(pU,3)-1.or.ipc<lbound(pV,1).or.ipc>ubound(pV,1)-1.or. &
                jpv<lbound(pV,2).or.jpv>ubound(pV,2)-1.or.kpw<lbound(pW,3).or.kpw>ubound(pW,3)-1) then
               print*,'Interpolation out of bounds',ipu,jpc,kpc,ipc,jpv,kpw
            end if
            ipu = max(lbound(pU,1), min(ubound(pU,1)-1, ipu))
            jpc = max(lbound(pU,2), min(ubound(pU,2)-1, jpc))
            kpc = max(lbound(pU,3), min(ubound(pU,3)-1, kpc))
            ipc = max(lbound(pV,1), min(ubound(pV,1)-1, ipc))
            jpv = max(lbound(pV,2), min(ubound(pV,2)-1, jpv))
            kpw = max(lbound(pW,3), min(ubound(pW,3)-1, kpw))
            ! Cell-centered weights
            wxc1 = (pos(1) - (this%amr%xlo + (real(ipc,WP)+0.5_WP)*dx)) * dxi
            wyc1 = (pos(2) - (this%amr%ylo + (real(jpc,WP)+0.5_WP)*dy)) * dyi
            wzc1 = (pos(3) - (this%amr%zlo + (real(kpc,WP)+0.5_WP)*dz)) * dzi
            wxc1 = max(0.0_WP, min(1.0_WP, wxc1)); wxc2 = 1.0_WP - wxc1
            wyc1 = max(0.0_WP, min(1.0_WP, wyc1)); wyc2 = 1.0_WP - wyc1
            wzc1 = max(0.0_WP, min(1.0_WP, wzc1)); wzc2 = 1.0_WP - wzc1
            ! Face-centered weights
            wxu1 = (pos(1) - (this%amr%xlo + real(ipu,WP)*dx)) * dxi
            wyv1 = (pos(2) - (this%amr%ylo + real(jpv,WP)*dy)) * dyi
            wzw1 = (pos(3) - (this%amr%zlo + real(kpw,WP)*dz)) * dzi
            wxu1 = max(0.0_WP, min(1.0_WP, wxu1)); wxu2 = 1.0_WP - wxu1
            wyv1 = max(0.0_WP, min(1.0_WP, wyv1)); wyv2 = 1.0_WP - wyv1
            wzw1 = max(0.0_WP, min(1.0_WP, wzw1)); wzw2 = 1.0_WP - wzw1
            ! U at x-faces: face-centered in x, cell-centered in y,z
            vel(1) = wzc1*(wyc1*(wxu1*pU(ipu+1,jpc+1,kpc+1,1)+wxu2*pU(ipu,jpc+1,kpc+1,1)) + &
            &              wyc2*(wxu1*pU(ipu+1,jpc  ,kpc+1,1)+wxu2*pU(ipu,jpc  ,kpc+1,1))) + &
            &        wzc2*(wyc1*(wxu1*pU(ipu+1,jpc+1,kpc  ,1)+wxu2*pU(ipu,jpc+1,kpc  ,1)) + &
            &              wyc2*(wxu1*pU(ipu+1,jpc  ,kpc  ,1)+wxu2*pU(ipu,jpc  ,kpc  ,1)))
            ! V at y-faces: cell-centered in x, face-centered in y, cell-centered in z
            vel(2) = wzc1*(wyv1*(wxc1*pV(ipc+1,jpv+1,kpc+1,1)+wxc2*pV(ipc,jpv+1,kpc+1,1)) + &
            &              wyv2*(wxc1*pV(ipc+1,jpv  ,kpc+1,1)+wxc2*pV(ipc,jpv  ,kpc+1,1))) + &
            &        wzc2*(wyv1*(wxc1*pV(ipc+1,jpv+1,kpc  ,1)+wxc2*pV(ipc,jpv+1,kpc  ,1)) + &
            &              wyv2*(wxc1*pV(ipc+1,jpv  ,kpc  ,1)+wxc2*pV(ipc,jpv  ,kpc  ,1)))
            ! W at z-faces: cell-centered in x,y, face-centered in z
            vel(3) = wzw1*(wyc1*(wxc1*pW(ipc+1,jpc+1,kpw+1,1)+wxc2*pW(ipc,jpc+1,kpw+1,1)) + &
            &              wyc2*(wxc1*pW(ipc+1,jpc  ,kpw+1,1)+wxc2*pW(ipc,jpc  ,kpw+1,1))) + &
            &        wzw2*(wyc1*(wxc1*pW(ipc+1,jpc+1,kpw  ,1)+wxc2*pW(ipc,jpc+1,kpw  ,1)) + &
            &              wyc2*(wxc1*pW(ipc+1,jpc  ,kpw  ,1)+wxc2*pW(ipc,jpc  ,kpw  ,1)))
         else
            ! All cell-centered
            ipc = floor((pos(1) - this%amr%xlo) * dxi - 0.5_WP)
            jpc = floor((pos(2) - this%amr%ylo) * dyi - 0.5_WP)
            kpc = floor((pos(3) - this%amr%zlo) * dzi - 0.5_WP)
            ! Clamp to array bounds
            ipc = max(lbound(pU,1), min(ubound(pU,1)-1, ipc))
            jpc = max(lbound(pU,2), min(ubound(pU,2)-1, jpc))
            kpc = max(lbound(pU,3), min(ubound(pU,3)-1, kpc))
            ! Cell-centered weights
            wxc1 = (pos(1) - (this%amr%xlo + (real(ipc,WP)+0.5_WP)*dx)) * dxi
            wyc1 = (pos(2) - (this%amr%ylo + (real(jpc,WP)+0.5_WP)*dy)) * dyi
            wzc1 = (pos(3) - (this%amr%zlo + (real(kpc,WP)+0.5_WP)*dz)) * dzi
            wxc1 = max(0.0_WP, min(1.0_WP, wxc1)); wxc2 = 1.0_WP - wxc1
            wyc1 = max(0.0_WP, min(1.0_WP, wyc1)); wyc2 = 1.0_WP - wyc1
            wzc1 = max(0.0_WP, min(1.0_WP, wzc1)); wzc2 = 1.0_WP - wzc1
            vel(1) = wzc1*(wyc1*(wxc1*pU(ipc+1,jpc+1,kpc+1,1)+wxc2*pU(ipc,jpc+1,kpc+1,1)) + &
            &              wyc2*(wxc1*pU(ipc+1,jpc  ,kpc+1,1)+wxc2*pU(ipc,jpc  ,kpc+1,1))) + &
            &        wzc2*(wyc1*(wxc1*pU(ipc+1,jpc+1,kpc  ,1)+wxc2*pU(ipc,jpc+1,kpc  ,1)) + &
            &              wyc2*(wxc1*pU(ipc+1,jpc  ,kpc  ,1)+wxc2*pU(ipc,jpc  ,kpc  ,1)))
            vel(2) = wzc1*(wyc1*(wxc1*pV(ipc+1,jpc+1,kpc+1,1)+wxc2*pV(ipc,jpc+1,kpc+1,1)) + &
            &              wyc2*(wxc1*pV(ipc+1,jpc  ,kpc+1,1)+wxc2*pV(ipc,jpc  ,kpc+1,1))) + &
            &        wzc2*(wyc1*(wxc1*pV(ipc+1,jpc+1,kpc  ,1)+wxc2*pV(ipc,jpc+1,kpc  ,1)) + &
            &              wyc2*(wxc1*pV(ipc+1,jpc  ,kpc  ,1)+wxc2*pV(ipc,jpc  ,kpc  ,1)))
            vel(3) = wzc1*(wyc1*(wxc1*pW(ipc+1,jpc+1,kpc+1,1)+wxc2*pW(ipc,jpc+1,kpc+1,1)) + &
            &              wyc2*(wxc1*pW(ipc+1,jpc  ,kpc+1,1)+wxc2*pW(ipc,jpc  ,kpc+1,1))) + &
            &        wzc2*(wyc1*(wxc1*pW(ipc+1,jpc+1,kpc  ,1)+wxc2*pW(ipc,jpc+1,kpc  ,1)) + &
            &              wyc2*(wxc1*pW(ipc+1,jpc  ,kpc  ,1)+wxc2*pW(ipc,jpc  ,kpc  ,1)))
         end if
      end function interp_velocity

   end subroutine advance_vof


   ! ============================================================================
   ! DEFERRED METHODS
   ! ============================================================================

   !> Get solver information
   subroutine get_info(this)
      use parallel, only: MPI_REAL_WP
      use mpi_f08
      class(amrvof), intent(inout) :: this
      integer :: lvl, ierr

      ! Initialize
      this%VFmin = huge(1.0_WP)
      this%VFmax = -huge(1.0_WP)
      this%VFint = 0.0_WP

      ! Loop over levels
       do lvl = 0, this%amr%clvl()
         this%VFmin = min(this%VFmin, this%VF%get_min(lvl=lvl))
         this%VFmax = max(this%VFmax, this%VF%get_max(lvl=lvl))
      end do

      ! Compute volume integral at level 0
      this%VFint = this%VF%get_sum(lvl=0) * this%amr%dx(0) * this%amr%dy(0) * this%amr%dz(0)

      ! Reduce across MPI ranks
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%VFint, 1, MPI_REAL_WP, MPI_SUM, this%amr%comm, ierr)
   end subroutine get_info

   !> Register checkpoint
   subroutine register_checkpoint(this, io)
      use amrio_class, only: amrio
      class(amrvof), intent(inout) :: this
      class(amrio), intent(inout) :: io
      call io%add_data(this%VF, 'VF')
      call io%add_data(this%Cliq, 'Cliq')
      call io%add_data(this%Cgas, 'Cgas')
      call io%add_data(this%PLIC, 'PLIC')
   end subroutine register_checkpoint

   !> Restore checkpoint
   subroutine restore_checkpoint(this, io, dirname)
      use amrio_class, only: amrio
      class(amrvof), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      call io%read_data(dirname, this%VF, 'VF')
      call io%read_data(dirname, this%Cliq, 'Cliq')
      call io%read_data(dirname, this%Cgas, 'Cgas')
      call io%read_data(dirname, this%PLIC, 'PLIC')
   end subroutine restore_checkpoint

   !> Compute advective CFL at finest level
   !> Takes external staggered velocity MultiFabs and returns max CFL
   subroutine get_cfl(this, U, V, W, dt, cfl)
      class(amrvof), intent(inout) :: this
      type(amrex_multifab), intent(in) :: U, V, W  !< Staggered velocity at clvl
      real(WP), intent(in) :: dt
      real(WP), intent(out) :: cfl
      real(WP) :: Umax, Vmax, Wmax, CFLx, CFLy, CFLz
      integer :: lvl
      ! Get finest level metrics
      lvl = this%amr%clvl()
      ! Get max velocity norms
      Umax = U%norm0()
      Vmax = V%norm0()
      Wmax = W%norm0()
      ! Compute directional CFLs
      CFLx = dt * Umax / this%amr%dx(lvl)
      CFLy = dt * Vmax / this%amr%dy(lvl)
      CFLz = dt * Wmax / this%amr%dz(lvl)
      ! Return max CFL
      cfl = max(CFLx, CFLy, CFLz)
   end subroutine get_cfl

   !> Print solver info to screen
   subroutine amrvof_print(this)
      use messager, only: log
      use string, only: str_long
      class(amrvof), intent(in) :: this
      character(len=str_long) :: message
      call log("VOF solver: "//trim(this%name))
      write(message,'("  VF range: [",ES12.5,", ",ES12.5,"]")') VFlo, VFhi
      call log(trim(message))
      call log("  Grid: "//trim(this%amr%name))
   end subroutine amrvof_print

end module amrvof_class
