!> AMR Incompressible solver class
!> Provides projection and basic operations for constant-density flow
module amrincomp_class
   use iso_c_binding,    only: c_ptr, c_null_ptr, c_loc, c_f_pointer
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata
   use amrsolver_class,  only: amrsolver
   use amrmg_class,      only: amrmg, amrmg_cstcoef
   use amrex_amr_module, only: amrex_multifab, amrex_mfiter, amrex_box, &
   &                           amrex_boxarray, amrex_distromap, amrex_geometry, &
   &                           amrex_interp_face_divfree
   implicit none
   private

   ! Expose type and dispatchers
   public :: amrincomp
   public :: amrincomp_on_init, amrincomp_on_coarse, amrincomp_on_remake
   public :: amrincomp_on_clear, amrincomp_tagging, amrincomp_postregrid

   !> AMR Incompressible solver type
   type, extends(amrsolver) :: amrincomp
      ! User-configurable callbacks
      procedure(incomp_init_iface), pointer, nopass :: user_init => null()
      procedure(incomp_tagging_iface), pointer, nopass :: user_tagging => null()
      procedure(incomp_dirichlet_iface), pointer, nopass :: user_dirichlet => null()

      ! Flow data (solver owns these)
      type(amrdata) :: U, V, W          !< Current face velocities
      type(amrdata) :: Uold, Vold, Wold !< Old face velocities
      type(amrdata) :: P                !< Pressure (cell-centered)
      type(amrdata) :: div              !< Divergence (cell-centered)

      ! Pressure solver
      type(amrmg) :: psolver

      ! Physical properties
      real(WP) :: rho  = 1.0_WP   !< Density (constant)
      real(WP) :: visc = 0.0_WP   !< Dynamic viscosity

      ! Monitoring quantities
      real(WP) :: Umax = 0.0_WP                 !< Max U velocity
      real(WP) :: Vmax = 0.0_WP                 !< Max V velocity
      real(WP) :: Wmax = 0.0_WP                 !< Max W velocity
      real(WP) :: Pmax = 0.0_WP                 !< Max pressure
      real(WP) :: divmax = 0.0_WP               !< Maximum divergence
      real(WP) :: rhoUint = 0.0_WP              !< Integral of rho*U
      real(WP) :: rhoVint = 0.0_WP              !< Integral of rho*V
      real(WP) :: rhoWint = 0.0_WP              !< Integral of rho*W
      real(WP) :: rhoKint = 0.0_WP              !< Integral of rho*K
      real(WP) :: CFLc_x = 0.0_WP, CFLc_y = 0.0_WP, CFLc_z = 0.0_WP  !< Convective CFL
      real(WP) :: CFLv_x = 0.0_WP, CFLv_y = 0.0_WP, CFLv_z = 0.0_WP  !< Viscous CFL
      real(WP) :: CFL = 0.0_WP                  !< Maximum CFL

      ! Velocity interpolation method (for C/F interface fills)
      integer :: interp_vel = amrex_interp_face_divfree

   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: get_div      !< Compute divergence (assumes velocity ghosts filled)
      procedure :: get_pgrad    !< Compute pressure gradient (assumes pressure ghosts filled)
      ! Override internal type-bound callbacks from amrsolver
      procedure :: on_init
      procedure :: on_coarse
      procedure :: on_remake
      procedure :: on_clear
      procedure :: post_regrid
      procedure :: fill_velocity_lvl         !< Fill velocity ghosts at single level
      procedure :: fill_velocity             !< Fill velocity ghosts on all levels
      procedure :: fill_velocity_from_coarse !< Fill velocity from coarse
      procedure :: sync_velocity_lvl         !< Sync velocity ghosts at single level
      procedure :: sync_velocity             !< Sync velocity ghosts on all levels
      procedure :: fill_velocity_mfab        !< Fill dest MultiFabs for regridding
      procedure :: average_down_velocity     !< Average down MAC velocity for C/F consistency
      procedure :: average_down_velocity_to  !< Average down MAC velocity for single level
      procedure :: average_down_pressure     !< Average down pressure for C/F consistency
      procedure :: average_down_pressure_to  !< Average down pressure for single level
      ! Deferred from amrsolver base class
      procedure :: get_info
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
      ! Physics procedures
      procedure :: get_dmomdt                 !< Compute momentum advection RHS
      procedure :: get_cfl                    !< Compute CFL numbers
      procedure :: correct_outflow            !< Correct outflow for global mass conservation
      procedure :: print => amrincomp_print   !< Print solver info
   end type amrincomp

   !> Abstract interface for user-overridable on_init callback
   abstract interface
      subroutine incomp_init_iface(solver, lvl, time, ba, dm)
         import :: amrincomp, WP, amrex_boxarray, amrex_distromap
         class(amrincomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine incomp_init_iface
   end interface

   !> Abstract interface for user-overridable tagging callback
   abstract interface
      subroutine incomp_tagging_iface(solver, lvl, tags, time)
         import :: amrincomp, c_ptr, WP
         class(amrincomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         type(c_ptr), intent(in) :: tags
         real(WP), intent(in) :: time
      end subroutine incomp_tagging_iface
   end interface

   !> Abstract interface for user-provided Dirichlet velocity BC values
   !> User receives a boundary box and fills it efficiently (e.g., Poiseuille profile)
   abstract interface
      subroutine incomp_dirichlet_iface(solver, bx, p, comp, face, time, geom)
         import :: amrincomp, amrex_box, amrex_geometry, WP
         class(amrincomp), intent(in) :: solver
         type(amrex_box), intent(in) :: bx          !< Boundary region to fill
         real(WP), dimension(:,:,:,:), pointer, intent(inout) :: p
         character(len=1), intent(in) :: comp       !< 'U', 'V', or 'W'
         integer, intent(in) :: face                !< 1=xlo,2=xhi,3=ylo,4=yhi,5=zlo,6=zhi
         real(WP), intent(in) :: time
         type(amrex_geometry), intent(in) :: geom
      end subroutine incomp_dirichlet_iface
   end interface

contains

   ! ============================================================================
   ! DISPATCHERS (module-level) - recover concrete amrincomp type
   ! ============================================================================

   !> Dispatch on_init: calls type-bound method then user callback
   subroutine amrincomp_on_init(ctx, lvl, time, ba, dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_init(lvl, time, ba, dm)
      if (associated(this%user_init)) call this%user_init(this, lvl, time, ba, dm)
   end subroutine amrincomp_on_init

   !> Dispatch on_coarse: calls type-bound method
   subroutine amrincomp_on_coarse(ctx, lvl, time, ba, dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_coarse(lvl, time, ba, dm)
   end subroutine amrincomp_on_coarse

   !> Dispatch on_remake: calls type-bound method
   subroutine amrincomp_on_remake(ctx, lvl, time, ba, dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_remake(lvl, time, ba, dm)
   end subroutine amrincomp_on_remake

   !> Dispatch on_clear: calls type-bound method
   subroutine amrincomp_on_clear(ctx, lvl)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_clear(lvl)
   end subroutine amrincomp_on_clear

   !> Dispatch tagging: calls user callback if set
   subroutine amrincomp_tagging(ctx, lvl, tags, time)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags
      real(WP), intent(in) :: time
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx, this)
      if (associated(this%user_tagging)) call this%user_tagging(this, lvl, tags, time)
   end subroutine amrincomp_tagging

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrincomp_postregrid(ctx, lbase, time)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx, this)
      call this%post_regrid(lbase, time)
   end subroutine amrincomp_postregrid

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize the incompressible solver
   subroutine initialize(this, amr, name)
      implicit none
      class(amrincomp), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name

      ! Set name
      if (present(name)) then
         this%name = trim(name)
      else
         this%name = 'UNNAMED_INCOMP'
      end if

      ! Store amrgrid pointer
      this%amr => amr

      ! Initialize staggered velocity (do we need ghosts?)
      call this%U%initialize(amr, name='U', ncomp=1, ng=1, nodal=[.true., .false., .false.])
      call this%V%initialize(amr, name='V', ncomp=1, ng=1, nodal=[.false., .true., .false.])
      call this%W%initialize(amr, name='W', ncomp=1, ng=1, nodal=[.false., .false., .true.])
      call this%Uold%initialize(amr, name='Uold', ncomp=1, ng=1, nodal=[.true., .false., .false.])
      call this%Vold%initialize(amr, name='Vold', ncomp=1, ng=1, nodal=[.false., .true., .false.])
      call this%Wold%initialize(amr, name='Wold', ncomp=1, ng=1, nodal=[.false., .false., .true.])

      ! Initialize pressure (cell-centered) (do we need ghosts?)
      call this%P%initialize(amr, name='P', ncomp=1, ng=1)

      ! Initialize divergence (cell-centered, no ghosts needed)
      call this%div%initialize(amr, name='div', ncomp=1, ng=0)

      ! Set parent pointers for BC context access
      this%U%parent => this; this%V%parent => this; this%W%parent => this
      this%Uold%parent => this; this%Vold%parent => this; this%Wold%parent => this
      this%P%parent => this
      this%div%parent => this

      ! Set velocity fillbc callbacks to shared handler
      this%U%fillbc => velocity_fillbc
      this%V%fillbc => velocity_fillbc
      this%W%fillbc => velocity_fillbc

      ! Initialize pressure solver
      call this%psolver%initialize(amr, type=amrmg_cstcoef)

      ! Register all 6 callbacks with amrgrid using concrete dispatchers
      select type (this)
       type is (amrincomp)
         call this%amr%add_on_init   (amrincomp_on_init,    c_loc(this))
         call this%amr%add_on_coarse (amrincomp_on_coarse,  c_loc(this))
         call this%amr%add_on_remake (amrincomp_on_remake,  c_loc(this))
         call this%amr%add_on_clear  (amrincomp_on_clear,   c_loc(this))
         call this%amr%add_tagging   (amrincomp_tagging,    c_loc(this))
         call this%amr%add_postregrid(amrincomp_postregrid, c_loc(this))
      end select

      ! Print solver info
      call this%print()

   end subroutine initialize

   !> Finalize the incompressible solver
   subroutine finalize(this)
      implicit none
      class(amrincomp), intent(inout) :: this
      call this%U%finalize()
      call this%V%finalize()
      call this%W%finalize()
      call this%Uold%finalize()
      call this%Vold%finalize()
      call this%Wold%finalize()
      call this%P%finalize()
      call this%div%finalize()
      call this%psolver%finalize()
      nullify(this%amr)
      nullify(this%user_init)
      nullify(this%user_tagging)
      nullify(this%user_dirichlet)
   end subroutine finalize

   ! ============================================================================
   ! INTERNAL CALLBACK OVERRIDES
   ! ============================================================================

   !> Override on_init: reset levels and set to zero
   subroutine on_init(this, lvl, time, ba, dm)
      implicit none
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Reset level layouts
      call this%U%reset_level(lvl, ba, dm)
      call this%V%reset_level(lvl, ba, dm)
      call this%W%reset_level(lvl, ba, dm)
      call this%Uold%reset_level(lvl, ba, dm)
      call this%Vold%reset_level(lvl, ba, dm)
      call this%Wold%reset_level(lvl, ba, dm)
      call this%P%reset_level(lvl, ba, dm)
      call this%div%reset_level(lvl, ba, dm)
      ! Set to zero
      call this%U%setval(val=0.0_WP, lvl=lvl)
      call this%V%setval(val=0.0_WP, lvl=lvl)
      call this%W%setval(val=0.0_WP, lvl=lvl)
      call this%Uold%setval(val=0.0_WP, lvl=lvl)
      call this%Vold%setval(val=0.0_WP, lvl=lvl)
      call this%Wold%setval(val=0.0_WP, lvl=lvl)
      call this%P%setval(val=0.0_WP, lvl=lvl)
      call this%div%setval(val=0.0_WP, lvl=lvl)
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse using divergence-free interpolation
   subroutine on_coarse(this, lvl, time, ba, dm)
      implicit none
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Velocity: allocate then fill with divergence-free interpolation
      call this%U%reset_level(lvl, ba, dm)
      call this%V%reset_level(lvl, ba, dm)
      call this%W%reset_level(lvl, ba, dm)
      call this%fill_velocity_from_coarse(lvl, time)
      ! Old velocity just needs geometry
      call this%Uold%reset_level(lvl, ba, dm)
      call this%Vold%reset_level(lvl, ba, dm)
      call this%Wold%reset_level(lvl, ba, dm)
      ! Pressure on coarse
      call this%P%on_coarse(this%P, lvl, time, ba, dm)
      ! Divergence just needs geometry
      call this%div%reset_level(lvl, ba, dm)
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid using divergence-free interpolation
   subroutine on_remake(this, lvl, time, ba, dm)
      use amrex_amr_module, only: amrex_multifab_build, amrex_multifab_destroy
      implicit none
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_multifab) :: Utmp, Vtmp, Wtmp
      ! Build temp MultiFabs with new layout (0 ghost cells for FillPatch)
      call amrex_multifab_build(Utmp, ba, dm, 1, 0, this%U%nodal)
      call amrex_multifab_build(Vtmp, ba, dm, 1, 0, this%V%nodal)
      call amrex_multifab_build(Wtmp, ba, dm, 1, 0, this%W%nodal)
      ! Fill temps from old data via coupled FillPatch
      call this%fill_velocity_mfab(Utmp, Vtmp, Wtmp, lvl, time)
      ! Reset levels and copy from temps
      call this%U%reset_level(lvl, ba, dm)
      call this%V%reset_level(lvl, ba, dm)
      call this%W%reset_level(lvl, ba, dm)
      call this%U%mf(lvl)%copy(Utmp, 1, 1, 1, 0)
      call this%V%mf(lvl)%copy(Vtmp, 1, 1, 1, 0)
      call this%W%mf(lvl)%copy(Wtmp, 1, 1, 1, 0)
      ! Destroy temps
      call amrex_multifab_destroy(Utmp)
      call amrex_multifab_destroy(Vtmp)
      call amrex_multifab_destroy(Wtmp)
      ! Old velocity just needs geometry
      call this%Uold%reset_level(lvl, ba, dm)
      call this%Vold%reset_level(lvl, ba, dm)
      call this%Wold%reset_level(lvl, ba, dm)
      ! Pressure remake
      call this%P%on_remake(this%P, lvl, time, ba, dm)
      ! Divergence just needs geometry
      call this%div%reset_level(lvl, ba, dm)
   end subroutine on_remake

   !> Override on_clear: delete level
   subroutine on_clear(this, lvl)
      implicit none
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%U%clear_level(lvl)
      call this%V%clear_level(lvl)
      call this%W%clear_level(lvl)
      call this%Uold%clear_level(lvl)
      call this%Vold%clear_level(lvl)
      call this%Wold%clear_level(lvl)
      call this%P%clear_level(lvl)
      call this%div%clear_level(lvl)
   end subroutine on_clear

   !> Override post_regrid: average down for C/F consistency
   subroutine post_regrid(this, lbase, time)
      implicit none
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      call this%average_down_velocity(lbase)
      call this%average_down_pressure(lbase)
      ! Rebuild pressure solver operators for new grid
      call this%psolver%setup()
   end subroutine post_regrid

   !> Average down MAC velocity for a single level (lvl+1 -> lvl)
   !> Uses amrdata infrastructure which handles face-centered averaging correctly
   subroutine average_down_velocity_to(this, lvl)
      implicit none
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%U%average_downto(lvl)
      call this%V%average_downto(lvl)
      call this%W%average_downto(lvl)
   end subroutine average_down_velocity_to

   !> Average down MAC velocity from finest to lbase
   !> Simply calls average_down_velocity_to in a loop
   !> @param lbase Optional: lowest level to average down to (default 0)
   subroutine average_down_velocity(this, lbase)
      implicit none
      class(amrincomp), intent(inout) :: this
      integer, intent(in), optional :: lbase
      integer :: lvl, lb
      lb = 0; if (present(lbase)) lb = lbase
      do lvl = this%amr%clvl()-1, lb, -1
         call this%average_down_velocity_to(lvl)
      end do
   end subroutine average_down_velocity

   !> Fill velocity ghost cells at a single level using divergence-free interpolation
   subroutine fill_velocity_lvl(this, lvl, time)
      use iso_c_binding, only: c_loc, c_funloc, c_funptr, c_ptr
      use amrex_interface, only: amrmfab_fillpatch_single, amrmfab_fillpatch_two_faces
      use amrdata_class, only: amrdata_fillbc
      implicit none
      class(amrincomp), target, intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr) :: ctx_u, ctx_v, ctx_w
      type(c_funptr) :: bc_dispatch
      integer :: rr, lo_bc(9), hi_bc(9)
      real(WP) :: t_old, t_new

      ! Get contexts for each velocity component
      ctx_u = c_loc(this%U); ctx_v = c_loc(this%V); ctx_w = c_loc(this%W)
      bc_dispatch = c_funloc(amrdata_fillbc)
      t_old = time - 1.0e200_WP
      t_new = time

      if (lvl .eq. 0) then
         ! Level 0: single-level fill (just physical BCs)
         call amrmfab_fillpatch_single(this%U%mf(0), time, this%U%mf(0), &
         &   time, this%U%mf(0), this%amr%geom(0), ctx_u, bc_dispatch, time, 1, 1, 1)
         call amrmfab_fillpatch_single(this%V%mf(0), time, this%V%mf(0), &
         &   time, this%V%mf(0), this%amr%geom(0), ctx_v, bc_dispatch, time, 1, 1, 1)
         call amrmfab_fillpatch_single(this%W%mf(0), time, this%W%mf(0), &
         &   time, this%W%mf(0), this%amr%geom(0), ctx_w, bc_dispatch, time, 1, 1, 1)
      else
         ! Build combined BC array: [U_x,U_y,U_z, V_x,V_y,V_z, W_x,W_y,W_z]
         lo_bc(1:3) = this%U%lo_bc(:,1)
         lo_bc(4:6) = this%V%lo_bc(:,1)
         lo_bc(7:9) = this%W%lo_bc(:,1)
         hi_bc(1:3) = this%U%hi_bc(:,1)
         hi_bc(4:6) = this%V%hi_bc(:,1)
         hi_bc(7:9) = this%W%hi_bc(:,1)
         rr = this%amr%rref(lvl-1)
         ! Call 3-component divfree FillPatch from two levels
         call amrmfab_fillpatch_two_faces( &
         &   this%U%mf(lvl), this%V%mf(lvl), this%W%mf(lvl), time, &
         &   t_old, this%U%mf(lvl-1), this%V%mf(lvl-1), this%W%mf(lvl-1), &
         &   t_new, this%U%mf(lvl-1), this%V%mf(lvl-1), this%W%mf(lvl-1), &
         &   this%amr%geom(lvl-1), &
         &   t_old, this%U%mf(lvl), this%V%mf(lvl), this%W%mf(lvl), &
         &   t_new, this%U%mf(lvl), this%V%mf(lvl), this%W%mf(lvl), &
         &   this%amr%geom(lvl), &
         &   ctx_u, ctx_v, ctx_w, bc_dispatch, bc_dispatch, bc_dispatch, &
         &   1, 1, 1, rr, this%interp_vel, lo_bc, hi_bc)
      end if
      ! Reconcile shared face values at box boundaries
      call this%U%mf(lvl)%override_sync(this%amr%geom(lvl))
      call this%V%mf(lvl)%override_sync(this%amr%geom(lvl))
      call this%W%mf(lvl)%override_sync(this%amr%geom(lvl))
   end subroutine fill_velocity_lvl

   !> Fill velocity ghost cells on all levels
   subroutine fill_velocity(this, time)
      implicit none
      class(amrincomp), intent(inout) :: this
      real(WP), intent(in) :: time
      integer :: lvl
      do lvl = 0, this%amr%clvl()
         call this%fill_velocity_lvl(lvl, time)
      end do
   end subroutine fill_velocity

   !> Fill new fine level velocity from coarse using divergence-free interpolation
   !> Used during MakeNewLevelFromCoarse (creation of new fine levels)
   subroutine fill_velocity_from_coarse(this, lvl, time)
      use iso_c_binding, only: c_loc, c_funloc, c_funptr, c_ptr
      use amrex_interface, only: amrmfab_fillcoarsepatch_faces
      use amrdata_class, only: amrdata_fillbc
      implicit none
      class(amrincomp), target, intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr) :: ctx_u, ctx_v, ctx_w
      type(c_funptr) :: bc_dispatch
      integer :: rr, lo_bc(9), hi_bc(9)
      ! Get contexts for each velocity component
      ctx_u = c_loc(this%U); ctx_v = c_loc(this%V); ctx_w = c_loc(this%W)
      bc_dispatch = c_funloc(amrdata_fillbc)
      ! Build combined BC array: [U_x,U_y,U_z, V_x,V_y,V_z, W_x,W_y,W_z]
      lo_bc(1:3) = this%U%lo_bc(:,1)
      lo_bc(4:6) = this%V%lo_bc(:,1)
      lo_bc(7:9) = this%W%lo_bc(:,1)
      hi_bc(1:3) = this%U%hi_bc(:,1)
      hi_bc(4:6) = this%V%hi_bc(:,1)
      hi_bc(7:9) = this%W%hi_bc(:,1)
      rr = this%amr%rref(lvl-1)
      ! Call 3-component divfree FillCoarsePatch
      call amrmfab_fillcoarsepatch_faces( &
      &   this%U%mf(lvl), this%V%mf(lvl), this%W%mf(lvl), time, &
      &   this%U%mf(lvl-1), this%V%mf(lvl-1), this%W%mf(lvl-1), &
      &   this%amr%geom(lvl-1), this%amr%geom(lvl), &
      &   ctx_u, ctx_v, ctx_w, bc_dispatch, bc_dispatch, bc_dispatch, &
      &   1, 1, 1, rr, this%interp_vel, lo_bc, hi_bc)
      ! Reconcile shared face values at box boundaries
      call this%U%mf(lvl)%override_sync(this%amr%geom(lvl))
      call this%V%mf(lvl)%override_sync(this%amr%geom(lvl))
      call this%W%mf(lvl)%override_sync(this%amr%geom(lvl))
   end subroutine fill_velocity_from_coarse

   !> Sync velocity ghost cells at a single level (lightweight, no C/F interpolation)
   subroutine sync_velocity_lvl(this, lvl)
      implicit none
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%U%sync_lvl(lvl)
      call this%V%sync_lvl(lvl)
      call this%W%sync_lvl(lvl)
   end subroutine sync_velocity_lvl

   !> Sync velocity ghost cells on all levels
   subroutine sync_velocity(this)
      implicit none
      class(amrincomp), intent(inout) :: this
      integer :: lvl
      do lvl = 0, this%amr%clvl()
         call this%sync_velocity_lvl(lvl)
      end do
   end subroutine sync_velocity

   !> Fill destination MultiFabs with velocity using divergence-free interpolation
   !> Used during regridding (on_remake) to fill new layout MultiFabs
   subroutine fill_velocity_mfab(this, Udest, Vdest, Wdest, lvl, time)
      use iso_c_binding, only: c_loc, c_funloc, c_funptr, c_ptr
      use amrex_interface, only: amrmfab_fillpatch_single, amrmfab_fillpatch_two_faces
      use amrdata_class, only: amrdata_fillbc
      implicit none
      class(amrincomp), target, intent(inout) :: this
      type(amrex_multifab), intent(inout) :: Udest, Vdest, Wdest
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr) :: ctx_u, ctx_v, ctx_w
      type(c_funptr) :: bc_dispatch
      integer :: rr, lo_bc(9), hi_bc(9)
      real(WP) :: t_old, t_new

      ! Get contexts for each velocity component
      ctx_u = c_loc(this%U); ctx_v = c_loc(this%V); ctx_w = c_loc(this%W)
      bc_dispatch = c_funloc(amrdata_fillbc)
      t_old = time - 1.0e200_WP
      t_new = time

      if (lvl .eq. 0) then
         ! Level 0: single-level fill (just physical BCs)
         call amrmfab_fillpatch_single(Udest, t_old, this%U%mf(0), &
         &   t_new, this%U%mf(0), this%amr%geom(0), ctx_u, bc_dispatch, time, 1, 1, 1)
         call amrmfab_fillpatch_single(Vdest, t_old, this%V%mf(0), &
         &   t_new, this%V%mf(0), this%amr%geom(0), ctx_v, bc_dispatch, time, 1, 1, 1)
         call amrmfab_fillpatch_single(Wdest, t_old, this%W%mf(0), &
         &   t_new, this%W%mf(0), this%amr%geom(0), ctx_w, bc_dispatch, time, 1, 1, 1)
      else
         ! Build combined BC array
         lo_bc(1:3) = this%U%lo_bc(:,1)
         lo_bc(4:6) = this%V%lo_bc(:,1)
         lo_bc(7:9) = this%W%lo_bc(:,1)
         hi_bc(1:3) = this%U%hi_bc(:,1)
         hi_bc(4:6) = this%V%hi_bc(:,1)
         hi_bc(7:9) = this%W%hi_bc(:,1)
         rr = this%amr%rref(lvl-1)
         ! Call 3-component divfree FillPatch
         call amrmfab_fillpatch_two_faces( &
         &   Udest, Vdest, Wdest, time, &
         &   t_old, this%U%mf(lvl-1), this%V%mf(lvl-1), this%W%mf(lvl-1), &
         &   t_new, this%U%mf(lvl-1), this%V%mf(lvl-1), this%W%mf(lvl-1), &
         &   this%amr%geom(lvl-1), &
         &   t_old, this%U%mf(lvl), this%V%mf(lvl), this%W%mf(lvl), &
         &   t_new, this%U%mf(lvl), this%V%mf(lvl), this%W%mf(lvl), &
         &   this%amr%geom(lvl), &
         &   ctx_u, ctx_v, ctx_w, bc_dispatch, bc_dispatch, bc_dispatch, &
         &   1, 1, 1, rr, this%interp_vel, lo_bc, hi_bc)
      end if
      ! Reconcile shared face values at box boundaries
      call Udest%override_sync(this%amr%geom(lvl))
      call Vdest%override_sync(this%amr%geom(lvl))
      call Wdest%override_sync(this%amr%geom(lvl))
   end subroutine fill_velocity_mfab

   !> Average down pressure for a single level (lvl+1 -> lvl)
   !> Uses amrdata infrastructure which handles cell-centered averaging correctly
   subroutine average_down_pressure_to(this, lvl)
      implicit none
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%P%average_downto(lvl)
   end subroutine average_down_pressure_to

   !> Average down pressure from finest to lbase
   !> Simply calls average_down_pressure_to in a loop
   !> @param lbase Optional: lowest level to average down to (default 0)
   subroutine average_down_pressure(this, lbase)
      implicit none
      class(amrincomp), intent(inout) :: this
      integer, intent(in), optional :: lbase
      integer :: lvl, lb
      lb = 0; if (present(lbase)) lb = lbase
      do lvl = this%amr%clvl()-1, lb, -1
         call this%average_down_pressure_to(lvl)
      end do
   end subroutine average_down_pressure

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Compute divergence of velocity into internal div field, update divmax
   subroutine get_div(this)
      use amrex_interface, only: amrmfab_compute_divergence
      implicit none
      class(amrincomp), intent(inout) :: this
      integer :: lvl
      ! Use our wrapper to amrex's 2nd order staggered divergence
      do lvl = 0, this%amr%clvl()
         call amrmfab_compute_divergence(this%div%mf(lvl), &
            this%U%mf(lvl), this%V%mf(lvl), this%W%mf(lvl), &
            this%amr%geom(lvl))
      end do
      ! Update divmax
      this%divmax = 0.0_WP
      do lvl = 0, this%amr%clvl()
         this%divmax = max(this%divmax, this%div%norm0(lvl=lvl))
      end do
   end subroutine get_div

   !> Compute pressure gradient into user-provided face amrdata
   subroutine get_pgrad(this, dPdx, dPdy, dPdz)
      implicit none
      class(amrincomp), intent(inout) :: this
      type(amrdata), intent(inout) :: dPdx, dPdy, dPdz
      integer :: lvl, i, j, k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pP, pGx, pGy, pGz
      real(WP) :: dxi, dyi, dzi

      do lvl = 0, this%amr%clvl()
         dxi = 1.0_WP / this%amr%dx(lvl)
         dyi = 1.0_WP / this%amr%dy(lvl)
         dzi = 1.0_WP / this%amr%dz(lvl)

         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            pP  => this%P%mf(lvl)%dataptr(mfi)
            pGx => dPdx%mf(lvl)%dataptr(mfi)
            pGy => dPdy%mf(lvl)%dataptr(mfi)
            pGz => dPdz%mf(lvl)%dataptr(mfi)

            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     pGx(i,j,k,1) = (pP(i,j,k,1) - pP(i-1,j,k,1)) * dxi
                     pGy(i,j,k,1) = (pP(i,j,k,1) - pP(i,j-1,k,1)) * dyi
                     pGz(i,j,k,1) = (pP(i,j,k,1) - pP(i,j,k-1,1)) * dzi
                  end do
               end do
            end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine get_pgrad

   !> Get solver information: min/max velocity, min/max pressure, divergence, momentum, TKE
   subroutine get_info(this)
      use amrex_amr_module, only: amrex_mfiter, amrex_box, amrex_imultifab, amrex_imultifab_build, amrex_imultifab_destroy
      use amrex_interface, only: amrmask_make_fine
      use parallel, only: MPI_REAL_WP
      use mpi_f08
      implicit none
      class(amrincomp), intent(inout) :: this
      integer :: lvl, i, j, k, ierr
      real(WP) :: dV, Uc, Vc, Wc
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_imultifab) :: mask
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU, pV, pW
      integer, dimension(:,:,:,:), contiguous, pointer :: pMask

      ! First compute divergence (this updates divmax)
      call this%get_div()

      ! Initialize min/max values
      this%Umax = -huge(1.0_WP)
      this%Vmax = -huge(1.0_WP)
      this%Wmax = -huge(1.0_WP)
      this%Pmax = -huge(1.0_WP)

      ! Loop over all levels for min/max
      do lvl = 0, this%amr%clvl()
         this%Umax = max(this%Umax, this%U%norm0(lvl=lvl))
         this%Vmax = max(this%Vmax, this%V%norm0(lvl=lvl))
         this%Wmax = max(this%Wmax, this%W%norm0(lvl=lvl))
         this%Pmax = max(this%Pmax, this%P%norm0(lvl=lvl))
      end do

      ! Momentum integrals (rho * U * dV, summed over cells at level 0)
      this%rhoUint = this%rho * this%U%get_sum(lvl=0) * this%amr%dx(0) * this%amr%dy(0) * this%amr%dz(0)
      this%rhoVint = this%rho * this%V%get_sum(lvl=0) * this%amr%dx(0) * this%amr%dy(0) * this%amr%dz(0)
      this%rhoWint = this%rho * this%W%get_sum(lvl=0) * this%amr%dx(0) * this%amr%dy(0) * this%amr%dz(0)

      ! Kinetic energy integral: 0.5 * rho * (Uc^2 + Vc^2 + Wc^2) * dV
      ! Uses composite integration with fine masking to avoid double-counting
      this%rhoKint = 0.0_WP
      do lvl = 0, this%amr%clvl()
         dV = this%amr%dx(lvl) * this%amr%dy(lvl) * this%amr%dz(lvl)

         ! Build fine mask for this level (if not finest)
         if (lvl .lt. this%amr%clvl()) then
            call amrex_imultifab_build(mask, this%amr%ba(lvl), this%amr%dm(lvl), 1, 0)
            call amrmask_make_fine(mask, this%amr%ba(lvl+1), [this%amr%rref(lvl), this%amr%rref(lvl), this%amr%rref(lvl)], 0, 1)
         end if

         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            pU => this%U%mf(lvl)%dataptr(mfi)
            pV => this%V%mf(lvl)%dataptr(mfi)
            pW => this%W%mf(lvl)%dataptr(mfi)
            if (lvl .lt. this%amr%clvl()) pMask => mask%dataptr(mfi)

            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     ! Skip cells covered by finer level
                     if (lvl .lt. this%amr%clvl()) then
                        if (pMask(i,j,k,1) .eq. 0) cycle
                     end if
                     ! Interpolate face velocities to cell center
                     Uc = 0.5_WP * (pU(i,j,k,1) + pU(i+1,j,k,1))
                     Vc = 0.5_WP * (pV(i,j,k,1) + pV(i,j+1,k,1))
                     Wc = 0.5_WP * (pW(i,j,k,1) + pW(i,j,k+1,1))
                     ! Accumulate kinetic energy
                     this%rhoKint = this%rhoKint + 0.5_WP * this%rho * (Uc**2 + Vc**2 + Wc**2) * dV
                  end do
               end do
            end do
         end do
         call this%amr%mfiter_destroy(mfi)

         if (lvl .lt. this%amr%clvl()) call amrex_imultifab_destroy(mask)
      end do

      ! Reduce across MPI ranks
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%rhoKint, 1, MPI_REAL_WP, MPI_SUM, this%amr%comm, ierr)
   end subroutine get_info

   !> Compute CFL numbers (convective and viscous)
   subroutine get_cfl(this, dt, cfl, cflc)
      implicit none
      class(amrincomp), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), intent(out) :: cfl
      real(WP), intent(out), optional :: cflc
      integer :: lvl
      real(WP) :: dx, Umax_lvl, Vmax_lvl, Wmax_lvl, visc_cfl
      ! Reset CFLs
      this%CFLc_x = 0.0_WP; this%CFLc_y = 0.0_WP; this%CFLc_z = 0.0_WP
      this%CFLv_x = 0.0_WP; this%CFLv_y = 0.0_WP; this%CFLv_z = 0.0_WP
      ! Compute CFL at each level (finest level determines dt)
      do lvl = 0, this%amr%clvl()
         Umax_lvl = this%U%norm0(lvl=lvl)
         Vmax_lvl = this%V%norm0(lvl=lvl)
         Wmax_lvl = this%W%norm0(lvl=lvl)
         ! Convective CFL
         this%CFLc_x = max(this%CFLc_x, dt * Umax_lvl / this%amr%dx(lvl))
         this%CFLc_y = max(this%CFLc_y, dt * Vmax_lvl / this%amr%dy(lvl))
         this%CFLc_z = max(this%CFLc_z, dt * Wmax_lvl / this%amr%dz(lvl))
         ! Viscous CFL (explicit stability: dt < dx^2 / (4*nu))
         dx = min(this%amr%dx(lvl), this%amr%dy(lvl), this%amr%dz(lvl))
         visc_cfl = 4.0_WP * this%visc * dt / (this%rho * dx**2)
         this%CFLv_x = max(this%CFLv_x, visc_cfl)
         this%CFLv_y = max(this%CFLv_y, visc_cfl)
         this%CFLv_z = max(this%CFLv_z, visc_cfl)
      end do
      ! Compute max overall CFL
      this%CFL = max(this%CFLc_x, this%CFLc_y, this%CFLc_z, this%CFLv_x, this%CFLv_y, this%CFLv_z)
      ! Return max overall CFL
      cfl = this%CFL
      ! Optionally return max convective CFL
      if (present(cflc)) cflc = max(this%CFLc_x, this%CFLc_y, this%CFLc_z)
   end subroutine get_cfl

   !> Print solver info to screen
   subroutine amrincomp_print(this)
      use messager, only: log
      use string, only: str_long
      implicit none
      class(amrincomp), intent(in) :: this
      character(len=str_long) :: message
      call log("Incompressible solver: "//trim(this%name))
      write(message,'("  rho = ",ES12.5,", visc = ",ES12.5)') this%rho, this%visc
      call log(trim(message))
      call log("  Grid: "//trim(this%amr%name))
   end subroutine amrincomp_print

   !> Register solver data for checkpoint
   subroutine register_checkpoint(this, io)
      use amrio_class, only: amrio
      implicit none
      class(amrincomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      call io%add_data(this%U, 'U')
      call io%add_data(this%V, 'V')
      call io%add_data(this%W, 'W')
      call io%add_data(this%P, 'P')
   end subroutine register_checkpoint

   !> Restore solver data from checkpoint
   subroutine restore_checkpoint(this, io, dirname)
      use amrio_class, only: amrio
      implicit none
      class(amrincomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      call io%read_data(dirname, this%U, 'U')
      call io%read_data(dirname, this%V, 'V')
      call io%read_data(dirname, this%W, 'W')
      call io%read_data(dirname, this%P, 'P')
   end subroutine restore_checkpoint

   !> Compute momentum advection RHS for all levels
   !> drhoUdt = -∇·(ρuu), drhoVdt = -∇·(ρuv), drhoWdt = -∇·(ρuw)
   !> Uses flux averaging at C/F interfaces for conservation
   subroutine get_dmomdt(this, U, V, W, drhoUdt, drhoVdt, drhoWdt, time)
      use amrex_amr_module, only: amrex_multifab, amrex_multifab_destroy, amrex_mfiter, amrex_box
      use amrex_interface,  only: amrmfab_average_down_cell, amrmfab_average_down_edge
      implicit none
      class(amrincomp), intent(inout) :: this
      class(amrdata), intent(inout) :: U, V, W                          !< Velocity state (face-centered)
      class(amrdata), intent(inout) :: drhoUdt, drhoVdt, drhoWdt        !< Output: momentum RHS (face-centered)
      real(WP), intent(in) :: time
      ! Flux MultiFabs (9 total: 3 CC, 2 xy-edge, 2 xz-edge, 2 yz-edge)
      type(amrex_multifab), dimension(0:this%amr%maxlvl) :: FUx, FUy, FUz
      type(amrex_multifab), dimension(0:this%amr%maxlvl) :: FVx, FVy, FVz
      type(amrex_multifab), dimension(0:this%amr%maxlvl) :: FWx, FWy, FWz
      type(amrex_multifab) :: Ufill, Vfill, Wfill
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      integer :: lvl, i, j, k
      real(WP) :: dxi, dyi, dzi
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU, pV, pW
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pFUx, pFUy, pFUz
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pFVx, pFVy, pFVz
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pFWx, pFWy, pFWz
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pdUdt, pdVdt, pdWdt

      ! ========================================================================
      ! Phase 1: Compute fluxes on all levels
      ! ========================================================================
      do lvl = 0, this%amr%clvl()
         dxi = 1.0_WP / this%amr%dx(lvl)
         dyi = 1.0_WP / this%amr%dy(lvl)
         dzi = 1.0_WP / this%amr%dz(lvl)

         ! Build momentum flux MultiFabs
         ! FUx, FVy, FWz: cell-centered (diagonal fluxes)
         call this%amr%mfab_build(lvl, FUx(lvl), ncomp=1, nover=1, atface=[.false.,.false.,.false.])
         call this%amr%mfab_build(lvl, FVy(lvl), ncomp=1, nover=1, atface=[.false.,.false.,.false.])
         call this%amr%mfab_build(lvl, FWz(lvl), ncomp=1, nover=1, atface=[.false.,.false.,.false.])
         ! FUy, FVx: xy-edge (cross-fluxes)
         call this%amr%mfab_build(lvl, FUy(lvl), ncomp=1, nover=0, atface=[.true. ,.true. ,.false.])
         call this%amr%mfab_build(lvl, FVx(lvl), ncomp=1, nover=0, atface=[.true. ,.true. ,.false.])
         ! FUz, FWx: xz-edge (cross-fluxes)
         call this%amr%mfab_build(lvl, FUz(lvl), ncomp=1, nover=0, atface=[.true. ,.false.,.true. ])
         call this%amr%mfab_build(lvl, FWx(lvl), ncomp=1, nover=0, atface=[.true. ,.false.,.true. ])
         ! FVz, FWy: yz-edge (cross-fluxes)
         call this%amr%mfab_build(lvl, FVz(lvl), ncomp=1, nover=0, atface=[.false.,.true. ,.true. ])
         call this%amr%mfab_build(lvl, FWy(lvl), ncomp=1, nover=0, atface=[.false.,.true. ,.true. ])

         ! Build velocity fills with ghost cells (need 1 ghost for interpolation)
         call this%amr%mfab_build(lvl, Ufill, ncomp=1, nover=1, atface=[.true. ,.false.,.false.])
         call this%amr%mfab_build(lvl, Vfill, ncomp=1, nover=1, atface=[.false.,.true. ,.false.])
         call this%amr%mfab_build(lvl, Wfill, ncomp=1, nover=1, atface=[.false.,.false.,.true. ])
         call this%fill_velocity_mfab(Udest=Ufill, Vdest=Vfill, Wdest=Wfill, lvl=lvl, time=time)

         ! MFIter loop: compute all 9 fluxes
         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()  ! Cell-centered tile

            ! Get pointers to velocity data
            pU => Ufill%dataptr(mfi)
            pV => Vfill%dataptr(mfi)
            pW => Wfill%dataptr(mfi)

            ! Get pointers to flux data
            pFUx => FUx(lvl)%dataptr(mfi)
            pFUy => FUy(lvl)%dataptr(mfi)
            pFUz => FUz(lvl)%dataptr(mfi)
            pFVx => FVx(lvl)%dataptr(mfi)
            pFVy => FVy(lvl)%dataptr(mfi)
            pFVz => FVz(lvl)%dataptr(mfi)
            pFWx => FWx(lvl)%dataptr(mfi)
            pFWy => FWy(lvl)%dataptr(mfi)
            pFWz => FWz(lvl)%dataptr(mfi)

            ! Diagonal fluxes
            do k = bx%lo(3)-1, bx%hi(3)+1; do j = bx%lo(2)-1, bx%hi(2)+1; do i = bx%lo(1)-1, bx%hi(1)+1
               pFUx(i,j,k,1) = - 0.25_WP * this%rho * sum(pU(i:i+1,j,k,1))**2 + 2.0_WP * this%visc * (pU(i+1,j,k,1) - pU(i,j,k,1)) * dxi
               pFVy(i,j,k,1) = - 0.25_WP * this%rho * sum(pV(i,j:j+1,k,1))**2 + 2.0_WP * this%visc * (pV(i,j+1,k,1) - pV(i,j,k,1)) * dyi
               pFWz(i,j,k,1) = - 0.25_WP * this%rho * sum(pW(i,j,k:k+1,1))**2 + 2.0_WP * this%visc * (pW(i,j,k+1,1) - pW(i,j,k,1)) * dzi
            end do; end do; end do

            ! Edge cross-fluxes with proper bounds for each type
            ! xy-edge (FUy, FVx): nodal in x,y; cell in z -> [lo,hi] in z; [lo,hi+1] in x,y
            do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2)+1; do i = bx%lo(1), bx%hi(1)+1
               pFUy(i,j,k,1) = - 0.25_WP * this%rho * sum(pV(i-1:i,j,k,1)) * sum(pU(i,j-1:j,k,1)) &
               &              + this%visc * ((pU(i,j,k,1) - pU(i,j-1,k,1)) * dyi + (pV(i,j,k,1) - pV(i-1,j,k,1)) * dxi)
               pFVx(i,j,k,1) = pFUy(i,j,k,1)
            end do; end do; end do
            ! yz-edge (FVz, FWy): nodal in y,z; cell in x -> [lo,hi] in x; [lo,hi+1] in y,z
            do k = bx%lo(3), bx%hi(3)+1; do j = bx%lo(2), bx%hi(2)+1; do i = bx%lo(1), bx%hi(1)
               pFVz(i,j,k,1) = - 0.25_WP * this%rho * sum(pW(i,j-1:j,k,1)) * sum(pV(i,j,k-1:k,1)) &
               &              + this%visc * ((pV(i,j,k,1) - pV(i,j,k-1,1)) * dzi + (pW(i,j,k,1) - pW(i,j-1,k,1)) * dyi)
               pFWy(i,j,k,1) = pFVz(i,j,k,1)
            end do; end do; end do
            ! zx-edge (FWx, FUz): nodal in z,x; cell in y -> [lo,hi] in y; [lo,hi+1] in z,x
            do k = bx%lo(3), bx%hi(3)+1; do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)+1
               pFWx(i,j,k,1) = - 0.25_WP * this%rho * sum(pU(i,j,k-1:k,1)) * sum(pW(i-1:i,j,k,1)) &
               &              + this%visc * ((pW(i,j,k,1) - pW(i-1,j,k,1)) * dxi + (pU(i,j,k,1) - pU(i,j,k-1,1)) * dzi)
               pFUz(i,j,k,1) = pFWx(i,j,k,1)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)

         ! Destroy velocity fills (no longer needed)
         call amrex_multifab_destroy(Ufill)
         call amrex_multifab_destroy(Vfill)
         call amrex_multifab_destroy(Wfill)
      end do

      ! ========================================================================
      ! Phase 2: Average down fluxes (fine -> coarse) for conservation
      ! ========================================================================
      do lvl = this%amr%clvl(), 1, -1
         ! Cell-centered fluxes
         call amrmfab_average_down_cell(fmf=FUx(lvl), cmf=FUx(lvl-1), rr=this%amr%rref(lvl-1), cgeom=this%amr%geom(lvl-1), ngcrse=0)
         call amrmfab_average_down_cell(fmf=FVy(lvl), cmf=FVy(lvl-1), rr=this%amr%rref(lvl-1), cgeom=this%amr%geom(lvl-1), ngcrse=0)
         call amrmfab_average_down_cell(fmf=FWz(lvl), cmf=FWz(lvl-1), rr=this%amr%rref(lvl-1), cgeom=this%amr%geom(lvl-1), ngcrse=0)
         ! Edge-centered fluxes
         call amrmfab_average_down_edge(fmf=FUy(lvl), cmf=FUy(lvl-1), rr=this%amr%rref(lvl-1), cgeom=this%amr%geom(lvl-1), ngcrse=0)
         call amrmfab_average_down_edge(fmf=FVx(lvl), cmf=FVx(lvl-1), rr=this%amr%rref(lvl-1), cgeom=this%amr%geom(lvl-1), ngcrse=0)
         call amrmfab_average_down_edge(fmf=FUz(lvl), cmf=FUz(lvl-1), rr=this%amr%rref(lvl-1), cgeom=this%amr%geom(lvl-1), ngcrse=0)
         call amrmfab_average_down_edge(fmf=FWx(lvl), cmf=FWx(lvl-1), rr=this%amr%rref(lvl-1), cgeom=this%amr%geom(lvl-1), ngcrse=0)
         call amrmfab_average_down_edge(fmf=FVz(lvl), cmf=FVz(lvl-1), rr=this%amr%rref(lvl-1), cgeom=this%amr%geom(lvl-1), ngcrse=0)
         call amrmfab_average_down_edge(fmf=FWy(lvl), cmf=FWy(lvl-1), rr=this%amr%rref(lvl-1), cgeom=this%amr%geom(lvl-1), ngcrse=0)
      end do

      ! ========================================================================
      ! Phase 3: Compute divergence to get momentum RHS
      ! ========================================================================
      do lvl = 0, this%amr%clvl()
         dxi = 1.0_WP / this%amr%dx(lvl)
         dyi = 1.0_WP / this%amr%dy(lvl)
         dzi = 1.0_WP / this%amr%dz(lvl)

         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            ! Get pointers to flux data
            pFUx => FUx(lvl)%dataptr(mfi)
            pFUy => FUy(lvl)%dataptr(mfi)
            pFUz => FUz(lvl)%dataptr(mfi)
            pFVx => FVx(lvl)%dataptr(mfi)
            pFVy => FVy(lvl)%dataptr(mfi)
            pFVz => FVz(lvl)%dataptr(mfi)
            pFWx => FWx(lvl)%dataptr(mfi)
            pFWy => FWy(lvl)%dataptr(mfi)
            pFWz => FWz(lvl)%dataptr(mfi)

            ! Get pointers to output RHS
            pdUdt => drhoUdt%mf(lvl)%dataptr(mfi)
            pdVdt => drhoVdt%mf(lvl)%dataptr(mfi)
            pdWdt => drhoWdt%mf(lvl)%dataptr(mfi)

            ! U-momentum RHS at x-faces: -d(FUx)/dx - d(FUy)/dy - d(FUz)/dz
            bx = mfi%nodaltilebox(1)
            do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
               pdUdt(i,j,k,1) = dxi * (pFUx(i,j,k,1) - pFUx(i-1,j,k,1)) &
               &              + dyi * (pFUy(i,j+1,k,1) - pFUy(i,j,k,1)) &
               &              + dzi * (pFUz(i,j,k+1,1) - pFUz(i,j,k,1))
            end do; end do; end do

            ! V-momentum RHS at y-faces: -d(FVx)/dx - d(FVy)/dy - d(FVz)/dz
            bx = mfi%nodaltilebox(2)
            do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
               pdVdt(i,j,k,1) = dxi * (pFVx(i+1,j,k,1) - pFVx(i,j,k,1)) &
               &              + dyi * (pFVy(i,j,k,1) - pFVy(i,j-1,k,1)) &
               &              + dzi * (pFVz(i,j,k+1,1) - pFVz(i,j,k,1))
            end do; end do; end do

            ! W-momentum RHS at z-faces: -d(FWx)/dx - d(FWy)/dy - d(FWz)/dz
            bx = mfi%nodaltilebox(3)
            do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
               pdWdt(i,j,k,1) = dxi * (pFWx(i+1,j,k,1) - pFWx(i,j,k,1)) &
               &              + dyi * (pFWy(i,j+1,k,1) - pFWy(i,j,k,1)) &
               &              + dzi * (pFWz(i,j,k,1) - pFWz(i,j,k-1,1))
            end do; end do; end do

         end do
         call this%amr%mfiter_destroy(mfi)
      end do

      ! ========================================================================
      ! Cleanup flux MultiFabs
      ! ========================================================================
      do lvl = 0, this%amr%clvl()
         call amrex_multifab_destroy(FUx(lvl))
         call amrex_multifab_destroy(FUy(lvl))
         call amrex_multifab_destroy(FUz(lvl))
         call amrex_multifab_destroy(FVx(lvl))
         call amrex_multifab_destroy(FVy(lvl))
         call amrex_multifab_destroy(FVz(lvl))
         call amrex_multifab_destroy(FWx(lvl))
         call amrex_multifab_destroy(FWy(lvl))
         call amrex_multifab_destroy(FWz(lvl))
      end do

   end subroutine get_dmomdt

   !> Velocity boundary condition callback (shared by U, V, W)
   !> Handles staggering-aware BC fills for face-centered velocity data
   !> - ext_dir: calls user_dirichlet callback for Dirichlet values
   !> - foextrap: copies from interior (Neumann, zero gradient)
   !> - hoextrap: higher-order extrapolation
   !> - reflect_even/odd: symmetry/anti-symmetry
   subroutine velocity_fillbc(this, mf, scomp, ncomp, time, geom)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, &
      &   amrex_bc_ext_dir, amrex_bc_foextrap, amrex_bc_hoextrap, amrex_bc_reflect_even, amrex_bc_reflect_odd
      implicit none
      class(amrdata), intent(inout) :: this
      type(amrex_multifab), intent(inout) :: mf
      integer, intent(in) :: scomp, ncomp
      real(WP), intent(in) :: time
      type(amrex_geometry), intent(in) :: geom
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      class(amrincomp), pointer :: solver
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer :: i, j, k, dlo(3), dhi(3), flo(3), fhi(3)
      integer :: ilo, ihi, jlo, jhi, klo, khi
      character(len=1) :: comp
      logical :: xper, yper, zper

      ! Access solver for periodicity and user callback
      select type (s => this%parent)
       class is (amrincomp)
         solver => s
         xper = solver%amr%xper; yper = solver%amr%yper; zper = solver%amr%zper
         if (xper .and. yper .and. zper) return
      end select

      ! Get domain bounds and component ID
      dlo = geom%domain%lo
      dhi = geom%domain%hi
      comp = this%name(1:1)  ! 'U', 'V', or 'W'

      ! Compute staggered domain bounds for this component
      ! Cell domain: [dlo, dhi]. Face domains extend by 1 in staggered direction.
      flo = dlo; fhi = dhi
      select case (comp)
       case ('U'); fhi(1) = dhi(1) + 1
       case ('V'); fhi(2) = dhi(2) + 1
       case ('W'); fhi(3) = dhi(3) + 1
      end select

      ! Loop over FABs
      call amrex_mfiter_build(mfi, mf, tiling=.false.)
      do while (mfi%next())
         p => mf%dataptr(mfi)
         ilo = lbound(p,1); ihi = ubound(p,1)
         jlo = lbound(p,2); jhi = ubound(p,2)
         klo = lbound(p,3); khi = ubound(p,3)

         ! Skip if FAB entirely within staggered domain
         if (ilo .ge. flo(1) .and. ihi .le. fhi(1) .and. &
         &   jlo .ge. flo(2) .and. jhi .le. fhi(2) .and. &
         &   klo .ge. flo(3) .and. khi .le. fhi(3)) cycle

         select type (solver => this%parent)
          class is (amrincomp)
            ! X-LOW BOUNDARY
            if (.not.xper .and. ilo .lt. flo(1)) call apply_vel_bc_lo(dir=1, bnd=flo(1), bctype=this%lo_bc(1,1), face=1)
            ! X-HIGH BOUNDARY
            if (.not.xper .and. ihi .gt. fhi(1)) call apply_vel_bc_hi(dir=1, bnd=fhi(1), bctype=this%hi_bc(1,1), face=2)
            ! Y-LOW BOUNDARY
            if (.not.yper .and. jlo .lt. flo(2)) call apply_vel_bc_lo(dir=2, bnd=flo(2), bctype=this%lo_bc(2,1), face=3)
            ! Y-HIGH BOUNDARY
            if (.not.yper .and. jhi .gt. fhi(2)) call apply_vel_bc_hi(dir=2, bnd=fhi(2), bctype=this%hi_bc(2,1), face=4)
            ! Z-LOW BOUNDARY
            if (.not.zper .and. klo .lt. flo(3)) call apply_vel_bc_lo(dir=3, bnd=flo(3), bctype=this%lo_bc(3,1), face=5)
            ! Z-HIGH BOUNDARY
            if (.not.zper .and. khi .gt. fhi(3)) call apply_vel_bc_hi(dir=3, bnd=fhi(3), bctype=this%hi_bc(3,1), face=6)
         end select
      end do
      call amrex_mfiter_destroy(mfi)

   contains

      !> Apply BC at low boundary in direction dir
      !> For NORMAL component (e.g., U in x): fills boundary face + ghosts
      !> For TANGENT component (e.g., V in x): fills ghosts only
      subroutine apply_vel_bc_lo(dir, bnd, bctype, face)
         implicit none
         integer, intent(in) :: dir, bnd, bctype, face
         type(amrex_box) :: bc_bx
         integer :: ii, jj, kk, fill_to, src_from
         logical :: is_normal

         ! Check if this is a normal component (U in x, V in y, W in z)
         is_normal = (comp .eq. 'U' .and. dir .eq. 1) .or. &
         &           (comp .eq. 'V' .and. dir .eq. 2) .or. &
         &           (comp .eq. 'W' .and. dir .eq. 3)

         ! For normal: fill up to bnd (includes boundary face), source from bnd+1
         ! For tangent: fill up to bnd-1 (ghosts only), source from bnd
         if (is_normal) then
            fill_to = bnd
            src_from = bnd + 1
         else
            fill_to = bnd - 1
            src_from = bnd
         end if

         select case (bctype)
          case (amrex_bc_ext_dir)
            ! Dirichlet: fill via user callback
            if (associated(solver%user_dirichlet)) then
               select case (dir)
                case (1); bc_bx = amrex_box([ilo, jlo, klo], [fill_to, jhi, khi])
                case (2); bc_bx = amrex_box([ilo, jlo, klo], [ihi, fill_to, khi])
                case (3); bc_bx = amrex_box([ilo, jlo, klo], [ihi, jhi, fill_to])
               end select
               call solver%user_dirichlet(solver, bc_bx, p, comp, face, time, geom)
            else
               ! Default to zero if no callback provided
               select case (dir)
                case (1); do kk=klo,khi; do jj=jlo,jhi; do ii=ilo,fill_to; p(ii,jj,kk,1)=0.0_WP; end do; end do; end do
                case (2); do kk=klo,khi; do jj=jlo,fill_to; do ii=ilo,ihi; p(ii,jj,kk,1)=0.0_WP; end do; end do; end do
                case (3); do kk=klo,fill_to; do jj=jlo,jhi; do ii=ilo,ihi; p(ii,jj,kk,1)=0.0_WP; end do; end do; end do
               end select
            end if

          case (amrex_bc_foextrap)
            ! Neumann: copy from interior
            select case (dir)
             case (1); do kk=klo,khi; do jj=jlo,jhi; do ii=ilo,fill_to; p(ii,jj,kk,1)=p(src_from,jj,kk,1); end do; end do; end do
             case (2); do kk=klo,khi; do jj=jlo,fill_to; do ii=ilo,ihi; p(ii,jj,kk,1)=p(ii,src_from,kk,1); end do; end do; end do
             case (3); do kk=klo,fill_to; do jj=jlo,jhi; do ii=ilo,ihi; p(ii,jj,kk,1)=p(ii,jj,src_from,1); end do; end do; end do
            end select

          case (amrex_bc_hoextrap)
            ! Higher-order extrapolation
            select case (dir)
             case (1); do kk=klo,khi; do jj=jlo,jhi; do ii=ilo,fill_to
                        p(ii,jj,kk,1) = 2.0_WP*p(src_from,jj,kk,1) - p(src_from+1,jj,kk,1)
                     end do; end do; end do
             case (2); do kk=klo,khi; do jj=jlo,fill_to; do ii=ilo,ihi
                        p(ii,jj,kk,1) = 2.0_WP*p(ii,src_from,kk,1) - p(ii,src_from+1,kk,1)
                     end do; end do; end do
             case (3); do kk=klo,fill_to; do jj=jlo,jhi; do ii=ilo,ihi
                        p(ii,jj,kk,1) = 2.0_WP*p(ii,jj,src_from,1) - p(ii,jj,src_from+1,1)
                     end do; end do; end do
            end select

          case (amrex_bc_reflect_even)
            ! Symmetry: F(-n) = F(n) - always ghosts only (bnd-1)
            select case (dir)
             case (1); do kk=klo,khi; do jj=jlo,jhi; do ii=ilo,bnd-1; p(ii,jj,kk,1)=p(2*bnd-ii-1,jj,kk,1); end do; end do; end do
             case (2); do kk=klo,khi; do jj=jlo,bnd-1; do ii=ilo,ihi; p(ii,jj,kk,1)=p(ii,2*bnd-jj-1,kk,1); end do; end do; end do
             case (3); do kk=klo,bnd-1; do jj=jlo,jhi; do ii=ilo,ihi; p(ii,jj,kk,1)=p(ii,jj,2*bnd-kk-1,1); end do; end do; end do
            end select

          case (amrex_bc_reflect_odd)
            ! Anti-symmetry: F(-n) = -F(n) - always ghosts only (bnd-1)
            select case (dir)
             case (1); do kk=klo,khi; do jj=jlo,jhi; do ii=ilo,bnd-1; p(ii,jj,kk,1)=-p(2*bnd-ii-1,jj,kk,1); end do; end do; end do
             case (2); do kk=klo,khi; do jj=jlo,bnd-1; do ii=ilo,ihi; p(ii,jj,kk,1)=-p(ii,2*bnd-jj-1,kk,1); end do; end do; end do
             case (3); do kk=klo,bnd-1; do jj=jlo,jhi; do ii=ilo,ihi; p(ii,jj,kk,1)=-p(ii,jj,2*bnd-kk-1,1); end do; end do; end do
            end select

         end select
      end subroutine apply_vel_bc_lo

      !> Apply BC at high boundary in direction dir
      !> For NORMAL component (e.g., U in x): fills boundary face + ghosts
      !> For TANGENT component (e.g., V in x): fills ghosts only
      subroutine apply_vel_bc_hi(dir, bnd, bctype, face)
         implicit none
         integer, intent(in) :: dir, bnd, bctype, face
         type(amrex_box) :: bc_bx
         integer :: ii, jj, kk, fill_from, src_from
         logical :: is_normal

         ! Check if this is a normal component (U in x, V in y, W in z)
         is_normal = (comp .eq. 'U' .and. dir .eq. 1) .or. &
         &           (comp .eq. 'V' .and. dir .eq. 2) .or. &
         &           (comp .eq. 'W' .and. dir .eq. 3)

         ! For normal: fill from bnd (includes boundary face), source from bnd-1
         ! For tangent: fill from bnd+1 (ghosts only), source from bnd
         if (is_normal) then
            fill_from = bnd
            src_from = bnd - 1
         else
            fill_from = bnd + 1
            src_from = bnd
         end if

         select case (bctype)
          case (amrex_bc_ext_dir)
            ! Dirichlet: fill via user callback
            if (associated(solver%user_dirichlet)) then
               select case (dir)
                case (1); bc_bx = amrex_box([fill_from, jlo, klo], [ihi, jhi, khi])
                case (2); bc_bx = amrex_box([ilo, fill_from, klo], [ihi, jhi, khi])
                case (3); bc_bx = amrex_box([ilo, jlo, fill_from], [ihi, jhi, khi])
               end select
               call solver%user_dirichlet(solver, bc_bx, p, comp, face, time, geom)
            else
               select case (dir)
                case (1); do kk=klo,khi; do jj=jlo,jhi; do ii=fill_from,ihi; p(ii,jj,kk,1)=0.0_WP; end do; end do; end do
                case (2); do kk=klo,khi; do jj=fill_from,jhi; do ii=ilo,ihi; p(ii,jj,kk,1)=0.0_WP; end do; end do; end do
                case (3); do kk=fill_from,khi; do jj=jlo,jhi; do ii=ilo,ihi; p(ii,jj,kk,1)=0.0_WP; end do; end do; end do
               end select
            end if

          case (amrex_bc_foextrap)
            ! Neumann: copy from interior
            select case (dir)
             case (1); do kk=klo,khi; do jj=jlo,jhi; do ii=fill_from,ihi; p(ii,jj,kk,1)=p(src_from,jj,kk,1); end do; end do; end do
             case (2); do kk=klo,khi; do jj=fill_from,jhi; do ii=ilo,ihi; p(ii,jj,kk,1)=p(ii,src_from,kk,1); end do; end do; end do
             case (3); do kk=fill_from,khi; do jj=jlo,jhi; do ii=ilo,ihi; p(ii,jj,kk,1)=p(ii,jj,src_from,1); end do; end do; end do
            end select

          case (amrex_bc_hoextrap)
            ! Higher-order extrapolation
            select case (dir)
             case (1); do kk=klo,khi; do jj=jlo,jhi; do ii=fill_from,ihi
                        p(ii,jj,kk,1) = 2.0_WP*p(src_from,jj,kk,1) - p(src_from-1,jj,kk,1)
                     end do; end do; end do
             case (2); do kk=klo,khi; do jj=fill_from,jhi; do ii=ilo,ihi
                        p(ii,jj,kk,1) = 2.0_WP*p(ii,src_from,kk,1) - p(ii,src_from-1,kk,1)
                     end do; end do; end do
             case (3); do kk=fill_from,khi; do jj=jlo,jhi; do ii=ilo,ihi
                        p(ii,jj,kk,1) = 2.0_WP*p(ii,jj,src_from,1) - p(ii,jj,src_from-1,1)
                     end do; end do; end do
            end select

          case (amrex_bc_reflect_even)
            ! Symmetry: F(-n) = F(n) - always ghosts only (bnd+1)
            select case (dir)
             case (1); do kk=klo,khi; do jj=jlo,jhi; do ii=bnd+1,ihi; p(ii,jj,kk,1)=p(2*bnd-ii+1,jj,kk,1); end do; end do; end do
             case (2); do kk=klo,khi; do jj=bnd+1,jhi; do ii=ilo,ihi; p(ii,jj,kk,1)=p(ii,2*bnd-jj+1,kk,1); end do; end do; end do
             case (3); do kk=bnd+1,khi; do jj=jlo,jhi; do ii=ilo,ihi; p(ii,jj,kk,1)=p(ii,jj,2*bnd-kk+1,1); end do; end do; end do
            end select

          case (amrex_bc_reflect_odd)
            ! Anti-symmetry: F(-n) = -F(n) - always ghosts only (bnd+1)
            select case (dir)
             case (1); do kk=klo,khi; do jj=jlo,jhi; do ii=bnd+1,ihi; p(ii,jj,kk,1)=-p(2*bnd-ii+1,jj,kk,1); end do; end do; end do
             case (2); do kk=klo,khi; do jj=bnd+1,jhi; do ii=ilo,ihi; p(ii,jj,kk,1)=-p(ii,2*bnd-jj+1,kk,1); end do; end do; end do
             case (3); do kk=bnd+1,khi; do jj=jlo,jhi; do ii=ilo,ihi; p(ii,jj,kk,1)=-p(ii,jj,2*bnd-kk+1,1); end do; end do; end do
            end select

         end select
      end subroutine apply_vel_bc_hi

   end subroutine velocity_fillbc

   !> Correct outflow velocity to ensure global mass conservation
   subroutine correct_outflow(this)
      use mpi_f08,  only: MPI_ALLREDUCE, MPI_SUM, MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(amrincomp), intent(inout) :: this
      real(WP) :: Qin, Qout, Ucorr, Aout
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU
      integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi
      integer :: dlo(3), dhi(3), ilo_bnd, ihi_bnd, ierr
      real(WP) :: dy, dz, dA
      type(amrex_mfiter) :: mfi

      if (this%amr%xper) return

      ! Level 0 only
      dlo = this%amr%geom(0)%domain%lo
      dhi = this%amr%geom(0)%domain%hi
      ilo_bnd = dlo(1)
      ihi_bnd = dhi(1) + 1
      dy = this%amr%dy(0)
      dz = this%amr%dz(0)
      dA = dy * dz

      Qin = 0.0_WP; Qout = 0.0_WP; Aout = 0.0_WP

      call this%amr%mfiter_build(0, mfi, tiling=.false.)
      do while (mfi%next())
         pU => this%U%mf(0)%dataptr(mfi)
         ilo = lbound(pU,1); ihi = ubound(pU,1)
         jlo = lbound(pU,2); jhi = ubound(pU,2)
         klo = lbound(pU,3); khi = ubound(pU,3)

         if (ilo .le. ilo_bnd .and. ilo_bnd .le. ihi) then
            do k = max(klo,dlo(3)), min(khi,dhi(3))
               do j = max(jlo,dlo(2)), min(jhi,dhi(2))
                  Qin = Qin + pU(ilo_bnd,j,k,1) * dA
               end do
            end do
         end if

         if (ilo .le. ihi_bnd .and. ihi_bnd .le. ihi) then
            do k = max(klo,dlo(3)), min(khi,dhi(3))
               do j = max(jlo,dlo(2)), min(jhi,dhi(2))
                  Qout = Qout + pU(ihi_bnd,j,k,1) * dA
                  Aout = Aout + dA
               end do
            end do
         end if
      end do
      call this%amr%mfiter_destroy(mfi)

      call MPI_ALLREDUCE(MPI_IN_PLACE, Qin,  1, MPI_REAL_WP, MPI_SUM, this%amr%comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, Qout, 1, MPI_REAL_WP, MPI_SUM, this%amr%comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, Aout, 1, MPI_REAL_WP, MPI_SUM, this%amr%comm, ierr)

      if (Aout .gt. 0.0_WP) then
         Ucorr = (Qin - Qout) / Aout
      else
         return
      end if

      call this%amr%mfiter_build(0, mfi, tiling=.false.)
      do while (mfi%next())
         pU => this%U%mf(0)%dataptr(mfi)
         ilo = lbound(pU,1); ihi = ubound(pU,1)
         jlo = lbound(pU,2); jhi = ubound(pU,2)
         klo = lbound(pU,3); khi = ubound(pU,3)

         if (ilo .le. ihi_bnd .and. ihi_bnd .le. ihi) then
            do k = max(klo,dlo(3)), min(khi,dhi(3))
               do j = max(jlo,dlo(2)), min(jhi,dhi(2))
                  pU(ihi_bnd,j,k,1) = pU(ihi_bnd,j,k,1) + Ucorr
               end do
            end do
         end if
      end do
      call this%amr%mfiter_destroy(mfi)
      
   end subroutine correct_outflow

end module amrincomp_class
