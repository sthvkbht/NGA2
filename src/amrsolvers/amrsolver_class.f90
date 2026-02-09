!> Abstract base class for AMR solvers
!> Establishes callback pattern for internal callbacks (on_coarse, on_remake, on_clear, post_regrid)
!> User-overridable callbacks (on_init, tagging) are defined in concrete solvers with concrete types
module amrsolver_class
   use iso_c_binding,    only: c_ptr, c_null_ptr, c_loc, c_f_pointer, c_int, c_double
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid
   use amrex_amr_module, only: amrex_boxarray, amrex_distromap
   use amrio_class,      only: amrio
   implicit none
   private

   public :: amrsolver
   public :: amrsolver_on_init, amrsolver_on_coarse, amrsolver_on_remake
   public :: amrsolver_on_clear, amrsolver_tagging, amrsolver_postregrid

   !> Abstract solver base - extended by concrete solvers (amrscalar, etc.)
   type, abstract :: amrsolver
      ! Pointer to amrgrid
      class(amrgrid), pointer :: amr => null()
      ! Solver name for logging
      character(len=str_medium) :: name = 'UNNAMED_SOLVER'

   contains
      ! Internal callbacks (type-bound, implemented in base, overridden by concrete)
      procedure :: on_coarse
      procedure :: on_remake
      procedure :: on_clear
      procedure :: post_regrid

      ! Deferred - concrete solvers must implement
      procedure(info_iface), deferred :: get_info
      procedure(chkpt_iface), deferred :: register_checkpoint
      procedure(restore_iface), deferred :: restore_checkpoint
   end type amrsolver

   !> Interface for get_info
   abstract interface
      subroutine info_iface(this)
         import :: amrsolver
         class(amrsolver), intent(inout) :: this
      end subroutine info_iface
   end interface

   !> Interface for register_checkpoint
   abstract interface
      subroutine chkpt_iface(this, io)
         import :: amrsolver, amrio
         class(amrsolver), intent(inout) :: this
         class(amrio), intent(inout) :: io
      end subroutine chkpt_iface
   end interface

   !> Interface for restore_checkpoint
   abstract interface
      subroutine restore_iface(this, io, dirname)
         import :: amrsolver, amrio
         class(amrsolver), intent(inout) :: this
         class(amrio), intent(inout) :: io
         character(len=*), intent(in) :: dirname
      end subroutine restore_iface
   end interface

contains

   !> Internal on_coarse: called when new fine level created from coarse
   subroutine on_coarse(this, lvl, time, ba, dm)
      class(amrsolver), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Base implementation does nothing - concrete solvers override
   end subroutine on_coarse

   !> Internal on_remake: called when level is remade during regrid
   subroutine on_remake(this, lvl, time, ba, dm)
      class(amrsolver), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Base implementation does nothing - concrete solvers override
   end subroutine on_remake

   !> Internal on_clear: called when level is deleted
   subroutine on_clear(this, lvl)
      class(amrsolver), intent(inout) :: this
      integer, intent(in) :: lvl
      ! Base implementation does nothing - concrete solvers override
   end subroutine on_clear

   !> Internal post_regrid: called after regrid completes
   subroutine post_regrid(this, lbase, time)
      class(amrsolver), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      ! Base implementation does nothing - concrete solvers override
   end subroutine post_regrid


   ! ============================================================================
   ! DISPATCHERS - placeholder stubs, concrete solvers provide real dispatchers
   ! ============================================================================

   subroutine amrsolver_on_init(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
   end subroutine amrsolver_on_init

   subroutine amrsolver_on_coarse(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
   end subroutine amrsolver_on_coarse

   subroutine amrsolver_on_remake(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
   end subroutine amrsolver_on_remake

   subroutine amrsolver_on_clear(ctx, lvl)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
   end subroutine amrsolver_on_clear

   subroutine amrsolver_tagging(ctx, lvl, tags, time)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags
      real(WP), intent(in) :: time
   end subroutine amrsolver_tagging

   subroutine amrsolver_postregrid(ctx, lbase, time)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
   end subroutine amrsolver_postregrid

end module amrsolver_class
