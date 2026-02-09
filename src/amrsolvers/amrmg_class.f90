!> Unified AMR Multigrid solver for elliptic problems
!> Supports both constant-coefficient (Poisson) and variable-coefficient problems.
!>   - cstcoef type: Uses amrex_poisson for ∇²ϕ = f
!>   - varcoef type: Uses amrex_abeclaplacian for (α·A - β·∇·(B∇))ϕ = f
module amrmg_class
   use precision,     only: WP
   use amrgrid_class, only: amrgrid
   use amrdata_class, only: amrdata
   use amrex_amr_module, only: amrex_multifab, amrex_geometry, amrex_boxarray, amrex_distromap
   use amrex_lo_bctypes_module
   use amrex_linear_solver_module
   implicit none
   private

   ! Solver type constants (required at init)
   integer, parameter, public :: amrmg_cstcoef = 1  !< Constant coefficient (amrex_poisson)
   integer, parameter, public :: amrmg_varcoef = 2  !< Variable coefficient (amrex_abeclaplacian)

   ! Bottom solver types
   integer, parameter, public :: amrmg_bottom_smoother = 0
   integer, parameter, public :: amrmg_bottom_bicgstab = 1
   integer, parameter, public :: amrmg_bottom_cg       = 2
   integer, parameter, public :: amrmg_bottom_hypre    = 3
   integer, parameter, public :: amrmg_bottom_petsc    = 4

   ! Expose AMReX boundary condition types for user convenience
   integer, parameter, public :: amrmg_bc_interior    = amrex_lo_bogus      !< Interior (not at domain boundary)
   integer, parameter, public :: amrmg_bc_dirichlet   = amrex_lo_dirichlet  !< Dirichlet BC
   integer, parameter, public :: amrmg_bc_neumann     = amrex_lo_neumann    !< Neumann BC
   integer, parameter, public :: amrmg_bc_reflect_odd = amrex_lo_reflect_odd
   integer, parameter, public :: amrmg_bc_periodic    = amrex_lo_periodic   !< Periodic BC

   !> Unified AMR Multigrid solver
   type, public :: amrmg
      ! Solver type (required at initialize)
      integer :: type = -1  !< amrmg_cstcoef or amrmg_varcoef (no default)

      ! Equation coefficients (varcoef type only)
      real(WP) :: alpha = 0.0_WP  !< Scalar for A term
      real(WP) :: beta  = 1.0_WP  !< Scalar for ∇·(B∇) term

      ! Solver parameters
      real(WP) :: tol_rel = 1.0e-10_WP
      real(WP) :: tol_abs = 0.0_WP
      integer  :: max_iter = 200
      integer  :: verbose = 0
      integer  :: bottom_solver = amrmg_bottom_bicgstab
      integer  :: maxorder = 2  !< BC/interpolation order

      ! Boundary conditions (set from grid periodicity at init)
      integer, dimension(3) :: lo_bc = amrmg_bc_interior
      integer, dimension(3) :: hi_bc = amrmg_bc_interior

      ! Output from solve
      real(WP) :: res   = 0.0_WP  !< Final residual
      integer  :: niter = 0       !< Number of iterations

      ! Private internals
      class(amrgrid), pointer, private :: amr => null()
      logical, private :: setup_done = .false.
      ! Stored AMReX objects for reuse
      type(amrex_poisson), private :: poisson
      type(amrex_abeclaplacian), private :: abeclap
      type(amrex_multigrid), private :: multigrid
   contains
      procedure :: initialize
      procedure :: setup
      procedure :: solve
      procedure :: solve_level
      procedure :: get_fluxes
      procedure :: destroy
      procedure :: print_short
      procedure :: print => print_long
      procedure :: log => log_solver
      procedure :: finalize
   end type amrmg

contains

   !> Initialize solver with grid and type
   subroutine initialize(this, amr, type)
      use messager, only: die, log
      class(amrmg), intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      integer, intent(in) :: type  !< amrmg_cstcoef or amrmg_varcoef (required)

      ! Validate type
      if (type .ne. amrmg_cstcoef .and. type .ne. amrmg_varcoef) then
         call die('[amrmg initialize] Invalid type: must be amrmg_cstcoef or amrmg_varcoef')
      end if

      this%amr => amr
      this%type = type
      this%setup_done = .false.

      ! Set smart BC defaults based on grid periodicity
      ! Periodic directions: amrmg_bc_periodic
      ! Non-periodic directions: amrmg_bc_neumann
      if (amr%xper) then
         this%lo_bc(1) = amrmg_bc_periodic
         this%hi_bc(1) = amrmg_bc_periodic
      else
         this%lo_bc(1) = amrmg_bc_neumann
         this%hi_bc(1) = amrmg_bc_neumann
      end if
      if (amr%yper) then
         this%lo_bc(2) = amrmg_bc_periodic
         this%hi_bc(2) = amrmg_bc_periodic
      else
         this%lo_bc(2) = amrmg_bc_neumann
         this%hi_bc(2) = amrmg_bc_neumann
      end if
      if (amr%zper) then
         this%lo_bc(3) = amrmg_bc_periodic
         this%hi_bc(3) = amrmg_bc_periodic
      else
         this%lo_bc(3) = amrmg_bc_neumann
         this%hi_bc(3) = amrmg_bc_neumann
      end if

      if (type .eq. amrmg_cstcoef) then
         call log('[amrmg] Initialized constant-coefficient solver')
      else
         call log('[amrmg] Initialized variable-coefficient solver')
      end if
   end subroutine initialize

   !> Setup solver - call after coefficients change or after regrid
   !> For varcoef type, bcoef fields must be provided
   subroutine setup(this, acoef, bcoef_x, bcoef_y, bcoef_z)
      use messager, only: die
      class(amrmg), intent(inout) :: this
      type(amrdata), intent(in), optional :: acoef
      type(amrdata), intent(in), optional :: bcoef_x, bcoef_y, bcoef_z

      type(amrex_geometry), dimension(:), allocatable :: geom
      type(amrex_boxarray), dimension(:), allocatable :: ba
      type(amrex_distromap), dimension(:), allocatable :: dm

      if (this%type .eq. -1) call die('[amrmg setup] Solver not initialized')
      if (this%setup_done) call this%destroy()

      ! Build AMReX wrapper arrays
      build_arrays: block
         integer :: lev
         allocate(geom(0:this%amr%clvl()))
         allocate(ba(0:this%amr%clvl()))
         allocate(dm(0:this%amr%clvl()))
         do lev = 0, this%amr%clvl()
            geom(lev) = this%amr%geom(lev)
            ba(lev) = this%amr%get_boxarray(lev)
            dm(lev) = this%amr%get_distromap(lev)
         end do
      end block build_arrays

      ! Build operator based on type
      select case (this%type)

       case (amrmg_cstcoef)
         poisson_setup: block
            call amrex_poisson_build(this%poisson, geom, ba, dm, &
               metric_term=.false., agglomeration=.true., consolidation=.true., &
               max_coarsening_level=30)
            call this%poisson%set_domain_bc(this%lo_bc, this%hi_bc)
            call amrex_multigrid_build(this%multigrid, this%poisson)
            call this%multigrid%set_verbose(this%verbose)
            call this%multigrid%set_max_iter(this%max_iter)
            call this%multigrid%set_bottom_solver(this%bottom_solver)
         end block poisson_setup

       case (amrmg_varcoef)
         abeclap_setup: block
            type(amrex_multifab) :: bcoef(3)
            integer :: lev
            call amrex_abeclaplacian_build(this%abeclap, geom, ba, dm, &
               metric_term=.false., agglomeration=.true., consolidation=.true., &
               max_coarsening_level=30)
            call this%abeclap%set_domain_bc(this%lo_bc, this%hi_bc)
            call this%abeclap%set_maxorder(this%maxorder)
            call this%abeclap%set_scalars(this%alpha, this%beta)
            ! Set A coefficient
            if (present(acoef)) then
               do lev = 0, this%amr%clvl()
                  call this%abeclap%set_acoeffs(lev, acoef%mf(lev))
               end do
            end if
            ! Set B coefficients
            if (present(bcoef_x) .and. present(bcoef_y) .and. present(bcoef_z)) then
               do lev = 0, this%amr%clvl()
                  bcoef(1) = bcoef_x%mf(lev)
                  bcoef(2) = bcoef_y%mf(lev)
                  bcoef(3) = bcoef_z%mf(lev)
                  call this%abeclap%set_bcoeffs(lev, bcoef)
               end do
            end if
            call amrex_multigrid_build(this%multigrid, this%abeclap)
            call this%multigrid%set_verbose(this%verbose)
            call this%multigrid%set_max_iter(this%max_iter)
            call this%multigrid%set_bottom_solver(this%bottom_solver)
         end block abeclap_setup

      end select

      this%setup_done = .true.
   end subroutine setup

   !> Solve - uses stored operator and multigrid
   !> @param phi Solution field (in: initial guess with BC in ghosts, out: solution)
   !> @param rhs Right-hand side field
   subroutine solve(this, phi, rhs)
      use messager, only: die
      class(amrmg), intent(inout) :: this
      type(amrdata), intent(inout) :: phi
      type(amrdata), intent(in) :: rhs

      type(amrex_multifab), dimension(:), allocatable :: sol, rhsmf

      if (this%type .eq. -1) call die('[amrmg solve] Solver not initialized')
      if (.not. this%setup_done) call die('[amrmg solve] Solver not setup')

      ! Build solution/rhs arrays
      build_arrays: block
         integer :: lev
         allocate(sol(0:this%amr%clvl()))
         allocate(rhsmf(0:this%amr%clvl()))
         do lev = 0, this%amr%clvl()
            sol(lev) = phi%mf(lev)
            rhsmf(lev) = rhs%mf(lev)
         end do
      end block build_arrays

      ! Set level BCs
      set_bcs: block
         integer :: lev
         select case (this%type)
          case (amrmg_cstcoef)
            do lev = 0, this%amr%clvl()
               call this%poisson%set_level_bc(lev, sol(lev))
            end do
          case (amrmg_varcoef)
            do lev = 0, this%amr%clvl()
               call this%abeclap%set_level_bc(lev, sol(lev))
            end do
         end select
      end block set_bcs

      ! Solve and get iteration count
      this%res = this%multigrid%solve(sol, rhsmf, this%tol_rel, this%tol_abs)
      get_niters: block
         use amrex_interface, only: amrmlmg_get_niters
         this%niter = amrmlmg_get_niters(this%multigrid%p)
      end block get_niters
   end subroutine solve

   !> Level-by-level solve (for subcycling)
   !> Builds single-level operator on-the-fly (not stored)
   !> @param lev Level to solve on (0-indexed, for geometry lookup)
   !> @param phi_mf Solution MultiFab (in: initial guess, out: solution)
   !> @param rhs_mf Right-hand side MultiFab
   !> @param phi_crse_mf Coarse level solution for C/F BC (required if lev>0)
   !> @param acoef_mf Optional cell-centered A coefficient (varcoef only)
   !> @param bcoef_x_mf,bcoef_y_mf,bcoef_z_mf Optional face-centered B coefficients (varcoef only)
   subroutine solve_level(this, lev, phi_mf, rhs_mf, phi_crse_mf, acoef_mf, bcoef_x_mf, bcoef_y_mf, bcoef_z_mf)
      use messager, only: die, log
      use string, only: str_long
      class(amrmg), intent(inout) :: this
      integer, intent(in) :: lev
      type(amrex_multifab), intent(inout) :: phi_mf
      type(amrex_multifab), intent(in) :: rhs_mf
      type(amrex_multifab), intent(in), optional :: phi_crse_mf
      type(amrex_multifab), intent(in), optional :: acoef_mf
      type(amrex_multifab), intent(in), optional :: bcoef_x_mf, bcoef_y_mf, bcoef_z_mf

      ! Single-level arrays for operator (always size 1, index 0)
      type(amrex_geometry) :: geom(0:0)
      type(amrex_boxarray) :: ba(0:0)
      type(amrex_distromap) :: dm(0:0)
      type(amrex_multifab) :: sol(0:0), rhsmf(0:0)
      integer :: rref

      ! Validate inputs
      if (this%type .eq. -1) call die('[amrmg solve_level] Solver not initialized')
      if (lev .gt. 0 .and. .not.present(phi_crse_mf)) call die('[amrmg solve_level] phi_crse_mf required for lev>0')

      ! Build single-level arrays
      setup_arrays: block
         geom(0) = this%amr%geom(lev)
         ba(0) = this%amr%get_boxarray(lev)
         dm(0) = this%amr%get_distromap(lev)
         sol(0) = phi_mf
         rhsmf(0) = rhs_mf
         if (lev .gt. 0) rref = this%amr%rref(lev-1)
      end block setup_arrays

      ! Solve based on operator type
      select case (this%type)

       case (amrmg_cstcoef)
         poisson_solve: block
            type(amrex_poisson) :: linop
            type(amrex_multigrid) :: mlmg
            ! Build operator
            call amrex_poisson_build(linop, geom, ba, dm, &
               metric_term=.false., agglomeration=.true., consolidation=.true., &
               max_coarsening_level=30)
            call linop%set_domain_bc(this%lo_bc, this%hi_bc)
            ! Set C/F BC if on refined level
            if (lev .gt. 0) call linop%set_coarse_fine_bc(phi_crse_mf, rref)
            call linop%set_level_bc(0, sol(0))
            ! Solve
            call amrex_multigrid_build(mlmg, linop)
            call mlmg%set_verbose(this%verbose)
            call mlmg%set_max_iter(this%max_iter)
            call mlmg%set_bottom_solver(this%bottom_solver)
            this%res = mlmg%solve(sol, rhsmf, this%tol_rel, this%tol_abs)
            ! Get iteration count before cleanup
            get_niters_p: block
               use amrex_interface, only: amrmlmg_get_niters
               this%niter = amrmlmg_get_niters(mlmg%p)
            end block get_niters_p
            call amrex_multigrid_destroy(mlmg)
            call amrex_poisson_destroy(linop)
         end block poisson_solve

       case (amrmg_varcoef)
         abeclap_solve: block
            type(amrex_abeclaplacian) :: linop
            type(amrex_multigrid) :: mlmg
            type(amrex_multifab) :: bcoef(3)
            ! Build operator
            call amrex_abeclaplacian_build(linop, geom, ba, dm, &
               metric_term=.false., agglomeration=.true., consolidation=.true., &
               max_coarsening_level=30)
            call linop%set_domain_bc(this%lo_bc, this%hi_bc)
            call linop%set_maxorder(this%maxorder)
            call linop%set_scalars(this%alpha, this%beta)
            ! Set coefficients
            if (present(acoef_mf)) call linop%set_acoeffs(0, acoef_mf)
            if (present(bcoef_x_mf) .and. present(bcoef_y_mf) .and. present(bcoef_z_mf)) then
               bcoef(1) = bcoef_x_mf
               bcoef(2) = bcoef_y_mf
               bcoef(3) = bcoef_z_mf
               call linop%set_bcoeffs(0, bcoef)
            end if
            ! Set C/F BC if on refined level
            if (lev .gt. 0) call linop%set_coarse_fine_bc(phi_crse_mf, rref)
            call linop%set_level_bc(0, sol(0))
            ! Solve
            call amrex_multigrid_build(mlmg, linop)
            call mlmg%set_verbose(this%verbose)
            call mlmg%set_max_iter(this%max_iter)
            call mlmg%set_bottom_solver(this%bottom_solver)
            this%res = mlmg%solve(sol, rhsmf, this%tol_rel, this%tol_abs)
            ! Get iteration count before cleanup
            get_niters_a: block
               use amrex_interface, only: amrmlmg_get_niters
               this%niter = amrmlmg_get_niters(mlmg%p)
            end block get_niters_a
            call amrex_multigrid_destroy(mlmg)
            call amrex_abeclaplacian_destroy(linop)
         end block abeclap_solve

      end select

      ! Log result
      log_result: block
         character(len=str_long) :: message
         if (this%amr%amRoot) then
            write(message,'("[amrmg solve_level] lev=",i0," res=",es12.5)') lev, this%res
            call log(message)
         end if
      end block log_result
   end subroutine solve_level

   !> Get face-centered fluxes using solver's C/F-consistent stencils
   !> For (alpha*A - beta*div(B*grad))phi = rhs, flux = -B*grad(phi)
   !> This uses the solver's internal quadratic interpolation which is
   !> consistent at C/F interfaces (unlike FillPatch's linear interpolation).
   !> Works for both cstcoef (Poisson) and varcoef (ABecLaplacian).
   !> @param phi Solution field
   !> @param flux_x X-face-centered flux output
   !> @param flux_y Y-face-centered flux output
   !> @param flux_z Z-face-centered flux output
   subroutine get_fluxes(this, phi, flux_x, flux_y, flux_z)
      use iso_c_binding, only: c_ptr
      use messager, only: die
      use amrex_interface, only: amrmlmg_get_fluxes
      class(amrmg), intent(in) :: this
      type(amrdata), intent(in) :: phi
      type(amrdata), intent(inout) :: flux_x, flux_y, flux_z

      type(c_ptr), dimension(:), allocatable :: sol_ptrs, fx_ptrs, fy_ptrs, fz_ptrs
      integer :: lev, nlevs

      if (.not. this%setup_done) call die('[amrmg get_fluxes] Solver not setup')

      nlevs = this%amr%clvl() + 1

      ! Build pointer arrays
      allocate(sol_ptrs(0:this%amr%clvl()))
      allocate(fx_ptrs(0:this%amr%clvl()))
      allocate(fy_ptrs(0:this%amr%clvl()))
      allocate(fz_ptrs(0:this%amr%clvl()))

      do lev = 0, this%amr%clvl()
         sol_ptrs(lev) = phi%mf(lev)%p
         fx_ptrs(lev) = flux_x%mf(lev)%p
         fy_ptrs(lev) = flux_y%mf(lev)%p
         fz_ptrs(lev) = flux_z%mf(lev)%p
      end do

      ! Call C++ wrapper - works for both Poisson and ABecLaplacian
      call amrmlmg_get_fluxes(this%multigrid%p, sol_ptrs, fx_ptrs, fy_ptrs, fz_ptrs, nlevs)

      ! Deallocate pointer arrays
      deallocate(sol_ptrs)
      deallocate(fx_ptrs)
      deallocate(fy_ptrs)
      deallocate(fz_ptrs)
   end subroutine get_fluxes

   !> Destroy solver internals - call before setup when operator changes
   subroutine destroy(this)
      class(amrmg), intent(inout) :: this
      if (this%setup_done) then
         call amrex_multigrid_destroy(this%multigrid)
         select case (this%type)
          case (amrmg_cstcoef)
            call amrex_poisson_destroy(this%poisson)
          case (amrmg_varcoef)
            call amrex_abeclaplacian_destroy(this%abeclap)
         end select
      end if
      this%setup_done = .false.
   end subroutine destroy

   !> Log solver info
   subroutine log_solver(this)
      use string,   only: str_long
      use messager, only: log
      class(amrmg), intent(in) :: this
      character(len=str_long) :: message
      if (this%amr%amRoot) then
         write(message,'("AMR Multigrid solver for grid [",a,"]")') trim(this%amr%name); call log(message)
         if (this%type .eq. amrmg_cstcoef) then
            write(message,'(" >     type = constant-coefficient (Poisson)")'); call log(message)
         else if (this%type .eq. amrmg_varcoef) then
            write(message,'(" >     type = variable-coefficient (ABecLaplacian)")'); call log(message)
            write(message,'(" >    alpha = ",es12.5)') this%alpha; call log(message)
            write(message,'(" >     beta = ",es12.5)') this%beta; call log(message)
         end if
         write(message,'(" >   it/max = ",i0,"/",i0)') this%niter,this%max_iter; call log(message)
         write(message,'(" >      res = ",es12.5)') this%res; call log(message)
         write(message,'(" >  tol_rel = ",es12.5)') this%tol_rel; call log(message)
      end if
   end subroutine log_solver

   !> Print solver info to screen
   subroutine print_long(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      class(amrmg), intent(in) :: this
      if (this%amr%amRoot) then
         write(output_unit,'("AMR Multigrid solver for grid [",a,"]")') trim(this%amr%name)
         if (this%type .eq. amrmg_cstcoef) then
            write(output_unit,'(" >     type = constant-coefficient (Poisson)")')
         else if (this%type .eq. amrmg_varcoef) then
            write(output_unit,'(" >     type = variable-coefficient (ABecLaplacian)")')
            write(output_unit,'(" >    alpha = ",es12.5)') this%alpha
            write(output_unit,'(" >     beta = ",es12.5)') this%beta
         end if
         write(output_unit,'(" >   it/max = ",i0,"/",i0)') this%niter,this%max_iter
         write(output_unit,'(" >      res = ",es12.5)') this%res
         write(output_unit,'(" >  tol_rel = ",es12.5)') this%tol_rel
      end if
   end subroutine print_long

   !> Short print of solver info
   subroutine print_short(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      class(amrmg), intent(in) :: this
      if (this%amr%amRoot) then
         write(output_unit,'("AMR MG solver [",a16,"] -> it/max = ",i3,"/",i3," res = ",es12.5)') &
            trim(this%amr%name),this%niter,this%max_iter,this%res
      end if
   end subroutine print_short

   !> Finalize and release resources
   subroutine finalize(this)
      class(amrmg), intent(inout) :: this
      if (this%setup_done) call this%destroy()
      nullify(this%amr)
      this%type = -1
   end subroutine finalize

end module amrmg_class
