!> AMR Compressible Multiphase solver class
!> Grown from amrcomp_class by embedding amrvof into it
module amrmpcomp_class
   use iso_c_binding,    only: c_ptr
   use precision,        only: WP
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata
   use amrsolver_class,  only: amrsolver
   use surfmesh_class,   only: surfmesh
   use amrex_amr_module, only: amrex_boxarray,amrex_distromap
   implicit none
   private

   ! VF clipping thresholds
   real(WP), parameter, public :: VFlo=1.0e-12_WP   !< Below this VF, liquid vanishes
   real(WP), parameter, public :: VFhi=1.0_WP-VFlo  !< Above this VF, gas vanishes
   real(WP), parameter, public :: vol_eps=1.0e-8_WP !< Volume epsilon for division by zero

   ! PLIC boundary condition types
   integer, parameter, public :: BC_LIQ     = 1  !< All liquid in ghost
   integer, parameter, public :: BC_GAS     = 2  !< All gas in ghost
   integer, parameter, public :: BC_REFLECT = 3  !< Symmetry (mirror across boundary)
   integer, parameter, public :: BC_USER    = 4  !< User-defined callback

   ! Expose type and dispatchers
   public :: amrmpcomp
   public :: amrmpcomp_on_init,amrmpcomp_on_coarse,amrmpcomp_on_remake
   public :: amrmpcomp_on_clear,amrmpcomp_tagging,amrmpcomp_postregrid

   !> AMR Compressible Multiphase solver type
   type, extends(amrsolver) :: amrmpcomp
      ! User-configurable callbacks
      procedure(comp_init_iface), pointer, nopass :: user_init=>null()
      procedure(comp_tagging_iface), pointer, nopass :: user_tagging=>null()
      procedure(comp_bc_iface), pointer, nopass :: user_bc=>null()
      procedure(vof_bc_iface), pointer, nopass :: user_vof_bc=>null()

      ! PLIC boundary conditions (per face, only used if direction is non-periodic)
      integer :: vof_lo_bc(3) = BC_REFLECT
      integer :: vof_hi_bc(3) = BC_REFLECT

      ! Liquid equation of state function pointers: PL=PL(rho,I), CL=CL(rho,P), TL=TL(rho,P)
      procedure(eos_P_iface), pointer, nopass :: getPL=>null()
      procedure(eos_C_iface), pointer, nopass :: getCL=>null()
      procedure(eos_T_iface), pointer, nopass :: getTL=>null()

      ! Gas equation of state function pointers: PG=PG(rho,I), CG=CG(rho,P), TG=TG(rho,P)
      procedure(eos_P_iface), pointer, nopass :: getPG=>null()
      procedure(eos_C_iface), pointer, nopass :: getCG=>null()
      procedure(eos_T_iface), pointer, nopass :: getTG=>null()

      ! Pointer to subroutine for mixture cell relaxation
      procedure(relax_iface), pointer, nopass :: relax=>null()

      ! Conserved variables (1: VF*rhoL, 2: (1-VF)*rhoG, 3: VF*rhoL*IL, 4: (1-VF)*rhoG*IG, 5: rhoU, 6: rhoV, 7: rhoW)
      type(amrdata) :: Q,Qold

      ! Volume fraction and barycenters (1st moments)
      type(amrdata) :: VF,VFold            !< Liquid volume fraction
      type(amrdata) :: Cliq,Cliqold        !< Liquid barycenter (3 components)
      type(amrdata) :: Cgas,Cgasold        !< Gas barycenter (3 components)

      ! PLIC interface (4 components: nx, ny, nz, d)
      type(amrdata) :: PLIC,PLICold

      ! Primitive variables: mixture velocity
      type(amrdata) :: U,V,W

      ! Phasic primitive variables
      type(amrdata) :: RHOL,RHOG           !< Phasic densities
      type(amrdata) :: IL,IG               !< Phasic internal energies
      type(amrdata) :: PL,PG               !< Phasic pressures
      type(amrdata) :: TL,TG               !< Phasic temperatures

      ! Mixture speed of sound
      type(amrdata) :: C

      ! Physical properties
      type(amrdata) :: visc              !< Dynamic viscosity
      type(amrdata) :: beta              !< Bulk viscosity
      type(amrdata) :: diff              !< Heat diffusivity

      ! CFL numbers
      real(WP) :: CFLc_x=0.0_WP,CFLc_y=0.0_WP,CFLc_z=0.0_WP  !< Convective
      real(WP) :: CFLa_x=0.0_WP,CFLa_y=0.0_WP,CFLa_z=0.0_WP  !< Acoustic
      real(WP) :: CFLv_x=0.0_WP,CFLv_y=0.0_WP,CFLv_z=0.0_WP  !< Viscous

      ! Monitoring quantities
      real(WP) :: Umax=0.0_WP,Vmax=0.0_WP,Wmax=0.0_WP
      real(WP) :: ILmin=0.0_WP,ILmax=0.0_WP,IGmin=0.0_WP,IGmax=0.0_WP
      real(WP) :: PLmin=0.0_WP,PLmax=0.0_WP,PGmin=0.0_WP,PGmax=0.0_WP
      real(WP) :: TLmin=0.0_WP,TLmax=0.0_WP,TGmin=0.0_WP,TGmax=0.0_WP
      real(WP) :: Cmin=0.0_WP,Cmax=0.0_WP
      real(WP), dimension(7) :: Qint=0.0_WP,Qmin=0.0_WP,Qmax=0.0_WP
      real(WP) :: rhoKint=0.0_WP
      real(WP) :: VFint=0.0_WP,VFmin=0.0_WP,VFmax=0.0_WP

      ! Minimum density for stability
      real(WP) :: rho_floor=1.0e-10_WP

      ! Number of overlap cells (2 for WENO3 stencil)
      integer :: nover=2

      ! Tagging parameter for VOF
      integer :: regrid_buffer=10  !< Number of cells to buffer around interface for tagging

      ! Surface mesh for visualization (finest level polygons)
      type(surfmesh) :: smesh

   contains
      procedure :: initialize
      procedure :: finalize
      ! Lifecycle callbacks
      procedure :: on_init
      procedure :: on_coarse
      procedure :: on_remake
      procedure :: on_clear
      procedure :: post_regrid
      ! Physics
      procedure :: get_primitive
      procedure :: get_conserved
      procedure :: get_dQdt
      procedure :: get_viscartif
      procedure :: get_cfl
      procedure :: build_plic
      procedure :: reset_moments
      procedure :: apply_relax
      ! VOF sync/fill/average utilities
      procedure :: fill_moments_lvl
      procedure :: sync_moments_lvl
      procedure :: sync_moments
      procedure :: sync_plic_lvl
      procedure :: sync_plic
      procedure :: fill_plic_lvl
      procedure :: vof_average_down
      ! Print
      procedure :: get_info
      procedure :: print=>amrmpcomp_print
      ! Checkpoint I/O (deferred from amrsolver)
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
   end type amrmpcomp

   !> Abstract interface for user init callback
   abstract interface
      subroutine comp_init_iface(solver,lvl,time,ba,dm)
         import :: amrmpcomp,WP,amrex_boxarray,amrex_distromap
         class(amrmpcomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine comp_init_iface
   end interface

   !> Abstract interface for tagging callback
   abstract interface
      subroutine comp_tagging_iface(solver,lvl,tags,time)
         import :: amrmpcomp,c_ptr,WP
         class(amrmpcomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         type(c_ptr), intent(in) :: tags
         real(WP), intent(in) :: time
      end subroutine comp_tagging_iface
   end interface

   !> Abstract interface for EoS: P=P(rho,I)
   abstract interface
      pure real(WP) function eos_P_iface(rho,I)
         import :: WP
         real(WP), intent(in) :: rho
         real(WP), intent(in) :: I
      end function eos_P_iface
   end interface

   !> Abstract interface for EoS: C=C(rho,P)
   abstract interface
      pure real(WP) function eos_C_iface(rho,P)
         import :: WP
         real(WP), intent(in) :: rho
         real(WP), intent(in) :: P
      end function eos_C_iface
   end interface

   !> Abstract interface for EoS: T=T(rho,P)
   abstract interface
      pure real(WP) function eos_T_iface(rho,P)
         import :: WP
         real(WP), intent(in) :: rho
         real(WP), intent(in) :: P
      end function eos_T_iface
   end interface

   !> Abstract interface for pressure relaxation callback
   abstract interface
      subroutine relax_iface(VF,Q)
         import :: WP
         real(WP), intent(inout) :: VF
         real(WP), dimension(:), intent(inout) :: Q
      end subroutine relax_iface
   end interface

   !> Abstract interface for user BC callback
   abstract interface
      subroutine comp_bc_iface(solver,pQ,bc_bx,face,time)
         use amrex_amr_module, only: amrex_box
         import :: amrmpcomp,WP
         class(amrmpcomp), intent(inout) :: solver
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
         type(amrex_box), intent(in) :: bc_bx
         integer, intent(in) :: face
         real(WP), intent(in) :: time
      end subroutine comp_bc_iface
   end interface

   !> Abstract interface for VOF boundary condition callback
   abstract interface
      subroutine vof_bc_iface(solver,bx,pVF,pCliq,pCgas,pPLIC,face,time,what)
         use amrex_amr_module, only: amrex_box
         import :: amrmpcomp,WP
         class(amrmpcomp), intent(inout) :: solver
         type(amrex_box), intent(in) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCliq,pCgas,pPLIC
         integer, intent(in) :: face      !< 1=xlo,2=xhi,3=ylo,4=yhi,5=zlo,6=zhi
         real(WP), intent(in) :: time
         integer, intent(in) :: what      !< 1=PLIC, 2=moments
      end subroutine vof_bc_iface
   end interface

contains

   ! ============================================================================
   ! DISPATCHERS
   ! ============================================================================

   !> Dispatch on_init: calls type-bound method then user callback
   subroutine amrmpcomp_on_init(ctx,lvl,time,ba,dm)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_init(lvl,time,ba,dm)
      if (associated(this%user_init)) call this%user_init(this,lvl,time,ba,dm)
   end subroutine amrmpcomp_on_init

   !> Dispatch on_coarse: calls type-bound method
   subroutine amrmpcomp_on_coarse(ctx,lvl,time,ba,dm)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_coarse(lvl,time,ba,dm)
   end subroutine amrmpcomp_on_coarse

   !> Dispatch on_remake: calls type-bound method
   subroutine amrmpcomp_on_remake(ctx,lvl,time,ba,dm)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_remake(lvl,time,ba,dm)
   end subroutine amrmpcomp_on_remake

   !> Dispatch on_clear: calls type-bound method
   subroutine amrmpcomp_on_clear(ctx,lvl)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_clear(lvl)
   end subroutine amrmpcomp_on_clear

   !> Dispatch tagging: VOF interface tagging + user callback
   subroutine amrmpcomp_tagging(ctx,lvl,tags,time)
      use amrex_amr_module, only: amrex_tagboxarray,amrex_mfiter,amrex_box,amrex_multifab
      use amrgrid_class, only: SETtag
      use iso_c_binding, only: c_f_pointer,c_char
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags
      real(WP), intent(in) :: time
      type(amrmpcomp), pointer :: this
      type(amrex_tagboxarray) :: tba
      type(amrex_multifab) :: band
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pBand
      integer :: i,j,k,dir,n,layer
      integer, dimension(3) :: ind
      integer :: eff_buffer
      
      call c_f_pointer(ctx,this)
      tba=tags
      
      ! Build band MultiFab with 1 ghost cell
      call this%amr%mfab_build(lvl=lvl,mfab=band,ncomp=1,nover=1)
      call band%setval(0.0_WP)
      
      ! Effective buffer (shrink at coarser levels)
      eff_buffer=max(1,this%regrid_buffer/(2**(this%amr%clvl()-lvl)))
      
      ! Pass 1: Mark interface cells (band=1)
      call this%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         bx=mfi%tilebox()
         pVF  =>this%VF%mf(lvl)%dataptr(mfi)
         pBand=>band%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Flag all obvious mixture cells
            if (pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi) then
               pBand(i,j,k,1)=1.0_WP
               cycle
            end if
            ! Check for implicit interfaces (pure cell adjacent to opposite phase)
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
      do layer=2,eff_buffer
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            bx=mfi%tilebox()
            pBand=>band%dataptr(mfi)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               if (pBand(i,j,k,1).eq.0.0_WP.and.any(pBand(i-1:i+1,j-1:j+1,k-1:k+1,1).eq.real(layer-1,WP))) pBand(i,j,k,1)=real(layer,WP)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         call band%fill_boundary(this%amr%geom(lvl))
      end do
      
      ! Pass 3: Set tags from band
      call this%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         bx=mfi%tilebox()
         tagarr=>tba%dataPtr(mfi)
         pBand =>band%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            if (pBand(i,j,k,1).gt.0.0_WP) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call this%amr%mfiter_destroy(mfi)
      
      ! Cleanup
      call this%amr%mfab_destroy(band)
      
      ! Call user tagging if provided
      if (associated(this%user_tagging)) call this%user_tagging(this,lvl,tags,time)
   end subroutine amrmpcomp_tagging

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrmpcomp_postregrid(ctx,lbase,time)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%post_regrid(lbase,time)
   end subroutine amrmpcomp_postregrid

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize the compressible solver
   subroutine initialize(this,amr,name)
      use iso_c_binding, only: c_loc
      use amrex_amr_module, only: amrex_bc_foextrap
      implicit none
      class(amrmpcomp), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name

      ! Set name
      if (present(name)) then
         this%name=trim(name)
      else
         this%name='UNNAMED_COMP'
      end if
      this%amr=>amr

      ! Initialize conserved variables Q (7 components for multiphase)
      call this%Q%initialize   (amr,name='Q'   ,ncomp=7,ng=this%nover); this%Q%parent   =>this
      call this%Qold%initialize(amr,name='Qold',ncomp=7,ng=this%nover); this%Qold%parent=>this

      ! Initialize volume fraction and barycenters
      call this%VF%initialize     (amr,name='VF'     ,ncomp=1,ng=this%nover); this%VF%parent     =>this
      call this%VFold%initialize  (amr,name='VFold'  ,ncomp=1,ng=this%nover); this%VFold%parent  =>this
      call this%Cliq%initialize   (amr,name='Cliq'   ,ncomp=3,ng=this%nover); this%Cliq%parent   =>this
      call this%Cliqold%initialize(amr,name='Cliqold',ncomp=3,ng=this%nover); this%Cliqold%parent=>this
      call this%Cgas%initialize   (amr,name='Cgas'   ,ncomp=3,ng=this%nover); this%Cgas%parent   =>this
      call this%Cgasold%initialize(amr,name='Cgasold',ncomp=3,ng=this%nover); this%Cgasold%parent=>this

      ! Initialize PLIC (4 components: nx, ny, nz, d)
      call this%PLIC%initialize   (amr,name='PLIC'   ,ncomp=4,ng=this%nover); this%PLIC%parent   =>this
      call this%PLICold%initialize(amr,name='PLICold',ncomp=4,ng=this%nover); this%PLICold%parent=>this

      ! Initialize mixture velocity
      call this%U%initialize(amr,name='U',ncomp=1,ng=this%nover); this%U%parent=>this
      call this%V%initialize(amr,name='V',ncomp=1,ng=this%nover); this%V%parent=>this
      call this%W%initialize(amr,name='W',ncomp=1,ng=this%nover); this%W%parent=>this

      ! Initialize phasic primitive variables
      call this%RHOL%initialize(amr,name='RHOL',ncomp=1,ng=this%nover); this%RHOL%parent=>this
      call this%RHOG%initialize(amr,name='RHOG',ncomp=1,ng=this%nover); this%RHOG%parent=>this
      call this%IL%initialize  (amr,name='IL'  ,ncomp=1,ng=this%nover); this%IL%parent  =>this
      call this%IG%initialize  (amr,name='IG'  ,ncomp=1,ng=this%nover); this%IG%parent  =>this
      call this%PL%initialize  (amr,name='PL'  ,ncomp=1,ng=this%nover); this%PL%parent  =>this
      call this%PG%initialize  (amr,name='PG'  ,ncomp=1,ng=this%nover); this%PG%parent  =>this
      call this%TL%initialize  (amr,name='TL'  ,ncomp=1,ng=this%nover); this%TL%parent  =>this
      call this%TG%initialize  (amr,name='TG'  ,ncomp=1,ng=this%nover); this%TG%parent  =>this
      call this%C%initialize   (amr,name='C'   ,ncomp=1,ng=this%nover); this%C%parent   =>this

      ! Initialize physical properties (Neumann BCs on those)
      call this%visc%initialize(amr,name='visc',ncomp=1,ng=this%nover); this%visc%parent=>this
      call this%beta%initialize(amr,name='beta',ncomp=1,ng=this%nover); this%beta%parent=>this
      call this%diff%initialize(amr,name='diff',ncomp=1,ng=this%nover); this%diff%parent=>this
      if (.not.amr%xper) then
         this%visc%lo_bc(1,1)=amrex_bc_foextrap; this%visc%hi_bc(1,1)=amrex_bc_foextrap
         this%beta%lo_bc(1,1)=amrex_bc_foextrap; this%beta%hi_bc(1,1)=amrex_bc_foextrap
         this%diff%lo_bc(1,1)=amrex_bc_foextrap; this%diff%hi_bc(1,1)=amrex_bc_foextrap
      end if
      if (.not.amr%yper) then
         this%visc%lo_bc(2,1)=amrex_bc_foextrap; this%visc%hi_bc(2,1)=amrex_bc_foextrap
         this%beta%lo_bc(2,1)=amrex_bc_foextrap; this%beta%hi_bc(2,1)=amrex_bc_foextrap
         this%diff%lo_bc(2,1)=amrex_bc_foextrap; this%diff%hi_bc(2,1)=amrex_bc_foextrap
      end if
      if (.not.amr%zper) then
         this%visc%lo_bc(3,1)=amrex_bc_foextrap; this%visc%hi_bc(3,1)=amrex_bc_foextrap
         this%beta%lo_bc(3,1)=amrex_bc_foextrap; this%beta%hi_bc(3,1)=amrex_bc_foextrap
         this%diff%lo_bc(3,1)=amrex_bc_foextrap; this%diff%hi_bc(3,1)=amrex_bc_foextrap
      end if

      ! Set Q fillbc callback to internal handler
      this%Q%fillbc=>Q_fillbc

      ! Register callbacks with amrgrid
      select type (this)
       type is (amrmpcomp)
         call this%amr%add_on_init   (amrmpcomp_on_init,   c_loc(this))
         call this%amr%add_on_coarse (amrmpcomp_on_coarse, c_loc(this))
         call this%amr%add_on_remake (amrmpcomp_on_remake, c_loc(this))
         call this%amr%add_on_clear  (amrmpcomp_on_clear,  c_loc(this))
         call this%amr%add_tagging   (amrmpcomp_tagging,   c_loc(this))
         call this%amr%add_postregrid(amrmpcomp_postregrid,c_loc(this))
      end select

      ! Print solver info
      call this%print()

   end subroutine initialize

   !> Finalize the compressible multiphase solver
   subroutine finalize(this)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      ! Conserved variables
      call this%Q%finalize(); call this%Qold%finalize()
      ! VOF moments
      call this%VF%finalize(); call this%VFold%finalize()
      call this%Cliq%finalize(); call this%Cliqold%finalize()
      call this%Cgas%finalize(); call this%Cgasold%finalize()
      ! PLIC
      call this%PLIC%finalize(); call this%PLICold%finalize()
      ! Velocity
      call this%U%finalize(); call this%V%finalize(); call this%W%finalize()
      ! Phasic primitives
      call this%RHOL%finalize(); call this%RHOG%finalize()
      call this%IL%finalize(); call this%IG%finalize()
      call this%PL%finalize(); call this%PG%finalize()
      call this%TL%finalize(); call this%TG%finalize()
      call this%C%finalize()
      ! Physical properties
      call this%visc%finalize(); call this%beta%finalize(); call this%diff%finalize()
      ! Nullify pointers
      nullify(this%amr); nullify(this%user_init); nullify(this%user_tagging)
      nullify(this%getPL); nullify(this%getCL); nullify(this%getTL)
      nullify(this%getPG); nullify(this%getCG); nullify(this%getTG)
      nullify(this%relax)
   end subroutine finalize

   ! ============================================================================
   ! INTERNAL CALLBACK OVERRIDES
   ! ============================================================================

   !> Override on_init: reset levels and set to zero
   subroutine on_init(this,lvl,time,ba,dm)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Reset level layouts for conserved variables
      call this%Q%reset_level(lvl,ba,dm); call this%Qold%reset_level(lvl,ba,dm)
      ! Reset VOF moments
      call this%VF%reset_level(lvl,ba,dm); call this%VFold%reset_level(lvl,ba,dm)
      call this%Cliq%reset_level(lvl,ba,dm); call this%Cliqold%reset_level(lvl,ba,dm)
      call this%Cgas%reset_level(lvl,ba,dm); call this%Cgasold%reset_level(lvl,ba,dm)
      ! Reset PLIC
      call this%PLIC%reset_level(lvl,ba,dm); call this%PLICold%reset_level(lvl,ba,dm)
      ! Reset mixture velocity
      call this%U%reset_level(lvl,ba,dm)
      call this%V%reset_level(lvl,ba,dm)
      call this%W%reset_level(lvl,ba,dm)
      ! Reset phasic primitives
      call this%RHOL%reset_level(lvl,ba,dm); call this%RHOG%reset_level(lvl,ba,dm)
      call this%IL%reset_level(lvl,ba,dm); call this%IG%reset_level(lvl,ba,dm)
      call this%PL%reset_level(lvl,ba,dm); call this%PG%reset_level(lvl,ba,dm)
      call this%TL%reset_level(lvl,ba,dm); call this%TG%reset_level(lvl,ba,dm)
      call this%C%reset_level(lvl,ba,dm)
      ! Reset physical properties
      call this%visc%reset_level(lvl,ba,dm)
      call this%beta%reset_level(lvl,ba,dm)
      call this%diff%reset_level(lvl,ba,dm)
      ! Zero out everything
      call this%Q%setval(val=0.0_WP,lvl=lvl); call this%Qold%setval(val=0.0_WP,lvl=lvl)
      call this%VF%setval(val=0.0_WP,lvl=lvl); call this%VFold%setval(val=0.0_WP,lvl=lvl)
      call this%Cliq%setval(val=0.0_WP,lvl=lvl); call this%Cliqold%setval(val=0.0_WP,lvl=lvl)
      call this%Cgas%setval(val=0.0_WP,lvl=lvl); call this%Cgasold%setval(val=0.0_WP,lvl=lvl)
      call this%PLIC%setval(val=0.0_WP,lvl=lvl); call this%PLICold%setval(val=0.0_WP,lvl=lvl)
      call this%U%setval(val=0.0_WP,lvl=lvl)
      call this%V%setval(val=0.0_WP,lvl=lvl)
      call this%W%setval(val=0.0_WP,lvl=lvl)
      call this%RHOL%setval(val=0.0_WP,lvl=lvl); call this%RHOG%setval(val=0.0_WP,lvl=lvl)
      call this%IL%setval(val=0.0_WP,lvl=lvl); call this%IG%setval(val=0.0_WP,lvl=lvl)
      call this%PL%setval(val=0.0_WP,lvl=lvl); call this%PG%setval(val=0.0_WP,lvl=lvl)
      call this%TL%setval(val=0.0_WP,lvl=lvl); call this%TG%setval(val=0.0_WP,lvl=lvl)
      call this%C%setval(val=0.0_WP,lvl=lvl)
      call this%visc%setval(val=0.0_WP,lvl=lvl)
      call this%beta%setval(val=0.0_WP,lvl=lvl)
      call this%diff%setval(val=0.0_WP,lvl=lvl)
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse using conservative interpolation
   subroutine on_coarse(this,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      
      ! --- Conserved variables Q ---
      call this%Q%on_coarse(this%Q,lvl,time,ba,dm)
      call this%Qold%reset_level(lvl,ba,dm)
      
      ! --- Auxiliary / derived quantities (just reset, will be recomputed) ---
      call this%U%reset_level(lvl,ba,dm)
      call this%V%reset_level(lvl,ba,dm)
      call this%W%reset_level(lvl,ba,dm)
      call this%RHOL%reset_level(lvl,ba,dm); call this%RHOG%reset_level(lvl,ba,dm)
      call this%IL%reset_level(lvl,ba,dm);   call this%IG%reset_level(lvl,ba,dm)
      call this%PL%reset_level(lvl,ba,dm);   call this%PG%reset_level(lvl,ba,dm)
      call this%TL%reset_level(lvl,ba,dm);   call this%TG%reset_level(lvl,ba,dm)
      call this%C%reset_level(lvl,ba,dm)
      call this%visc%reset_level(lvl,ba,dm)
      call this%beta%reset_level(lvl,ba,dm)
      call this%diff%reset_level(lvl,ba,dm)
      
      ! --- VOF moments: interpolate from coarse ---
      call this%VF%on_coarse(this%VF,lvl,time,ba,dm)
      call this%Cliq%on_coarse(this%Cliq,lvl,time,ba,dm)
      call this%Cgas%on_coarse(this%Cgas,lvl,time,ba,dm)
      call this%VFold%reset_level(lvl,ba,dm)
      call this%Cliqold%reset_level(lvl,ba,dm)
      call this%Cgasold%reset_level(lvl,ba,dm)
      
      ! --- PLIC: allocate and set trivial planes for pure cells ---
      call this%PLIC%reset_level(lvl,ba,dm)
      call this%PLICold%reset_level(lvl,ba,dm)
      call set_trivial_plic()
      
      ! --- Fill VOF ghosts (sync + periodic corrections + physical BC) ---
      call this%fill_moments_lvl(lvl,time)
      call this%fill_plic_lvl(lvl,time)
      
   contains
      !> Set PLIC to trivial planes based on VF
      subroutine set_trivial_plic()
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pPLIC
         integer :: i,j,k
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            bx=mfi%tilebox()
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            pPLIC=>this%PLIC%mf(lvl)%dataptr(mfi)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pPLIC(i,j,k,1:3)=0.0_WP
               pPLIC(i,j,k,4)=sign(1.0e10_WP,pVF(i,j,k,1)-0.5_WP)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end subroutine set_trivial_plic
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid using conservative interpolation
   subroutine on_remake(this,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      
      ! --- Conserved variables Q ---
      call this%Q%on_remake(this%Q,lvl,time,ba,dm)
      call this%Qold%reset_level(lvl,ba,dm)
      
      ! --- Auxiliary / derived quantities (just reset, will be recomputed) ---
      call this%U%reset_level(lvl,ba,dm)
      call this%V%reset_level(lvl,ba,dm)
      call this%W%reset_level(lvl,ba,dm)
      call this%RHOL%reset_level(lvl,ba,dm); call this%RHOG%reset_level(lvl,ba,dm)
      call this%IL%reset_level(lvl,ba,dm);   call this%IG%reset_level(lvl,ba,dm)
      call this%PL%reset_level(lvl,ba,dm);   call this%PG%reset_level(lvl,ba,dm)
      call this%TL%reset_level(lvl,ba,dm);   call this%TG%reset_level(lvl,ba,dm)
      call this%C%reset_level(lvl,ba,dm)
      call this%visc%reset_level(lvl,ba,dm)
      call this%beta%reset_level(lvl,ba,dm)
      call this%diff%reset_level(lvl,ba,dm)
      
      ! --- VOF moments: copy existing + fill new areas from coarse ---
      call this%VF%on_remake(this%VF,lvl,time,ba,dm)
      call this%Cliq%on_remake(this%Cliq,lvl,time,ba,dm)
      call this%Cgas%on_remake(this%Cgas,lvl,time,ba,dm)
      call this%VFold%reset_level(lvl,ba,dm)
      call this%Cliqold%reset_level(lvl,ba,dm)
      call this%Cgasold%reset_level(lvl,ba,dm)
      
      ! --- PLIC: copy existing + fill new areas, then fix pure cells ---
      call this%PLIC%on_remake(this%PLIC,lvl,time,ba,dm)
      call this%PLICold%reset_level(lvl,ba,dm)
      call fix_pure_plic()
      
      ! --- Fill VOF ghosts (sync + periodic corrections + physical BC) ---
      call this%fill_moments_lvl(lvl,time)
      call this%fill_plic_lvl(lvl,time)
      
   contains
      !> Fix PLIC for pure cells (new cells from coarse have wrong interpolated PLIC)
      subroutine fix_pure_plic()
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pPLIC
         integer :: i,j,k
         real(WP) :: vf
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            bx=mfi%tilebox()
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            pPLIC=>this%PLIC%mf(lvl)%dataptr(mfi)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               vf=pVF(i,j,k,1)
               if (vf.lt.VFlo.or.vf.gt.VFhi) then
                  pPLIC(i,j,k,1:3)=0.0_WP
                  pPLIC(i,j,k,4)=sign(1.0e10_WP,vf-0.5_WP)
               end if
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end subroutine fix_pure_plic
   end subroutine on_remake

   !> Override on_clear: delete level
   subroutine on_clear(this,lvl)
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%Q%clear_level(lvl); call this%Qold%clear_level(lvl)
      call this%VF%clear_level(lvl); call this%VFold%clear_level(lvl)
      call this%Cliq%clear_level(lvl); call this%Cliqold%clear_level(lvl)
      call this%Cgas%clear_level(lvl); call this%Cgasold%clear_level(lvl)
      call this%PLIC%clear_level(lvl); call this%PLICold%clear_level(lvl)
      call this%U%clear_level(lvl); call this%V%clear_level(lvl); call this%W%clear_level(lvl)
      call this%RHOL%clear_level(lvl); call this%RHOG%clear_level(lvl)
      call this%IL%clear_level(lvl); call this%IG%clear_level(lvl)
      call this%PL%clear_level(lvl); call this%PG%clear_level(lvl)
      call this%TL%clear_level(lvl); call this%TG%clear_level(lvl)
      call this%C%clear_level(lvl)
      call this%visc%clear_level(lvl)
      call this%beta%clear_level(lvl)
      call this%diff%clear_level(lvl)
   end subroutine on_clear

   !> Override post_regrid: average down for C/F consistency
   subroutine post_regrid(this,lbase,time)
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      integer :: lvl
      ! Average down VOF (VF/Cliq/Cgas + trivial PLIC cleanup at coarse levels)
      call this%vof_average_down(lbase)
      ! Average down conserved variables Q for C/F consistency
      do lvl=this%amr%clvl()-1,lbase,-1
         call this%Q%average_downto(lvl)
      end do
      ! Fill Q ghosts and rebuild primitives
      call this%Q%fill(time)
      call this%get_primitive(this%Q)
   end subroutine post_regrid

   !> Internal fillbc for Q - calls default_fillbc first, then user_bc for ext_dir faces
   subroutine Q_fillbc(this,mf,scomp,ncomp,time,geom)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,&
      &                           amrex_box,amrex_geometry,amrex_multifab,amrex_bc_ext_dir
      use amrdata_class, only: default_fillbc
      implicit none
      class(amrdata), intent(inout) :: this
      type(amrex_multifab), intent(inout) :: mf
      integer, intent(in) :: scomp, ncomp
      real(WP), intent(in) :: time
      type(amrex_geometry), intent(in) :: geom
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bc_bx
      class(amrmpcomp), pointer :: solver
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer :: ilo,ihi,jlo,jhi,klo,khi
      integer :: dlo(3),dhi(3)
      
      ! First apply default BC handling (foextrap,hoextrap,reflect,etc.)
      call default_fillbc(this,mf,scomp,ncomp,time,geom)
      
      ! Access parent solver
      select type (s=>this%parent)
       class is (amrmpcomp)
         solver=>s
      end select
      
      ! Check if user callback exists
      if (.not.associated(solver%user_bc)) return

      ! Get domain bounds
      dlo=geom%domain%lo
      dhi=geom%domain%hi

      ! Loop over FABs and apply user_bc for ext_dir faces
      call amrex_mfiter_build(mfi,mf,tiling=.false.)
      do while (mfi%next())
         p=>mf%dataptr(mfi)
         ilo=lbound(p,1); ihi=ubound(p,1)
         jlo=lbound(p,2); jhi=ubound(p,2)
         klo=lbound(p,3); khi=ubound(p,3)

         ! X-LOW (face=1)
         if (this%lo_bc(1,1).eq.amrex_bc_ext_dir .and. ilo.lt.dlo(1)) then
            bc_bx=amrex_box([ilo,jlo,klo],[dlo(1)-1,jhi,khi])
            call solver%user_bc(solver,p,bc_bx,1,time)
         end if

         ! X-HIGH (face=2)
         if (this%hi_bc(1,1).eq.amrex_bc_ext_dir .and. ihi.gt.dhi(1)) then
            bc_bx=amrex_box([dhi(1)+1,jlo,klo],[ihi,jhi,khi])
            call solver%user_bc(solver,p,bc_bx,2,time)
         end if

         ! Y-LOW (face=3)
         if (this%lo_bc(2,1).eq.amrex_bc_ext_dir .and. jlo.lt.dlo(2)) then
            bc_bx=amrex_box([ilo,jlo,klo],[ihi,dlo(2)-1,khi])
            call solver%user_bc(solver,p,bc_bx,3,time)
         end if

         ! Y-HIGH (face=4)
         if (this%hi_bc(2,1).eq.amrex_bc_ext_dir .and. jhi.gt.dhi(2)) then
            bc_bx=amrex_box([ilo,dhi(2)+1,klo],[ihi,jhi,khi])
            call solver%user_bc(solver,p,bc_bx,4,time)
         end if

         ! Z-LOW (face=5)
         if (this%lo_bc(3,1).eq.amrex_bc_ext_dir .and. klo.lt.dlo(3)) then
            bc_bx=amrex_box([ilo,jlo,klo],[ihi,jhi,dlo(3)-1])
            call solver%user_bc(solver,p,bc_bx,5,time)
         end if

         ! Z-HIGH (face=6)
         if (this%hi_bc(3,1).eq.amrex_bc_ext_dir .and. khi.gt.dhi(3)) then
            bc_bx=amrex_box([ilo,jlo,dhi(3)+1],[ihi,jhi,khi])
            call solver%user_bc(solver,p,bc_bx,6,time)
         end if

      end do
      call amrex_mfiter_destroy(mfi)

   end subroutine Q_fillbc

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Calculate primitive variables from conserved variables
   !> Q layout: (1) VF*rhoL, (2) (1-VF)*rhoG, (3) VF*rhoL*IL, (4) (1-VF)*rhoG*IG, (5) rhoU, (6) rhoV, (7) rhoW
   subroutine get_primitive(this,Q)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use messager, only: die
      implicit none
      class(amrmpcomp), intent(inout) :: this
      type(amrdata), intent(in) :: Q
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pU,pV,pW
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pRHOL,pRHOG,pIL,pIG,pPL,pPG,pTL,pTG,pC
      real(WP) :: rho_inv,CL,CG
      ! Check passed Q is as expected
      if (Q%ncomp.ne.7) call die('[amrmpcomp get_primitive] Q must have 7 components')
      if (Q%ng.lt.this%nover) call die('[amrmpcomp get_primitive] Q must have at least nover ghost cells')
      ! Check EoS functions are set
      if (.not.associated(this%getPL)) call die('[amrmpcomp get_primitive] getPL not set')
      if (.not.associated(this%getCL)) call die('[amrmpcomp get_primitive] getCL not set')
      if (.not.associated(this%getTL)) call die('[amrmpcomp get_primitive] getTL not set')
      if (.not.associated(this%getPG)) call die('[amrmpcomp get_primitive] getPG not set')
      if (.not.associated(this%getCG)) call die('[amrmpcomp get_primitive] getCG not set')
      if (.not.associated(this%getTG)) call die('[amrmpcomp get_primitive] getTG not set')
      ! Loop over levels
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pQ   =>Q%mf(lvl)%dataptr(mfi)
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            pU   =>this%U%mf(lvl)%dataptr(mfi)
            pV   =>this%V%mf(lvl)%dataptr(mfi)
            pW   =>this%W%mf(lvl)%dataptr(mfi)
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
            pIL  =>this%IL%mf(lvl)%dataptr(mfi)
            pIG  =>this%IG%mf(lvl)%dataptr(mfi)
            pPL  =>this%PL%mf(lvl)%dataptr(mfi)
            pPG  =>this%PG%mf(lvl)%dataptr(mfi)
            pTL  =>this%TL%mf(lvl)%dataptr(mfi)
            pTG  =>this%TG%mf(lvl)%dataptr(mfi)
            pC   =>this%C%mf(lvl)%dataptr(mfi)
            ! Loop over grown tiles
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Compute mixture velocity from momentum
               rho_inv=1.0_WP/max(pQ(i,j,k,1)+pQ(i,j,k,2),this%rho_floor)
               pU(i,j,k,1)=pQ(i,j,k,5)*rho_inv
               pV(i,j,k,1)=pQ(i,j,k,6)*rho_inv
               pW(i,j,k,1)=pQ(i,j,k,7)*rho_inv
               ! Get liquid primitive variables
               if (pVF(i,j,k,1).ge.VFlo) then
                  pRHOL(i,j,k,1)=pQ(i,j,k,1)/pVF(i,j,k,1)
                  pIL  (i,j,k,1)=pQ(i,j,k,3)/pQ(i,j,k,1)
                  pPL  (i,j,k,1)=this%getPL(pRHOL(i,j,k,1),pIL(i,j,k,1))
                  pTL  (i,j,k,1)=this%getTL(pRHOL(i,j,k,1),pPL(i,j,k,1))
                  CL            =this%getCL(pRHOL(i,j,k,1),pPL(i,j,k,1))
               else
                  pRHOL(i,j,k,1)=0.0_WP
                  pIL  (i,j,k,1)=0.0_WP
                  pPL  (i,j,k,1)=0.0_WP
                  pTL  (i,j,k,1)=0.0_WP
                  CL            =0.0_WP
               end if
               ! Get gas primitive variables
               if (pVF(i,j,k,1).le.VFhi) then
                  pRHOG(i,j,k,1)=pQ(i,j,k,2)/(1.0_WP-pVF(i,j,k,1))
                  pIG  (i,j,k,1)=pQ(i,j,k,4)/pQ(i,j,k,2)
                  pPG  (i,j,k,1)=this%getPG(pRHOG(i,j,k,1),pIG(i,j,k,1))
                  pTG  (i,j,k,1)=this%getTG(pRHOG(i,j,k,1),pPG(i,j,k,1))
                  CG            =this%getCG(pRHOG(i,j,k,1),pPG(i,j,k,1))
               else
                  pRHOG(i,j,k,1)=0.0_WP
                  pIG  (i,j,k,1)=0.0_WP
                  pPG  (i,j,k,1)=0.0_WP
                  pTG  (i,j,k,1)=0.0_WP
                  CG            =0.0_WP
               end if
               ! Get mixture speed of sound
               pC(i,j,k,1)=sqrt((pQ(i,j,k,1)*CL**2+pQ(i,j,k,2)*CG**2)*rho_inv)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine get_primitive

   !> Calculate conserved variables from primitive variables
   !> Rebuilds Q from VF, phasic densities, phasic energies, and mixture velocity
   subroutine get_conserved(this)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pU,pV,pW,pRHOL,pRHOG,pIL,pIG
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pQ   =>this%Q%mf(lvl)%dataptr(mfi)
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            pU   =>this%U%mf(lvl)%dataptr(mfi)
            pV   =>this%V%mf(lvl)%dataptr(mfi)
            pW   =>this%W%mf(lvl)%dataptr(mfi)
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
            pIL  =>this%IL%mf(lvl)%dataptr(mfi)
            pIG  =>this%IG%mf(lvl)%dataptr(mfi)
            ! Loop over grown tiles
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pQ(i,j,k,1)=(       pVF(i,j,k,1))*pRHOL(i,j,k,1)
               pQ(i,j,k,2)=(1.0_WP-pVF(i,j,k,1))*pRHOG(i,j,k,1)
               pQ(i,j,k,3)=pQ(i,j,k,1)*pIL(i,j,k,1)
               pQ(i,j,k,4)=pQ(i,j,k,2)*pIG(i,j,k,1)
               pQ(i,j,k,5)=(pQ(i,j,k,1)+pQ(i,j,k,2))*pU(i,j,k,1)
               pQ(i,j,k,6)=(pQ(i,j,k,1)+pQ(i,j,k,2))*pV(i,j,k,1)
               pQ(i,j,k,7)=(pQ(i,j,k,1)+pQ(i,j,k,2))*pW(i,j,k,1)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine get_conserved

   !> Calculate dQdt from passed Q
   subroutine get_dQdt(this,Q,dQdt,dt,time)
      use amrex_amr_module, only: amrex_multifab,amrex_mfiter,amrex_box
      implicit none
      class(amrmpcomp), intent(inout) :: this
      type(amrdata), intent(inout) :: Q
      type(amrdata), intent(inout) :: dQdt
      real(WP), intent(in) :: dt,time
      type(amrex_multifab), dimension(0:this%amr%maxlvl) :: Fx,Fy,Fz
      type(amrex_multifab) :: Vx,Vy,Vz
      type(amrex_multifab) :: band

      ! Shared variables for internal functions
      real(WP) :: dx,dy,dz,dxi,dyi,dzi                              ! Needed for SL transport
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW ! Velocity used for project
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pPLICold ! PLICold used in tet2flux_plic
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQold    ! Qold used in tet2flux_plic
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVFold   ! VFold used in tet2flux_plic

      ! First build primitive variables from Q
      call this%get_primitive(Q)
      
      ! Build transport band at finest level to localize SL computation
      build_band: block
         integer :: dir,n,i,j,k,lvl
         integer, dimension(3) :: ind
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pBand
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         lvl=this%amr%clvl()
         ! Build temporary MultiFab with 1 ghost cell
         call this%amr%mfab_build(lvl=lvl,mfab=band,ncomp=1,nover=1)
         ! Pass 1: Mark interface cells (band=1)
         call band%setval(0.0_WP)
         call this%amr%mfiter_build(lvl=lvl,mfi=mfi)
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
         ! Nullify VFold
         nullify(pVFold)
         ! Synchronize within level
         call band%fill_boundary(this%amr%geom(lvl))
         ! Pass 2: Extend by 1 layer (band=2)
         call this%amr%mfiter_build(lvl=lvl,mfi=mfi)
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

      ! Allocate all fluxes
      define_fluxes: block
         integer :: lvl
         ! Face-centered conserved variable fluxes (7 components)
         do lvl=0,this%amr%clvl()
            call this%amr%mfab_build(lvl=lvl,mfab=Fx(lvl),ncomp=7,nover=0,atface=[.true. ,.false.,.false.]); call Fx(lvl)%setval(0.0_WP)
            call this%amr%mfab_build(lvl=lvl,mfab=Fy(lvl),ncomp=7,nover=0,atface=[.false.,.true. ,.false.]); call Fy(lvl)%setval(0.0_WP)
            call this%amr%mfab_build(lvl=lvl,mfab=Fz(lvl),ncomp=7,nover=0,atface=[.false.,.false.,.true. ]); call Fz(lvl)%setval(0.0_WP)
         end do
         ! Volume moment fluxes at finest level (8 components: Lvol,Gvol,Lbar,Gbar)
         call this%amr%mfab_build(lvl=this%amr%clvl(),mfab=Vx,ncomp=8,nover=0,atface=[.true. ,.false.,.false.]); call Vx%setval(0.0_WP)
         call this%amr%mfab_build(lvl=this%amr%clvl(),mfab=Vy,ncomp=8,nover=0,atface=[.false.,.true. ,.false.]); call Vy%setval(0.0_WP)
         call this%amr%mfab_build(lvl=this%amr%clvl(),mfab=Vz,ncomp=8,nover=0,atface=[.false.,.false.,.true. ]); call Vz%setval(0.0_WP)
      end block define_fluxes

      ! Phase 1a: Semi-Lagrangian fluxes at finest level
      semilagrangian_fluxes: block
         use amrvof_geometry, only: tet_sign,tet_map,correct_flux_poly
         integer :: lvl,i,j,k,n,nn
         real(WP), dimension(3,9) :: face
         real(WP), dimension(3,4) :: tet
         integer , dimension(3,4) :: ijk
         real(WP), dimension(8) :: Vflux
         real(WP), dimension(7) :: Qflux
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pBand,pVx,pVy,pVz,pFx,pFy,pFz
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: fbx
         ! Get finest level info
         lvl=this%amr%clvl()
         dx=this%amr%dx(lvl); dxi=1.0_WP/this%amr%dx(lvl)
         dy=this%amr%dy(lvl); dyi=1.0_WP/this%amr%dy(lvl)
         dz=this%amr%dz(lvl); dzi=1.0_WP/this%amr%dz(lvl)
         ! Loop over finest level tiles
         call this%amr%mfiter_build(lvl=lvl,mfi=mfi)
         do while (mfi%next())
            ! Get data pointers
            pPLICold=>this%PLICold%mf(lvl)%dataptr(mfi)
            pQold   =>this%Qold%mf(lvl)%dataptr(mfi)
            pVFold  =>this%VFold%mf(lvl)%dataptr(mfi)
            pBand   =>band%dataptr(mfi)
            pU      =>this%U%mf(lvl)%dataptr(mfi)
            pV      =>this%V%mf(lvl)%dataptr(mfi)
            pW      =>this%W%mf(lvl)%dataptr(mfi)
            pVx     =>Vx%dataptr(mfi)
            pVy     =>Vy%dataptr(mfi)
            pVz     =>Vz%dataptr(mfi)
            pFx     =>Fx(lvl)%dataptr(mfi)
            pFy     =>Fy(lvl)%dataptr(mfi)
            pFz     =>Fz(lvl)%dataptr(mfi)
            ! X-fluxes
            fbx=mfi%nodaltilebox(1)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               if (maxval(pBand(i-1:i,j,k,1)).eq.0.0_WP) cycle
               ! Build flux polyhedron
               face(:,1)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,5)=project(face(:,1),-dt)
               face(:,2)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,6)=project(face(:,2),-dt)
               face(:,3)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,7)=project(face(:,3),-dt)
               face(:,4)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,8)=project(face(:,4),-dt)
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               call correct_flux_poly(poly=face,target_volume=dt*dy*dz*0.5_WP*(pU(i-1,j,k,1)+pU(i,j,k,1)))
               ! Decompose into tets, cut, and accumulate
               pVx(i,j,k,1:8)=0.0_WP
               pFx(i,j,k,1:7)=0.0_WP
               do n=1,8
                  do nn=1,4
                     tet(:,nn)=face(:,tet_map(nn,n))
                     ijk(:,nn)=floor([(tet(1,nn)-this%amr%xlo)*dxi,(tet(2,nn)-this%amr%ylo)*dyi,(tet(3,nn)-this%amr%zlo)*dzi])
                  end do
                  call tet2flux(tet,ijk,Vflux,Qflux)
                  pVx(i,j,k,1:8)=pVx(i,j,k,1:8)+tet_sign(tet)*Vflux
                  pFx(i,j,k,1:7)=pFx(i,j,k,1:7)+tet_sign(tet)*Qflux
               end do
               ! Convert to flux rate
               pFx(i,j,k,1:7)=pFx(i,j,k,1:7)/(dt*dy*dz)
            end do; end do; end do
            ! Y-fluxes
            fbx=mfi%nodaltilebox(2)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               if (maxval(pBand(i,j-1:j,k,1)).eq.0.0_WP) cycle
               face(:,1)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,5)=project(face(:,1),-dt)
               face(:,2)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,6)=project(face(:,2),-dt)
               face(:,3)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,7)=project(face(:,3),-dt)
               face(:,4)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,8)=project(face(:,4),-dt)
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               call correct_flux_poly(poly=face,target_volume=dt*dz*dx*0.5_WP*(pV(i,j-1,k,1)+pV(i,j,k,1)))
               pVy(i,j,k,1:8)=0.0_WP
               pFy(i,j,k,1:7)=0.0_WP
               do n=1,8
                  do nn=1,4
                     tet(:,nn)=face(:,tet_map(nn,n))
                     ijk(:,nn)=floor([(tet(1,nn)-this%amr%xlo)*dxi,(tet(2,nn)-this%amr%ylo)*dyi,(tet(3,nn)-this%amr%zlo)*dzi])
                  end do
                  call tet2flux(tet,ijk,Vflux,Qflux)
                  pVy(i,j,k,1:8)=pVy(i,j,k,1:8)+tet_sign(tet)*Vflux
                  pFy(i,j,k,1:7)=pFy(i,j,k,1:7)+tet_sign(tet)*Qflux
               end do
               pFy(i,j,k,1:7)=pFy(i,j,k,1:7)/(dt*dz*dx)
            end do; end do; end do
            ! Z-fluxes
            fbx=mfi%nodaltilebox(3)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               if (maxval(pBand(i,j,k-1:k,1)).eq.0.0_WP) cycle
               face(:,1)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,5)=project(face(:,1),-dt)
               face(:,2)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,6)=project(face(:,2),-dt)
               face(:,3)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,7)=project(face(:,3),-dt)
               face(:,4)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,8)=project(face(:,4),-dt)
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               call correct_flux_poly(poly=face,target_volume=dt*dx*dy*0.5_WP*(pW(i,j,k-1,1)+pW(i,j,k,1)))
               pVz(i,j,k,1:8)=0.0_WP
               pFz(i,j,k,1:7)=0.0_WP
               do n=1,8
                  do nn=1,4
                     tet(:,nn)=face(:,tet_map(nn,n))
                     ijk(:,nn)=floor([(tet(1,nn)-this%amr%xlo)*dxi,(tet(2,nn)-this%amr%ylo)*dyi,(tet(3,nn)-this%amr%zlo)*dzi])
                  end do
                  call tet2flux(tet,ijk,Vflux,Qflux)
                  pVz(i,j,k,1:8)=pVz(i,j,k,1:8)+tet_sign(tet)*Vflux
                  pFz(i,j,k,1:7)=pFz(i,j,k,1:7)+tet_sign(tet)*Qflux
               end do
               pFz(i,j,k,1:7)=pFz(i,j,k,1:7)/(dt*dx*dy)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         ! Nullify pointers
         nullify(pU,pV,pW,pPLICold,pQold,pVFold)
      end block semilagrangian_fluxes
      
      ! Phase 1b: Finite volume fluxes for all levels (skip band cells at finest level)
      finitevolume_fluxes: block
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pFx,pFy,pFz,pBand
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pPL,pPG,pTL,pTG,pIL,pIG,pVF
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVisc,pBeta,pDiff
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW ! Intentional masking
         real(WP), dimension(-2: 0) :: wenop
         real(WP), dimension(-1:+1) :: wenom
         real(WP), dimension(1:3,1:3) :: gradU
         real(WP) :: w,div,vel,mflux
         real(WP), parameter :: eps=1.0e-15_WP
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: fbx
         integer :: lvl,i,j,k
         do lvl=0,this%amr%clvl()
            ! Grid spacings for this level
            dx=this%amr%dx(lvl); dxi=1.0_WP/this%amr%dx(lvl)
            dy=this%amr%dy(lvl); dyi=1.0_WP/this%amr%dy(lvl)
            dz=this%amr%dz(lvl); dzi=1.0_WP/this%amr%dz(lvl)
            ! Loop over tiles
            call this%amr%mfiter_build(lvl=lvl,mfi=mfi)
            do while (mfi%next())
               ! Get data pointers
               pQ   =>Q%mf(lvl)%dataptr(mfi)
               pU   =>this%U%mf(lvl)%dataptr(mfi)
               pV   =>this%V%mf(lvl)%dataptr(mfi)
               pW   =>this%W%mf(lvl)%dataptr(mfi)
               pVF  =>this%VF%mf(lvl)%dataptr(mfi)
               pPL  =>this%PL%mf(lvl)%dataptr(mfi)
               pPG  =>this%PG%mf(lvl)%dataptr(mfi)
               pTL  =>this%TL%mf(lvl)%dataptr(mfi)
               pTG  =>this%TG%mf(lvl)%dataptr(mfi)
               pIL  =>this%IL%mf(lvl)%dataptr(mfi)
               pIG  =>this%IG%mf(lvl)%dataptr(mfi)
               pVisc=>this%visc%mf(lvl)%dataptr(mfi)
               pBeta=>this%beta%mf(lvl)%dataptr(mfi)
               pDiff=>this%diff%mf(lvl)%dataptr(mfi)
               pFx  =>Fx(lvl)%dataptr(mfi)
               pFy  =>Fy(lvl)%dataptr(mfi)
               pFz  =>Fz(lvl)%dataptr(mfi)
               if (lvl.eq.this%amr%clvl()) pBand=>band%dataptr(mfi)
               ! X-fluxes
               fbx=mfi%nodaltilebox(1)
               do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  ! Skip band cells at finest level (already filled with SL fluxes)
                  if (lvl.eq.this%amr%clvl()) then; if (maxval(pBand(i-1:i,j,k,1)).gt.0.0_WP) cycle; end if
                  ! Face velocity
                  vel=0.5_WP*sum(pU(i-1:i,j,k,1))
                  ! WENO phasic mass fluxes: Q(1)=rhoL*VF, Q(2)=rhoG*(1-VF)
                  w=weno_weight((abs(pQ(i-1,j,k,1)-pQ(i-2,j,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i-1,j,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pQ(i+1,j,k,1)-pQ(i  ,j,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i-1,j,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFx(i,j,k,1)=-0.5_WP*(vel+abs(vel))*sum(wenop*pQ(i-2:i  ,j,k,1)) &
                  &            -0.5_WP*(vel-abs(vel))*sum(wenom*pQ(i-1:i+1,j,k,1))
                  w=weno_weight((abs(pQ(i-1,j,k,2)-pQ(i-2,j,k,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i-1,j,k,2))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pQ(i+1,j,k,2)-pQ(i  ,j,k,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i-1,j,k,2))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFx(i,j,k,2)=-0.5_WP*(vel+abs(vel))*sum(wenop*pQ(i-2:i  ,j,k,2)) &
                  &            -0.5_WP*(vel-abs(vel))*sum(wenom*pQ(i-1:i+1,j,k,2))
                  ! WENO phasic energy fluxes: F(3)=F(1)*WENO(IL), F(4)=F(2)*WENO(IG)
                  w=weno_weight((abs(pIL(i-1,j,k,1)-pIL(i-2,j,k,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i-1,j,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pIL(i+1,j,k,1)-pIL(i  ,j,k,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i-1,j,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFx(i,j,k,3)=0.5_WP*(pFx(i,j,k,1)-abs(pFx(i,j,k,1)))*sum(wenop*pIL(i-2:i  ,j,k,1)) &
                  &           +0.5_WP*(pFx(i,j,k,1)+abs(pFx(i,j,k,1)))*sum(wenom*pIL(i-1:i+1,j,k,1))
                  w=weno_weight((abs(pIG(i-1,j,k,1)-pIG(i-2,j,k,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i-1,j,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pIG(i+1,j,k,1)-pIG(i  ,j,k,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i-1,j,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFx(i,j,k,4)=0.5_WP*(pFx(i,j,k,2)-abs(pFx(i,j,k,2)))*sum(wenop*pIG(i-2:i  ,j,k,1)) &
                  &           +0.5_WP*(pFx(i,j,k,2)+abs(pFx(i,j,k,2)))*sum(wenom*pIG(i-1:i+1,j,k,1))
                  ! Total mass flux for momentum
                  mflux=pFx(i,j,k,1)+pFx(i,j,k,2)
                  ! Momentum fluxes with pressure stress: P_face = VF*PL + (1-VF)*PG
                  pFx(i,j,k,5)=mflux*0.5_WP*sum(pU(i-1:i,j,k,1))-0.5_WP*sum(pVF(i-1:i,j,k,1)*pPL(i-1:i,j,k,1)+(1.0_WP-pVF(i-1:i,j,k,1))*pPG(i-1:i,j,k,1))
                  pFx(i,j,k,6)=mflux*0.5_WP*sum(pV(i-1:i,j,k,1))
                  pFx(i,j,k,7)=mflux*0.5_WP*sum(pW(i-1:i,j,k,1))
                  ! Velocity gradients at x-face
                  gradU(1,1)=dxi*(pU(i,j,k,1)-pU(i-1,j,k,1))
                  gradU(2,1)=0.25_WP*dyi*(pU(i-1,j+1,k,1)-pU(i-1,j-1,k,1)+pU(i,j+1,k,1)-pU(i,j-1,k,1))
                  gradU(3,1)=0.25_WP*dzi*(pU(i-1,j,k+1,1)-pU(i-1,j,k-1,1)+pU(i,j,k+1,1)-pU(i,j,k-1,1))
                  gradU(1,2)=dxi*(pV(i,j,k,1)-pV(i-1,j,k,1))
                  gradU(2,2)=0.25_WP*dyi*(pV(i-1,j+1,k,1)-pV(i-1,j-1,k,1)+pV(i,j+1,k,1)-pV(i,j-1,k,1))
                  gradU(3,2)=0.25_WP*dzi*(pV(i-1,j,k+1,1)-pV(i-1,j,k-1,1)+pV(i,j,k+1,1)-pV(i,j,k-1,1))
                  gradU(1,3)=dxi*(pW(i,j,k,1)-pW(i-1,j,k,1))
                  gradU(2,3)=0.25_WP*dyi*(pW(i-1,j+1,k,1)-pW(i-1,j-1,k,1)+pW(i,j+1,k,1)-pW(i,j-1,k,1))
                  gradU(3,3)=0.25_WP*dzi*(pW(i-1,j,k+1,1)-pW(i-1,j,k-1,1)+pW(i,j,k+1,1)-pW(i,j,k-1,1))
                  div=gradU(1,1)+gradU(2,2)+gradU(3,3)
                  ! Viscous stress at x-face (added to momentum fluxes)
                  pFx(i,j,k,5)=pFx(i,j,k,5)+0.5_WP*sum(pVisc(i-1:i,j,k,1))*(gradU(1,1)+gradU(1,1))+0.5_WP*(sum(pBeta(i-1:i,j,k,1))-2.0_WP/3.0_WP*sum(pVisc(i-1:i,j,k,1)))*div
                  pFx(i,j,k,6)=pFx(i,j,k,6)+0.5_WP*sum(pVisc(i-1:i,j,k,1))*(gradU(2,1)+gradU(1,2))
                  pFx(i,j,k,7)=pFx(i,j,k,7)+0.5_WP*sum(pVisc(i-1:i,j,k,1))*(gradU(3,1)+gradU(1,3))
                  ! VF-weighted heat diffusion flux
                  pFx(i,j,k,3)=pFx(i,j,k,3)+0.5_WP*sum(       pVF(i-1:i,j,k,1))*0.5_WP*sum(pDiff(i-1:i,j,k,1))*dxi*(pTL(i,j,k,1)-pTL(i-1,j,k,1))
                  pFx(i,j,k,4)=pFx(i,j,k,4)+0.5_WP*sum(1.0_WP-pVF(i-1:i,j,k,1))*0.5_WP*sum(pDiff(i-1:i,j,k,1))*dxi*(pTG(i,j,k,1)-pTG(i-1,j,k,1))
               end do; end do; end do
               ! Y-fluxes
               fbx=mfi%nodaltilebox(2)
               do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  ! Skip band cells at finest level (already filled with SL fluxes)
                  if (lvl.eq.this%amr%clvl()) then; if (maxval(pBand(i,j-1:j,k,1)).gt.0.0_WP) cycle; end if
                  ! Face velocity
                  vel=0.5_WP*sum(pV(i,j-1:j,k,1))
                  ! WENO phasic mass fluxes: Q(1)=rhoL*VF, Q(2)=rhoG*(1-VF)
                  w=weno_weight((abs(pQ(i,j-1,k,1)-pQ(i,j-2,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j-1,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pQ(i,j+1,k,1)-pQ(i,j  ,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j-1,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFy(i,j,k,1)=-0.5_WP*(vel+abs(vel))*sum(wenop*pQ(i,j-2:j  ,k,1)) &
                  &            -0.5_WP*(vel-abs(vel))*sum(wenom*pQ(i,j-1:j+1,k,1))
                  w=weno_weight((abs(pQ(i,j-1,k,2)-pQ(i,j-2,k,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i,j-1,k,2))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pQ(i,j+1,k,2)-pQ(i,j  ,k,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i,j-1,k,2))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFy(i,j,k,2)=-0.5_WP*(vel+abs(vel))*sum(wenop*pQ(i,j-2:j  ,k,2)) &
                  &            -0.5_WP*(vel-abs(vel))*sum(wenom*pQ(i,j-1:j+1,k,2))
                  ! WENO phasic energy fluxes: F(3)=F(1)*WENO(IL), F(4)=F(2)*WENO(IG)
                  w=weno_weight((abs(pIL(i,j-1,k,1)-pIL(i,j-2,k,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i,j-1,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pIL(i,j+1,k,1)-pIL(i,j  ,k,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i,j-1,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFy(i,j,k,3)=0.5_WP*(pFy(i,j,k,1)-abs(pFy(i,j,k,1)))*sum(wenop*pIL(i,j-2:j  ,k,1)) &
                  &           +0.5_WP*(pFy(i,j,k,1)+abs(pFy(i,j,k,1)))*sum(wenom*pIL(i,j-1:j+1,k,1))
                  w=weno_weight((abs(pIG(i,j-1,k,1)-pIG(i,j-2,k,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i,j-1,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pIG(i,j+1,k,1)-pIG(i,j  ,k,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i,j-1,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFy(i,j,k,4)=0.5_WP*(pFy(i,j,k,2)-abs(pFy(i,j,k,2)))*sum(wenop*pIG(i,j-2:j  ,k,1)) &
                  &           +0.5_WP*(pFy(i,j,k,2)+abs(pFy(i,j,k,2)))*sum(wenom*pIG(i,j-1:j+1,k,1))
                  mflux=pFy(i,j,k,1)+pFy(i,j,k,2)
                  ! Momentum fluxes with pressure stress
                  pFy(i,j,k,5)=mflux*0.5_WP*sum(pU(i,j-1:j,k,1))
                  pFy(i,j,k,6)=mflux*0.5_WP*sum(pV(i,j-1:j,k,1))-0.5_WP*sum(pVF(i,j-1:j,k,1)*pPL(i,j-1:j,k,1)+(1.0_WP-pVF(i,j-1:j,k,1))*pPG(i,j-1:j,k,1))
                  pFy(i,j,k,7)=mflux*0.5_WP*sum(pW(i,j-1:j,k,1))
                  ! Velocity gradients at y-face
                  gradU(1,1)=0.25_WP*dxi*(pU(i+1,j-1,k,1)-pU(i-1,j-1,k,1)+pU(i+1,j,k,1)-pU(i-1,j,k,1))
                  gradU(2,1)=dyi*(pU(i,j,k,1)-pU(i,j-1,k,1))
                  gradU(3,1)=0.25_WP*dzi*(pU(i,j-1,k+1,1)-pU(i,j-1,k-1,1)+pU(i,j,k+1,1)-pU(i,j,k-1,1))
                  gradU(1,2)=0.25_WP*dxi*(pV(i+1,j-1,k,1)-pV(i-1,j-1,k,1)+pV(i+1,j,k,1)-pV(i-1,j,k,1))
                  gradU(2,2)=dyi*(pV(i,j,k,1)-pV(i,j-1,k,1))
                  gradU(3,2)=0.25_WP*dzi*(pV(i,j-1,k+1,1)-pV(i,j-1,k-1,1)+pV(i,j,k+1,1)-pV(i,j,k-1,1))
                  gradU(1,3)=0.25_WP*dxi*(pW(i+1,j-1,k,1)-pW(i-1,j-1,k,1)+pW(i+1,j,k,1)-pW(i-1,j,k,1))
                  gradU(2,3)=dyi*(pW(i,j,k,1)-pW(i,j-1,k,1))
                  gradU(3,3)=0.25_WP*dzi*(pW(i,j-1,k+1,1)-pW(i,j-1,k-1,1)+pW(i,j,k+1,1)-pW(i,j,k-1,1))
                  div=gradU(1,1)+gradU(2,2)+gradU(3,3)
                  ! Viscous stress at y-face
                  pFy(i,j,k,5)=pFy(i,j,k,5)+0.5_WP*sum(pVisc(i,j-1:j,k,1))*(gradU(1,2)+gradU(2,1))
                  pFy(i,j,k,6)=pFy(i,j,k,6)+0.5_WP*sum(pVisc(i,j-1:j,k,1))*(gradU(2,2)+gradU(2,2))+0.5_WP*(sum(pBeta(i,j-1:j,k,1))-2.0_WP/3.0_WP*sum(pVisc(i,j-1:j,k,1)))*div
                  pFy(i,j,k,7)=pFy(i,j,k,7)+0.5_WP*sum(pVisc(i,j-1:j,k,1))*(gradU(3,2)+gradU(2,3))
                  ! VF-weighted heat diffusion flux
                  pFy(i,j,k,3)=pFy(i,j,k,3)+0.5_WP*sum(       pVF(i,j-1:j,k,1))*0.5_WP*sum(pDiff(i,j-1:j,k,1))*dyi*(pTL(i,j,k,1)-pTL(i,j-1,k,1))
                  pFy(i,j,k,4)=pFy(i,j,k,4)+0.5_WP*sum(1.0_WP-pVF(i,j-1:j,k,1))*0.5_WP*sum(pDiff(i,j-1:j,k,1))*dyi*(pTG(i,j,k,1)-pTG(i,j-1,k,1))
               end do; end do; end do
               ! Z-fluxes
               fbx=mfi%nodaltilebox(3)
               do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  ! Skip band cells at finest level (already filled with SL fluxes)
                  if (lvl.eq.this%amr%clvl()) then; if (maxval(pBand(i,j,k-1:k,1)).gt.0.0_WP) cycle; end if
                  ! Face velocity
                  vel=0.5_WP*sum(pW(i,j,k-1:k,1))
                  ! WENO phasic mass fluxes: Q(1)=rhoL*VF, Q(2)=rhoG*(1-VF)
                  w=weno_weight((abs(pQ(i,j,k-1,1)-pQ(i,j,k-2,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j,k-1,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pQ(i,j,k+1,1)-pQ(i,j,k  ,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j,k-1,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFz(i,j,k,1)=-0.5_WP*(vel+abs(vel))*sum(wenop*pQ(i,j,k-2:k  ,1)) &
                  &            -0.5_WP*(vel-abs(vel))*sum(wenom*pQ(i,j,k-1:k+1,1))
                  w=weno_weight((abs(pQ(i,j,k-1,2)-pQ(i,j,k-2,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i,j,k-1,2))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pQ(i,j,k+1,2)-pQ(i,j,k  ,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i,j,k-1,2))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFz(i,j,k,2)=-0.5_WP*(vel+abs(vel))*sum(wenop*pQ(i,j,k-2:k  ,2)) &
                  &            -0.5_WP*(vel-abs(vel))*sum(wenom*pQ(i,j,k-1:k+1,2))
                  ! WENO phasic energy fluxes: F(3)=F(1)*WENO(IL), F(4)=F(2)*WENO(IG)
                  w=weno_weight((abs(pIL(i,j,k-1,1)-pIL(i,j,k-2,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i,j,k-1,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pIL(i,j,k+1,1)-pIL(i,j,k  ,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i,j,k-1,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFz(i,j,k,3)=0.5_WP*(pFz(i,j,k,1)-abs(pFz(i,j,k,1)))*sum(wenop*pIL(i,j,k-2:k  ,1)) &
                  &           +0.5_WP*(pFz(i,j,k,1)+abs(pFz(i,j,k,1)))*sum(wenom*pIL(i,j,k-1:k+1,1))
                  w=weno_weight((abs(pIG(i,j,k-1,1)-pIG(i,j,k-2,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i,j,k-1,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pIG(i,j,k+1,1)-pIG(i,j,k  ,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i,j,k-1,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFz(i,j,k,4)=0.5_WP*(pFz(i,j,k,2)-abs(pFz(i,j,k,2)))*sum(wenop*pIG(i,j,k-2:k  ,1)) &
                  &           +0.5_WP*(pFz(i,j,k,2)+abs(pFz(i,j,k,2)))*sum(wenom*pIG(i,j,k-1:k+1,1))
                  mflux=pFz(i,j,k,1)+pFz(i,j,k,2)
                  ! Momentum fluxes with pressure stress
                  pFz(i,j,k,5)=mflux*0.5_WP*sum(pU(i,j,k-1:k,1))
                  pFz(i,j,k,6)=mflux*0.5_WP*sum(pV(i,j,k-1:k,1))
                  pFz(i,j,k,7)=mflux*0.5_WP*sum(pW(i,j,k-1:k,1))-0.5_WP*sum(pVF(i,j,k-1:k,1)*pPL(i,j,k-1:k,1)+(1.0_WP-pVF(i,j,k-1:k,1))*pPG(i,j,k-1:k,1))
                  ! Velocity gradients at z-face
                  gradU(1,1)=0.25_WP*dxi*(pU(i+1,j,k-1,1)-pU(i-1,j,k-1,1)+pU(i+1,j,k,1)-pU(i-1,j,k,1))
                  gradU(2,1)=0.25_WP*dyi*(pU(i,j+1,k-1,1)-pU(i,j-1,k-1,1)+pU(i,j+1,k,1)-pU(i,j-1,k,1))
                  gradU(3,1)=dzi*(pU(i,j,k,1)-pU(i,j,k-1,1))
                  gradU(1,2)=0.25_WP*dxi*(pV(i+1,j,k-1,1)-pV(i-1,j,k-1,1)+pV(i+1,j,k,1)-pV(i-1,j,k,1))
                  gradU(2,2)=0.25_WP*dyi*(pV(i,j+1,k-1,1)-pV(i,j-1,k-1,1)+pV(i,j+1,k,1)-pV(i,j-1,k,1))
                  gradU(3,2)=dzi*(pV(i,j,k,1)-pV(i,j,k-1,1))
                  gradU(1,3)=0.25_WP*dxi*(pW(i+1,j,k-1,1)-pW(i-1,j,k-1,1)+pW(i+1,j,k,1)-pW(i-1,j,k,1))
                  gradU(2,3)=0.25_WP*dyi*(pW(i,j+1,k-1,1)-pW(i,j-1,k-1,1)+pW(i,j+1,k,1)-pW(i,j-1,k,1))
                  gradU(3,3)=dzi*(pW(i,j,k,1)-pW(i,j,k-1,1))
                  div=gradU(1,1)+gradU(2,2)+gradU(3,3)
                  ! Viscous stress at z-face
                  pFz(i,j,k,5)=pFz(i,j,k,5)+0.5_WP*sum(pVisc(i,j,k-1:k,1))*(gradU(1,3)+gradU(3,1))
                  pFz(i,j,k,6)=pFz(i,j,k,6)+0.5_WP*sum(pVisc(i,j,k-1:k,1))*(gradU(2,3)+gradU(3,2))
                  pFz(i,j,k,7)=pFz(i,j,k,7)+0.5_WP*sum(pVisc(i,j,k-1:k,1))*(gradU(3,3)+gradU(3,3))+0.5_WP*(sum(pBeta(i,j,k-1:k,1))-2.0_WP/3.0_WP*sum(pVisc(i,j,k-1:k,1)))*div
                  ! VF-weighted heat diffusion flux
                  pFz(i,j,k,3)=pFz(i,j,k,3)+0.5_WP*sum(       pVF(i,j,k-1:k,1))*0.5_WP*sum(pDiff(i,j,k-1:k,1))*dzi*(pTL(i,j,k,1)-pTL(i,j,k-1,1))
                  pFz(i,j,k,4)=pFz(i,j,k,4)+0.5_WP*sum(1.0_WP-pVF(i,j,k-1:k,1))*0.5_WP*sum(pDiff(i,j,k-1:k,1))*dzi*(pTG(i,j,k,1)-pTG(i,j,k-1,1))
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
      end block finitevolume_fluxes

      ! Phase 2: Average down all fluxes for C/F conservation
      c_f_consistency: block
         use amrex_interface, only: amrmfab_average_down_face
         integer :: lvl
         do lvl=this%amr%clvl(),1,-1
            call amrmfab_average_down_face(fmf=Fx(lvl),cmf=Fx(lvl-1),rr=this%amr%rref(lvl-1),cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=Fy(lvl),cmf=Fy(lvl-1),rr=this%amr%rref(lvl-1),cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=Fz(lvl),cmf=Fz(lvl-1),rr=this%amr%rref(lvl-1),cgeom=this%amr%geom(lvl-1))
         end do
      end block c_f_consistency
      
      ! Phase 3: Compute divergence and source terms for all levels, update VF/bary at band
      divergence_and_sources: block
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         integer :: lvl,i,j,k
         real(WP), dimension(:,:,:,:), contiguous, pointer :: rhs,pFx,pFy,pFz,pBand
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVFold  ! Intentional masking
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pPL,pPG,pVisc,pBeta
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVx,pVy,pVz
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pCliq,pCgas,pCliqold,pCgasold
         real(WP), dimension(1:3,1:3) :: gradU
         real(WP) :: div,vol
         real(WP) :: Lvol_old,Lvol_new,Lvol_flux
         real(WP) :: Gvol_old,Gvol_new,Gvol_flux
         real(WP), dimension(3) :: Lbar_old,Lbar_new,Lbar_flux
         real(WP), dimension(3) :: Gbar_old,Gbar_new,Gbar_flux
         do lvl=0,this%amr%clvl()
            ! Grid spacings for this level
            dx=this%amr%dx(lvl); dxi=1.0_WP/dx
            dy=this%amr%dy(lvl); dyi=1.0_WP/dy
            dz=this%amr%dz(lvl); dzi=1.0_WP/dz
            vol=dx*dy*dz
            ! Loop over tiles
            call this%amr%mfiter_build(lvl=lvl,mfi=mfi)
            do while (mfi%next())
               ! Get data pointers
               rhs  =>dQdt%mf(lvl)%dataptr(mfi)
               pFx  =>Fx(lvl)%dataptr(mfi)
               pFy  =>Fy(lvl)%dataptr(mfi)
               pFz  =>Fz(lvl)%dataptr(mfi)
               pU   =>this%U%mf(lvl)%dataptr(mfi)
               pV   =>this%V%mf(lvl)%dataptr(mfi)
               pW   =>this%W%mf(lvl)%dataptr(mfi)
               pVF  =>this%VF%mf(lvl)%dataptr(mfi)
               pPL  =>this%PL%mf(lvl)%dataptr(mfi)
               pPG  =>this%PG%mf(lvl)%dataptr(mfi)
               pVisc=>this%visc%mf(lvl)%dataptr(mfi)
               pBeta=>this%beta%mf(lvl)%dataptr(mfi)
               ! Extra pointers at finest level
               if (lvl.eq.this%amr%clvl()) then
                  pBand   =>band%dataptr(mfi)
                  pVx     =>Vx%dataptr(mfi)
                  pVy     =>Vy%dataptr(mfi)
                  pVz     =>Vz%dataptr(mfi)
                  pVFold  =>this%VFold%mf(lvl)%dataptr(mfi)
                  pCliq   =>this%Cliq%mf(lvl)%dataptr(mfi)
                  pCgas   =>this%Cgas%mf(lvl)%dataptr(mfi)
                  pCliqold=>this%Cliqold%mf(lvl)%dataptr(mfi)
                  pCgasold=>this%Cgasold%mf(lvl)%dataptr(mfi)
               end if
               ! Loop over interior
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Divergence of conserved variable fluxes (7 components)
                  rhs(i,j,k,:)=dxi*(pFx(i+1,j,k,:)-pFx(i,j,k,:))+dyi*(pFy(i,j+1,k,:)-pFy(i,j,k,:))+dzi*(pFz(i,j,k+1,:)-pFz(i,j,k,:))
                  ! Velocity gradients at cell center
                  gradU(1,1)=0.5_WP*dxi*(pU(i+1,j,k,1)-pU(i-1,j,k,1))
                  gradU(2,1)=0.5_WP*dyi*(pU(i,j+1,k,1)-pU(i,j-1,k,1))
                  gradU(3,1)=0.5_WP*dzi*(pU(i,j,k+1,1)-pU(i,j,k-1,1))
                  gradU(1,2)=0.5_WP*dxi*(pV(i+1,j,k,1)-pV(i-1,j,k,1))
                  gradU(2,2)=0.5_WP*dyi*(pV(i,j+1,k,1)-pV(i,j-1,k,1))
                  gradU(3,2)=0.5_WP*dzi*(pV(i,j,k+1,1)-pV(i,j,k-1,1))
                  gradU(1,3)=0.5_WP*dxi*(pW(i+1,j,k,1)-pW(i-1,j,k,1))
                  gradU(2,3)=0.5_WP*dyi*(pW(i,j+1,k,1)-pW(i,j-1,k,1))
                  gradU(3,3)=0.5_WP*dzi*(pW(i,j,k+1,1)-pW(i,j,k-1,1))
                  div=gradU(1,1)+gradU(2,2)+gradU(3,3)
                  ! Pressure dilatation: split by VF between phasic energies
                  rhs(i,j,k,3)=rhs(i,j,k,3)-(       pVF(i,j,k,1))*pPL(i,j,k,1)*div
                  rhs(i,j,k,4)=rhs(i,j,k,4)-(1.0_WP-pVF(i,j,k,1))*pPG(i,j,k,1)*div
                  ! Viscous heating: τ:∇U, split by VF between phasic energies
                  rhs(i,j,k,3)=rhs(i,j,k,3)+(       pVF(i,j,k,1))*( &
                  & (2.0_WP*pVisc(i,j,k,1)*gradU(1,1)+(pBeta(i,j,k,1)-2.0_WP/3.0_WP*pVisc(i,j,k,1))*div)*gradU(1,1) &
                  &+(2.0_WP*pVisc(i,j,k,1)*gradU(2,2)+(pBeta(i,j,k,1)-2.0_WP/3.0_WP*pVisc(i,j,k,1))*div)*gradU(2,2) &
                  &+(2.0_WP*pVisc(i,j,k,1)*gradU(3,3)+(pBeta(i,j,k,1)-2.0_WP/3.0_WP*pVisc(i,j,k,1))*div)*gradU(3,3) &
                  &+pVisc(i,j,k,1)*(gradU(2,1)+gradU(1,2))*(gradU(2,1)+gradU(1,2)) &
                  &+pVisc(i,j,k,1)*(gradU(3,1)+gradU(1,3))*(gradU(3,1)+gradU(1,3)) &
                  &+pVisc(i,j,k,1)*(gradU(3,2)+gradU(2,3))*(gradU(3,2)+gradU(2,3)))
                  rhs(i,j,k,4)=rhs(i,j,k,4)+(1.0_WP-pVF(i,j,k,1))*( &
                  & (2.0_WP*pVisc(i,j,k,1)*gradU(1,1)+(pBeta(i,j,k,1)-2.0_WP/3.0_WP*pVisc(i,j,k,1))*div)*gradU(1,1) &
                  &+(2.0_WP*pVisc(i,j,k,1)*gradU(2,2)+(pBeta(i,j,k,1)-2.0_WP/3.0_WP*pVisc(i,j,k,1))*div)*gradU(2,2) &
                  &+(2.0_WP*pVisc(i,j,k,1)*gradU(3,3)+(pBeta(i,j,k,1)-2.0_WP/3.0_WP*pVisc(i,j,k,1))*div)*gradU(3,3) &
                  &+pVisc(i,j,k,1)*(gradU(2,1)+gradU(1,2))*(gradU(2,1)+gradU(1,2)) &
                  &+pVisc(i,j,k,1)*(gradU(3,1)+gradU(1,3))*(gradU(3,1)+gradU(1,3)) &
                  &+pVisc(i,j,k,1)*(gradU(3,2)+gradU(2,3))*(gradU(3,2)+gradU(2,3)))
                  ! VF/barycenter update at band cells (finest level only)
                  if (lvl.eq.this%amr%clvl()) then
                     ! Skip if cell not in band
                     if (pBand(i,j,k,1).eq.0.0_WP) cycle
                     ! Old phasic moments
                     Lvol_old=(       pVFold(i,j,k,1))*vol
                     Gvol_old=(1.0_WP-pVFold(i,j,k,1))*vol
                     Lbar_old=pCliqold(i,j,k,1:3)
                     Gbar_old=pCgasold(i,j,k,1:3)
                     ! Net volume flux (outflow positive) from SL volume moments
                     Lvol_flux=pVx(i+1,j,k, 1 )-pVx(i,j,k, 1 )+pVy(i,j+1,k, 1 )-pVy(i,j,k, 1 )+pVz(i,j,k+1, 1 )-pVz(i,j,k, 1 )
                     Gvol_flux=pVx(i+1,j,k, 2 )-pVx(i,j,k, 2 )+pVy(i,j+1,k, 2 )-pVy(i,j,k, 2 )+pVz(i,j,k+1, 2 )-pVz(i,j,k, 2 )
                     Lbar_flux=pVx(i+1,j,k,3:5)-pVx(i,j,k,3:5)+pVy(i,j+1,k,3:5)-pVy(i,j,k,3:5)+pVz(i,j,k+1,3:5)-pVz(i,j,k,3:5)
                     Gbar_flux=pVx(i+1,j,k,6:8)-pVx(i,j,k,6:8)+pVy(i,j+1,k,6:8)-pVy(i,j,k,6:8)+pVz(i,j,k+1,6:8)-pVz(i,j,k,6:8)
                     ! New phasic volumes
                     Lvol_new=Lvol_old-Lvol_flux
                     Gvol_new=Gvol_old-Gvol_flux
                     ! New VF and default barycenters
                     pVF(i,j,k,1)=Lvol_new/(Lvol_new+Gvol_new)
                     pCliq(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                     pCgas(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
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
                  end if
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
      end block divergence_and_sources

      ! Cleanup temporary mfabs
      cleanup: block
         integer :: lvl
         call this%amr%mfab_destroy(band)
         call this%amr%mfab_destroy(Vx)
         call this%amr%mfab_destroy(Vy)
         call this%amr%mfab_destroy(Vz)
         do lvl=0,this%amr%clvl()
            call this%amr%mfab_destroy(Fx(lvl))
            call this%amr%mfab_destroy(Fy(lvl))
            call this%amr%mfab_destroy(Fz(lvl))
         end do
      end block cleanup

      ! Sync and apply BC
      call this%fill_moments_lvl(this%amr%clvl(),time)
      
   contains

      !> WENO switch function
      real(WP) function weno_weight(ratio)
         implicit none
         real(WP), intent(in) :: ratio
         real(WP), parameter :: lambda=0.13_WP
         real(WP), parameter :: delta=0.01_WP
         weno_weight=(1.0_WP-tanh((ratio-lambda)/delta))/3.0_WP+(1.0_WP-tanh((ratio-1.0_WP/lambda)/delta))/6.0_WP
      end function weno_weight

      !> Recursive subroutine that cuts a tet by grid planes to compute volume and Q fluxes
      recursive subroutine tet2flux(mytet,myind,myVflux,myQflux)
         use amrvof_geometry, only: cut_side,cut_v1,cut_v2,cut_vtet,cut_ntets,cut_nvert
         real(WP), dimension(3,4), intent(in) :: mytet
         integer,  dimension(3,4), intent(in) :: myind
         real(WP), dimension(8),  intent(out) :: myVflux
         real(WP), dimension(7),  intent(out) :: myQflux
         integer :: dir,cut_ind,icase,n1,n2,v1,v2
         real(WP), dimension(4) :: dd
         real(WP), dimension(3,8) :: vert
         integer,  dimension(3,8,2) :: vert_ind
         real(WP) :: mu,my_vol
         real(WP), dimension(3,4) :: newtet
         integer,  dimension(3,4) :: newind
         real(WP), dimension(3) :: a,b,c
         real(WP), dimension(8) :: subVflux
         real(WP), dimension(7) :: subQflux
         real(WP) :: xcut,ycut,zcut
         
         myVflux=0.0_WP
         myQflux=0.0_WP
         
         ! Determine if tet spans multiple cells and needs cutting
         if (maxval(myind(1,:))-minval(myind(1,:)).gt.0) then
            dir=1; cut_ind=maxval(myind(1,:))
            xcut=this%amr%xlo+real(cut_ind,WP)*dx
            dd(:)=mytet(1,:)-xcut
         else if (maxval(myind(2,:))-minval(myind(2,:)).gt.0) then
            dir=2; cut_ind=maxval(myind(2,:))
            ycut=this%amr%ylo+real(cut_ind,WP)*dy
            dd(:)=mytet(2,:)-ycut
         else if (maxval(myind(3,:))-minval(myind(3,:)).gt.0) then
            dir=3; cut_ind=maxval(myind(3,:))
            zcut=this%amr%zlo+real(cut_ind,WP)*dz
            dd(:)=mytet(3,:)-zcut
         else
            ! All vertices in same cell - cut by PLIC and return
            call tet2flux_plic(mytet,myind(1,1),myind(2,1),myind(3,1),myVflux,myQflux)
            return
         end if
         
         ! Find cut case (1-indexed: 1-16)
         icase=1+int(0.5_WP+sign(0.5_WP,dd(1))) &
              +2*int(0.5_WP+sign(0.5_WP,dd(2))) &
              +4*int(0.5_WP+sign(0.5_WP,dd(3))) &
              +8*int(0.5_WP+sign(0.5_WP,dd(4)))
         
         ! Copy vertices and indices
         do n1=1,4
            vert(:,n1)=mytet(:,n1)
            vert_ind(:,n1,1)=myind(:,n1)
            vert_ind(:,n1,2)=myind(:,n1)
            vert_ind(dir,n1,1)=min(vert_ind(dir,n1,1),cut_ind-1)
            vert_ind(dir,n1,2)=max(vert_ind(dir,n1,1),cut_ind)
         end do
         
         ! Create interpolated vertices on cut plane
         do n1=1,cut_nvert(icase)
            v1=cut_v1(n1,icase); v2=cut_v2(n1,icase)
            mu=min(1.0_WP,max(0.0_WP,-dd(v1)/(sign(abs(dd(v2)-dd(v1))+epsilon(1.0_WP),dd(v2)-dd(v1)))))
            vert(:,4+n1)=(1.0_WP-mu)*vert(:,v1)+mu*vert(:,v2)
            vert_ind(1,4+n1,1)=floor((vert(1,4+n1)-this%amr%xlo)*dxi)
            vert_ind(2,4+n1,1)=floor((vert(2,4+n1)-this%amr%ylo)*dyi)
            vert_ind(3,4+n1,1)=floor((vert(3,4+n1)-this%amr%zlo)*dzi)
            vert_ind(:,4+n1,1)=max(vert_ind(:,4+n1,1),min(vert_ind(:,v1,1),vert_ind(:,v2,1)))
            vert_ind(:,4+n1,1)=min(vert_ind(:,4+n1,1),max(vert_ind(:,v1,1),vert_ind(:,v2,1)))
            vert_ind(:,4+n1,2)=vert_ind(:,4+n1,1)
            vert_ind(dir,4+n1,1)=cut_ind-1
            vert_ind(dir,4+n1,2)=cut_ind
         end do
         
         ! Create and process sub-tets
         do n1=1,cut_ntets(icase)
            do n2=1,4
               newtet(:,n2)=vert(:,cut_vtet(n2,n1,icase))
               newind(:,n2)=vert_ind(:,cut_vtet(n2,n1,icase),cut_side(n1,icase))
            end do
            a=newtet(:,1)-newtet(:,4)
            b=newtet(:,2)-newtet(:,4)
            c=newtet(:,3)-newtet(:,4)
            my_vol=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
            if (my_vol.lt.1.0e-15_WP*dx*dy*dz) cycle
            call tet2flux(newtet,newind,subVflux,subQflux)
            myVflux=myVflux+subVflux
            myQflux=myQflux+subQflux
         end do
         
      end subroutine tet2flux

      !> Cut tet by PLIC and compute volume + conserved variable fluxes
      subroutine tet2flux_plic(mytet,i0,j0,k0,myVflux,myQflux)
         use amrvof_geometry, only: cut_v1,cut_v2,cut_vtet,cut_ntets,cut_nvert,cut_nntet,tet_vol
         use messager, only: die
         real(WP), dimension(3,4), intent(in) :: mytet
         integer,  intent(in) :: i0,j0,k0
         real(WP), dimension(8),  intent(out) :: myVflux
         real(WP), dimension(7),  intent(out) :: myQflux
         integer :: icase,n1,v1,v2
         real(WP), dimension(4) :: dd
         real(WP), dimension(3,8) :: vert
         real(WP), dimension(3) :: a,b,c,bary,normal
         real(WP) :: mu,my_vol,dist,VF0
         
         myVflux=0.0_WP
         myQflux=0.0_WP

         ! Check indices are within PLICold bounds
         if (i0.lt.lbound(pPLICold,1).or.i0.gt.ubound(pPLICold,1).or. &
             j0.lt.lbound(pPLICold,2).or.j0.gt.ubound(pPLICold,2).or. &
             k0.lt.lbound(pPLICold,3).or.k0.gt.ubound(pPLICold,3)) then
            call die('[tet2flux_plic] Index out of bounds - check CFL or ghost cells')
         end if
         
         ! Get old VF for this cell
         VF0=pVFold(i0,j0,k0,1)
         
         ! Pure cell shortcut
         if (pPLICold(i0,j0,k0,4).gt.+1.0e9_WP) then
            ! Pure liquid
            my_vol=abs(tet_vol(mytet))
            bary=0.25_WP*(mytet(:,1)+mytet(:,2)+mytet(:,3)+mytet(:,4))
            myVflux( 1 )=my_vol
            myVflux(3:5)=my_vol*bary
            ! Q flux: all mass is liquid
            myQflux=my_vol*pQold(i0,j0,k0,:)
            return
         else if (pPLICold(i0,j0,k0,4).lt.-1.0e9_WP) then
            ! Pure gas
            my_vol=abs(tet_vol(mytet))
            bary=0.25_WP*(mytet(:,1)+mytet(:,2)+mytet(:,3)+mytet(:,4))
            myVflux( 2 )=my_vol
            myVflux(6:8)=my_vol*bary
            ! Q flux: all mass is gas
            myQflux=my_vol*pQold(i0,j0,k0,:)
            return
         end if
         
         ! Get PLIC from this cell
         normal=pPLICold(i0,j0,k0,1:3)
         dist=pPLICold(i0,j0,k0,4)
         
         ! Compute signed distance to plane for each vertex
         dd(1)=normal(1)*mytet(1,1)+normal(2)*mytet(2,1)+normal(3)*mytet(3,1)-dist
         dd(2)=normal(1)*mytet(1,2)+normal(2)*mytet(2,2)+normal(3)*mytet(3,2)-dist
         dd(3)=normal(1)*mytet(1,3)+normal(2)*mytet(2,3)+normal(3)*mytet(3,3)-dist
         dd(4)=normal(1)*mytet(1,4)+normal(2)*mytet(2,4)+normal(3)*mytet(3,4)-dist
         
         ! Find cut case
         icase=1+int(0.5_WP+sign(0.5_WP,dd(1))) &
         &    +2*int(0.5_WP+sign(0.5_WP,dd(2))) &
         &    +4*int(0.5_WP+sign(0.5_WP,dd(3))) &
         &    +8*int(0.5_WP+sign(0.5_WP,dd(4)))
         
         ! Copy vertices
         vert(:,1:4)=mytet(:,1:4)
         
         ! Create interpolated vertices on cut plane
         do n1=1,cut_nvert(icase)
            v1=cut_v1(n1,icase); v2=cut_v2(n1,icase)
            mu=min(1.0_WP,max(0.0_WP,-dd(v1)/(sign(abs(dd(v2)-dd(v1))+epsilon(1.0_WP),dd(v2)-dd(v1)))))
            vert(:,4+n1)=(1.0_WP-mu)*vert(:,v1)+mu*vert(:,v2)
         end do
         
         ! Gas tets: from 1 to cut_nntet-1
         do n1=1,cut_nntet(icase)-1
            a=vert(:,cut_vtet(1,n1,icase))-vert(:,cut_vtet(4,n1,icase))
            b=vert(:,cut_vtet(2,n1,icase))-vert(:,cut_vtet(4,n1,icase))
            c=vert(:,cut_vtet(3,n1,icase))-vert(:,cut_vtet(4,n1,icase))
            my_vol=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
            bary=0.25_WP*(vert(:,cut_vtet(1,n1,icase))+vert(:,cut_vtet(2,n1,icase)) &
            &            +vert(:,cut_vtet(3,n1,icase))+vert(:,cut_vtet(4,n1,icase)))
            myVflux( 2 )=myVflux( 2 )+my_vol
            myVflux(6:8)=myVflux(6:8)+my_vol*bary
         end do
         
         ! Liquid tets: from cut_ntets down to cut_nntet
         do n1=cut_ntets(icase),cut_nntet(icase),-1
            a=vert(:,cut_vtet(1,n1,icase))-vert(:,cut_vtet(4,n1,icase))
            b=vert(:,cut_vtet(2,n1,icase))-vert(:,cut_vtet(4,n1,icase))
            c=vert(:,cut_vtet(3,n1,icase))-vert(:,cut_vtet(4,n1,icase))
            my_vol=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
            bary=0.25_WP*(vert(:,cut_vtet(1,n1,icase))+vert(:,cut_vtet(2,n1,icase)) &
            &            +vert(:,cut_vtet(3,n1,icase))+vert(:,cut_vtet(4,n1,icase)))
            myVflux( 1 )=myVflux( 1 )+my_vol
            myVflux(3:5)=myVflux(3:5)+my_vol*bary
         end do

         ! Compute Q flux from Qold
         myQflux(1)=myVflux(1)*pQold(i0,j0,k0,1)/(       VF0)
         myQflux(2)=myVflux(2)*pQold(i0,j0,k0,2)/(1.0_WP-VF0)
         myQflux(3)=myVflux(1)*pQold(i0,j0,k0,3)/(       VF0)
         myQflux(4)=myVflux(2)*pQold(i0,j0,k0,4)/(1.0_WP-VF0)
         myQflux(5:7)=sum(myQflux(1:2))*pQold(i0,j0,k0,5:7)/max(sum(pQold(i0,j0,k0,1:2)),this%rho_floor)
         
      end subroutine tet2flux_plic

      !> RK2 vertex projection back in time
      function project(p1,mydt) result(p2)
         implicit none
         real(WP), dimension(3), intent(in) :: p1
         real(WP), dimension(3)             :: p2
         real(WP),               intent(in) :: mydt
         p2=p1+mydt*interp_velocity(        p1    )
         p2=p1+mydt*interp_velocity(0.5_WP*(p1+p2))
      end function project

      !> Trilinear interpolation of collocated velocity - uses pU,pV,pW
      function interp_velocity(pos) result(vel)
         implicit none
         real(WP), dimension(3), intent(in) :: pos
         real(WP), dimension(3) :: vel
         integer  :: ipc,jpc,kpc
         real(WP) :: wxc1,wyc1,wzc1,wxc2,wyc2,wzc2
         ! All cell-centered
         ipc=floor((pos(1)-this%amr%xlo)*dxi-0.5_WP)
         jpc=floor((pos(2)-this%amr%ylo)*dyi-0.5_WP)
         kpc=floor((pos(3)-this%amr%zlo)*dzi-0.5_WP)
         ! Clamp to array bounds
         ipc=max(lbound(pU,1),min(ubound(pU,1)-1,ipc))
         jpc=max(lbound(pU,2),min(ubound(pU,2)-1,jpc))
         kpc=max(lbound(pU,3),min(ubound(pU,3)-1,kpc))
         ! Cell-centered weights
         wxc1=(pos(1)-(this%amr%xlo+(real(ipc,WP)+0.5_WP)*dx))*dxi
         wyc1=(pos(2)-(this%amr%ylo+(real(jpc,WP)+0.5_WP)*dy))*dyi
         wzc1=(pos(3)-(this%amr%zlo+(real(kpc,WP)+0.5_WP)*dz))*dzi
         wxc1=max(0.0_WP,min(1.0_WP,wxc1)); wxc2=1.0_WP-wxc1
         wyc1=max(0.0_WP,min(1.0_WP,wyc1)); wyc2=1.0_WP-wyc1
         wzc1=max(0.0_WP,min(1.0_WP,wzc1)); wzc2=1.0_WP-wzc1
         vel(1)=wzc1*(wyc1*(wxc1*pU(ipc+1,jpc+1,kpc+1,1)+wxc2*pU(ipc,jpc+1,kpc+1,1))+ &
         &            wyc2*(wxc1*pU(ipc+1,jpc  ,kpc+1,1)+wxc2*pU(ipc,jpc  ,kpc+1,1)))+&
         &      wzc2*(wyc1*(wxc1*pU(ipc+1,jpc+1,kpc  ,1)+wxc2*pU(ipc,jpc+1,kpc  ,1))+ &
         &            wyc2*(wxc1*pU(ipc+1,jpc  ,kpc  ,1)+wxc2*pU(ipc,jpc  ,kpc  ,1)))
         vel(2)=wzc1*(wyc1*(wxc1*pV(ipc+1,jpc+1,kpc+1,1)+wxc2*pV(ipc,jpc+1,kpc+1,1))+ &
         &            wyc2*(wxc1*pV(ipc+1,jpc  ,kpc+1,1)+wxc2*pV(ipc,jpc  ,kpc+1,1)))+&
         &      wzc2*(wyc1*(wxc1*pV(ipc+1,jpc+1,kpc  ,1)+wxc2*pV(ipc,jpc+1,kpc  ,1))+ &
         &            wyc2*(wxc1*pV(ipc+1,jpc  ,kpc  ,1)+wxc2*pV(ipc,jpc  ,kpc  ,1)))
         vel(3)=wzc1*(wyc1*(wxc1*pW(ipc+1,jpc+1,kpc+1,1)+wxc2*pW(ipc,jpc+1,kpc+1,1))+ &
         &            wyc2*(wxc1*pW(ipc+1,jpc  ,kpc+1,1)+wxc2*pW(ipc,jpc  ,kpc+1,1)))+&
         &      wzc2*(wyc1*(wxc1*pW(ipc+1,jpc+1,kpc  ,1)+wxc2*pW(ipc,jpc+1,kpc  ,1))+ &
         &            wyc2*(wxc1*pW(ipc+1,jpc  ,kpc  ,1)+wxc2*pW(ipc,jpc  ,kpc  ,1)))
      end function interp_velocity

   end subroutine get_dQdt

   !> Build PLIC reconstruction from VF and barycenters using PLICnet
   subroutine build_plic(this,time)
      use plicnet, only: get_normal, reflect_moments
      use mathtools, only: normalize
      use amrvof_geometry, only: get_plane_dist, cut_hex_polygon
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: time
      integer :: lvl
      real(WP) :: dx, dy, dz, dxi, dyi, dzi

      ! Only build at finest level
      lvl = this%amr%clvl()

      ! Get cell size at this level
      dx = this%amr%dx(lvl); dxi = 1.0_WP / dx
      dy = this%amr%dy(lvl); dyi = 1.0_WP / dy
      dz = this%amr%dz(lvl); dzi = 1.0_WP / dz

      ! Compute PLIC planes
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
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pPLIC2
         real(WP), dimension(3) :: lo2, hi2
         real(WP), dimension(4) :: plane2
         real(WP), dimension(3,8) :: hex2
         type(amrex_mfiter) :: mfi2
         type(amrex_box) :: bx2, gbx2
         
         ! Per-FAB polygon storage (allocatable, indexed by cell)
         real(WP), dimension(:,:,:,:,:), allocatable :: polygon_local  ! (3, 6, ilo:ihi, jlo:jhi, klo:khi)
         integer, dimension(:,:,:), allocatable :: poly_nv_local       ! (ilo:ihi, jlo:jhi, klo:khi)
         real(WP), dimension(3,6) :: poly_verts
         integer :: poly_nv
         integer :: glo(3), ghi(3)
         
         call this%amr%mfiter_build(lvl, mfi2)
         do while (mfi2%next())
            bx2 = mfi2%tilebox()
            gbx2 = mfi2%growntilebox(2)  ! Grown by 2 for 5x5x5 stencil
            pPLIC2 => this%PLIC%mf(lvl)%dataptr(mfi2)
            
            ! Get bounds for per-FAB allocation
            glo = gbx2%lo
            ghi = gbx2%hi
            
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
                     if (abs(pPLIC2(ii,jj,kk,4)) .gt. 1.0e+9_WP) cycle
                     
                     ! Build hex and plane
                     lo2 = [this%amr%xlo + real(ii  ,WP)*dx, this%amr%ylo + real(jj  ,WP)*dy, this%amr%zlo + real(kk  ,WP)*dz]
                     hi2 = [this%amr%xlo + real(ii+1,WP)*dx, this%amr%ylo + real(jj+1,WP)*dy, this%amr%zlo + real(kk+1,WP)*dz]
                     plane2 = [pPLIC2(ii,jj,kk,1), pPLIC2(ii,jj,kk,2), pPLIC2(ii,jj,kk,3), pPLIC2(ii,jj,kk,4)]
                     hex2(:,1) = [hi2(1), lo2(2), lo2(3)]
                     hex2(:,2) = [hi2(1), hi2(2), lo2(3)]
                     hex2(:,3) = [hi2(1), hi2(2), hi2(3)]
                     hex2(:,4) = [hi2(1), lo2(2), hi2(3)]
                     hex2(:,5) = [lo2(1), lo2(2), lo2(3)]
                     hex2(:,6) = [lo2(1), hi2(2), lo2(3)]
                     hex2(:,7) = [lo2(1), hi2(2), hi2(3)]
                     hex2(:,8) = [lo2(1), lo2(2), hi2(3)]
                     
                     call cut_hex_polygon(hex2, plane2, poly_nv, poly_verts)
                     
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
            do k = bx2%lo(3), bx2%hi(3)
               do j = bx2%lo(2), bx2%hi(2)
                  do i = bx2%lo(1), bx2%hi(1)
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
         call this%amr%mfiter_destroy(mfi2)
      end block polygon_extraction
      
   end subroutine build_plic

   !> Reset VF and barycenters from PLIC plane to ensure consistency
   !> Computes in valid + ghost cells from PLIC (which is already filled)
   subroutine reset_moments(this)
      use amrvof_geometry, only: cut_hex_vol
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      class(amrmpcomp), intent(inout) :: this
      integer :: lvl,i,j,k
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCliq,pCgas,pPLIC,pQ
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
         pQ   =>this%Q%mf(lvl)%dataptr(mfi)
         ! Loop over grown tiles
         bx=mfi%growntilebox(this%nover)
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
                  pQ(i,j,k,2)=0.0_WP
                  pQ(i,j,k,4)=0.0_WP
               else
                  pVF(i,j,k,1)=0.0_WP
                  pQ(i,j,k,1)=0.0_WP
                  pQ(i,j,k,3)=0.0_WP
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
               pQ(i,j,k,1)=0.0_WP
               pQ(i,j,k,3)=0.0_WP
            end if
            if (pVF(i,j,k,1).gt.VFhi) then
               pVF(i,j,k,1)=1.0_WP
               pCliq(i,j,k,1:3)=cell_center
               pCgas(i,j,k,1:3)=cell_center
               pQ(i,j,k,2)=0.0_WP
               pQ(i,j,k,4)=0.0_WP
            end if
         end do; end do; end do
      end do
      call this%amr%mfiter_destroy(mfi)
      
      ! Average down VF/Cliq/Cgas to coarse levels + clean up PLIC
      call this%vof_average_down()
      
   end subroutine reset_moments

   ! ============================================================================
   ! VOF SYNC / FILL / AVERAGE_DOWN METHODS (from amrvof)
   ! ============================================================================

   !> Sync VF/Cliq/Cgas ghosts on all levels + fix periodic barycenters
   subroutine sync_moments(this)
      class(amrmpcomp), intent(inout) :: this
      integer :: lvl
      do lvl=0,this%amr%clvl()
         call this%sync_moments_lvl(lvl)
      end do
   end subroutine sync_moments

   !> Sync VF/Cliq/Cgas ghosts at level + fix periodic barycenters
   subroutine sync_moments_lvl(this, lvl)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_box, amrex_geometry
      class(amrmpcomp), intent(inout) :: this
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

   !> Sync PLIC ghosts on all levels + fix periodic plane distance
   subroutine sync_plic(this)
      class(amrmpcomp), intent(inout) :: this
      integer :: lvl
      do lvl = 0, this%amr%clvl()
         call this%sync_plic_lvl(lvl)
      end do
   end subroutine sync_plic

   !> Sync PLIC ghosts at level + fix periodic plane distance
   !> d <- d +/- n*L where L is domain length and n is normal component
   subroutine sync_plic_lvl(this, lvl)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_geometry
      implicit none
      class(amrmpcomp), intent(inout) :: this
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

   !> Fill PLIC ghosts at a level (fill + physical BC)
   subroutine fill_plic_lvl(this, lvl, time)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_geometry, amrex_box
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_mfiter) :: mfi
      type(amrex_geometry) :: geom
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas, pPLIC
      integer :: ig, jg, kg, ilo, ihi, jlo, jhi, klo, khi
      integer :: dlo(3), dhi(3)
      real(WP) :: xL, yL, zL, dx, dy, dz
      
      ! Sync ghosts (periodic + MPI exchange)
      call this%PLIC%fill_lvl(lvl, time)
      
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
      subroutine apply_bc_face(dir, side, bc_type, i1, i2, j1, j2, k1, k2, bnd, x_bnd)
         integer, intent(in) :: dir, side, bc_type, i1, i2, j1, j2, k1, k2, bnd
         real(WP), intent(in) :: x_bnd
         integer :: ig2, jg2, kg2, isrc, jsrc, ksrc, face
         type(amrex_box) :: bc_bx
         
         select case (bc_type)
         
          case (BC_LIQ)
            ! Trivial PLIC: full liquid
            do kg2 = k1, k2; do jg2 = j1, j2; do ig2 = i1, i2
               pPLIC(ig2,jg2,kg2,1:3) = 0.0_WP
               pPLIC(ig2,jg2,kg2,4) = 1.0e10_WP
            end do; end do; end do
            
          case (BC_GAS)
            ! Trivial PLIC: full gas
            do kg2 = k1, k2; do jg2 = j1, j2; do ig2 = i1, i2
               pPLIC(ig2,jg2,kg2,1:3) = 0.0_WP
               pPLIC(ig2,jg2,kg2,4) = -1.0e10_WP
            end do; end do; end do
            
          case (BC_REFLECT)
            ! Mirror PLIC from interior + flip normal component
            do kg2 = k1, k2; do jg2 = j1, j2; do ig2 = i1, i2
               isrc = ig2; jsrc = jg2; ksrc = kg2
               if (dir.eq.1) isrc = 2*bnd - ig2 - side
               if (dir.eq.2) jsrc = 2*bnd - jg2 - side
               if (dir.eq.3) ksrc = 2*bnd - kg2 - side
               ! Copy plane
               pPLIC(ig2,jg2,kg2,1:4) = pPLIC(isrc,jsrc,ksrc,1:4)
               ! Flip normal component
               pPLIC(ig2,jg2,kg2,dir) = -pPLIC(ig2,jg2,kg2,dir)
               ! Correct plane distance
               pPLIC(ig2,jg2,kg2,4) = pPLIC(ig2,jg2,kg2,4) - 2.0_WP*pPLIC(isrc,jsrc,ksrc,dir)*x_bnd
            end do; end do; end do
            
          case (BC_USER)
            ! User callback sets PLIC
            if (associated(this%user_vof_bc)) then
               bc_bx = amrex_box([i1, j1, k1], [i2, j2, k2])
               face = 2*dir - 1 + (1+side)/2
               call this%user_vof_bc(this, bc_bx, pVF, pCliq, pCgas, pPLIC, face, time, 1)
            end if
            
          case default
            ! Do nothing
            
         end select
         
      end subroutine apply_bc_face
      
   end subroutine fill_plic_lvl

   !> Fill moments ghosts at a level (fill + physical BC)
   subroutine fill_moments_lvl(this, lvl, time)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_box, amrex_geometry
      class(amrmpcomp), intent(inout) :: this
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
      subroutine apply_bc_face(dir, side, bc_type, i1, i2, j1, j2, k1, k2, bnd, x_bnd)
         integer, intent(in) :: dir, side, bc_type, i1, i2, j1, j2, k1, k2, bnd
         real(WP), intent(in) :: x_bnd
         integer :: ig2, jg2, kg2, isrc, jsrc, ksrc, face
         real(WP), dimension(3) :: center
         type(amrex_box) :: bc_bx
         
         select case (bc_type)
         
          case (BC_LIQ)
            ! Full liquid: VF=1, barycenters at cell center
            do kg2 = k1, k2; do jg2 = j1, j2; do ig2 = i1, i2
               center = [this%amr%xlo + (real(ig2,WP)+0.5_WP)*dx, &
               &         this%amr%ylo + (real(jg2,WP)+0.5_WP)*dy, &
               &         this%amr%zlo + (real(kg2,WP)+0.5_WP)*dz]
               pVF(ig2,jg2,kg2,1) = 1.0_WP
               pCliq(ig2,jg2,kg2,1:3) = center
               pCgas(ig2,jg2,kg2,1:3) = center
            end do; end do; end do
            
          case (BC_GAS)
            ! Full gas: VF=0, barycenters at cell center
            do kg2 = k1, k2; do jg2 = j1, j2; do ig2 = i1, i2
               center = [this%amr%xlo + (real(ig2,WP)+0.5_WP)*dx, &
               &         this%amr%ylo + (real(jg2,WP)+0.5_WP)*dy, &
               &         this%amr%zlo + (real(kg2,WP)+0.5_WP)*dz]
               pVF(ig2,jg2,kg2,1) = 0.0_WP
               pCliq(ig2,jg2,kg2,1:3) = center
               pCgas(ig2,jg2,kg2,1:3) = center
            end do; end do; end do
            
          case (BC_REFLECT)
            ! Mirror VF/Cliq/Cgas from interior + reflect barycenters
            do kg2 = k1, k2; do jg2 = j1, j2; do ig2 = i1, i2
               isrc = ig2; jsrc = jg2; ksrc = kg2
               if (dir.eq.1) isrc = 2*bnd - ig2 - side
               if (dir.eq.2) jsrc = 2*bnd - jg2 - side
               if (dir.eq.3) ksrc = 2*bnd - kg2 - side
               ! Copy VF
               pVF(ig2,jg2,kg2,1) = pVF(isrc,jsrc,ksrc,1)
               ! Copy and reflect barycenters
               pCliq(ig2,jg2,kg2,1:3) = pCliq(isrc,jsrc,ksrc,1:3)
               pCgas(ig2,jg2,kg2,1:3) = pCgas(isrc,jsrc,ksrc,1:3)
               pCliq(ig2,jg2,kg2,dir) = 2.0_WP*x_bnd - pCliq(isrc,jsrc,ksrc,dir)
               pCgas(ig2,jg2,kg2,dir) = 2.0_WP*x_bnd - pCgas(isrc,jsrc,ksrc,dir)
            end do; end do; end do
            
          case (BC_USER)
            ! User callback sets moments
            if (associated(this%user_vof_bc)) then
               bc_bx = amrex_box([i1, j1, k1], [i2, j2, k2])
               face = 2*dir - 1 + (1+side)/2
               call this%user_vof_bc(this, bc_bx, pVF, pCliq, pCgas, pPLIC, face, time, 2)
            end if
            
          case default
            ! Do nothing
            
         end select
         
      end subroutine apply_bc_face
      
   end subroutine fill_moments_lvl

   !> Average down VF/Cliq/Cgas from finest to lbase, then sync ghost cells
   !> Clean up PLIC at coarse levels and sync ghost cells
   subroutine vof_average_down(this, lbase)
      use amrex_interface, only: amrmfab_average_down_cell
      class(amrmpcomp), intent(inout) :: this
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
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF2, pPLIC2
         integer :: i, j, k
         call amrex_mfiter_build(mfi, this%PLIC%mf(lvl), tiling=.false.)
         do while (mfi%next())
            bx = mfi%tilebox()
            pVF2 => this%VF%mf(lvl)%dataptr(mfi)
            pPLIC2 => this%PLIC%mf(lvl)%dataptr(mfi)
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     pPLIC2(i,j,k,1:3) = 0.0_WP
                     pPLIC2(i,j,k,4) = sign(1.0e10_WP, pVF2(i,j,k,1) - 0.5_WP)
                  end do
               end do
            end do
         end do
         call amrex_mfiter_destroy(mfi)
      end subroutine set_trivial_plic
   end subroutine vof_average_down

   !> Apply pressure relaxation to mixture cells
   subroutine apply_relax(this)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use amrvof_geometry,  only: get_plane_dist,cut_hex_vol
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pQ,pPLIC,pCliq,pCgas
      real(WP), dimension(3) :: lo,hi,cell_center,bary_liq,bary_gas
      real(WP), dimension(3,8) :: hex
      real(WP), dimension(4) :: plane
      real(WP) :: dx,dy,dz,vol_liq,vol_gas
      ! If no relaxation model was provided, return
      if (.not.associated(this%relax)) return
      ! Apply relaxation on finest level only (mixture cells are always at finest)
      lvl=this%amr%clvl()
      dx=this%amr%dx(lvl); dy=this%amr%dy(lvl); dz=this%amr%dz(lvl)
      call this%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         pVF  =>this%VF%mf(lvl)%dataptr(mfi)
         pQ   =>this%Q%mf(lvl)%dataptr(mfi)
         pPLIC=>this%PLIC%mf(lvl)%dataptr(mfi)
         pCliq=>this%Cliq%mf(lvl)%dataptr(mfi)
         pCgas=>this%Cgas%mf(lvl)%dataptr(mfi)
         ! Loop over grown tiles
         bx=mfi%growntilebox(this%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Only relax mixture cells
            if (pVF(i,j,k,1).lt.VFlo.or.pVF(i,j,k,1).gt.VFhi) cycle
            ! Apply user-provided relaxation model (modifies VF and Q)
            call this%relax(pVF(i,j,k,1),pQ(i,j,k,:))
            ! Adjust PLIC plane to match new VF
            lo=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hi=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            cell_center=0.5_WP*(lo+hi)
            ! Reposition plane: keep normal, adjust distance for new VF
            pPLIC(i,j,k,4)=get_plane_dist(pPLIC(i,j,k,1:3),lo,hi,pVF(i,j,k,1))
            ! Recompute barycenters from adjusted PLIC
            hex(:,1)=[lo(1),lo(2),lo(3)]; hex(:,2)=[hi(1),lo(2),lo(3)]
            hex(:,3)=[hi(1),hi(2),lo(3)]; hex(:,4)=[lo(1),hi(2),lo(3)]
            hex(:,5)=[lo(1),lo(2),hi(3)]; hex(:,6)=[hi(1),lo(2),hi(3)]
            hex(:,7)=[hi(1),hi(2),hi(3)]; hex(:,8)=[lo(1),hi(2),hi(3)]
            plane=pPLIC(i,j,k,1:4)
            call cut_hex_vol(hex,plane,vol_liq,vol_gas,bary_liq,bary_gas)
            pVF(i,j,k,1)=vol_liq/(dx*dy*dz)
            pCliq(i,j,k,1:3)=bary_liq
            pCgas(i,j,k,1:3)=bary_gas
            ! Handle newly-pure cells from relaxation
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
   end subroutine apply_relax
   
   !> Compute artificial bulk viscosity
   subroutine get_viscartif(this,dt,beta)
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_multifab,amrex_multifab_destroy
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      class(amrdata), intent(inout) :: beta
      ! Local variables
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab) :: div
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pBeta,pDiv,pU,pV,pW,pC
      real(WP) :: dxi,dyi,dzi,dx,dy,dz,max_beta
      real(WP) :: dudy,dudz,dvdx,dvdz,dwdx,dwdy,vort,grad_div
      integer :: lvl,i,j,k,si,sj,sk,n
      ! Parameters
      real(WP), parameter :: Cartif=2.0_WP
      real(WP), parameter :: max_cfl=0.5_WP
      integer, parameter :: nfilter=2
      real(WP), dimension(-1:+1), parameter :: filter=[1.0_WP/6.0_WP,2.0_WP/3.0_WP,1.0_WP/6.0_WP]
      
      ! Zero out beta
      call beta%setval(val=0.0_WP)
      
      ! Loop over levels
      do lvl=0,this%amr%clvl()
         
         ! Grid spacings
         dx=this%amr%dx(lvl); dxi=1.0_WP/dx
         dy=this%amr%dy(lvl); dyi=1.0_WP/dy
         dz=this%amr%dz(lvl); dzi=1.0_WP/dz
         
         ! Max beta from CFL
         max_beta=max_cfl*min(dx**2,dy**2,dz**2)/(4.0_WP*dt)
         
         ! Build temp multifab for divergence at this level (1 ghost cell)
         call this%amr%mfab_build(lvl=lvl,mfab=div,ncomp=1,nover=1); call div%setval(0.0_WP)
         
         ! Phase 1: Compute divergence
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            ! Get data pointers
            pU=>this%U%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi)
            pDiv=>div%dataptr(mfi)
            ! Loop over tiles grown by 1 (safe as velocity has 2 ghost cells)
            bx=mfi%growntilebox(1)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pDiv(i,j,k,1)=0.5_WP*(dxi*(pU(i+1,j,k,1)-pU(i-1,j,k,1))+dyi*(pV(i,j+1,k,1)-pV(i,j-1,k,1))+dzi*(pW(i,j,k+1,1)-pW(i,j,k-1,1)))
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         
         ! Phase 2: Compute beta
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            ! Get data pointers
            pU=>this%U%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi)
            pC=>this%C%mf(lvl)%dataptr(mfi)
            pDiv=>div%dataptr(mfi)
            pBeta=>beta%mf(lvl)%dataptr(mfi)
            ! Loop over interior tiles
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Only work in compression regions
               if (pDiv(i,j,k,1).ge.0.0_WP) cycle
               ! Compute local vorticity
               dudy=0.5_WP*dyi*(pU(i,j+1,k,1)-pU(i,j-1,k,1))
               dudz=0.5_WP*dzi*(pU(i,j,k+1,1)-pU(i,j,k-1,1))
               dvdx=0.5_WP*dxi*(pV(i+1,j,k,1)-pV(i-1,j,k,1))
               dvdz=0.5_WP*dzi*(pV(i,j,k+1,1)-pV(i,j,k-1,1))
               dwdx=0.5_WP*dxi*(pW(i+1,j,k,1)-pW(i-1,j,k,1))
               dwdy=0.5_WP*dyi*(pW(i,j+1,k,1)-pW(i,j-1,k,1))
               vort=(dwdy-dvdz)**2+(dudz-dwdx)**2+(dvdx-dudy)**2
               ! Compute |grad(div)|
               grad_div=max(abs(pDiv(i+1,j,k,1)-pDiv(i,j,k,1)),abs(pDiv(i,j,k,1)-pDiv(i-1,j,k,1)))*dx**2 &
               &       +max(abs(pDiv(i,j+1,k,1)-pDiv(i,j,k,1)),abs(pDiv(i,j,k,1)-pDiv(i,j-1,k,1)))*dy**2 &
               &       +max(abs(pDiv(i,j,k+1,1)-pDiv(i,j,k,1)),abs(pDiv(i,j,k,1)-pDiv(i,j,k-1,1)))*dz**2
               ! Floor vorticity with sound speed
               vort=max(vort,(0.05_WP*pC(i,j,k,1)/min(dx,dy,dz))**2)
               ! Compute beta
               pBeta(i,j,k,1)=Cartif*grad_div*min(4.0_WP/3.0_WP*pDiv(i,j,k,1)**2/(pDiv(i,j,k,1)**2+vort+1.0e-15_WP),1.0_WP)
               ! Clip to max
               pBeta(i,j,k,1)=min(pBeta(i,j,k,1),max_beta)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         
         ! Phase 3: Filter beta using div as scratch
         do n=1,nfilter
            ! Reset div to zero
            call div%setval(0.0_WP)
            ! Copy beta to div
            call div%copy(srcmf=beta%mf(lvl),srccomp=1,dstcomp=1,nc=1,ng=0)
            ! Fill internal/periodic ghosts
            call div%fill_boundary(this%amr%geom(lvl))
            ! Apply Neumann BCs
            call this%amr%mfab_foextrap(lvl=lvl,mfab=div)
            ! Apply filter
            call this%amr%mfiter_build(lvl,mfi)
            do while(mfi%next())
               ! Get data pointers
               pDiv=>div%dataptr(mfi)
               pBeta=>beta%mf(lvl)%dataptr(mfi)
               ! Loop over interior tiles
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pBeta(i,j,k,1)=0.0_WP
                  do sk=-1,+1; do sj=-1,+1; do si=-1,+1
                     pBeta(i,j,k,1)=pBeta(i,j,k,1)+filter(si)*filter(sj)*filter(sk)*pDiv(i+si,j+sj,k+sk,1)
                  end do; end do; end do
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do

         ! Destroy temp divergence multifab
         call amrex_multifab_destroy(div)

         ! Fill beta's ghosts at this level
         call beta%mf(lvl)%fill_boundary(this%amr%geom(lvl))
         call this%amr%mfab_foextrap(lvl=lvl,mfab=beta%mf(lvl))

      end do

   end subroutine get_viscartif

   !> Calculate CFL numbers
   subroutine get_cfl(this,dt,cfl)
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(out) :: cfl
      integer :: lvl
      real(WP) :: Umax,Vmax,Wmax,Cmax,viscmax
      ! Reset CFLs
      this%CFLc_x=0.0_WP; this%CFLc_y=0.0_WP; this%CFLc_z=0.0_WP
      this%CFLa_x=0.0_WP; this%CFLa_y=0.0_WP; this%CFLa_z=0.0_WP
      this%CFLv_x=0.0_WP; this%CFLv_y=0.0_WP; this%CFLv_z=0.0_WP
      ! Compute CFL at each level (finest level determines dt)
      do lvl=0,this%amr%clvl()
         ! Max velocity
         Umax=this%U%norm0(lvl=lvl)
         Vmax=this%V%norm0(lvl=lvl)
         Wmax=this%W%norm0(lvl=lvl)
         ! Max speed of sound
         Cmax=this%C%norm0(lvl=lvl)
         ! Max viscosities (TODO: use phasic viscosities viscL/RHOL, viscG/RHOG instead of mixture)
         get_viscmax: block
            use amrex_amr_module, only: amrex_mfiter,amrex_box
            use parallel, only: MPI_REAL_WP
            use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_MAX
            type(amrex_mfiter) :: mfi
            type(amrex_box) :: bx
            real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVisc,pBeta,pDiff
            integer :: i,j,k,ierr
            real(WP) :: rho
            viscmax=0.0_WP
            call this%amr%mfiter_build(lvl,mfi)
            do while(mfi%next())
               ! Get data pointers
               pQ=>this%Q%mf(lvl)%dataptr(mfi)
               pVisc=>this%visc%mf(lvl)%dataptr(mfi)
               pBeta=>this%beta%mf(lvl)%dataptr(mfi)
               pDiff=>this%diff%mf(lvl)%dataptr(mfi)
               ! Loop over interior tiles
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  rho=max(pQ(i,j,k,1)+pQ(i,j,k,2),this%rho_floor)
                  viscmax=max(viscmax,pVisc(i,j,k,1)/rho,pBeta(i,j,k,1)/rho,pDiff(i,j,k,1)/rho)
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
            call MPI_ALLREDUCE(MPI_IN_PLACE,viscmax,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         end block get_viscmax
         ! Convective+acoustic
         this%CFLc_x=max(this%CFLc_x,(Umax+Cmax)*dt/this%amr%dx(lvl))
         this%CFLc_y=max(this%CFLc_y,(Vmax+Cmax)*dt/this%amr%dy(lvl))
         this%CFLc_z=max(this%CFLc_z,(Wmax+Cmax)*dt/this%amr%dz(lvl))
         ! Acoustic
         this%CFLa_x=max(this%CFLa_x,Cmax*dt/this%amr%dx(lvl))
         this%CFLa_y=max(this%CFLa_y,Cmax*dt/this%amr%dy(lvl))
         this%CFLa_z=max(this%CFLa_z,Cmax*dt/this%amr%dz(lvl))
         ! Viscous
         this%CFLv_x=max(this%CFLv_x,4.0_WP*viscmax*dt/this%amr%dx(lvl)**2)
         this%CFLv_y=max(this%CFLv_y,4.0_WP*viscmax*dt/this%amr%dy(lvl)**2)
         this%CFLv_z=max(this%CFLv_z,4.0_WP*viscmax*dt/this%amr%dz(lvl)**2)
      end do
      ! Return max CFL
      cfl=max(this%CFLc_x,this%CFLc_y,this%CFLc_z,this%CFLa_x,this%CFLa_y,this%CFLa_z,this%CFLv_x,this%CFLv_y,this%CFLv_z)
   end subroutine get_cfl

   !> Calculate monitoring info
   subroutine get_info(this)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer :: lvl,n
      real(WP) :: dV

      ! Compute extrema across all levels
      this%Umax=0.0_WP; this%Vmax=0.0_WP; this%Wmax=0.0_WP
      this%ILmin=huge(1.0_WP); this%ILmax=-huge(1.0_WP); this%IGmin=huge(1.0_WP); this%IGmax=-huge(1.0_WP)
      this%PLmin=huge(1.0_WP); this%PLmax=-huge(1.0_WP); this%PGmin=huge(1.0_WP); this%PGmax=-huge(1.0_WP)
      this%TLmin=huge(1.0_WP); this%TLmax=-huge(1.0_WP); this%TGmin=huge(1.0_WP); this%TGmax=-huge(1.0_WP)
      this%Cmin=huge(1.0_WP); this%Cmax=-huge(1.0_WP)
      this%VFmin=huge(1.0_WP); this%VFmax=-huge(1.0_WP); this%VFint=0.0_WP
      this%Qmin=huge(1.0_WP); this%Qmax=-huge(1.0_WP)
      do lvl=0,this%amr%clvl()
         ! Velocity norm 0
         this%Umax=max(this%Umax,this%U%norm0(lvl=lvl))
         this%Vmax=max(this%Vmax,this%V%norm0(lvl=lvl))
         this%Wmax=max(this%Wmax,this%W%norm0(lvl=lvl))
         ! Extrema of phasic internal energy, pressure, and temperature
         this%ILmin=min(this%ILmin,this%IL%get_min(lvl=lvl)); this%ILmax=max(this%ILmax,this%IL%get_max(lvl=lvl))
         this%IGmin=min(this%IGmin,this%IG%get_min(lvl=lvl)); this%IGmax=max(this%IGmax,this%IG%get_max(lvl=lvl))
         this%PLmin=min(this%PLmin,this%PL%get_min(lvl=lvl)); this%PLmax=max(this%PLmax,this%PL%get_max(lvl=lvl))
         this%PGmin=min(this%PGmin,this%PG%get_min(lvl=lvl)); this%PGmax=max(this%PGmax,this%PG%get_max(lvl=lvl))
         this%TLmin=min(this%TLmin,this%TL%get_min(lvl=lvl)); this%TLmax=max(this%TLmax,this%TL%get_max(lvl=lvl))
         this%TGmin=min(this%TGmin,this%TG%get_min(lvl=lvl)); this%TGmax=max(this%TGmax,this%TG%get_max(lvl=lvl))
         ! Extrema of mixture speed of sound
         this%Cmin=min(this%Cmin,this%C%get_min(lvl=lvl)); this%Cmax=max(this%Cmax,this%C%get_max(lvl=lvl))
         ! Extrema of volume fraction
         this%VFmin=min(this%VFmin,this%VF%get_min(lvl=lvl)); this%VFmax=max(this%VFmax,this%VF%get_max(lvl=lvl))
         ! Extrema of conserved variables
         do n=1,this%Q%ncomp
            this%Qmin(n)=min(this%Qmin(n),this%Q%get_min(lvl=lvl,comp=n))
            this%Qmax(n)=max(this%Qmax(n),this%Q%get_max(lvl=lvl,comp=n))
         end do
      end do

      ! Conserved integrals at base level
      dV=this%amr%dx(0)*this%amr%dy(0)*this%amr%dz(0)
      do n=1,this%Q%ncomp
         this%Qint(n)=this%Q%get_sum(lvl=0,comp=n)*dV
      end do

      ! Kinetic energy integral: 0.5 * rho * (U^2 + V^2 + W^2) * dV
      ! Uses composite integration with fine masking to avoid double-counting
      get_rhoKint: block
         use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface, only: amrmask_make_fine
         use parallel, only: MPI_REAL_WP
         use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         type(amrex_imultifab) :: mask
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pU,pV,pW
         integer, dimension(:,:,:,:), contiguous, pointer :: pMask
         integer :: i,j,k,ierr
         this%rhoKint=0.0_WP
         do lvl=0,this%amr%clvl()
            ! Get cell volume
            dV=this%amr%dx(lvl)*this%amr%dy(lvl)*this%amr%dz(lvl)
            ! Build fine mask for this level (if not finest)
            if (lvl.lt.this%amr%clvl()) then
               call amrex_imultifab_build(mask,this%amr%ba(lvl),this%amr%dm(lvl),1,0)
               call amrmask_make_fine(mask,this%amr%ba(lvl+1),[this%amr%rref(lvl),this%amr%rref(lvl),this%amr%rref(lvl)],0,1)
            end if
            ! Loop over tiles
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               bx=mfi%tilebox()
               ! Get pointers to data
               pQ=>this%Q%mf(lvl)%dataptr(mfi)
               pU=>this%U%mf(lvl)%dataptr(mfi)
               pV=>this%V%mf(lvl)%dataptr(mfi)
               pW=>this%W%mf(lvl)%dataptr(mfi)
               ! Get pointer to fine mask
               if (lvl.lt.this%amr%clvl()) pMask=>mask%dataptr(mfi)
               ! Loop over cells
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Skip cells covered by finer level
                  if (lvl.lt.this%amr%clvl()) then; if (pMask(i,j,k,1).eq.0) cycle; end if
                  ! Accumulate kinetic energy (rho = Q(1)+Q(2) = VF*rhoL + (1-VF)*rhoG)
                  this%rhoKint=this%rhoKint+0.5_WP*(pQ(i,j,k,1)+pQ(i,j,k,2))*(pU(i,j,k,1)**2+pV(i,j,k,1)**2+pW(i,j,k,1)**2)*dV
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
            if (lvl.lt.this%amr%clvl()) call amrex_imultifab_destroy(mask)
         end do
         ! Reduce across MPI ranks
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%rhoKint,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
      end block get_rhoKint
      
   end subroutine get_info

   !> Print solver info to screen
   subroutine amrmpcomp_print(this)
      use messager, only: log
      implicit none
      class(amrmpcomp), intent(in) :: this
      call log("Compressible Multiphase solver: "//trim(this%name))
      call log("  Grid: "//trim(this%amr%name))
   end subroutine amrmpcomp_print

   !> Register solver data for checkpoint
   subroutine register_checkpoint(this,io)
      use amrio_class, only: amrio
      implicit none
      class(amrmpcomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      call io%add_data(this%Q,'Q')
      call io%add_data(this%VF,'VF')
      call io%add_data(this%Cliq,'Cliq')
      call io%add_data(this%Cgas,'Cgas')
   end subroutine register_checkpoint

   !> Restore solver data from checkpoint
   subroutine restore_checkpoint(this,io,dirname)
      use amrio_class, only: amrio
      implicit none
      class(amrmpcomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      call io%read_data(dirname,this%Q,'Q')
      call io%read_data(dirname,this%VF,'VF')
      call io%read_data(dirname,this%Cliq,'Cliq')
      call io%read_data(dirname,this%Cgas,'Cgas')
   end subroutine restore_checkpoint
end module amrmpcomp_class
