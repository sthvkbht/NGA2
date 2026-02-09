!> AMR compressible drop test case
module simulation
   use precision,         only: WP
   use amrgrid_class,     only: amrgrid
   use amrmpcomp_class,   only: amrmpcomp
   use amrviz_class,      only: amrviz
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> AMR grid
   type(amrgrid), target :: amr

   !> Timetracker and compressible multiphase solver
   type(timetracker) :: time
   type(amrmpcomp), target :: fs
   type(amrdata) :: dQdt,rho_mix,Umag,Mach
   
   !> Visualization
   type(event) :: viz_evt
   type(amrviz) :: viz

   ! Regrid parameters
   type(event) :: regrid_evt
   
   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile
   
   !> Stiffened gas EOS parameters (liquid and gas)
   real(WP) :: GammaL,PinfL,CvL
   real(WP) :: GammaG,PinfG,CvG

   !> Flow parameters
   real(WP) :: rhoG1,pG1,u1           !< Pre-shock gas state
   real(WP) :: rhoG2,pG2,u2           !< Post-shock gas state
   real(WP) :: rhoL1,pL1              !< Initial liquid state
   real(WP) :: M2,Xs                  !< Post-shock Mach and shock location
   real(WP) :: Ms                     !< Shock Mach number
   real(WP) :: density_ratio          !< rhoL1/rhoG1
   real(WP) :: sound_ratio            !< cL1/cG1
   real(WP) :: Reynolds,visc_ratio    !< Viscosity 
   real(WP) :: Prandtl ,diff_ratio    !< Heat diffusivity
   
   !> Sutherland viscosity parameters: mu_g = (1+Suth_T)*T^Suth_n / (Re*(T+Suth_T))
   real(WP) :: Suth_n=1.5_WP          !< Sutherland exponent (1.0 for constant)
   real(WP) :: Suth_T=0.4042_WP       !< Sutherland temperature (0.0 for constant)

   !> Sponge parameters
   real(WP) :: R_spg=3.0_WP
   real(WP) :: L_spg=1.0_WP
   real(WP) :: nu_spg=0.01_WP

   !> Tagging parameters
   real(WP) :: Rec_tag=huge(1.0_WP)
   real(WP) :: Res_tag=huge(1.0_WP)
   
contains

   !> Smooth Heaviside function
   real(WP) function Hshock(x,delta)
      real(WP), intent(in) :: x,delta
      Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock

   !> Levelset function for sphere (centered at origin)
   function sphere_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=0.5_WP-sqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)
   end function sphere_levelset

   !> Liquid EOS: P=f(RHO,I) - Stiffened gas
   pure real(WP) function get_PL(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_PL=RHO*I*(GammaL-1.0_WP)-GammaL*PinfL
   end function get_PL
   !> Liquid EOS: T=f(RHO,P)
   pure real(WP) function get_TL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_TL=(P+PinfL)/(CvL*RHO*(GammaL-1.0_WP))
   end function get_TL
   !> Liquid EOS: C=f(RHO,P)
   pure real(WP) function get_CL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_CL=sqrt(GammaL*(P+PinfL)/RHO)
   end function get_CL
   !> Liquid EOS: I=f(RHO,P) (used for initialization)
   pure real(WP) function get_IL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_IL=(P+GammaL*PinfL)/(RHO*(GammaL-1.0_WP))
   end function get_IL

   !> Gas EOS: P=f(RHO,I) - Ideal gas
   pure real(WP) function get_PG(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_PG=RHO*I*(GammaG-1.0_WP)-GammaG*PinfG
   end function get_PG
   !> Gas EOS: T=f(RHO,P)
   pure real(WP) function get_TG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_TG=(P+PinfG)/(CvG*RHO*(GammaG-1.0_WP))
   end function get_TG
   !> Gas EOS: C=f(RHO,P)
   pure real(WP) function get_CG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_CG=sqrt(GammaG*(P+PinfG)/RHO)
   end function get_CG
   !> Gas EOS: I=f(RHO,P) (used for initialization)
   pure real(WP) function get_IG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_IG=(P+GammaG*PinfG)/(RHO*(GammaG-1.0_WP))
   end function get_IG

   !> Implicit mechanical relaxation for stiffened gas EOS pair
   !> Solves quadratic for equilibrium pressure Peq where PL=PG=Peq,
   !> then adjusts VF and internal energies via p*dV work exchange.
   !> Conserves: phasic masses Q(1:2), total internal energy Q(3)+Q(4), momentum Q(5:7)
   subroutine P_relax_implicit(VF,Q)
      implicit none
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP) :: a,b,d,Peq,VFeq
      real(WP) :: invG1G,invG1L,d0,d1,facG,facL
      real(WP), parameter :: RHOGmin=1.0e-3_WP
      ! Skip near-pure-liquid cells (gas density too low)
      if (Q(2)/(1.0_WP-VF).lt.RHOGmin) return
      ! Precompute EOS constants
      invG1G=1.0_WP/(GammaG-1.0_WP)
      invG1L=1.0_WP/(GammaL-1.0_WP)
      d0=GammaL*PinfL*invG1L
      d1=1.0_WP+invG1L
      facG=GammaG*PinfG*invG1G
      facL=invG1G+VF
      ! Quadratic coefficients: a*Peq^2 + b*Peq + d = 0
      a=d1*facL-VF*(invG1G+1.0_WP)
      b=d1*(facG-Q(4))-VF*facG+d0*facL-Q(3)*(invG1G+1.0_WP)
      d=d0*(facG-Q(4))-Q(3)*facG
      ! Solve for equilibrium pressure (positive root)
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Bail if pressure is unphysical
      if (Peq.le.max(-PinfG,-PinfL)) return
      ! Equilibrium volume fraction from liquid energy constraint
      VFeq=(VF*Peq+Q(3))/(d1*Peq+d0)
      ! Update internal energies via p*dV work exchange
      Q(3)=Q(3)-Peq*(VFeq-VF)
      Q(4)=Q(4)+Peq*(VFeq-VF)
      VF=VFeq
   end subroutine P_relax_implicit

   !> Compute viscosity: Sutherland for gas, VF-weighted blend with liquid
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pTG,pVF,pQ,pVisc,pBeta,pDiff
      real(WP) :: r_cyl,blend,mu_g,mu_l,k_g,k_l
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pTG=>fs%TG%mf(lvl)%dataptr(mfi)
            pVF=>fs%VF%mf(lvl)%dataptr(mfi)
            pQ=>fs%Q%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta=>fs%beta%mf(lvl)%dataptr(mfi)
            pDiff=>fs%diff%mf(lvl)%dataptr(mfi)
            ! Get tilebox with overlap
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Gas viscosity from Sutherland
               mu_g=(1.0_WP+Suth_T)*pTG(i,j,k,1)**Suth_n/(Reynolds*(pTG(i,j,k,1)+Suth_T))
               ! Liquid viscosity from ratio
               mu_l=visc_ratio*Reynolds**(-1.0_WP)
               ! Mixture viscosity
               !pVisc(i,j,k,1)=pVF(i,j,k,1)*mu_l+(1.0_WP-pVF(i,j,k,1))*mu_g
               pVisc(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/mu_l+(1.0_WP-pVF(i,j,k,1))/mu_g)
               ! Zero bulk viscosity
               pBeta(i,j,k,1)=0.0_WP
               ! Gas heat diffusivity: k=Cv*Gamma*mu/Pr
               k_g=GammaG*CvG*mu_g/Prandtl
               ! Liquid heat diffusivity from ratio
               k_l=diff_ratio*GammaG*CvG/(Reynolds*Prandtl)
               ! Mixture diffusivity
               !pDiff(i,j,k,1)=pVF(i,j,k,1)*k_l+(1.0_WP-pVF(i,j,k,1))*k_g
               pDiff(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/k_l+(1.0_WP-pVF(i,j,k,1))/k_g)
               ! Apply sponge layer viscosity
               r_cyl=sqrt((amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl))**2+(amr%zlo+(real(k,WP)+0.5_WP)*amr%dz(lvl))**2)
               if (r_cyl.gt.R_spg) then
                  blend=min((r_cyl-R_spg)/L_spg,1.0_WP)**2
                  pVisc(i,j,k,1)=max(pVisc(i,j,k,1),blend*nu_spg*sum(pQ(i,j,k,1:2)))
                  pDiff(i,j,k,1)=max(pDiff(i,j,k,1),blend*nu_spg*sum(pQ(i,j,k,1:2)))
               end if
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities
   
   !> User init callback - set Q and VF/barycenters for a drop at rest with a shock
   subroutine shockdrop_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      use mms_geom, only: initialize_volume_moments
      use amrmpcomp_class, only: VFlo
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pCliq,pCgas
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz,myVF,IEL,x_cc,rhoG,pG,uG,H
      integer :: i,j,k
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      IEL=get_IL(rhoL1,pL1)
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         ! Get pointers to data
         pQ   =>solver%Q%mf(lvl)%dataptr(mfi)
         pVF  =>solver%VF%mf(lvl)%dataptr(mfi)
         pCliq=>solver%Cliq%mf(lvl)%dataptr(mfi)
         pCgas=>solver%Cgas%mf(lvl)%dataptr(mfi)
         ! Loop over grown tilebox
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Compute VF and barycenters from levelset
            call initialize_volume_moments(lo=[solver%amr%xlo+real(i  ,WP)*dx,solver%amr%ylo+real(j  ,WP)*dy,solver%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[solver%amr%xlo+real(i+1,WP)*dx,solver%amr%ylo+real(j+1,WP)*dy,solver%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=sphere_levelset,time=time,level=4,VFlo=VFlo,VF=myVF,BL=BL,BG=BG)
            ! Store volume moments
            pVF(i,j,k,1)=myVF
            pCliq(i,j,k,1:3)=BL
            pCgas(i,j,k,1:3)=BG
            ! Compute local gas state from shock profile
            x_cc=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
            H=Hshock(x=Xs-x_cc,delta=0.5_WP*dx)
            rhoG=rhoG1+(rhoG2-rhoG1)*H
            pG  =pG1  +(pG2  -pG1  )*H
            uG  =u1   +(u2   -u1   )*H
            ! Set conserved variables: Q=(VF*rhoL, (1-VF)*rhoG, VF*rhoL*IL, (1-VF)*rhoG*IG, rho_mix*U, 0, 0)
            pQ(i,j,k,1)=(       myVF)*rhoL1
            pQ(i,j,k,2)=(1.0_WP-myVF)*rhoG
            pQ(i,j,k,3)=pQ(i,j,k,1)*IEL
            pQ(i,j,k,4)=pQ(i,j,k,2)*get_IG(rhoG,pG)
            pQ(i,j,k,5)=(pQ(i,j,k,1)+pQ(i,j,k,2))*uG
            pQ(i,j,k,6)=0.0_WP
            pQ(i,j,k,7)=0.0_WP
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine shockdrop_init

   !> Apply inflow BC at low-x (face=1)
   subroutine shock_dirichlet(solver,pQ,bc_bx,face,time)
      use amrex_amr_module, only: amrex_box
      class(amrmpcomp), intent(inout) :: solver
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      type(amrex_box), intent(in) :: bc_bx
      integer, intent(in) :: face
      real(WP), intent(in) :: time
      integer :: i,j,k
      select case (face)
       case (1)  ! X-LOW: Dirichlet inflow with post-shock (gas only, no liquid)
         do k=bc_bx%lo(3),bc_bx%hi(3); do j=bc_bx%lo(2),bc_bx%hi(2); do i=bc_bx%lo(1),bc_bx%hi(1)
            pQ(i,j,k,1)=0.0_WP                  ! No liquid
            pQ(i,j,k,2)=rhoG2                   ! Gas density
            pQ(i,j,k,3)=0.0_WP                  ! No liquid energy
            pQ(i,j,k,4)=rhoG2*get_IG(rhoG2,pG2) ! Gas internal energy
            pQ(i,j,k,5)=rhoG2*u2                ! X-momentum
            pQ(i,j,k,6)=0.0_WP
            pQ(i,j,k,7)=0.0_WP
         end do; end do; end do
      end select
   end subroutine shock_dirichlet

   !> Tagger based on vorticity, divergence, and VF gradient
   subroutine my_tagger(solver,lvl,tags_ptr,time)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags_ptr
      real(WP), intent(in) :: time
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVisc
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
      real(WP) :: vort_mag,div_neg,rho,mu,Rec,Res
      integer :: i,j,k
      dx=solver%amr%dx(lvl); dxi=1.0_WP/dx
      dy=solver%amr%dy(lvl); dyi=1.0_WP/dy
      dz=solver%amr%dz(lvl); dzi=1.0_WP/dz
      tags=tags_ptr
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         bx=mfi%tilebox()
         tagarr=>tags%dataPtr(mfi)
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         pVisc=>solver%visc%mf(lvl)%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Get rho and mu
            rho=sum(pQ(i,j,k,1:2))
            mu=pVisc(i,j,k,1)
            if (mu.le.0.0_WP) mu=1.0_WP/Reynolds
            ! Get velocity gradient (Q(5:7) are momentum components)
            dudx=0.5_WP*dxi*((pQ(i+1,j,k,5)/max(sum(pQ(i+1,j,k,1:2)),solver%rho_floor))-(pQ(i-1,j,k,5)/max(sum(pQ(i-1,j,k,1:2)),solver%rho_floor)))
            dvdx=0.5_WP*dxi*((pQ(i+1,j,k,6)/max(sum(pQ(i+1,j,k,1:2)),solver%rho_floor))-(pQ(i-1,j,k,6)/max(sum(pQ(i-1,j,k,1:2)),solver%rho_floor)))
            dwdx=0.5_WP*dxi*((pQ(i+1,j,k,7)/max(sum(pQ(i+1,j,k,1:2)),solver%rho_floor))-(pQ(i-1,j,k,7)/max(sum(pQ(i-1,j,k,1:2)),solver%rho_floor)))
            dudy=0.5_WP*dyi*((pQ(i,j+1,k,5)/max(sum(pQ(i,j+1,k,1:2)),solver%rho_floor))-(pQ(i,j-1,k,5)/max(sum(pQ(i,j-1,k,1:2)),solver%rho_floor)))
            dvdy=0.5_WP*dyi*((pQ(i,j+1,k,6)/max(sum(pQ(i,j+1,k,1:2)),solver%rho_floor))-(pQ(i,j-1,k,6)/max(sum(pQ(i,j-1,k,1:2)),solver%rho_floor)))
            dwdy=0.5_WP*dyi*((pQ(i,j+1,k,7)/max(sum(pQ(i,j+1,k,1:2)),solver%rho_floor))-(pQ(i,j-1,k,7)/max(sum(pQ(i,j-1,k,1:2)),solver%rho_floor)))
            dudz=0.5_WP*dzi*((pQ(i,j,k+1,5)/max(sum(pQ(i,j,k+1,1:2)),solver%rho_floor))-(pQ(i,j,k-1,5)/max(sum(pQ(i,j,k-1,1:2)),solver%rho_floor)))
            dvdz=0.5_WP*dzi*((pQ(i,j,k+1,6)/max(sum(pQ(i,j,k+1,1:2)),solver%rho_floor))-(pQ(i,j,k-1,6)/max(sum(pQ(i,j,k-1,1:2)),solver%rho_floor)))
            dwdz=0.5_WP*dzi*((pQ(i,j,k+1,7)/max(sum(pQ(i,j,k+1,1:2)),solver%rho_floor))-(pQ(i,j,k-1,7)/max(sum(pQ(i,j,k-1,1:2)),solver%rho_floor)))
            ! Get vorticity magnitude
            vort_mag=sqrt((dwdy-dvdz)**2+(dudz-dwdx)**2+(dvdx-dudy)**2)
            ! Get dilatation
            div_neg=min(dudx+dvdy+dwdz,0.0_WP)
            ! Tag based on cell Reynolds numbers
            Rec=rho*vort_mag*min(dx,dy,dz)**2/mu
            if (Rec.gt.Rec_tag) tagarr(i,j,k,1)=SETtag
            ! Also tag based on cell shock Reynolds number
            Res=rho*abs(div_neg)*min(dx,dy,dz)**2/mu
            if (Res.gt.Res_tag) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Read EoS and flow parameters
      init_eos_and_flow: block
         use messager, only: log,die
         use string,   only: str_long
         character(len=str_long) :: message
         real(WP) :: A,B,C,c1,cL0
         ! Gas EoS parameters (ideal gas = stiffened gas with Pinf=0)
         call param_read('GammaG',GammaG)
         PinfG=0.0_WP
         ! Liquid EoS: gamma only, PinfL is computed below
         call param_read('GammaL',GammaL)
         ! Shock parameters (gas phase, uses GammaG)
         call param_read('Mach number',M2)
         call param_read('Shock location',Xs)
         ! Post-shock normalization: rhoG2=1, Deltau=1, T2=1
         rhoG2=1.0_WP
         pG2=1.0_WP/(GammaG*M2**2)
         ! Quadratic for rhoG1: A*rhoG1^2 - B*rhoG1 + C = 0
         A=2.0_WP*GammaG*pG2+(GammaG-1.0_WP)
         B=4.0_WP*GammaG*pG2+(GammaG+1.0_WP)
         C=2.0_WP*GammaG*pG2
         rhoG1=(B-sqrt(B**2-4.0_WP*A*C))/(2.0_WP*A)  ! smaller root for compression
         ! Shock-fixed frame velocities and pressure
         u1=1.0_WP/(1.0_WP-rhoG1)
         u2=u1-1.0_WP
         pG1=pG2-rhoG1/(1.0_WP-rhoG1)
         if (pG1.le.0.0_WP) call die('[simulation_init] Cannot achieve requested Mach number - negative pre-shock pressure')
         ! Shock Mach number
         Ms=u1/sqrt(GammaG*pG1/rhoG1)
         ! Shift to lab frame: pre-shock stationary
         u2=1.0_WP
         u1=0.0_WP
         ! CvG from T2=1
         CvG=pG2/(rhoG2*(GammaG-1.0_WP))
         ! Liquid state from density and sound-speed ratios (pre-shock, AIAA nondim)
         call param_read('Density ratio',density_ratio)
         call param_read('Sound speed ratio',sound_ratio)
         c1=sqrt(GammaG*pG1/rhoG1)
         cL0=sound_ratio*c1
         rhoL1=density_ratio*rhoG1
         PinfL=pG1*(density_ratio*GammaG/GammaL*sound_ratio**2-1.0_WP)
         pL1=pG1                                                        ! Force pressure equilibrium
         CvL=(pL1+PinfL)/(rhoL1*(GammaL-1.0_WP))
         ! Viscous parameters
         call param_read('Reynolds number',Reynolds)
         call param_read('Prandtl number',Prandtl)
         call param_read('Viscosity ratio',visc_ratio)
         call param_read('Diffusivity ratio',diff_ratio)
         call param_read('Sutherland exponent',Suth_n)
         call param_read('Sutherland temperature',Suth_T)
         ! Log
         write(message,'("[Post-shock Mach] M2=",es12.5)') M2; call log(message)
         write(message,'("[Shock Mach]      Ms=",es12.5)') Ms; call log(message)
         write(message,'("[Pre-shock]  rhoG1=",es12.5," pG1=",es12.5)') rhoG1,pG1; call log(message)
         write(message,'("[Post-shock] rhoG2=",es12.5," pG2=",es12.5)') rhoG2,pG2; call log(message)
         write(message,'("[Liquid] rhoL1=",es12.5," pL1=",es12.5," cL0=",es12.5)') rhoL1,pL1,cL0; call log(message)
         write(message,'("[Liquid] GammaL=",es12.5," PinfL=",es12.5," CvL=",es12.5)') GammaL,PinfL,CvL; call log(message)
         write(message,'("[Gas]    GammaG=",es12.5," PinfG=",es12.5," CvG=",es12.5)') GammaG,PinfG,CvG; call log(message)
         write(message,'("[Visc]   Re=",es12.5," mu*=",es12.5," Suth_n=",es12.5," Suth_T=",es12.5)') Reynolds,visc_ratio,Suth_n,Suth_T; call log(message)
      end block init_eos_and_flow
      
      ! Initialize AMR grid
      create_amrgrid: block
         amr%name='amrcomp_drop'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=-05.0_WP; amr%xhi=+15.0_WP
         amr%ylo=-10.0_WP; amr%yhi=+10.0_WP
         amr%zlo=-10.0_WP; amr%zhi=+10.0_WP
         amr%xper=.false.; amr%yper=.true.; amr%zper=.true.
         call param_read('Max levels',amr%maxlvl)
         call amr%initialize()
      end block create_amrgrid
      
      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt=time%dtmax
      end block initialize_timetracker

      ! Initialize compressible multiphase solver
      create_solver: block
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap
         use amrmpcomp_class,  only: BC_GAS
         ! Create flow solver
         call fs%initialize(amr=amr,name='drop')
         ! Provide thermodynamic model (6 EOS pointers)
         fs%getPL=>get_PL; fs%getCL=>get_CL; fs%getTL=>get_TL
         fs%getPG=>get_PG; fs%getCG=>get_CG; fs%getTG=>get_TG
         ! Provide pressure relaxation model
         fs%relax=>P_relax_implicit
         ! Set initial conditions and BCs
         fs%vof_lo_bc(1)=BC_GAS
         fs%Q%lo_bc(1,:)=amrex_bc_ext_dir
         fs%Q%hi_bc(1,:)=amrex_bc_foextrap
         fs%user_init=>shockdrop_init
         fs%user_bc=>shock_dirichlet
      end block create_solver
      
      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: amrex_interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=7,ng=0,interp=amrex_interp_none); call dQdt%register()
         call rho_mix%initialize(amr,name='rho_mix',ncomp=1,ng=fs%nover,interp=amrex_interp_none); call rho_mix%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=amrex_interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1,ng=0,interp=amrex_interp_none); call Mach%register()
      end block create_workspace
      
      ! Initialize regridding
      init_regridding: block
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_tagging=>my_tagger
         call param_read('Tagging Rec',Rec_tag)
         call param_read('Tagging Res',Res_tag)
         ! Create initial grid
         call amr%init_from_scratch(time=time%t)
         ! Build PLIC and reset moments
         call fs%build_plic(time%t)
         call fs%reset_moments()
         ! Compute viscosities
         call get_viscosities()
         ! Add artificial bulk viscosity
         call fs%get_viscartif(dt=time%dt,beta=fs%beta)
         call rho_mix%copy(src=fs%Q,srccomp=1,ncomp=1); call rho_mix%add(src=fs%Q,srccomp=2,ncomp=1); call fs%beta%multiply(src=rho_mix,srccomp=1)
         ! Compute Umag and Mach number
         call Umag%get_magnitude(fs%U,fs%V,fs%W)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
      end block init_regridding
      
      ! Initialize visualization
      create_viz: block
         ! Create visualization object
         call viz%initialize(amr,'drop')
         call viz%add_scalar(fs%VF,1,'VF')
         call viz%add_scalar(fs%RHOL,1,'RHOL')
         call viz%add_scalar(fs%RHOG,1,'RHOG')
         call viz%add_scalar(fs%PL,1,'PL')
         call viz%add_scalar(fs%PG,1,'PG')
         call viz%add_scalar(fs%U,1,'U')
         call viz%add_scalar(fs%V,1,'V')
         call viz%add_scalar(fs%W,1,'W')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         call viz%add_surfmesh(fs%smesh,'plic')
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) call viz%write(time=time%t)
      end block create_viz
      
      ! Create monitors
      create_monitors: block
         ! Get solver info and cfl
         call fs%get_info()
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         ! Create simulation monitor
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Qmin(1),'rhoLmin')
         call mfile%add_column(fs%Qmax(1),'rhoLmax')
         call mfile%add_column(fs%PLmin,'PLmin')
         call mfile%add_column(fs%PLmax,'PLmax')
         call mfile%add_column(fs%Qmin(2),'rhoGmin')
         call mfile%add_column(fs%Qmax(2),'rhoGmax')
         call mfile%add_column(fs%PGmin,'PGmin')
         call mfile%add_column(fs%PGmax,'PGmax')
         call mfile%add_column(fs%VFmin,'VFmin')
         call mfile%add_column(fs%VFmax,'VFmax')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLc_x,'CFLc_x')
         call cflfile%add_column(fs%CFLc_y,'CFLc_y')
         call cflfile%add_column(fs%CFLc_z,'CFLc_z')
         call cflfile%add_column(fs%CFLa_x,'CFLa_x')
         call cflfile%add_column(fs%CFLa_y,'CFLa_y')
         call cflfile%add_column(fs%CFLa_z,'CFLa_z')
         call cflfile%add_column(fs%CFLv_x,'CFLv_x')
         call cflfile%add_column(fs%CFLv_y,'CFLv_y')
         call cflfile%add_column(fs%CFLv_z,'CFLv_z')
         call cflfile%write()
         ! Create conservation monitor
         consfile=monitor(amRoot=amr%amRoot,name='conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(fs%Qint(1),'Liquid Mass')
         call consfile%add_column(fs%Qint(2),'Gas Mass')
         call consfile%add_column(fs%Qint(3),'Liquid IntEnergy')
         call consfile%add_column(fs%Qint(4),'Gas IntEnergy')
         call consfile%add_column(fs%Qint(5),'U Momentum')
         call consfile%add_column(fs%Qint(6),'V Momentum')
         call consfile%add_column(fs%Qint(7),'W Momentum')
         call consfile%add_column(fs%rhoKint,'Kinetic energy')
         call consfile%write()
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
      end block create_monitors

   end subroutine simulation_init
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old state
         call fs%Qold%copy(src=fs%Q)
         call fs%VFold%copy(src=fs%VF)
         call fs%Cliqold%copy(src=fs%Cliq)
         call fs%Cgasold%copy(src=fs%Cgas)
         call fs%PLICold%copy(src=fs%PLIC)
         
         ! ===== RK2 Stage 1: dQdt = f(t, Q) =====
         call fs%get_dQdt(Q=fs%Q,dQdt=dQdt,dt=0.5_WP*time%dt,time=time%t)
         
         ! ===== RK2 Stage 2: Q* = Qold + dt/2*dQdt, dQdt* = f(t+dt/2, Q*) =====
         call fs%Q%copy(src=fs%Qold); call fs%Q%saxpy(a=0.5_WP*time%dt,src=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t+0.5_WP*time%dt)
         call fs%apply_relax()
         call fs%get_dQdt(Q=fs%Q,dQdt=dQdt,dt=time%dt,time=time%t+0.5_WP*time%dt)

         ! ===== RK2 Final: Q = Qold + dt*dQdt* =====
         call fs%Q%copy(src=fs%Qold)
         call fs%Q%saxpy(a=time%dt,src=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call fs%apply_relax()

         ! Rebuild PLIC and reset moments
         call fs%build_plic(time%t); call fs%reset_moments()
         
         ! Recompute primitive variables
         call fs%get_primitive(fs%Q)
         
         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Compute viscosities
         call get_viscosities()

         ! Add artificial bulk viscosity
         call fs%get_viscartif(dt=time%dt,beta=fs%beta)
         call rho_mix%copy(src=fs%Q,srccomp=1,ncomp=1); call rho_mix%add(src=fs%Q,srccomp=2,ncomp=1); call fs%beta%multiply(src=rho_mix,srccomp=1)

         ! Compute Umag and Mach number
         call Umag%get_magnitude(fs%U,fs%V,fs%W)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time%t)

         ! Perform and output monitoring
         call fs%get_info()
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         
      end do
      
   end subroutine simulation_run
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Finalize time
      call time%finalize()
      ! Finalize grid
      call amr%finalize()
      call regrid_evt%finalize()
      ! Finalize solver
      call fs%finalize()
      call dQdt%finalize()
      call rho_mix%finalize()
      call Umag%finalize()
      call Mach%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call gridfile%finalize()
   end subroutine simulation_final
   
end module simulation
