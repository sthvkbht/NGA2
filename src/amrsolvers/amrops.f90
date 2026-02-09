!> AMR discrete operators module
!> Provides stateless utility subroutines for stencil-based operations
!> on AMReX MultiFabs. Inspired by AMReX-Hydro's HydroUtils pattern.
module amrops
   use precision, only: WP
   use amrex_amr_module, only: amrex_multifab, amrex_mfiter, amrex_mfiter_build, &
      amrex_mfiter_destroy, amrex_box
   implicit none
   private

   public :: div_face_to_cell
   public :: grad_cell_to_face

contains

   !> Compute cell-centered divergence from face-centered velocities
   !> div(i,j,k) = (U(i+1)-U(i))*dxi + (V(j+1)-V(j))*dyi + (W(k+1)-W(k))*dzi
   !> @param dxi Inverse mesh spacing in x
   !> @param dyi Inverse mesh spacing in y
   !> @param dzi Inverse mesh spacing in z
   !> @param U_mf Face-centered velocity in x (on x-faces)
   !> @param V_mf Face-centered velocity in y (on y-faces)
   !> @param W_mf Face-centered velocity in z (on z-faces)
   !> @param div_mf Cell-centered divergence (output)
   !> @param tiling Optional: enable tiling (default: false)
   subroutine div_face_to_cell(dxi, dyi, dzi, U_mf, V_mf, W_mf, div_mf, tiling)
      real(WP), intent(in) :: dxi, dyi, dzi
      type(amrex_multifab), intent(in) :: U_mf, V_mf, W_mf
      type(amrex_multifab), intent(inout) :: div_mf
      logical, intent(in), optional :: tiling
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU, pV, pW, pDiv
      logical :: do_tiling
      integer :: i, j, k

      do_tiling = .false.; if (present(tiling)) do_tiling = tiling

      call amrex_mfiter_build(mfi, div_mf, tiling=do_tiling)
      do while (mfi%next())
         bx = mfi%tilebox()
         pU   => U_mf%dataptr(mfi)
         pV   => V_mf%dataptr(mfi)
         pW   => W_mf%dataptr(mfi)
         pDiv => div_mf%dataptr(mfi)

         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)
                  pDiv(i,j,k,1) = (pU(i+1,j,k,1) - pU(i,j,k,1)) * dxi &
                     + (pV(i,j+1,k,1) - pV(i,j,k,1)) * dyi &
                     + (pW(i,j,k+1,1) - pW(i,j,k,1)) * dzi
               end do
            end do
         end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine div_face_to_cell

   !> Compute face-centered gradient from cell-centered pressure
   !> gradx(i,j,k) = (P(i)-P(i-1))*dxi  on x-face
   !> grady(i,j,k) = (P(j)-P(j-1))*dyi  on y-face
   !> gradz(i,j,k) = (P(k)-P(k-1))*dzi  on z-face
   !> @param dxi Inverse mesh spacing in x
   !> @param dyi Inverse mesh spacing in y
   !> @param dzi Inverse mesh spacing in z
   !> @param P_mf Cell-centered pressure
   !> @param gradx_mf Gradient on x-faces (output)
   !> @param grady_mf Gradient on y-faces (output)
   !> @param gradz_mf Gradient on z-faces (output)
   !> @param tiling Optional: enable tiling (default: false)
   subroutine grad_cell_to_face(dxi, dyi, dzi, P_mf, gradx_mf, grady_mf, gradz_mf, tiling)
      real(WP), intent(in) :: dxi, dyi, dzi
      type(amrex_multifab), intent(in) :: P_mf
      type(amrex_multifab), intent(inout) :: gradx_mf, grady_mf, gradz_mf
      logical, intent(in), optional :: tiling
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pP, pGx, pGy, pGz
      logical :: do_tiling
      integer :: i, j, k

      do_tiling = .false.; if (present(tiling)) do_tiling = tiling

      call amrex_mfiter_build(mfi, gradx_mf, tiling=do_tiling)
      do while (mfi%next())
         bx = mfi%tilebox()
         pP  => P_mf%dataptr(mfi)
         pGx => gradx_mf%dataptr(mfi)
         pGy => grady_mf%dataptr(mfi)
         pGz => gradz_mf%dataptr(mfi)

         ! X-gradient on x-faces
         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)+1
                  pGx(i,j,k,1) = (pP(i,j,k,1) - pP(i-1,j,k,1)) * dxi
               end do
            end do
         end do

         ! Y-gradient on y-faces
         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)+1
               do i = bx%lo(1), bx%hi(1)
                  pGy(i,j,k,1) = (pP(i,j,k,1) - pP(i,j-1,k,1)) * dyi
               end do
            end do
         end do

         ! Z-gradient on z-faces
         do k = bx%lo(3), bx%hi(3)+1
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)
                  pGz(i,j,k,1) = (pP(i,j,k,1) - pP(i,j,k-1,1)) * dzi
               end do
            end do
         end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine grad_cell_to_face

end module amrops
