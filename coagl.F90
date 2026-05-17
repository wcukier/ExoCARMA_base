#include "carma_globaer.h"

!! This routine calculates coagulation loss rates <coaglgg>.
!! See [Jacobson, et al., Atmos. Env., 28, 1327, 1994] for details
!! on the coagulation algorithm.
!!
!! The loss rates for all particle elements in a particle group are equal.
!!
!! @author Eric Jensen
!! @version Oct-1995
subroutine coagl(carma, cstate, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_planet_mod
  use carma_condensate_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local Variables
  integer :: ig
  integer :: jg
  integer :: je
  integer :: igrp
  integer :: iz
  integer :: i
  integer :: j
  real(kind=f) :: volx_val  ! hoisted volx(igrp,ig,jg,i,j) — iz-invariant


  ! Loop over particle groups for which coagulation loss is being
  ! calculated.
  do ig = 1,NGROUP

    ! Loop over particle groups that particle in group ig might
    ! collide with.
    do jg = 1,NGROUP
  
      ! Element corresponding to particle number concentration
      je = ienconc(jg)
  
      ! Particle resulting from coagulation between groups <ig> and <jg> goes
      ! to group <igrp>
      igrp = icoag(ig,jg)
  
      ! Resulting particle is in same group as particle under consideration --
      ! partial loss (muliplies <volx>).
      if( igrp .eq. ig )then

        ! Reordered (i, j, iz) so the inner iz loop is contiguous on
        ! ckernel(:,i,j,...), pcl(:,j,je), coaglg(:,i,ig) and the iz-invariant
        ! volx lookup is done once per (i,j) instead of NZ times. The
        ! per-iz pconmax guard is preserved. Floating-point summation order
        ! over j for each (iz,i) is unchanged, so output is bit-identical.
        do i = 1, NBIN-1
          do j = 1, NBIN
            volx_val = volx(igrp,ig,jg,i,j)
            do iz = 1, NZ
              if( pconmax(iz,jg) .gt. FEW_PC .and. &
                  pconmax(iz,ig) .gt. FEW_PC )then
                coaglg(iz,i,ig) = coaglg(iz,i,ig) &
                          + ckernel(iz,i,j,ig,jg) * &
                          pcl(iz,j,je) * volx_val
              endif
            enddo
          enddo
        enddo
  
      !  Resulting particle is in a different group -- complete loss (no <volx>).
      else if( igrp .ne. ig .and. igrp .ne. 0 )then

        ! Reordered (i, j, iz) for the same reason as the partial-loss branch
        ! above. No volx here (different group → complete loss).
        do i = 1, NBIN
          do j = 1, NBIN
            do iz = 1, NZ
              if( pconmax(iz,jg) .gt. FEW_PC .and. &
                  pconmax(iz,ig) .gt. FEW_PC )then
                coaglg(iz,i,ig) = coaglg(iz,i,ig) &
                      + ckernel(iz,i,j,ig,jg) * &
                      pcl(iz,j,je)
              endif
            enddo
          enddo
        enddo
      endif  ! igrp .eq. ig ?
    enddo  ! jg
  enddo  ! ig

  ! Boundary condition: Particles from bin <NBIN> are only lost by
  ! coagulating into other elements. (This is taken care of by <NBIN>-1
  ! limit above)

  !  Return to caller with particle loss rates due to coagulation evaluated.
  return
end
