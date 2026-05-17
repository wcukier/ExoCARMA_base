! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine zeroes the fast microphysics sinks and sources, 
!! at one spatial point per call.
!!
!! @author Andy Ackerman
!! @version Oct-1997
subroutine zeromicro(carma, cstate, iz, rc)

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
  integer, intent(in)                  :: iz      !! vertical index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure
  

  ! Set production terms and loss rates due to nucleation, growth,
  ! and evaporation to zero.  Also set index of smallest bin nuceleated
  ! during time step equal to <NBIN> first time through spatial loop.
  
  if (do_grow) then

    phprod(iz)     = 0._f
    rlprod(iz)     = 0._f
    dtpart(iz,:,:) = 0._f

    if (NGAS > 0) gasprod(:,iz) = 0._f

    rhompe(:,:,iz)  = 0._f
    rnucpe(:,:,iz)   = 0._f
    rnucpeup(:,:,iz) = 0._f		!PETER
    growpe(:,:,iz)   = 0._f
    evappe(:,:,iz)   = 0._f
    rnuclg(:,:,:,iz) = 0._f
    rnuclgsum(:,:,iz) = 0._f
    growlg(:,:,iz)   = 0._f
    evaplg(:,:,iz)   = 0._f

    coreavg(iz)   = 0._f
    evdrop(iz)    = 0._f
    coresig(iz)   = 0._f
    evcore(:,iz)  = 0._f
    too_small(iz) = .false.
    too_big(iz)   = .false.
    nuc_small(iz) = .false.

  end if

  ! Return to caller with fast microphysics sinks and sources zeroed.
  return
end
