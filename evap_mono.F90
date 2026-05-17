! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine calculates particle source terms <evappe> due to total
!! evaporation from bin <ibin> group <ig> into a monodisperse
!! distribution.
!!
!! Distinct evaporation of cores has not been treated.
!!
!! @author Andy Ackerman
!! @version Aug-2001
subroutine evap_mono(carma,cstate,iz,ibin,ig,iavg,ieto,igto,rc)

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
  integer, intent(in)                  :: iz      !! z index
  integer, intent(in)                  :: ibin    !! bin index
  integer, intent(in)                  :: ig      !! group index
  integer, intent(in)                  :: iavg
  integer, intent(in)                  :: ieto
  integer, intent(in)                  :: igto
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer                              :: ic  
  integer                              :: iecore
  integer                              :: ie2cn
  integer                              :: jbin
  logical                              :: conserve_mass
  real(kind=f)                         :: factor
  real(kind=f)                         :: fracmass

 ! write(*,*) 'mono', iz, ibin, ig, iavg, ieto, igto

  ! Define option to conserve mass or number when a choice must be made
  ! during monodisperse total evaporation beyond CN grid -- should be done in setupaer()
  conserve_mass = .true.

  ! Set automatic flag for total evaporation used in gasexchange()
  totevap(ibin,ig,iz) = .true.

  ! Possibly put all of core mass into largest, smallest, or
  ! smallest nucelated CN bin 
  if( too_big(iz) .or. too_small(iz) .or. nuc_small(iz) )then

    if( too_big(iz) )then
      jbin = NBIN
    elseif( too_small(iz) )then
      jbin = 1
    else
      jbin = 1
    endif

    if( conserve_mass )then
      factor = coreavg(iz)/rmass(jbin,igto)
    else
      factor = ONE
    endif

    ! First the CN number concentration element
    evappe(jbin,ieto,iz) = evappe(jbin,ieto,iz) + factor*evdrop(iz)

    ! Now the CN cores
    do ic = 2, ncore(ig)
      iecore = icorelem(ic,ig)
      ie2cn  = ievp2elem(iecore)
      evappe(jbin,ie2cn,iz) = evappe(jbin,ie2cn,iz) + &
        factor*evcore(ic,iz)*rmass(jbin,igto)
    enddo
  else

    ! Partition core mass between two CN bins, conserving total core mass
    ! and number.  The number will be subdivided into bins <iavg> and <iavg>-1.
    if( iavg .le. 1 .or. iavg .gt. NBIN )then
      if (do_print) write(LUNOPRT, *) "evap_mono: bad iavg = , ", iavg
      rc = RC_ERROR
      return
    endif

    fracmass = ( rmass(iavg,igto) - coreavg(iz) ) / diffmass(iavg,igto,iavg-1,igto)
!    fracmass = max( 0._f, min( ONE, fracmass ) )

    ! First the CN number concentration element
    evappe(iavg-1,ieto,iz) = evappe(iavg-1,ieto,iz) + evdrop(iz)*fracmass
    evappe(iavg,ieto,iz) = evappe(iavg,ieto,iz) + evdrop(iz)*( ONE - fracmass )

    ! Now the cores
    do ic = 2, ncore(ig)
      iecore = icorelem(ic,ig)
      ie2cn  = ievp2elem(iecore)
      evappe(iavg-1,ie2cn,iz) = evappe(iavg-1,ie2cn,iz) + &
          rmass(iavg-1,igto)*evcore(ic,iz)*fracmass
      evappe(iavg,ie2cn,iz) = evappe(iavg,ie2cn,iz) + &
          rmass(iavg,igto)*evcore(ic,iz)*( ONE - fracmass )
    enddo
  endif

  return
end
