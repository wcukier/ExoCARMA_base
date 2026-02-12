! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine evaluates supersaturations <supsatl> and <supsati> for all gases.
!!
!! @author Andy Ackerman, Chuck Bardeen
!! @version Dec-1995, Aug-2010
subroutine supersat(carma, cstate, iz, igroup, rc)

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
  integer, intent(in)                  :: igroup    !! group index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  real(kind=f)  :: rvap
  real(kind=f)  :: gc_cgs
  real(kind=f)  :: alpha
  integer                  :: ielem, igas, stofact    !! gas index

  ielem = ienconc(igroup)     ! element of particle number concentration
  igas = igrowgas(ielem) 

  ! Calculate vapor pressures.
  call vaporp(carma, cstate, iz, igroup, rc)

  ! Define gas constant for this gas
  rvap = RGAS / gwtmol_dif(igas)

  gc_cgs = gc(iz,igas) / (zmet(iz)*xmet(iz)*ymet(iz))
  

  ! Add in reaction saturation ratio correction for type III reactions (Helling
  ! and Woitke 2006)
  if ( carma%f_group(igroup)%f_is_type3 .eq. 1 ) then ! WC
    stofact = carma%f_group(igroup)%f_stofact
    ! note this calculates S_r = S^(1/stofact) - 1, the below calculation is for S-1
    supsatl(iz,igroup) = (gc_cgs * rvap * t(iz) / pvapl(iz,igroup))**(1/stofact) - 1._f 
    supsati(iz,igroup) = (gc_cgs * rvap * t(iz) / pvapi(iz,igroup))**(1/stofact) - 1._f
  else
    supsatl(iz,igroup) = (gc_cgs * rvap * t(iz) - pvapl(iz,igroup)) / pvapl(iz,igroup)
    supsati(iz,igroup) = (gc_cgs * rvap * t(iz) - pvapi(iz,igroup)) / pvapi(iz,igroup)
  endif



  ! For subgrid scale clouds, the supersaturation needs to be increased be scaled
  ! based upon cloud fraction. This approach is similar to Wilson and Ballard (1999),
  ! except that only the water vapor (no liquid water) is used to determine the available
  ! water.
  !
  ! NOTE: This assumes that the cloud is an ice cloud.
  if (do_incloud) then
    alpha = rhcrit(iz) * (1._f - cldfrc(iz)) + cldfrc(iz)
    
    supsatl(iz,igroup) = (gc_cgs * rvap * t(iz) - alpha * pvapl(iz,igroup)) / pvapl(iz,igroup)
    supsati(iz,igroup) = (gc_cgs * rvap * t(iz) - alpha * pvapi(iz,igroup)) / pvapi(iz,igroup)
    
    ! Limit supersaturation to liquid saturation.
    supsatl(iz,igroup) = min(supsatl(iz,igroup), 0._f)
    supsati(iz,igroup) = min(supsati(iz,igroup), (pvapl(iz,igroup) &
	- alpha * pvapi(iz,igroup)) / pvapi(iz,igroup))        
  end if


  ! supsatl close to 0 causes carma to crash
  if (abs(supsatl(iz,igroup)) < 1e-16_f) then
    supsatl(iz,igroup) = 1e-16_f
  end if

  return


end


!! This routine evaluates supersaturations <supsatl> and <supsati> for all gases, but
!! thus version of the routine does not scale the supersaturation based on the cloud
!! fraction. It also assumes that vaporp has already been called.
!!
!! @author Andy Ackerman, Chuck Bardeen
!! @version Dec-1995, Aug-2010
subroutine supersat_nocldf(carma, cstate, iz, igroup, ssi, ssl, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(in)                  :: iz      !! z index
  integer, intent(in)                  :: igroup    !! gas index
  real(kind=f), intent(out)            :: ssl
  real(kind=f), intent(out)            :: ssi
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  real(kind=f)  :: rvap
  real(kind=f)  :: gc_cgs
  real(kind=f)  :: alpha
  integer                  :: ielem, igas, stofact    !! gas index

  ielem = ienconc(igroup)     ! element of particle number concentration
  igas = igrowgas(ielem) 

  ! Calculate vapor pressures.
  call vaporp(carma, cstate, iz, igroup, rc)

  ! Define gas constant for this gas
  rvap = RGAS / gwtmol_dif(igas)

  gc_cgs = gc(iz,igas) / (zmet(iz)*xmet(iz)*ymet(iz))

  ssl = (gc_cgs * rvap * t(iz) - pvapl(iz,igroup)) / pvapl(iz,igroup)
  ssi = (gc_cgs * rvap * t(iz) - pvapi(iz,igroup)) / pvapi(iz,igroup)

  return
end
