! LLM-GENERATED — UNVERIFIED (tier: physics)
! Remove this banner only after completing the checks in Assistant/POLICY.md §2.
! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine calculates the total amount of condensate associated with
!! each gas, in units of the gas mass consumed to form it, and the total
!! latent heat content of the condensate.
!!
!! @author Chuck Bardeen
!! @version Nov-2009
subroutine totalcondensate(carma, cstate, iz, total_ice, total_liquid, total_rlheat, rc)

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
  real(kind=f), intent(out)            :: total_ice(NGAS)      !! total ice [g gas equivalent / x/y/z]
  real(kind=f), intent(out)            :: total_liquid(NGAS)   !! total liquid [g gas equivalent / x/y/z]
  real(kind=f), intent(out)            :: total_rlheat         !! total latent heat of the condensate [erg / x/y/z]
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer                              :: igroup  ! group index
  integer                              :: icore   ! core index
  integer                              :: igas    ! gas index
  integer                              :: igascore    ! gas index of volatile core
  integer                              :: ibin    ! bin index
  integer                              :: ielem   ! element index
  integer                              :: i
  integer                              :: jgas    ! depleted gas index
  integer                              :: jgroup  ! group of the volatile core
  real(kind=f)                         :: coremass
  real(kind=f)                         :: volatilemass
  real(kind=f)                         :: sto_vol      ! g gas per g condensate for this group
  real(kind=f)                         :: sto_core     ! g gas per g condensate for the volatile core


  ! Initialize local variables for keeping track of gas changes due
  ! to nucleation and growth in each particle group.
  total_ice(:)    = 0._f
  total_liquid(:) = 0._f
  total_rlheat    = 0._f

  ! Iterate over each particle type and total up that ones that interact
  ! with the gases.
  !
  ! This code assumes that all changes in condensate are associated with
  ! growth in a particular gas. This doesn't handle all possible changes
  ! associated with nucleation, if the group do not also participate in
  ! growth.
  do igroup = 1,NGROUP

    ielem = ienconc(igroup)     ! element of particle number concentration

    igas = igrowgas(ielem)      ! condensing gas

    igascore = 0
    icore = 0

    if ((itype(ielem) == I_VOLATILE) .and. (igas /= 0)) then

      ! Mass of key gas consumed per mass of condensate, for the latent heat
      sto_vol = carma%f_dgas_dcond(igroup, igas)

      do ibin = 1, NBIN
   
        ! If this group has core masses, then determine the involatile component.
        coremass = 0._f
      
        do i = 1, ncore(igroup)
          icore = icorelem(i, igroup)
          coremass = coremass + pc(iz, ibin, icore)
        end do

	! Assume ncore = 1 or 0 always
	if (icore > 0) then
	  if ( itype(ievp2elem(icore)) == I_VOLATILE ) then
	    igascore = igrowgas(ievp2elem(icore))
	    jgroup   = igelem(ievp2elem(icore))
	    sto_core = carma%f_dgas_dcond(jgroup, igascore)
	  endif
	endif
        
        volatilemass = (pc(iz, ibin, ielem) * rmass(ibin, igroup)) - coremass
        
        ! There seem to be times when the coremass becomes larger than the total
        ! mass. This shouldn't happen, but check for it here.
        !
        ! NOTE: This can be caused by advection in the parent model or sedimentation
        ! in this model.
        if (volatilemass > 0._f) then
          ! Charge every gas the group depletes; count the latent heat once,
          ! per mass of key gas.
          if (is_grp_ice(igroup)) then
            do jgas = 1, NGAS
              if (carma%f_dgas_dcond(igroup, jgas) == 0._f) cycle
              total_ice(jgas) = total_ice(jgas) + volatilemass * carma%f_dgas_dcond(igroup, jgas)
            end do
            total_rlheat = total_rlheat + volatilemass * sto_vol &
                         * (rlhe(iz,igroup) + rlhm(iz,igroup))
          else
            do jgas = 1, NGAS
              if (carma%f_dgas_dcond(igroup, jgas) == 0._f) cycle
              total_liquid(jgas) = total_liquid(jgas) + volatilemass * carma%f_dgas_dcond(igroup, jgas)
            end do
            total_rlheat = total_rlheat + volatilemass * sto_vol * rlhe(iz,igroup)
          end if
        end if
	if (igascore > 0) then
	  if (is_grp_ice(jgroup)) then
            do jgas = 1, NGAS
              if (carma%f_dgas_dcond(jgroup, jgas) == 0._f) cycle
              total_ice(jgas) = total_ice(jgas) + coremass * carma%f_dgas_dcond(jgroup, jgas)
            end do
            total_rlheat = total_rlheat + coremass * sto_core &
                         * (rlhe(iz,jgroup) + rlhm(iz,jgroup))
          else
            do jgas = 1, NGAS
              if (carma%f_dgas_dcond(jgroup, jgas) == 0._f) cycle
              total_liquid(jgas) = total_liquid(jgas) + coremass * carma%f_dgas_dcond(jgroup, jgas)
            end do
            total_rlheat = total_rlheat + coremass * sto_core * rlhe(iz,jgroup)
          end if
	end if
      end do
    end if
  end do
      
  return
end
