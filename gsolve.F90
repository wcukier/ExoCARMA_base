! LLM-GENERATED — UNVERIFIED (tier: physics)
! Remove this banner only after completing the checks in Assistant/POLICY.md §2.
! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine calculates new gas concentrations.
!!
!! @author Andy Ackerman, Bill McKie, Chuck Bardeen
!! @version Dec-1995, Sep-1997, Nov-2009
subroutine gsolve(carma, cstate, iz, previous_ice, previous_liquid, previous_rlheat, rc, dtime_in, nretries_in, dt_suggest, &
                  gc_base, gas_check)

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
  real(kind=f), intent(in)             :: previous_ice(NGAS)      !! total ice at the start of substep [g gas equivalent]
  real(kind=f), intent(in)             :: previous_liquid(NGAS)   !! total liquid at the start of substep [g gas equivalent]
  real(kind=f), intent(in)             :: previous_rlheat         !! latent heat of condensate at the start of substep [erg]
  integer, intent(inout)               :: rc      !! return code, negative indicates failure
  real(kind=f), intent(in)             :: dtime_in    !! thread-local dtime for this substep
  real(kind=f), intent(in)             :: nretries_in !! thread-local retry count
  real(kind=f), intent(inout)          :: dt_suggest  !! substep size the failing rates could tolerate [s];
                                                      !! only reduced, and set only on RC_WARNING_RETRY
  real(kind=f), intent(in)             :: gc_base(NGAS) !! gas at the start of the substep [g/x/y/z]; the condensate
                                                        !! change is applied to this, not to the gc the rates saw
  integer, intent(in)                  :: gas_check   !! 0 = no positivity/threshold checks, 1 = checks that request a
                                                      !! retry but never clamp, 2 = checks with the clamp at maxretries

  ! Local Variables
  integer                              :: igas    !! gas index
  integer                              :: jgroup  !! a group condensing this gas (diagnostics)
  integer                              :: ig
  real(kind=f)                         :: gc_cgs
  real(kind=f)                         :: gc_old
  real(kind=f)                         :: rvap
  real(kind=f)                         :: total_ice(NGAS)      ! total ice [g gas equivalent]
  real(kind=f)                         :: total_liquid(NGAS)   ! total liquid [g gas equivalent]
  real(kind=f)                         :: total_rlheat         ! latent heat of condensate [erg]
  real(kind=f)                         :: gc_threshold         ! threshold for changes to gc
  real(kind=f)                         :: net_prod             ! net gas production rate [g/cm3/s]
  real(kind=f), parameter              :: DT_SAFETY = 0.5_f    ! margin on rate-informed substep suggestions
  
  
  1 format(/,'gsolve::ERROR - negtive gas concentration for ',a,' : iz=',i4,',lat=', &
              f7.2,',lon=',f7.2,',gc=',e10.3,',gasprod=',e10.3,',supsati=',e10.3, &
              ',supsatl=',e10.3,',t=',f6.2)
  2 format('gsolve::ERROR - conditions at beginning of the step : gc=',e10.3,',supsati=',e17.10, &
              ',supsatl=',e17.10,',t=',f6.2,',d_gc=',e10.3,',d_t=',f6.2)
  3 format(/,'microfast::WARNING - gas concentration change exceeds threshold: ',a,' : iz=',i4,',lat=', &
              f7.2,',lon=',f7.2, ', (gc-gcl)/gcl=', e10.3)
  4 format(/,'gsolve::WARNING: negtive gas concentration reached for ',a,' : iz=',i4,',lat=', &
              f7.2,',lon=',f7.2,',gc=',e10.3,',gasprod=',e10.3,',supsati=',e10.3, &
              ',supsatl=',e10.3,',t=',f6.2, ', setting gc to 1e-50*rhoa(iz)')
  

  ! Determine the total amount of condensate for each gas, in units of the
  ! gas mass it consumed.
  call totalcondensate(carma, cstate, iz, total_ice, total_liquid, total_rlheat, rc)

  ! Latent heating from the net change in condensate across all groups.
  rlprod(iz) = rlprod(iz) - (previous_rlheat - total_rlheat) / (CP * rhoa(iz) * dtime_in)

  do igas = 1,NGAS

    ! We do not seem to be conserving mass and energy, so rather than relying upon gasprod
    ! and rlheat, recalculate the total change in condensate to determine the change
    ! in gas and energy.
    !
    ! This is because in the old scheme, the particles were solved for implicitly, but the
    ! gas and latent heat were solved for explicitly using the same rates.
    gasprod(igas,iz) = ((previous_ice(igas) - total_ice(igas)) + &
		                  (previous_liquid(igas) - total_liquid(igas))) / dtime_in

    gc_old = gc_base(igas)

    gc(iz,igas) = gc_base(igas) + dtime_in * (gasprod(igas,iz) + phochemprod_gas(iz,igas))

    if (gas_check == 0) cycle

    ! Don't let the gas concentration go negative.
    if (gc(iz,igas) < 0.0_f) then
#ifdef CARMA_TEST_CHECKS
      ! Test-build diagnostic of the failing gas
      write(91,'(a,2i4,f6.0,4es15.6e3,i2,*(es13.4e3))') 'F', iz, igas, nretries_in, dtime_in, gc_old, &
        gasprod(igas,iz), phochemprod_gas(iz,igas), merge(1, 0, nretries_in == maxretries .and. gas_check == 2), &
        (supsatl(iz,ig), ig = 1, NGROUP)
#endif
      if (do_substep) then
        if (nretries_in == maxretries .and. gas_check == 2) then 
          !if (do_print) write(LUNOPRT,1) trim(gasname(igas)), iz, &
          !      lat, lon, gc(iz,igas), gasprod(igas), &
          !      supsati(iz,igas), supsatl(iz,igas), t(iz)
          ! if (do_print) write(LUNOPRT,4) trim(gasname(igas)), iz, &
          !       lat, lon, gc(iz,igas), gasprod(igas), &
          !       supsati(iz,igas), supsatl(iz,igas), t(iz)
          ! if (do_print) write(LUNOPRT,2) gcl(iz,igas), supsatiold(iz,igas), &
          !       supsatlold(iz,igas), told(iz), d_gc(iz,igas), d_t(iz)
          gc(iz,igas) = 1e-50_f * rhoa(iz)
          rc=RC_WARNING
        else
          ! Suggest the substep the depletion rate could tolerate: the time to
          ! consume the gas on hand, with margin.
          net_prod = abs(gasprod(igas,iz) + phochemprod_gas(iz,igas))
          if (net_prod > 0._f .and. gc_old > 0._f) then
            dt_suggest = min(dt_suggest, DT_SAFETY * gc_old / net_prod)
          end if
          rc = RC_WARNING_RETRY
        endif
      else
        ! Report the saturation of the first condensate keyed to this gas.
        jgroup = 0
        do ig = 1, NGROUP
          if (igrowgas(ienconc(ig)) == igas) then
            jgroup = ig
            exit
          end if
        end do
        if (do_print .and. jgroup > 0) write(LUNOPRT,1) trim(gasname(igas)), iz, &
                lat, lon, gc(iz,igas), gasprod(igas,iz), &
                supsati(iz, jgroup), supsatl(iz,jgroup), t(iz)
        rc = RC_WARNING_RETRY
      end if
    end if
    
    ! If gas changes by too much, then retry the calculation.
    gc_threshold = dgc_threshold(igas)
    if (do_incloud) gc_threshold = gc_threshold / cldfrc(iz)
    !
    ! NOTE: If doing incloud calculations, then the threshold needs to be scaled by the
    ! cloud fraction since the cloud mass has been scaled by the cloud fraction.
    
    if (gc_threshold /= 0._f) then
      if ((dtime_in * gasprod(igas,iz) / gc(iz,igas)) > gc_threshold) then
        if (do_substep) then
          if (nretries_in == maxretries) then
            if (do_print) write(LUNOPRT,3) trim(gasname(igas)), iz, &
              lat, lon, dtime_in * gasprod(igas,iz) / gc(iz,igas)
            jgroup = 0
            do ig = 1, NGROUP
              if (igrowgas(ienconc(ig)) == igas) then
                jgroup = ig
                exit
              end if
            end do
            if (do_print .and. jgroup > 0) write(LUNOPRT,2) gcl(iz,igas), supsatiold(iz,jgroup), &
              supsatlold(iz,jgroup), told(iz), d_gc(iz,igas), d_t(iz)
          end if
        else
          if (do_print) write(LUNOPRT,3) trim(gasname(igas)), iz, &
		        lat, lon, dtime_in * gasprod(igas,iz) / gc(iz,igas)
        end if

        ! Suggest the substep that would keep the relative change at the
        ! threshold, with margin.
        if (abs(gasprod(igas,iz)) > 0._f .and. gc(iz,igas) > 0._f) then
          dt_suggest = min(dt_suggest, &
                           DT_SAFETY * gc_threshold * gc(iz,igas) / abs(gasprod(igas,iz)))
        end if

        rc = RC_WARNING_RETRY
      end if
    end if
  end do

  ! Return to caller with new gas concentrations.
  return
end
