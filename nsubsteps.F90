! LLM-GENERATED — UNVERIFIED (tier: physics)
! Remove this banner only after completing the checks in Assistant/POLICY.md §2.
! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine calculates the number of sub-timesteps <ntsubsteps>
!!  for the current model spatial point.
!!
!! @author Eric Jensen
!! @version Apr-2000
subroutine nsubsteps(carma, cstate, iz, dtime_save, ntsubsteps, rc)

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
  integer, intent(in)                  :: iz          !! z index
  real(kind=f), intent(in)             :: dtime_save  !! original (not substepped) dtime
  integer, intent(inout)               :: ntsubsteps  !! suggested number of substeps
  integer, intent(inout)               :: rc          !! return code, negative indicates failure

  ! Local declarations
  integer                        :: ig      ! group index
  integer                        :: igas    ! gas index
  integer                        :: jg      ! group of the nucleation target element
  integer                        :: ibin    ! bin index
  integer                        :: iepart
  integer                        :: inuc
  integer                        :: ienucto
  integer                        :: ibin_small(NGROUP)
  integer                        :: j       ! depleted-gas list index
  integer                        :: jgas    ! depleted gas
  real(kind=f)                   :: g0
  real(kind=f)                   :: g1
  real(kind=f)                   :: dmdt
  real(kind=f)                   :: dt_adv
  real(kind=f)                   :: dt_exh  ! substep limit from gas exhaustion [s]
  real(kind=f)                   :: demand  ! gas consumption rate of a group [g/cm3/s]
  real(kind=f), parameter        :: DT_SAFETY = 0.5_f  ! margin on the exhaustion estimate
  real(kind=f)                   :: ss
  real(kind=f)                   :: ssold
  real(kind=f)                   :: pvap
  real(kind=f)                   :: vf_max


  ! If substepping is disabled, then use one substep
  if (.not. do_substep) then
    ntsubsteps = 1
  else 
    ! Set default values
    ntsubsteps = minsubsteps
    
    ! Find the bin number of the smallest particle bin that
    ! contains a significant number of particles.
    ! Also check for significant activation of water droplets.
  
    if( ntsubsteps .lt. maxsubsteps )then
  
      do ig = 1, NGROUP
  
        if( pconmax(iz,ig) .gt. FEW_PC) then
  
          ibin_small(ig) = NBIN

          ! element of particle number concentration  
          iepart = ienconc(ig)
          
          if( (itype(iepart) .eq. I_INVOLATILE) .or. &
	    (itype(iepart) .eq. I_VOLATILE) ) then
  
            do inuc = 1,nnuc2elem(iepart)
              ienucto = inuc2elem(inuc,iepart)

              if( (inucproc(iepart,ienucto) .eq. I_DROPACT) .or. &
	  (inucproc(iepart,ienucto) .eq. I_HETGEN) ) then

                ! Supersaturation of the condensate that nucleates onto this
                ! seed: the target element's group.
                jg = igelem(ienucto)
                ss = max( supsatl(iz,jg), supsatlold(iz,jg) )

                do ibin = 1, NBIN
                  if( pc(iz,ibin,iepart) / xmet(iz) / ymet(iz) / zmet(iz) &
		.gt. conmax * pconmax(iz,ig) .and. &
                      ss .gt. scrit(iz,ibin,jg) )then
                    ntsubsteps = maxsubsteps
                  endif
                enddo
              endif
            enddo
  
          elseif( itype(iepart) .eq. I_VOLATILE ) then
          
            do ibin = NBIN-1, 1, -1
              if( pc(iz,ibin,iepart) / xmet(iz) / ymet(iz) / zmet(iz) &
		.gt. conmax * pconmax(iz,ig) )then
                ibin_small(ig) = ibin
              endif
            enddo
  
          endif
        endif
      enddo
    endif
  
    ! Calculate the growth rate of a particle with the mode radius for
    ! each volatile group.  The maximum time-step to use is then the
    ! mass growth rate divided by the mass bin width / 2.
    if( ntsubsteps .lt. maxsubsteps )then
  
      dt_adv = dtime_save
      dt_exh = dtime_save
      do ig = 1, NGROUP
      
        ! element of particle number concentration
        iepart = ienconc(ig)
        
        ! condensing gas
        igas = igrowgas(iepart)

        if (igas /= 0) then
  
          if( pconmax(iz,ig) .gt. FEW_PC ) then
    
            if( itype(iepart) .eq. I_VOLATILE ) then
    
              if( is_grp_ice(ig) )then
                ss = supsati(iz,ig)
                pvap = pvapi(iz,ig)
              else
                ss = supsatl(iz,ig)
                pvap = pvapl(iz,ig)
              endif
    
              g0 = gro(iz,ibin_small(ig),ig)
              g1 = gro1(iz,ibin_small(ig),ig)
              dmdt = abs( pvap * ss * g0 / ( 1._f + g0*g1*pvap ) )
              
              if (dmdt /= 0._f) then
                dt_adv = min( dt_adv, dm(ibin_small(ig),ig)/dmdt )

                ! Gas exhaustion: this group's growth draws on every gas in its
                ! depletion list; don't let one substep consume more than a
                ! fraction of any of them.  Only condensation (ss > 0) can
                ! exhaust a gas -- evaporation returns it.
                if (ss > 0._f) then
                do j = 1, carma%f_ndepgas(ig)
                  jgas = carma%f_idepgas(j, ig)
                  demand = dmdt * pconmax(iz,ig) * carma%f_dgas_dcond(ig, jgas)
                  if (demand > 0._f .and. gc(iz,jgas) > 0._f) then
                    dt_exh = min( dt_exh, DT_SAFETY * gc(iz,jgas) / demand )
                  end if
                end do
                end if
              end if
            endif
          endif
        endif
      enddo

      ntsubsteps = nint(min(real(maxsubsteps, kind=f), real(dtime_save, kind=f) / dt_adv))
      ntsubsteps = max( minsubsteps, ntsubsteps )

      ! The exhaustion criterion obeys the same maxsubsteps bound as the other
      ! predictor criteria: the predictor prescribes on every step, so an
      ! aggressive estimate here taxes the whole run.  A layer that truly
      ! needs more than maxsubsteps fails its first attempt and the
      ! rate-informed retry in newstate jumps straight to the required count
      ! (bounded there by maxntsubsteps_cap).  Real arithmetic so a tiny
      ! dt_exh cannot overflow an integer.
      if (dt_exh < dtime_save) then
        ntsubsteps = max( ntsubsteps, &
             int(min(real(maxsubsteps, kind=f), dtime_save / dt_exh + 1._f)) )
      end if
    endif

    ! If the ice supersaturation is large enough for homogeneous freezing
    ! of sulfate aerosols, then use maximum number of substeps
    if( ntsubsteps .lt. (maxsubsteps) )then
      do ig = 1, NGROUP
      
        ! element of particle number concentration  
        iepart = ienconc(ig)

        do inuc = 1,nnuc2elem(iepart)
          ienucto = inuc2elem(inuc,iepart)

          if (iand(inucproc(iepart,ienucto), I_AERFREEZE) .ne. 0) then
            if( (supsati(iz,igelem(ienucto)) .gt. 0.4_f) .and. (t(iz) .lt. 233.16_f) ) then
              ntsubsteps = maxsubsteps
            endif
          endif
        enddo
      enddo
    endif
  endif

  ! Return to caller with number of sub-timesteps evaluated.
  return
end
