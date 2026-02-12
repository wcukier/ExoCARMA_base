! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine defines time-independent parameters used to calculate
!! condensational growth/evaporation.
!!
!! The parameters defined for each gas are 
!1>
!!   gwtmol:   molecular weight [g/mol]
!!   diffus:   diffusivity      [cm^2/s]
!!   rlhe  :   latent heat of evaporation [cm^2/s^2]
!!   rlhm  :   latent heat of melting [cm^2/s^2]
!!<
!! Time-independent parameters that depend on particle radius are
!! defined in setupgkern.f.
!!
!! This routine requires that vertical profiles of temperature <T>,
!! and pressure <p> are defined.
!!
!! @author Andy Ackerman
!! @version Dec-1995
subroutine setupgrow(carma, cstate, rc)

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
  integer, intent(inout)               :: rc       !! return code, negative indicates failure
  
  ! Local Variable
  integer                        :: ielem    !! element index
  integer                        :: k        !! z index
  integer                        :: i, igas, igroup
  real(kind=f)                   :: rhoa_cgs, aden, wtpctf, wtmol_dif, coldia
  ! Define formats
  1 format(a,':  ',12i6)
  2 format(a,':  ',i6)
  3 format( /' id  gwtmol   gasname',( /,i3,3x,f5.1,3x,a))
  5 format( /,'Particle growth mapping arrays (setupgrow):')


  !-----Check that values are valid------------------------------------------
  do ielem = 1, NELEM
    if( igrowgas(ielem) .gt. NGAS )then
      if (do_print) write(LUNOPRT,*) 'setupgrow::ERROR - component of igrowgas > NGAS'
      rc = -1
      return
    endif
  enddo

  ! Define parameters with weak time-dependence to be used in
  ! growth equation.
  do k = 1, NZ

    ! Diffusivity of water vapor in air from Pruppacher & Klett (eq. 13-3);
    ! units are [cm^2/s].
    do igroup = 1, NGROUP
      ielem = ienconc(igroup)     ! element of particle number concentration
      igas = igrowgas(ielem) 


      if (igas .eq. igash2o) then 
        rhoa_cgs = rhoa(k) / (xmet(k) * ymet(k) * zmet(k))
        ! diffus(k, igash2o) = 0.211_f * (1.01325e+6_f / p(k)) * (t(k) / 273.15_f )**1.94_f
        diffus(k, igroup) = 5._f / (16._f * AVG * COLDIA_H2O**2_f * rhoa_cgs * COLINT) * &
        sqrt(RGAS * t(k) * wtmol_air(k) * (WTMOL_H2O + wtmol_air(k)) / (2._f * PI * WTMOL_H2O)) 
    
        ! Latent heat of evaporation for water; units are [cm^2/s^2]
        if (do_cnst_rlh) then
          rlhe(k, igroup) = RLHE_CNST
        else
          ! from Stull
          rlhe(k, igroup) = (2.5_f - .00239_f * (t(k) - 273.16_f)) * 1.e10_f      
        end if
    
        ! Latent heat of ice melting; units are [cm^2/s^2]
        if (do_cnst_rlh) then
          rlhm(k, igroup) = RLHM_CNST
        else
    
          ! from Pruppacher & Klett (eq. 4-85b)
          !
          ! NOTE: This expression yields negative values for rlmh at mesospheric
          ! temperatures.
          rlhm(k, igroup) = (79.7_f + 0.485_f * (t(k) - 273.16_f) - 2.5e-3_f * &
            ((t(k) - 273.16_f)**2)) * 4.186e7_f
        end if

        ! Properties for H2SO4
      else if (igas .eq. igash2so4) then
        ! Diffusivity
        rhoa_cgs = rhoa(k) / (xmet(k) * ymet(k) * zmet(k))
        aden     = rhoa_cgs * AVG / wtmol_air(k)
        diffus(k,igroup) = 5._f / (16._f * AVG * COLDIA_H2SO4**2_f * rhoa_cgs * COLINT) * &
        sqrt(RGAS * t(k) * wtmol_air(k) * (WTMOL_H2SO4 + wtmol_air(k)) / (2._f * PI * WTMOL_H2SO4)) 
    
        wtpctf = wtpct(k)/100._f
        ! HACK: make H2SO4 latent heats same as water
        !      rlhe(k,igash2so4) = rlhe(k, igash2o)
        !      rlhm(k,igash2so4) = rlhe(k, igash2o)
        ! From Jang et al. 2006, Transactions of the Korean Nuclear Society Autumn Meeting
        rlhe(k,igroup) = 1364.93*wtpctf**3._f - 1226.46*wtpctf**2._f + 382.23*wtpctf + 540.52
        rlhe(k,igroup) = rlhe(k,igroup) * 4.184e7_f
        rlhm(k,igroup) = rlhe(k,igroup)

      else ! WC
        coldia = carma%f_group(igroup)%f_coldia
        wtmol_dif = carma%f_gas(igas)%f_wtmol_dif
        ! Diffusivity

        diffus(k,igroup) = 1._f/carma%f_group(igroup)%f_stofact * 5._f / (16._f * AVG * coldia**2_f * rhoa_cgs * COLINT) * &
          sqrt(RGAS * t(k) * wtmol_air(k) * (wtmol_dif + wtmol_air(k)) / (2._f * PI * wtmol_dif)) 
        
        if (carma%f_group(igroup)%f_lat_heat_e .gt. 0) then
          rlhe(k,igroup) = carma%f_group(igroup)%f_lat_heat_e 
        else 
          rlhe(k,igroup) =  carma%f_group(igroup)%f_vp_tcoeff * log(10._f) * RGAS / wtmol_dif
        endif
        rlhm(k,igroup) = rlhe(k, igroup)
      end if
    enddo

    
  enddo

#ifdef DEBUG
  ! Report some initialization values
  if (do_print_init) then
    write(LUNOPRT,5)
    write(LUNOPRT,2) 'NGAS    ',NGAS
    write(LUNOPRT,1) 'igrowgas',(igrowgas(i),i=1,NELEM)
    ! write(LUNOPRT,3) (i,gwtmol(i),gasname(i),i=1,NGAS)
  endif
#endif

  ! Return to caller with particle growth mapping arrays and time-dependent
  ! parameters initialized.
  return
end
