!! Rayleigh scattering by the gas, for the shortwave.
!!
!! The correlated-k tables are absorption only, so Rayleigh scattering is added
!! separately in the shortwave. Polarizability form:
!!
!! ```
!!   sigma_s(lambda) = (128 pi^5 / 3) * alpha_s^2 / lambda^4 * F_s
!! ```
!!
!! with ``alpha`` the polarizability volume [cm^3] and ``F`` the King
!! correction. A mixture is the number-fraction weighted sum over species.
!!
!! Holding ``alpha`` constant with wavelength gives, against the Dalgarno &
!! Williams (1962) H2 series:
!!
!! ```
!!   0.3 um   -14%      1.0 um   +1.6%
!!   0.4 um   -6.7%     2.0 um   +2.8%
!!   0.5 um   -3.2%
!!   0.7 um   -0.1%
!! ```
!!
!! The Sonora bands start at 0.268 um, so the largest errors are inside the
!! solved range; a result sensitive to the near-UV albedo needs a dispersing
!! ``alpha``.
!!
!! @author Wolf Cukier
module carma_rayleigh

  use carma_precision_mod

  implicit none

  private
  public :: RAY_NSPEC, RAY_H2, RAY_HE
  public :: RAY_SOLAR_H2HE
  public :: ray_sigma, ray_sigma_mix

  !! Number of species the cross-section table knows about.
  integer, parameter :: RAY_NSPEC = 2

  !! Indices into the species table.
  integer, parameter :: RAY_H2 = 1
  integer, parameter :: RAY_HE = 2

  !! Polarizability volume [cm^3].
  real(kind=f), parameter :: RAY_ALPHA(RAY_NSPEC) = &
      (/ 0.802e-24_f, &      ! H2
         0.205e-24_f /)      ! He

  !! King correction factor [-]. Exactly 1 for He; within a percent of 1 for H2.
  real(kind=f), parameter :: RAY_KING(RAY_NSPEC) = &
      (/ 1.0_f, &            ! H2
         1.0_f /)            ! He

  !! Solar-composition H2/He number fractions.
  real(kind=f), parameter :: RAY_SOLAR_H2HE(RAY_NSPEC) = &
      (/ 0.86_f, &           ! H2
         0.14_f /)           ! He

  real(kind=f), parameter :: PI_F = 3.14159265358979323846_f

  !! 128 pi^5 / 3, the prefactor of the polarizability form.
  real(kind=f), parameter :: RAY_PREFAC = 128._f * PI_F**5 / 3._f

contains

  !! Rayleigh cross-section of one species [cm^2 / molecule].
  !!
  !! Returns zero for a non-positive wavelength or unknown species.
  pure function ray_sigma(ispec, wavelen_cm) result(sigma)

    implicit none

    integer, intent(in)      :: ispec        !! RAY_H2, RAY_HE, ...
    real(kind=f), intent(in) :: wavelen_cm   !! [cm]
    real(kind=f)             :: sigma        !! [cm^2/molecule]

    if (wavelen_cm <= 0._f .or. ispec < 1 .or. ispec > RAY_NSPEC) then
      sigma = 0._f
      return
    end if

    sigma = RAY_PREFAC * RAY_ALPHA(ispec)**2 * RAY_KING(ispec) &
            / wavelen_cm**4

    return
  end function ray_sigma


  !! Number-fraction weighted cross-section of a mixture [cm^2 / molecule].
  !!
  !! ``x`` is used as given, not normalised.
  pure function ray_sigma_mix(x, wavelen_cm) result(sigma)

    implicit none

    real(kind=f), intent(in) :: x(RAY_NSPEC)  !! number fractions [-]
    real(kind=f), intent(in) :: wavelen_cm    !! [cm]
    real(kind=f)             :: sigma         !! [cm^2/molecule]

    integer :: is

    sigma = 0._f
    do is = 1, RAY_NSPEC
      sigma = sigma + x(is) * ray_sigma(is, wavelen_cm)
    end do

    return
  end function ray_sigma_mix

end module carma_rayleigh
