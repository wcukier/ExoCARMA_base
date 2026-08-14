!! Rayleigh scattering by the gas, for the shortwave.
!!
!! The correlated-k tables carry **absorption only**, which is the right
!! approximation in the thermal infrared -- molecular scattering there is
!! negligible against absorption and against the clouds. In the shortwave it is
!! not: Rayleigh scattering by H2 rises as ``lambda^-4`` and dominates the
!! clear-sky blue, so leaving it out makes a cloud-free column far too
!! absorbing and biases the albedo the clouds are supposed to set.
!!
!! The cross-section is the standard polarizability form,
!!
!! ```
!!   sigma_s(lambda) = (128 pi^5 / 3) * alpha_s^2 / lambda^4 * F_s
!! ```
!!
!! with ``alpha`` the polarizability volume in cm^3 and ``F`` the King
!! correction for molecular anisotropy. A mixture is the number-fraction
!! weighted sum over species.
!!
!! **Why this form rather than a fitted series.** Dalgarno & Williams (1962)
!! give H2 as a three-term series in ``1/lambda``, which is more accurate to
!! the blue, but it is specific to H2. This form takes one polarizability and
!! one King factor per species, so adding CO2, H2O or N2 for a
!! higher-metallicity atmosphere is a row in the table below rather than a new
!! function.
!!
!! **The accuracy that buys, measured against Dalgarno & Williams for H2:**
!!
!! ```
!!   0.3 um   -14%      1.0 um   +1.6%
!!   0.4 um   -6.7%     2.0 um   +2.8%
!!   0.5 um   -3.2%
!!   0.7 um   -0.1%
!! ```
!!
!! The error comes from holding ``alpha`` constant with wavelength; real
!! polarizabilities disperse upward to the blue. Note the Sonora band set
!! starts at 0.268 um, so the worst of this is inside the range actually
!! solved -- it is not safely off the end of the grid. It is accepted because
!! those bands carry little flux for a cool object and are optically thick in
!! Rayleigh anyway, but a run whose answer depends on the near-UV albedo wants
!! a dispersing ``alpha`` here first.
!!
!! Composition is a **number-fraction vector**, not hard-coded H/He, even
!! though H/He is the only mixture shipped. The physics has no business
!! assuming the atmosphere.
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

  !! King correction factor for molecular anisotropy [-]. Monatomic helium is
  !! isotropic and so exactly 1; H2's is within a percent of it across the
  !! optical, which is inside the error of holding alpha constant anyway.
  real(kind=f), parameter :: RAY_KING(RAY_NSPEC) = &
      (/ 1.0_f, &            ! H2
         1.0_f /)            ! He

  !! Solar-composition number fractions, the only mixture shipped. Phase 8
  !! turns this into an input; until then it is the default a caller gets by
  !! passing it explicitly, rather than something the routines assume.
  real(kind=f), parameter :: RAY_SOLAR_H2HE(RAY_NSPEC) = &
      (/ 0.86_f, &           ! H2
         0.14_f /)           ! He

  real(kind=f), parameter :: PI_F = 3.14159265358979323846_f

  !! 128 pi^5 / 3, the prefactor of the polarizability form.
  real(kind=f), parameter :: RAY_PREFAC = 128._f * PI_F**5 / 3._f

contains

  !! Rayleigh cross-section of one species [cm^2 / molecule].
  !!
  !! ``wavelen_cm`` is a wavelength in cm, matching the units the polarizability
  !! is tabulated in. Returns zero for a non-positive wavelength rather than
  !! dividing by it.
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
  !! ``x`` need not sum to one; it is used as given, so a caller carrying only
  !! the scattering species of a larger mixture gets the answer for those. That
  !! is deliberate -- normalising here would silently reinterpret such a call.
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
