!! Band-integrated Planck intensity, for the radiatively coupled version of
!! CARMA.
!!
!! This is a direct port of pyHARP's ``bbflux_wavenumber(wn1, wn2, temp)``
!! (``pyharp/src/radiation/bbflux.cpp``), which is itself the classic Wiscombe
!! PLKAVG algorithm shipped with DISORT. It integrates the Planck function over
!! a finite wavenumber interval, which is what a correlated-k band needs: every
!! g-point in a Sonora band carries the blackbody intensity of the *whole* band,
!! and the Gauss weights (which sum to one per band) then average within it.
!!
!! Note this is a different quantity from the existing ``planck.F90``, which
!! gives the monochromatic intensity in cgs. This module works in **SI**
!! (W/m^2/sr), matching pyHARP, so that the port can be validated against it
!! directly. Conversion to cgs happens at the boundary of ``carma_rtsolve``.
!!
!! @author Wolf Cukier
module carma_planck

  use carma_precision_mod

  implicit none

  private
  public :: bbflux_wavenumber, bbflux_wavenumber_col

  ! Constants, matching bbflux.cpp exactly. Any drift here shows up as a
  ! multiplicative error on every flux, so they are deliberately not shared
  ! with carma_constants_mod (whose values are cgs and differ in the last
  ! digits).
  real(kind=f), parameter :: C2     = 1.438786_f     ! h*c/k [cm K]
  real(kind=f), parameter :: SIGMA  = 5.67032e-8_f   ! Stefan-Boltzmann [W/m^2/K^4]
  real(kind=f), parameter :: VCUT   = 1.5_f          ! power- vs exponential-series cutoff
  real(kind=f), parameter :: PI_F   = 3.14159265358979323846_f
  real(kind=f), parameter :: SIGDPI = SIGMA / PI_F
  real(kind=f), parameter :: CONC   = 15._f / PI_F**4
  real(kind=f), parameter :: C1     = 1.1911e-18_f   ! h*c^2 [W/(m^2 sr cm^-4)]

  real(kind=f), parameter :: A1 =  1._f / 3._f
  real(kind=f), parameter :: A2 = -1._f / 8._f
  real(kind=f), parameter :: A3 =  1._f / 60._f
  real(kind=f), parameter :: A4 = -1._f / 5040._f
  real(kind=f), parameter :: A5 =  1._f / 272160._f
  real(kind=f), parameter :: A6 = -1._f / 13305600._f

contains

  !! Planck intensity integrated from wavenumber wn1 to wn2 at temperature
  !! temp, in W/m^2/sr.
  !!
  !! Both series expansions below are in the dimensionless variable
  !! v = C2 * wn / T. For small v the integral from 0 to v is a power series;
  !! for large v it is an exponential series. Which of the two endpoints falls
  !! on which side of VCUT selects the combination.
  function bbflux_wavenumber(wn1, wn2, temp)

    implicit none

    real(kind=f), intent(in) :: wn1   !! lower wavenumber [cm^-1]
    real(kind=f), intent(in) :: wn2   !! upper wavenumber [cm^-1], >= wn1
    real(kind=f), intent(in) :: temp  !! temperature [K], > 0
    real(kind=f)             :: bbflux_wavenumber  !! [W/m^2/sr]

    real(kind=f) :: v(2), p(2), d(2)
    real(kind=f) :: vsq, ex, exm, mv, arg, ans
    integer      :: i, m, smallv

    ! A degenerate interval means monochromatic intensity, not an integral.
    if (wn1 == wn2) then
      arg = exp(-C2 * wn1 / temp)
      bbflux_wavenumber = C1 * wn1**3 * arg / (1._f - arg)
      return
    end if

    v(1) = C2 * wn1 / temp
    v(2) = C2 * wn2 / temp

    smallv = 0
    p(:)   = 0._f
    d(:)   = 0._f

    do i = 1, 2
      if (v(i) < VCUT) then
        smallv = smallv + 1

        vsq  = v(i) * v(i)
        p(i) = CONC * vsq * v(i) * &
               (A1 + v(i) * (A2 + v(i) * (A3 + vsq * (A4 + vsq * (A5 + vsq * A6)))))
      else
        ! Exponential series, truncated at 6 terms as in DISORT.
        ex  = exp(-v(i))
        exm = 1._f

        do m = 1, 6
          mv   = real(m, kind=f) * v(i)
          exm  = exm * ex
          d(i) = d(i) + exm * (6._f + mv * (6._f + mv * (3._f + mv))) &
                        / real(m, kind=f)**4
        end do

        d(i) = d(i) * CONC
      end if
    end do

    if (smallv == 2) then
      ! Both endpoints small: difference of two power series.
      ans = p(2) - p(1)
    else if (smallv == 1) then
      ! Straddles the cutoff: total minus the tail below wn1 and above wn2.
      ans = 1._f - p(1) - d(2)
    else
      ! Both endpoints large: difference of two exponential series.
      ans = d(1) - d(2)
    end if

    bbflux_wavenumber = ans * SIGDPI * temp**4

    return
  end function bbflux_wavenumber


  !! Band-integrated Planck intensity for a whole column of temperatures.
  !!
  !! Convenience wrapper; the level Planck source in the Toon solver needs one
  !! of these per band per level.
  subroutine bbflux_wavenumber_col(wn1, wn2, temp, be)

    implicit none

    real(kind=f), intent(in)  :: wn1        !! lower wavenumber [cm^-1]
    real(kind=f), intent(in)  :: wn2        !! upper wavenumber [cm^-1]
    real(kind=f), intent(in)  :: temp(:)    !! temperatures [K]
    real(kind=f), intent(out) :: be(:)      !! intensity [W/m^2/sr]

    integer :: i

    do i = 1, size(temp)
      be(i) = bbflux_wavenumber(wn1, wn2, temp(i))
    end do

    return
  end subroutine bbflux_wavenumber_col

end module carma_planck
