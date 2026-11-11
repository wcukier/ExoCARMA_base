!! Two-stream shortwave radiative transfer for the radiatively coupled version
!! of CARMA.
!!
!! Toon, McKay & Ackerman (1989), following pyHARP's ``toon_mckay89_shortwave``
!! (MIT licensed). The source is a collimated beam entering the top.
!!
!! Conventions are those of ``carma_rtsolve``: arrays bottom-to-top, SI units.
!!
!!  * Delta-Eddington scaling is always applied (forward-peaked cloud scattering).
!!  * A column with no scattering and no surface reflection is pure
!!    Beer-Lambert; the general coefficients are singular in that limit.
!!
!! @author Wolf Cukier
module carma_swsolve

  use carma_precision_mod
  use carma_rtsolve, only : dtridgl

  implicit none

  private
  public :: toon_sw_column

  !! Cap on the exponent argument, so optically thick layers do not overflow exp().
  real(kind=f), parameter :: EXPMAX = 35._f

  !! Below this single-scattering albedo a layer is treated as purely absorbing.
  real(kind=f), parameter :: WMIN = 1.e-12_f

  !! Substituted for a coefficient that would otherwise be exactly zero and divide.
  real(kind=f), parameter :: TINY_SUB = 1.e-10_f

contains

  !! Solve the shortwave two-stream problem for a single column at a single
  !! spectral point.
  !!
  !! All input and output arrays are bottom-to-top: index 1 is the deepest
  !! layer/level, index nlay (or nlay+1) is the top of the atmosphere.
  !!
  !! ``f0`` is the beam flux through a surface normal to it; the flux entering
  !! the column is ``f0 * mu``.
  subroutine toon_sw_column(nlay, f0, mu_in, dtau_in, w0_in, g_in, w_surf, &
                            f_up, f_dn, f_up_mid, f_dn_mid)

    implicit none

    integer, intent(in)      :: nlay             !! number of layers
    real(kind=f), intent(in) :: f0               !! incident beam flux [W/m^2]
    real(kind=f), intent(in) :: mu_in(nlay+1)    !! cos(incidence) per level [-]
    real(kind=f), intent(in) :: dtau_in(nlay)    !! layer optical thickness
    real(kind=f), intent(in) :: w0_in(nlay)      !! single scattering albedo
    real(kind=f), intent(in) :: g_in(nlay)       !! asymmetry parameter
    real(kind=f), intent(in) :: w_surf           !! surface albedo (0 = gas giant)

    real(kind=f), intent(out) :: f_up(nlay+1)    !! upward flux [W/m^2]
    real(kind=f), intent(out) :: f_dn(nlay+1)    !! downward flux [W/m^2]
    real(kind=f), intent(out), optional :: f_up_mid(nlay) !! layer-midpoint upward
    real(kind=f), intent(out), optional :: f_dn_mid(nlay) !! layer-midpoint downward

    integer :: nlev, l, k, n, i

    ! top-down working copies
    real(kind=f) :: dtau_raw(nlay), w_raw(nlay), g_raw(nlay), mu(nlay+1)
    real(kind=f) :: dtau(nlay), w0(nlay), hg(nlay), mu_zm(nlay)
    real(kind=f) :: tau_in(nlay+1), tau(nlay+1), cum_trans(nlay+1)
    real(kind=f) :: dir(nlay+1)
    real(kind=f) :: gam(nlay), Ap(nlay), Am(nlay)
    real(kind=f) :: Cpm1(nlay), Cmm1(nlay), Cp(nlay), Cm(nlay)
    real(kind=f) :: Ep(nlay), Em(nlay), E1(nlay), E2(nlay), E3(nlay), E4(nlay)

    real(kind=f) :: Af(2*nlay), Bf(2*nlay), Cf(2*nlay), Df(2*nlay), xkk(2*nlay)
    real(kind=f) :: Cwrk(2*nlay), Dwrk(2*nlay)
    real(kind=f) :: xk1(nlay), xk2(nlay)

    real(kind=f) :: up_td(nlay+1), dn_td(nlay+1)
    real(kind=f) :: up_mid(nlay), dn_mid(nlay)

    real(kind=f) :: g_sq, g1, g2, g3, g4, lam, denom, exptrm
    real(kind=f) :: bsurf, dir_k, dir_kp1, dir_mid
    real(kind=f) :: sqrt3, sqrt3d2, btop
    logical :: all_zero_w, uniform_mu

    nlev = nlay + 1
    l    = 2 * nlay

    sqrt3   = sqrt(3._f)
    sqrt3d2 = sqrt3 / 2._f
    btop    = 0._f

    ! ---- reverse to top-down -----------------------------------------------
    do k = 1, nlay
      i = nlay - k + 1
      dtau_raw(k) = dtau_in(i)
      w_raw(k)    = w0_in(i)
      g_raw(k)    = g_in(i)
    end do

    do k = 1, nlev
      mu(k) = mu_in(nlev - k + 1)
    end do

    ! mu(1) is the top of the column, mu(nlev) the bottom.
    uniform_mu = (mu(nlev) == mu(1))

    ! ---- unscaled cumulative optical depth ---------------------------------
    ! Used by the purely absorbing branch, before delta-Eddington rescaling.
    tau_in(1) = 0._f
    do k = 1, nlay
      tau_in(k+1) = tau_in(k) + dtau_raw(k)
    end do

    all_zero_w = .true.
    do k = 1, nlay
      if (w_raw(k) > WMIN) then
        all_zero_w = .false.
        exit
      end if
    end do

    ! A reflecting surface feeds the diffuse streams, so needs the full solver.
    if (w_surf > WMIN) all_zero_w = .false.

    if (all_zero_w) then

      ! ---- purely absorbing: Beer-Lambert on the beam ----------------------
      if (uniform_mu) then
        do k = 1, nlev
          dn_td(k) = f0 * mu(nlev) * exp(-tau_in(k) / mu(nlev))
        end do
      else
        cum_trans(1) = tau_in(1) / mu(1)
        do k = 1, nlev - 1
          cum_trans(k+1) = cum_trans(k) + dtau_raw(k) / mu(k+1)
        end do
        do k = 1, nlev
          dn_td(k) = f0 * mu(nlev) * exp(-cum_trans(k))
        end do
      end if

      dn_td(nlev) = dn_td(nlev) * (1._f - w_surf)
      up_td(:)    = 0._f

    else

      ! ---- delta-Eddington rescaling ---------------------------------------
      do k = 1, nlay
        g_sq    = g_raw(k) * g_raw(k)
        w0(k)   = ((1._f - g_sq) * w_raw(k)) / (1._f - w_raw(k) * g_sq)
        dtau(k) = (1._f - w_raw(k) * g_sq) * dtau_raw(k)
        hg(k)   = g_raw(k) / (1._f + g_raw(k))
      end do

      tau(1) = 0._f
      do k = 1, nlay
        tau(k+1) = tau(k) + dtau(k)
      end do

      ! ---- direct beam ------------------------------------------------------
      if (uniform_mu) then
        do k = 1, nlev
          dir(k) = f0 * mu(nlev) * exp(-tau(k) / mu(nlev))
        end do
        do k = 1, nlay
          mu_zm(k) = mu(nlev)
        end do
      else
        cum_trans(1) = tau(1) / mu(1)
        do k = 1, nlev - 1
          cum_trans(k+1) = cum_trans(k) + (tau(k+1) - tau(k)) / mu(k+1)
        end do
        do k = 1, nlev
          dir(k) = f0 * mu(nlev) * exp(-cum_trans(k))
        end do
        do k = 1, nlay
          mu_zm(k) = 0.5_f * (mu(k) + mu(k+1))
        end do
      end if

      ! ---- per-layer two-stream coefficients --------------------------------
      do k = 1, nlay
        g1 = sqrt3d2 * (2._f - w0(k) * (1._f + hg(k)))
        g2 = (sqrt3d2 * w0(k)) * (1._f - hg(k))
        if (g2 == 0._f) g2 = TINY_SUB
        g3 = (1._f - sqrt3 * hg(k) * mu_zm(k)) / 2._f
        g4 = 1._f - g3

        lam    = sqrt(g1 * g1 - g2 * g2)
        gam(k) = (g1 - lam) / g2

        denom = (lam * lam) - 1._f / (mu_zm(k) * mu_zm(k))
        if (denom == 0._f) denom = TINY_SUB

        Ap(k) = f0 * w0(k) * (g3 * (g1 - 1._f / mu_zm(k)) + g2 * g4) / denom
        Am(k) = f0 * w0(k) * (g4 * (g1 + 1._f / mu_zm(k)) + g2 * g3) / denom

        Cpm1(k) = Ap(k) * exp(-tau(k) / mu_zm(k))
        Cmm1(k) = Am(k) * exp(-tau(k) / mu_zm(k))
        Cp(k)   = Ap(k) * exp(-tau(k+1) / mu_zm(k))
        Cm(k)   = Am(k) * exp(-tau(k+1) / mu_zm(k))

        exptrm = min(lam * dtau(k), EXPMAX)
        Ep(k)  = exp(exptrm)
        Em(k)  = 1._f / Ep(k)

        E1(k) = Ep(k) + gam(k) * Em(k)
        E2(k) = Ep(k) - gam(k) * Em(k)
        E3(k) = gam(k) * Ep(k) + Em(k)
        E4(k) = gam(k) * Ep(k) - Em(k)
      end do

      ! ---- tridiagonal system ----------------------------------------------
      Af(1) = 0._f
      Bf(1) = gam(1) + 1._f
      Cf(1) = gam(1) - 1._f
      Df(1) = btop - Cmm1(1)

      n = 1
      do i = 2, l - 2, 2
        Af(i) = (E1(n) + E3(n)) * (gam(n+1) - 1._f)
        Bf(i) = (E2(n) + E4(n)) * (gam(n+1) - 1._f)
        Cf(i) = 2._f * (1._f - gam(n+1) * gam(n+1))
        Df(i) = (gam(n+1) - 1._f) * (Cpm1(n+1) - Cp(n)) &
                + (1._f - gam(n+1)) * (Cm(n) - Cmm1(n+1))
        n = n + 1
      end do

      n = 1
      do i = 3, l - 1, 2
        Af(i) = 2._f * (1._f - gam(n) * gam(n))
        Bf(i) = (E1(n) - E3(n)) * (1._f + gam(n+1))
        Cf(i) = (E1(n) + E3(n)) * (gam(n+1) - 1._f)
        Df(i) = E3(n) * (Cpm1(n+1) - Cp(n)) + E1(n) * (Cm(n) - Cmm1(n+1))
        n = n + 1
      end do

      ! Bottom boundary: the surface reflects the diffuse streams and the
      ! direct beam (Toon89 eqn 37).
      bsurf = w_surf * dir(nlev)
      Af(l) = E1(nlay) - w_surf * E3(nlay)
      Bf(l) = E2(nlay) - w_surf * E4(nlay)
      Cf(l) = 0._f
      Df(l) = bsurf - Cp(nlay) + w_surf * Cm(nlay)

      Cwrk(:) = Cf(:)
      Dwrk(:) = Df(:)
      call dtridgl(l, Af, Bf, Cwrk, Dwrk, xkk)

      ! ---- fluxes from the solution vector ----------------------------------
      do n = 1, nlay
        xk1(n) = xkk(2*n - 1) + xkk(2*n)
        xk2(n) = xkk(2*n - 1) - xkk(2*n)

        ! Catastrophic cancellation guard.
        if (abs(xk2(n)) < 1.e-30_f * abs(xkk(2*n))) xk2(n) = 0._f

        up_td(n) = xk1(n) + gam(n) * xk2(n) + Cpm1(n)
        dn_td(n) = xk1(n) * gam(n) + xk2(n) + Cmm1(n)
      end do

      up_td(nlev) = xk1(nlay) * Ep(nlay) &
                    + gam(nlay) * xk2(nlay) * Em(nlay) + Cp(nlay)
      dn_td(nlev) = xk1(nlay) * Ep(nlay) * gam(nlay) &
                    + xk2(nlay) * Em(nlay) + Cm(nlay)

      ! f_dn is diffuse plus direct beam.
      do k = 1, nlev
        dn_td(k) = dn_td(k) + dir(k)
      end do

    end if

    ! ---- layer midpoints ----------------------------------------------------
    ! Diffuse streams take the upper level value (O(dtau)); the direct beam
    ! takes its geometric mean across the layer, exact for an exponential.
    if (present(f_up_mid) .or. present(f_dn_mid)) then
      do k = 1, nlay
        up_mid(k) = up_td(k)

        if (all_zero_w) then
          dir_k   = dn_td(k)       ! with no scattering the flux *is* the beam
          dir_kp1 = dn_td(k+1)
        else
          dir_k   = dir(k)
          dir_kp1 = dir(k+1)
        end if

        dir_mid   = sqrt(max(dir_k * dir_kp1, 0._f))
        dn_mid(k) = dn_td(k) - dir_k + dir_mid
      end do
    end if

    ! ---- back to bottom-to-top ordering -------------------------------------
    do k = 1, nlev
      f_up(k) = up_td(nlev - k + 1)
      f_dn(k) = dn_td(nlev - k + 1)
    end do

    if (present(f_up_mid)) then
      do k = 1, nlay
        f_up_mid(k) = up_mid(nlay - k + 1)
      end do
    end if

    if (present(f_dn_mid)) then
      do k = 1, nlay
        f_dn_mid(k) = dn_mid(nlay - k + 1)
      end do
    end if

    return
  end subroutine toon_sw_column

end module carma_swsolve
