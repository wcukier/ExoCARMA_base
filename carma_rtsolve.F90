!! Two-stream longwave radiative transfer for the radiatively coupled version
!! of CARMA.
!!
!! This is a transliteration of pyHARP's ``toon_mckay89_longwave``
!! (``pyharp/src/rtsolver/toon_mckay89_longwave_impl.h``, MIT licensed), after
!! Toon, McKay & Ackerman (1989), which pyHARP in turn derived from Elsie Lee's
!! ``Exo-FMS_column_ck``. Keeping it a faithful line-for-line port is what makes
!! the Fortran verifiable against pyHARP directly.
!!
!! Why not DISORT: ``pydisort`` wraps ``cdisort`` (a C rewrite of DISORT 2.1)
!! which is GPL-3, incompatible with CARMApy's Apache-2.0 licence. The solver
!! sits behind a deliberately narrow interface -- optical properties and a
!! Planck source in, level fluxes out -- so a permissively licensed DISORT can
!! be added later as an alternative backend without touching callers.
!!
!! Conventions:
!!  * All arrays passed in and out are **bottom-to-top**, matching CARMA's
!!    ``p(1)`` = deepest. The Toon recursion runs top-down, so the reversal is
!!    done explicitly on entry and undone on exit; nothing outside this module
!!    sees the reversed ordering.
!!  * Units are **SI** (W/m^2), matching pyHARP, so results compare directly.
!!    The caller converts to cgs.
!!
!! @author Wolf Cukier
module carma_rtsolve

  use carma_precision_mod

  implicit none

  private
  public :: toon_lw_column
  !! Shared with ``carma_swsolve``: both two-stream solvers close on the same
  !! tridiagonal system and the reference has one implementation of it, so
  !! duplicating it here would be two copies to keep in step for no gain.
  public :: dtridgl

  real(kind=f), parameter :: PI_F  = 3.14159265358979323846_f
  real(kind=f), parameter :: TWOPI = 2._f * PI_F
  real(kind=f), parameter :: UBARI = 0.5_f

  !! Cap on the exponent argument, as in the reference. Without it, an
  !! optically thick layer overflows exp() before the tridiagonal solve can
  !! cancel the large terms against each other.
  real(kind=f), parameter :: EXPMAX = 35._f

  integer, parameter :: NMU = 5

  !! Gauss quadrature points and weights for the angular integration.
  real(kind=f), parameter :: UARR(NMU) = &
      (/ 0.0985350858_f, 0.3045357266_f, 0.5620251898_f, &
         0.8019865821_f, 0.9601901429_f /)
  real(kind=f), parameter :: WUARR(NMU) = &
      (/ 0.0157479145_f, 0.0739088701_f, 0.1463869871_f, &
         0.1671746381_f, 0.0967815902_f /)

contains

  !! Thomas algorithm for a tridiagonal system, matching pyHARP's ``dtridgl``.
  !!
  !! ``c`` and ``d`` are modified in place, exactly as in the reference (the
  !! caller passes copies).
  pure subroutine dtridgl(n, a, b, c, d, x)

    implicit none

    integer, intent(in)         :: n
    real(kind=f), intent(in)    :: a(n), b(n)
    real(kind=f), intent(inout) :: c(n), d(n)
    real(kind=f), intent(out)   :: x(n)

    integer      :: i
    real(kind=f) :: denom

    c(1) = c(1) / b(1)
    d(1) = d(1) / b(1)

    do i = 2, n
      denom = b(i) - a(i) * c(i-1)
      if (denom == 0._f) denom = 1.e-12_f   ! guard, as in the reference

      if (i < n) then
        c(i) = c(i) / denom
      else
        c(i) = 0._f
      end if

      d(i) = (d(i) - a(i) * d(i-1)) / denom
    end do

    x(n) = d(n)

    do i = n-1, 1, -1
      x(i) = d(i) - c(i) * x(i+1)
    end do

    return
  end subroutine dtridgl


  !! Solve the longwave two-stream problem for a single column at a single
  !! spectral point.
  !!
  !! All input and output arrays are bottom-to-top: index 1 is the deepest
  !! layer/level, index nlay (or nlay+1) is the top of the atmosphere.
  subroutine toon_lw_column(nlay, dtau_in, w0_in, g_in, be_in, a_surf, &
                            top_emission_flag, btop_factor, hard_surface, &
                            delta_eddington_lw, f_up, f_dn, &
                            f_up_mid, f_dn_mid, be_corr_in)

    implicit none

    integer, intent(in)      :: nlay               !! number of layers
    real(kind=f), intent(in) :: dtau_in(nlay)      !! layer optical thickness
    real(kind=f), intent(in) :: w0_in(nlay)        !! single scattering albedo
    real(kind=f), intent(in) :: g_in(nlay)         !! asymmetry parameter
    real(kind=f), intent(in) :: be_in(nlay+1)      !! level Planck intensity [W/m^2/sr]
    real(kind=f), intent(in) :: a_surf             !! surface albedo
    integer, intent(in)      :: top_emission_flag  !! <0 auto, 0 none, 1 full Planck
    real(kind=f), intent(in) :: btop_factor        !! scales tau_top when auto
    logical, intent(in)      :: hard_surface       !! .false. for a gas giant
    logical, intent(in)      :: delta_eddington_lw !! rescale w0/dtau/g

    real(kind=f), intent(out) :: f_up(nlay+1)      !! upward flux [W/m^2]
    real(kind=f), intent(out) :: f_dn(nlay+1)      !! downward flux [W/m^2]
    real(kind=f), intent(out), optional :: f_up_mid(nlay)  !! layer-midpoint upward flux
    real(kind=f), intent(out), optional :: f_dn_mid(nlay)  !! layer-midpoint downward flux

    !! The part of each layer's emission its own bounding levels cannot see
    !! [W/m^2/sr]: the Planck intensity at the layer's temperature less the
    !! Planck intensity the level interpolation implies there.
    !!
    !! The levels are interpolated from the layer centres, so for any profile
    !! linear in ln p they reproduce the centre exactly and this is zero -- the
    !! source below is then the one built from levels alone, unchanged. It
    !! departs from zero only where a layer's temperature differs from what its
    !! neighbours imply, which an interpolation cannot represent: a
    !! layer-to-layer oscillation lives entirely in the part an average throws
    !! away. Without this term such a layer emits at its neighbours'
    !! temperature rather than its own, so one colder than its surroundings
    !! radiates as though it were warm and cools further -- a positive feedback
    !! at the grid scale, which manufactures the oscillation it feeds on.
    real(kind=f), intent(in), optional :: be_corr_in(nlay)

    integer :: nlev, l, k, n, m, i

    ! top-down working copies
    real(kind=f) :: dtau(nlay), w0(nlay), hg(nlay), be(nlay+1), be_c(nlay)
    logical      :: have_corr
    real(kind=f) :: alp(nlay), lam(nlay), gam(nlay), term(nlay)
    real(kind=f) :: B0(nlay), B1(nlay)
    real(kind=f) :: Cpm1(nlay), Cmm1(nlay), Cp(nlay), Cm(nlay)
    real(kind=f) :: Ep(nlay), Em(nlay), E1(nlay), E2(nlay), E3(nlay), E4(nlay)
    real(kind=f) :: em1(nlay)

    real(kind=f) :: Af(2*nlay), Bf(2*nlay), Cf(2*nlay), Df(2*nlay), xkk(2*nlay)
    real(kind=f) :: Cwrk(2*nlay), Dwrk(2*nlay)
    real(kind=f) :: xk1(nlay), xk2(nlay)

    real(kind=f) :: g(nlay), h(nlay), xj(nlay), xk(nlay)
    real(kind=f) :: alpha1(nlay), alpha2(nlay), sigma1(nlay), sigma2(nlay)

    real(kind=f) :: up_g(nlay+1), dn_g(nlay+1)
    real(kind=f) :: up_lev(nlay+1), dn_lev(nlay+1)
    real(kind=f) :: up_mid(nlay), dn_mid(nlay)

    real(kind=f) :: g_sq, exptrm, Btop, Btop_g, Bsurf, bsurf_flux
    real(kind=f) :: u, em2, em2_mid, Ep_mid, em1_mid, em3
    real(kind=f) :: l_u_p1, l_u_m1, common_den, term_val, mid_val

    nlev = nlay + 1
    l    = 2 * nlay

    ! ---- reverse to top-down and apply optional delta-Eddington scaling -----
    ! When disabled (the default) the raw values are used, as in PICASO.
    do k = 1, nlay
      i = nlay - k + 1

      if (delta_eddington_lw) then
        g_sq    = g_in(i) * g_in(i)
        w0(k)   = ((1._f - g_sq) * w0_in(i)) / (1._f - w0_in(i) * g_sq)
        dtau(k) = (1._f - w0_in(i) * g_sq) * dtau_in(i)
        hg(k)   = g_in(i) / (1._f + g_in(i))
      else
        w0(k)   = w0_in(i)
        dtau(k) = dtau_in(i)
        hg(k)   = g_in(i)
      end if
    end do

    do k = 1, nlev
      be(k) = be_in(nlev - k + 1)
    end do

    have_corr = present(be_corr_in)

    be_c(:) = 0._f
    if (have_corr) then
      do k = 1, nlay
        be_c(k) = be_corr_in(nlay - k + 1)
      end do
    end if

    ! ---- per-layer two-stream coefficients ----------------------------------
    do k = 1, nlay
      alp(k)  = sqrt((1._f - w0(k)) / (1._f - w0(k) * hg(k)))
      lam(k)  = alp(k) * (1._f - w0(k) * hg(k)) / UBARI
      gam(k)  = (1._f - alp(k)) / (1._f + alp(k))
      term(k) = UBARI / (1._f - w0(k) * hg(k))

      ! An optically negligible layer has no meaningful Planck gradient; using
      ! the midpoint value avoids dividing by ~0.
      !
      ! The source is linear in tau across the layer, B(tau) = B0 + B1*tau: the
      ! levels set the slope, and `be_c` adds back the part of the layer's own
      ! emission they cannot represent. It is zero for a resolved profile, so
      ! this is the level-only source wherever the column is smooth.
      if (dtau(k) <= 1.e-6_f) then
        B1(k) = 0._f
        B0(k) = 0.5_f * (be(k+1) + be(k)) + be_c(k)
      else
        B1(k) = (be(k+1) - be(k)) / dtau(k)
        B0(k) = be(k) + be_c(k)
      end if

      Cpm1(k) = B0(k) + B1(k) * term(k)
      Cmm1(k) = B0(k) - B1(k) * term(k)
      Cp(k)   = B0(k) + B1(k) * dtau(k) + B1(k) * term(k)
      Cm(k)   = B0(k) + B1(k) * dtau(k) - B1(k) * term(k)

      exptrm = min(lam(k) * dtau(k), EXPMAX)
      Ep(k)  = exp(exptrm)
      Em(k)  = 1._f / Ep(k)
      E1(k)  = Ep(k) + gam(k) * Em(k)
      E2(k)  = Ep(k) - gam(k) * Em(k)
      E3(k)  = gam(k) * Ep(k) + Em(k)
      E4(k)  = gam(k) * Ep(k) - Em(k)
    end do

    ! ---- boundary conditions -------------------------------------------------
    if (top_emission_flag < 0) then
      ! Auto: an effective optical depth above the model top, so a finite
      ! downward flux enters at TOA. btop_factor = p(1)/(p(2)-p(1)) reproduces
      ! PICASO's get_thermal_1d.
      Btop = (1._f - exp(-(dtau(1) * btop_factor) / UBARI)) * be(1)
    else
      Btop = real(top_emission_flag, kind=f) * be(1)
    end if

    Bsurf = be(nlev)

    if (hard_surface) then
      ! Terrestrial: surface emits as a blackbody reduced by emissivity.
      bsurf_flux = (1._f - a_surf) * Bsurf
    else
      ! Gas giant: Planck plus its gradient extrapolated below the deepest
      ! level. This is what lets the internal heat flux emerge, and is the
      ! branch a brown dwarf needs.
      bsurf_flux = Bsurf + B1(nlay) * UBARI
    end if

    ! ---- tridiagonal system --------------------------------------------------
    Af(1) = 0._f
    Bf(1) = gam(1) + 1._f
    Cf(1) = gam(1) - 1._f
    Df(1) = Btop - Cmm1(1)

    do n = 1, nlay - 1
      Af(2*n) = (E1(n) + E3(n)) * (gam(n+1) - 1._f)
      Bf(2*n) = (E2(n) + E4(n)) * (gam(n+1) - 1._f)
      Cf(2*n) = 2._f * (1._f - gam(n+1) * gam(n+1))
      Df(2*n) = (gam(n+1) - 1._f) * (Cpm1(n+1) - Cp(n)) &
              + (1._f - gam(n+1)) * (Cm(n) - Cmm1(n+1))

      Af(2*n+1) = 2._f * (1._f - gam(n) * gam(n))
      Bf(2*n+1) = (E1(n) - E3(n)) * (1._f + gam(n+1))
      Cf(2*n+1) = (E1(n) + E3(n)) * (gam(n+1) - 1._f)
      Df(2*n+1) = E3(n) * (Cpm1(n+1) - Cp(n)) &
                + E1(n) * (Cm(n) - Cmm1(n+1))
    end do

    Af(l) = E1(nlay) - a_surf * E3(nlay)
    Bf(l) = E2(nlay) - a_surf * E4(nlay)
    Cf(l) = 0._f
    Df(l) = bsurf_flux - Cp(nlay) + a_surf * Cm(nlay)

    Cwrk(:) = Cf(:)
    Dwrk(:) = Df(:)

    call dtridgl(l, Af, Bf, Cwrk, Dwrk, xkk)

    ! ---- source function terms ----------------------------------------------
    do n = 1, nlay
      xk1(n) = xkk(2*n-1) + xkk(2*n)
      xk2(n) = xkk(2*n-1) - xkk(2*n)

      ! Catastrophic cancellation guard, as in the reference.
      if (abs(xk2(n)) < 1.e-30_f * abs(xkk(2*n))) xk2(n) = 0._f

      if (w0(n) <= 1.e-4_f) then
        ! Effectively pure absorption: the scattering source vanishes.
        g(n)      = 0._f
        h(n)      = 0._f
        xj(n)     = 0._f
        xk(n)     = 0._f
        alpha1(n) = TWOPI * B0(n)
        alpha2(n) = TWOPI * B1(n)
        sigma1(n) = alpha1(n)
        sigma2(n) = alpha2(n)
      else
        common_den = 1._f + alp(n)
        g(n)  = TWOPI * w0(n) * xk1(n) * (1._f + hg(n) * alp(n)) / common_den
        h(n)  = TWOPI * w0(n) * xk2(n) * (1._f - hg(n) * alp(n)) / common_den
        xj(n) = TWOPI * w0(n) * xk1(n) * (1._f - hg(n) * alp(n)) / common_den
        xk(n) = TWOPI * w0(n) * xk2(n) * (1._f + hg(n) * alp(n)) / common_den

        term_val  = UBARI * w0(n) * hg(n) / (1._f - w0(n) * hg(n))
        alpha1(n) = TWOPI * (B0(n) + B1(n) * term_val)
        alpha2(n) = TWOPI * B1(n)
        sigma1(n) = TWOPI * (B0(n) - B1(n) * term_val)
        sigma2(n) = alpha2(n)
      end if

      em1(n) = 1._f / exp(min(lam(n) * dtau(n), EXPMAX))
    end do

    ! ---- angular quadrature --------------------------------------------------
    up_lev(:) = 0._f
    dn_lev(:) = 0._f
    up_mid(:) = 0._f
    dn_mid(:) = 0._f

    do m = 1, NMU
      u = UARR(m)

      ! downward sweep, TOA to bottom
      if (top_emission_flag < 0) then
        Btop_g = (1._f - exp(-(dtau(1) * btop_factor) / u)) * be(1)
      else
        Btop_g = real(top_emission_flag, kind=f) * be(1)
      end if

      dn_g(1) = TWOPI * Btop_g

      do k = 1, nlay
        em2    = exp(-dtau(k) / u)
        l_u_p1 = lam(k) * u + 1._f
        l_u_m1 = lam(k) * u - 1._f

        dn_g(k+1) = dn_g(k) * em2 &
                  + (xj(k) / l_u_p1) * (Ep(k) - em2) &
                  + (xk(k) / l_u_m1) * (em2 - em1(k)) &
                  + sigma1(k) * (1._f - em2) &
                  + sigma2(k) * (u * em2 + dtau(k) - u)

        ! Same source-function integral, stopped at half the layer's optical
        ! depth (PICASO's flux_minus_mdpt).
        em2_mid = exp(-0.5_f * dtau(k) / u)
        Ep_mid  = exp(min(0.5_f * lam(k) * dtau(k), EXPMAX))
        em1_mid = 1._f / Ep_mid

        mid_val = dn_g(k) * em2_mid &
                + (xj(k) / l_u_p1) * (Ep_mid - em2_mid) &
                + (xk(k) / l_u_m1) * (em2_mid - em1_mid) &
                + sigma1(k) * (1._f - em2_mid) &
                + sigma2(k) * (u * em2_mid + 0.5_f * dtau(k) - u)

        dn_mid(k) = dn_mid(k) + mid_val * WUARR(m)
      end do

      ! upward sweep, bottom to TOA
      if (hard_surface) then
        up_g(nlev) = TWOPI * (1._f - a_surf) * Bsurf
      else
        up_g(nlev) = TWOPI * (Bsurf + B1(nlay) * u)
      end if

      do k = nlay, 1, -1
        em2    = exp(-dtau(k) / u)
        em3    = em1(k) * em2
        l_u_m1 = lam(k) * u - 1._f
        l_u_p1 = lam(k) * u + 1._f

        up_g(k) = up_g(k+1) * em2 &
                + (g(k) / l_u_m1) * (Ep(k) * em2 - 1._f) &
                + (h(k) / l_u_p1) * (1._f - em3) &
                + alpha1(k) * (1._f - em2) &
                + alpha2(k) * (u - (dtau(k) + u) * em2)

        em2_mid = exp(-0.5_f * dtau(k) / u)
        Ep_mid  = exp(min(0.5_f * lam(k) * dtau(k), EXPMAX))
        em1_mid = 1._f / Ep_mid

        mid_val = up_g(k+1) * em2_mid &
                + (g(k) / l_u_m1) * (Ep(k) * em2_mid - Ep_mid) &
                + (h(k) / l_u_p1) * (em1_mid - em1(k) * em2_mid) &
                + alpha1(k) * (1._f - em2_mid) &
                + alpha2(k) * (u + 0.5_f * dtau(k) - (dtau(k) + u) * em2_mid)

        up_mid(k) = up_mid(k) + mid_val * WUARR(m)
      end do

      do k = 1, nlev
        dn_lev(k) = dn_lev(k) + dn_g(k) * WUARR(m)
        up_lev(k) = up_lev(k) + up_g(k) * WUARR(m)
      end do
    end do

    ! ---- back to bottom-to-top ordering --------------------------------------
    do k = 1, nlev
      f_up(k) = up_lev(nlev - k + 1)
      f_dn(k) = dn_lev(nlev - k + 1)
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
  end subroutine toon_lw_column

end module carma_rtsolve
