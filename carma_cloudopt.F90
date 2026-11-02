!! Cloud optical properties from CARMA's own particle state, for the
!! radiatively coupled version of CARMA.
!!
!! The Mie tables ``qext``, ``ssa`` and ``asym`` are generated on the run's bin
!! grid and on the correlated-k band centres (``carmapy.radiation.gen_band_mie_table``),
!! so the only per-step work is the number-density-weighted sum over bins. The
!! caller broadcasts each band's values across its g-points. All units are cgs.
!!
!! Approximations: Mie properties are evaluated at the dry bin radius, ignoring
!! ``cstate%f_r_wet``; and a heterogeneous group uses the mantle species'
!! refractive index regardless of the core/mantle mass ratio.
!!
!! @author Wolf Cukier
module carma_cloudopt

  use carma_precision_mod

  implicit none

  private
  public :: cloud_optics_column

  real(kind=f), parameter :: PI_F = 3.14159265358979323846_f

contains

  !! Column cloud optical depth, single scattering albedo and asymmetry factor,
  !! summed over every bin of every group.
  !!
  !!   beta_ext(iw,iz) = sum_g sum_b numden(iz,b,g) * pi * r(b,g)^2 * qext(iw,b,g)
  !!
  !! with ``beta_sca`` likewise using ``qext*ssa``, and the asymmetry factor
  !! accumulated as a scattering-weighted mean. The wavelength index is
  !! innermost because it is the contiguous dimension of the optics arrays.
  subroutine cloud_optics_column(nz, nbin, nwave, ngroup, numden, radius, dz, &
                                 qext, ssa, asym, tau, w0, gasym)

    integer, intent(in)       :: nz                              !! number of layers
    integer, intent(in)       :: nbin                            !! number of size bins
    integer, intent(in)       :: nwave                           !! number of bands
    integer, intent(in)       :: ngroup                          !! number of particle groups
    real(kind=f), intent(in)  :: numden(nz, nbin, ngroup)        !! particle number density [#/cm^3]
    real(kind=f), intent(in)  :: radius(nbin, ngroup)            !! dry bin radius [cm]
    real(kind=f), intent(in)  :: dz(nz)                          !! layer thickness [cm]
    real(kind=f), intent(in)  :: qext(nwave, nbin, ngroup)       !! extinction efficiency
    real(kind=f), intent(in)  :: ssa(nwave, nbin, ngroup)        !! single scattering albedo
    real(kind=f), intent(in)  :: asym(nwave, nbin, ngroup)       !! asymmetry factor
    real(kind=f), intent(out) :: tau(nwave, nz)                  !! cloud optical depth
    real(kind=f), intent(out) :: w0(nwave, nz)                   !! cloud single scattering albedo
    real(kind=f), intent(out) :: gasym(nwave, nz)                !! cloud asymmetry factor

    integer      :: iz, ibin, igroup, iw
    real(kind=f) :: xsec, qs
    real(kind=f) :: bext(nwave), bsca(nwave), bg(nwave)

    ! Layers are independent; the band accumulators are per-thread.
    !$OMP PARALLEL DO &
    !$OMP& PRIVATE(iz, ibin, igroup, iw, xsec, qs, bext, bsca, bg) &
    !$OMP& SCHEDULE(static)
    do iz = 1, nz

      bext(:) = 0._f
      bsca(:) = 0._f
      bg(:)   = 0._f

      do igroup = 1, ngroup
        do ibin = 1, nbin

          ! Geometric cross section per unit volume [1/cm]; skip empty bins.
          xsec = numden(iz, ibin, igroup) * PI_F * radius(ibin, igroup)**2
          if (xsec <= 0._f) cycle

          do iw = 1, nwave
            qs = qext(iw, ibin, igroup) * ssa(iw, ibin, igroup)
            bext(iw) = bext(iw) + xsec * qext(iw, ibin, igroup)
            bsca(iw) = bsca(iw) + xsec * qs
            bg(iw)   = bg(iw)   + xsec * qs * asym(iw, ibin, igroup)
          end do

        end do
      end do

      do iw = 1, nwave
        tau(iw, iz) = bext(iw) * dz(iz)

        ! A cloud-free layer gets tau = w0 = g = 0.
        if (bext(iw) > 0._f) then
          w0(iw, iz) = bsca(iw) / bext(iw)
        else
          w0(iw, iz) = 0._f
        end if

        if (bsca(iw) > 0._f) then
          gasym(iw, iz) = bg(iw) / bsca(iw)
        else
          gasym(iw, iz) = 0._f
        end if
      end do

    end do
    !$OMP END PARALLEL DO

  end subroutine cloud_optics_column

end module carma_cloudopt
