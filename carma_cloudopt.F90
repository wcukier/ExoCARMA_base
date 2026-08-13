!! Cloud optical properties from CARMA's own particle state, for the
!! radiatively coupled version of CARMA.
!!
!! This is the Fortran counterpart of ``_get_cloud_opacities`` in
!! ``src/carmapy/results.py``, and deliberately reproduces it rather than
!! redesigning it. The one structural difference is where the size dependence
!! is resolved: the Python path interpolates a Mie table onto the bin radii with
!! a spline on every call, whereas here the table has already been generated on
!! *exactly* the run's bin grid and band set (see
!! ``carmapy.radiation.gen_band_mie_table``). Nothing is interpolated at run
!! time; ``qext``, ``ssa`` and ``asym`` are static ``(NWAVE, NBIN)`` arrays and
!! the only per-step work is the number-density-weighted sum over bins.
!!
!! Optical properties are evaluated on the correlated-k **band centres**, not on
!! the g-points. Mie properties vary smoothly in wavelength and the g-points
!! within a band share a wavelength range, so this is an 8x saving on the
!! dominant cost of the sum and introduces no error the band discretization does
!! not already have. The caller broadcasts each band's values across its
!! g-points when assembling the optical depth handed to ``carma_rtsolve``.
!!
!! Everything here is cgs, matching the rest of CARMA. Optical depth is
!! dimensionless, so nothing crosses the SI boundary that ``carma_rtsolve``
!! guards.
!!
!! Two approximations are inherited from the Python implementation and are
!! stated rather than fixed: Mie properties are evaluated at the **dry** bin
!! radius, ignoring ``cstate%f_r_wet``; and a heterogeneous group uses a single
!! refractive index (that of the mantle species) regardless of the core/mantle
!! mass ratio, which varies with both height and time.
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
  !! The per-bin extinction cross section is the geometric cross section times
  !! the Mie efficiency, so the bin sum is
  !!
  !!   beta_ext(iw,iz) = sum_g sum_b numden(iz,b,g) * pi * r(b,g)^2 * qext(iw,b,g)
  !!
  !! with ``beta_sca`` likewise using ``qext*ssa``, and the asymmetry factor
  !! accumulated as a scattering-weighted mean. This matches results.py's
  !! ``weighted_qext`` / ``weighted_qsca`` / ``weighted_g`` sums exactly.
  !!
  !! The wavelength index is innermost and is the contiguous dimension of the
  !! optics arrays, making the bin sum a GEMV per layer. Written the other way
  !! round -- wavelength outermost -- this loop nest dominates the timestep.
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

    do iz = 1, nz

      bext(:) = 0._f
      bsca(:) = 0._f
      bg(:)   = 0._f

      do igroup = 1, ngroup
        do ibin = 1, nbin

          ! Geometric cross section per unit volume for this bin. Empty bins
          ! are the common case high in the column, and skipping them here is
          ! what keeps the nwave loop off the critical path.
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

        ! A cloud-free layer has no scattering to weight with; leaving w0 and g
        ! at zero makes the layer purely absorbing with zero optical depth,
        ! which is the correct no-op for the solver.
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

  end subroutine cloud_optics_column

end module carma_cloudopt
