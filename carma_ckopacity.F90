!! Correlated-k gas opacity for the radiatively coupled version of CARMA.
!!
!! Reads the flat binary written by ``carmapy.radiation.export_ck_table`` --
!! Sonora 2020 correlated-k data, originally a TorchScript container produced by
!! pyHARP -- and interpolates it onto the model's (pressure, temperature)
!! profile.
!!
!! Contract, taken from ``pyharp/python/sonora/sonora.py`` and
!! ``pyharp/src/opacity/multiband.cpp``:
!!
!!  * Spectral points are band-major: ``iwave = (iband-1)*ng + ig``.
!!  * ``kappa`` is stored as **ln(cm^2 / molecule)**; pyHARP converts to an
!!    attenuation coefficient with ``1e-4 * N_A * exp(kappa) * conc``, where
!!    ``conc`` is a mole concentration in mol/m^3. That gives 1/m.
!!  * The temperature axis is **absolute** temperature.
!!  * pyHARP interpolates with ``extrapolate=false``, so values are **clamped**
!!    at the table edges rather than extrapolated. This matters for brown
!!    dwarfs: ``extend_atmosphere`` reaches well past the table's deepest
!!    pressure, and those layers are convective anyway, so holding kappa at the
!!    edge value is both faithful to the reference and physically harmless.
!!
!! @author Wolf Cukier
module carma_ckopacity

  use carma_precision_mod

  implicit none

  private
  public :: ck_table_type, ck_load, ck_destroy, ck_kappa_column

  ! Must match _CK_MAGIC / _CK_VERSION in carmapy/radiation.py.
  integer(kind=8), parameter :: CK_MAGIC   = int(z'43524D41434B3031', kind=8)
  integer(kind=8), parameter :: CK_VERSION = 1_8

  !! Avogadro's number, matching harp::constants::Avogadro.
  real(kind=f), parameter, public :: CK_AVOGADRO = 6.02214076e23_f

  !! A loaded correlated-k table.
  type ck_table_type
    integer :: nband = 0   !! number of spectral bands
    integer :: ng    = 0   !! g-points per band
    integer :: nwave = 0   !! nband * ng
    integer :: npres = 0   !! pressure axis length
    integer :: ntemp = 0   !! temperature axis length

    real(kind=f), allocatable :: wmin(:)        !! (nband) band lower edge [cm^-1]
    real(kind=f), allocatable :: wmax(:)        !! (nband) band upper edge [cm^-1]
    real(kind=f), allocatable :: gauss_pts(:)   !! (ng)
    real(kind=f), allocatable :: gauss_wts(:)   !! (ng), sums to 1 per band
    real(kind=f), allocatable :: wavenumber(:)  !! (nwave) [cm^-1]
    real(kind=f), allocatable :: weights(:)     !! (nwave) spectral integration weights
    real(kind=f), allocatable :: lnp(:)         !! (npres) ln(pressure [Pa])
    real(kind=f), allocatable :: temp(:)        !! (ntemp) [K]
    real(kind=f), allocatable :: kappa(:,:,:)   !! (nwave, npres, ntemp) ln(cm^2/molecule)
  end type ck_table_type

contains

  !! Load a correlated-k table from the binary written by
  !! ``carmapy.radiation.export_ck_table``.
  !!
  !! The file is raw stream-access float64/int64, written in Fortran
  !! (column-major) order, so it maps straight onto the arrays below.
  subroutine ck_load(path, tbl, rc)

    implicit none

    character(len=*), intent(in)       :: path  !! path to the .ck file
    type(ck_table_type), intent(inout) :: tbl   !! table to fill
    integer, intent(inout)             :: rc    !! negative on failure

    integer(kind=8) :: magic, version, nband8, ng8, npres8, ntemp8
    integer         :: lun, ios, i
    real(kind=f), allocatable :: pres(:)

    rc = 0

    open(newunit=lun, file=trim(path), access='stream', form='unformatted', &
         status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ck_load::ERROR - cannot open ck table: ', trim(path)
      rc = -1
      return
    end if

    read(lun, iostat=ios) magic
    if (ios /= 0 .or. magic /= CK_MAGIC) then
      write(*,*) 'ck_load::ERROR - not a carmapy ck table, or written on a ', &
                 'machine of the opposite endianness: ', trim(path)
      close(lun)
      rc = -1
      return
    end if

    read(lun, iostat=ios) version, nband8, ng8, npres8, ntemp8
    if (ios /= 0) then
      write(*,*) 'ck_load::ERROR - truncated header: ', trim(path)
      close(lun)
      rc = -1
      return
    end if

    if (version /= CK_VERSION) then
      write(*,*) 'ck_load::ERROR - ck table version', version, &
                 'but this build reads version', CK_VERSION
      close(lun)
      rc = -1
      return
    end if

    tbl%nband = int(nband8)
    tbl%ng    = int(ng8)
    tbl%npres = int(npres8)
    tbl%ntemp = int(ntemp8)
    tbl%nwave = tbl%nband * tbl%ng

    call ck_destroy(tbl)

    allocate(tbl%wmin(tbl%nband), tbl%wmax(tbl%nband), &
             tbl%gauss_pts(tbl%ng), tbl%gauss_wts(tbl%ng), &
             tbl%wavenumber(tbl%nwave), tbl%weights(tbl%nwave), &
             tbl%lnp(tbl%npres), tbl%temp(tbl%ntemp), &
             tbl%kappa(tbl%nwave, tbl%npres, tbl%ntemp), &
             pres(tbl%npres), stat=ios)
    if (ios /= 0) then
      write(*,*) 'ck_load::ERROR - allocation failed for ck table'
      close(lun)
      rc = -1
      return
    end if

    read(lun, iostat=ios) tbl%wmin, tbl%wmax, tbl%gauss_pts, tbl%gauss_wts, &
                          tbl%wavenumber, tbl%weights, pres, tbl%temp, tbl%kappa
    close(lun)

    if (ios /= 0) then
      write(*,*) 'ck_load::ERROR - truncated ck table: ', trim(path)
      rc = -1
      return
    end if

    ! pyHARP interpolates in ln(pressure); do the log once here.
    do i = 1, tbl%npres
      tbl%lnp(i) = log(pres(i))
    end do

    deallocate(pres)

    return
  end subroutine ck_load


  !! Release a table's storage.
  subroutine ck_destroy(tbl)

    implicit none

    type(ck_table_type), intent(inout) :: tbl

    if (allocated(tbl%wmin))       deallocate(tbl%wmin)
    if (allocated(tbl%wmax))       deallocate(tbl%wmax)
    if (allocated(tbl%gauss_pts))  deallocate(tbl%gauss_pts)
    if (allocated(tbl%gauss_wts))  deallocate(tbl%gauss_wts)
    if (allocated(tbl%wavenumber)) deallocate(tbl%wavenumber)
    if (allocated(tbl%weights))    deallocate(tbl%weights)
    if (allocated(tbl%lnp))        deallocate(tbl%lnp)
    if (allocated(tbl%temp))       deallocate(tbl%temp)
    if (allocated(tbl%kappa))      deallocate(tbl%kappa)

    return
  end subroutine ck_destroy


  !! Locate x in an ascending grid, returning the lower bracket index i and the
  !! interpolation fraction w such that the value is
  !! ``(1-w)*v(i) + w*v(i+1)``.
  !!
  !! Outside the grid, w is clamped to 0 or 1 -- matching pyHARP's
  !! ``interpn(..., extrapolate=false)``.
  pure subroutine locate_clamped(grid, x, i, w)

    implicit none

    real(kind=f), intent(in)  :: grid(:)
    real(kind=f), intent(in)  :: x
    integer, intent(out)      :: i
    real(kind=f), intent(out) :: w

    integer :: n, lo, hi, mid

    n = size(grid)

    if (n == 1) then
      i = 1
      w = 0._f
      return
    end if

    if (x <= grid(1)) then
      i = 1
      w = 0._f
      return
    end if

    if (x >= grid(n)) then
      i = n - 1
      w = 1._f
      return
    end if

    ! Binary search for the bracketing interval.
    lo = 1
    hi = n

    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (x >= grid(mid)) then
        lo = mid
      else
        hi = mid
      end if
    end do

    i = lo
    w = (x - grid(lo)) / (grid(lo+1) - grid(lo))

    return
  end subroutine locate_clamped


  !! Gas attenuation coefficient for every spectral point over a whole column.
  !!
  !! Bilinear in (ln p, T) at each g-point. The wavenumber axis of pyHARP's
  !! trilinear ``interpn`` is a no-op here because we always evaluate at the
  !! table's own spectral points, so it is omitted.
  !!
  !! Returns the attenuation coefficient in **1/m**, i.e. pyHARP's
  !! ``1e-4 * N_A * exp(kappa) * conc`` with ``conc`` the gas mole
  !! concentration in mol/m^3.
  subroutine ck_kappa_column(tbl, pres_pa, temp_k, conc, beta)

    implicit none

    type(ck_table_type), intent(in) :: tbl
    real(kind=f), intent(in)  :: pres_pa(:)   !! (nz) layer pressure [Pa]
    real(kind=f), intent(in)  :: temp_k(:)    !! (nz) layer temperature [K]
    real(kind=f), intent(in)  :: conc(:)      !! (nz) gas mole concentration [mol/m^3]
    real(kind=f), intent(out) :: beta(:,:)    !! (nwave, nz) attenuation [1/m]

    integer      :: nz, iz, iw, ip, it
    real(kind=f) :: wp, wt, k00, k01, k10, k11, lnk, scale

    nz = size(pres_pa)

    do iz = 1, nz
      call locate_clamped(tbl%lnp,  log(pres_pa(iz)), ip, wp)
      call locate_clamped(tbl%temp, temp_k(iz),       it, wt)

      ! ln(cm^2/molecule) -> 1/m, folding in the column's concentration.
      scale = 1.e-4_f * CK_AVOGADRO * conc(iz)

      do iw = 1, tbl%nwave
        k00 = tbl%kappa(iw, ip,   it)
        k10 = tbl%kappa(iw, ip+1, it)
        k01 = tbl%kappa(iw, ip,   it+1)
        k11 = tbl%kappa(iw, ip+1, it+1)

        ! Interpolate in log space, as pyHARP does -- it stores ln(kappa) and
        ! exponentiates only after interpolating.
        lnk = (1._f - wp) * (1._f - wt) * k00 + wp * (1._f - wt) * k10 &
            + (1._f - wp) *         wt  * k01 + wp *         wt  * k11

        beta(iw, iz) = scale * exp(lnk)
      end do
    end do

    return
  end subroutine ck_kappa_column

end module carma_ckopacity
