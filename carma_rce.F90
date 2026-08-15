!! Radiative-convective temperature evolution for an isolated brown dwarf.
!!
!! Turns the radiative transfer of ``carma_rtsolve`` into a temperature
!! tendency and applies it to the column, so the P-T profile evolves with the
!! cloud field instead of being held fixed for the whole run.
!!
!! Every layer follows the radiative flux divergence,
!! ``dT/dt = -dF_net/dz / (rho*Cp)``, and a dry convective adjustment then
!! removes whatever superadiabatic gradients that leaves behind, mixing each
!! unstable pair of layers at fixed enthalpy onto the adiabat selected by
!! ``adiabat`` -- either the Parmentier (2015) analytic fit or the tabulated
!! H/He gradient. See ``rce_conv_adj``.
!!
!! The convective zones are therefore an *output*: they relocate themselves as
!! the cloud field reshapes the profile, and there may be more than one, since
!! a cloud deck can drive a detached convective layer aloft over a stable
!! region. ``is_conv`` carries the per-layer answer.
!!
!! There is no incident radiation: the internal flux is imposed as the net flux
!! through the base of the column, so radiative-convective equilibrium is the
!! state in which ``F_TOA == sigma*t_int**4``. For an isolated object that is
!! also the effective temperature, but the two are distinct quantities and only
!! the internal one is an input here.
!!
!! **Cadence.** A full solve costs roughly as much as the microphysics it sits
!! next to, and the cloud field barely moves in one timestep, so the solve runs
!! on an adaptive cadence while the resulting ``dTdt`` is held and applied on
!! *every* step. Temperature therefore evolves continuously. A fresh solve is
!! triggered by accumulated temperature drift, by a change in the column's
!! cloud cross-section, or by a hard ceiling on elapsed steps.
!!
!! **Units.** This module is SI throughout -- Pa, m, K, W/m^2 -- matching the
!! driver's ``p``/``zc``/``zl`` and the radiative transfer modules. The two
!! exceptions are documented at their conversion sites: ``Cp`` arrives in
!! erg/g/K and gravity in cm/s^2 from the namelist, and
!! ``cloud_optics_column`` is cgs.
!!
!! **Level temperatures.** The driver carries layer-centre temperatures, but
!! the Planck source and the hydrostatic regrid both want level values. They
!! are reconstructed by linear interpolation in ``ln p``, with the two
!! endpoints extrapolated. The reconstruction depends only on ``t(:)``, so a
!! restarted run rebuilds the same grid it stopped on; it does not reproduce
!! the level temperatures Python supplied at step 1, which differ from the
!! reconstruction by a sub-percent effect on layer thickness.
module carma_rce

  use carma_precision_mod
  use carma_enums_mod, only : I_CART
  use carma_planck,    only : bbflux_wavenumber, bbflux_wavenumber_col
  use carma_ckopacity, only : ck_table_type, ck_load, ck_destroy, ck_kappa_column, &
                              CK_AVOGADRO
  use carma_cloudopt,  only : cloud_optics_column
  use carma_rtsolve,   only : toon_lw_column
  use carma_swsolve,   only : toon_sw_column
  use carma_linalg,    only : lu_factor, lu_solve
  use carma_rayleigh,  only : RAY_NSPEC, RAY_SOLAR_H2HE, ray_sigma_mix

  implicit none

  private
  public :: rce_type, rce_init, rce_update, rce_destroy
  public :: rce_adiabat, rce_level_temps, rce_regrid_z
  public :: rce_grad_ad, rce_integrate_adiabat, rce_load_adiabat
  public :: rce_conv_adj
  !! Exposed so the Jacobian's bandwidth can be measured from outside: the
  !! banded truncation below is only legitimate if a narrow band captures the
  !! operator, and that is a property of this discretisation, not of radiative
  !! transfer, so it has to be measured rather than assumed.
  public :: rce_fluxes, rce_jacobian, rce_solve
  public :: I_RCE_EQUILIBRIUM, I_RCE_PHYSICAL
  public :: I_ADIABAT_PARMENTIER, I_ADIABAT_TABLE

  !! Heating is accelerated relative to the microphysical clock, to reach
  !! radiative-convective equilibrium in a tractable number of steps.
  integer, parameter :: I_RCE_EQUILIBRIUM = 0
  !! Heating runs on the microphysical clock, so the time history is physical.
  integer, parameter :: I_RCE_PHYSICAL    = 1

  !! Stefan-Boltzmann constant [ W / m^2 / K^4 ]
  real(kind=f), parameter :: SIGMA_SB = 5.670374419e-8_f
  !! Universal gas constant [ J / mol / K ]
  real(kind=f), parameter :: RGAS_SI  = 8.31446261815324_f
  !! Boltzmann's constant [ J / K ]
  real(kind=f), parameter :: KB_SI    = 1.380649e-23_f
  !! Proton mass [ kg ]
  real(kind=f), parameter :: MPROT_SI = 1.67262192369e-27_f

  !! Parmentier et al. (2015) adiabat coefficients: the gradient of the fit is
  !! d ln T / d ln P = PARM_A - PARM_B * T.
  real(kind=f), parameter :: PARM_A = 0.32_f
  real(kind=f), parameter :: PARM_B = 1._f / 30000._f

  !! The analytic fit of Parmentier et al. (2015).
  integer, parameter :: I_ADIABAT_PARMENTIER = 0
  !! The tabulated H/He gradient, integrated from the anchor.
  integer, parameter :: I_ADIABAT_TABLE      = 1

  !! Ceiling on the RK4 step used to integrate the tabulated adiabat, in
  !! ln p. The table's own pressure spacing is 0.2 in log10 p (0.46 in ln p),
  !! so this puts several steps in every table cell.
  real(kind=f), parameter :: ADIABAT_DLNP = 0.05_f

  !! 1 bar in Pa. The table's pressure axis is log10 bar; this module is SI.
  real(kind=f), parameter :: BAR2PA = 1.e5_f

  !! Tolerance on the convective adjustment's stability test, as a fraction of
  !! the pair's temperature: a pair this close to the adiabat counts as
  !! neutral. Relative rather than the reference's absolute 1e-6 K, which at a
  !! brown dwarf's 3000 K interior would be 3e-10 -- within roundoff of a chain
  !! of sixty layers. It is what decides ``is_conv``; in a convecting column
  !! CONV_ITERMAX, not this, is what ends the sweeps.
  real(kind=f), parameter :: CONV_TOL     = 1.e-9_f
  !! Cap on the adjustment's sweeps. Raised from the reference's 10 because
  !! ``pfact`` there is a constant, which makes an unstable pair repairable in
  !! a single pass, whereas here it depends on the temperatures being adjusted.
  !!
  !! Reaching an exactly neutral column takes far more sweeps than this -- the
  !! pairwise mixing relaxes the profile diffusively, so a sixty-layer column
  !! converges at the rate a Gauss-Seidel sweep does. Stopping short is by
  !! design and is how the reference uses the scheme too: the adjustment runs
  !! every timestep against a profile the previous step already left nearly
  !! neutral, so the residual relaxes across steps rather than within one call.
  !! Every sweep conserves enthalpy exactly, so a partial adjustment is a
  !! smaller correction, never a wrong one. ``resid`` reports what is left.
  integer, parameter      :: CONV_ITERMAX = 30

  !! Fewest layers a convective zone may span. Flags are set pairwise, so two
  !! is the shortest run the identification can produce, and a lone pair is the
  !! signature of a grid-scale oscillation rather than of convection. Three
  !! demands two adjacent neutral interfaces, which a 2*dz mode cannot fake.
  integer, parameter      :: CONV_MIN_LAYERS = 3

  !! Floor on the adiabatic gradient used by the adjustment. The Parmentier
  !! fit crosses zero at PARM_A/PARM_B = 9600 K and is negative above it,
  !! which is outside the domain it was fitted on and would make the mixing
  !! below drive the pair apart instead of onto an adiabat. The tabulated
  !! gradient is positive everywhere, so this binds only on the fit.
  real(kind=f), parameter :: GRAD_AD_MIN = 1.e-3_f

  !! Temperature perturbation used to measure the Jacobian [K].
  !!
  !! Large enough that the flux difference clears round-off in a sum over 1568
  !! spectral points, and small enough to stay inside one cell of the opacity
  !! table: the Sonora ck grid is ~50 K wide near 900 K and only piecewise
  !! linear, so a difference straddling a node measures the kink in the table
  !! rather than the derivative. A different table needs this re-chosen.
  real(kind=f), parameter :: JAC_DT = 2._f

  !! Levenberg damping added to the implicit system's diagonal when a solve
  !! fails, and the factor it decays by when one succeeds. Zero damping is the
  !! pure Newton step; large damping approaches the diagonal update.
  real(kind=f), parameter :: LAMBDA_MIN   = 1.e-3_f
  real(kind=f), parameter :: LAMBDA_DECAY = 0.5_f

  !! Solves between full Jacobian re-measurements, and the temperature drift
  !! that forces an early one. A Jacobian costs 2*nz flux evaluations, so it
  !! cannot ride on the solve trigger -- in a coupled run that fires every
  !! step.
  integer, parameter      :: JAC_GAP_DEFAULT = 200
  real(kind=f), parameter :: JAC_DRIFT_MAX   = 100._f

  !! Floor on the interval between re-measurements, whatever the drift says.
  !!
  !! Without it the drift trigger defeats the whole cadence. During spin-up a
  !! column far from equilibrium moves by tens of kelvin *per step* -- measured
  !! at 57 K/step on the 196-band case -- so any drift threshold worth having
  !! fires every step, and each firing costs 2*nz flux evaluations. That took
  !! the run from minutes to an estimated forty hours. The Jacobian being
  !! somewhat stale costs accuracy in the step direction, which the trust
  !! region and Levenberg damping already handle; rebuilding it every step
  !! costs the run.
  integer, parameter      :: JAC_GAP_MIN = 50

  !! Trust radius on the implicit step [K]: the largest |dT| any layer may take
  !! in one step before the whole vector is scaled down.
  !!
  !! This is *not* the same thing as ``dt_max``, and cannot be folded into it.
  !! ``dt_max`` is a user-facing backstop that a run may legitimately disable;
  !! the trust radius is what keeps a Newton step honest when the Jacobian is
  !! only an approximation to a nonlinear problem. Far from equilibrium the
  !! linear model happily predicts a step of thousands of kelvin, and taking it
  !! sends the column through zero temperature and straight to NaN -- which is
  !! exactly what a real 196-band column did within fifty steps when the two
  !! were conflated and the backstop was off.
  real(kind=f), parameter :: DT_TRUST = 20._f

  !! erg/g/K -> J/kg/K
  real(kind=f), parameter :: CP_CGS2SI = 1.e-4_f
  !! cm/s^2 -> m/s^2, and m -> cm
  real(kind=f), parameter :: ACC_CGS2SI = 1.e-2_f
  real(kind=f), parameter :: M2CM       = 1.e2_f

  !! Configuration and held state for the radiative-convective update.
  type rce_type

    ! ---- configuration ----------------------------------------------------
    integer      :: nz      = 0     !! layers
    integer      :: nbin    = 0     !! size bins
    integer      :: ngroup  = 0     !! particle groups
    integer      :: igridv  = 0     !! vertical grid type
    integer      :: mode    = I_RCE_EQUILIBRIUM

    real(kind=f) :: t_int   = 0._f  !! internal temperature [K]
    real(kind=f) :: cp      = 0._f  !! specific heat [J/kg/K]
    real(kind=f) :: grav    = 0._f  !! gravity [m/s^2]
    real(kind=f) :: wtmol   = 0._f  !! mean molecular weight [g/mol]

    ! ---- incident shortwave -----------------------------------------------
    !! Irradiation temperature [K]: the incident flux normal to the beam is
    !! ``sigma * t_irr**4``. Zero means an isolated object and switches the
    !! whole shortwave path off.
    real(kind=f) :: t_irr   = 0._f
    !! Temperature of the blackbody the incident spectrum is *shaped* like.
    !! It carries no energy scale -- that is t_irr's job -- so the two are
    !! independent and a hot star at weak irradiation is expressible.
    real(kind=f) :: t_star  = 0._f
    !! Cosine of the incidence angle. Also where redistribution lives: 0.5 is
    !! a dayside average, 0.25 full redistribution.
    real(kind=f) :: mu0     = 0.5_f
    !! Surface albedo. Zero for a gas giant, which has no surface.
    real(kind=f) :: w_surf  = 0._f
    !! True when t_irr > 0, i.e. when there is anything to solve.
    logical      :: sw_on   = .false.

    !! Per-band incident flux [W/m^2], summing to sigma*t_irr**4.
    real(kind=f), allocatable :: f0(:)
    !! Per-band Rayleigh attenuation factor: multiply by concentration
    !! [mol/m^3] and path [m] to get an optical depth.
    real(kind=f), allocatable :: ray_fac(:)
    !! Per-level incidence cosine, carried as an array so a sphericity
    !! correction can be added without changing the solver's interface.
    real(kind=f), allocatable :: mu_lev(:)
    !! Shortwave flux workspace.
    real(kind=f), allocatable :: sw_up(:), sw_dn(:)

    !! Diagnostics, refreshed on every solve.
    real(kind=f) :: absorbed_sw = 0._f  !! shortwave absorbed by the column [W/m^2]
    real(kind=f) :: reflected_sw = 0._f !! shortwave reflected to space [W/m^2]

    !! Which adiabat the convective interior follows.
    integer      :: adiabat = I_ADIABAT_PARMENTIER
    !! The tabulated gradient, only allocated for I_ADIABAT_TABLE.
    !! Axes are log10 T [K] and log10 p [bar]; grad is (nt, np).
    real(kind=f), allocatable :: ad_t(:)
    real(kind=f), allocatable :: ad_p(:)
    real(kind=f), allocatable :: ad_grad(:,:)

    real(kind=f) :: accel   = 1._f  !! heating acceleration factor
    real(kind=f) :: dt_max  = 0.5_f !! per-step |dT| clamp [K]
    real(kind=f) :: dt_tol  = 1._f  !! drift before re-solving [K]
    real(kind=f) :: dtau_tol = 0.02_f !! cloud cross-section change before re-solving
    integer      :: gap_max = 100   !! hard ceiling on steps between solves

    ! ---- the loaded opacity table -----------------------------------------
    type(ck_table_type) :: ck
    integer      :: nband = 0
    integer      :: ng    = 0
    integer      :: nwave = 0

    ! ---- held state -------------------------------------------------------
    !! Layer heating rate from the last solve [K/s]. Applied every step.
    real(kind=f), allocatable :: dtdt(:)
    !! Net upward flux at levels from the last solve [W/m^2].
    real(kind=f), allocatable :: fnet(:)
    !! Local radiative time constant per layer [s], used to stabilise the
    !! temperature update. See ``rce_update``.
    real(kind=f), allocatable :: tau_rad(:)
    !! Which layers the adjustment left convective. The primary description of
    !! the column's stability, since there can be more than one zone.
    logical, allocatable :: is_conv(:)
    !! Top layer of the convective zone reaching the base of the column, or 0
    !! when the deepest layer is not convective. A convenience derived from
    !! ``is_conv``; nothing branches on it.
    integer      :: nz_rcb = 0
    !! How many contiguous convective zones ``is_conv`` holds.
    integer      :: nzone  = 0
    !! Largest fractional superadiabaticity the last adjustment left behind.
    real(kind=f) :: conv_resid = 0._f

    integer      :: steps_since_solve = 0
    integer      :: nsolve            = 0    !! solves so far, for diagnostics
    integer      :: nclamp            = 0    !! how often the |dT| clamp fired
    !! Steps that fell back to the per-layer damping because the tridiagonal
    !! solve was singular or produced a non-finite step.
    integer      :: nfallback         = 0
    real(kind=f) :: dt_since_solve    = 0._f !! model time since last solve [s]
    real(kind=f) :: drift             = 0._f !! accumulated max|dT| since last solve [K]
    real(kind=f) :: xsec_ref         = -1._f !! cloud cross-section at last solve

    ! ---- work arrays, allocated once --------------------------------------
    real(kind=f), allocatable :: beta(:,:)      !! (nwave, nz) gas attenuation [1/m]
    real(kind=f), allocatable :: tau_c(:,:)     !! (nband, nz) cloud optical depth
    real(kind=f), allocatable :: w0_c(:,:)      !! (nband, nz)
    real(kind=f), allocatable :: g_c(:,:)       !! (nband, nz)
    real(kind=f), allocatable :: be(:,:)        !! (nband, nz+1) band Planck [W/m^2/sr]
    !! (nband, nz) band Planck at the layer centres [W/m^2/sr]. The levels
    !! above are interpolated from these, and an interpolation cannot see a
    !! layer-to-layer oscillation, so the solver is given both: the centres
    !! set how brightly each layer emits, the levels how that varies across it.
    real(kind=f), allocatable :: be_mid(:,:)
    real(kind=f), allocatable :: numden_grp(:,:,:) !! (nz, nbin, ngroup) [#/cm^3]
    real(kind=f), allocatable :: dtau(:), w0(:), gasym(:)
    real(kind=f), allocatable :: f_up(:), f_dn(:)
    real(kind=f), allocatable :: tl(:)          !! (nz+1) level temperature [K]
    real(kind=f), allocatable :: dz(:)          !! (nz) layer thickness [m]

    !! Planck-weighted optical depth, accumulated in the spectral loop and
    !! consumed by the radiative time constant.
    real(kind=f), allocatable :: tau_num(:), tau_den(:)

    !! Jacobian of the heating rate, ``jac(i,k) = d(dTdt_i)/dT_k`` [1/s].
    !!
    !! Dense, because it was measured to be. A tridiagonal truncation captures
    !! only 47-58% of the operator at every depth on this discretisation, and
    !! gets the sign of the terms it does keep backwards: perturbing one layer
    !! raises the levels it shares with its neighbours, so they emit more and
    !! cool, where genuine exchange would warm them. See ``rce_jacobian``.
    real(kind=f), allocatable :: jac(:,:)
    !! LU factors of the last implicit system, and its pivots.
    real(kind=f), allocatable :: lu(:,:)
    integer, allocatable      :: piv(:)
    !! Whether jac holds anything yet.
    logical :: have_jac = .false.
    !! The profile jac was measured at, for the drift trigger.
    real(kind=f), allocatable :: t_jac(:)
    !! Solves since jac was last measured, and the ceiling on that.
    integer :: jac_age = 0
    integer :: jac_gap = JAC_GAP_DEFAULT
    !! Levenberg damping on the implicit system, raised on a rejected step.
    real(kind=f) :: lambda = 0._f
    !! Scratch for the Jacobian's perturbed evaluations.
    real(kind=f), allocatable :: t_pert(:), fnet_pert(:), fnet_minus(:)
    real(kind=f), allocatable :: tnum_pert(:), tden_pert(:)

  end type rce_type

contains

  !! Set up the radiative-convective update and load the opacity table.
  !!
  !! ``cp_cgs`` and ``grav_cgs`` are taken in the units the ``physical_params``
  !! namelist uses (erg/g/K and cm/s^2) and converted here, so callers do not
  !! have to.
  subroutine rce_init(rce, nz, nbin, ngroup, nband, igridv, ck_path, &
                      t_int, t_irr, t_star, mu0, w_surf, &
                      cp_cgs, grav_cgs, wtmol, &
                      adiabat, adiabat_path, &
                      mode, accel, dt_max, dt_tol, dtau_tol, gap_max, rc)

    implicit none

    type(rce_type), intent(inout) :: rce
    integer, intent(in)      :: nz, nbin, ngroup, nband, igridv
    character(len=*), intent(in) :: ck_path
    integer, intent(in)      :: adiabat     !! I_ADIABAT_PARMENTIER or _TABLE
    character(len=*), intent(in) :: adiabat_path !! only read for _TABLE
    real(kind=f), intent(in) :: t_int       !! internal temperature [K]
    real(kind=f), intent(in) :: t_irr       !! irradiation temperature [K], 0 = isolated
    real(kind=f), intent(in) :: t_star      !! shape of the incident spectrum [K]
    real(kind=f), intent(in) :: mu0         !! cos(incidence), carries redistribution
    real(kind=f), intent(in) :: w_surf      !! surface albedo, 0 for a gas giant
    real(kind=f), intent(in) :: cp_cgs      !! [erg/g/K]
    real(kind=f), intent(in) :: grav_cgs    !! [cm/s^2]
    real(kind=f), intent(in) :: wtmol       !! [g/mol]
    integer, intent(in)      :: mode
    real(kind=f), intent(in) :: accel, dt_max, dt_tol, dtau_tol
    integer, intent(in)      :: gap_max
    integer, intent(inout)   :: rc

    rc = 0

    call ck_load(ck_path, rce%ck, rc)
    if (rc < 0) return

    ! Cloud optics are tabulated on band centres, so the two spectral grids
    ! have to be the same set of bands. A mismatch means the Mie tables were
    ! generated against a different ck table.
    if (rce%ck%nband /= nband) then
      write(*,*) 'rce_init::ERROR - ck table has', rce%ck%nband, &
                 'bands but the cloud optics were generated for', nband
      rc = -1
      return
    end if

    rce%nz     = nz
    rce%nbin   = nbin
    rce%ngroup = ngroup
    rce%igridv = igridv
    rce%nband  = rce%ck%nband
    rce%ng     = rce%ck%ng
    rce%nwave  = rce%ck%nwave

    rce%t_int  = t_int
    rce%t_irr  = t_irr
    rce%t_star = t_star
    rce%mu0    = mu0
    rce%w_surf = w_surf

    ! A shaped spectrum needs a temperature to be shaped by. Falling back on
    ! the irradiation temperature reproduces the common case where the beam is
    ! a blackbody carrying exactly its own flux.
    if (rce%t_star <= 0._f) rce%t_star = t_irr
    rce%cp     = cp_cgs   * CP_CGS2SI
    rce%grav   = grav_cgs * ACC_CGS2SI
    rce%wtmol  = wtmol

    rce%adiabat = adiabat
    if (rce%adiabat == I_ADIABAT_TABLE) then
      call rce_load_adiabat(adiabat_path, rce, rc)
      if (rc < 0) return
    end if

    rce%mode     = mode
    rce%accel    = accel
    rce%dt_max   = dt_max
    rce%dt_tol   = dt_tol
    rce%dtau_tol = dtau_tol
    rce%gap_max  = gap_max

    allocate(rce%dtdt(nz), rce%fnet(nz+1), rce%tau_rad(nz))
    allocate(rce%beta(rce%nwave, nz))
    allocate(rce%tau_c(rce%nband, nz), rce%w0_c(rce%nband, nz), &
             rce%g_c(rce%nband, nz))
    allocate(rce%be(rce%nband, nz+1))
    allocate(rce%be_mid(rce%nband, nz))
    allocate(rce%numden_grp(nz, nbin, ngroup))
    allocate(rce%dtau(nz), rce%w0(nz), rce%gasym(nz))
    allocate(rce%f_up(nz+1), rce%f_dn(nz+1))
    allocate(rce%tl(nz+1), rce%dz(nz))
    allocate(rce%tau_num(nz), rce%tau_den(nz))
    allocate(rce%jac(nz,nz), rce%lu(nz,nz), rce%piv(nz), rce%t_jac(nz))
    allocate(rce%t_pert(nz), rce%fnet_pert(nz+1), rce%fnet_minus(nz+1))
    allocate(rce%tnum_pert(nz), rce%tden_pert(nz))
    allocate(rce%is_conv(nz))

    allocate(rce%f0(rce%nband), rce%ray_fac(rce%nband))
    allocate(rce%mu_lev(nz+1), rce%sw_up(nz+1), rce%sw_dn(nz+1))

    rce%dtdt(:)    = 0._f
    rce%fnet(:)    = 0._f
    rce%tau_rad(:) = huge(0._f)
    rce%is_conv(:) = .false.

    call rce_init_shortwave(rce)

    return
  end subroutine rce_init


  !! Build the incident spectrum and the Rayleigh cross-sections, once.
  !!
  !! Both depend only on the band grid and on ``t_star`` / ``t_irr``, none of
  !! which move during a run, so this is startup work rather than per-solve
  !! work.
  subroutine rce_init_shortwave(rce)

    implicit none

    type(rce_type), intent(inout) :: rce

    integer      :: iband
    real(kind=f) :: total, lam_cm, wn_c, sig
    real(kind=f) :: x(RAY_NSPEC)

    rce%sw_on = (rce%t_irr > 0._f)

    rce%f0(:)      = 0._f
    rce%ray_fac(:) = 0._f
    rce%mu_lev(:)  = rce%mu0

    if (.not. rce%sw_on) return

    ! ---- incident spectrum -------------------------------------------------
    ! The shape is a blackbody at t_star; the magnitude is set by t_irr. They
    ! are separated because a companion's spectrum and the amount of it that
    ! arrives are independent facts, and folding them together would make a
    ! weakly irradiated hot primary inexpressible.
    do iband = 1, rce%nband
      rce%f0(iband) = bbflux_wavenumber(rce%ck%wmin(iband), &
                                        rce%ck%wmax(iband), rce%t_star)
    end do

    total = sum(rce%f0(:))

    if (total > 0._f) then
      rce%f0(:) = rce%f0(:) * (SIGMA_SB * rce%t_irr**4) / total
    else
      ! t_star so low that the band set captures no flux at all. Rather than
      ! divide by zero, put the whole budget in the longest-wavelength band,
      ! which is where it would have gone.
      rce%f0(rce%nband) = SIGMA_SB * rce%t_irr**4
    end if

    ! ---- Rayleigh ----------------------------------------------------------
    ! Evaluated at the band centre, matching how the cloud optics are
    ! tabulated. For a lambda^-4 law across a wide band that slightly
    ! underestimates the band mean; the bands are narrow enough that this is
    ! well inside the error already accepted in the cross-section itself.
    x(:) = RAY_SOLAR_H2HE(:)

    do iband = 1, rce%nband
      wn_c = 0.5_f * (rce%ck%wmin(iband) + rce%ck%wmax(iband))
      if (wn_c <= 0._f) cycle

      lam_cm = 1._f / wn_c
      sig    = ray_sigma_mix(x, lam_cm)

      ! Same conversion the gas opacity uses: cm^2/molecule -> 1/m, given a
      ! concentration in mol/m^3 supplied later.
      rce%ray_fac(iband) = 1.e-4_f * CK_AVOGADRO * sig
    end do

    return
  end subroutine rce_init_shortwave


  !! Release the opacity table and all work arrays.
  subroutine rce_destroy(rce)

    implicit none

    type(rce_type), intent(inout) :: rce

    call ck_destroy(rce%ck)

    if (allocated(rce%dtdt))       deallocate(rce%dtdt)
    if (allocated(rce%fnet))       deallocate(rce%fnet)
    if (allocated(rce%tau_rad))    deallocate(rce%tau_rad)
    if (allocated(rce%beta))       deallocate(rce%beta)
    if (allocated(rce%tau_c))      deallocate(rce%tau_c)
    if (allocated(rce%w0_c))       deallocate(rce%w0_c)
    if (allocated(rce%g_c))        deallocate(rce%g_c)
    if (allocated(rce%be))         deallocate(rce%be)
    if (allocated(rce%be_mid))     deallocate(rce%be_mid)
    if (allocated(rce%numden_grp)) deallocate(rce%numden_grp)
    if (allocated(rce%dtau))       deallocate(rce%dtau)
    if (allocated(rce%w0))         deallocate(rce%w0)
    if (allocated(rce%gasym))      deallocate(rce%gasym)
    if (allocated(rce%f_up))       deallocate(rce%f_up)
    if (allocated(rce%f_dn))       deallocate(rce%f_dn)
    if (allocated(rce%tl))         deallocate(rce%tl)
    if (allocated(rce%dz))         deallocate(rce%dz)
    if (allocated(rce%tau_num))    deallocate(rce%tau_num)
    if (allocated(rce%tau_den))    deallocate(rce%tau_den)
    if (allocated(rce%jac))        deallocate(rce%jac)
    if (allocated(rce%lu))         deallocate(rce%lu)
    if (allocated(rce%piv))        deallocate(rce%piv)
    if (allocated(rce%t_jac))      deallocate(rce%t_jac)
    if (allocated(rce%t_pert))     deallocate(rce%t_pert)
    if (allocated(rce%fnet_pert))  deallocate(rce%fnet_pert)
    if (allocated(rce%fnet_minus)) deallocate(rce%fnet_minus)
    if (allocated(rce%tnum_pert))  deallocate(rce%tnum_pert)
    if (allocated(rce%tden_pert))  deallocate(rce%tden_pert)
    if (allocated(rce%is_conv))    deallocate(rce%is_conv)
    if (allocated(rce%ad_t))       deallocate(rce%ad_t)
    if (allocated(rce%ad_p))       deallocate(rce%ad_p)
    if (allocated(rce%ad_grad))    deallocate(rce%ad_grad)
    if (allocated(rce%f0))         deallocate(rce%f0)
    if (allocated(rce%ray_fac))    deallocate(rce%ray_fac)
    if (allocated(rce%mu_lev))     deallocate(rce%mu_lev)
    if (allocated(rce%sw_up))      deallocate(rce%sw_up)
    if (allocated(rce%sw_dn))      deallocate(rce%sw_dn)

    return
  end subroutine rce_destroy


  !! The Parmentier (2015) adiabat, anchored at ``(t0, p0)``.
  !!
  !! Mirrors ``Carma.extend_atmosphere``'s ``new_T``, and is the closed-form
  !! counterpart of the gradient the convective adjustment uses, so a profile
  !! this module leaves neutrally stable is on the adiabat Python extended it
  !! along. At ``p == p0`` this returns ``t0`` exactly.
  pure function rce_adiabat(pres, t0, p0) result(temp)

    implicit none

    real(kind=f), intent(in) :: pres   !! pressure to evaluate at [Pa]
    real(kind=f), intent(in) :: t0     !! anchor temperature [K]
    real(kind=f), intent(in) :: p0     !! anchor pressure [Pa]
    real(kind=f)             :: temp

    real(kind=f) :: kappa

    kappa = t0 / (PARM_A - PARM_B * t0) * (pres / p0) ** PARM_A
    temp  = PARM_A * kappa / (1._f + PARM_B * kappa)

    return
  end function rce_adiabat


  !! Read the tabulated adiabatic gradient shipped with carmapy.
  !!
  !! The file is the one ``carmapy.adiabat`` reads, passed by path rather than
  !! copied into the run directory: '#' comment lines, then ``nt np``, the
  !! log10 T axis, the log10 p axis, and the gradient as nt rows of np values.
  !! A specific heat block follows it in the file and is not read.
  subroutine rce_load_adiabat(path, rce, rc)

    implicit none

    character(len=*), intent(in)  :: path
    type(rce_type), intent(inout) :: rce
    integer, intent(inout)        :: rc

    integer            :: lun, ios, nt, np, i, j
    character(len=512) :: line

    rc = 0

    open(newunit=lun, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'rce_load_adiabat::ERROR - cannot open adiabat table: ', &
                 trim(path)
      rc = -1
      return
    end if

    ! The dimensions are the first line that is neither blank nor a comment.
    do
      read(lun, '(A)', iostat=ios) line
      if (ios /= 0) then
        write(*,*) 'rce_load_adiabat::ERROR - no dimension line in: ', &
                   trim(path)
        close(lun)
        rc = -1
        return
      end if
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      exit
    end do

    read(line, *, iostat=ios) nt, np
    if (ios /= 0 .or. nt < 2 .or. np < 2) then
      write(*,*) 'rce_load_adiabat::ERROR - bad dimensions in: ', trim(path)
      close(lun)
      rc = -1
      return
    end if

    allocate(rce%ad_t(nt), rce%ad_p(np), rce%ad_grad(nt, np))

    read(lun, *, iostat=ios) (rce%ad_t(i), i = 1, nt)
    if (ios == 0) read(lun, *, iostat=ios) (rce%ad_p(j), j = 1, np)
    if (ios == 0) read(lun, *, iostat=ios) ((rce%ad_grad(i,j), j = 1, np), &
                                            i = 1, nt)
    close(lun)

    if (ios /= 0) then
      write(*,*) 'rce_load_adiabat::ERROR - truncated table: ', trim(path)
      rc = -1
      return
    end if

    return
  end subroutine rce_load_adiabat


  !! Tabulated adiabatic gradient d ln T / d ln p at ``(temp, pres)``.
  !!
  !! Bilinear in (log10 T, log10 p), clamped at the table edges -- the same
  !! treatment ``carmapy.adiabat.grad_ad`` gives it, so the two agree wherever
  !! they are both asked. Clamping is silent here; the Python side warns when
  !! it builds the profile, which is where a user can act on it.
  pure function rce_grad_ad(rce, temp, pres) result(grad)

    implicit none

    type(rce_type), intent(in) :: rce
    real(kind=f), intent(in)   :: temp   !! [K]
    real(kind=f), intent(in)   :: pres   !! [Pa]
    real(kind=f)               :: grad

    integer      :: it, ip, nt, np
    real(kind=f) :: lt, lp, ft, fp

    nt = size(rce%ad_t)
    np = size(rce%ad_p)

    lt = min(max(log10(temp), rce%ad_t(1)), rce%ad_t(nt))
    lp = min(max(log10(pres / BAR2PA), rce%ad_p(1)), rce%ad_p(np))

    it = 1
    do while (it < nt - 1 .and. rce%ad_t(it + 1) < lt)
      it = it + 1
    end do

    ip = 1
    do while (ip < np - 1 .and. rce%ad_p(ip + 1) < lp)
      ip = ip + 1
    end do

    ft = (lt - rce%ad_t(it)) / (rce%ad_t(it + 1) - rce%ad_t(it))
    fp = (lp - rce%ad_p(ip)) / (rce%ad_p(ip + 1) - rce%ad_p(ip))

    grad = (1._f - ft) * (1._f - fp) * rce%ad_grad(it,     ip    ) &
         +          ft  * (1._f - fp) * rce%ad_grad(it + 1, ip    ) &
         +          ft  *          fp  * rce%ad_grad(it + 1, ip + 1) &
         + (1._f - ft) *          fp  * rce%ad_grad(it,     ip + 1)

    return
  end function rce_grad_ad


  !! Integrate ln T along the tabulated adiabat, from ``x0`` to ``x1`` in
  !! ln p. RK4 with a capped step, so the result does not depend on how far
  !! apart the caller's levels happen to be.
  pure function rce_integrate_adiabat(rce, y0, x0, x1) result(y)

    implicit none

    type(rce_type), intent(in) :: rce
    real(kind=f), intent(in)   :: y0   !! ln T at x0
    real(kind=f), intent(in)   :: x0   !! ln p, start [Pa]
    real(kind=f), intent(in)   :: x1   !! ln p, end [Pa]
    real(kind=f)               :: y

    integer      :: nstep, i
    real(kind=f) :: h, x, k1, k2, k3, k4

    nstep = max(1, ceiling(abs(x1 - x0) / ADIABAT_DLNP))
    h     = (x1 - x0) / real(nstep, f)

    y = y0
    x = x0

    do i = 1, nstep
      k1 = rce_grad_ad(rce, exp(y),                    exp(x))
      k2 = rce_grad_ad(rce, exp(y + 0.5_f * h * k1),   exp(x + 0.5_f * h))
      k3 = rce_grad_ad(rce, exp(y + 0.5_f * h * k2),   exp(x + 0.5_f * h))
      k4 = rce_grad_ad(rce, exp(y + h * k3),           exp(x + h))

      y = y + h * (k1 + 2._f * k2 + 2._f * k3 + k4) / 6._f
      x = x + h
    end do

    return
  end function rce_integrate_adiabat


  !! The adiabatic temperature ratio across the pair of layers ``(iz, iz+1)``.
  !!
  !! ``exp(integral of grad_ad dlnp)`` across the pair, evaluated by the
  !! midpoint rule: one gradient at the pair's mean temperature and the
  !! geometric mean of its pressures. Exact where ``grad_ad`` is constant over
  !! a layer spacing, second order otherwise, and it costs one interpolation
  !! rather than an RK4 sweep -- which matters because the adjustment runs on
  !! every timestep, not only on radiative solves.
  !!
  !! Layer ``iz`` is the deeper of the two, so the result is >= 1.
  pure function rce_pfact(rce, tbar, p_lo, p_hi) result(pfact)

    implicit none

    type(rce_type), intent(in) :: rce
    real(kind=f), intent(in)   :: tbar   !! mean temperature of the pair [K]
    real(kind=f), intent(in)   :: p_lo   !! pressure of the deeper layer [Pa]
    real(kind=f), intent(in)   :: p_hi   !! pressure of the shallower layer [Pa]
    real(kind=f)               :: pfact

    real(kind=f) :: grad

    if (rce%adiabat == I_ADIABAT_TABLE) then
      grad = rce_grad_ad(rce, tbar, sqrt(p_lo * p_hi))
    else
      grad = PARM_A - PARM_B * tbar
    end if

    pfact = (p_lo / p_hi) ** max(grad, GRAD_AD_MIN)

    return
  end function rce_pfact


  !! Dry convective adjustment: mix away any superadiabatic gradient.
  !!
  !! After Ray Pierrehumbert's scheme as implemented in Elspeth Lee's
  !! ``Exo-FMS_column_ck`` (``src/dry_conv_adj_mod.f90``, ``Ray_dry_adj``),
  !! generalised from that code's constant ``kappa`` to the run's own adiabat
  !! so the H2-dissociation flattening the tabulated gradient carries is not
  !! thrown away in the deep atmosphere.
  !!
  !! A pair of adjacent layers is unstable when the deeper one is hotter than
  !! the adiabat through the shallower one allows. An unstable pair is mixed at
  !! **fixed enthalpy** and placed exactly on the adiabat: the two assignments
  !! below solve ``dp*T`` summed over the pair being unchanged together with
  !! ``t(iz) = pfact*t(iz+1)``, which holds for any ``pfact`` and so survives
  !! ``pfact`` depending on the temperatures being adjusted. Sweeping down and
  !! then up propagates that through the column.
  !!
  !! Arrays are bottom-to-top: layer ``iz`` sits below layer ``iz+1``, and
  !! ``p(1)`` is the deepest.
  pure subroutine rce_conv_adj(rce, p, pl, t, is_conv, nz_rcb, nzone, resid)

    implicit none

    type(rce_type), intent(in)  :: rce
    real(kind=f), intent(in)    :: p(rce%nz)      !! layer pressure [Pa]
    real(kind=f), intent(in)    :: pl(rce%nz+1)   !! level pressure [Pa]
    real(kind=f), intent(inout) :: t(rce%nz)      !! layer temperature [K]
    logical, intent(out)        :: is_conv(rce%nz) !! left neutrally stable
    integer, intent(out)        :: nz_rcb  !! top of the zone reaching the base
    integer, intent(out)        :: nzone   !! contiguous convective zones
    !! Largest fractional superadiabaticity left in the column, 0 if none.
    real(kind=f), intent(out)   :: resid

    integer      :: nz, iz, iter
    logical      :: did_adj
    logical      :: neutral(rce%nz)
    real(kind=f) :: dp(rce%nz), tbar, pfact

    nz = rce%nz

    do iz = 1, nz
      dp(iz) = pl(iz) - pl(iz+1)
    end do

    do iter = 1, CONV_ITERMAX
      did_adj = .false.

      do iz = nz - 1, 1, -1
        tbar  = (dp(iz) * t(iz) + dp(iz+1) * t(iz+1)) / (dp(iz) + dp(iz+1))
        pfact = rce_pfact(rce, tbar, p(iz), p(iz+1))

        if (t(iz) > pfact * t(iz+1) * (1._f + CONV_TOL)) then
          t(iz+1) = (dp(iz) + dp(iz+1)) * tbar / (dp(iz+1) + pfact * dp(iz))
          t(iz)   = pfact * t(iz+1)
          did_adj = .true.
        end if
      end do

      do iz = 1, nz - 1
        tbar  = (dp(iz) * t(iz) + dp(iz+1) * t(iz+1)) / (dp(iz) + dp(iz+1))
        pfact = rce_pfact(rce, tbar, p(iz), p(iz+1))

        if (t(iz) > pfact * t(iz+1) * (1._f + CONV_TOL)) then
          t(iz+1) = (dp(iz) + dp(iz+1)) * tbar / (dp(iz+1) + pfact * dp(iz))
          t(iz)   = pfact * t(iz+1)
          did_adj = .true.
        end if
      end do

      if (.not. did_adj) exit
    end do

    ! A neutrally stable zone no longer adjusts, so it cannot be recognised
    ! from the sweeps. Both layers of every interface sitting on or above the
    ! adiabat belong to a convective zone -- "above" because a pair the sweeps
    ! did not have the budget to finish is still a convecting one.
    is_conv(:) = .false.
    neutral(:) = .false.
    resid      = 0._f

    do iz = 1, nz - 1
      tbar  = (dp(iz) * t(iz) + dp(iz+1) * t(iz+1)) / (dp(iz) + dp(iz+1))
      pfact = rce_pfact(rce, tbar, p(iz), p(iz+1))

      if (t(iz) >= pfact * t(iz+1) * (1._f - CONV_TOL)) then
        neutral(iz)   = .true.
        is_conv(iz)   = .true.
        is_conv(iz+1) = .true.
      end if

      resid = max(resid, t(iz) / (pfact * t(iz+1)) - 1._f)
    end do

    ! ---- discard zones too thin to be convection --------------------------
    ! Flags are set an interface at a time, so the shortest run expressible is
    ! two layers -- which is also exactly what a grid-scale (2*dz) oscillation
    ! in t manufactures: one interface made momentarily neutral by an
    ! anomalously warm layer, with both its neighbours left strongly
    ! subadiabatic. Measured on a real column: an isolated pair at 0.005 bar
    ! sat at S = -4e-6 between interfaces at S = -6e-2, while genuine zones ran
    ! 13 and 14 layers deep.
    !
    ! Keeping such a pair is not harmless bookkeeping. A layer marked
    ! convective is taken out of the implicit temperature solve and given the
    ! per-layer diagonal update instead, so the flag switches off the one
    ! mechanism able to damp a coupled mode at the very layers carrying it, and
    ! the oscillation that produced the flag is then sustained by it.
    !
    ! What joins two layers into one zone is therefore the *interface* between
    ! them, not the fact that both carry a flag. A sawtooth flags every layer
    ! it spans -- each one sits on a neutral interface with one neighbour --
    ! but the interfaces alternate neutral and strongly stable, so by
    ! interfaces it is a chain of two-layer islands and the purge below sees
    ! it for what it is. Reading contiguity off the layer flags instead welds
    ! those islands into a single deep zone that clears any length test.
    call rce_purge_thin_zones(nz, neutral, is_conv)

    nzone = 0
    do iz = 1, nz
      if (.not. is_conv(iz)) cycle
      if (iz == 1) then
        nzone = nzone + 1
      else if (.not. is_conv(iz-1)) then
        nzone = nzone + 1
      else if (.not. neutral(iz-1)) then
        nzone = nzone + 1
      end if
    end do

    nz_rcb = 0
    do iz = 1, nz
      if (.not. is_conv(iz)) exit
      nz_rcb = iz
      if (iz < nz) then
        if (.not. neutral(iz)) exit
      end if
    end do

    return
  end subroutine rce_conv_adj


  !! Clear the convective flag from any run shorter than CONV_MIN_LAYERS.
  !!
  !! A run extends from one layer to the next only across a neutral interface,
  !! ``neutral(iz)`` describing the interface between layers ``iz`` and
  !! ``iz+1``. Two flagged layers with a stably stratified interface between
  !! them belong to different zones, however adjacent they are: convection
  !! does not cross a stable interface, and a run that appears to is a
  !! grid-scale oscillation wearing the flags of the neutral interfaces it
  !! manufactured.
  !!
  !! Only the flags move. The temperatures are left exactly as the adjustment
  !! sweeps left them, which is the right division of labour: a pair the sweeps
  !! did not touch was never superadiabatic enough to mix, and unflagging it
  !! only returns it to the radiative solve that should have been handling it.
  pure subroutine rce_purge_thin_zones(nz, neutral, is_conv)

    implicit none

    integer, intent(in)    :: nz
    logical, intent(in)    :: neutral(nz)
    logical, intent(inout) :: is_conv(nz)

    integer :: iz, ibeg

    iz = 1
    do while (iz <= nz)

      if (.not. is_conv(iz)) then
        iz = iz + 1
        cycle
      end if

      ibeg = iz
      do while (iz < nz)
        if (.not. is_conv(iz+1)) exit
        if (.not. neutral(iz)) exit
        iz = iz + 1
      end do

      if (iz - ibeg + 1 < CONV_MIN_LAYERS) is_conv(ibeg:iz) = .false.

      iz = iz + 1

    end do

    return
  end subroutine rce_purge_thin_zones


  !! Reconstruct level temperatures from layer-centre values.
  !!
  !! Linear in ``ln p``, with the top and bottom levels extrapolated from the
  !! two nearest layers. Both arrays are bottom-to-top.
  pure subroutine rce_level_temps(nz, p, pl, t, tl)

    implicit none

    integer, intent(in)       :: nz
    real(kind=f), intent(in)  :: p(nz)     !! layer pressure [Pa]
    real(kind=f), intent(in)  :: pl(nz+1)  !! level pressure [Pa]
    real(kind=f), intent(in)  :: t(nz)     !! layer temperature [K]
    real(kind=f), intent(out) :: tl(nz+1)  !! level temperature [K]

    integer      :: iz
    real(kind=f) :: w

    if (nz == 1) then
      tl(:) = t(1)
      return
    end if

    do iz = 2, nz
      w = (log(pl(iz)) - log(p(iz))) / (log(p(iz-1)) - log(p(iz)))
      tl(iz) = t(iz) + w * (t(iz-1) - t(iz))
    end do

    ! Bottom and top levels lie outside the layer centres, so extrapolate on
    ! the same slope as the adjacent pair.
    w = (log(pl(1)) - log(p(1))) / (log(p(1)) - log(p(2)))
    tl(1) = t(1) + w * (t(1) - t(2))

    w = (log(pl(nz+1)) - log(p(nz))) / (log(p(nz)) - log(p(nz-1)))
    tl(nz+1) = t(nz) + w * (t(nz) - t(nz-1))

    return
  end subroutine rce_level_temps


  !! Recompute the altitude grid hydrostatically for the current temperature.
  !!
  !! Ports ``Carma.calculate_z``. For ``I_LOGP`` the grid is a function of
  !! pressure and the bottom scale height only, so it does not move as
  !! temperature evolves and this routine leaves it untouched.
  subroutine rce_regrid_z(nz, igridv, pl, tl, wtmol, grav, zc, zl)

    implicit none

    integer, intent(in)      :: nz
    integer, intent(in)      :: igridv
    real(kind=f), intent(in) :: pl(nz+1)   !! level pressure [Pa]
    real(kind=f), intent(in) :: tl(nz+1)   !! level temperature [K]
    real(kind=f), intent(in) :: wtmol      !! mean molecular weight [g/mol]
    real(kind=f), intent(in) :: grav       !! gravity [m/s^2]

    real(kind=f), intent(inout) :: zc(nz)    !! layer centre altitude [m]
    real(kind=f), intent(inout) :: zl(nz+1)  !! level altitude [m]

    integer      :: iz
    real(kind=f) :: scale_h

    if (igridv /= I_CART) return

    zl(1) = 0._f
    do iz = 2, nz + 1
      scale_h = KB_SI * tl(iz) / (wtmol * MPROT_SI * grav)
      zl(iz)  = zl(iz-1) + scale_h * log(pl(iz-1) / pl(iz))
    end do

    do iz = 1, nz
      zc(iz) = 0.5_f * (zl(iz) + zl(iz+1))
    end do

    return
  end subroutine rce_regrid_z


  !! Total cloud geometric cross-section per unit area of the column,
  !! ``sum_z sum_bin sum_group numden * pi * r^2 * dz``.
  !!
  !! A cheap proxy for how much the cloud field's opacity has moved since the
  !! last solve. Absolute normalisation is irrelevant -- only the relative
  !! change is used.
  pure function cloud_xsec(nz, nbin, ngroup, numden, radius, dz) result(xsec)

    implicit none

    integer, intent(in)      :: nz, nbin, ngroup
    real(kind=f), intent(in) :: numden(nz, nbin, ngroup)  !! [#/cm^3]
    real(kind=f), intent(in) :: radius(nbin, ngroup)      !! [cm]
    real(kind=f), intent(in) :: dz(nz)                    !! [cm]
    real(kind=f)             :: xsec

    integer      :: iz, ibin, igroup
    real(kind=f) :: acc

    xsec = 0._f
    do igroup = 1, ngroup
      do ibin = 1, nbin
        acc = 0._f
        do iz = 1, nz
          acc = acc + numden(iz, ibin, igroup) * dz(iz)
        end do
        xsec = xsec + acc * radius(ibin, igroup) ** 2
      end do
    end do

    return
  end function cloud_xsec


  !! Full radiative transfer solve: fills ``rce%fnet`` and ``rce%dtdt``.
  !!
  !! Gas and cloud optical properties are combined per spectral point, the
  !! two-stream problem is solved at each of the ``nband*ng`` points, and the
  !! results are summed with the correlated-k weights. The heating rate then
  !! follows from the flux divergence across each layer.
  subroutine rce_solve(rce, p, pl, t, radius, qext, ssa, asym)

    implicit none

    type(rce_type), intent(inout) :: rce
    real(kind=f), intent(in) :: p(rce%nz)      !! layer pressure [Pa]
    real(kind=f), intent(in) :: pl(rce%nz+1)   !! level pressure [Pa]
    real(kind=f), intent(in) :: t(rce%nz)      !! layer temperature [K]
    real(kind=f), intent(in) :: radius(rce%nbin, rce%ngroup)  !! [cm]
    real(kind=f), intent(in) :: qext(rce%nband, rce%nbin, rce%ngroup)
    real(kind=f), intent(in) :: ssa(rce%nband, rce%nbin, rce%ngroup)
    real(kind=f), intent(in) :: asym(rce%nband, rce%nbin, rce%ngroup)

    integer      :: nz, iz
    real(kind=f) :: dmass, tau_mean, emiss

    nz = rce%nz

    call rce_fluxes(rce, p, pl, t, radius, qext, ssa, asym, &
                    rce%fnet, rce%tau_num, rce%tau_den)

    ! The internal flux is a boundary condition, not something to be inferred:
    ! the net flux through the base of the column is sigma*t_int^4 by
    ! definition. toon_lw_column's gas-giant lower boundary extrapolates it
    ! from the deep Planck gradient instead, which is only right once the deep
    ! profile already is. Imposing it makes F_TOA == sigma*t_int^4 the fixed
    ! point of the column's energy budget, which is what the convective
    ! adjustment then transports the flux to satisfy.
    rce%fnet(1) = SIGMA_SB * rce%t_int ** 4

    ! ---- flux divergence -> heating rate ----------------------------------
    do iz = 1, nz
      ! The layer's mass per unit area. Hydrostatic balance makes this exact,
      ! where rho*dz is only exact in the limit of thin layers -- they differ
      ! by (dlnp)^2/24, which is 0.7% at a layer spacing of 0.4 in ln p. It
      ! has to be the exact one, because the convective adjustment redistributes
      ! enthalpy weighted by dp and any mismatch between the two measures
      ! shows up as the column gaining or losing energy in the handoff.
      dmass = (pl(iz) - pl(iz+1)) / rce%grav

      rce%dtdt(iz) = -(rce%fnet(iz+1) - rce%fnet(iz)) / (dmass * rce%cp)

      ! Local radiative time constant: the layer's heat capacity per unit area
      ! over how fast its own emission responds to a temperature change. Used
      ! only to stabilise the update in rce_update, never as a physical
      ! result, so the effective emissivity is the usual grey estimate.
      if (rce%tau_den(iz) > 0._f) then
        tau_mean = rce%tau_num(iz) / rce%tau_den(iz)
      else
        tau_mean = 0._f
      end if
      emiss = max(1._f - exp(-tau_mean), 1.e-30_f)

      rce%tau_rad(iz) = dmass * rce%cp &
                        / (4._f * SIGMA_SB * t(iz) ** 3 * emiss)
    end do

    rce%nsolve = rce%nsolve + 1

    return

    rce%nsolve = rce%nsolve + 1

    return
  end subroutine rce_solve


  !! Net flux at every level for a given temperature profile.
  !!
  !! Split out of ``rce_solve`` so ``rce_jacobian`` can evaluate it on a
  !! perturbed profile without disturbing any held state. Everything that
  !! depends on temperature is recomputed here from ``t``: the level
  !! temperatures the Planck source is built from, and the gas opacity.
  !! ``dz`` deliberately is not -- the altitude grid is regridded once per
  !! step, after the update, so within a step it is a constant and a
  !! temperature perturbation must not move it.
  subroutine rce_fluxes(rce, p, pl, t, radius, qext, ssa, asym, &
                        fnet, tau_num, tau_den)

    implicit none

    type(rce_type), intent(inout) :: rce
    real(kind=f), intent(in) :: p(rce%nz)      !! layer pressure [Pa]
    real(kind=f), intent(in) :: pl(rce%nz+1)   !! level pressure [Pa]
    real(kind=f), intent(in) :: t(rce%nz)      !! layer temperature [K]
    real(kind=f), intent(in) :: radius(rce%nbin, rce%ngroup)  !! [cm]
    real(kind=f), intent(in) :: qext(rce%nband, rce%nbin, rce%ngroup)
    real(kind=f), intent(in) :: ssa(rce%nband, rce%nbin, rce%ngroup)
    real(kind=f), intent(in) :: asym(rce%nband, rce%nbin, rce%ngroup)

    real(kind=f), intent(out) :: fnet(rce%nz+1)   !! net upward flux [W/m^2]
    real(kind=f), intent(out) :: tau_num(rce%nz)  !! Planck-weighted tau, numerator
    real(kind=f), intent(out) :: tau_den(rce%nz)  !! ... and its denominator

    integer      :: nz, nlev, iz, iband, ig, iw
    real(kind=f) :: conc(rce%nz), dz_cm(rce%nz), tl(rce%nz+1)
    !! (nband, nz) the part of each layer's emission its bounding levels
    !! cannot see: the layer-centre Planck less the level mean.
    real(kind=f) :: be_corr(rce%nband, rce%nz)
    real(kind=f) :: btop_factor, tau_gas, tau_tot, wt, be_lay
    real(kind=f) :: tau_ray, tau_cld, tau_sca
    real(kind=f) :: sw_top, sw_bot, sw_ref

    call rce_level_temps(rce%nz, p, pl, t, tl)

    nz   = rce%nz
    nlev = nz + 1

    ! ---- gas attenuation --------------------------------------------------
    ! The ck tables are per molecule of the gas mixture, so the concentration
    ! is that of the whole atmosphere.
    do iz = 1, nz
      conc(iz)  = p(iz) / (RGAS_SI * t(iz))
      dz_cm(iz) = rce%dz(iz) * M2CM
    end do

    call ck_kappa_column(rce%ck, p, t, conc, rce%beta)

    ! ---- cloud optics on band centres -------------------------------------
    call cloud_optics_column(nz, rce%nbin, rce%nband, rce%ngroup, &
                             rce%numden_grp, radius, dz_cm, &
                             qext, ssa, asym, &
                             rce%tau_c, rce%w0_c, rce%g_c)

    ! ---- band-integrated Planck at every level, and every layer -----------
    ! The layer values are what each layer emits from; the level values give
    ! the solver the gradient across it. Both are needed because `tl` is an
    ! interpolation of `t` and so is blind to a layer-to-layer oscillation.
    do iband = 1, rce%nband
      do iz = 1, nlev
        rce%be(iband, iz) = bbflux_wavenumber(rce%ck%wmin(iband), &
                                              rce%ck%wmax(iband), tl(iz))
      end do

      call bbflux_wavenumber_col(rce%ck%wmin(iband), rce%ck%wmax(iband), &
                                 t(1:nz), rce%be_mid(iband, 1:nz))

      do iz = 1, nz
        be_corr(iband, iz) = rce%be_mid(iband, iz) &
                             - 0.5_f * (rce%be(iband, iz) + rce%be(iband, iz+1))
      end do
    end do

    ! No incident radiation: the top boundary is the auto-emission branch with
    ! the factor that reproduces PICASO's, in top-down pressure terms.
    btop_factor = pl(nlev) / (pl(nlev-1) - pl(nlev))

    ! ---- solve every spectral point ---------------------------------------
    fnet(:)    = 0._f
    tau_num(:)  = 0._f
    tau_den(:)  = 0._f

    sw_top = 0._f
    sw_bot = 0._f
    sw_ref = 0._f

    do iband = 1, rce%nband
      do ig = 1, rce%ng
        iw = (iband - 1) * rce%ng + ig

        do iz = 1, nz
          tau_gas = rce%beta(iw, iz) * rce%dz(iz)

          ! Rayleigh scattering by the gas. Only assembled when there is a
          ! beam to scatter: it falls as lambda^-4, so in the thermal infrared
          ! it is five orders of magnitude below the optical and contributes
          ! nothing, while leaving it out of an unirradiated run keeps such a
          ! run bit-for-bit what it was before the shortwave existed.
          if (rce%sw_on) then
            tau_ray = rce%ray_fac(iband) * conc(iz) * rce%dz(iz)
          else
            tau_ray = 0._f
          end if

          tau_cld = rce%tau_c(iband, iz)
          tau_tot = tau_gas + tau_cld + tau_ray

          rce%dtau(iz) = tau_tot

          ! Scattering optical depth: the cloud's scattering part plus all of
          ! the Rayleigh, which is conservative.
          tau_sca = tau_cld * rce%w0_c(iband, iz) + tau_ray

          if (tau_tot > 0._f) then
            ! Gas absorption does not scatter, so the scatterers are the cloud
            ! and the Rayleigh term.
            rce%w0(iz) = tau_sca / tau_tot
          else
            rce%w0(iz) = 0._f
          end if
          ! The two-stream coefficients are singular at w0 == 1 exactly, which
          ! a non-absorbing grain in a transparent band can reach.
          rce%w0(iz) = min(rce%w0(iz), 1._f - 1.e-12_f)

          ! Asymmetry is a scattering-weighted mean. Rayleigh is isotropic and
          ! so enters with g = 0, dragging the column's asymmetry down wherever
          ! it competes with the forward-peaked cloud.
          if (tau_sca > 0._f) then
            rce%gasym(iz) = tau_cld * rce%w0_c(iband, iz) &
                            * rce%g_c(iband, iz) / tau_sca
          else
            rce%gasym(iz) = 0._f
          end if

          ! Planck-weighted mean optical depth, accumulated here so the
          ! stabilisation timescale below costs nothing extra.
          be_lay = rce%be_mid(iband, iz)
          wt = rce%ck%weights(iw)
          tau_num(iz) = tau_num(iz) + wt * be_lay * tau_tot
          tau_den(iz) = tau_den(iz) + wt * be_lay
        end do

        call toon_lw_column(nz, rce%dtau, rce%w0, rce%gasym, &
                            rce%be(iband, :), 0._f, -1, btop_factor, &
                            .false., .false., rce%f_up, rce%f_dn, &
                            be_corr_in=be_corr(iband, :))

        wt = rce%ck%weights(iw)
        do iz = 1, nlev
          fnet(iz) = fnet(iz) + wt * (rce%f_up(iz) - rce%f_dn(iz))
        end do

        ! ---- the incident beam, into the same net flux --------------------
        ! Accumulating shortwave into `fnet` rather than tracking it apart is
        ! what keeps the rest of this module unchanged. In equilibrium the
        ! longwave carries sigma*t_int^4 + absorbed upward at the top while the
        ! shortwave net there is -absorbed, so the two cancel and
        ! `fnet(top) == sigma*t_int^4` remains the condition to converge on,
        ! exactly as for an isolated object. The heating rate picks up
        ! shortwave absorption for free, being a divergence of this same array.
        if (rce%sw_on) then
          if (rce%f0(iband) > 0._f) then
            call toon_sw_column(nz, rce%f0(iband), rce%mu_lev, &
                                rce%dtau, rce%w0, rce%gasym, rce%w_surf, &
                                rce%sw_up, rce%sw_dn)

            do iz = 1, nlev
              fnet(iz) = fnet(iz) + wt * (rce%sw_up(iz) - rce%sw_dn(iz))
            end do

            sw_top = sw_top + wt * (rce%sw_dn(nlev) - rce%sw_up(nlev))
            sw_bot = sw_bot + wt * (rce%sw_dn(1) - rce%sw_up(1))
            sw_ref = sw_ref + wt * rce%sw_up(nlev)
          end if
        end if
      end do
    end do

    ! What the column kept, and what it sent back. Both are diagnostics: the
    ! energy is already in `fnet`.
    rce%absorbed_sw  = sw_top - sw_bot
    rce%reflected_sw = sw_ref

    return
  end subroutine rce_fluxes


  !! Jacobian of the heating rate, ``jac(i,k) = d(dTdt_i)/dT_k`` [1/s].
  !!
  !! A per-layer damping can only see the diagonal -- how a layer responds to
  !! its own temperature. That makes the update a Jacobi sweep: every layer
  !! jumps to the equilibrium it would have if its neighbours held still, and
  !! since they all move at once they overshoot together. The off-diagonals are
  !! the coupling that makes them not independent.
  !!
  !! **Measured, not derived.** Differentiating the two-stream solution
  !! analytically means committing to a linearisation of the whole chain --
  !! opacity, Planck source, boundary terms -- and a wrong one is not obviously
  !! wrong. Finite differences on the real solver cannot fail in that way.
  !!
  !! **Dense, because it was measured to be.** The fraction of
  !! ``sum_j |dH_j/dT_k|`` captured by a band of half-width w on this column:
  !!
  !!     p [bar]    w=1     w=2     w=4     w=8    w=16
  !!       22.5    0.517   0.999   1.000   1.000   1.000
  !!       0.111   0.518   0.670   0.790   0.884   0.939
  !!       0.0058  0.528   0.669   0.773   0.860   0.936
  !!
  !! A tridiagonal keeps about half the operator at every depth and never
  !! saturates. It also gets the sign of what it keeps backwards: perturbing a
  !! layer raises the levels it shares with its neighbours, so their Planck
  !! source rises and they *cool*, where genuine exchange would warm them. That
  !! is a property of this layer-to-level reconstruction, not of radiative
  !! transfer -- a different vertical scheme needs the measurement redone. See
  !! ``tests/integration/test_rce.py::test_jacobian_bandwidth``.
  !!
  !! Cost is ``2*nz`` flux evaluations, so this runs on its own slow cadence
  !! rather than on every solve; see ``rce_update``.
  subroutine rce_jacobian(rce, p, pl, t, radius, qext, ssa, asym)

    implicit none

    type(rce_type), intent(inout) :: rce
    real(kind=f), intent(in) :: p(rce%nz)
    real(kind=f), intent(in) :: pl(rce%nz+1)
    real(kind=f), intent(in) :: t(rce%nz)
    real(kind=f), intent(in) :: radius(rce%nbin, rce%ngroup)
    real(kind=f), intent(in) :: qext(rce%nband, rce%nbin, rce%ngroup)
    real(kind=f), intent(in) :: ssa(rce%nband, rce%nbin, rce%ngroup)
    real(kind=f), intent(in) :: asym(rce%nband, rce%nbin, rce%ngroup)

    integer      :: nz, iz, k
    real(kind=f) :: dmass, hp, hm

    nz = rce%nz

    do k = 1, nz

      ! Centred, not one-sided: the heating goes as T^4, so a one-sided
      ! difference carries a first-order bias that a centred one cancels.
      rce%t_pert(:) = t(:)
      rce%t_pert(k) = t(k) + JAC_DT
      call rce_fluxes(rce, p, pl, rce%t_pert, radius, qext, ssa, asym, &
                      rce%fnet_pert, rce%tnum_pert, rce%tden_pert)

      rce%t_pert(k) = t(k) - JAC_DT
      call rce_fluxes(rce, p, pl, rce%t_pert, radius, qext, ssa, asym, &
                      rce%fnet_minus, rce%tnum_pert, rce%tden_pert)

      ! The base flux is imposed, not solved for, so it does not respond to a
      ! temperature perturbation and neither does layer 1's heating through it.
      rce%fnet_pert(1)  = rce%fnet(1)
      rce%fnet_minus(1) = rce%fnet(1)

      do iz = 1, nz
        dmass = (pl(iz) - pl(iz+1)) / rce%grav

        hp = -(rce%fnet_pert(iz+1)  - rce%fnet_pert(iz))  / (dmass * rce%cp)
        hm = -(rce%fnet_minus(iz+1) - rce%fnet_minus(iz)) / (dmass * rce%cp)

        rce%jac(iz, k) = (hp - hm) / (2._f * JAC_DT)
      end do

    end do

    rce%t_jac(:)  = t(:)
    rce%jac_age   = 0
    rce%have_jac  = .true.

    return
  end subroutine rce_jacobian


  !! Thomas algorithm for a tridiagonal system, with the sub- and
  !! super-diagonals given as full-length arrays (``a(1)`` and ``c(n)`` unused).
  !!
  !! Kept local rather than shared with ``carma_rtsolve``'s copy so the
  !! validated solver is not touched.
  pure subroutine rce_tridiag(n, a, b, c, d, x, ok)

    implicit none

    integer, intent(in)       :: n
    real(kind=f), intent(in)  :: a(n), b(n), c(n), d(n)
    real(kind=f), intent(out) :: x(n)
    logical, intent(out)      :: ok

    integer      :: i
    real(kind=f) :: cp(n), dp(n), denom

    ok = .true.

    if (b(1) == 0._f) then
      ok = .false.
      x(:) = 0._f
      return
    end if

    cp(1) = c(1) / b(1)
    dp(1) = d(1) / b(1)

    do i = 2, n
      denom = b(i) - a(i) * cp(i-1)
      if (denom == 0._f) then
        ok = .false.
        x(:) = 0._f
        return
      end if
      cp(i) = c(i) / denom
      dp(i) = (d(i) - a(i) * dp(i-1)) / denom
    end do

    x(n) = dp(n)
    do i = n-1, 1, -1
      x(i) = dp(i) - cp(i) * x(i+1)
    end do

    return
  end subroutine rce_tridiag



  !! Advance the column's temperature and altitude grid by one timestep.
  !!
  !! Called once per microphysical step, before ``CARMASTATE_Create``. A full
  !! radiative transfer solve happens only when one of the cadence triggers
  !! fires; the heating rate from the most recent solve is applied on every
  !! call, so temperature evolves continuously. The convective adjustment runs
  !! on every call too -- it is cheap, and leaving a superadiabatic profile
  !! standing between solves would feed the microphysics a column that is not
  !! one the model believes in.
  !!
  !! ``numden`` is indexed by element, as the driver holds it. Only the number
  !! element of each group carries a number concentration -- core-mass elements
  !! are a mass concentration in the same array -- so ``elem_is_number``
  !! selects the right one.
  subroutine rce_update(rce, dtime, p, pl, t, zc, zl, &
                        nelem, numden, elem2group, elem_is_number, &
                        radius, qext, ssa, asym)

    implicit none

    type(rce_type), intent(inout) :: rce
    real(kind=f), intent(in) :: dtime          !! microphysical timestep [s]
    real(kind=f), intent(in) :: p(rce%nz)      !! layer pressure [Pa]
    real(kind=f), intent(in) :: pl(rce%nz+1)   !! level pressure [Pa]

    real(kind=f), intent(inout) :: t(rce%nz)     !! layer temperature [K]
    real(kind=f), intent(inout) :: zc(rce%nz)    !! layer centre altitude [m]
    real(kind=f), intent(inout) :: zl(rce%nz+1)  !! level altitude [m]

    integer, intent(in)      :: nelem
    real(kind=f), intent(in) :: numden(rce%nz, nelem, rce%nbin)  !! [#/cm^3]
    integer, intent(in)      :: elem2group(nelem)
    logical, intent(in)      :: elem_is_number(nelem)

    real(kind=f), intent(in) :: radius(rce%nbin, rce%ngroup)  !! [cm]
    real(kind=f), intent(in) :: qext(rce%nband, rce%nbin, rce%ngroup)
    real(kind=f), intent(in) :: ssa(rce%nband, rce%nbin, rce%ngroup)
    real(kind=f), intent(in) :: asym(rce%nband, rce%nbin, rce%ngroup)

    integer      :: nz, iz, iz_top, ibin, ielem, igroup
    real(kind=f) :: dz_cm(rce%nz), damp(rce%nz)
    real(kind=f) :: xsec, dt_eff, dt0, dt_lay
    real(kind=f) :: dt_col(rce%nz), rhs(rce%nz), sol(rce%nz), dt_cap
    integer      :: rad_idx(rce%nz), nrad, i, j
    logical      :: dense_ok, do_jac
    logical      :: do_solve

    nz = rce%nz

    ! ---- geometry and level temperatures for the current profile -----------
    do iz = 1, nz
      rce%dz(iz) = zl(iz+1) - zl(iz)
      dz_cm(iz)  = rce%dz(iz) * M2CM
    end do

    call rce_level_temps(nz, p, pl, t, rce%tl)

    ! ---- repack number density from elements to groups --------------------
    rce%numden_grp(:,:,:) = 0._f
    do ielem = 1, nelem
      if (.not. elem_is_number(ielem)) cycle
      igroup = elem2group(ielem)
      do ibin = 1, rce%nbin
        do iz = 1, nz
          rce%numden_grp(iz, ibin, igroup) = numden(iz, ielem, ibin)
        end do
      end do
    end do

    ! ---- decide whether to re-solve ---------------------------------------
    xsec = cloud_xsec(nz, rce%nbin, rce%ngroup, rce%numden_grp, radius, dz_cm)

    do_solve = (rce%nsolve == 0)
    if (.not. do_solve) do_solve = rce%steps_since_solve >= rce%gap_max
    if (.not. do_solve) do_solve = rce%drift >= rce%dt_tol
    if (.not. do_solve .and. rce%xsec_ref > 0._f) then
      do_solve = abs(xsec - rce%xsec_ref) > rce%dtau_tol * rce%xsec_ref
    end if

    if (do_solve) then
      call rce_solve(rce, p, pl, t, radius, qext, ssa, asym)

      ! The Jacobian costs 2*nz flux evaluations against the solve's one, so it
      ! gets its own, much slower trigger: it describes how the column responds
      ! to temperature, which changes far more slowly than the heating rate
      ! itself. Re-measured when it has never been measured, when it is stale by
      ! count, or when the profile has moved far enough that it no longer
      ! describes this column.
      do_jac = .not. rce%have_jac
      if (.not. do_jac) do_jac = rce%jac_age >= rce%jac_gap
      if (.not. do_jac .and. rce%jac_age >= JAC_GAP_MIN) then
        do_jac = maxval(abs(t - rce%t_jac)) > JAC_DRIFT_MAX
      end if

      if (do_jac) then
        call rce_jacobian(rce, p, pl, t, radius, qext, ssa, asym)
      else
        rce%jac_age = rce%jac_age + 1
      end if

      rce%steps_since_solve = 0
      rce%dt_since_solve    = 0._f
      rce%drift             = 0._f
      rce%xsec_ref          = xsec
    end if

    ! ---- apply the held heating rate to every layer -----------------------
    if (rce%mode == I_RCE_EQUILIBRIUM) then
      dt_eff = dtime * rce%accel
    else
      dt_eff = dtime
    end if

    ! Backward Euler on the coupled column:
    !
    !   (C_i/dt) dT_i - sum_k J_ik dT_k = C_i * dtdt_i
    !
    ! with J the dense Jacobian of the heating rate. Dividing through by the
    ! layer's heat capacity leaves the system below in units of 1/s.
    !
    ! With the off-diagonals zeroed this reduces algebraically to the per-layer
    ! damping it replaces, dT = dtdt*dt_eff/(1 + dt_eff/tau_rad), since the
    ! diagonal of J is what tau_rad estimates. The off-diagonals are the whole
    ! point: without them each layer relaxes to its own equilibrium holding its
    ! neighbours fixed, which is a Jacobi sweep and oscillates.
    !
    ! As dt_eff grows this becomes dT = -J^-1 dtdt, which is Newton's method
    ! for the steady state. accel is therefore the continuation parameter
    ! between a near-explicit step and a full Newton step, not a separate
    ! acceleration hack.
    !
    ! In a radiative layer the fixed point is untouched -- dT = 0 requires
    ! dtdt = 0 whatever J is -- so this changes only the approach to it.
    !
    ! The system is solved over the radiative layers only, rows and columns
    ! both. A convective zone is handled separately below, so its rows must not
    ! appear in the right-hand side here: the coupling would carry the zone's
    ! dtdt into the radiative layers beside it, and overwriting the zone's own
    ! rows afterwards would leave that behind -- measured at 2.6e-4 in the
    ! emergent flux, and it does not shrink with more steps. Dropping the zone's
    ! rows and columns treats it as a fixed boundary, which is what handling it
    ! separately means.
    dense_ok  = rce%have_jac
    dt_col(:) = 0._f

    if (dense_ok) then
      nrad = 0
      do iz = 1, nz
        if (.not. rce%is_conv(iz)) then
          nrad = nrad + 1
          rad_idx(nrad) = iz
        end if
      end do

      if (nrad > 0) then
        do j = 1, nrad
          do i = 1, nrad
            rce%lu(i,j) = -dt_eff * rce%jac(rad_idx(i), rad_idx(j))
          end do
          ! Levenberg damping. The measured diagonal is genuinely positive in
          ! the upper atmosphere, where a layer carrying far more flux than it
          ! emits absorbs more of the through-beam as it warms, so it is not
          ! sign-clipped -- lambda is what keeps such a row solvable.
          rce%lu(j,j) = rce%lu(j,j) + 1._f + rce%lambda
          rhs(j)      = dt_eff * rce%dtdt(rad_idx(j))
        end do

        call lu_factor(nrad, rce%lu(1:nrad,1:nrad), rce%piv(1:nrad), dense_ok)

        if (dense_ok) then
          call lu_solve(nrad, rce%lu(1:nrad,1:nrad), rce%piv(1:nrad), &
                        rhs(1:nrad), sol(1:nrad))

          ! A measured Jacobian can be poor where the column is optically thin
          ! enough that a layer barely talks to anything.
          if (any(sol(1:nrad) /= sol(1:nrad))) dense_ok = .false.
        end if

        if (dense_ok) then
          do j = 1, nrad
            dt_col(rad_idx(j)) = sol(j)
          end do
        end if
      end if
    end if

    if (.not. dense_ok) then
      do iz = 1, nz
        dt_col(iz) = rce%dtdt(iz) * dt_eff &
                     / (1._f + dt_eff / rce%tau_rad(iz))
      end do
      rce%nfallback = rce%nfallback + 1
      ! Raise the damping so the next attempt is closer to the diagonal step
      ! that just had to be substituted anyway.
      rce%lambda = max(2._f * rce%lambda, LAMBDA_MIN)
    else
      rce%lambda = rce%lambda * LAMBDA_DECAY
      if (rce%lambda < LAMBDA_MIN) rce%lambda = 0._f
    end if

    ! Convective layers keep the per-layer damping they had before, made
    ! uniform within each zone; only radiative layers take the implicit step.
    !
    ! Not conservatism: the implicit step cannot be rescaled to sit inside a
    ! zone. A radiative layer reaches equilibrium on its own, at dtdt == 0, so
    ! any positive factor leaves that fixed point alone. A convective zone does
    ! not -- it equilibrates collectively, when the flux entering its base
    ! equals the flux leaving its top, with individual layers still heating and
    ! cooling against the mixing. That balance is a dp-weighted sum over the
    ! zone, so a factor varying from layer to layer reweights it and moves the
    ! equilibrium rather than only the approach to it: measured at 2.6e-4 in
    ! the emergent flux, and it does not shrink with more steps.
    !
    ! Collapsing the zone onto a single factor is what fixes that, and it needs
    ! the step to be a positive multiple of dtdt. The coupled step is not --
    ! the off-diagonals can give a layer a dT opposing its own dtdt, which is
    ! the whole point of them -- so there is no factor to take. The zone gets
    ! the diagonal update instead, which is what it was validated with.
    !
    ! The oscillation the implicit step exists to cure lives in the radiative
    ! upper atmosphere, where the layers are optically thin and exchange
    ! strongly, so this gives it up exactly where it was not needed.
    do iz = 1, nz
      damp(iz) = 1._f / (1._f + dt_eff / rce%tau_rad(iz))
    end do

    iz = 1
    do while (iz <= nz)
      if (.not. rce%is_conv(iz)) then
        iz = iz + 1
        cycle
      end if

      iz_top = iz
      do while (iz_top < nz)
        if (.not. rce%is_conv(iz_top + 1)) exit
        iz_top = iz_top + 1
      end do

      damp(iz:iz_top) = minval(damp(iz:iz_top))
      iz = iz_top + 1
    end do

    do iz = 1, nz
      if (rce%is_conv(iz)) then
        dt_col(iz) = rce%dtdt(iz) * dt_eff * damp(iz)
      end if
    end do

    ! Trust region: if the step is too big, scale the *whole vector* rather
    ! than clipping layer by layer. A Newton step is a direction; clipping
    ! individual entries changes that direction into one nothing computed, and
    ! silently discards the difference -- which is how the old per-layer clamp
    ! came to destroy 823 kW/m^2 while reporting a converged-looking flux. A
    ! uniform scale is still a descent step, just a shorter one.
    !
    ! The radius is the tighter of the trust radius, which always applies, and
    ! the user's dt_max, which a run may disable. The linear model is only
    ! local: far from equilibrium it will confidently predict a step of
    ! thousands of kelvin.
    dt_lay = 0._f
    do iz = 1, nz
      dt_lay = max(dt_lay, abs(dt_col(iz)))
    end do

    dt_cap = min(DT_TRUST, rce%dt_max)

    if (dt_lay > dt_cap .and. dt_lay > 0._f) then
      dt_col(:) = dt_col(:) * (dt_cap / dt_lay)
      dt_lay    = dt_cap
      rce%nclamp = rce%nclamp + 1
      ! A step that had to be cut back means the linear model overreached, so
      ! lean towards the diagonal next time.
      rce%lambda = max(2._f * rce%lambda, LAMBDA_MIN)
    end if

    do iz = 1, nz
      t(iz) = t(iz) + dt_col(iz)
    end do

    ! ---- mix away whatever the heating left superadiabatic ----------------
    call rce_conv_adj(rce, p, pl, t, rce%is_conv, rce%nz_rcb, rce%nzone, &
                      rce%conv_resid)

    ! ---- move the altitude grid to match the new temperature --------------
    call rce_level_temps(nz, p, pl, t, rce%tl)
    call rce_regrid_z(nz, rce%igridv, pl, rce%tl, rce%wtmol, rce%grav, zc, zl)

    rce%steps_since_solve = rce%steps_since_solve + 1
    rce%dt_since_solve    = rce%dt_since_solve + dtime
    rce%drift             = rce%drift + dt_lay

    return
  end subroutine rce_update

end module carma_rce
