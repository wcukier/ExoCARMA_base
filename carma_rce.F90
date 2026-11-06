!! Radiative-convective temperature evolution for a brown dwarf column.
!!
!! Each layer is heated by the radiative flux divergence,
!! ``dT/dt = -dF_net/dz / (rho*Cp)``, and a dry convective adjustment then
!! mixes every superadiabatic pair of layers at fixed enthalpy onto the
!! adiabat selected by ``adiabat`` (the Parmentier 2015 fit or the tabulated
!! H/He gradient). Convective zones are diagnosed per layer in ``is_conv``;
!! there may be more than one.
!!
!! The net flux through the base of the column is imposed as sigma*t_int**4,
!! so equilibrium is ``F_TOA == sigma*t_int**4``. An incident beam
!! (``t_irr > 0``) is added into the same net flux.
!!
!! A full solve runs on an adaptive cadence, triggered by accumulated
!! temperature drift, a change in the cloud cross-section, or a ceiling on
!! steps since the last solve; the held ``dTdt`` is applied on every step.
!! ``rce_kzz`` optionally recomputes the mixing-length Kzz from the same solve.
!!
!! Units are SI (Pa, m, K, W/m^2), except that Cp [erg/g/K] and gravity
!! [cm/s^2] arrive in cgs from the namelist, and ``cloud_optics_column`` and
!! ``rce_kzz`` work in cgs.
!!
!! Level temperatures are reconstructed from the layer-centre ``t`` by linear
!! interpolation in ``ln p``, extrapolated at the two ends.
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
  public :: rce_grad_ad, rce_grad, rce_integrate_adiabat, rce_load_adiabat
  public :: rce_conv_adj, rce_kzz
  public :: I_KZZ_STATIC, I_KZZ_MIXING_LENGTH
  public :: rce_fluxes, rce_jacobian, rce_solve, rce_optics_prep
  public :: I_ADIABAT_PARMENTIER, I_ADIABAT_TABLE

  !! Hold the input Kzz profile for the whole run.
  integer, parameter :: I_KZZ_STATIC        = 0
  !! Recompute Kzz from mixing length theory on each solve.
  integer, parameter :: I_KZZ_MIXING_LENGTH = 1

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

  !! Largest RK4 step in ln p when integrating the tabulated adiabat. The
  !! table spacing is 0.46 in ln p.
  real(kind=f), parameter :: ADIABAT_DLNP = 0.05_f

  !! 1 bar in Pa. The table's pressure axis is log10 bar.
  real(kind=f), parameter :: BAR2PA = 1.e5_f

  !! Relative tolerance on the convective stability test: a pair this close to
  !! the adiabat counts as neutral. Relative because an absolute 1e-6 K (as in
  !! Exo-FMS) is near roundoff at a 3000 K interior.
  real(kind=f), parameter :: CONV_TOL     = 1.e-9_f
  !! Cap on the adjustment's sweeps per call. Each sweep conserves enthalpy, so
  !! stopping early leaves a residual that relaxes over later steps; ``resid``
  !! reports it.
  integer, parameter      :: CONV_ITERMAX = 30

  !! Fewest layers a convective zone may span. A two-layer zone is the
  !! signature of a grid-scale (2*dz) oscillation rather than convection.
  integer, parameter      :: CONV_MIN_LAYERS = 3

  !! Floor on the adiabatic gradient. The Parmentier fit is negative above
  !! PARM_A/PARM_B = 9600 K, outside its fitted domain.
  real(kind=f), parameter :: GRAD_AD_MIN = 1.e-3_f

  !! ---- mixing length eddy diffusion (PICASO get_kzz), see rce_kzz -------
  !! Scale factor on the mixing length velocity.
  real(kind=f), parameter :: KZZ_SCALEF   = 1._f / 3._f
  !! Floor on the mixing length / scale height; keeps Kzz finite in stable
  !! layers.
  real(kind=f), parameter :: KZZ_MIXL_MIN = 0.1_f
  !! Overshoot: the convective flux may fall with height by at most this
  !! fraction of the pressure ratio across a pair of layers. Arbitrary.
  real(kind=f), parameter :: KZZ_CHF_FALL = 1._f / 3._f
  !! Floor on the convective flux, sigma*(KZZ_TEFF_MIN*t_int)**4.
  real(kind=f), parameter :: KZZ_TEFF_MIN = 0.05_f
  !! Cp/R used by the Kzz formula (diatomic ideal gas). Not rce%cp: the
  !! parameterisation was calibrated with this value.
  real(kind=f), parameter :: KZZ_CP_OVER_R = 7._f / 2._f
  !! Universal gas constant [ erg / mol / K ], PICASO's value.
  real(kind=f), parameter :: RGAS_CGS = 8.3143e7_f
  !! Stefan-Boltzmann constant [ erg / cm^2 / s / K^4 ].
  real(kind=f), parameter :: SIGMA_SB_CGS = 0.56687e-4_f
  !! W/m^2 -> erg/cm^2/s
  real(kind=f), parameter :: FLUX_SI2CGS = 1.e3_f
  !! Pa -> barye
  real(kind=f), parameter :: PA2BARYE = 1.e1_f

  !! Temperature perturbation used to measure the Jacobian [K]. Large enough to
  !! clear roundoff in the spectral sum, small enough to stay within one cell
  !! of the opacity table (the Sonora ck grid is ~50 K wide near 900 K).
  real(kind=f), parameter :: JAC_DT = 2._f

  !! Levenberg damping added to the implicit system's diagonal when a solve
  !! fails, and the factor it decays by when one succeeds.
  real(kind=f), parameter :: LAMBDA_MIN   = 1.e-3_f
  real(kind=f), parameter :: LAMBDA_DECAY = 0.5_f

  !! Solves between Jacobian re-measurements (each costs 2*nz flux
  !! evaluations), and the temperature drift [K] that forces an early one.
  integer, parameter      :: JAC_GAP         = 200
  real(kind=f), parameter :: JAC_DRIFT_MAX   = 100._f

  !! Fewest solves between re-measurements, so the drift trigger cannot fire
  !! on every solve during spin-up.
  integer, parameter      :: JAC_GAP_MIN = 50

  !! Trust radius on the implicit step [K]: the largest |dT| any layer may take
  !! before the whole step is scaled down. Far from equilibrium the linear
  !! model predicts steps of thousands of kelvin.
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

    real(kind=f) :: t_int   = 0._f  !! internal temperature [K]
    real(kind=f) :: cp      = 0._f  !! specific heat [J/kg/K]
    real(kind=f) :: grav    = 0._f  !! gravity [m/s^2]
    real(kind=f) :: wtmol   = 0._f  !! mean molecular weight [g/mol]

    ! ---- incident shortwave -----------------------------------------------
    !! Irradiation temperature [K]: the incident flux normal to the beam is
    !! ``sigma * t_irr**4``. Zero switches the shortwave off.
    real(kind=f) :: t_irr   = 0._f
    !! Blackbody temperature setting the shape of the incident spectrum [K];
    !! t_irr sets its magnitude.
    real(kind=f) :: t_star  = 0._f
    !! Cosine of the incidence angle, which also carries redistribution: 0.5 is
    !! a dayside average, 0.25 full redistribution.
    real(kind=f) :: mu0     = 0.5_f
    !! Surface albedo; zero for a gas giant.
    real(kind=f) :: w_surf  = 0._f
    !! True when t_irr > 0.
    logical      :: sw_on   = .false.

    !! Per-band incident flux [W/m^2], summing to sigma*t_irr**4.
    real(kind=f), allocatable :: f0(:)
    !! Per-band Rayleigh attenuation factor: multiply by concentration
    !! [mol/m^3] and path [m] to get an optical depth.
    real(kind=f), allocatable :: ray_fac(:)
    !! Per-level incidence cosine.
    real(kind=f), allocatable :: mu_lev(:)

    !! Diagnostics, refreshed on every solve.
    real(kind=f) :: absorbed_sw = 0._f  !! shortwave absorbed by the column [W/m^2]
    real(kind=f) :: reflected_sw = 0._f !! shortwave reflected to space [W/m^2]

    !! I_KZZ_STATIC or I_KZZ_MIXING_LENGTH.
    integer      :: kzz_mode = I_KZZ_STATIC

    !! Multiplier on the mixing length; 1 is PICASO's parameterisation. Kzz
    !! goes as L^(4/3), so this scales the profile by its 4/3 power.
    real(kind=f) :: kzz_mixl_scale = 1._f

    !! Which adiabat the convective interior follows.
    integer      :: adiabat = I_ADIABAT_PARMENTIER
    !! The tabulated gradient, only allocated for I_ADIABAT_TABLE.
    !! Axes are log10 T [K] and log10 p [bar]; grad is (nt, np).
    real(kind=f), allocatable :: ad_t(:)
    real(kind=f), allocatable :: ad_p(:)
    real(kind=f), allocatable :: ad_grad(:,:)

    !! Heating acceleration factor on the radiative timestep. 1 keeps the
    !! radiative and microphysical clocks together.
    real(kind=f) :: accel   = 1._f
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
    !! Eddy diffusion at levels from the last solve [cm^2/s]. Only refreshed
    !! for I_KZZ_MIXING_LENGTH.
    real(kind=f), allocatable :: kzz(:)
    !! Factor the convective flux was scaled by to carry the internal flux at
    !! the base of the column; 1 in equilibrium.
    real(kind=f) :: kzz_rescale = 1._f
    !! Local radiative time constant per layer [s], used to damp the
    !! temperature update.
    real(kind=f), allocatable :: tau_rad(:)
    !! Which layers the adjustment left convective.
    logical, allocatable :: is_conv(:)
    !! Top layer of the convective zone reaching the base of the column, or 0
    !! when the deepest layer is not convective.
    integer      :: nz_rcb = 0
    !! How many contiguous convective zones ``is_conv`` holds.
    integer      :: nzone  = 0
    !! Largest fractional superadiabaticity the last adjustment left behind.
    real(kind=f) :: conv_resid = 0._f

    integer      :: steps_since_solve = 0
    integer      :: nsolve            = 0    !! solves so far, for diagnostics
    integer      :: nclamp            = 0    !! how often the |dT| clamp fired
    !! Steps that fell back to the per-layer damping because the implicit
    !! solve was singular or produced a non-finite step.
    integer      :: nfallback         = 0
    real(kind=f) :: dt_since_solve    = 0._f !! model time since last solve [s]
    real(kind=f) :: drift             = 0._f !! accumulated max|dT| since last solve [K]
    real(kind=f) :: xsec_ref         = -1._f !! cloud cross-section at last solve
    !! Cumulative energy budget of the update [J/m^2]: heating computed by the
    !! solves (sum of C*dtdt*dt), heating applied to t, and enthalpy change of
    !! the convective adjustment.
    !! LLM-GENERATED — UNVERIFIED (tier: physics)
    real(kind=f) :: e_comp = 0._f
    real(kind=f) :: e_appl = 0._f
    real(kind=f) :: e_adj  = 0._f
    !! Computed and applied, per layer [J/m^2].
    real(kind=f), allocatable :: e_comp_z(:), e_appl_z(:)
    !! Computed and applied, split into radiative and convective layers [J/m^2].
    real(kind=f) :: e_appl_rad = 0._f, e_comp_rad = 0._f
    real(kind=f) :: e_appl_cnv = 0._f, e_comp_cnv = 0._f

    ! ---- work arrays, allocated once --------------------------------------
    real(kind=f), allocatable :: beta(:,:)      !! (nwave, nz) gas attenuation [1/m]
    real(kind=f), allocatable :: tau_c(:,:)     !! (nband, nz) cloud optical depth
    real(kind=f), allocatable :: w0_c(:,:)      !! (nband, nz)
    real(kind=f), allocatable :: g_c(:,:)       !! (nband, nz)
    real(kind=f), allocatable :: be(:,:)        !! (nz+1, nband) band Planck [W/m^2/sr]
    !! (nz, nband) band Planck at the layer centres [W/m^2/sr]. Passed to the
    !! solver alongside the levels, which are interpolated and cannot see a
    !! layer-to-layer oscillation.
    real(kind=f), allocatable :: be_mid(:,:)
    real(kind=f), allocatable :: numden_grp(:,:,:) !! (nz, nbin, ngroup) [#/cm^3]
    !! Per-spectral-point results, summed in a fixed order after the threaded
    !! spectral loop so the total does not depend on the thread layout. The
    !! operands are stored rather than wt*d, so each accumulation is one
    !! fused multiply-add.
    real(kind=f), allocatable :: dfl_w(:,:)      !! (nz+1, nwave) f_up - f_dn
    real(kind=f), allocatable :: dfs_w(:,:)      !! (nz+1, nwave) sw_up - sw_dn
    real(kind=f), allocatable :: tautot_w(:,:)   !! (nz, nwave)
    real(kind=f), allocatable :: swt_w(:), swb_w(:), swr_w(:)
    logical, allocatable      :: sw_hit(:)       !! did this point carry a beam
    real(kind=f), allocatable :: tl(:)          !! (nz+1) level temperature [K]
    real(kind=f), allocatable :: dz(:)          !! (nz) layer thickness [m]

    !! Planck-weighted optical depth, accumulated in the spectral loop and
    !! consumed by the radiative time constant.
    real(kind=f), allocatable :: tau_num(:), tau_den(:)

    !! Dense Jacobian of the heating rate, ``jac(i,k) = d(dTdt_i)/dT_k`` [1/s].
    !! See ``rce_jacobian``.
    real(kind=f), allocatable :: jac(:,:)
    !! LU factors of the last implicit system, and its pivots.
    real(kind=f), allocatable :: lu(:,:)
    integer, allocatable      :: piv(:)
    !! Whether jac holds anything yet.
    logical :: have_jac = .false.
    !! The profile jac was measured at, for the drift trigger.
    real(kind=f), allocatable :: t_jac(:)
    !! Solves since jac was last measured.
    integer :: jac_age = 0
    !! Levenberg damping on the implicit system, raised on a rejected step.
    real(kind=f) :: lambda = 0._f
    !! Scratch for the Jacobian's perturbed evaluations.
    real(kind=f), allocatable :: t_pert(:), fnet_pert(:), fnet_minus(:)
    real(kind=f), allocatable :: tnum_pert(:), tden_pert(:)

  end type rce_type

contains

  !! Set up the radiative-convective update and load the opacity table.
  !! ``cp_cgs`` and ``grav_cgs`` are in the namelist's cgs units.
  subroutine rce_init(rce, nz, nbin, ngroup, nband, igridv, ck_path, &
                      t_int, t_irr, t_star, mu0, w_surf, &
                      cp_cgs, grav_cgs, wtmol, &
                      adiabat, adiabat_path, kzz_mode, kzz_mixl_scale, &
                      accel, dt_tol, dtau_tol, gap_max, rc)

    implicit none

    type(rce_type), intent(inout) :: rce
    integer, intent(in)      :: nz, nbin, ngroup, nband, igridv
    character(len=*), intent(in) :: ck_path
    integer, intent(in)      :: adiabat     !! I_ADIABAT_PARMENTIER or _TABLE
    integer, intent(in)      :: kzz_mode    !! I_KZZ_STATIC or _MIXING_LENGTH
    real(kind=f), intent(in) :: kzz_mixl_scale !! multiplier on the mixing length
    character(len=*), intent(in) :: adiabat_path !! only read for _TABLE
    real(kind=f), intent(in) :: t_int       !! internal temperature [K]
    real(kind=f), intent(in) :: t_irr       !! irradiation temperature [K], 0 = isolated
    real(kind=f), intent(in) :: t_star      !! shape of the incident spectrum [K]
    real(kind=f), intent(in) :: mu0         !! cos(incidence), carries redistribution
    real(kind=f), intent(in) :: w_surf      !! surface albedo, 0 for a gas giant
    real(kind=f), intent(in) :: cp_cgs      !! [erg/g/K]
    real(kind=f), intent(in) :: grav_cgs    !! [cm/s^2]
    real(kind=f), intent(in) :: wtmol       !! [g/mol]
    real(kind=f), intent(in) :: accel, dt_tol, dtau_tol
    integer, intent(in)      :: gap_max
    integer, intent(inout)   :: rc

    rc = 0

    call ck_load(ck_path, rce%ck, rc)
    if (rc < 0) return

    ! The cloud optics must be on the ck table's bands.
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

    ! Default to a blackbody beam at the irradiation temperature.
    if (rce%t_star <= 0._f) rce%t_star = t_irr

    rce%cp     = cp_cgs   * CP_CGS2SI
    rce%grav   = grav_cgs * ACC_CGS2SI
    rce%wtmol  = wtmol

    rce%adiabat = adiabat
    if (rce%adiabat == I_ADIABAT_TABLE) then
      call rce_load_adiabat(adiabat_path, rce, rc)
      if (rc < 0) return
    end if

    rce%kzz_mode = kzz_mode

    ! Checked here because rce_kzz is pure and cannot report an error.
    if (kzz_mixl_scale <= 0._f) then
      write(*,*) 'rce_init::ERROR - kzz_mixl_scale must be positive, got', &
                 kzz_mixl_scale
      rc = -1
      return
    end if
    rce%kzz_mixl_scale = kzz_mixl_scale

    rce%accel    = accel
    rce%dt_tol   = dt_tol
    rce%dtau_tol = dtau_tol
    rce%gap_max  = gap_max

    allocate(rce%dtdt(nz), rce%fnet(nz+1), rce%kzz(nz+1), rce%tau_rad(nz))
    allocate(rce%beta(rce%nwave, nz))
    allocate(rce%tau_c(rce%nband, nz), rce%w0_c(rce%nband, nz), &
             rce%g_c(rce%nband, nz))
    allocate(rce%be(nz+1, rce%nband))
    allocate(rce%be_mid(nz, rce%nband))
    allocate(rce%numden_grp(nz, nbin, ngroup))
    allocate(rce%dfl_w(nz+1, rce%nwave), rce%dfs_w(nz+1, rce%nwave))
    allocate(rce%tautot_w(nz, rce%nwave))
    allocate(rce%swt_w(rce%nwave), rce%swb_w(rce%nwave), rce%swr_w(rce%nwave))
    allocate(rce%sw_hit(rce%nwave))
    allocate(rce%tl(nz+1), rce%dz(nz))
    allocate(rce%tau_num(nz), rce%tau_den(nz))
    allocate(rce%jac(nz,nz), rce%lu(nz,nz), rce%piv(nz), rce%t_jac(nz))
    allocate(rce%t_pert(nz), rce%fnet_pert(nz+1), rce%fnet_minus(nz+1))
    allocate(rce%tnum_pert(nz), rce%tden_pert(nz))
    allocate(rce%is_conv(nz))

    allocate(rce%f0(rce%nband), rce%ray_fac(rce%nband))
    allocate(rce%mu_lev(nz+1))

    rce%dtdt(:)    = 0._f
    rce%fnet(:)    = 0._f
    rce%kzz(:)     = 0._f
    rce%tau_rad(:) = huge(0._f)
    rce%is_conv(:) = .false.

    call rce_init_shortwave(rce)

    return
  end subroutine rce_init


  !! Build the incident spectrum and the Rayleigh cross-sections, which are
  !! fixed for the run.
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
    ! Shape is a blackbody at t_star, magnitude sigma*t_irr**4.
    do iband = 1, rce%nband
      rce%f0(iband) = bbflux_wavenumber(rce%ck%wmin(iband), &
                                        rce%ck%wmax(iband), rce%t_star)
    end do

    total = sum(rce%f0(:))

    if (total > 0._f) then
      rce%f0(:) = rce%f0(:) * (SIGMA_SB * rce%t_irr**4) / total
    else
      ! t_star too low for the bands to capture any flux: put it all in the
      ! longest-wavelength band.
      rce%f0(rce%nband) = SIGMA_SB * rce%t_irr**4
    end if

    ! ---- Rayleigh ----------------------------------------------------------
    ! Evaluated at the band centre, like the cloud optics.
    x(:) = RAY_SOLAR_H2HE(:)

    do iband = 1, rce%nband
      wn_c = 0.5_f * (rce%ck%wmin(iband) + rce%ck%wmax(iband))
      if (wn_c <= 0._f) cycle

      lam_cm = 1._f / wn_c
      sig    = ray_sigma_mix(x, lam_cm)

      ! cm^2/molecule -> m^2/mol, so that times mol/m^3 gives 1/m.
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
    if (allocated(rce%kzz))        deallocate(rce%kzz)
    if (allocated(rce%tau_rad))    deallocate(rce%tau_rad)
    if (allocated(rce%beta))       deallocate(rce%beta)
    if (allocated(rce%tau_c))      deallocate(rce%tau_c)
    if (allocated(rce%w0_c))       deallocate(rce%w0_c)
    if (allocated(rce%g_c))        deallocate(rce%g_c)
    if (allocated(rce%be))         deallocate(rce%be)
    if (allocated(rce%be_mid))     deallocate(rce%be_mid)
    if (allocated(rce%numden_grp)) deallocate(rce%numden_grp)
    if (allocated(rce%dfl_w))      deallocate(rce%dfl_w)
    if (allocated(rce%dfs_w))      deallocate(rce%dfs_w)
    if (allocated(rce%tautot_w))   deallocate(rce%tautot_w)
    if (allocated(rce%swt_w))      deallocate(rce%swt_w)
    if (allocated(rce%swb_w))      deallocate(rce%swb_w)
    if (allocated(rce%swr_w))      deallocate(rce%swr_w)
    if (allocated(rce%sw_hit))     deallocate(rce%sw_hit)
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

    return
  end subroutine rce_destroy


  !! The Parmentier (2015) adiabat through ``(t0, p0)``: the closed form of
  !! the gradient PARM_A - PARM_B*T. Returns ``t0`` at ``p == p0``.
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
  !! Format: '#' comment lines, then ``nt np``, the log10 T axis, the log10 p
  !! [bar] axis, and the gradient as nt rows of np values. The specific heat
  !! block that follows is not read.
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
  !! Bilinear in (log10 T, log10 p), silently clamped at the table edges.
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


  !! The adiabatic temperature ratio across the pair of layers ``(iz, iz+1)``,
  !! ``exp(integral of grad_ad dlnp)`` by the midpoint rule: one gradient at
  !! the pair's mean temperature and geometric-mean pressure.
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

    grad = rce_grad(rce, tbar, sqrt(p_lo * p_hi))

    pfact = (p_lo / p_hi) ** max(grad, GRAD_AD_MIN)

    return
  end function rce_pfact


  !! The adiabatic gradient ``d ln T / d ln p`` of the selected adiabat.
  pure function rce_grad(rce, temp, pres) result(grad)

    implicit none

    type(rce_type), intent(in) :: rce
    real(kind=f), intent(in)   :: temp   !! [K]
    real(kind=f), intent(in)   :: pres   !! [Pa]
    real(kind=f)               :: grad

    if (rce%adiabat == I_ADIABAT_TABLE) then
      grad = rce_grad_ad(rce, temp, pres)
    else
      grad = PARM_A - PARM_B * temp
    end if

    return
  end function rce_grad


  !! Dry convective adjustment: mix away any superadiabatic gradient.
  !!
  !! After Ray Pierrehumbert's scheme as implemented in Elspeth Lee's
  !! ``Exo-FMS_column_ck`` (``src/dry_conv_adj_mod.f90``, ``Ray_dry_adj``),
  !! with that code's constant ``kappa`` replaced by the run's adiabat.
  !!
  !! A pair is unstable when the deeper layer is hotter than the adiabat
  !! through the shallower one allows. It is mixed at fixed enthalpy (sum of
  !! ``dp*T`` over the pair) onto ``t(iz) = pfact*t(iz+1)``. Sweeps run down
  !! and then up the column.
  !!
  !! Arrays are bottom-to-top: ``p(1)`` is the deepest layer.
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

    ! Both layers of every interface on or above the adiabat are convective;
    ! "above" covers pairs the sweeps did not finish.
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
    ! A grid-scale (2*dz) oscillation makes isolated neutral interfaces. Left
    ! flagged, those layers would leave the implicit solve, the one update
    ! that damps such a mode. Zones are joined across neutral interfaces, not
    ! flagged layers, so a sawtooth reads as a chain of two-layer islands.
    call rce_purge_thin_zones(nz, neutral, is_conv)

    ! Count zones, again joined only across neutral interfaces.
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
  !! ``neutral(iz)`` being the interface between layers ``iz`` and ``iz+1``.
  !! Only the flags change; temperatures are left as the adjustment left them.
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


  !! Eddy diffusion from mixing length theory [cm^2/s].
  !!
  !! Follows PICASO's ``get_kzz`` (``picaso/climate.py``):
  !!
  !!   Kzz = scalef * H * (l/H)^(4/3) * ( R * F_conv / (rho * Cp) )^(1/3)
  !!
  !! with ``l = kzz_mixl_scale * max(0.1, min(1, dlnT/dlnp / grad_ad)) * H``.
  !! Evaluated in every layer; the floor on ``l/H`` leaves stable layers a
  !! small residual Kzz.
  !!
  !! The convective flux is the emergent flux less the layer's net radiative
  !! flux, then limited by an overshoot floor, rescaled to put the internal
  !! flux through the base, and floored at ``sigma*(0.05*t_int)^4``.
  !!
  !! cgs internally. Arrays are bottom-to-top; ``kzz`` is returned on levels.
  pure subroutine rce_kzz(rce, p, pl, t, fnet, kzz, rescale)

    implicit none

    type(rce_type), intent(in)  :: rce
    real(kind=f), intent(in)    :: p(rce%nz)      !! layer pressure [Pa]
    real(kind=f), intent(in)    :: pl(rce%nz+1)   !! level pressure [Pa]
    real(kind=f), intent(in)    :: t(rce%nz)      !! layer temperature [K]
    real(kind=f), intent(in)    :: fnet(rce%nz+1) !! net upward flux [W/m^2]
    real(kind=f), intent(out)   :: kzz(rce%nz+1)  !! eddy diffusion [cm^2/s]
    !! What the convective flux had to be scaled by at the base of the column.
    real(kind=f), intent(out), optional :: rescale

    integer      :: nz, iz
    real(kind=f) :: chf(rce%nz), kz_lay(rce%nz), p_cgs(rce%nz)
    real(kind=f) :: dtdp(rce%nz)
    real(kind=f) :: r_atmos, cp_kzz, grav_cgs, f_sum, f_target, flx_min
    real(kind=f) :: grad_ad, lapse, mixl, scale_h, rho, fscale, w

    nz = rce%nz

    ! Specific gas constant, and the Cp of KZZ_CP_OVER_R (not rce%cp).
    r_atmos  = RGAS_CGS / rce%wtmol
    cp_kzz   = KZZ_CP_OVER_R * r_atmos
    grav_cgs = rce%grav / ACC_CGS2SI

    do iz = 1, nz
      p_cgs(iz) = p(iz) * PA2BARYE
    end do

    ! Emergent flux, and the internal flux it equals in equilibrium.
    f_sum    = fnet(nz+1) * FLUX_SI2CGS
    f_target = SIGMA_SB_CGS * rce%t_int**4
    flx_min  = SIGMA_SB_CGS * (KZZ_TEFF_MIN * rce%t_int)**4

    ! ---- convective flux ---------------------------------------------------
    ! What the radiation is not carrying, with the level net flux averaged
    ! onto layers. The deepest layer carries the whole flux convectively.
    ! fnet includes any shortwave; PICASO uses the infrared net flux alone.
    chf(1) = f_sum
    do iz = 2, nz
      chf(iz) = f_sum - 0.5_f * (fnet(iz) + fnet(iz+1)) * FLUX_SI2CGS
    end do

    ! Convective overshoot floor; see KZZ_CHF_FALL.
    do iz = 2, nz
      chf(iz) = max(chf(iz), &
                    KZZ_CHF_FALL * (p_cgs(iz) / p_cgs(iz-1)) * chf(iz-1))
    end do

    ! Rescale every layer so the base carries the imposed internal flux.
    if (chf(1) > 0._f) then
      fscale = f_target / chf(1)
    else
      fscale = 1._f
    end if

    if (present(rescale)) rescale = fscale

    do iz = 1, nz
      chf(iz) = max(chf(iz) * fscale, flx_min)
    end do

    ! ---- mixing length -----------------------------------------------------
    ! Lapse rate of the current profile, centred where there are neighbours on
    ! both sides and one-sided at the two ends.
    do iz = 1, nz
      if (nz == 1) then
        dtdp(iz) = 0._f
      else if (iz == 1) then
        dtdp(iz) = log(t(2) / t(1)) / log(p_cgs(2) / p_cgs(1))
      else if (iz == nz) then
        dtdp(iz) = log(t(nz) / t(nz-1)) / log(p_cgs(nz) / p_cgs(nz-1))
      else
        dtdp(iz) = log(t(iz+1) / t(iz-1)) / log(p_cgs(iz+1) / p_cgs(iz-1))
      end if
    end do

    do iz = 1, nz
      grad_ad = max(rce_grad(rce, t(iz), p(iz)), GRAD_AD_MIN)
      lapse   = min(1._f, dtdp(iz) / grad_ad)

      scale_h = r_atmos * t(iz) / grav_cgs
      ! kzz_mixl_scale applies after the floor, to stable layers too.
      mixl    = rce%kzz_mixl_scale * max(KZZ_MIXL_MIN, lapse) * scale_h
      rho     = p_cgs(iz) / (r_atmos * t(iz))

      kz_lay(iz) = KZZ_SCALEF * scale_h * (mixl / scale_h)**(4._f/3._f) &
                   * (r_atmos * chf(iz) / (rho * cp_kzz))**(1._f/3._f)
    end do

    ! ---- onto the levels ---------------------------------------------------
    ! Linear in (ln Kzz, ln p); the boundary levels take the end layers'
    ! values rather than extrapolating.
    kzz(1)    = kz_lay(1)
    kzz(nz+1) = kz_lay(nz)

    do iz = 2, nz
      w = log((pl(iz) * PA2BARYE) / p_cgs(iz-1)) &
          / log(p_cgs(iz) / p_cgs(iz-1))
      kzz(iz) = exp((1._f - w) * log(kz_lay(iz-1)) + w * log(kz_lay(iz)))
    end do

    return
  end subroutine rce_kzz


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


  !! Recompute the altitude grid hydrostatically for the current temperature,
  !! as ``Carma.calculate_z`` does. Only ``I_CART`` grids depend on
  !! temperature; any other grid is left untouched.
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
  !! A proxy for the change in cloud opacity since the last solve; only its
  !! relative change is used, so the factor of pi is omitted.
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
  subroutine rce_solve(rce, p, pl, t)

    implicit none

    type(rce_type), intent(inout) :: rce
    real(kind=f), intent(in) :: p(rce%nz)      !! layer pressure [Pa]
    real(kind=f), intent(in) :: pl(rce%nz+1)   !! level pressure [Pa]
    real(kind=f), intent(in) :: t(rce%nz)      !! layer temperature [K]

    integer      :: nz, iz
    real(kind=f) :: dmass, tau_mean, emiss

    nz = rce%nz

    call rce_fluxes(rce, p, pl, t, rce%fnet, rce%tau_num, rce%tau_den)

    ! Impose the internal flux at the base, overriding toon_lw_column's
    ! extrapolated lower boundary.
    rce%fnet(1) = SIGMA_SB * rce%t_int ** 4

    ! ---- flux divergence -> heating rate ----------------------------------
    do iz = 1, nz
      ! Mass per unit area from dp/g rather than rho*dz, so the heating is
      ! consistent with the dp-weighted enthalpy of the convective adjustment.
      dmass = (pl(iz) - pl(iz+1)) / rce%grav

      rce%dtdt(iz) = -(rce%fnet(iz+1) - rce%fnet(iz)) / (dmass * rce%cp)

      ! Local radiative time constant, heat capacity over 4*sigma*T^3*emiss,
      ! with a grey emissivity from the Planck-weighted optical depth. Used
      ! only to damp the update.
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


  !! Cloud optical depth, single-scattering albedo and asymmetry on the band
  !! centres, for the cloud field currently in ``rce%numden_grp``.
  !!
  !! Independent of temperature, so it runs once per cloud field rather than
  !! inside ``rce_fluxes``. Must run before any flux evaluation that reads
  !! ``rce%tau_c``, ``rce%w0_c`` or ``rce%g_c``.
  subroutine rce_optics_prep(rce, radius, qext, ssa, asym)

    implicit none

    type(rce_type), intent(inout) :: rce
    real(kind=f), intent(in) :: radius(rce%nbin, rce%ngroup)  !! [cm]
    real(kind=f), intent(in) :: qext(rce%nband, rce%nbin, rce%ngroup)
    real(kind=f), intent(in) :: ssa(rce%nband, rce%nbin, rce%ngroup)
    real(kind=f), intent(in) :: asym(rce%nband, rce%nbin, rce%ngroup)

    integer      :: iz
    real(kind=f) :: dz_cm(rce%nz)

    do iz = 1, rce%nz
      dz_cm(iz) = rce%dz(iz) * M2CM
    end do

    call cloud_optics_column(rce%nz, rce%nbin, rce%nband, rce%ngroup, &
                             rce%numden_grp, radius, dz_cm, &
                             qext, ssa, asym, &
                             rce%tau_c, rce%w0_c, rce%g_c)

    return
  end subroutine rce_optics_prep


  !! Net flux at every level for a given temperature profile. Called by
  !! ``rce_solve``, and by ``rce_jacobian`` on perturbed profiles. Level
  !! temperatures and gas opacity are recomputed from ``t``; ``dz`` is held
  !! fixed within a step.
  subroutine rce_fluxes(rce, p, pl, t, fnet, tau_num, tau_den)

    implicit none

    type(rce_type), intent(inout) :: rce
    real(kind=f), intent(in) :: p(rce%nz)      !! layer pressure [Pa]
    real(kind=f), intent(in) :: pl(rce%nz+1)   !! level pressure [Pa]
    real(kind=f), intent(in) :: t(rce%nz)      !! layer temperature [K]

    real(kind=f), intent(out) :: fnet(rce%nz+1)   !! net upward flux [W/m^2]
    real(kind=f), intent(out) :: tau_num(rce%nz)  !! Planck-weighted tau, numerator
    real(kind=f), intent(out) :: tau_den(rce%nz)  !! ... and its denominator

    integer      :: nz, nlev, iz, iband, ig, iw
    !! Per-point column workspace, private to each thread.
    real(kind=f) :: dtau_l(rce%nz), w0_l(rce%nz), g_l(rce%nz)
    real(kind=f) :: fup_l(rce%nz+1), fdn_l(rce%nz+1)
    real(kind=f) :: swup_l(rce%nz+1), swdn_l(rce%nz+1)
    real(kind=f) :: conc(rce%nz), tl(rce%nz+1)
    !! (nz, nband) the part of each layer's emission its bounding levels
    !! cannot see: the layer-centre Planck less the level mean.
    real(kind=f) :: be_corr(rce%nz, rce%nband)
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
      conc(iz) = p(iz) / (RGAS_SI * t(iz))
    end do

    call ck_kappa_column(rce%ck, p, t, conc, rce%beta)

    ! ---- band-integrated Planck at every level, and every layer -----------
    !$OMP PARALLEL DO PRIVATE(iband, iz) SCHEDULE(static)
    do iband = 1, rce%nband
      do iz = 1, nlev
        rce%be(iz, iband) = bbflux_wavenumber(rce%ck%wmin(iband), &
                                              rce%ck%wmax(iband), tl(iz))
      end do

      call bbflux_wavenumber_col(rce%ck%wmin(iband), rce%ck%wmax(iband), &
                                 t(1:nz), rce%be_mid(1:nz, iband))

      do iz = 1, nz
        be_corr(iz, iband) = rce%be_mid(iz, iband) &
                             - 0.5_f * (rce%be(iz, iband) + rce%be(iz+1, iband))
      end do
    end do
    !$OMP END PARALLEL DO

    ! No incident longwave at the top: the auto-emission boundary with
    ! PICASO's factor, in top-down pressure terms.
    btop_factor = pl(nlev) / (pl(nlev-1) - pl(nlev))

    ! ---- solve every spectral point ---------------------------------------
    fnet(:)    = 0._f
    tau_num(:)  = 0._f
    tau_den(:)  = 0._f

    sw_top = 0._f
    sw_bot = 0._f
    sw_ref = 0._f

    rce%sw_hit(:) = .false.

    ! One iteration per spectral point (iband, ig flattened). Each writes only
    ! its own column of the per-point arrays. Static schedule, since every
    ! point costs the same; OMP_SCHEDULE is left to the loop in newstate.F90.
    !$OMP PARALLEL DO &
    !$OMP& PRIVATE(iband, ig, iz, tau_gas, tau_ray, tau_cld, tau_tot, &
    !$OMP&         tau_sca, dtau_l, w0_l, g_l, fup_l, fdn_l, &
    !$OMP&         swup_l, swdn_l) &
    !$OMP& SCHEDULE(static)
    do iw = 1, rce%nwave
        iband = (iw - 1) / rce%ng + 1
        ig    = iw - (iband - 1) * rce%ng

        do iz = 1, nz
          tau_gas = rce%beta(iw, iz) * rce%dz(iz)

          ! Rayleigh scattering, only with an incident beam; it is negligible
          ! in the thermal infrared.
          if (rce%sw_on) then
            tau_ray = rce%ray_fac(iband) * conc(iz) * rce%dz(iz)
          else
            tau_ray = 0._f
          end if

          tau_cld = rce%tau_c(iband, iz)
          tau_tot = tau_gas + tau_cld + tau_ray

          dtau_l(iz) = tau_tot

          ! Scattering optical depth: cloud scattering plus all of Rayleigh.
          tau_sca = tau_cld * rce%w0_c(iband, iz) + tau_ray

          if (tau_tot > 0._f) then
            w0_l(iz) = tau_sca / tau_tot
          else
            w0_l(iz) = 0._f
          end if

          ! The two-stream coefficients are singular at w0 == 1.
          w0_l(iz) = min(w0_l(iz), 1._f - 1.e-12_f)

          ! Scattering-weighted asymmetry; Rayleigh enters with g = 0.
          if (tau_sca > 0._f) then
            g_l(iz) = tau_cld * rce%w0_c(iband, iz) &
                      * rce%g_c(iband, iz) / tau_sca
          else
            g_l(iz) = 0._f
          end if

          ! For the Planck-weighted mean optical depth.
          rce%tautot_w(iz, iw) = tau_tot
        end do

        call toon_lw_column(nz, dtau_l, w0_l, g_l, &
                            rce%be(:, iband), 0._f, -1, btop_factor, &
                            .false., .false., fup_l, fdn_l, &
                            be_corr_in=be_corr(:, iband))

        do iz = 1, nlev
          rce%dfl_w(iz, iw) = fup_l(iz) - fdn_l(iz)
        end do

        ! ---- the incident beam, into the same net flux --------------------
        ! In equilibrium the longwave carries sigma*t_int^4 + absorbed upward at
        ! the top and the shortwave net is -absorbed, so fnet(top) ==
        ! sigma*t_int^4 remains the equilibrium condition.
        if (rce%sw_on) then
          if (rce%f0(iband) > 0._f) then
            call toon_sw_column(nz, rce%f0(iband), rce%mu_lev, &
                                dtau_l, w0_l, g_l, rce%w_surf, &
                                swup_l, swdn_l)

            do iz = 1, nlev
              rce%dfs_w(iz, iw) = swup_l(iz) - swdn_l(iz)
            end do

            rce%swt_w(iw) = swdn_l(nlev) - swup_l(nlev)
            rce%swb_w(iw) = swdn_l(1) - swup_l(1)
            rce%swr_w(iw) = swup_l(nlev)
            rce%sw_hit(iw) = .true.
          end if
        end if
    end do
    !$OMP END PARALLEL DO

    ! Sum in spectral-point order, independent of the thread layout.
    do iw = 1, rce%nwave
      iband = (iw - 1) / rce%ng + 1
      wt    = rce%ck%weights(iw)

      do iz = 1, nz
        be_lay = rce%be_mid(iz, iband)
        tau_num(iz) = tau_num(iz) + wt * be_lay * rce%tautot_w(iz, iw)
        tau_den(iz) = tau_den(iz) + wt * be_lay
      end do

      do iz = 1, nlev
        fnet(iz) = fnet(iz) + wt * rce%dfl_w(iz, iw)
      end do

      if (rce%sw_hit(iw)) then
        do iz = 1, nlev
          fnet(iz) = fnet(iz) + wt * rce%dfs_w(iz, iw)
        end do
        sw_top = sw_top + wt * rce%swt_w(iw)
        sw_bot = sw_bot + wt * rce%swb_w(iw)
        sw_ref = sw_ref + wt * rce%swr_w(iw)
      end if
    end do

    ! Diagnostics only; the energy is already in fnet.
    rce%absorbed_sw  = sw_top - sw_bot
    rce%reflected_sw = sw_ref

    return
  end subroutine rce_fluxes


  !! Jacobian of the heating rate, ``jac(i,k) = d(dTdt_i)/dT_k`` [1/s], by
  !! centred finite differences on the full flux solve.
  !!
  !! Dense: in the upper atmosphere a band of half-width 1 holds only about
  !! half of ``sum_j |dH_j/dT_k|``, and half-width 16 about 94%. With level
  !! temperatures interpolated from layers, perturbing a layer also warms the
  !! levels it shares with its neighbours, so the nearest off-diagonals have
  !! the opposite sign to radiative exchange.
  !!
  !! Costs ``2*nz`` flux evaluations, so it runs on its own cadence.
  subroutine rce_jacobian(rce, p, pl, t)

    implicit none

    type(rce_type), intent(inout) :: rce
    real(kind=f), intent(in) :: p(rce%nz)
    real(kind=f), intent(in) :: pl(rce%nz+1)
    real(kind=f), intent(in) :: t(rce%nz)

    integer      :: nz, iz, k
    real(kind=f) :: dmass, hp, hm

    nz = rce%nz

    do k = 1, nz

      ! Centred difference, since the heating goes as T^4.
      rce%t_pert(:) = t(:)
      rce%t_pert(k) = t(k) + JAC_DT
      call rce_fluxes(rce, p, pl, rce%t_pert, &
                      rce%fnet_pert, rce%tnum_pert, rce%tden_pert)

      rce%t_pert(k) = t(k) - JAC_DT
      call rce_fluxes(rce, p, pl, rce%t_pert, &
                      rce%fnet_minus, rce%tnum_pert, rce%tden_pert)

      ! The base flux is imposed and does not respond to temperature.
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
  !! radiative solve runs only when a cadence trigger fires; the held heating
  !! rate and the convective adjustment are applied on every call.
  !!
  !! ``numden`` is indexed by element; ``elem_is_number`` selects the number
  !! element of each group, since core-mass elements hold a mass
  !! concentration.
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
    real(kind=f) :: dt_col(rce%nz), rhs(rce%nz), sol(rce%nz)
    real(kind=f) :: kzz_new(rce%nz+1), rescale_new
    real(kind=f) :: hcap(rce%nz), e_pre
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
      ! Once per cloud field: every flux evaluation below reads the result.
      call rce_optics_prep(rce, radius, qext, ssa, asym)

      call rce_solve(rce, p, pl, t)

      ! Through locals, since rce_kzz takes rce as intent(in).
      if (rce%kzz_mode == I_KZZ_MIXING_LENGTH) then
        call rce_kzz(rce, p, pl, t, rce%fnet, kzz_new, rescale_new)
        rce%kzz(:)     = kzz_new(:)
        rce%kzz_rescale = rescale_new
      end if

      ! Re-measure the Jacobian if it has never been measured, is older than
      ! JAC_GAP solves, or the profile has drifted by more than JAC_DRIFT_MAX.
      do_jac = .not. rce%have_jac
      if (.not. do_jac) do_jac = rce%jac_age >= JAC_GAP
      if (.not. do_jac .and. rce%jac_age >= JAC_GAP_MIN) then
        do_jac = maxval(abs(t - rce%t_jac)) > JAC_DRIFT_MAX
      end if

      if (do_jac) then
        call rce_jacobian(rce, p, pl, t)
      else
        rce%jac_age = rce%jac_age + 1
      end if

      rce%steps_since_solve = 0
      rce%dt_since_solve    = 0._f
      rce%drift             = 0._f
      rce%xsec_ref          = xsec
    end if

    ! ---- apply the held heating rate to every layer -----------------------
    dt_eff = dtime * rce%accel

    ! Backward Euler on the coupled column,
    !
    !   (C_i/dt) dT_i - sum_k J_ik dT_k = C_i * dtdt_i,
    !
    ! with J the dense Jacobian, divided through by C_i. With J diagonal this
    ! is the per-layer damping dT = dtdt*dt_eff/(1 + dt_eff/tau_rad); as
    ! dt_eff grows it tends to the Newton step dT = -J^-1 dtdt. The fixed point
    ! dtdt = 0 is unchanged.
    !
    ! Solved over the radiative layers only, rows and columns both, so the
    ! convective zones stepped below act as a fixed boundary and their dtdt is
    ! not carried into neighbouring radiative layers.
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
          ! Levenberg damping. The diagonal of J can be positive in the upper
          ! atmosphere and is not sign-clipped.
          rce%lu(j,j) = rce%lu(j,j) + 1._f + rce%lambda
          rhs(j)      = dt_eff * rce%dtdt(rad_idx(j))
        end do

        call lu_factor(nrad, rce%lu(1:nrad,1:nrad), rce%piv(1:nrad), dense_ok)

        if (dense_ok) then
          call lu_solve(nrad, rce%lu(1:nrad,1:nrad), rce%piv(1:nrad), &
                        rhs(1:nrad), sol(1:nrad))

          ! Reject a NaN solution.
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

      ! Raise the damping for the next attempt.
      rce%lambda = max(2._f * rce%lambda, LAMBDA_MIN)
    else
      rce%lambda = rce%lambda * LAMBDA_DECAY
      if (rce%lambda < LAMBDA_MIN) rce%lambda = 0._f
    end if

    ! Convective zones take a diagonal update with one damping factor per
    ! zone, the fastest layer's. A zone equilibrates collectively, when its
    ! dp-weighted net heating is zero, so a factor varying across its layers
    ! would move that equilibrium rather than only the approach to it.
    do iz = 1, nz
      damp(iz) = 1._f / (1._f + dt_eff / rce%tau_rad(iz))
      hcap(iz) = (pl(iz) - pl(iz+1)) / rce%grav * rce%cp
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

    ! Trust region: scale the whole step vector, rather than clipping layer by
    ! layer, so the step keeps its direction.
    dt_lay = 0._f
    do iz = 1, nz
      dt_lay = max(dt_lay, abs(dt_col(iz)))
    end do

    if (dt_lay > DT_TRUST) then
      dt_col(:) = dt_col(:) * (DT_TRUST / dt_lay)
      dt_lay    = DT_TRUST
      rce%nclamp = rce%nclamp + 1

      ! The linear model overreached; raise the damping.
      rce%lambda = max(2._f * rce%lambda, LAMBDA_MIN)
    end if

    ! ---- energy budget diagnostics ----------------------------------------
    rce%e_comp = rce%e_comp + dt_eff * sum(hcap(:) * rce%dtdt(:))
    rce%e_appl = rce%e_appl + sum(hcap(:) * dt_col(:))

    if (.not. allocated(rce%e_comp_z)) then
      allocate(rce%e_comp_z(nz), rce%e_appl_z(nz))
      rce%e_comp_z(:) = 0._f
      rce%e_appl_z(:) = 0._f
    end if

    rce%e_comp_z(:) = rce%e_comp_z(:) + dt_eff * hcap(:) * rce%dtdt(:)
    rce%e_appl_z(:) = rce%e_appl_z(:) + hcap(:) * dt_col(:)

    do iz = 1, nz
      if (rce%is_conv(iz)) then
        rce%e_comp_cnv = rce%e_comp_cnv + dt_eff * hcap(iz) * rce%dtdt(iz)
        rce%e_appl_cnv = rce%e_appl_cnv + hcap(iz) * dt_col(iz)
      else
        rce%e_comp_rad = rce%e_comp_rad + dt_eff * hcap(iz) * rce%dtdt(iz)
        rce%e_appl_rad = rce%e_appl_rad + hcap(iz) * dt_col(iz)
      end if
    end do

    do iz = 1, nz
      t(iz) = t(iz) + dt_col(iz)
    end do

    ! ---- mix away whatever the heating left superadiabatic ----------------
    e_pre = sum(hcap(:) * t(:))
    call rce_conv_adj(rce, p, pl, t, rce%is_conv, rce%nz_rcb, rce%nzone, &
                      rce%conv_resid)
    rce%e_adj = rce%e_adj + (sum(hcap(:) * t(:)) - e_pre)

    ! ---- move the altitude grid to match the new temperature --------------
    call rce_level_temps(nz, p, pl, t, rce%tl)
    call rce_regrid_z(nz, rce%igridv, pl, rce%tl, rce%wtmol, rce%grav, zc, zl)

    rce%steps_since_solve = rce%steps_since_solve + 1
    rce%dt_since_solve    = rce%dt_since_solve + dtime
    rce%drift             = rce%drift + dt_lay

    return
  end subroutine rce_update

end module carma_rce
