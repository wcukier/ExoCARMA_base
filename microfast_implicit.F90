! LLM-GENERATED — UNVERIFIED (tier: physics)
! Remove this banner only after completing the checks in Assistant/POLICY.md §2.
! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! Fast microphysics for one layer and substep with the gas depletion solved
!! implicitly (backward Euler in the gas).  Without do_implicit_gas this is a
!! single explicit microfast pass.
!!
!! Finds the gas g* for which the rates evaluated at g* return g* at the end of
!! the substep, so every process sees the same end-of-substep gas and no
!! reservoir is overdrawn.  The final gas is always the start-of-substep gas less
!! the condensate formed, so mass is conserved whatever g* is; convergence only
!! sets how consistent the rates are with that final gas.
!!
!! The explicit pass is kept when no gas changes by more than implicit_trigger.
!! If the iteration does not converge a retry is requested; on the last allowed
!! retry the best iterate is accepted.  Temperatures and latent heating stay
!! explicit: each pass starts from the start-of-substep temperature.
subroutine microfast_implicit(carma, cstate, iz, rc, maxrate, dtime_in, nretries_in, dt_suggest)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_planet_mod
  use carma_condensate_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(inout)      :: carma       !! the carma object
  type(carmastate_type), intent(inout) :: cstate      !! the carma state object
  integer, intent(in)                  :: iz          !! z index
  integer, intent(inout)               :: rc          !! return code, negative indicates failure
  real(kind=f), intent(in)             :: maxrate     !! PETER
  real(kind=f), intent(in)             :: dtime_in    !! thread-local dtime for this substep
  real(kind=f), intent(in)             :: nretries_in !! thread-local retry count
  real(kind=f), intent(inout)          :: dt_suggest  !! substep size the failing rates could tolerate [s]

  ! Local declarations
  integer, parameter        :: MAXBT    = 8          ! backtracking halvings per Newton step
  real(kind=f), parameter   :: FD_STEP  = 1.e-5_f    ! finite-difference step in ln(g)
  real(kind=f), parameter   :: MAX_DLN  = 2._f       ! largest change of ln(g) per Newton step
  real(kind=f), parameter   :: ACT_FRAC = 1.e-3_f    ! a gas joins the solve when it changes by more than
                                                     ! ACT_FRAC*implicit_trigger in the explicit pass
  real(kind=f), parameter   :: ATOL_FRAC = 1.e-12_f  ! absolute tolerance, as a fraction of the starting gas

  real(kind=f)              :: pc0(NBIN,NELEM)       ! start-of-substep state of this layer
  real(kind=f)              :: pcnucl0(NBIN,NELEM)
  real(kind=f)              :: gc0(NGAS)
  real(kind=f)              :: t_save
  real(kind=f)              :: rlheat0
  real(kind=f)              :: partheat0

  integer                   :: iact(NGAS)            ! gases solved for
  integer                   :: nact
  integer                   :: igas
  integer                   :: i, j, k, iter, ibt, ipiv
  integer                   :: npass                 ! microfast passes spent on this substep
  integer                   :: method                ! 0 = Newton, 1 = bracketed fallback, 2 = bracketed, at a switch point
  logical                   :: at_switch             ! the bracketed solve stopped on a jump in the residual
  integer                   :: rc_pass
  logical                   :: converged
  logical                   :: accepted
  real(kind=f)              :: y(NGAS)               ! ln of the trial gas
  real(kind=f)              :: fres(NGAS)            ! residual (g_new - g*) / gc0
  real(kind=f)              :: ytry(NGAS)
  real(kind=f)              :: ftry(NGAS)
  real(kind=f)              :: ybest(NGAS)
  real(kind=f)              :: phibest
  real(kind=f)              :: gnew(NGAS)
  real(kind=f)              :: jac(NGAS,NGAS)
  real(kind=f)              :: dlng(NGAS)
  real(kind=f)              :: phi, phitry, lambda, dmax, piv, tmp
  real(kind=f)              :: dt_dummy
  real(kind=f)              :: gstar(NGAS)
  integer                   :: isolve                ! Newton solves this substep (a key swap forces a second)
  integer                   :: iexp                  ! passes over the active set (a gas left out may join)
  logical                   :: added
  logical                   :: swapped
  real(kind=f), parameter   :: KEY_HYST = 0.1_f      ! a reactant takes over as key only once its supply is
                                                     ! this fraction below the current key's
  integer, parameter        :: MAXSWEEP = 10         ! sweeps over the solved gases in the bracketed fallback
  integer, parameter        :: MAXBRK   = 60         ! bracket-shrinking steps per gas in the bracketed fallback
  integer, parameter        :: NILLIN   = 10         ! Illinois steps before the bracketed fallback bisects
  integer, parameter        :: MAXEXPND = 40         ! doublings when searching upward for a bracket
  real(kind=f), parameter   :: LO_FRAC  = 1.e-30_f   ! lower bracket end, as a fraction of the starting gas

  ! Snapshot everything a microfast pass changes for this layer, so that each
  ! pass starts from the start-of-substep state.
  pc0(:,:)     = pc(iz,:,:)
  pcnucl0(:,:) = pc_nucl(iz,:,:)
  gc0(:)       = gc(iz,:)
  t_save       = t(iz)
  rlheat0      = rlheat(iz)
  partheat0    = partheat(iz)

  ! Key reactants are chosen from the start-of-substep gas and held fixed while
  ! the rates are iterated, since a key change is a discontinuity in the rates.
  call select_keys(gc0, swapped)

  if (.not. do_implicit_gas) then
    call microfast(carma, cstate, iz, rc, maxrate, dtime_in, nretries_in, dt_suggest, gc0, 2)
    return
  end if

  ! Explicit pass: rates at the start-of-substep gas.  Kept if no gas moves far.
  rc_pass  = RC_OK
  dt_dummy = dt_suggest
  call microfast(carma, cstate, iz, rc_pass, maxrate, dtime_in, nretries_in, dt_dummy, gc0, 1)
  if (rc_pass < RC_OK) then
    rc = rc_pass
    return
  end if
  npass = 1

  gnew(:) = gc(iz,:)
  if (rc_pass == RC_OK) then
    if (all(abs(gnew(:) - gc0(:)) <= implicit_trigger * gc0(:))) then
      dt_suggest = dt_dummy
      return
    end if
  end if

  ! Gases the explicit pass moved appreciably are solved for; the others keep the
  ! start-of-substep value for their rates.
  nact = 0
  do igas = 1, NGAS
    if (gc0(igas) <= 0._f) cycle
    if ((gnew(igas) < 0._f) .or. &
        (abs(gnew(igas) - gc0(igas)) > ACT_FRAC * implicit_trigger * gc0(igas))) then
      nact = nact + 1
      iact(nact) = igas
    end if
  end do

  ! Newton iteration in ln(g*) from the start-of-substep gas, whose residual the
  ! explicit pass already gives.  If the converged gas changes a condensate's key
  ! reactant, the solve is repeated once with the new keys.  A gas outside the
  ! solve that the converged state moves appreciably joins it, and the solve is
  ! repeated from the converged gas.
  do i = 1, nact
    y(i) = log(gc0(iact(i)))
  end do

  do iexp = 1, NGAS
  do isolve = 1, 2
  if ((isolve == 1) .and. (iexp == 1)) then
    do i = 1, nact
      fres(i) = (gnew(iact(i)) - gc0(iact(i))) / gc0(iact(i))
    end do
  else
    call residual(y, fres, rc_pass)
    if (rc_pass < RC_OK) then
      rc = rc_pass
      return
    end if
    gnew(:) = gc(iz,:)
  end if

  phi       = sum(fres(1:nact)**2)
  ybest(:)  = y(:)
  phibest   = phi
  converged = is_converged(y, gnew)

  iter = 0
  do while ((.not. converged) .and. (iter < implicit_maxit))
    iter = iter + 1

    ! Finite-difference Jacobian of the residual with respect to ln(g*).
    do j = 1, nact
      ytry(1:nact) = y(1:nact)
      ytry(j)      = y(j) + FD_STEP
      call residual(ytry, ftry, rc_pass)
      if (rc_pass < RC_OK) then
        rc = rc_pass
        return
      end if
      jac(1:nact, j) = (ftry(1:nact) - fres(1:nact)) / FD_STEP
    end do

    ! Solve jac * dy = -fres by Gaussian elimination with partial pivoting.
    dlng(1:nact) = -fres(1:nact)
    do k = 1, nact
      ipiv = k
      do i = k+1, nact
        if (abs(jac(i,k)) > abs(jac(ipiv,k))) ipiv = i
      end do
      if (jac(ipiv,k) == 0._f) exit
      if (ipiv /= k) then
        do j = 1, nact
          tmp = jac(k,j); jac(k,j) = jac(ipiv,j); jac(ipiv,j) = tmp
        end do
        tmp = dlng(k); dlng(k) = dlng(ipiv); dlng(ipiv) = tmp
      end if
      do i = k+1, nact
        piv = jac(i,k) / jac(k,k)
        jac(i,k:nact) = jac(i,k:nact) - piv * jac(k,k:nact)
        dlng(i) = dlng(i) - piv * dlng(k)
      end do
    end do
    if (any(jac(1:nact,1:nact) /= jac(1:nact,1:nact))) exit

    do k = nact, 1, -1
      if (jac(k,k) == 0._f) exit
      dlng(k) = (dlng(k) - sum(jac(k,k+1:nact) * dlng(k+1:nact))) / jac(k,k)
    end do
    if (k > 0) exit

    ! Bound the step in ln(g), then backtrack until the residual norm drops.
    dmax = maxval(abs(dlng(1:nact)))
    if (dmax > MAX_DLN) dlng(1:nact) = dlng(1:nact) * (MAX_DLN / dmax)

    lambda   = 1._f
    accepted = .false.
    do ibt = 0, MAXBT
      ytry(1:nact) = y(1:nact) + lambda * dlng(1:nact)
      call residual(ytry, ftry, rc_pass)
      if (rc_pass < RC_OK) then
        rc = rc_pass
        return
      end if
      phitry = sum(ftry(1:nact)**2)
      if (phitry < (1._f - 1.e-4_f * lambda) * phi) then
        accepted = .true.
        exit
      end if
      lambda = 0.5_f * lambda
    end do
    if (.not. accepted) exit

    y(1:nact)    = ytry(1:nact)
    fres(1:nact) = ftry(1:nact)
    phi          = phitry
    gnew(:)      = gc(iz,:)
    if (phi < phibest) then
      ybest(1:nact) = y(1:nact)
      phibest       = phi
    end if
    converged = is_converged(y, gnew)
  end do

  ! Newton needs a residual smooth enough for its slope to point at the root;
  ! near an exhausted reservoir it is not (nucleation switching on and off,
  ! bins evaporating away), and the line search stalls.  The bracketed solve
  ! needs only the sign change of the residual.
  method = 0
  if (.not. converged) then
    call bracket_solve(converged, at_switch)
    if (rc < RC_OK) return
    if (converged) method = merge(2, 1, at_switch)
  end if

  if (.not. (converged .and. isolve == 1)) exit

  gstar(:) = gc0(:)
  do i = 1, nact
    gstar(iact(i)) = exp(y(i))
  end do
  call select_keys(gstar, swapped)
  if (.not. swapped) exit
  end do

  if (.not. converged) exit

  added = .false.
  do igas = 1, NGAS
    if (gc0(igas) <= 0._f) cycle
    if (any(iact(1:nact) == igas)) cycle
    if ((gnew(igas) < 0._f) .or. &
        (abs(gnew(igas) - gc0(igas)) > ACT_FRAC * implicit_trigger * gc0(igas))) then
      nact       = nact + 1
      iact(nact) = igas
      y(nact)    = log(gc0(igas))
      added      = .true.
    end if
  end do
  if (.not. added) exit
  end do

  ! Diagnostic: layer, number of solved gases, Newton iterations, passes,
  ! converged flag, best residual, method and the solved gases.
  if (do_printdiag) write(lundiag,'(a,i4,2i3,i5,i2,es13.4e3,i2,*(i3))') 'I', iz, nact, iter, npass, &
    merge(1, 0, converged), phibest, method, (iact(i), i = 1, nact)

  ! Not converged: ask for a shorter substep while retries remain; on the last
  ! one, accept the best iterate (still mass conserving).
  if (.not. converged) then
    if (nretries_in < maxretries) then
      call restore()
      dt_suggest = min(dt_suggest, 0.5_f * dtime_in)
      rc = RC_WARNING_RETRY
      return
    end if
    y(1:nact) = ybest(1:nact)
  end if

  ! Final pass at g*, with the usual positivity and threshold checks.
  call restore()
  do i = 1, nact
    gc(iz,iact(i)) = exp(y(i))
  end do
  call microfast(carma, cstate, iz, rc, maxrate, dtime_in, nretries_in, dt_suggest, gc0, 2)

  return

contains

  !! Bracketed fallback: each solved gas in turn is found in ln(g) with the
  !! others held, sweeping until the joint fixed point converges; Illinois
  !! (safeguarded regula falsi) steps, then bisection.  For a single gas a
  !! bracket exists: with no gas the residual is positive, with ample gas
  !! negative.  A process that switches on at some abundance makes the residual
  !! jump across zero; the bracket then closes on the switch point and its lower
  !! end is accepted once narrower than implicit_rtol.  ok = .false. when no
  !! bracket is found.
  subroutine bracket_solve(ok, jumped_any)
    logical, intent(out) :: ok
    logical, intent(out) :: jumped_any
    integer              :: isweep, jj, ib, ibr, iside
    real(kind=f)         :: ulo, uhi, uc, flo, fhi, fc, gs
    real(kind=f)         :: flo0, fhi0
    real(kind=f)         :: ytmp(NGAS), ftmp(NGAS)
    logical              :: gas_ok
    logical              :: jumped(NGAS)

    ok = .false.
    jumped_any = .false.
    do isweep = 1, MAXSWEEP
      jumped(:) = .false.
      do jj = 1, nact
        ytmp(1:nact) = y(1:nact)
        call residual(ytmp, ftmp, rc_pass)
        if (rc_pass < RC_OK) then
          rc = rc_pass
          return
        end if
        if (ftmp(jj) == 0._f) cycle

        ! Bracket the root of this gas's residual: ulo has f > 0, uhi has f < 0.
        if (ftmp(jj) < 0._f) then
          uhi = ytmp(jj)
          fhi = ftmp(jj)
          ulo = log(LO_FRAC * gc0(iact(jj)))
          ytmp(jj) = ulo
          call residual(ytmp, ftmp, rc_pass)
          if (rc_pass < RC_OK) then
            rc = rc_pass
            return
          end if
          flo = ftmp(jj)
          if (flo <= 0._f) return
        else
          ulo = ytmp(jj)
          flo = ftmp(jj)
          uhi = ulo
          fhi = flo
          do ib = 1, MAXEXPND
            uhi = uhi + log(2._f)
            ytmp(jj) = uhi
            call residual(ytmp, ftmp, rc_pass)
            if (rc_pass < RC_OK) then
              rc = rc_pass
              return
            end if
            fhi = ftmp(jj)
            if (fhi < 0._f) exit
            ulo = uhi
            flo = fhi
          end do
          if (fhi >= 0._f) return
        end if

        flo0 = flo
        fhi0 = fhi

        ! Illinois iteration (regula falsi, halving the retained end's residual
        ! whenever the same end is kept twice), then bisection, which halves
        ! the bracket every step even across a jump.
        iside  = 0
        gas_ok = .false.
        do ibr = 1, MAXBRK
          if (abs(uhi - ulo) <= implicit_rtol) then
            ! Closed on a jump: take the lower end, whose residual is positive,
            ! so the end-of-step gas exceeds the trial gas and stays positive.
            if (do_printdiag) write(lundiag,'(a,i4,i3,4es13.4e3)') 'J', iz, iact(jj), exp(ulo), exp(uhi), flo0, fhi0
            uc = ulo
            gas_ok     = .true.
            jumped(jj) = .true.
            jumped_any = .true.
            exit
          end if

          if (ibr <= NILLIN) then
            uc = (ulo * fhi - uhi * flo) / (fhi - flo)
          else
            uc = 0.5_f * (ulo + uhi)
          end if
          ytmp(jj) = uc
          call residual(ytmp, ftmp, rc_pass)
          if (rc_pass < RC_OK) then
            rc = rc_pass
            return
          end if
          fc = ftmp(jj)
          gs = exp(uc)
          if (abs(fc) * gc0(iact(jj)) <= implicit_rtol * gs + ATOL_FRAC * gc0(iact(jj))) then
            gas_ok = .true.
            exit
          end if

          ! flo0 and fhi0 keep the true residuals at the bracket ends; flo and
          ! fhi carry the Illinois down-weighting.
          if (fc > 0._f) then
            ulo  = uc
            flo  = fc
            flo0 = fc
            if (iside == 1) fhi = 0.5_f * fhi
            iside = 1
          else
            uhi  = uc
            fhi  = fc
            fhi0 = fc
            if (iside == -1) flo = 0.5_f * flo
            iside = -1
          end if
        end do
        if (.not. gas_ok) return
        y(jj) = uc
      end do

      ! Joint check at the updated gases.
      call residual(y, fres, rc_pass)
      if (rc_pass < RC_OK) then
        rc = rc_pass
        return
      end if
      gnew(:) = gc(iz,:)
      phi = sum(fres(1:nact)**2)
      if (phi < phibest) then
        ybest(1:nact) = y(1:nact)
        phibest       = phi
      end if
      if (converged_or_switched(y, gnew, jumped)) then
        ok = .true.
        return
      end if
      if (nact == 1) return
    end do
  end subroutine bracket_solve

  !! As is_converged, but a gas whose bracket closed on a switch point counts
  !! as converged while its end-of-substep gas is not below its trial gas: the
  !! other gases may have moved since its bracket closed, and a negative
  !! residual there would drive it below zero.
  logical function converged_or_switched(yt, gn, jmp)
    real(kind=f), intent(in) :: yt(NGAS)
    real(kind=f), intent(in) :: gn(NGAS)
    logical, intent(in)      :: jmp(NGAS)
    integer                  :: ii
    real(kind=f)             :: gsw, tol

    converged_or_switched = .true.
    do ii = 1, nact
      gsw = exp(yt(ii))
      tol = implicit_rtol * gsw + ATOL_FRAC * gc0(iact(ii))
      if (jmp(ii)) then
        if (gn(iact(ii)) < gsw - tol) then
          converged_or_switched = .false.
          return
        end if
      else if (abs(gn(iact(ii)) - gsw) > tol) then
        converged_or_switched = .false.
        return
      end if
    end do
  end function converged_or_switched

  !! Make each multi-reactant condensate's key the reactant with the smallest
  !! impingement supply per formula unit, n_j v_th,j / nu_j (proportional to
  !! gc_j / M_j**1.5 / nu_j at a common temperature), switching only when the
  !! newcomer undercuts the current key by KEY_HYST.
  subroutine select_keys(gk, changed)
    real(kind=f), intent(in) :: gk(NGAS)
    logical, intent(out)     :: changed
    integer                  :: ig, jg, ifit, icur, ibest
    real(kind=f)             :: scur, sbest, sj

    changed = .false.
    do ig = 1, NGROUP
      ifit = igrowgas(ienconc(ig))
      if (ifit == 0) cycle
      if (count(carma%f_dgas_dcond(ig, :) /= 0._f) < 2) cycle
      icur  = ikey(iz,ig)
      scur  = supply(ig, icur, gk)
      ibest = icur
      sbest = scur
      do jg = 1, NGAS
        if (carma%f_dgas_dcond(ig, jg) == 0._f) cycle
        if ((jg /= ifit) .and. (xref(ig, jg) <= 0._f)) cycle
        sj = supply(ig, jg, gk)
        if (sj < sbest) then
          ibest = jg
          sbest = sj
        end if
      end do
      if ((ibest /= icur) .and. (sbest < (1._f - KEY_HYST) * scur)) then
        if (do_printdiag) write(lundiag,'(a,4i4)') 'K', iz, ig, icur, ibest
        ikey(iz,ig) = ibest
        changed = .true.
      end if
    end do
  end subroutine select_keys

  !! Impingement supply of reactant jg per formula unit of condensate ig.
  real(kind=f) function supply(ig, jg, gk)
    integer, intent(in)      :: ig, jg
    real(kind=f), intent(in) :: gk(NGAS)
    supply = max(gk(jg), 0._f) / (gwtmol_dif(jg) * sqrt(gwtmol_dif(jg))) &
             / real(carma%f_nu_dep(ig, jg), f)
  end function supply

  !! Restore the start-of-substep state of the layer.
  subroutine restore()
    pc(iz,:,:)      = pc0(:,:)
    pc_nucl(iz,:,:) = pcnucl0(:,:)
    gc(iz,:)        = gc0(:)
    t(iz)           = t_save
    rlheat(iz)      = rlheat0
    partheat(iz)    = partheat0
  end subroutine restore

  !! Residual (g_new - g*) / gc0 of a pass whose rates see the trial gas exp(yt).
  subroutine residual(yt, fout, rcp)
    real(kind=f), intent(in)  :: yt(NGAS)
    real(kind=f), intent(out) :: fout(NGAS)
    integer, intent(out)      :: rcp
    integer                   :: ii
    real(kind=f)              :: dts

    call restore()
    do ii = 1, nact
      gc(iz,iact(ii)) = exp(yt(ii))
    end do
    rcp = RC_OK
    dts = huge(1._f)
    call microfast(carma, cstate, iz, rcp, maxrate, dtime_in, nretries_in, dts, gc0, 0)
    npass = npass + 1
    if (rcp >= RC_OK) rcp = RC_OK
    do ii = 1, nact
      fout(ii) = (gc(iz,iact(ii)) - exp(yt(ii))) / gc0(iact(ii))
    end do
  end subroutine residual

  !! The trial gas reproduces itself to within implicit_rtol (relative) or
  !! ATOL_FRAC of the starting gas (absolute) for every solved gas.
  logical function is_converged(yt, gn)
    real(kind=f), intent(in) :: yt(NGAS)
    real(kind=f), intent(in) :: gn(NGAS)
    integer                  :: ii
    real(kind=f)             :: gs

    is_converged = .true.
    do ii = 1, nact
      gs = exp(yt(ii))
      if (abs(gn(iact(ii)) - gs) > implicit_rtol * gs + ATOL_FRAC * gc0(iact(ii))) then
        is_converged = .false.
        return
      end if
    end do
  end function is_converged

end subroutine microfast_implicit
