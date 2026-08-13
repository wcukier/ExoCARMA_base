!! Dense linear algebra for the radiatively coupled version of CARMA.
!!
!! LU factorisation with partial pivoting, and the corresponding triangular
!! solve. Written out rather than called from LAPACK because this build links
!! no external libraries at all -- there is no `-l` anywhere in
!! `src/CARMA/Makefile` -- and `carmapy.exe` is meant to stay self-contained.
!!
!! The only consumer is ``carma_rce``, whose implicit temperature step needs a
!! dense system solved: the Jacobian of the heating rate is measured to be
!! dense on this vertical discretisation, with a tridiagonal capturing only
!! about half of it at every depth. See the Phase 6 notes in
!! ``claude/RCE_STATUS.md``.
!!
!! Kept in its own module, depending on nothing but the precision kind, so its
!! correctness can be checked on its own against a reference implementation
!! rather than only through the physics that uses it.
!!
!! @author Wolf Cukier
module carma_linalg

  use carma_precision_mod

  implicit none

  private
  public :: lu_factor, lu_solve

contains

  !! Factor ``a`` into ``L*U`` in place, with partial pivoting by row.
  !!
  !! On exit the unit-diagonal ``L`` occupies the strict lower triangle and
  !! ``U`` the upper triangle, in the usual compact form; ``piv(k)`` is the row
  !! swapped with row ``k`` at step ``k``. ``ok`` is false if the matrix is
  !! singular to working precision, in which case ``a`` is left partly factored
  !! and must not be used.
  pure subroutine lu_factor(n, a, piv, ok)

    implicit none

    integer, intent(in)         :: n
    real(kind=f), intent(inout) :: a(n,n)   !! in: matrix; out: its factors
    integer, intent(out)        :: piv(n)   !! pivot row chosen at each step
    logical, intent(out)        :: ok

    integer      :: k, i, j, ip
    real(kind=f) :: amax, swap, mult

    ok = .true.

    do k = 1, n - 1

      ! Partial pivoting: without it a zero on the diagonal stops the
      ! factorisation even when the matrix is perfectly well conditioned, and
      ! a small one costs digits.
      ip   = k
      amax = abs(a(k,k))
      do i = k + 1, n
        if (abs(a(i,k)) > amax) then
          amax = abs(a(i,k))
          ip   = i
        end if
      end do

      piv(k) = ip

      if (amax == 0._f) then
        ok = .false.
        return
      end if

      if (ip /= k) then
        do j = 1, n
          swap    = a(k,j)
          a(k,j)  = a(ip,j)
          a(ip,j) = swap
        end do
      end if

      do i = k + 1, n
        mult   = a(i,k) / a(k,k)
        a(i,k) = mult
        do j = k + 1, n
          a(i,j) = a(i,j) - mult * a(k,j)
        end do
      end do

    end do

    piv(n) = n
    if (a(n,n) == 0._f) ok = .false.

    return
  end subroutine lu_factor


  !! Solve ``A*x = b`` from the factors produced by ``lu_factor``.
  !!
  !! ``a`` and ``piv`` are untouched, so one factorisation can be reused for
  !! many right-hand sides -- which is the point of splitting the two, since
  !! the factorisation is the expensive half.
  pure subroutine lu_solve(n, a, piv, b, x)

    implicit none

    integer, intent(in)       :: n
    real(kind=f), intent(in)  :: a(n,n)   !! factors from lu_factor
    integer, intent(in)       :: piv(n)
    real(kind=f), intent(in)  :: b(n)
    real(kind=f), intent(out) :: x(n)

    integer      :: k, i, j, ip
    real(kind=f) :: swap, s

    x(:) = b(:)

    ! Apply the same row interchanges the factorisation made.
    do k = 1, n - 1
      ip = piv(k)
      if (ip /= k) then
        swap  = x(k)
        x(k)  = x(ip)
        x(ip) = swap
      end if
    end do

    ! Forward substitution; L has an implicit unit diagonal.
    do i = 2, n
      s = x(i)
      do j = 1, i - 1
        s = s - a(i,j) * x(j)
      end do
      x(i) = s
    end do

    ! Back substitution.
    x(n) = x(n) / a(n,n)
    do i = n - 1, 1, -1
      s = x(i)
      do j = i + 1, n
        s = s - a(i,j) * x(j)
      end do
      x(i) = s / a(i,i)
    end do

    return
  end subroutine lu_solve

end module carma_linalg
