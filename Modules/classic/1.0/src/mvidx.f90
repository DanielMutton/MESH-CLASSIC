!> \file
!! Utility to determine location of closest value in an array
!      * June 18/2019 - M. Fortier
!      *   - significantly refactored the code, removing the goto statement
!      * May 12/2004 - L. Solheim
!      *
!>     * Given a monotonic vector v of length n and a value x,
!!     * return the index mvidx such that x is between
!!     * v(mvidx) and v(mvidx+1).
!!     *
!!     * v must be monotonic, either increasing of decreasing.
!!     * there is no check on whether or not this vector is
!!     * monotonic.
!!     *
!!     * This function returns 1 or n-1 if x is out of range.
!!     *
!!     * input:
!!     *
!!     *   real    :: v(n) ...monitonic vector (increasing or decreasing)
!!     *   integer :: n    ...size of v
!!     *   real    :: x    ...single real :: value
!!     *
!!     * output:
!!     *
!!     *   v(mvidx) .le/>= x <=/>= v(mvidx+1)
!!
integer function MVIDX (V, N, X)
  !-----------------------------------------------------------------------
  implicit none

  real, intent(in) :: X, V(N)
  integer, intent(in) :: N
  integer :: JL, JM, JU
  !-----------------------------------------------------------------------

  if (X == V(1)) then
    MVIDX = 1
    return
  else if (X == V(N)) then
    MVIDX = N - 1
    return
  end if
  JL = 1
  JU = N
  do while (JU - JL > 1)
    JM = (JU + JL) / 2
    if ((V(N) > V(1)).eqv.(X > V(JM))) then
      JL = JM
    else
      JU = JM
    end if
  end do
  MVIDX = JL
  return
end function MVIDX
