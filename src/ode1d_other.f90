module ode1d_other

! General utilities for solving 1D ODEs, not used directly by dftatom. They are
! available here for reuse.

use types, only: dp
use constants, only: pi
use utils, only: stop_error
use ode1d, only: get_midpoints

implicit none
private
public rk4_integrate_inward2, rk4_integrate_inward, rk4_integrate_inward4, &
        adams_extrapolation_inward, adams_interp_inward, &
        integrate_simpson_direct, eulerstep, rk4step, rk4step2, &
        adams_interp_inward_8, adams_interp_outward_8, &
        adams_extrapolation_inward_8, adams_extrapolation_outward_8, &
        adams_extrapolation_inward_6, adams_interp_outward_6, &
        adams_interp_inward_6, adams_extrapolation_outward_6


contains

subroutine rk4_integrate_inward2(R, y0, C1, C1mid, max_val, y1, y2, imin)
! Integrates the following set of equations inwards:
! dy1/dx =      y2
! dy2/dx = C1 * y1
real(dp), intent(in) :: R(:) ! Grid
real(dp), intent(in) :: y0(:) ! Initial condition
! Coefficients C1 at grid points and midpoints:
real(dp), intent(in) :: C1(:), C1mid(:)
! Maximum value (if y1 > max_val, the integration stops)
real(dp), intent(in) :: max_val
! Solution y1 and y2
real(dp), intent(out) :: y1(:), y2(:)
! The integration stops at R(imax)
integer, intent(out) :: imin

integer :: i
integer :: n
real(dp), dimension(size(y0)) :: dym, dyt, yt, dydx, y
real(dp) :: h

n = size(R)
y = y0
y1(n) = y(1)
y2(n) = y(2)
do i = n-1, 1, -1
    ! rk4 step size
    h = R(i)-R(i+1)

    ! Do rk4 step:
    dydx(1) =           y(2)
    dydx(2) = C1(i+1) * y(1)
    yt = y + h/2 * dydx
    dyt(1) =            yt(2)
    dyt(2) = C1mid(i) * yt(1)
    yt = y + h/2 * dyt
    dym(1) =            yt(2)
    dym(2) = C1mid(i) * yt(1)
    yt = y + h * dym
    dym = dyt + dym
    dyt(1) =         yt(2)
    dyt(2) = C1(i) * yt(1)
    y = y + h/6 * (dydx + dyt + 2*dym)

    y1(i) = y(1)
    y2(i) = y(2)
    if (abs(y(1)) >= max_val) then
        imin = i
        return
    end if
end do
imin = 1
end subroutine

real(dp) function adams_extrapolation_outward_6(y, i) result(r)
! Adams extrapolation formula
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = +(4277 * y(i  ) &
    -7923  * y(i-1) &
    +9982  * y(i-2) &
    -7298  * y(i-3) &
    +2877  * y(i-4) &
    -475   * y(i-5) &
    ) / 1440
end function

real(dp) function adams_extrapolation_inward_6(y, i) result(r)
! Adams extrapolation formula
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = -(4277 * y(i  ) &
    -7923  * y(i+1) &
    +9982  * y(i+2) &
    -7298  * y(i+3) &
    +2877  * y(i+4) &
    -475   * y(i+5) &
    ) / 1440
end function

real(dp) function adams_interp_outward_6(y, i) result(r)
! Adams interpolation formula
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = +(475 * y(i+1) &
    +1427 * y(i  ) &
    -798  * y(i-1) &
    +482  * y(i-2) &
    -173  * y(i-3) &
    +27   * y(i-4) &
    ) / 1440
end function

real(dp) function adams_interp_inward_6(y, i) result(r)
! Adams interpolation formula
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = -(475 * y(i-1) &
    +1427 * y(i  ) &
    -798  * y(i+1) &
    +482  * y(i+2) &
    -173  * y(i+3) &
    +27   * y(i+4) &
    ) / 1440
end function

real(dp) function adams_extrapolation_outward_8(y, i) result(r)
! Adams extrapolation formula
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = +(434241  * y(i  ) &
    -1152169  * y(i-1) &
    +2183877  * y(i-2) &
    -2664477  * y(i-3) &
    +2102243  * y(i-4) &
    -1041723  * y(i-5) &
    +295767   * y(i-6) &
    -36799    * y(i-7) &
    ) / 120960
end function

real(dp) function adams_extrapolation_inward_8(y, i) result(r)
! Adams extrapolation formula
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = -(434241  * y(i  ) &
    -1152169  * y(i+1) &
    +2183877  * y(i+2) &
    -2664477  * y(i+3) &
    +2102243  * y(i+4) &
    -1041723  * y(i+5) &
    +295767   * y(i+6) &
    -36799    * y(i+7) &
    ) / 120960
end function

real(dp) function adams_interp_outward_8(y, i) result(r)
! Adams interpolation formula
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = +(36799 * y(i+1) &
    +139849 * y(i  ) &
    -121797 * y(i-1) &
    +123133 * y(i-2) &
    -88547  * y(i-3) &
    +41499  * y(i-4) &
    -11351  * y(i-5) &
    +1375   * y(i-6) &
    ) / 120960
end function

real(dp) function adams_interp_inward_8(y, i) result(r)
! Adams interpolation formula
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = -(36799 * y(i-1) &
    +139849 * y(i  ) &
    -121797 * y(i+1) &
    +123133 * y(i+2) &
    -88547  * y(i+3) &
    +41499  * y(i+4) &
    -11351  * y(i+5) &
    +1375   * y(i+6) &
    ) / 120960
end function

real(dp) function adams_interp_inward(y, i) result(r)
! Adams interpolation formula
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = -(9*y(i-1) + 19*y(i) - 5*y(i+1) + y(i+2)) / 24
end function

real(dp) function adams_extrapolation_inward(y, i) result(r)
! Adams extrapolation formula
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = -(55*y(i) - 59*y(i+1) + 37*y(i+2) - 9*y(i+3)) / 24
end function

subroutine rk4step(y, x, h, derivs)
!advances the solution of
!y_i'(x)=f_i(x,y_1,y_2,...,y_N)       for i=1,2,...,N; f_i are known
!from x to x+h.
!parameters:
!  y ... y(x) on the input and y(x+h) on the output.
!  derivs(x,y,dydx) .... returns a vector of f_i(x,y1,y2,...)
integer, parameter :: N=2
real(dp), intent(in) :: h, x
real(dp), intent(inout) :: y(N)
interface
    subroutine derivs(x, y, dydx)
        use types
        integer, parameter :: N=2
        real(dp), intent(in) :: x, y(N)
        real(dp), intent(out) :: dydx(N)
    end subroutine
end interface
real(dp) :: dym(N), dyt(N), yt(N), dydx(N)
call derivs(x, y, dydx)
yt = y + h/2 * dydx
call derivs(x+h/2, yt, dyt)
yt = y + h/2 * dyt
call derivs(x+h/2, yt, dym)
yt = y + h * dym
dym = dyt + dym
call derivs(x+h, yt, dyt)
y = y + h/6 * (dydx + dyt + 2*dym)
end subroutine

subroutine rk4step2(y, x, h, derivs)
!advances the solution of
!y_i'(x)=f_i(x,y_1,y_2,...,y_N)       for i=1,2,...,N; f_i are known
!from x to x+h.
!parameters:
!  y ... y(x) on the input and y(x+h) on the output.
!  derivs(x,y,dydx) .... returns a vector of f_i(x,y1,y2,...)
integer, parameter :: N=2
real(dp), intent(in) :: h, x
real(dp), intent(inout) :: y(N)
interface
    subroutine derivs(x, y, x_index, dydx)
        use types
        integer, parameter :: N=2
        integer, intent(in) :: x_index ! 1, 2, 3
        real(dp), intent(in) :: x, y(N)
        real(dp), intent(out) :: dydx(N)
    end subroutine
end interface
real(dp) :: dym(N), dyt(N), yt(N), dydx(N)
call derivs(x, y, 1, dydx)
yt = y + h/2 * dydx
call derivs(x+h/2, yt, 2, dyt)
yt = y + h/2 * dyt
call derivs(x+h/2, yt, 2, dym)
yt = y + h * dym
dym = dyt + dym
call derivs(x+h, yt, 3, dyt)
y = y + h/6 * (dydx + dyt + 2*dym)
end subroutine

subroutine eulerstep(y, x, h, derivs)
integer, parameter :: N=2
real(dp), intent(in) :: h, x
real(dp), intent(inout) :: y(N)
interface
    subroutine derivs(x, y, dydx)
        use types
        integer, parameter :: N=2
        real(dp), intent(in) :: x, y(N)
        real(dp), intent(out) :: dydx(N)
    end subroutine
end interface
real(dp) :: dydx(N)
call derivs(x, y, dydx)
y = y + h * dydx
end subroutine

real(dp) function integrate_simpson_direct(x, f) result(s)
! Integrates the function f(x) using Simpson rule
!
! Returns: the integral
real(dp), intent(in), dimension(:) :: x ! the grid
real(dp), intent(in), dimension(:) :: f ! function defined on the grid x

real(dp) :: dx, f_a, f_mid(size(x)-1), f_b
integer :: i
s = 0
f_mid = get_midpoints(x, f)
do i = 1, size(x)-1
    dx = x(i+1) - x(i)
    f_a = f(i)
    f_b = f(i+1)
    s = s + dx/6 * (f_a + 4*f_mid(i) + f_b)
end do
end function

subroutine rk4_integrate_inward(R, y0, C1, C2, C1mid, C2mid, max_val, &
    y1, y2, imin)
! Integrates the following set of equations inwards:
! dy1/dx =                y2
! dy2/dx = C1 * y1 + C2 * y2
real(dp), intent(in) :: R(:) ! Grid
real(dp), intent(in) :: y0(:) ! Initial condition
! Coefficients C1 and C2 at grid points and midpoints:
real(dp), intent(in) :: C1(:), C2(:), C1mid(:), C2mid(:)
! Maximum value (if y1 > max_val, the integration stops)
real(dp), intent(in) :: max_val
! Solution y1 and y2
real(dp), intent(out) :: y1(:), y2(:)
! The integration stops at R(imax)
integer, intent(out) :: imin

integer :: i
integer :: n
real(dp), dimension(size(y0)) :: dym, dyt, yt, dydx, y
real(dp) :: h

n = size(R)
y = y0
y1(n) = y(1)
y2(n) = y(2)
do i = n-1, 1, -1
    ! rk4 step size
    h = R(i)-R(i+1)

    ! Do rk4 step:
    dydx(1) =                            y(2)
    dydx(2) = C1(i+1) * y(1) + C2(i+1) * y(2)
    yt = y + h/2 * dydx
    dyt(1) =                               yt(2)
    dyt(2) = C1mid(i) * yt(1) + C2mid(i) * yt(2)
    yt = y + h/2 * dyt
    dym(1) =                               yt(2)
    dym(2) = C1mid(i) * yt(1) + C2mid(i) * yt(2)
    yt = y + h * dym
    dym = dyt + dym
    dyt(1) =                         yt(2)
    dyt(2) = C1(i) * yt(1) + C2(i) * yt(2)
    y = y + h/6 * (dydx + dyt + 2*dym)

    y1(i) = y(1)
    y2(i) = y(2)
    if (abs(y(1)) >= max_val) then
        imin = i
        return
    end if
end do
imin = 1
end subroutine

subroutine rk4_integrate_inward4(R, y0, C, Cmid, max_val, y1, y2, imin)
! Integrates the following set of equations inwards:
! dy1/dx = C(:, 1, 1) * y1 + C(:, 1, 2) * y2
! dy2/dx = C(:, 2, 1) * y1 + C(:, 2, 2) * y2
real(dp), intent(in) :: R(:) ! Grid
real(dp), intent(in) :: y0(:) ! Initial condition
! Coefficients C1 and C2 at grid points and midpoints:
real(dp), intent(in) :: C(:, :, :), Cmid(:, :, :)
! Maximum value (if y1 > max_val, the integration stops)
real(dp), intent(in) :: max_val
! Solution y1 and y2
real(dp), intent(out) :: y1(:), y2(:)
! The integration stops at R(imax)
integer, intent(out) :: imin

integer :: i
integer :: n
real(dp), dimension(size(y0)) :: dym, dyt, yt, dydx, y
real(dp) :: h

n = size(R)
y = y0
y1(n) = y(1)
y2(n) = y(2)
do i = n-1, 1, -1
    ! rk4 step size
    h = R(i)-R(i+1)

    ! Do rk4 step:
    dydx(1) = C(i+1, 1, 1) * y(1) + C(i+1, 1, 2) * y(2)
    dydx(2) = C(i+1, 2, 1) * y(1) + C(i+1, 2, 2) * y(2)
    !dydx = matmul(C(i-1, :, :), y)
    yt = y + h/2 * dydx
    dyt(1) = Cmid(i, 1, 1) * yt(1) + Cmid(i, 1, 2) * yt(2)
    dyt(2) = Cmid(i, 2, 1) * yt(1) + Cmid(i, 2, 2) * yt(2)
    !dyt = matmul(Cmid(i-1, :, :), yt)
    yt = y + h/2 * dyt
    dym(1) = Cmid(i, 1, 1) * yt(1) + Cmid(i, 1, 2) * yt(2)
    dym(2) = Cmid(i, 2, 1) * yt(1) + Cmid(i, 2, 2) * yt(2)
    !dym = matmul(Cmid(i-1, :, :), yt)
    yt = y + h * dym
    dym = dyt + dym
    dyt(1) = C(i, 1, 1) * yt(1) + C(i, 1, 2) * yt(2)
    dyt(2) = C(i, 2, 1) * yt(1) + C(i, 2, 2) * yt(2)
    !dyt = matmul(C(i, :, :), yt)
    y = y + h/6 * (dydx + dyt + 2*dym)

    y1(i) = y(1)
    y2(i) = y(2)
    if (abs(y(1)) >= max_val .or. abs(y(2)) >= max_val) then
        imin = i
        return
    end if
end do
imin = 1
end subroutine

end module
