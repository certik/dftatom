module ode1d

! General utilities for solving 1D ODEs. the Adams and rk4 subroutines
! are used by Schroedinger, Dirac and Poisson solvers. The integrate
! function is used at other places in dftatom to calculate integrals of the
! radial density/orbitals.

use types, only: dp
use utils, only: stop_error

implicit none
private
public integrate, normalize, parsefunction, get_n_nodes, get_min_idx, &
        adams_interp_outward, adams_extrapolation_outward, &
        adams_interp_outward_implicit, &
        adams_interp_inward_implicit, &
        get_midpoints, rk4_integrate3, rk4_integrate4, rk4_integrate, &
        integrate_trapz_1, integrate_trapz_3, integrate_trapz_5, &
        integrate_trapz_7, integrate_simpson, integrate_adams

contains

real(dp) function integrate(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)

! Choose one from the integration rules below:
!s = integrate_trapz_1(Rp, f)
!s = integrate_trapz_3(Rp, f)
!s = integrate_trapz_5(Rp, f)
s = integrate_trapz_7(Rp, f)
!s = integrate_simpson(Rp, f)
!s = integrate_adams(Rp, f)
end function

real(dp) function integrate_trapz_1(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)

real(dp) :: g(size(Rp))
integer :: N
N = size(Rp)
g = f * Rp
s = (g(1) + g(N)) / 2
s = s + sum(g(2:N-1))
end function

real(dp) function integrate_trapz_3(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)

real(dp) :: g(size(Rp))
integer :: N
N = size(Rp)
g = f * Rp
s = (9 * (g(1) + g(N)) + 28 * (g(2) + g(N-1)) + 23 * (g(3) + g(N-2))) / 24
s = s + sum(g(4:N-3))
end function

real(dp) function integrate_trapz_5(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)

real(dp) :: g(size(Rp))
integer :: N
N = size(Rp)
g = f * Rp
s = (  475 * (g(1) + g(N  )) &
    + 1902 * (g(2) + g(N-1)) &
    + 1104 * (g(3) + g(N-2)) &
    + 1586 * (g(4) + g(N-3)) &
    + 1413 * (g(5) + g(N-4)) &
    ) / 1440
s = s + sum(g(6:N-5))
end function

real(dp) function integrate_trapz_7(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)

real(dp) :: g(size(Rp))
integer :: N
N = size(Rp)
g = f * Rp
s = (  36799 * (g(1) + g(N  )) &
    + 176648 * (g(2) + g(N-1)) &
    +  54851 * (g(3) + g(N-2)) &
    + 177984 * (g(4) + g(N-3)) &
    +  89437 * (g(5) + g(N-4)) &
    + 130936 * (g(6) + g(N-5)) &
    + 119585 * (g(7) + g(N-6)) &
    ) / 120960
s = s + sum(g(8:N-7))
end function

real(dp) function integrate_simpson(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)

real(dp) :: g(size(Rp))
integer :: i, N
N = size(Rp)
g = f * Rp
s = 0
do i = 2, N-1, 2
    s = s + g(i-1) + 4*g(i) + g(i+1)
end do
s = s / 3
if (modulo(N, 2) == 0) then
    ! If N is even, add the last slice separately
    s = s + (5*g(N) + 8*g(N-1) - g(N-2)) / 12
end if
end function

real(dp) function integrate_adams(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)

real(dp) :: g(size(Rp))
integer :: i
s = integrate_trapz_1(Rp(:4), f(:4))
g = f * Rp
do i = 4, size(Rp)-1
    s = s + adams_interp_outward(g, i)
end do
end function

subroutine normalize(Rp, Y)
!     normalizes Y inplace on the grid R
!     I.e. int |Y|^2 dr = 1
!     input parameters:
!     Y (array, in, out) .... the unnormalized wavefunction
!     R (array, in) .... the radial grid
!     output parameters:
!     Y is normalized inplace
real(dp), intent(in) :: Rp(:)
real(dp), intent(inout) :: Y(size(Rp))

real(dp) :: S
S = integrate(Rp, Y**2)
S = sqrt(abs(S))
if (S > 0) then
    Y = Y / S
else
    call stop_error("normalize: zero function")
end if
end subroutine

subroutine parsefunction(y, nodes, minidx, positive)
! parses the function y, returns:
! nodes: the number of intersection with the x-axis, not counting
!   the beginning and infinity
! minidx: the index of the last minimum, i.e. the place, where the
!   function was last closest to the x-axis (in absolute value), so for
!   indexes > minidx, the function goes to +- infinity, or possibly to
!   zero in infinity
! positive: true or false, depending if the y is approaching the x-axis
!   in the infinity from above (positive) or below (negative).
! This information is is used in the routine checke to determine, if the
! function lies below or above a certain energy.
real(dp), intent(in) :: y(:)
integer, intent(out) :: nodes, minidx
logical, intent(out) :: positive

nodes = get_n_nodes(y)
minidx = get_min_idx(y)
positive = y(size(y)) > 0
end subroutine

integer function get_n_nodes(y) result(nodes)
! Returns the number of nodes of the function 'y'
real(dp), intent(in) :: y(:)

integer :: last_sign, last_i, i, isy
nodes = 0
last_sign = int(sign(1.0_dp, y(1)))
last_i = -1

do i = 2, size(y)
  isy = int(sign(1.0_dp, y(i)))
  if (isy == -last_sign) then
      last_sign = isy
      last_i = i - 1
      nodes = nodes + 1
  end if
end do
end function

integer function get_min_idx(y) result(k)
! Returns the index of the last minimum of the function 'y'
real(dp), intent(in) :: y(:)
k = size(y)
do while (abs(y(k-1)) < abs(y(k)))
  k = k-1
  if (k == 1) exit
end do
k = k - 1
end function

real(dp) function adams_extrapolation_outward(y, i) result(r)
! Adams extrapolation formula
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = +(55*y(i) - 59*y(i-1) + 37*y(i-2) - 9*y(i-3)) / 24
end function

real(dp) function adams_interp_outward(y, i) result(r)
! Adams interpolation formula
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = +(9*y(i+1) + 19*y(i) - 5*y(i-1) + y(i-2)) / 24
end function

real(dp) function adams_interp_outward_implicit(y, i) result(r)
! Adams interpolation formula for implicit method
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = +(19*y(i) - 5*y(i-1) + y(i-2)) / 24
end function

real(dp) function adams_interp_inward_implicit(y, i) result(r)
! Adams interpolation formula for implicit method
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = -(19*y(i) - 5*y(i+1) + y(i+2)) / 24
end function

subroutine interp(t, X, n, Y, val)
!     Interpolates the x, y points to estimate the value at the point "t".
!     val (out) ...... the estimated function value at the point "t".
!     t (in).......... the "x" value at which to estimate the value
!     X (array, in) .. the x-values of the points
!     Y (array, in) .. the y-values of the points
!     n (in) ......... the length of the arrays X and Y
real(dp), intent(in) :: t
integer, intent(in) :: n
real(dp), intent(in) :: X(n), Y(n)
real(dp), intent(out) :: val
real(dp) :: f, denum
integer :: i, j
val = 0
do j = 1, n
    f = 1
    denum = 1
    do i=1, n
        if (i == j) cycle
        f = f * (t-X(i))
        denum = denum * (X(j)-X(i))
    end do
    val = val + Y(j) * f / denum
end do
end subroutine


subroutine get_val(f, x, R, n, V, i)
! returns the interpolated value of f(x), where f is defined on the
! grid R
! input parameters
! R ... grid
! n ... length of the grid
! f ... function defined on the grid R
! x ... point at which to interpolate
! i ... the integer number of the radial point at which "x" lies
! output parameters
! V ... the interpolated function value
integer, intent(in) :: n
real(dp), intent(in) :: f(n), x, R(n)
real(dp), intent(out) :: V
integer, intent(in) :: i

integer :: j1, j2, n1, n2
if (n < 4) call stop_error("get_val: n >= 4 required")
j1 = i-1
j2 = j1+1

n1 = j1-1
n2 = j2+1
if (n1 < 1) then
    n2 = n2-n1+1
    n1 = 1
end if
if (n2 > n) then
    n1 = n1-(n2-n)
    n2 = n
end if
call interp(x, r(n1:n2), n2-n1+1, f(n1:n2), V)
end subroutine

function get_midpoints(R, V) result(Vmid)
real(dp), intent(in) :: R(:), V(:)
real(dp) :: Vmid(size(R)-1)
integer :: i
if (.not.(size(R) == size(V))) then
    call stop_error("get_midpoints: incorrect array sizes")
end if
do i = 1, size(R) - 1
    call get_val(V, (R(i) + R(i+1))/2, R, size(R), Vmid(i), i+1)
end do
end function

subroutine rk4_integrate3(R, y0, C1, C2, C1mid, C2mid, max_val, y1, y2, imax)
! Integrates the following set of equations outwards:
! dy1/dx =           y2
! dy2/dx = C1 + C2 * y2
! The above system can be equivalently written as a second order ODE:
! y1'' - C2 * y1' = C1
! The coefficients C1 and C2 are passed in as arrays, so they can have any
! dependence on R. For example for a Poisson equation of the form:
! V'' + 2 * V' / R = -4*pi*n
! we would get C1 = -4*pi*n, C2 = -2/R
real(dp), intent(in) :: R(:) ! Grid
real(dp), intent(in) :: y0(:) ! Initial condition
! Coefficients C1 and C2 at grid points and midpoints:
real(dp), intent(in) :: C1(:), C2(:), C1mid(:), C2mid(:)
! Maximum value (if y1 > max_val, the integration stops)
real(dp), intent(in) :: max_val
! Solution y1 and y2
real(dp), intent(out) :: y1(:), y2(:)
! The integration stops at R(imax)
integer, intent(out) :: imax

integer :: i
integer :: n
real(dp), dimension(size(y0)) :: dym, dyt, yt, dydx, y
real(dp) :: h

n = size(R)
y = y0
y1(1) = y(1)
y2(1) = y(2)
do i = 2, n
    ! rk4 step size
    h = R(i)-R(i-1)

    ! Do rk4 step:
    dydx(1) =                      y(2)
    dydx(2) = C1(i-1)  + C2(i-1) * y(2)
    yt = y + h/2 * dydx
    dyt(1) =                           yt(2)
    dyt(2) = C1mid(i-1) + C2mid(i-1) * yt(2)
    yt = y + h/2 * dyt
    dym(1) =                           yt(2)
    dym(2) = C1mid(i-1) + C2mid(i-1) * yt(2)
    yt = y + h * dym
    dym = dyt + dym
    dyt(1) =                 yt(2)
    dyt(2) = C1(i) + C2(i) * yt(2)
    y = y + h/6 * (dydx + dyt + 2*dym)

    y1(i) = y(1)
    y2(i) = y(2)
    if (abs(y(1)) >= max_val) then
        imax = i
        return
    end if
end do
imax = n
end subroutine

subroutine rk4_integrate4(R, y0, C, Cmid, max_val, y1, y2, imax)
! Integrates the following set of equations outwards:
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
integer, intent(out) :: imax

integer :: i
integer :: n
real(dp), dimension(size(y0)) :: dym, dyt, yt, dydx, y
real(dp) :: h

n = size(R)
y = y0
y1(1) = y(1)
y2(1) = y(2)
do i = 2, n
    ! rk4 step size
    h = R(i)-R(i-1)

    ! Do rk4 step:
    dydx(1) = C(i-1, 1, 1) * y(1) + C(i-1, 1, 2) * y(2)
    dydx(2) = C(i-1, 2, 1) * y(1) + C(i-1, 2, 2) * y(2)
    !dydx = matmul(C(i-1, :, :), y)
    yt = y + h/2 * dydx
    dyt(1) = Cmid(i-1, 1, 1) * yt(1) + Cmid(i-1, 1, 2) * yt(2)
    dyt(2) = Cmid(i-1, 2, 1) * yt(1) + Cmid(i-1, 2, 2) * yt(2)
    !dyt = matmul(Cmid(i-1, :, :), yt)
    yt = y + h/2 * dyt
    dym(1) = Cmid(i-1, 1, 1) * yt(1) + Cmid(i-1, 1, 2) * yt(2)
    dym(2) = Cmid(i-1, 2, 1) * yt(1) + Cmid(i-1, 2, 2) * yt(2)
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
        imax = i
        return
    end if
end do
imax = n
end subroutine

subroutine rk4_integrate(R, y0, C1, C2, C1mid, C2mid, max_val, y1, y2, imax)
! Integrates the following set of equations outwards:
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
integer, intent(out) :: imax

integer :: i
integer :: n
real(dp), dimension(size(y0)) :: dym, dyt, yt, dydx, y
real(dp) :: h

n = size(R)
y = y0
y1(1) = y(1)
y2(1) = y(2)
do i = 2, n
    ! rk4 step size
    h = R(i)-R(i-1)

    ! Do rk4 step:
    dydx(1) =                            y(2)
    dydx(2) = C1(i-1) * y(1) + C2(i-1) * y(2)
    yt = y + h/2 * dydx
    dyt(1) =                                   yt(2)
    dyt(2) = C1mid(i-1) * yt(1) + C2mid(i-1) * yt(2)
    yt = y + h/2 * dyt
    dym(1) =                                   yt(2)
    dym(2) = C1mid(i-1) * yt(1) + C2mid(i-1) * yt(2)
    yt = y + h * dym
    dym = dyt + dym
    dyt(1) =                         yt(2)
    dyt(2) = C1(i) * yt(1) + C2(i) * yt(2)
    y = y + h/6 * (dydx + dyt + 2*dym)

    y1(i) = y(1)
    y2(i) = y(2)
    if (abs(y(1)) >= max_val) then
        imax = i
        return
    end if
end do
imax = n
end subroutine

end module
