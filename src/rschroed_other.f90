module rschroed_other

! Other Schroedinger integrators, not directly used by dftatom, but available
! for reuse. This module contains various rk4 integrators and Adams
! predictor-corrector integrators (both for outward and inward integration).

use types, only: dp
use ode1d_other, only: rk4_integrate_inward2
use utils, only: stop_error
use ode1d, only: adams_extrapolation_outward, adams_interp_outward, &
    rk4_integrate
use ode1d_other, only: adams_extrapolation_inward, adams_interp_inward

implicit none
private
public schroed_inward_adams_pr, schroed_outward_adams_pr, &
    integrate_rschroed_u_inward_rk4, integrate_rschroed_u_rk4

contains

subroutine integrate_rschroed_u_rk4(l, E, R, V, Vmid, P, Q, imax)
! Integrates the Schrodinger equation outward using Runge-Kutta 4th order method
!
! It integrates the Schroedinger equation in the P(r), Q(r) form outwards:
! P' = Q
! Q' = P'' = (2*(V-E) + l*(l+1)/R**2) * P
!
! Returns P(r), Q(r).
integer, intent(in) :: l
real(dp), intent(in) :: E
real(dp), intent(in) :: R(:)
real(dp), intent(in) :: V(:), Vmid(:)
real(dp), intent(out) :: P(:), Q(:)
integer, intent(out) :: imax
real(dp), parameter :: max_val = 1e300_dp

real(dp) :: y(2)
integer :: nr

real(dp), dimension(size(R)) :: C1, C2
real(dp), dimension(size(R)-1) :: C1mid, C2mid
real(dp), dimension(size(R)-1) :: Rmid

! y(:) are values of the components (2) of the equation, in our case:
! y(1) = P, y(2) = Q
! dydx(:) are the derivatives of the components (2) of the equation

nr = size(R)

y(1) = R(1)**(l+1)
y(2) = (l+1)*R(1)**l

C1 = 2*(V-E) + l*(l+1)/R**2
C2 = 0
Rmid = (R(:size(R)-1) + R(2:)) / 2
C1mid = 2*(Vmid-E) + l*(l+1)/Rmid**2
C2mid = 0
call rk4_integrate(R, y, C1, C2, C1mid, C2mid, max_val, P, Q, imax)
end subroutine

subroutine integrate_rschroed_u_inward_rk4(n, l, E, R, V, Vmid, &
        P, Q, imin)
! Integrates the Schrodinger equation inward using Runge-Kutta 4th order method
!
! It integrates the Schroedinger equation in the P(r), Q(r) form inwards:
! P' = Q
! Q' = P'' = (2*(V-E) + l*(l+1)/R**2) * P
!
! Returns P(r), Q(r).
integer, intent(in) :: l, n
real(dp), intent(in) :: E
real(dp), intent(in) :: R(:)
real(dp), intent(in) :: V(:), Vmid(:)
real(dp), intent(out) :: P(:), Q(:)
integer, intent(out) :: imin

real(dp) :: y(2), lambda, zeta, sigma
integer :: nr
integer :: i, i_max
real(dp), parameter :: max_val = 1e300_dp
integer, parameter :: asymptotic = 4
real(dp), dimension(size(R)) :: C1
real(dp), dimension(size(R)-1) :: C1mid
real(dp), dimension(size(R)-1) :: Rmid

! y(:) are values of the components (2) of the equation, in our case:
! y(1) = P, y(2) = Q
! dydx(:) are the derivatives of the components (2) of the equation

nr = size(R)

if (E > 0) then
    print *, "E =", E
    call stop_error("E < 0 required")
end if
if (nr < 4) call stop_error("integrate_rschroed_u_inward: n >= 4 required")
i_max = nr
lambda = sqrt(-2*E)
! One possible condition:
!do while (lambda*R(i_max) > 40)
! Start integrating at the 100x classical turning point
do while (R(i_max) > 100 * R(1))
    if (i_max == 2) then
        print *, "E =", E, "lambda =", lambda
        call stop_error("Can't start the inward integration")
    end if
    i_max = i_max - 1
end do
! The simplest asymptotic seems to be working just fine, but just in case, one
! can try more precise ones by changing the 'asymptotic' parameter above
if (asymptotic == 1) then
    ! Simplest asymptotic (simplified derivative)
    y(1) = R(i_max)**n * exp(-lambda*R(i_max))
    y(2) = - lambda * R(i_max)**n * exp(-lambda*R(i_max))
else if (asymptotic == 2) then
    ! Simple asymptotic (correct derivative)
    y(1) = R(i_max)**n * exp(-lambda*R(i_max))
    y(2) = (n * R(i_max)**(n-1) - lambda * R(i_max)**n) * exp(-lambda*R(i_max))
else if (asymptotic == 3) then
    ! More precise asymptotic:
    zeta = -V(i_max)*R(i_max)
    sigma = zeta/lambda
    y(1) = 0
    y(2) = 0
    ! Use 4 terms (we could use more if needed):
    do i = 0, 4
        y(1) = y(1) + R(i_max)**(sigma-i) * exp(-lambda*R(i_max)) * &
            a_k(i, lambda, sigma, l)
        y(2) = y(2) + R(i_max)**(sigma-i) * exp(-lambda*R(i_max)) * &
            b_k(i, lambda, sigma, l)
    end do
else if (asymptotic == 4) then
    lambda = sqrt(l*(l+1)/R(i_max)**2 + 2 * (V(i_max) - E))
    y(1) = 1
    y(2) = - lambda * R(i_max) * y(1)
else
    call stop_error("Unknown 'asymptotic' parameter")
end if


P(i_max+1:) = 0
Q(i_max+1:) = 0

C1 = 2*(V-E) + l*(l+1)/R**2
Rmid = (R(:size(R)-1) + R(2:)) / 2
C1mid = 2*(Vmid-E) + l*(l+1)/Rmid**2
call rk4_integrate_inward2(R(:i_max), y, C1(:i_max), C1mid(:i_max-1), &
        max_val, P(:i_max), Q(:i_max), imin)
end subroutine

real(dp) recursive function a_k(k, lambda, sigma, l) result(r)
integer, intent(in) :: k, l
real(dp), intent(in) :: lambda, sigma
if (k == 0) then
    r = 1
else
    r = (l*(l+1) - (sigma-k)*(sigma-k+1))/(2*k*lambda)
    r = r * a_k(k-1, lambda, sigma, l)
end if
end function

real(dp) recursive function b_k(k, lambda, sigma, l) result(r)
integer, intent(in) :: k, l
real(dp), intent(in) :: lambda, sigma
if (k == 0) then
    r = -lambda
else
    r = (-l*(l+1) + (sigma+k)*(sigma-k+1))/(2*k)
    r = r * a_k(k-1, lambda, sigma, l)
end if
end function

subroutine schroed_outward_adams_pr(l, E, R, Rp, V, P, Q, imax)
! Integrates the Schrodinger equation outward using Adams 4th order method
!
! It integrates the Schroedinger equation in the P(r), Q(r) form outwards using
! predictor-corrector:
! P' = Q
! Q' = P'' = C * P
! where C = 2*(V-E) + l*(l+1)/R**2
!
! Returns P(r), Q(r).
!
! It tranforms the problem to a uniform mesh 1 <= t <= N + 1, defined by
! R(t) and Rp(t) = dR/dt and Rpp(t) = d^2R/dt^2:
!
! u(t)   = P(R(t))
! up(t)  = du/dt  = dP/dR * dR/dt             = Q * Rp
! upp(t) = dup/dt = dQ/dR * Rp^2 + Q * dRp/dt = C*P*Rp^2 + Q * Rpp
!
! So
!
! up  = u * Rp
! upp = u * C * Rp^2 + up * Rpp / Rp
!
! For example if Rp = al*R, Rpp = al^2 * R, we get:
!
! up  = u * al * R
! upp = u * C * al^2 * R^2 + up * al

integer, intent(in) :: l
real(dp), intent(in) :: E
real(dp), intent(in) :: R(:), Rp(:)
real(dp), intent(in) :: V(:)
real(dp), intent(out) :: P(:), Q(:)
real(dp), parameter :: max_val = 1e6_dp
integer, intent(out) :: imax ! The integration was carried to R(imax)

real(dp), dimension(size(R)) :: C, u1, u2, u1p, u2p
integer :: i, it

C = (l*(l+1)/R**2 + 2 * (V-E))
u1(1:4) = R(1:4) ** (l+1)
u2(1:4) = (l+1) * R(1:4) ** l
u1p(1:4) = Rp(1:4)          * u2(1:4)
u2p(1:4) = Rp(1:4) * C(1:4) * u1(1:4)

do i = 4, size(R)-1
    u1(i+1)  = u1(i)  + adams_extrapolation_outward(u1p, i)
    u2(i+1)  = u2(i)  + adams_extrapolation_outward(u2p, i)
    do it = 1, 2
        u1p(i+1) = Rp(i+1)          * u2(i+1)
        u2p(i+1) = Rp(i+1) * C(i+1) * u1(i+1)
        u1(i+1)  = u1(i) + adams_interp_outward(u1p, i)
        u2(i+1)  = u2(i) + adams_interp_outward(u2p, i)
    end do
    if (abs(u1(i+1)) >= max_val .or. abs(u2(i+1)) >= max_val) then
        P = u1
        Q = u2
        imax = i
        return
    end if
end do
P = u1
Q = u2
imax = size(R)
end subroutine

subroutine schroed_inward_adams_pr(l, E, R, Rp, V, P, Q, imin)
! Integrates the Schrodinger equation inward using Adams 4th order method
!
! It integrates the Schroedinger equation in the P(r), Q(r) form inwards using
! predictor-corrector:
! P' = Q
! Q' = P'' = (2*(V-E) + l*(l+1)/R**2) * P
!
! Returns P(r), Q(r).
!
integer, intent(in) :: l
real(dp), intent(in) :: E
real(dp), intent(in) :: R(:), Rp(:)
real(dp), intent(in) :: V(:)
real(dp), intent(out) :: P(:), Q(:)
real(dp), parameter :: max_val = 1e300_dp
integer, intent(out) :: imin

real(dp), dimension(size(R)) :: C, u1, u2, u1p, u2p
real(dp) :: xkap
integer :: i, it, i_max

C = (l*(l+1)/R**2 + 2 * (V-E))

i_max = size(R)-4
do while (R(i_max) > 100 * R(1))
    if (i_max == 2) then
        call stop_error("Can't start the inward integration")
    end if
    i_max = i_max - 1
end do

xkap = sqrt(l*(l+1) / r(i_max)**2 + 2 * (V(i_max) - E))
u1(i_max:i_max+4) = exp(-xkap * (R(i_max:i_max+4) - R(i_max)))
u2(i_max:i_max+4) = - xkap * u1(i_max:i_max+4)
u1p(i_max:i_max+4) = Rp(i_max:i_max+4)                * u2(i_max:i_max+4)
u2p(i_max:i_max+4) = Rp(i_max:i_max+4) * C(i_max:i_max+4) * u1(i_max:i_max+4)

do i = i_max, 2, -1
    u1(i-1)  = u1(i)  + adams_extrapolation_inward(u1p, i)
    u2(i-1)  = u2(i)  + adams_extrapolation_inward(u2p, i)
    do it = 1, 2
        u1p(i-1) = Rp(i-1)          * u2(i-1)
        u2p(i-1) = Rp(i-1) * C(i-1) * u1(i-1)
        u1(i-1)  = u1(i) + adams_interp_inward(u1p, i)
        u2(i-1)  = u2(i) + adams_interp_inward(u2p, i)
    end do
    if (abs(u1(i-1)) >= max_val .or. abs(u2(i-1)) >= max_val) then
        P = u1
        Q = u2
        P(i_max+1:) = 0
        Q(i_max+1:) = 0
        imin = i
        return
    end if
end do
P = u1
Q = u2
P(i_max+1:) = 0
Q(i_max+1:) = 0
imin = 1
end subroutine

end module
