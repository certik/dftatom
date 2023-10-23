module rschroed

! Routines in this module solve the radial Schroedinger equation outward and
! inward using the implicit Adams method.

use types, only: dp
use ode1d, only: adams_interp_outward_implicit, &
    adams_interp_inward_implicit, get_midpoints, rk4_integrate
use utils, only: stop_error

implicit none
private
public schroed_outward_adams, schroed_inward_adams

contains

subroutine schroed_outward_adams(l, Z, E, R, Rp, V, P, Q, imax)
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
integer, intent(in) :: Z
real(dp), intent(in) :: E
real(dp), intent(in) :: R(:), Rp(:)
real(dp), intent(in) :: V(:)
real(dp), intent(out) :: P(:), Q(:)
real(dp), parameter :: max_val = 1e6_dp
integer, intent(out) :: imax ! The integration was carried to R(imax)

real(dp), dimension(size(R)) :: C, u1, u2, u1p, u2p, Vmid
integer :: i
real(dp) :: lambda, Delta, M(2, 2), u1_tmp, u2_tmp

if (size(R) < 5) call stop_error("size(R) < 5")
if (.not. (size(R) == size(Rp) .and. size(R) == size(V) .and. &
    size(R) == size(P) .and. size(R) == size(P) .and. size(R) == size(Q))) then
    call stop_error("Array sizes mismatch")
end if
C = (l*(l+1)/R**2 + 2 * (V-E))
Vmid(:3) = get_midpoints(R(:4), V(:4))
call integrate_rschroed_rk4(l, Z, E, R(:4), V(:4), Vmid(:3), &
    u1(:4), u2(:4), imax)
!u1(1:4) = R(1:4) ** (l+1)
!u2(1:4) = (l+1) * R(1:4) ** l
u1p(1:4) = Rp(1:4)          * u2(1:4)
u2p(1:4) = Rp(1:4) * C(1:4) * u1(1:4)

do i = 4, size(R)-1
    u1p(i) = Rp(i)        * u2(i)
    u2p(i) = Rp(i) * C(i) * u1(i)
    u1_tmp  = u1(i) + adams_interp_outward_implicit(u1p, i)
    u2_tmp  = u2(i) + adams_interp_outward_implicit(u2p, i)

    lambda = 9.0_dp / 24
    Delta = 1 - lambda**2 * C(i+1) * Rp(i+1)**2
    M(1, 1) = 1 / Delta
    M(2, 1) = lambda * C(i+1) * Rp(i+1) / Delta
    M(1, 2) = lambda * Rp(i+1) / Delta
    M(2, 2) = 1 / Delta

    u1(i+1) = M(1, 1) * u1_tmp + M(1, 2) * u2_tmp
    u2(i+1) = M(2, 1) * u1_tmp + M(2, 2) * u2_tmp
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

subroutine schroed_inward_adams(l, E, R, Rp, V, P, Q, imin)
! Integrates the Schrodinger equation inward using Adams 4th order method
!
! It integrates the Schroedinger equation in the P(r), Q(r) form inwards using
! predictor-corrector:
! P' = Q
! Q' = P'' = (2*(V-E) + l*(l+1)/R**2) * P
!
! Returns P(r), Q(r).
!
! Important note: This only works for exponential meshes, where the mesh
! marameter (as defined by mesh_exp()) is equal to "a = (rmax/rmin)**((N-1)/N)".
! Otherwise it will produce wrong answer.
integer, intent(in) :: l
real(dp), intent(in) :: E
real(dp), intent(in) :: R(:), Rp(:)
real(dp), intent(in) :: V(:)
real(dp), intent(out) :: P(:), Q(:)
real(dp), parameter :: max_val = 1e300_dp
integer, intent(out) :: imin

real(dp), dimension(size(R)) :: C, u1, u2, u1p, u2p
integer :: i, i_max
real(dp) :: lambda, Delta, M(2, 2), u1_tmp, u2_tmp
real(dp) :: R_max

C = (l*(l+1)/R**2 + 2 * (V-E))

i_max = size(R)-4
if (i_max < 4) call stop_error("imax < 4")
lambda = sqrt(-2*E)
! We require that exp(-lambda*(R-R(1)) ~ epsilon(1.0_dp),
! if we start further from
! the origin, it might sometimes blow up, if we start closer, we might not get
! as precise asymptotic.
! It follows that R ~ R(1)-log(epsilon(1.0_dp)) / lambda
R_max = R(1)-log(epsilon(1.0_dp))/lambda
do while (R(i_max) > R_max)
    if (i_max == 2) then
        call stop_error("Can't start the inward integration")
    end if
    i_max = i_max - 1
end do

u1(i_max:i_max+4) = exp(-lambda * (R(i_max:i_max+4)-R(1)))
u2(i_max:i_max+4) = - lambda * u1(i_max:i_max+4)
u1p(i_max:i_max+4) = Rp(i_max:i_max+4)                * u2(i_max:i_max+4)
u2p(i_max:i_max+4) = Rp(i_max:i_max+4) * C(i_max:i_max+4) * u1(i_max:i_max+4)

do i = i_max, 2, -1
    u1p(i) = Rp(i)        * u2(i)
    u2p(i) = Rp(i) * C(i) * u1(i)
    u1_tmp  = u1(i) + adams_interp_inward_implicit(u1p, i)
    u2_tmp  = u2(i) + adams_interp_inward_implicit(u2p, i)

    lambda = -9.0_dp / 24
    Delta = 1 - lambda**2 * C(i-1) * Rp(i-1)**2
    M(1, 1) = 1 / Delta
    M(2, 1) = lambda * C(i-1) * Rp(i-1) / Delta
    M(1, 2) = lambda * Rp(i-1) / Delta
    M(2, 2) = 1 / Delta

    u1(i-1) = M(1, 1) * u1_tmp + M(1, 2) * u2_tmp
    u2(i-1) = M(2, 1) * u1_tmp + M(2, 2) * u2_tmp
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
P(i_max+4:) = 0
Q(i_max+4:) = 0
imin = 1
end subroutine

subroutine integrate_rschroed_rk4(l, Z, E, R, V, Vmid, P, Q, imax)
! Integrates the Schrodinger equation outward using Runge-Kutta 4th order method
!
! It integrates the Schroedinger equation in the R(r) form outwards:
! R'' = -2/R * R' + (2*(V-E) + l*(l+1)/R**2) * R
!
! Returns P(r), Q(r), where these are defined as:
! P(r) = r*R(r)
! Q(r) = P'(r) = r * R'(r) + R(r)
! where R(r) is the radial solution.
integer, intent(in) :: l
integer, intent(in) :: Z
real(dp), intent(in) :: E
real(dp), intent(in) :: R(:)
real(dp), intent(in) :: V(:), Vmid(:)
real(dp), intent(out) :: P(:), Q(:)
integer, intent(out) :: imax ! The integration was carried to R(imax)

real(dp), parameter :: max_val = 1e6_dp
real(dp) :: y0(2), y1(size(R)), y2(size(R))
real(dp), dimension(size(R)) :: C1, C2
real(dp), dimension(size(R)-1) :: C1mid, C2mid
real(dp), dimension(size(R)-1) :: Rmid

! y(:) are values of the components (2) of the equation, in our case:
! y(1) = R, y(2) = R'
! dydx(:) are the derivatives of the components (2) of the equation

if (l == 0) then
    y0(1) = 1-Z*R(1)
    y0(2) = -Z
else
    y0(1) = R(1)**l
    y0(2) = l*R(1)**(l-1)
end if
if (size(V) /= size(Vmid) + 1) call stop_error("Vmid size is wrong")

C1 = 2*(V-E) + l*(l+1)/R**2
C2 = -2/R
Rmid = (R(:size(R)-1) + R(2:)) / 2
C1mid = 2*(Vmid-E) + l*(l+1)/Rmid**2
C2mid = -2/Rmid
call rk4_integrate(R, y0, C1, C2, C1mid, C2mid, max_val, y1, y2, imax)
P(:imax) = y1(:imax)*R(:imax) ! P(r) = r * R(r)
Q(:imax) = y2(:imax)*R(:imax) + y1(:imax) ! Q(r) = P'(r) = r * R'(r) + R(r)
end subroutine

end module
