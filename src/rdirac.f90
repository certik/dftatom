module rdirac

! Routines in this module solve the radial Dirac equation outward and
! inward using the implicit Adams method.

use types, only: dp
use ode1d, only: adams_interp_outward_implicit, &
    adams_interp_inward_implicit, get_midpoints, rk4_integrate4
use utils, only: stop_error

implicit none
private
public dirac_outward_adams, dirac_inward_adams


contains

subroutine dirac_outward_adams(c, kappa, Z, E, R, Rp, V, P, Q, imax)
!integrates the Dirac eq, returns r*R, where R is the radial
!solution.
! input parameters:
! kappa .... the kappa in the Dirac eq.
! Z .... the nucleus charge in Hartree atomic units
! E .... the energy at which to integration the equation
! R .... radial grid
! V .... potential on the radial grid
! c .... speed of light
! output parameters:
! Q .... f-component of the radial dirac wave function
! P .... g-component of the radial dirac wave function
real(dp), intent(in) :: c
integer, intent(in) :: kappa
integer, intent(in) :: Z
real(dp), intent(in) :: E
real(dp), intent(in) :: R(:), Rp(:)
real(dp), intent(in) :: V(:)
real(dp), intent(out) :: Q(:)
real(dp), intent(out) :: P(:)
integer, intent(out) :: imax

real(dp), parameter :: max_val = 1e6_dp
real(dp), dimension(size(R), 2, 2) :: Ctot

real(dp), dimension(size(R)) :: u1, u2, u1p, u2p
integer :: i
real(dp) :: lambda, Delta, M(2, 2), u1_tmp, u2_tmp
real(dp) :: Vmid(3)

if (size(R) < 4) call stop_error("size(R) <= 4")
Vmid = get_midpoints(R(:4), V(:4))
call integrate_radial_dirac_r_rk4(c, kappa, Z, E, R(:4), V(:4), Vmid, &
    u1(:4), u2(:4), imax)
if (imax /= 4) call stop_error("rk4 failed")

Ctot(:, 1, 1) = -kappa / R
Ctot(:, 2, 2) = +kappa / R
Ctot(:, 1, 2) = +(E-V)/c + 2*c
Ctot(:, 2, 1) = -(E-V)/c

u1p(:4) = Rp(:4) * (Ctot(:4, 1, 1) * u1(:4) + Ctot(:4, 1, 2) * u2(:4))
u2p(:4) = Rp(:4) * (Ctot(:4, 2, 1) * u1(:4) + Ctot(:4, 2, 2) * u2(:4))

do i = 4, size(R)-1
    u1p(i) = Rp(i) * (Ctot(i, 1, 1)*u1(i) + Ctot(i, 1, 2)*u2(i))
    u2p(i) = Rp(i) * (Ctot(i, 2, 1)*u1(i) + Ctot(i, 2, 2)*u2(i))
    u1_tmp  = u1(i) + adams_interp_outward_implicit(u1p, i)
    u2_tmp  = u2(i) + adams_interp_outward_implicit(u2p, i)

    lambda = 9.0_dp / 24
    Delta = 1 - lambda**2 * Rp(i+1)**2 * (Ctot(i+1, 1, 2) * Ctot(i+1, 2, 1) &
        -Ctot(i+1, 1, 1) * Ctot(i+1, 2, 2))
    M(1, 1) = (1 - lambda * Rp(i+1) * Ctot(i+1, 2, 2)) / Delta
    M(2, 1) = lambda * Rp(i+1) * Ctot(i+1, 2, 1) / Delta
    M(1, 2) = lambda * Rp(i+1) * Ctot(i+1, 1, 2) / Delta
    M(2, 2) = (1 - lambda * Rp(i+1) * Ctot(i+1, 1, 1)) / Delta

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

subroutine dirac_inward_adams(c, kappa, E, R, Rp, V, P, Q, imin)
!integrates the Dirac eq. inwards, returns r*R, where R is the radial
!solution.
! input parameters:
! kappa .... the kappa in the Dirac eq.
! E .... the energy at which to integration the equation
! R .... radial grid
! V .... potential on the radial grid
! c .... speed of light
! output parameters:
! f .... f-component of the radial dirac wave function
! P .... g-component of the radial dirac wave function
real(dp), intent(in) :: c
integer, intent(in) :: kappa
real(dp), intent(in) :: E
real(dp), intent(in) :: R(:), Rp(:)
real(dp), intent(in) :: V(:)
real(dp), intent(out) :: Q(:)
real(dp), intent(out) :: P(:)
integer, intent(out) :: imin

integer :: nr
integer :: i_max
real(dp), parameter :: max_val = 1e20_dp
real(dp) :: lambda
real(dp), dimension(size(R), 2, 2) :: Ctot
real(dp), dimension(size(R)) :: u1, u2, u1p, u2p
integer :: i
real(dp) :: Delta, M(2, 2), u1_tmp, u2_tmp
real(dp) :: R_max

nr = size(R)

if (E > 0) call stop_error("E < 0 required")

i_max = nr-4
if (i_max < 2) call stop_error("size(R) too small to start inward integraion")
lambda = sqrt(-2*E-E**2/c**2)
! We require that exp(-lambda*(R-R(1)) ~ epsilon(1.0_dp),
! if we start further from
! the origin, it might sometimes blow up, if we start closer, we might not get
! as precise asymptotic.
! It follows that R ~ R(1)-log(epsilon(1.0_dp)) / lambda
R_max = R(1)-log(epsilon(1.0_dp))/lambda
do while (R(i_max) > R_max)
    if (i_max == 2) then
        print *, "E =", E, "lambda =", lambda
        call stop_error("Can't start the inward integration")
    end if
    i_max = i_max - 1
end do

u1(i_max+4:) = 0
u2(i_max+4:) = 0

u1(i_max:i_max+4) = exp(-lambda * (R(i_max:i_max+4)-R(1))) / sqrt(-E/(E+2*c**2))
u2(i_max:i_max+4) = -exp(-lambda * (R(i_max:i_max+4)-R(1)))

Ctot(:, 1, 1) = -kappa / R
Ctot(:, 2, 2) = +kappa / R
Ctot(:, 1, 2) = +(E-V)/c + 2*c
Ctot(:, 2, 1) = -(E-V)/c

u1p(i_max:i_max+4) = Rp(i_max:i_max+4) * &
    (Ctot(i_max:i_max+4, 1, 1)*u1(i_max:i_max+4) &
        + Ctot(i_max:i_max+4, 1, 2) * u2(i_max:i_max+4))
u2p(i_max:i_max+4) = Rp(i_max:i_max+4) * &
    (Ctot(i_max:i_max+4, 2, 1)*u1(i_max:i_max+4) &
        + Ctot(i_max:i_max+4, 2, 2) * u2(i_max:i_max+4))

do i = i_max, 2, -1
    u1p(i) = Rp(i) * (Ctot(i, 1, 1)*u1(i) + Ctot(i, 1, 2)*u2(i))
    u2p(i) = Rp(i) * (Ctot(i, 2, 1)*u1(i) + Ctot(i, 2, 2)*u2(i))
    u1_tmp  = u1(i) + adams_interp_inward_implicit(u1p, i)
    u2_tmp  = u2(i) + adams_interp_inward_implicit(u2p, i)

    lambda = -9.0_dp / 24
    Delta = 1 - lambda**2 * Rp(i-1)**2 * (Ctot(i-1, 1, 2) * Ctot(i-1, 2, 1) &
        -Ctot(i-1, 1, 1) * Ctot(i-1, 2, 2))
    M(1, 1) = (1 - lambda * Rp(i-1) * Ctot(i-1, 2, 2)) / Delta
    M(2, 1) = lambda * Rp(i-1) * Ctot(i-1, 2, 1) / Delta
    M(1, 2) = lambda * Rp(i-1) * Ctot(i-1, 1, 2) / Delta
    M(2, 2) = (1 - lambda * Rp(i-1) * Ctot(i-1, 1, 1)) / Delta

    u1(i-1) = M(1, 1) * u1_tmp + M(1, 2) * u2_tmp
    u2(i-1) = M(2, 1) * u1_tmp + M(2, 2) * u2_tmp
    if (abs(u1(i-1)) >= max_val .or. abs(u2(i-1)) >= max_val) then
        P = u1
        Q = u2
        imin = i
        return
    end if
end do
P = u1
Q = u2
imin = 1
end subroutine

real(dp) recursive function a_k(n, lambda, sigma, zeta, kappa, c) result(r)
integer, intent(in) :: n, kappa
real(dp), intent(in) :: lambda, sigma, c, zeta
if (n < 1) then
    r = 0
    call stop_error("a_k: n >= 1 required")
else
    r = c/(n*lambda) * (kappa + (n-sigma)*sigma*lambda/zeta - zeta*lambda/c**2)
    r = r * b_k(n, lambda, sigma, zeta, kappa, c)
end if
end function

real(dp) recursive function b_k(k, lambda, sigma, zeta, kappa, c) result(r)
integer, intent(in) :: k, kappa
real(dp), intent(in) :: lambda, sigma, c, zeta
integer :: n
if (k < 1) then
    r = 0
    call stop_error("b_k: n >=1 required")
else if (k == 1) then
    r = 1/(2*c) * (kappa + zeta/lambda)
else
    n = k-1
    r = 1/(2*n*lambda) * (kappa**2 - (n-sigma)**2 - zeta**2/c**2) * &
        b_k(n, lambda, sigma, zeta, kappa, c)
end if
end function

function get_asymptotic(r, E, c, kappa, zeta, n_terms, norm)
! Calculates the Dirac asymptotic with 'n_terms' terms
! n_terms ... 0, 1, 2, 3, ...
! With n_terms=0, we get the simplest asymptotic, with high n_terms, we get a
! very precise value
!
! Example (calculate 20 terms):
!    zeta = -V(i_max)*R(i_max)
!    Y = get_asymptotic(R(i_max), E, c, kappa, zeta, 20, .true.)
real(dp), intent(in) :: r, E, c, zeta
integer, intent(in) :: kappa, n_terms
logical, intent(in) :: norm ! .true. ... normalize the asymptotics
real(dp) :: get_asymptotic(2)

real(dp) :: lambda, sigma, a_term, b_term
real(dp), parameter :: norm_constant = 1e-12_dp
integer :: i
lambda = sqrt(c**2 - E**2/c**2)
sigma = E*zeta/(c**2*lambda)
a_term = 1
b_term = 0
do i = 1, n_terms
    a_term = a_term + a_k(i, lambda, sigma, zeta, kappa, c) / r**i
    b_term = b_term + b_k(i, lambda, sigma, zeta, kappa, c) / r**i
end do
a_term = a_term * sqrt((c**2+E)/(2*c**2))
b_term = b_term * sqrt((c**2-E)/(2*c**2))
get_asymptotic(1) = (a_term + b_term)
get_asymptotic(2) = (a_term - b_term)
if (norm) then
    ! Normalized, so that P(r) = norm_constant:
    get_asymptotic = get_asymptotic / get_asymptotic(1) * norm_constant
else
    ! Correct asymptotic form:
    get_asymptotic = get_asymptotic * r**sigma * exp(-lambda*r)
end if
end function

subroutine integrate_radial_dirac_r_rk4(c, kappa, Z, E, R, V, Vmid, P, Q, imax)
!integrates the Dirac eq, returns r*R, where R is the radial
!solution.
! input parameters:
! kappa .... the kappa in the Dirac eq.
! Z .... the nucleus charge in Hartree atomic units
! E .... the energy at which to integration the equation
! R .... radial grid
! V .... potential on the radial grid
! c .... speed of light
! output parameters:
! Q .... f-component of the radial dirac wave function
! P .... g-component of the radial dirac wave function
real(dp), intent(in) :: c
integer, intent(in) :: kappa
integer, intent(in) :: Z
real(dp), intent(in) :: E
real(dp), intent(in) :: R(:)
real(dp), intent(in) :: V(:), Vmid(:)
real(dp), intent(out) :: Q(:)
real(dp), intent(out) :: P(:)
integer, intent(out) :: imax

real(dp), parameter :: max_val = 1e6_dp
real(dp) :: y(2)
real(dp), dimension(size(R), 2, 2) :: Ctot
real(dp), dimension(size(R)-1, 2, 2) :: Cmid
real(dp), dimension(size(R)-1) :: Rmid

real(dp) :: beta, Z1
integer :: l

beta = sqrt(kappa**2-(Z/c)**2)
if (Z /= 0) then
    y(1) = R(1)**beta
    y(2) = R(1)**beta * c * (beta + kappa) / Z
else
    Z1 = V(1)
    if (kappa < 0) then
        l = -kappa-1
        y(1) = +R(1)**(l+1)
        y(2) = +R(1)**(l+2) * (E + Z1) / (c*(2*l+3))
    else
        l = kappa
        y(1) = -R(1)**(l+2) * (E + Z1) / (c*(2*l+1))
        y(2) = +R(1)**(l+1)
    end if
end if

Ctot(:, 1, 1) = -kappa / R
Ctot(:, 2, 2) = +kappa / R
Ctot(:, 1, 2) = +(E-V)/c + 2*c
Ctot(:, 2, 1) = -(E-V)/c

Rmid = (R(:size(R)-1) + R(2:)) / 2
Cmid(:, 1, 1) = -kappa / Rmid
Cmid(:, 2, 2) = +kappa / Rmid
Cmid(:, 1, 2) = +(E-Vmid)/c + 2*c
Cmid(:, 2, 1) = -(E-Vmid)/c
call rk4_integrate4(R, y, Ctot, Cmid, max_val, P, Q, imax)
end subroutine


end module
