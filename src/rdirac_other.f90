module rdirac_other

! Other Dirac integrators, not directly used by dftatom, but available
! for reuse. This module contains various Adams predictor-corrector integrators
! (both for outward and inward integration) and functions to calculate analytic
! asymptotic.

use types, only: dp
use ode1d, only: adams_extrapolation_outward, adams_interp_outward
use ode1d_other, only: adams_extrapolation_inward, adams_interp_inward
use utils, only: stop_error

implicit none
private
public dirac_outward_adams_pr, dirac_inward_adams_pr

contains

subroutine dirac_outward_adams_pr(c, kappa, Z, E, R, Rp, V, P, Q, imax)
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

real(dp) :: beta
real(dp), dimension(size(R)) :: u1, u2, u1p, u2p
integer :: i, it

beta = sqrt(kappa**2-(Z/c)**2)
u1(:4) = R(:4)**(beta)
u2(:4) = R(:4)**(beta-1)*(beta+kappa)/((E-V(:4))/c+2*c)

Ctot(:, 1, 1) = -kappa / R
Ctot(:, 2, 2) = +kappa / R
Ctot(:, 1, 2) = +(E-V)/c + 2*c
Ctot(:, 2, 1) = -(E-V)/c

u1p(:4) = Rp(:4) * (Ctot(:4, 1, 1) * u1(:4) + Ctot(:4, 1, 2) * u2(:4))
u2p(:4) = Rp(:4) * (Ctot(:4, 2, 1) * u1(:4) + Ctot(:4, 2, 2) * u2(:4))

do i = 4, size(R)-1
    u1(i+1)  = u1(i)  + adams_extrapolation_outward(u1p, i)
    u2(i+1)  = u2(i)  + adams_extrapolation_outward(u2p, i)
    do it = 1, 2
        u1p(i+1) = Rp(i+1) * (Ctot(i+1, 1, 1)*u1(i+1) + Ctot(i+1, 1, 2)*u2(i+1))
        u2p(i+1) = Rp(i+1) * (Ctot(i+1, 2, 1)*u1(i+1) + Ctot(i+1, 2, 2)*u2(i+1))
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

subroutine dirac_inward_adams_pr(c, kappa, E, R, Rp, V, P, Q, imin)
!integrates the Dirac eq. inwards, returns r*R, where R is the radial
!solution.
! input parameters:
! kappa .... the kappa in the Dirac eq.
! Z .... the nucleus charge in Hartree atomic units
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
integer :: i, it

nr = size(R)

if (E > 0) call stop_error("E < 0 required")

! Find the starting point 'i_max' for the inward integration just like we do
! for the Schroedinger equation, as it works excellent for the Dirac equation
! as well:
i_max = nr-4
if (i_max < 2) call stop_error("size(R) too small to start inward integraion")
lambda = sqrt(-2*E)
do while (lambda*R(i_max) > 40)
    if (i_max == 2) then
        print *, "E =", E, "lambda =", lambda
        call stop_error("Can't start the inward integration")
    end if
    i_max = i_max - 1
end do

u1(i_max+1:) = 0
u2(i_max+1:) = 0

! This simple asymptotic gives accurate energies. One can use the
! get_asymptotic() function for a very precise asymptotic, but it doesn't seem
! to make any difference in final energies.
u1(i_max:i_max+4) = 1e-12_dp
u2(i_max:i_max+4) = 1e-12_dp

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
    u1(i-1)  = u1(i)  + adams_extrapolation_inward(u1p, i)
    u2(i-1)  = u2(i)  + adams_extrapolation_inward(u2p, i)
    do it = 1, 2
        u1p(i-1) = Rp(i-1) * (Ctot(i-1, 1, 1)*u1(i-1) + Ctot(i-1, 1, 2)*u2(i-1))
        u2p(i-1) = Rp(i-1) * (Ctot(i-1, 2, 1)*u1(i-1) + Ctot(i-1, 2, 2)*u2(i-1))
        u1(i-1)  = u1(i) + adams_interp_inward(u1p, i)
        u2(i-1)  = u2(i) + adams_interp_inward(u2p, i)
    end do
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

end module
