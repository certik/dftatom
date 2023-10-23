module rpoisson

! Routines in this module solve the radial Poisson equation outward using
! the predictor-corrector method (with Adams extrapolation/interpolation).

use types, only: dp
use utils, only: stop_error
use constants, only: pi
use ode1d, only: adams_extrapolation_outward, adams_interp_outward
use ode1d, only: integrate, get_midpoints, rk4_integrate3

implicit none

private
public rpoisson_outward_pc

contains

function rpoisson_outward_pc(R, Rp, rho) result(V)
! Solves the equation V''(r) + 2/r*V'(r) = -4*pi*rho
!
! Uses predictor corrector Adams method.
!
! It rewrites it to the equivalent system of first order ODEs on a uniform
! grid:
!   u1 = V
!   u2 = V'
!   u1p = u2 * Rp
!   u2p = -(4*pi*rho + 2*u2/r) * Rp
! and integrates outward using Adams method. The initial conditions are:
!   V (R(1)) = u1(1) = 4*pi * \int r * rho(r) dr
!   V'(R(1)) = u2(1) = 0
real(dp), intent(in) :: R(:), Rp(:), rho(:)
real(dp) :: V(size(R))

real(dp), dimension(size(R)) :: u1, u2, u1p, u2p
integer :: N, i, it
integer, parameter :: max_it = 2
real(dp) :: rho_mid(3)

N = size(R)
rho_mid = get_midpoints(R(:4), rho(:4))
call rpoisson_outward_rk4(rho(:4), rho_mid, R(:4), &
    4*pi*integrate(Rp, rho*R), &
    0.0_dp, &
    u1(:4), u2(:4))

u1p(:4) = u2(:4) * Rp(:4)
u2p(:4) = -(4*pi*rho(:4) + 2*u2(:4)/R(:4)) * Rp(:4)

do i = 4, N-1
    u1(i+1)  = u1(i)  + adams_extrapolation_outward(u1p, i)
    u2(i+1)  = u2(i)  + adams_extrapolation_outward(u2p, i)
    do it = 1, max_it
        u1p(i+1) = +Rp(i+1) * u2(i+1)
        u2p(i+1) = -Rp(i+1) * (4*pi*rho(i+1) + 2*u2(i+1)/R(i+1))
        u1(i+1)  = u1(i)  + adams_interp_outward(u1p, i)
        u2(i+1)  = u2(i)  + adams_interp_outward(u2p, i)
    end do
end do
V = u1
end function

subroutine rpoisson_outward_rk4(density, density_mid, R, V0, V0d, V, Vd)
! Solves V''(r) + 2/r*V'(r) = -4*pi*density (in this form) with initial
! conditions V(R(1)) = V0 and V'(R(1)) = V0d using 4th order Runge-Kutta method
! by rewriting it to the equivalent first order ODEs (see the documentation
! of rk4_integrate3() for more details). Returns V (value) and Vd (derivative).
real(dp), intent(in), dimension(:) :: density, density_mid ! density at the
    ! mesh points R and at midpoints
real(dp), intent(in), dimension(:) :: R ! radial grid
real(dp), intent(in) :: V0  ! value of V(r) at r=R(1)
real(dp), intent(in) :: V0d ! value of V'(r) at r=R(1)
real(dp), intent(out), dimension(:) :: V, Vd  ! V(r) and V'(r) on the grid R

real(dp) :: y(2)
real(dp), dimension(size(R)) :: C1, C2
real(dp), dimension(size(R)-1) :: C1mid, C2mid
real(dp), dimension(size(R)-1) :: Rmid
integer :: imax

y(1) = V0
y(2) = V0d

C1 = -4*pi*density
C2 = -2 / R
Rmid = (R(:size(R)-1) + R(2:)) / 2
C1mid = -4*pi*density_mid
C2mid = -2 / Rmid
call rk4_integrate3(R, y, C1, C2, C1mid, C2mid, 1e10_dp, V, Vd, imax)
if (imax /= size(R)) call stop_error("Poisson solver diverged.")
end subroutine

end module
