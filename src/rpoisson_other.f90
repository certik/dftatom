module rpoisson_other

! Other Poisson integrators, not directly used by dftatom, but available
! for reuse. This module contains various Adams predictor-corrector integrators
! (both for outward and inward integration).

use types, only: dp
use constants, only: pi
use ode1d_other, only: adams_extrapolation_inward, adams_interp_inward
use ode1d, only: adams_extrapolation_outward, adams_interp_outward

implicit none

private
public rpoisson_inward, rpoisson_outward_pc_rV, rpoisson_inward_pc, &
    rpoisson_inward_pc_rV

contains

function rpoisson_inward(R, Rp, rho, Z) result(V)
! Solves the equation V''(r) + 2/r*V'(r) = -4*pi*rho
!
! Uses explicit Adams method.
!
!
! It rewrites it to (r*V)'' = -4*pi*rho*r, converts to uniform grid:
!   u1 = r*V
!   u2 = (r*V)'
!   u1p = u2 * Rp
!   u2p = -4*pi*rho*R * Rp
! and integrates inward using Adams method. This gives u1=r*V and
! at the end we divide by "r": V = u1/r
!
! The boundary condition is such, that V -> Z/r for large r. So the boundary
! condition for u1=r*V is Z. We just set the initial 4 points of u1 for
! the Adams method equal to Z (and zero derivative u2) at the right hand side
! of the domain.
real(dp), intent(in) :: R(:), Rp(:), rho(:)
integer, intent(in) :: Z
real(dp) :: V(size(R))

real(dp), dimension(size(R)) :: u1, u1p, u2p
real(dp) :: u2
integer :: N, i

N = size(R)
u1(N-3:N) = Z
u2 = 0
u1p(N-3:N) = 0
u2p = -4*pi*rho * R * Rp

do i = N-3, 2, -1
    u1p(i) = Rp(i) * u2
    u1(i-1)  = u1(i)  + adams_extrapolation_inward(u1p, i)
    u2  = u2  + adams_extrapolation_inward(u2p, i)
end do
V = u1/R
end function

function rpoisson_inward_pc_rV(R, Rp, rho, Z) result(V)
! Solves the equation V''(r) + 2/r*V'(r) = -4*pi*rho
!
! Uses predictor corrector Adams method.
!
! It rewrites it to (r*V)'' = -4*pi*rho*r, converts to uniform grid:
!   u1 = r*V
!   u2 = (r*V)'
!   u1p = u2 * Rp
!   u2p = -4*pi*rho*R * Rp
! and integrates inward using Adams method. This gives u1=r*V and
! at the end we divide by "r": V = u1/r
!
! The boundary condition is such, that V -> Z/r for large r. So the boundary
! condition for u1=r*V is Z. We just set the initial 4 points of u1 for
! the Adams method equal to Z (and zero derivative u2) at the right hand side
! of the domain.
real(dp), intent(in) :: R(:), Rp(:), rho(:)
integer, intent(in) :: Z
real(dp) :: V(size(R))

real(dp), dimension(size(R)) :: u1, u2, u1p, u2p
integer :: N, i, it
integer, parameter :: max_it = 4
!real(dp), parameter :: alpha = 32.388_dp
!real(dp), parameter :: alpha = 1.0_dp/25

N = size(R)
!u1(N-3:N) = Z * erf(alpha * R(N-3:N))
!u2(N-3:N) = Z * alpha * 2 / sqrt(pi) * exp(-alpha**2 * R(N-3:N)**2)
u1 = Z
u2 = 0
!u1 = Z * erf(alpha * R)
!u2 = Z * alpha * 2 / sqrt(pi) * exp(-alpha**2 * R**2)
!u1p(N-3:N) = u2(N-3:N) * Rp(N-3:N)
u1p = u2 * Rp
u2p = -4*pi*rho * R * Rp

do i = N-3, 2, -1
    !u1p(i) = Rp(i) * u2(i)
    u1(i-1)  = u1(i)  + adams_extrapolation_inward(u1p, i)
    u2(i-1)  = u2(i)  + adams_extrapolation_inward(u2p, i)
    do it = 1, max_it
        u1p(i-1) = Rp(i-1) * u2(i-1)
        u1(i-1)  = u1(i)  + adams_interp_inward(u1p, i)
        u2(i-1)  = u2(i)  + adams_interp_inward(u2p, i)
    end do
end do
V = u1/R
end function

function rpoisson_inward_pc(R, Rp, rho, Z) result(V)
! Solves the equation V''(r) + 2/r*V'(r) = -4*pi*rho
!
! Uses predictor corrector Adams method.
!
! It rewrites it to (r*V)'' = -4*pi*rho*r, converts to uniform grid:
!   u1 = r*V
!   u2 = (r*V)'
!   u1p = u2 * Rp
!   u2p = -4*pi*rho*R * Rp
! and integrates inward using Adams method. This gives u1=r*V and
! at the end we divide by "r": V = u1/r
!
! The boundary condition is such, that V -> Z/r for large r. So the boundary
! condition for u1=r*V is Z. We just set the initial 4 points of u1 for
! the Adams method equal to Z (and zero derivative u2) at the right hand side
! of the domain.
real(dp), intent(in) :: R(:), Rp(:), rho(:)
integer, intent(in) :: Z
real(dp) :: V(size(R))

real(dp), dimension(size(R)) :: u1, u2, u1p, u2p
integer :: N, i, it
integer, parameter :: max_it = 4
!real(dp), parameter :: alpha = 32.388_dp
!real(dp), parameter :: alpha = 1.0_dp/25

N = size(R)
u1 = Z
u2 = 0
!u1(N-3:N) = Z * erf(alpha * R(N-3:N))
!u2(N-3:N) = Z * alpha * 2 / sqrt(pi) * exp(-alpha**2 * R(N-3:N)**2)
!u1 = Z * erf(alpha * R) / R
!u2 = -Z * erf(alpha * R) / R**2 + &
!            Z * alpha * 2 / sqrt(pi) * exp(-alpha**2 * R**2) / R
u1p = u2 * Rp
u2p = -(4*pi*rho + 2*u2/R) * Rp

do i = N-3, 2, -1
    u1(i-1)  = u1(i)  + adams_extrapolation_inward(u1p, i)
    u2(i-1)  = u2(i)  + adams_extrapolation_inward(u2p, i)
    do it = 1, max_it
        u1p(i-1) = +Rp(i-1) * u2(i-1)
        u2p(i-1) = -Rp(i-1) * (4*pi*rho(i-1) + 2*u2(i-1)/R(i-1))
        u1(i-1)  = u1(i)  + adams_interp_inward(u1p, i)
        u2(i-1)  = u2(i)  + adams_interp_inward(u2p, i)
    end do
end do
V = u1
end function

function rpoisson_outward_pc_rV(R, Rp, rho, Z) result(V)
! Solves the equation V''(r) + 2/r*V'(r) = -4*pi*rho
!
! Uses predictor corrector Adams method.
!
! It rewrites it to (r*V)'' = -4*pi*rho*r, converts to uniform grid:
!   u1 = r*V
!   u2 = (r*V)'
!   u1p = u2 * Rp
!   u2p = -4*pi*rho*R * Rp
! and integrates outward using Adams method. This gives u1=r*V and
! at the end we divide by "r": V = u1/r
real(dp), intent(in) :: R(:), Rp(:), rho(:)
integer, intent(in) :: Z
real(dp) :: V(size(R))

real(dp), dimension(size(R)) :: u1, u2, u1p, u2p
integer :: N, i, it
integer, parameter :: max_it = 2
!real(dp), parameter :: alpha = 32.388_dp
!real(dp), parameter :: alpha = 1.0_dp / 25

N = size(R)
!u1(:4) = Z * erf(alpha * R(:4))
!u2(:4) = Z * alpha * 2 / sqrt(pi) * exp(-alpha**2 * R(:4)**2)
!u1(:4) = Z * 2 * alpha * R(:4) / sqrt(pi)
!u2(:4) = Z * 2 * alpha / sqrt(pi)
u1(:4) = Z
u2(:4) = 0
!print *, u1(:4)
!print *, u2(:4)
!stop "OK"
u1p(:4) = u2(:4) * Rp(:4)
u2p = -4*pi*rho * R * Rp

do i = 4, N-1
    u1p(i) = Rp(i) * u2(i)
    u1(i+1)  = u1(i)  + adams_extrapolation_outward(u1p, i)
    u2(i+1)  = u2(i)  + adams_extrapolation_outward(u2p, i)
    do it = 1, max_it
        u1p(i+1) = Rp(i+1) * u2(i+1)
        u1(i+1)  = u1(i)  + adams_interp_outward(u1p, i)
        u2(i+1)  = u2(i)  + adams_interp_outward(u2p, i)
    end do
end do
V = u1/R
!V = V - V(N) + Z/R(N)
end function

end module
