module mixings

! This module contains SCF mixing algorithms.

use types, only: dp
use constants, only: pi
use dft_data, only: dft_data_t
use utils, only: stop_error
use ode1d, only: integrate
implicit none
private
public mixing_linear, mixing_linear_adapt, mixing_anderson

contains

function mixing_linear(R, x0, max_iter, alpha, d)
! Finds "x" so that R(x) = 0, uses x0 as the initial estimate
real(dp), intent(in) :: x0(:), alpha
integer, intent(in) :: max_iter
type(dft_data_t), intent(inout) :: d ! Data passed to "F"
real(dp) :: mixing_linear(size(x0))
interface
    function R(x, i, d)
    use types
    use dft_data
    implicit none
    real(dp), intent(in) :: x(:) ! "x"
    integer, intent(in) :: i ! iteration #
    type(dft_data_t), intent(inout) :: d ! F's data
    real(dp) :: R(size(x))
    end function
end interface

real(dp), dimension(size(x0)) :: x_i
integer :: i
x_i = x0
do i = 1, max_iter
    x_i = x_i + alpha * R(x_i, i, d)
end do
mixing_linear = x_i
end function

function mixing_linear_adapt(R, x0, max_iter, alpha, d)
! Finds "x" so that F(x) = 0, uses x0 as the initial estimate
real(dp), intent(in) :: x0(:)
type(dft_data_t), intent(inout) :: d ! Data passed to "F"
integer, intent(in) :: max_iter
real(dp) :: mixing_linear_adapt(size(x0)), alpha
interface
    function R(x, i, d)
    use types
    use dft_data
    implicit none
    real(dp), intent(in) :: x(:) ! "x"
    integer, intent(in) :: i ! iteration #
    type(dft_data_t), intent(inout) :: d ! F's data
    real(dp) :: R(size(x))
    end function
end interface

real(dp), parameter :: alpha_max = 1.0_dp ! mixing parameter
real(dp), dimension(size(x0)) :: x_m, R_m, R_mm1, beta
integer :: i, j
x_m = x0
beta(:) = alpha
do i = 1, max_iter
    R_m = R(x_m, i, d)
    x_m = x_m + beta * R_m
    if (i > 1) then
        do j = 1, size(beta)
            if (R_mm1(j) * R_m(j) > 0) then
                beta(j) = beta(j) + alpha
                if (beta(j) > alpha_max) beta(j) = alpha_max
            else
                beta(j) = alpha
            end if
        end do
    end if
    R_mm1 = R_m
end do
mixing_linear_adapt = x_m
end function

function mixing_anderson(R, x0, max_iter, energy_crit, d, alpha, eps)
! Finds "x" so that R(x) = 0, uses x0 as the initial estimate
real(dp), intent(in) :: x0(:)
integer, intent(in) :: max_iter
logical, intent(in) :: energy_crit
type(dft_data_t), intent(inout) :: d ! Data passed to "F"
real(dp), intent(in) :: alpha
real(dp), intent(in) :: eps
real(dp) :: mixing_anderson(size(x0))
interface
    function R(x, i, d)
    use types
    use dft_data
    implicit none
    real(dp), intent(in) :: x(:) ! "x"
    integer, intent(in) :: i ! iteration #
    type(dft_data_t), intent(inout) :: d ! F's data
    real(dp) :: R(size(x))
    end function
end interface

real(dp), dimension(size(x0)) :: x_i, x1_i, R_i, R1_i, delta_R, delta_x
real(dp) :: beta
real(dp) :: sn, sd
real(dp) :: ks_energies(size(d%ks_energies))
real(dp) :: x_i_norm, R_i_norm
real(dp) :: err_old, err, L2_err
integer :: i
x_i = x0
if (energy_crit) then
    ks_energies = d%ks_energies
    err_old = 1e12_dp
end if
do i = 1, max_iter
    R_i = R(x_i, i, d)
    if (energy_crit) then
        ! L2 norm of the "input" potential:
        x_i_norm = sqrt(4*pi*integrate(d%Rp, x_i**2 * d%R**2))
        ! L2 norm of the "output-input" potential:
        R_i_norm = sqrt(4*pi*integrate(d%Rp, R_i**2 * d%R**2))
        if (x_i_norm < 1e-12_dp) x_i_norm = 1e-12_dp
        L2_err = R_i_norm / x_i_norm
        err = maxval(abs(ks_energies - d%ks_energies))
        !print *, i, L2_err, err
        ! Do at least 3 iterations
        if (i >= 3 .and. L2_err < 5e-5_dp) then
            if (err < eps .and. err_old < eps) then
                mixing_anderson = x_i
                return
            end if
        end if
        ks_energies = d%ks_energies
        err_old = err
    end if

    if (i > 1) then
        delta_x = x_i - x1_i
        delta_R = R_i - R1_i
    end if
    x1_i = x_i
    R1_i = R_i
    x_i = x_i + alpha * R_i
    if (i > 1) then
        sn = sum(R_i * delta_R * d%R**2)
        sd = sum(delta_R**2 * d%R**2)
        beta = sn / sd
        x_i = x_i - beta * (delta_x + alpha * delta_R)
    end if
end do
mixing_anderson = x_i
if (energy_crit) call stop_error("SCF didn't converge")
end function

end module
