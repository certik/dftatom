module mixing_r

use types
use dft_data
use mixings
implicit none
private
public R1, R2, R3

contains

function R1(x, i, tmp)
real(dp), intent(in) :: x(:)
integer, intent(in) :: i
type(dft_data_t), intent(inout) :: tmp
real(dp) :: R1(size(x))

real(dp) :: d(size(x)), c
print *, i, associated(tmp%ks_energies)

d = (/ 3.0_dp, 2.0_dp, 1.5_dp, 1.0_dp, 0.5_dp /)
c = 0.01_dp
R1 = -d*x - c*sum(x*x)*x
end function

function R2(x, i, tmp)
real(dp), intent(in) :: x(:)
integer, intent(in) :: i
type(dft_data_t), intent(inout) :: tmp
real(dp) :: R2(size(x))
print *, i, associated(tmp%ks_energies)

R2 = cos(x) - x
end function

function R3(x, i, tmp)
real(dp), intent(in) :: x(:)
integer, intent(in) :: i
type(dft_data_t), intent(inout) :: tmp
real(dp) :: R3(size(x))

real(dp) :: a, b
real(dp), parameter :: c = 2.0_dp
print *, i, associated(tmp%ks_energies)
a = x(1)
b = x(2)
R3 = a**2 + b**2 - c**2
end function

end module

program mixings_test
! Tests all mixing schemes, that they work as excpected.
use types
use mixings
use dft_data
use mixing_r
use utils
implicit none

type(dft_data_t) :: d
real(dp) :: x0(5), x(5), y(50), y0(50), z(2), z0(2)
real(dp), target :: tmp1(5), tmp2(50), tmp3(2)
integer :: i

x0(:) = 1.0_dp

! For anderson mixing:
tmp1(:) = 1.0_dp
d%R => tmp1

x = mixing_linear(R1, x0, 60, 0.5_dp, d)
call assert(sqrt(sum(x**2)) < 1e-7_dp)
call assert(sqrt(sum(R1(x, 1, d)**2)) < 1e-7_dp)

x = mixing_linear_adapt(R1, x0, 20, 0.5_dp, d)
call assert(sqrt(sum(x**2)) < 1e-5_dp)
call assert(sqrt(sum(R1(x, 1, d)**2)) < 1e-5_dp)

x = mixing_anderson(R1, x0, 35, .false., d, alpha=0.5_dp, eps=0.0_dp)
call assert(sqrt(sum(x**2)) < 1e-5_dp)
call assert(sqrt(sum(R1(x, 1, d)**2)) < 1e-5_dp)

! ***************************************************

do i = 1, size(y)
    y0(i) = i
end do

! For anderson mixing:
tmp2(:) = 1.0_dp
d%R => tmp2

y = mixing_linear(R2, y0, 50, 1.0_dp, d)
call assert(sum((y-0.73908513321516067)**2) < 1e-10_dp)

y = mixing_linear_adapt(R2, y0, 20, 0.1_dp, d)
call assert(sum((y-0.73908513321516067)**2) < 1e-10_dp)

y = mixing_anderson(R2, y0, 9, .false., d, alpha=0.5_dp, eps=0.0_dp)
call assert(sum((y-0.73908513321516067)**2) < 1e-10_dp)

! ***************************************************

z0 = (/ 0.5, 0.5 /)

! For anderson mixing:
tmp3(:) = 1.0_dp
d%R => tmp3

z = mixing_linear(R3, z0, 40, 0.1_dp, d)
call assert(abs(z(1)**2 + z(2)**2 - 2.0_dp**2) < 1e-10_dp)

z = mixing_linear_adapt(R3, z0, 30, 0.1_dp, d)
call assert(abs(z(1)**2 + z(2)**2 - 2.0_dp**2) < 1e-10_dp)

z = mixing_anderson(R3, z0, 8, .false., d, alpha=0.5_dp, eps=0.0_dp)
call assert(abs(z(1)**2 + z(2)**2 - 2.0_dp**2) < 1e-10_dp)

end program
