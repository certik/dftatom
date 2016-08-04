program test_poisson_kernels
use utils, only: assert
use types, only: dp
use mesh, only: mesh_exp, mesh_exp_deriv
use rpoisson, only: rpoisson_kernel1, rpoisson_kernel2
implicit none

! Mesh parameters:
real(dp), parameter :: r_min = 1e-7_dp, r_max = 10, a = 2.7e6_dp
integer, parameter :: NN = 5000

real(dp) :: R(NN+1)
real(dp), dimension(size(R)) :: Rp, rho
real(dp), dimension(size(R), 2) :: u1, u2, u1p, u2p

R = mesh_exp(r_min, r_max, a, NN)
Rp = mesh_exp_deriv(r_min, r_max, a, NN)
rho = exp(-R)
u1 = 1
u2 = 1
u1p = 0
u2p = 0
call rpoisson_kernel1(R, Rp, rho, u1(:, 1), u2(:, 1), u1p(:, 1), u2p(:, 1))
call rpoisson_kernel2(R, Rp, rho, u1(:, 2), u2(:, 2), u1p(:, 2), u2p(:, 2))
print *, "Errors:"
print *, maxval(abs(u1(:,1)-u1(:,2)))
print *, maxval(abs(u2(:,1)-u2(:,2)))
print *, maxval(abs(u1p(:,1)-u1p(:,2)))
print *, maxval(abs(u2p(:,1)-u2p(:,2)))
call assert(maxval(abs(u1(:,1)-u1(:,2))) < 1e-10_dp)
call assert(maxval(abs(u2(:,1)-u2(:,2))) < 1e-15_dp)
call assert(maxval(abs(u1p(:,1)-u1p(:,2))) < 1e-9_dp)
call assert(maxval(abs(u2p(:,1)-u2p(:,2))) < 1e-16_dp)
end program
