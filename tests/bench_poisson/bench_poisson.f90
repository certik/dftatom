program bench_poisson
use types, only: dp
use mesh, only: mesh_exp, mesh_exp_deriv
use rpoisson, only: rpoisson_outward_pc
implicit none

! Mesh parameters:
real(dp), parameter :: r_min = 1e-7_dp, r_max = 10, a = 2.7e6_dp
integer, parameter :: NN = 5000

real(dp) :: R(NN+1)
real(dp), dimension(size(R)) :: Rp, u1, u2, u1p, u2p, rho

R = mesh_exp(r_min, r_max, a, NN)
Rp = mesh_exp_deriv(r_min, r_max, a, NN)
rho = exp(-R)
u1 = rpoisson_outward_pc(R, Rp, rho)
end program
