program well
use special, only: spherical_bessel_jn_zeros
use dftatom
implicit none

! For potential well, we need to set Z=0:
integer, parameter :: Z = 0
! Mesh parameters:
real(dp), parameter :: r_min = 1e-7_dp, r_max = 10, a = 2.7e6_dp
integer, parameter :: NN = 5000

real(dp), parameter :: c = 137.035999037_dp, eps = 1e-12_dp
integer :: n, l, relat, converged
real(dp) :: r(NN+1), u(size(r)), Ein, E, E_exact, error, P(size(r)), Q(size(r))
real(dp) :: Rp(NN+1)
real(dp) :: Emin_init, Emax_init
real(dp), allocatable :: zeros(:, :)

R = mesh_exp(r_min, r_max, a, NN)
Rp = mesh_exp_deriv(r_min, r_max, a, NN)
u = 0

allocate(zeros(10, 0:7))
zeros = spherical_bessel_jn_zeros(7, 10, 1e-13_dp)
print *, "Infinite spherical potential well"
print *, "Mesh parameters (r_min, r_max, a, N):"
print "(ES10.2, F10.2, ES10.2, I10)", r_min, r_max, a, NN
print *
print "(A3, A3, A15, A15, A10)", "n", "l", "E", "E_exact", "Error"
print *
do n = 1, 7
    do l = 0, n-1
        relat = 0
        E_exact = (zeros(n-l, l) / r_max)**2 / 2
        Ein = 10
        Emin_init = 0
        Emax_init = 100
        call solve_radial_eigenproblem(n, l, Ein, eps, 100, r, rp, u, &
            Z, c, relat, .false., Emin_init, Emax_init, converged, E, P, Q)
        error = abs(E - E_exact)
        if (converged /= 0) call stop_error("Not converged")
        print "(I3, I3, F15.6, F15.6, ES10.2)", n, l, E, E_exact, error
        if (error > 1e-6_dp) call stop_error("Error is higher than 1e-6")
    end do
end do
end program
