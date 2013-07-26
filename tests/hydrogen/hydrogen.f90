program Hydrogen
use dftatom, only: mesh_exp, mesh_exp_deriv, dp, &
    solve_radial_eigenproblem, stop_error
implicit none

! Atomic number:
integer, parameter :: Z = 92
! Mesh parameters:
real(dp), parameter :: r_min = 1e-8_dp, r_max = 50.0_dp, a = 1e6_dp
integer, parameter :: NN = 3000

real(dp), parameter :: c = 137.035999037_dp, eps = 1e-9_dp
integer :: n, l, relat, converged
real(dp) :: R(NN+1), u(size(r)), Ein, E, E_exact, error
real(dp) :: P(size(R)), Q(size(R)), Rp(NN+1)
logical, parameter :: perturb = .true.


R = mesh_exp(r_min, r_max, a, NN)
Rp = mesh_exp_deriv(r_min, r_max, a, NN)
u(:) = -Z/r

print *, "Hydrogen like energies for Z=92 (U)"
print *, "Mesh parameters (r_min, r_max, a, N):"
print "(ES10.2, F10.2, ES10.2, I10)", r_min, r_max, a, NN
print *
print "(A3, A3, A15, A15, A10)", "n", "l", "E", "E_exact", "Error"
print *
do n = 1, 7
    do l = 0, n-1
        E_exact = - Z**2 / (2.0_dp * n**2)
        Ein = -1000
        relat = 0
        call solve_radial_eigenproblem(n, l, Ein, eps, 100, R, Rp, u, Z, c, &
            relat, perturb, -5000._dp, 0._dp, converged, E, P, Q)
        error = abs(E -E_exact)
        if (converged /= 0) call stop_error("Not converged")
        print "(I3, I3, F15.6, F15.6, ES10.2)", n, l, E, E_exact, error
    end do
end do
end program
