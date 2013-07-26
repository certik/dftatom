program vloc_nonrel
use dftatom
implicit none

! For non-singular potential, we need to set Z=0:
integer, parameter :: Z = 0
! Mesh parameters:
real(dp), parameter :: r_min = 1e-7_dp, r_max = 50.0_dp, a = 2.7e6_dp
integer, parameter :: NN = 5000

real(dp), parameter :: c = 137.035999037_dp, eps = 1e-12_dp
integer :: n, l, relat, converged
real(dp) :: r(NN+1), u(size(r)), Ein, E, P(size(r)), Q(size(r))
real(dp) :: Rp(NN+1)
real(dp) :: Emin_init, Emax_init
real(dp) :: D(5)

R = mesh_exp(r_min, r_max, a, NN)
Rp = mesh_exp_deriv(r_min, r_max, a, NN)
D = [0.65435_dp, 2.45106_dp, -1.536643785333E-01_dp, 1.153664378533E+00_dp, &
    5.0000_dp]
u = -D(5)/R * (D(3)*erf(sqrt(D(1)) * R) + D(4)*erf(sqrt(D(2)) * R))

print *, "Local part of a pseudopotential for N (5 electrons)"
print *, "Mesh parameters (r_min, r_max, a, N):"
print "(ES10.2, F10.2, ES10.2, I10)", r_min, r_max, a, NN
print *
print "(A3, A3, A15)", "n", "l", "E"
print *
do n = 1, 7
    do l = 0, n-1
        relat = 0
        Ein = 10
        Emin_init = -100
        Emax_init = 100
        call solve_radial_eigenproblem(n, l, Ein, eps, 100, r, rp, u, &
            Z, c, relat, .false., Emin_init, Emax_init, converged, E, P, Q)
        if (converged /= 0) call stop_error("Not converged")
        print "(I3, I3, F15.6)", n, l, E
    end do
end do
end program
