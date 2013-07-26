program uranium
! Calculates Thomas-Fermi nonrelativistic energies for Z=92 (U)
!
use dftatom
implicit none

! Atomic number:
integer, parameter :: Z = 92
! Mesh parameters:
real(dp), parameter :: r_min = 1e-8_dp, r_max = 50.0_dp, a = 1e6_dp
integer, parameter :: NN = 3000

real(dp), parameter :: c = 137.035999037_dp, eps = 1e-6_dp
integer :: converged, n, l, Z_eff
real(dp) :: R(NN+1), V(size(r)), E, P(size(r)), Q(size(r)), E_initial
real(dp) :: Rp(NN+1)

integer, parameter :: relat = 0

R = mesh_exp(r_min, r_max, a, NN)
Rp = mesh_exp_deriv(r_min, r_max, a, NN)
V = thomas_fermi_potential(R, Z)

print *, "Thomas Fermi nonrelativistic energies for Z=92 (U)"
print *, "Mesh parameters (r_min, r_max, a, N):"
print "(ES10.2, F10.2, ES10.2, I10)", r_min, r_max, a, NN
print *
print "(A3, A3, A15, A15, A10)", "n", "l", "E", "E_initial"
print *
do n = 1, 7
    do l = 0, n-1
        Z_eff = Z - 2*((n-1)**2 + 2*l+1)
        if (Z_eff <= 0) Z_eff = 1
        E_initial = -Z_eff**2 * 1.0_dp/(2*n**2)
        call solve_radial_eigenproblem(n, l, E_initial, eps, 100, R, Rp, &
            V, Z, c, relat, .true., -10000._dp, 0._dp, converged, E, P, Q)
        print "(I3, I3, F15.6, F15.6)", n, l, E, E_initial
    end do
end do
end program
