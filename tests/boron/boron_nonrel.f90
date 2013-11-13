program boron
! Calculates Hydrogen like nonrelativistic energies for Z=5 (B)
!
! The purpose of this test is to check that:
! a) all U energies can converge to 1e-6
use dftatom
implicit none

! Atomic number:
integer, parameter :: Z = 5
! Mesh parameters:
real(dp), parameter :: r_min = 1e-8_dp, r_max = 50.0_dp, a = 1e6_dp
integer, parameter :: NN = 3000

real(dp), parameter :: c = 137.035999037_dp, eps = 1e-6_dp
integer :: n, l, relat, converged, i
real(dp) :: r(NN+1), u(size(r)), Ein, E, E_exact, error, P(size(r)), Q(size(r))
real(dp) :: Rp(NN+1)
integer :: n_orb
integer, pointer :: no(:), lo(:)
real(dp), pointer :: fo(:)
logical, parameter :: perturb = .true.

call get_atomic_states_nonrel(Z, no, lo, fo)
n_orb = size(no)

r = mesh_exp(r_min, r_max, a, NN)
Rp = mesh_exp_deriv(r_min, r_max, a, NN)
u(:) = -Z/r

print *, "Hydrogen like nonrelativistic energies for Z=5 (B)"
print *, "Mesh parameters (r_min, r_max, a, N):"
print "(ES10.2, F10.2, ES10.2, I10)", r_min, r_max, a, NN
print *
print "(A3, A3, A15, A15, A10)", "n", "l", "E", "E_exact", "Error"
print *
do i = 1, n_orb
    n = no(i)
    l = lo(i)
    relat = 0
    E_exact = E_nl(c, n, l, Z, relat)
    Ein = -100
    call solve_radial_eigenproblem(n, l, Ein, eps, 100, r, Rp, u, Z, c, &
        relat, perturb, -10000._dp, 0._dp, converged, E, P, Q)
    error = abs(E - E_exact)
    if (converged /= 0) call stop_error("Not converged")
    print "(I3, I3, F15.6, F15.6, ES10.2)", n, l, E, E_exact, error
    if (error > eps) call stop_error("Error is higher than 1e-6")
end do
deallocate(no, lo, fo)
end program
