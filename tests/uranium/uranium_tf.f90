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
integer :: converged, n, l, i
real(dp) :: R(NN+1), V(size(r)), E, P(size(r)), Q(size(r)), Rp(NN+1)

integer, parameter :: relat = 0
integer :: n_orb
integer, dimension(:), pointer :: no, lo
real(dp), dimension(:), pointer :: fo
real(dp), dimension(:), allocatable :: E_initial

call get_atomic_states_nonrel(Z, no, lo, fo)
n_orb = size(no)
allocate(E_initial(n_orb))

E_initial = get_tf_energies(Z, no, fo)

R = mesh_exp(r_min, r_max, a, NN)
Rp = mesh_exp_deriv(r_min, r_max, a, NN)
V = thomas_fermi_potential(R, Z)

print *, "Thomas Fermi nonrelativistic energies for Z=92 (U)"
print *, "Mesh parameters (r_min, r_max, a, N):"
print "(ES10.2, F10.2, ES10.2, I10)", r_min, r_max, a, NN
print *
print "(A3, A3, A15, A15, A10)", "n", "l", "E", "E_exact", "Error"
print *
do i = 1, n_orb
    n = no(i)
    l = lo(i)
    call solve_radial_eigenproblem(n, l, E_initial(i), eps, 100, &
        R, Rp, V, Z, c, relat, .true., -10000._dp, 0._dp, converged, E, P, Q)
    print "(I3, I3, F15.6, F15.6, ES10.2)", n, l, E, E_initial(i), 0.0_dp
end do
deallocate(no, lo, fo)
end program
