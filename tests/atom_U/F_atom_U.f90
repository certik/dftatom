program atom_U
use dftatom, only: dp, stop_error, get_atom_orb, &
    atom_lda
implicit none

! Atomic number:
integer :: Z
! Mesh parameters:
real(dp), parameter :: r_min = 1e-7_dp, r_max = 50.0_dp, a = 2.7e6_dp
integer :: NN


integer :: i, n_orb
character, parameter :: l_names(0:3) = (/ "s", "p", "d", "f" /)
real(dp) :: E_tot
real(dp), parameter :: reigen_eps = 1e-10_dp
real(dp), parameter :: mixing_eps = 5e-9_dp
integer, allocatable, dimension(:) :: no, lo
real(dp), allocatable, dimension(:) :: fo, ks_energies
real(dp), allocatable :: orbitals(:, :)
real(dp), allocatable :: R(:), Rp(:), V_tot(:), density(:)

Z = 92
NN = 5500
n_orb = get_atom_orb(Z)
allocate(ks_energies(n_orb), no(n_orb), &
    lo(n_orb), fo(n_orb), orbitals(NN+1, n_orb), R(NN+1), V_tot(NN+1), &
    density(NN+1), Rp(NN+1))
call atom_lda(Z, r_min, r_max, a, NN, no, lo, fo, ks_energies, E_tot, &
    R, Rp, V_tot, density, orbitals, reigen_eps, 100, mixing_eps, &
    0.5_dp, 200, .true.)

print *, "Z=", Z, "N=", NN
print '("E_tot=", F18.6)', E_tot
print *, "state    E            occupancy"
do i = 1, size(ks_energies)
    print "(I1, A, ' ', F18.6, '   ', F6.3)", no(i), l_names(lo(i)), &
        ks_energies(i), fo(i)
end do
print *, "Print the first 10 values of the 1st and 2nd orbitals:"
print *, orbitals(:10, 1)
print *, orbitals(:10, 2)

end program
