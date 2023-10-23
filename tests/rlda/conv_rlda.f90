program conv_rlda
! Calculates the relativistic DFT energies for Z=92 (U)
!
! The purpose of this test is to check that:
! a) SCF works
! b) all SCF U energies can converge to the given accuracy
! c) it converges to the correct energies (XC and other DFT things are correct)
use dftatom, only: atom_rlda
use types, only: dp
use utils, only: str, stop_error, newunit
implicit none

! Atomic number:
integer :: Z
! Mesh parameters:
real(dp), parameter :: r_min = 1e-8_dp, r_max = 50.0_dp, a = 6.2e7_dp
integer :: NN

! 1986 CODATA  -4223.41902095
real(dp), parameter :: c = 137.0359895_dp

integer :: i, n_orb
character, parameter :: l_names(0:3) = (/ "s", "p", "d", "f" /)
real(dp) :: err, E_tot, E_tot_exact
real(dp), parameter :: reigen_eps = 1e-10_dp
real(dp), parameter :: mixing_eps = 5e-9_dp
integer, allocatable, dimension(:) :: no, lo, so
real(dp), allocatable, dimension(:) :: fo, ks_energies
real(dp), pointer, dimension(:) :: ks_energies_exact
real(dp), allocatable :: orbitals(:, :)
real(dp), allocatable :: R(:), Rp(:), V_tot(:), density(:)
real(dp) :: eps
integer :: p

do p = 3, 8
    eps = 10.0_dp**(-p)
    eps = eps * 1.2_dp ! Allow numerical differences across compilers/platforms
    print *, "Test eps:", eps
    do Z = 92, 1, -1
        call get_LDA_energies(Z, ks_energies_exact, E_tot_exact)
        n_orb = size(ks_energies_exact)
        NN = get_N(Z, p)
        allocate(ks_energies(n_orb), no(n_orb), &
            lo(n_orb), fo(n_orb), orbitals(NN+1, n_orb), R(NN+1), V_tot(NN+1), &
            density(NN+1), so(n_orb), Rp(NN+1))
        call atom_rlda(Z, r_min, r_max, a, NN, c, no, lo, so, fo, &
            ks_energies, E_tot, &
            R, Rp, V_tot, density, orbitals, &
            reigen_eps, 100, mixing_eps, 0.5_dp, 200, .true.)

        print *, "Z=", Z
        print *, "N=", NN
        err = abs(E_tot - E_tot_exact)
        print '("E_tot=", F16.8, " E_tot_exact=", F16.8, " error:", ES10.2)', &
                E_tot, E_tot_exact, err
        if (err > eps) call error(err, eps)
        print *, "state    E            E_exact          error     occupancy"
        do i = 1, size(ks_energies)
            err = (ks_energies_exact(i) - ks_energies(i))
            print "(I1, A, ' ', F16.8, F16.8, ES10.2, '   ', F6.3)", no(i), &
                    l_names(lo(i)), ks_energies(i), ks_energies_exact(i), err, fo(i)
            if (err > eps) call error(err, eps)
        end do
        deallocate(ks_energies, ks_energies_exact, no, lo, fo, so, orbitals, &
            R, Rp, V_tot, density)
    end do
end do

contains

subroutine get_LDA_energies(Z, ks_energies, E_tot)
integer, intent(in) :: Z
real(dp), pointer :: ks_energies(:)
real(dp), intent(out) :: E_tot
integer :: i, fZ, n, u
open(newunit(u), file="rel_energies.dat", status="old")
do i = 1, 5
    read(u, *) ! header
end do
read(u, *) fZ
do while (fZ /= Z)
    read(u, *) ! E_tot
    read(u, *) ! n, ks_energies(n)
    read(u, *) fZ
end do
read(u, *) E_tot
read(u, *) n
allocate(ks_energies(n))
backspace(u)
read(u, *) n, ks_energies
close(u)
end subroutine

integer function get_N(Z, p) result(N)
integer, intent(in) :: Z, p
integer :: fZ, u
open(newunit(u), file="data" // str(p) // ".dat", status="old")
read(u, *) fZ, N
do while (fZ /= Z)
    read(u, *) fZ, N
end do
close(u)
end function

subroutine error(err, eps)
real(dp), intent(in) :: err, eps
print "('Test failed: error = ', es10.2, '   > ', es10.2, ' specified.')", &
        err, eps
call stop_error("Aborting...")
end subroutine

end program
