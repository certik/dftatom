program conv_lda
! Calculates the DFT nonrelativistic energies for Z=92 (U)
!
! The purpose of this test is to check that:
! a) SCF works
! b) all SCF U energies can converge to the given accuracy
! c) it converges to the correct energies (XC and other DFT things are correct)
use dftatom, only: dp, stop_error, get_atom_orb, &
    atom_lda
use utils, only: str, newunit
implicit none

! Atomic number:
integer :: Z
! Mesh parameters:
real(dp), parameter :: r_min = 1e-7_dp, r_max = 50.0_dp, a = 2.7e6_dp
integer :: NN


integer :: i, n_orb
character, parameter :: l_names(0:3) = (/ "s", "p", "d", "f" /)
real(dp) :: err, E_tot, E_tot_exact
real(dp), parameter :: reigen_eps = 1e-10_dp
real(dp), parameter :: mixing_eps = 5e-9_dp
integer, allocatable, dimension(:) :: no, lo
real(dp), allocatable, dimension(:) :: fo, ks_energies
real(dp), pointer, dimension(:) :: ks_energies_exact
real(dp), allocatable :: orbitals(:, :)
real(dp) :: eps
real(dp), allocatable :: R(:), Rp(:), V_tot(:), density(:)
integer :: p

do p = 3, 8
    eps = 10.0_dp**(-p)
    eps = eps * 1.2_dp ! Allow numerical differences across compilers/platforms
    print *, "Test eps:", eps
    do Z = 92, 1, -1
        n_orb = get_atom_orb(Z)
        NN = get_N(Z, p)
        allocate(ks_energies(n_orb), no(n_orb), &
            lo(n_orb), fo(n_orb), orbitals(NN+1, n_orb), R(NN+1), V_tot(NN+1), &
            density(NN+1), Rp(NN+1))
        call atom_lda(Z, r_min, r_max, a, NN, no, lo, fo, ks_energies, E_tot, &
            R, Rp, V_tot, density, orbitals, reigen_eps, 100, mixing_eps, &
            0.5_dp, 200, .true.)
        call get_LDA_energies(Z, ks_energies_exact, E_tot_exact)

        print *, "Z=", Z
        print *, "N=", NN
        err = abs(E_tot - E_tot_exact)
        print '("E_tot=", F18.10, " E_tot_exact=", F18.10, " error:", ES10.2)', &
                E_tot, E_tot_exact, err
        if (err > eps) call error(err, eps)
        print *, "state    E            E_exact          error     occupancy"
        do i = 1, size(ks_energies)
            err = (ks_energies_exact(i) - ks_energies(i))
            print "(I1, A, ' ', F18.10, F18.10, ES10.2, '   ', F6.3)", no(i), &
                    l_names(lo(i)), ks_energies(i), ks_energies_exact(i), err, fo(i)
            if (err > eps) call error(err, eps)
        end do
        deallocate(ks_energies, ks_energies_exact, no, lo, fo, orbitals, R, V_tot, &
            density, Rp)
    end do
end do

contains

subroutine get_LDA_energies(Z, ks_energies, E_tot)
integer, intent(in) :: Z
real(dp), pointer :: ks_energies(:)
real(dp), intent(out) :: E_tot
integer :: i, fZ, n, u
open(newunit(u), file="nonrel_energies.dat", status="old")
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
