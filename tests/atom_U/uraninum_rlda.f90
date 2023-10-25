program uranium_rlda
! Calculates the relativistic DFT energies for Z=92 (U) to 1e-6 a.u. accuracy.
!
! This program can be used to benchmark the performance of dftatom.
use dftatom, only: atom_rlda
use types, only: dp
use utils, only: stop_error
implicit none

! Atomic number:
integer :: Z = 92 ! Uraninum
! Mesh parameters:
real(dp), parameter :: r_min = 1e-8_dp, r_max = 50.0_dp, a = 6.2e7_dp
integer, parameter :: NN = 5269

! 1986 CODATA  -4223.41902095
real(dp), parameter :: c = 137.0359895_dp

integer :: i, n_orb
character, parameter :: l_names(0:3) = (/ "s", "p", "d", "f" /)
real(dp) :: err, E_tot
real(dp), parameter :: E_tot_exact = -28001.1323254868_dp
real(dp), parameter :: reigen_eps = 1e-10_dp
real(dp), parameter :: mixing_eps = 5e-9_dp
integer, allocatable, dimension(:) :: no, lo, so
real(dp), allocatable, dimension(:) :: fo, ks_energies
real(dp), parameter :: ks_energies_exact(29) = [ &
                -4223.4190204552_dp, &
                -789.4897823303_dp, &
                -761.3744759730_dp, &
                -622.8480945649_dp, &
                -199.4298056450_dp, &
                -186.6637131249_dp, &
                -154.7010266741_dp, &
                -134.5411802896_dp, &
                -128.0166573820_dp, &
                -50.7889480646_dp, &
                -45.0371712884_dp, &
                -36.6886104859_dp, &
                -27.5293062430_dp, &
                -25.9854289064_dp, &
                -13.8895142333_dp, &
                -13.4854696912_dp, &
                -11.2955870987_dp, &
                -9.0579642498_dp, &
                -7.0692956350_dp, &
                -3.7974162278_dp, &
                -3.5012171832_dp, &
                -0.1467883850_dp, &
                -0.1160471651_dp, &
                -1.7480399541_dp, &
                -1.1011189998_dp, &
                -0.7757841787_dp, &
                -0.1030408153_dp, &
                -0.0848020246_dp, &
                -0.1609472826_dp ]
real(dp), allocatable :: orbitals(:, :)
real(dp), allocatable :: R(:), Rp(:), V_tot(:), density(:)
real(dp) :: eps
integer :: p

p = 6
eps = 10.0_dp**(-p)
eps = eps * 1.2_dp ! Allow numerical differences across compilers/platforms
print *, "Test eps:", eps
Z = 92
n_orb = size(ks_energies_exact)
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
if (err > eps) call stop_error("err > eps")
print *, "state    E            E_exact          error     occupancy"
do i = 1, size(ks_energies)
    err = (ks_energies_exact(i) - ks_energies(i))
    print "(I1, A, ' ', F16.8, F16.8, ES10.2, '   ', F6.3)", no(i), &
            l_names(lo(i)), ks_energies(i), ks_energies_exact(i), err, fo(i)
    if (err > eps) call stop_error("err > eps")
end do

end program
