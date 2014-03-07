program nonrel_dft_externalpot
use dftatom, only: dp, mesh_exp, mesh_exp_deriv, thomas_fermi_potential, E_nl, &
        get_hydrogen_energies
use dft_data, only: dft_data_t
use dft, only: KS_step
use mixings, only: mixing_anderson
use interpolation, only: spline3, loadtxt
use utils, only: assert, newunit
implicit none

! Mesh parameters:
real(dp), parameter :: r_min = 1e-7_dp, r_max = 50.0_dp, a = 2.7e6_dp
integer, parameter :: N = 5000
real(dp), target :: R(N+1), Rp(N+1)
integer :: Z
!integer, parameter :: n_orb = 11
integer, parameter :: n_orb = 3
integer, target :: no(n_orb), lo(n_orb)
real(dp), target :: fo(n_orb)
real(dp), target :: ks_energies(n_orb)
real(dp), target :: Emin_init(n_orb), Emax_init(n_orb)
real(dp), target :: V_tot(N+1)
real(dp), target :: density(N+1)
real(dp), target :: orbitals(N+1, n_orb)
real(dp), parameter :: reigen_eps = 1e-10_dp
real(dp), parameter :: mixing_eps = 5e-9_dp
real(dp), parameter :: mixing_alpha = 0.5_dp
integer, parameter :: mixing_max_iter = 200, reigen_max_iter = 100
logical :: perturb = .false.
real(dp), dimension(size(R)), target :: V_h, V_xc, e_xc, V_coulomb, tmp
real(dp), dimension(size(R), 0:2), target :: V_l
type(dft_data_t) :: d
integer :: i, u
character, parameter :: l_names(0:3) = (/ "s", "p", "d", "f" /)
real(dp), allocatable :: data(:, :)

Z = 14
! Configuration for Z=50:
!no = [1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5]
!lo = [0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1]
!fo = [2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 2]
! You can modify the occupation here:
!fo(10) = 1.5_dp
!fo(11) = 2.0_dp

no = [1, 2, 3]
lo = [0, 1, 2]
fo = [2, 2, 10]

!no = [1]
!lo = [0]
!fo = [2]

R = mesh_exp(r_min, r_max, a, N)
Rp = mesh_exp_deriv(r_min, r_max, a, N)

!ks_energies = [-2, -1] ! Initial guess for energies
!ks_energies = get_hydrogen_energies(Z, no)
ks_energies(1) = -0.738977_dp
ks_energies(2) = -0.738977_dp
ks_energies(3) = -0.738977_dp
V_tot = thomas_fermi_potential(R, Z) ! Initial guess for the potential

! R = data(1, :)
! V0 = data(2, :)
! V1 = data(3, :)
! V2 = data(4, :)
call loadtxt("sn-pseudo.txt", data)
! Coulomb potential:
!     V_coulomb = -Z/R
! Uncomment this to use the potential loaded from the file:
V_l(:, 0) = spline3(data(1, :), data(2, :), R)
V_l(:, 1) = spline3(data(1, :), data(3, :), R)
V_l(:, 2) = spline3(data(1, :), data(4, :), R)
V_coulomb = -Z*erf(R)/R
forall(i=0:2) V_l(:, i) = V_l(:, i) - V_coulomb
V_tot = V_coulomb

! We allow a few unbounded states
Emax_init = 10
! Use Hydrogen Schroedinger energies for the lower limit:
do i = 1, size(no)
    Emin_init(i) = E_nl(0._dp, no(i), lo(i), Z, 0)
end do
Emin_init(1) = -0.738977_dp
Emin_init(1) = -30
Emin_init(2) = -30
Emin_init(3) = -30
! For robustness, decrease Emin by 10%:
Emin_init = 1.5_dp * Emin_init

d%Z = 0  ! The boundary condition for non-singular potential is Z=0.
d%R => R
d%Rp => Rp
d%rho => density
d%V_h => V_h
! Here V_coulomb means V_loc, the local part of a pseudopotential
d%V_coulomb => V_coulomb
d%V_l => V_l
d%V_xc => V_xc
d%e_xc => e_xc
d%V_tot => V_tot
d%orbitals => orbitals
d%alpha = mixing_alpha
d%pseudopot = .true.
d%xc_type = 2
d%dirac = .false.
d%perturb = perturb
d%reigen_eps = reigen_eps
d%reigen_max_iter = reigen_max_iter
d%no => no
d%lo => lo
d%fo => fo
d%ks_energies => ks_energies
d%Emax_init => Emax_init
d%Emin_init => Emin_init
tmp = mixing_anderson(KS_step, d%V_tot - d%V_coulomb, mixing_max_iter, &
    .true., d, mixing_alpha, mixing_eps)
! Prints the energies:
print *, "Ekin: ", d%Ekin
print *, "Ecoul:", d%Ecoul
print *, "Eenuc:", d%Eenuc
print *, "Exc:  ", d%Exc
print *, "E_tot:", d%Etot
print *, "state      E [a.u.]             E [Ry]      occupancy"
do i = 1, size(ks_energies)
    print "(I1, A, ' ', F18.6, '   ', F18.6, '   ', F6.3)", no(i), &
        l_names(lo(i)), ks_energies(i), ks_energies(i) * 2, fo(i)
end do
print *
print *, "The first 10 values of the radial grid:"
print *, R(:10)
print *
print *, "The first 10 values of the 1st and 2nd orbitals:"
print *, orbitals(:10, 1)
print *, orbitals(:10, 2)
print *
print *, "The first 10 values of the radial charge density:"
print *, density(:10)

print *, "The first 10 values of the Hartree potential (V_h):"
print *, V_h(:10)

! Save the radial grid, density, V_h
open(newunit(u), file="density.txt", status="replace")
write(u, "((es23.16, ' ', es23.16, ' ', es23.16))") (R(i), density(i), V_h(i), &
        i=1, size(R))
close(u)
end program
