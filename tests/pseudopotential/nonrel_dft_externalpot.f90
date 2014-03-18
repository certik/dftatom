program nonrel_dft_externalpot
use dftatom, only: dp, mesh_exp, mesh_exp_deriv, thomas_fermi_potential, E_nl, &
        get_hydrogen_energies, atom_lda_pseudo
use interpolation, only: spline3, loadtxt
use utils, only: assert, newunit
use ode1d, only: integrate
use constants, only: pi
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
logical :: perturb = .true.
real(dp), dimension(size(R)), target :: V_loc
real(dp), dimension(size(R), 0:2), target :: V_l
integer :: i
character, parameter :: l_names(0:3) = (/ "s", "p", "d", "f" /)
real(dp), allocatable :: data(:, :)
real(dp) :: Ekin, Eee, Een, Exc, Etot, Enl
integer :: xc_type

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
! Uncomment this to use the potential loaded from the file:
V_l(:, 0) = spline3(data(1, :), data(2, :), R)
V_l(:, 1) = spline3(data(1, :), data(3, :), R)
V_l(:, 2) = spline3(data(1, :), data(4, :), R)
! Use the s-channel for the V_loc potential:
V_loc = V_l(:, 0)
forall(i=0:2) V_l(:, i) = V_l(:, i) - V_loc
V_tot = V_loc ! Initial guess for the total potential

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

xc_type = 2
call atom_lda_pseudo(no, lo, fo, Emin_init, Emax_init, ks_energies, &
    R, Rp, V_loc, V_l, V_tot, density, orbitals, Ekin, Eee, Een, Exc, Etot, &
    Enl, &
    reigen_eps, reigen_max_iter, mixing_eps, mixing_alpha, mixing_max_iter, &
    perturb, xc_type)
! Prints the energies:
print *, "Ekin: ", Ekin
print *, "Ecoul:", Eee
print *, "Eenuc:", Een
print *, "Exc:  ", Exc
print *, "Etot: ", Etot
print *, "Enl:  ", Enl, "a.u. =", 2*Enl, "Ry"
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

end program
