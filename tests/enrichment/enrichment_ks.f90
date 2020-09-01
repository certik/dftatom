program enrichment_ks
use dftatom, only: dp, mesh_exp, mesh_exp_deriv, thomas_fermi_potential, E_nl
use dft_data, only: dft_data_t
use dft, only: KS_step
use mixings, only: mixing_anderson
use utils, only: assert, newunit
implicit none

! Mesh parameters:
real(dp), parameter :: r_min = 1e-3_dp, r_max = 30.0_dp, a = 1e6_dp
integer, parameter :: N = 150
real(dp), target :: R(N+1), Rp(N+1)
integer :: Z
integer, parameter :: n_orb = 1
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
integer, parameter :: mixing_max_iter = 200, reigen_max_iter = 40
logical :: perturb = .true.
real(dp), dimension(size(R)), target :: V_h, V_xc, e_xc, V_coulomb, tmp
type(dft_data_t) :: d
integer :: i, u
character, parameter :: l_names(0:3) = (/ "s", "p", "d", "f" /)
real(dp) :: Zion, rloc, C1, C2 !, C3, C4

Z = 1
no = [1]
lo = [0]
fo = [1]

R = mesh_exp(r_min, r_max, a, N)
Rp = mesh_exp_deriv(r_min, r_max, a, N)

ks_energies = [-2] ! Initial guess for energies
V_tot = thomas_fermi_potential(R, Z) ! Initial guess for the potential

rloc = 2.0_dp
C1 = -4.180237_dp
C2 =  0.725075_dp
!C3 = 0
!C4 = 0
Zion = 1
!V_coulomb = -Z/R
V_coulomb = - Zion/r * erf(r/(sqrt(2._dp)*rloc)) &
             + exp(-1._dp/2*(r/rloc)**2) * (C1 + C2*(r/rloc)**2)


! We allow a few unbounded states
Emax_init = 10
Emin_init = -10

d%Z = 0  ! The boundary condition for non-singular potential is Z=0.
d%R => R
d%Rp => Rp
d%rho => density
d%V_h => V_h
d%V_coulomb => V_coulomb
d%V_xc => V_xc
d%e_xc => e_xc
d%V_tot => V_tot
d%orbitals => orbitals
d%alpha = mixing_alpha
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
print *, "state    E            occupancy"
do i = 1, size(ks_energies)
    print "(I1, A, ' ', F18.6, '   ', F6.3)", no(i), l_names(lo(i)), &
        ks_energies(i), fo(i)
end do
print *
print *, "The first 10 values of the radial grid:"
print *, R(:10)
print *
print *, "The first 10 values of the 1st and 2nd orbitals:"
print *, orbitals(:10, 1)
!print *, orbitals(:10, 2)
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
