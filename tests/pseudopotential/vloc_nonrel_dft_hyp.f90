program vloc_nonrel_dft_hyp
! The same as vloc_nonrel_dft, but with a hyperbolic mesh. The purpose of this
! file is to show how to use a user defined mesh.
use dftatom, only: dp, thomas_fermi_potential, E_nl
use dft_data, only: dft_data_t
use dft, only: KS_step
use mixings, only: mixing_anderson
use utils, only: assert, stop_error, newunit
implicit none

! Mesh parameters:
real(dp), parameter :: r_min = 1e-7_dp, r_max = 50.0_dp, a = 1.003_dp
integer, parameter :: N = 5000
real(dp), target :: R(N+1), Rp(N+1)
real(dp) :: P(5)
integer :: Z
integer, parameter :: n_orb = 2
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
integer :: i!, u
character, parameter :: l_names(0:3) = (/ "s", "p", "d", "f" /)

Z = 5
no = [1, 2]
lo = [0, 1]
fo = [2, 3]

R = mesh_hyp(r_min, r_max, a, N)
Rp = mesh_hyp_deriv(r_min, r_max, a, N)

ks_energies = [-2, -1] ! Initial guess for energies
V_tot = thomas_fermi_potential(R, Z) ! Initial guess for the potential

!V_coulomb = -Z/R
P = [0.65435_dp, 2.45106_dp, -1.536643785333E-01_dp, 1.153664378533E+00_dp, &
    5.0000_dp]
V_coulomb = -P(5)/R * (P(3)*erf(sqrt(P(1)) * R) + P(4)*erf(sqrt(P(2)) * R))

! We allow a few unbounded states
Emax_init = 10
! Use Hydrogen Schroedinger energies for the lower limit:
do i = 1, size(no)
    Emin_init(i) = E_nl(0._dp, no(i), lo(i), Z, 0)
end do
! For robustness, decrease Emin by 10%:
Emin_init = 1.1_dp * Emin_init

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
print *, orbitals(:10, 2)
print *
print *, "The first 10 values of the radial charge density:"
print *, density(:10)

! Save the radial grid and density
!open(newunit(u), file="density.txt", status="replace")
!write(u, "((es23.16, ' ', es23.16))") (R(i), density(i), i=1, size(R))
!close(u)

call assert(abs(ks_energies(1) - (-1.613029_dp)) < 1e-6_dp)
call assert(abs(ks_energies(2) - (-0.181696_dp)) < 1e-6_dp)
call assert(abs(d%Etot - (-10.682131_dp)) < 1e-6_dp)

contains

function mesh_hyp(r_min, r_max, a, N) result(R)
! Generates hyperbolic mesh of N elements on [r_min, r_max]
real(dp), intent(in) :: r_min
real(dp), intent(in) :: r_max
real(dp), intent(in) :: a
integer, intent(in) :: N
real(dp) :: R(N+1)

integer :: i
if (N < 1) call stop_error("mesh_hyp() requires N >= 1")
do i = 0, N
    R(i+1) = i * (a - 1) / (a*N - i)  * (r_max - r_min) + r_min
end do
end function

function mesh_hyp_deriv(r_min, r_max, a, N) result(Rp)
! Generates the derivative of a hyperbolic mesh of N elements on [r_min, r_max]
real(dp), intent(in) :: r_min
real(dp), intent(in) :: r_max
real(dp), intent(in) :: a
integer, intent(in) :: N
real(dp) :: Rp(N+1)

integer :: i
if (N < 1) call stop_error("mesh_hyp_deriv() requires N >= 1")
do i = 0, N
    Rp(i+1) = a * N * (a - 1) / (a*N - i)**2  * (r_max - r_min)
end do
end function

end program
