module drivers

! This module contains high level drivers for atomic SCF calculations
!
! The idea is to use these drivers to do most frequent calculations with
! an exponential mesh and to get an idea how things work. They can be used as a
! starting point/template to write a custom solver for a particular problem,
! or to use a different mesh.

use types, only: dp
use constants, only: pi
use ode1d, only: integrate
use states, only: get_atomic_states_nonrel, get_atomic_states_rel
use dft, only: KS_step, rho2V
use mesh, only: mesh_exp, mesh_exp_deriv
use energies, only: get_hydrogen_energies, thomas_fermi_potential, E_nl
use dft_data, only: dft_data_t
use mixings, only: mixing_anderson
implicit none

private
public get_atom_orb, atom_lda, atom_rlda, get_atom_orb_rel

contains

integer function get_atom_orb(Z) result(n_orb)
! Returns the number of orbitals for the given atom
integer, intent(in) :: Z

integer, pointer, dimension(:) :: no, lo
real(dp), pointer, dimension(:) :: fo

call get_atomic_states_nonrel(Z, no, lo, fo)
n_orb = size(no)
deallocate(no, lo, fo)
end function

integer function get_atom_orb_rel(Z) result(n_orb)
! Returns the number of orbitals for the given atom
integer, intent(in) :: Z

integer, pointer, dimension(:) :: no, lo, so
real(dp), pointer, dimension(:) :: fo

call get_atomic_states_rel(Z, no, lo, so, fo)
n_orb = size(no)
deallocate(no, lo, so, fo)
end function

subroutine atom_lda(Z, r_min, r_max, a, N, no, lo, fo, ks_energies, E_tot, &
    R, Rp, V_tot, density, orbitals, reigen_eps, reigen_max_iter, &
    mixing_eps, mixing_alpha, &
    mixing_max_iter, perturb)
integer, intent(in) :: Z
real(dp), intent(in) :: r_min, r_max, a
integer, intent(in) :: N
integer, intent(out), target :: no(:), lo(:)
real(dp), intent(out), target :: fo(:)
real(dp), intent(out), target :: ks_energies(:)
real(dp), intent(out) :: E_tot
real(dp), intent(out), target :: R(:), Rp(:)
real(dp), intent(out), target :: V_tot(:)
real(dp), intent(out), target :: density(:)
real(dp), intent(out), target :: orbitals(:, :)
real(dp), intent(in) :: reigen_eps, mixing_eps, mixing_alpha
integer, intent(in) :: mixing_max_iter, reigen_max_iter
logical, intent(in) :: perturb

integer, pointer :: no_a(:), lo_a(:)
real(dp), pointer :: fo_a(:)
integer :: i
real(dp), dimension(size(ks_energies)), target :: Emin_init, Emax_init

real(dp), dimension(size(R)), target :: V_h, V_xc, e_xc, V_coulomb, tmp
type(dft_data_t) :: d

call get_atomic_states_nonrel(Z, no_a, lo_a, fo_a)
no = no_a
lo = lo_a
fo = fo_a
deallocate(no_a, lo_a, fo_a)

R = mesh_exp(r_min, r_max, a, N)
Rp = mesh_exp_deriv(r_min, r_max, a, N)

!V_tot = -Z / R
ks_energies = get_hydrogen_energies(Z, no)
V_tot = thomas_fermi_potential(R, Z)
!ks_energies = get_tf_energies(Z, no, fo)

V_coulomb = -Z/R

! We allow a few unbounded states
Emax_init = 10
! Use Hydrogen Schroedinger energies for the lower limit:
do i = 1, size(no)
    Emin_init(i) = E_nl(0._dp, no(i), lo(i), Z, 0)
end do
! For robustness, decrease Emin by 10%:
Emin_init = 1.1_dp * Emin_init


d%Z = Z
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
! Start from an initial density instead:
!d%rho = 1 / cosh(d%R)**2
!d%rho = d%rho * d%Z / integrate(d%Rp, 4*pi*d%rho*d%R**2)
!call rho2V(d)
tmp = mixing_anderson(KS_step, d%V_tot - d%V_coulomb, mixing_max_iter, &
    .true., d, mixing_alpha, mixing_eps)
E_tot = d%Etot
! Prints the energies:
!print *, "Ekin: ", d%Ekin
!print *, "Ecoul:", d%Ecoul
!print *, "Eenuc:", d%Eenuc
!print *, "Exc:  ", d%Exc
!print *, "E_tot:", d%Ekin + d%Ecoul + d%Eenuc + d%Exc
end subroutine

subroutine atom_rlda(Z, r_min, r_max, a, N, c, no, lo, so, fo, ks_energies, &
    E_tot, &
    R, Rp, V_tot, density, orbitals, reigen_eps, reigen_max_iter, &
    mixing_eps, mixing_alpha, &
    mixing_max_iter, perturb)
integer, intent(in) :: Z
real(dp), intent(in) :: r_min, r_max, a, c
integer, intent(in) :: N
integer, intent(out), target :: no(:), lo(:), so(:)
real(dp), intent(out), target :: fo(:)
real(dp), intent(out), target :: ks_energies(:)
real(dp), intent(out) :: E_tot
real(dp), intent(out), target :: R(:), Rp(:)
real(dp), intent(out), target :: V_tot(:)
real(dp), intent(out), target :: density(:)
real(dp), intent(out), target :: orbitals(:, :)
real(dp), intent(in) :: reigen_eps, mixing_eps, mixing_alpha
integer, intent(in) :: mixing_max_iter, reigen_max_iter
logical, intent(in) :: perturb

integer, pointer :: no_a(:), lo_a(:), so_a(:)
real(dp), pointer :: fo_a(:)
real(dp), dimension(size(ks_energies)), target :: Emin_init, Emax_init
integer :: i, relat

real(dp), dimension(size(R)), target :: V_h, V_xc, e_xc, V_coulomb, tmp
type(dft_data_t) :: d

call get_atomic_states_rel(Z, no_a, lo_a, so_a, fo_a)
no = no_a
lo = lo_a
so = so_a
fo = fo_a
deallocate(no_a, lo_a, so_a, fo_a)

R = mesh_exp(r_min, r_max, a, N)
Rp = mesh_exp_deriv(r_min, r_max, a, N)

!V_tot = -Z / R
ks_energies = get_hydrogen_energies(Z, no)
V_tot = thomas_fermi_potential(R, Z)
!ks_energies = get_tf_energies(Z, no, fo)

V_coulomb = -Z/R

! We allow a few unbounded states
Emax_init = 10
! Use Hydrogen Dirac energies for the lower limit:
do i = 1, size(no)
    if (so(i) == 1) then
        relat = 2
    else
        relat = 3
    end if
    Emin_init(i) = E_nl(c, no(i), lo(i), Z, relat)
end do
! For robustness, decrease Emin by 10%:
Emin_init = 1.1_dp * Emin_init


d%Z = Z
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
d%dirac = .true.
d%perturb = perturb
d%reigen_eps = reigen_eps
d%reigen_max_iter = reigen_max_iter
d%no => no
d%lo => lo
d%fo => fo
d%ks_energies => ks_energies
d%Emax_init => Emax_init
d%Emin_init => Emin_init

! Relativistic parameters:
d%c = c
d%so => so

tmp = mixing_anderson(KS_step, d%V_tot - d%V_coulomb, mixing_max_iter, &
    .true., d, mixing_alpha, mixing_eps)
E_tot = d%Etot
end subroutine

end module
