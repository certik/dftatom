module c_dftatom

! C interface to dftatom. When this file is modified, run the
! 'utils/generate` script to generate the c_dftatom.h header file
! (and the Cython include file).

use iso_c_binding, only: c_int, c_double, c_bool
use dftatom, only: get_Vxc, dp, solve_radial_eigenproblem, integrate, &
        rpoisson_outward_pc, get_atom_orb, atom_lda, mesh_exp, &
        mesh_exp_deriv, atom_rlda, get_atom_orb_rel, &
        integrate_rproblem_outward, get_Vh
implicit none

contains

subroutine dftatom_get_vxc(n, R, rho, relat, c, V) bind(c)
integer(c_int), intent(in) :: n
real(c_double), intent(in) :: R(n), rho(n), c
logical(c_bool), intent(in) :: relat
real(c_double), intent(out) :: V(n)
! We need to retype this, so that it compiles:
logical :: relat2
real(dp) :: exc(n)
relat2 = relat
call get_Vxc(R, rho, relat2, c, exc, V)
end subroutine

subroutine dftatom_integrate(n, x, f, s) bind(c)
integer(c_int), intent(in) :: n
real(c_double), intent(in), dimension(n) :: x
real(c_double), intent(in), dimension(n) :: f
real(c_double), intent(out) :: s
s = integrate(x, f)
end subroutine

subroutine dftatom_integrate_radial_poisson(n, R, Rp, rho, V) bind(c)
integer(c_int), intent(in) :: n
real(c_double), intent(in) :: R(n), Rp(n), rho(n)
real(c_double), intent(out) :: V(n)
V = rpoisson_outward_pc(R, Rp, rho)
end subroutine

subroutine dftatom_solve_radial_eigenproblem(n_array, n, l, Ein, eps, &
    max_iter, R, Rp, V, Z, c, relat, perturb, Emin_init, Emax_init, &
    converged, E, P, Q) bind(c)
integer(c_int), intent(in) :: n_array, n, l, relat, Z, max_iter
real(c_double), intent(in) :: R(n_array), Rp(n_array), V(n_array), eps, Ein, c
logical(c_bool), intent(in) :: perturb
real(c_double), intent(in) :: Emin_init, Emax_init
integer(c_int), intent(out) :: converged
real(c_double), intent(out) :: P(n_array), Q(n_array), E

logical :: perturb2
perturb2 = perturb

call solve_radial_eigenproblem(n, l, Ein, eps, max_iter, &
    R, Rp, V, Z, c, &
    relat, perturb2, Emin_init, Emax_init, converged, E, P, Q)
end subroutine

subroutine dftatom_integrate_rproblem_outward(N, l, E, R, Rp, V, Z, c, relat, &
    P, Q, imax) bind(c)
integer(c_int), intent(in) :: N, l, relat, Z
real(c_double), intent(in) :: E, c, R(N), Rp(N), V(N)
real(c_double), intent(out) :: P(N), Q(N)
integer(c_int), intent(out) :: imax
call integrate_rproblem_outward(l, E, R, Rp, V, Z, c, relat, P, Q, imax)
end subroutine

subroutine dftatom_get_atom_orb(Z, n_orb) bind(c)
integer(c_int), intent(in) :: Z
integer(c_int), intent(out) :: n_orb
n_orb = get_atom_orb(Z)
end subroutine

subroutine dftatom_get_atom_orb_rel(Z, n_orb) bind(c)
integer(c_int), intent(in) :: Z
integer(c_int), intent(out) :: n_orb
n_orb = get_atom_orb_rel(Z)
end subroutine

subroutine dftatom_atom_lda(Z, r_min, r_max, a, N, n_orb, &
    no, lo, fo, ks_energies, E_tot, R, Rp, V_tot, density, orbitals, &
    reigen_eps, reigen_max_iter, mixing_eps, mixing_alpha, &
    mixing_max_iter, perturb) bind(c)
integer(c_int), intent(in) :: Z
real(c_double), intent(in) :: r_min, r_max, a
integer(c_int), intent(in) :: N, n_orb
integer(c_int), intent(out) :: no(n_orb), lo(n_orb)
real(c_double), intent(out) :: fo(n_orb)
real(c_double), intent(out) :: ks_energies(n_orb)
real(c_double), intent(out) :: E_tot
real(c_double), intent(out) :: R(N+1), Rp(N+1)
real(c_double), intent(out) :: V_tot(N+1)
real(c_double), intent(out) :: density(N+1)
real(c_double), intent(out) :: orbitals(N+1, n_orb)
real(c_double), intent(in) :: reigen_eps, mixing_eps, mixing_alpha
integer(c_int), intent(in) :: mixing_max_iter, reigen_max_iter
logical(c_bool), intent(in) :: perturb

logical :: perturb2
perturb2 = perturb

call atom_lda(Z, r_min, r_max, a, N, no, lo, fo, ks_energies, E_tot, R, Rp, &
    V_tot, density, orbitals, reigen_eps, reigen_max_iter, mixing_eps, &
    mixing_alpha, mixing_max_iter, perturb2)
end subroutine

subroutine dftatom_atom_lda2(n, n_orb, Z, R, Rp, no, lo, fo, V_tot, &
    ks_energies, &
    E_tot, &
    density, orbitals, reigen_eps, reigen_max_iter, &
    mixing_eps, mixing_alpha, &
    mixing_max_iter, perturb) bind(c)
use dft_data, only: dft_data_t
use dft, only: KS_step
use mixings, only: mixing_anderson
integer(c_int), intent(in), target :: n, n_orb, Z
integer(c_int), intent(in), target :: no(n_orb), lo(n_orb)
real(c_double), intent(in), target :: fo(n_orb)
real(c_double), intent(inout), target :: ks_energies(n_orb)
real(c_double), intent(out) :: E_tot
real(c_double), intent(in), target :: R(n), Rp(n)
real(c_double), intent(inout), target :: V_tot(n)
real(c_double), intent(out), target :: density(n)
real(c_double), intent(out), target :: orbitals(n, n_orb)
real(c_double), intent(in) :: reigen_eps, mixing_eps, mixing_alpha
integer(c_int), intent(in) :: mixing_max_iter, reigen_max_iter
logical(c_bool), intent(in) :: perturb

real(dp), dimension(size(R)), target :: V_h, V_xc, e_xc, V_coulomb, tmp
type(dft_data_t) :: d

V_coulomb = -Z/R


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
tmp = mixing_anderson(KS_step, d%V_tot - d%V_coulomb, mixing_max_iter, &
    .true., d, mixing_alpha, mixing_eps)
E_tot = d%Etot
end subroutine

subroutine dftatom_atom_rlda(Z, r_min, r_max, a, N, c, n_orb, &
    no, lo, so, fo, ks_energies, E_tot, R, Rp, V_tot, density, orbitals, &
    reigen_eps, reigen_max_iter, mixing_eps, mixing_alpha, &
    mixing_max_iter, perturb) bind(c)
integer(c_int), intent(in) :: Z
real(c_double), intent(in) :: r_min, r_max, a, c
integer(c_int), intent(in) :: N, n_orb
integer(c_int), intent(out) :: no(n_orb), lo(n_orb), so(n_orb)
real(c_double), intent(out) :: fo(n_orb)
real(c_double), intent(out) :: ks_energies(n_orb)
real(c_double), intent(out) :: E_tot
real(c_double), intent(out) :: R(N+1), Rp(N+1)
real(c_double), intent(out) :: V_tot(N+1)
real(c_double), intent(out) :: density(N+1)
real(c_double), intent(out) :: orbitals(N+1, n_orb)
real(c_double), intent(in) :: reigen_eps, mixing_eps, mixing_alpha
integer(c_int), intent(in) :: mixing_max_iter, reigen_max_iter
logical(c_bool), intent(in) :: perturb

logical :: perturb2
perturb2 = perturb
call atom_rlda(Z, r_min, r_max, a, N, c, no, lo, so, fo, ks_energies, E_tot, &
    R, Rp, V_tot, &
    density, orbitals, reigen_eps, reigen_max_iter, mixing_eps, &
    mixing_alpha, mixing_max_iter, perturb2)
end subroutine

subroutine dftatom_mesh_exp(r_min, r_max, a, N, R) bind(c)
real(c_double), intent(in) :: r_min
real(c_double), intent(in) :: r_max
real(c_double), intent(in) :: a
integer(c_int), intent(in) :: N
real(c_double), intent(out) :: R(N+1)
R = mesh_exp(r_min, r_max, a, N)
end subroutine

subroutine dftatom_mesh_exp_deriv(r_min, r_max, a, N, R) bind(c)
real(c_double), intent(in) :: r_min
real(c_double), intent(in) :: r_max
real(c_double), intent(in) :: a
integer(c_int), intent(in) :: N
real(c_double), intent(out) :: R(N+1)
R = mesh_exp_deriv(r_min, r_max, a, N)
end subroutine

subroutine dftatom_get_vh(n, R, Rp, rho, V) bind(c)
integer(c_int), intent(in) :: n
real(c_double), intent(in) :: R(n), Rp(n), rho(n)
real(c_double), intent(out) :: V(n)
V = get_Vh(R, Rp, rho)
end subroutine

end module
