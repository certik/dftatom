from libcpp cimport bool

from numpy cimport import_array, ndarray
from numpy import empty

cimport c_dftatom

class ConvergeError(Exception):
    pass

def getvxc(ndarray[double, mode="c"] R not None,
        ndarray[double, mode="c"] rho not None,
        bool relat, double c):
    cdef int n = len(R)
    assert len(rho) == n
    cdef ndarray[double, mode="c"] V = empty(n)
    c_dftatom.dftatom_get_vxc(&n, &R[0], &rho[0], &relat, &c, &V[0])
    return V

def integrate(ndarray[double, mode="c"] x not None,
        ndarray[double, mode="c"] f not None):
    cdef int n = len(x)
    assert len(f) == n
    cdef double r
    c_dftatom.dftatom_integrate(&n, &x[0], &f[0], &r)
    return r

def integrate_radial_poisson(ndarray[double, mode="c"] density not None,
            ndarray[double, mode="c"] R not None,
            ndarray[double, mode="c"] Rp not None):
    cdef int n = len(density)
    assert len(R) == n
    cdef ndarray[double, mode="c"] V = empty(n)
    c_dftatom.dftatom_integrate_radial_poisson(&n, &R[0], &Rp[0], &density[0],
            &V[0])
    return V

def solve_radial_eigenproblem(double c, int n, int l, double Ein,
        double eps, int max_iter,
        ndarray[double, mode="c"] u not None,
        ndarray[double, mode="c"] R not None,
        ndarray[double, mode="c"] Rp not None,
        int Z, int relat, bool perturb, double Emin_init, double Emax_init):
    cdef int N = len(u)
    assert len(R) == N
    assert len(Rp) == N
    cdef double E
    cdef int converged
    cdef ndarray[double, mode="c"] P = empty(N)
    cdef ndarray[double, mode="c"] Q = empty(N)
    c_dftatom.dftatom_solve_radial_eigenproblem(&N, &n, &l, &Ein, &eps,
            &max_iter,
            &R[0], &Rp[0], &u[0], &Z, &c, &relat, &perturb,
            &Emin_init, &Emax_init,
            &converged, &E, &P[0], &Q[0])
    if converged == 0:
        return E, P, Q
    else:
        raise ConvergeError("Radial eigenproblem didn't converge (%d)" % \
                converged)

def integrate_radial_problem_outward(int l, double E,
        ndarray[double, mode="c"] R not None,
        ndarray[double, mode="c"] Rp not None,
        ndarray[double, mode="c"] V not None,
        int Z, double c, int relat):
    cdef int N = len(V)
    assert len(R) == N
    assert len(Rp) == N
    cdef int imax
    cdef ndarray[double, mode="c"] P = empty(N)
    cdef ndarray[double, mode="c"] Q = empty(N)
    c_dftatom.dftatom_integrate_rproblem_outward(&N, &l, &E, &R[0], &Rp[0],
            &V[0], &Z, &c, &relat, &P[0], &Q[0], &imax)
    return P[:imax], Q[:imax]

def atom_lda(int Z, double r_min, double r_max, double a, int N,
        double reigen_eps, int reigen_max_iter, double mixing_eps,
        double mixing_alpha,
        int mixing_max_iter, bool perturb):
    cdef int n_orb
    cdef double E_tot
    c_dftatom.dftatom_get_atom_orb(&Z, &n_orb)
    cdef ndarray[int, mode="c"] no = empty(n_orb, dtype="int32")
    cdef ndarray[int, mode="c"] lo = empty(n_orb, dtype="int32")
    cdef ndarray[double, mode="c"] fo = empty(n_orb)
    cdef ndarray[double, mode="c"] ks_energies = empty(n_orb)
    cdef ndarray[double, mode="c"] R = empty(N+1)
    cdef ndarray[double, mode="c"] Rp = empty(N+1)
    cdef ndarray[double, mode="c"] V_tot = empty(N+1)
    cdef ndarray[double, mode="c"] density = empty(N+1)
    cdef ndarray[double, ndim=2, mode="fortran"] orbitals = empty((N+1, n_orb),
            order="F")
    c_dftatom.dftatom_atom_lda(&Z, &r_min, &r_max, &a, &N, &n_orb,
            &no[0], &lo[0], &fo[0], &ks_energies[0], &E_tot,
            &R[0], &Rp[0], &V_tot[0], &density[0], &orbitals[0, 0],
            &reigen_eps, &reigen_max_iter, &mixing_eps, &mixing_alpha,
            &mixing_max_iter, &perturb)
    return E_tot, ks_energies, no, lo, fo, R, Rp, V_tot, density, orbitals

def atom_lda2(int Z,
        ndarray[double, mode="c"] R not None,
        ndarray[double, mode="c"] Rp not None,
        ndarray[int, mode="c"] no not None,
        ndarray[int, mode="c"] lo not None,
        ndarray[double, mode="c"] fo not None,
        ndarray[double, mode="c"] V_tot not None,
        ndarray[double, mode="c"] ks_energies not None,
        double reigen_eps, int reigen_max_iter, double mixing_eps,
        double mixing_alpha,
        int mixing_max_iter, bool perturb):
    cdef int n = len(R)
    cdef int n_orb = len(no)
    assert len(Rp) == n
    assert len(V_tot) == n
    assert len(lo) == n_orb
    assert len(fo) == n_orb
    assert len(ks_energies) == n_orb
    cdef double E_tot
    cdef ndarray[double, mode="c"] density = empty(n)
    cdef ndarray[double, ndim=2, mode="fortran"] orbitals = empty((n, n_orb),
            order="F")

    c_dftatom.dftatom_atom_lda2(&n, &n_orb, &Z, &R[0], &Rp[0],
            &no[0], &lo[0], &fo[0], &V_tot[0], &ks_energies[0], &E_tot,
            &density[0], &orbitals[0, 0],
            &reigen_eps, &reigen_max_iter, &mixing_eps, &mixing_alpha,
            &mixing_max_iter, &perturb)
    return E_tot, density, orbitals

def atom_rlda(int Z, double r_min, double r_max, double a, int N, double c,
        double reigen_eps, int reigen_max_iter,
        double mixing_eps, double mixing_alpha,
        int mixing_max_iter, bool perturb):
    cdef int n_orb
    cdef double E_tot
    c_dftatom.dftatom_get_atom_orb_rel(&Z, &n_orb)
    cdef ndarray[int, mode="c"] no = empty(n_orb, dtype="int32")
    cdef ndarray[int, mode="c"] lo = empty(n_orb, dtype="int32")
    cdef ndarray[int, mode="c"] so = empty(n_orb, dtype="int32")
    cdef ndarray[double, mode="c"] fo = empty(n_orb)
    cdef ndarray[double, mode="c"] ks_energies = empty(n_orb)
    cdef ndarray[double, mode="c"] R = empty(N+1)
    cdef ndarray[double, mode="c"] Rp = empty(N+1)
    cdef ndarray[double, mode="c"] V_tot = empty(N+1)
    cdef ndarray[double, mode="c"] density = empty(N+1)
    cdef ndarray[double, ndim=2, mode="fortran"] orbitals = empty((N+1, n_orb),
            order="F")
    c_dftatom.dftatom_atom_rlda(&Z, &r_min, &r_max, &a, &N, &c, &n_orb,
            &no[0], &lo[0], &so[0], &fo[0], &ks_energies[0], &E_tot,
            &R[0], &Rp[0], &V_tot[0], &density[0], &orbitals[0, 0], &reigen_eps,
            &reigen_max_iter,
            &mixing_eps, &mixing_alpha, &mixing_max_iter, &perturb)
    return E_tot, ks_energies, no, lo, so, fo, R, Rp, V_tot, density, orbitals

def mesh_exp(double r_min, double r_max, double a, int N):
    cdef ndarray[double, mode="c"] R = empty(N+1)
    c_dftatom.dftatom_mesh_exp(&r_min, &r_max, &a, &N, &R[0])
    return R

def mesh_exp_deriv(double r_min, double r_max, double a, int N):
    cdef ndarray[double, mode="c"] Rp = empty(N+1)
    c_dftatom.dftatom_mesh_exp_deriv(&r_min, &r_max, &a, &N, &Rp[0])
    return Rp

def integrate(ndarray[double, mode="c"] Rp not None,
        ndarray[double, mode="c"] f not None):
    cdef double s
    cdef int n = len(Rp)
    c_dftatom.dftatom_integrate(&n, &Rp[0], &f[0], &s)
    return s

def get_Vh(ndarray[double, mode="c"] R not None,
        ndarray[double, mode="c"] Rp not None,
        ndarray[double, mode="c"] rho not None):
    cdef ndarray[double, mode="c"] V = empty(len(R))
    cdef int n = len(R)
    c_dftatom.dftatom_get_vh(&n, &R[0], &Rp[0], &rho[0], &V[0])
    return V

import_array()
