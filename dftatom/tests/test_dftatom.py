# These tests test the Python interface (that all the types work as expected).
# The accuracy and correctness of the actual solver is tested in the Fortran
# tests.
#
# Another purpose of these tests is to show how to use the Python API.

from dftatom import (mesh_exp, mesh_exp_deriv, solve_radial_eigenproblem,
        atom_lda, atom_rlda)

def test_nonrelat():
    a = 2.7e6
    rmin = 1e-7
    rmax = 50
    N = 1000
    Z = 92
    relat = 0
    R = mesh_exp(rmin, rmax, a, N)
    Rp = mesh_exp_deriv(rmin, rmax, a, N)
    V_tot = -Z / R
    energies = []
    for n in range(8):
        for l in range(n):
            E, P, Q = solve_radial_eigenproblem(0.0, n, l, -1, 1e-11, 100,
                    V_tot, R, Rp, Z, relat, False, -10000, 0)
            energies.append(E)

def test_relat():
    a = 6.2e7
    rmin = 1e-8
    rmax = 50
    N = 1000
    Z = 92
    relat = 0
    c = 137.0359895
    R = mesh_exp(rmin, rmax, a, N)
    Rp = mesh_exp_deriv(rmin, rmax, a, N)
    V_tot = -Z / R
    energies = []
    for n in range(8):
        for l in range(n):
            relats = [2]
            if l > 0:
                relats.append(3)
            for relat in relats:
                E, P, Q = solve_radial_eigenproblem(c, n, l, relat, 1e-11, 100,
                        V_tot, R, Rp, Z, relat, False, -10000, 0)
                energies.append(E)

def test_lda():
    a = 2.7e6
    rmin = 1e-7
    rmax = 50
    Z = 92
    N = 1000
    E_tot, ks_energies, n, l, f, R, Rp, V_tot, density, orbitals = \
            atom_lda(Z, rmin, rmax, a, N, 1e-11, 100, 1e-10, 0.35, 100, True)

def test_rlda():
    a = 6.2e7
    rmin = 1e-8
    rmax = 50
    Z = 92
    N = 1000
    c = 137.0359895
    E_tot, ks_energies, n, l, s, f, R, Rp, V_tot, density, orbitals = \
            atom_rlda(Z, rmin, rmax, a, N, c, 1e-11, 100, 1e-10, 0.35, 100,
                    True)
