# This example shows how to calculate nonrelativistic and relativistic Uranium
# To 1e-6 Ha accuracy in total energy and eigenvalues.

from __future__ import print_function
from dftatom import atom_lda, atom_rlda

a = 2.7e6
rmin = 1e-7
rmax = 50
Z = 92
N = 5500
E_tot, ks_energies, n, l, f, R, Rp, V_tot, density, orbitals = \
        atom_lda(Z, rmin, rmax, a, N, 1e-11, 100, 1e-10, 0.35, 100, True)

print("Schroedinger LDA:")
print("E_tot = %13.6f" % E_tot)
print("n l   f       E")
for n_, l_, f_, E in zip(n, l, f, ks_energies):
    print("%d %d %5.2f %13.6f" % (n_, l_, f_, E))
print()

a = 6.2e7
rmin = 1e-8
c = 137.0359895
E_tot, ks_energies, n, l, s, f, R, Rp, V_tot, density, orbitals = \
        atom_rlda(Z, rmin, rmax, a, N, c, 1e-11, 100, 1e-10, 0.35, 100,
                True)

print("Dirac RLDA:")
print("E_tot = %13.6f" % E_tot)
print("n l s   f       E")
for n_, l_, s_, f_, E in zip(n, l, s, f, ks_energies):
    print("%d %d %d %5.2f %13.6f" % (n_, l_, s_, f_, E))
