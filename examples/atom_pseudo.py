# First create 'sn.param', then run opium as follows:
# opium sn.param sn.log all rpt
# Then run this script, which reads "sn.vi_plt" below

from math import pi
from numpy import array, empty, fromstring, size
from scipy.interpolate import splrep, splev
from dftatom import (atom_lda_pseudo, mesh_exp, mesh_exp_deriv, integrate,
        atom_lda)

def read_plt(filename):
    t = open(filename).read().split("@")
    V = []
    for n in range(4):
        V.append(fromstring(t[n], sep=" ").reshape([-1, 3]))
    V = array(V)
    # For the potential number n=0, 1, 2, 3:
    # R    = V[n, :, 0]
    # V(R) = V[n, :, 1]
    # V[n, :, 2] is radius

    # Convert Ry to Ha (atomic units):
    for n in range(4):
        V[n, :, 1] /= 2

    R = V[0, :, 0]

    return array([R, V[0, :, 1], V[1, :, 1], V[2, :, 1]])


a = 2.7e6
rmin = 1e-7
rmax = 50
N = 5000
R = mesh_exp(rmin, rmax, a, N)
Rp = mesh_exp_deriv(rmin, rmax, a, N)

no = array([1, 2, 3], dtype="int32")
lo = array([0, 1, 2], dtype="int32")
fo = array([2, 2, 10], dtype="double")
Emax_init = empty(3)
Emax_init[:] = 10
Emin_init = empty(3)
Emin_init[:] = -30
ks_energies = array([-2, -1, -1], dtype="double")
data = read_plt("sn.vi_plt")
V_l = empty((N+1, 3), order='F')
for i in range(3):
    tck = splrep(data[0, :], data[i+1, :])
    V_l[:, i] = splev(R, tck)
V_loc = V_l[:, 0].copy()
for i in range(3):
    V_l[:, i] -= V_loc
V_tot = V_loc.copy()
xc_type = 2
reigen_eps = 1e-10
mixing_eps = 5e-9
mixing_alpha = 0.5
mixing_max_iter = 200
reigen_max_iter = 100
perturb = True


Ekin, Eee, Een, Exc, Etot, Enl, density, orbitals = \
        atom_lda_pseudo(no, lo, fo, Emin_init, Emax_init, ks_energies, R, Rp,
                V_loc, V_l, V_tot, reigen_eps, reigen_max_iter, mixing_eps,
                mixing_alpha, mixing_max_iter, perturb, xc_type)

print "Pseudo:"
print "E_tot = %13.6f" % Etot
print "E_nl  = %13.6f a.u. = %13.6f Ry" % (Enl, Enl*2)
print "n l   f       E"
for n_, l_, f_, E in zip(no, lo, fo, ks_energies):
    print "%d %d %5.2f %13.6f" % (n_, l_, f_, E)

print
print "The first 10 values of the radial grid:"
print R[:10]

print
print "The first 10 values of the 1st and 2nd orbitals:"
print orbitals[:10, 0]
print orbitals[:10, 1]

print
print "The first 10 values of the radial charge density:"
print density[:10]

print
print "Total charge:", integrate(Rp, 4*pi*density*R**2)
idx = size(R)-1
while R[idx] > 2.0: idx -= 1
print "Rcut =", R[idx], "a.u. = R(idx); idx =", idx
print "Total charge within cut-off:", integrate(Rp[:idx],
        4*pi*density[:idx]*R[:idx]**2)

print
print "ALL ELECTRON CALCULATION:"
Z = 50
E_tot, ks_energies, n, l, f, R, Rp, V_tot, density, orbitals = \
        atom_lda(Z, rmin, rmax, a, N, 1e-11, 100, 1e-10, 0.35, 100, True)
valence_density = 0
print "Only treating the following valence states in charge calculation:"
for i in [8, 9, 10]: # valence states, starting from 0
    print i, n[i], l[i], f[i]
    # Construct the density corresponding to the orbital 'i':
    valence_density += f[i]*orbitals[:, i]**2/(4*pi)

density = valence_density
print "Total charge:", integrate(Rp, 4*pi*density*R**2)
print "Total charge within cut-off:", integrate(Rp[:idx],
        4*pi*density[:idx]*R[:idx]**2)
