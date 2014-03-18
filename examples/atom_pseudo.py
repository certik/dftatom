from numpy import array, empty, loadtxt
from scipy.interpolate import splrep, splev
from dftatom import atom_lda_pseudo, mesh_exp, mesh_exp_deriv

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
data = loadtxt("sn-pseudo.txt")
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
