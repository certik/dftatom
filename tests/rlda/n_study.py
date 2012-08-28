from dftatom import atom_rlda, atom_lda
from dftatom.plots import plot_rmin_conv, plot_N_conv
from numpy import array, loadtxt

def N_study(Z, E_tot_ref, ks_energies_ref):
    c = 137.0359895
    rmax = 50
    rmin = 1e-8
    a = 6.2e7
    Nmax = 20000
    Nmin = 100
    graph = open("graph2.dat", "w")
    for N in range(16278, 20000, 100):
        E_tot, ks_energies, n, l, s, f, R, V_tot, density, orbitals = \
                atom_rlda(Z, rmin, rmax, a, N, c, 1e-10, 1e-11, 0.35, 100)
        E_tot_err = abs(E_tot - E_tot_ref)
        ks_energies_err = max(abs(ks_energies-ks_energies_ref))
        graph.write("%d %.17e %.17e\n" % (N, E_tot_err, ks_energies_err))
        print N, E_tot_err, ks_energies_err

def load_ref(filename):
    f = open(filename)
    data = {}
    line = f.readline()
    # skip headers
    while line.startswith("#"):
        line = f.readline()
    for Z in range(1, 93):
        assert Z == int(line)
        line = f.readline()
        E_tot = float(line)
        line = f.readline()
        ks_energies = array([float(x) for x in line.split()[1:]])
        data[Z] = {"E_tot": E_tot, "ks_energies": ks_energies}
        line = f.readline()
    return data

ref = load_ref("rel_energies.dat")
Z = 90
E_tot_ref = ref[Z]["E_tot"]
ks_energies_ref = ref[Z]["ks_energies"]
N_study(Z, E_tot_ref, ks_energies_ref)
