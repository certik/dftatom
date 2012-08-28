from dftatom import atom_rlda, atom_lda
from dftatom.plots import plot_rmin_conv, plot_N_conv
from numpy import array, loadtxt

def get_energies_N(Z, E_tot_ref, ks_energies_ref, eps):
    c = 137.0359895
    rmax = 50
    rmin = 1e-8
    a = 6.2e7
    Nmax = 20000
    Nmin = 100
    while 1:
        N = (Nmax + Nmin) / 2
        E_tot, ks_energies, n, l, s, f, R, V_tot, density, orbitals = \
                atom_rlda(Z, rmin, rmax, a, N, c, 1e-10, 100, 5e-9, 0.5, 100)
        E_tot_err = abs(E_tot - E_tot_ref)
        ks_energies_err = max(abs(ks_energies-ks_energies_ref))
        print Nmin, Nmax, E_tot_err, ks_energies_err
        if E_tot_err < eps and ks_energies_err < eps:
            Nmax = N
            E_tot_ok = E_tot
            ks_energies_ok = ks_energies
        else:
            Nmin = N
        if (Nmax - Nmin <= 1):
            break
    return E_tot_ok, ks_energies_ok, Nmax

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
for p in [3, 4, 5, 6, 7, 8]:
    f = open("data%d.dat" % p, "w")
    r = range(1, 93)
    r.reverse()
    for Z in r:
        E_tot_ref = ref[Z]["E_tot"]
        ks_energies_ref = ref[Z]["ks_energies"]
        # Do the precision a little bit below the limit, to take into account
        # machine accuracy problems:
        eps = 0.9 * 10**(-p)
        E_tot, ks_energies, N = get_energies_N(Z, E_tot_ref, ks_energies_ref,
                eps)
        f.write("%d %d\n" % (Z, N))
        print "Z = %2d  N = %4d  E_tot_err = %e" %(Z, N, abs(E_tot-E_tot_ref))
        print "    max KS error: %e" % max(abs(ks_energies - ks_energies_ref))
