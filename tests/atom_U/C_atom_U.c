#include <stdio.h>

#include "c_dftatom.h"

int main()
{
    int Z, NN, n_orb;
    int iter1=100, iter2=200;
    double alpha=0.5;
    bool perturb=true;
    double r_min, r_max, a, E_tot, reigen_eps, mixing_eps;
    char l_names[4] = {'s', 'p', 'd', 'f'};
    int i;
    Z = 92;
    NN = 5500;
    dftatom_get_atom_orb(&Z, &n_orb);
    double ks_energies[n_orb], fo[n_orb],
           orbitals[n_orb][NN+1], R[NN+1], V_tot[NN+1], density[NN+1],
           Rp[NN+1];
    int no[n_orb], lo[n_orb];
    r_min = 1e-7;
    r_max = 50.0;
    a = 2.7e6;
    reigen_eps = 1e-10;
    mixing_eps = 5e-9;
    dftatom_atom_lda(&Z, &r_min, &r_max, &a, &NN, &n_orb, no, lo, fo,
            ks_energies, &E_tot, R, Rp, V_tot, density, &orbitals[0][0],
            &reigen_eps, &iter1, &mixing_eps, &alpha, &iter2, &perturb);
    printf("Z=%d, N=%d\n", Z, NN);
    printf("E_tot=%18.6f\n", E_tot);
    printf(" state    E            occupancy\n");
    for (i=0; i < n_orb; i++) {
        printf("%d%c %18.6f   %6.3f\n", no[i], l_names[lo[i]], ks_energies[i],
                fo[i]);
    }
    printf("\nPrint the first 10 values of the 1st and 2nd orbitals:\n");
    for (i=0; i < 10; i++) printf("%f   ", orbitals[0][i]);
    printf("\n");
    for (i=0; i < 10; i++) printf("%f   ", orbitals[1][i]);
    printf("\n");
    return 0;
}
