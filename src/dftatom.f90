module dftatom

! This module contains the high level public API for dftatom. One should only
! be using this module from external programs (as long as only the high level
! functionality is needed). For a low level usage, one can always call the
! individual modules directly.

use types
use utils
use ode1d
use states
use mesh
use reigen
use dft
use rpoisson
use drivers
use energies
implicit none
private
public dp, E_nl, get_atomic_states_nonrel, mesh_exp, stop_error, &
        solve_radial_eigenproblem, get_tf_energies, thomas_fermi_potential, &
        get_atomic_states_rel, integrate, get_Vxc, str, mesh_exp_deriv, &
        get_atom_orb, atom_lda, atom_rlda, rpoisson_outward_pc, &
        get_atom_orb_rel, integrate_rproblem_outward, get_Vh

end module
