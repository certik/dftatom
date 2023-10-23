module energies

! Calculates Hydrogen nonrelativistic and relativistic energies (exact),
! Thomas-Fermi (TF) energies (only very approximate), TF potential and charge
! density (very accurate).

use types, only: dp
use utils, only: stop_error
use constants, only: pi
implicit none
private
public E_nl, thomas_fermi_potential, get_tf_energies, get_hydrogen_energies, &
    thomas_fermi_density

contains

real(dp) function E_nl(c, n, l, Z, relat)
! Calculates exact energy for the radial Schroedinger/Dirac equations
real(dp), intent(in) :: c ! speed of light in atomic units
integer, intent(in) :: n, l, Z, relat
! quantum numbers (n, l), atomic number (z)
! relat == 0 ... Schroedinger equation
! relat == 2 ... Dirac equation, spin up
! relat == 3 ... Dirac equation, spin down

integer :: kappa
real(dp) :: beta
if (.not. (l >= 0)) call stop_error("'l' must be positive or zero")
if (.not. (n > l)) call stop_error("'n' must be greater than 'l'")
if (l == 0 .and. relat == 3) call stop_error("Spin must be up for l==0.")
if (relat == 0) then
    E_nl = - Z**2 / (2.0_dp * n**2)
else
    if (relat == 2) then
        kappa = -l - 1
    else
        kappa = l
    end if
    beta = sqrt(kappa**2 - (Z/c)**2)
    E_nl = c**2/sqrt(1 + (Z/c)**2/(n - abs(kappa) + beta)**2) - c**2
end if
end function

function thomas_fermi_potential(R, Z, cut) result(V)
! Generalized Thomas-Fermi atomic potential
real(dp), intent(in) :: R(:) ! Radial grid
integer, intent(in) :: Z     ! Atomic number
logical, intent(in), optional :: cut ! Cut the potential, default .true.
real(dp) :: x(size(R)), Z_eff(size(R)), V(size(R))
real(dp) :: alpha, beta, gamma

x = R * (128*Z/(9*pi**2)) ** (1.0_dp/3)
! Z_eff(x) = Z * phi(x), where phi(x) satisfies the Thomas-Fermi equation:
!   phi'' = phi**(3/2) / sqrt(x)
! with boundary conditions:
!   phi(0)  = 1
!   phi(oo) = 0
! There is no analytic solution, but one can solve this approximately. We use:
! http://arxiv.org/abs/physics/0511017
alpha = 0.7280642371_dp
beta = -0.5430794693_dp
gamma = 0.3612163121_dp
Z_eff = Z * (1 + alpha*sqrt(x) + beta*x*exp(-gamma*sqrt(x)))**2 * &
    exp(-2*alpha*sqrt(x))
! This keeps all the eigenvalues of the radial problem negative:
if (.not. present(cut)) where (Z_eff < 1) Z_eff = 1
V = -Z_eff / r
end function

function thomas_fermi_density(R, Z) result(rho)
! Generalized Thomas-Fermi atomic potential
real(dp), intent(in) :: R(:) ! Radial grid
integer, intent(in) :: Z     ! Atomic number
real(dp) :: V(size(R)), rho(size(R))
V = thomas_fermi_potential(R, Z, .false.)
rho = -1 / (3*pi**2) * (-2*V)**(3._dp/2)
end function

function get_tf_energies(Z, no, fo) result(E)
integer, intent(in) :: Z, no(:)
real(dp), intent(in) :: fo(:)
real(dp) :: E(size(no))

integer :: Zeff, i
Zeff = Z + 1
do i = 1, size(no)
    if (i > 1) then
        if (no(i) < no(i-1)) call stop_error("State order wrong")
    end if
    Zeff = Zeff - int(fo(i))
    if (Zeff <= 0) call stop_error("Negative ions not allowed")
    E(i) = -(1.0_dp * Zeff / no(i))**2/2
end do
end function

function get_hydrogen_energies(Z, no) result(E)
integer, intent(in) :: Z, no(:)
real(dp) :: E(size(no))
integer :: i
do i = 1, size(no)
    E(i) = -1.0_dp * Z**2 / (2*no(i)**2)
end do
end function

end module
