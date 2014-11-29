module dft

! Calculates the exchange and correlation potential, Hartree potential,
! and the full (single) Kohn Sham iteration.

use types, only: dp
use constants, only: pi
use dft_data, only: dft_data_t
use ode1d, only: integrate
use reigen, only: solve_radial_eigenproblem
use utils, only: stop_error
use energies, only: get_tf_energies, thomas_fermi_potential, &
    thomas_fermi_density
use rpoisson, only: rpoisson_outward_pc

use mixings, only: mixing_anderson
implicit none
private
public get_Vxc, get_Vh, KS_step, rho2V

contains

subroutine rho2V(d)
! Calculates V_xc, V_h and V_tot from rho.
! Assumes that d%rho is normalized.
type(dft_data_t), intent(inout) :: d

call get_Vxc(d%R, d%rho, d%dirac, d%c, d%e_xc, d%V_xc)
d%V_h = get_Vh(d%R, d%Rp, d%rho)
call total_energy(d%fo, d%ks_energies, d%V_tot, d%V_h, d%V_coulomb, d%e_xc, &
    d%R, d%Rp, d%rho, d%Ekin, d%Ecoul, d%Eenuc, d%Exc, d%Etot)
d%V_tot = d%V_coulomb + d%V_xc + d%V_h
end subroutine

subroutine V2rho(d)
! Calculates rho from V by solving Kohn-Sham equations
type(dft_data_t), intent(inout) :: d

real(dp), dimension(size(d%R)) :: P, Q, Y
integer :: converged, i, n, l, relat
real(dp) :: Ein, Emin_init, Emax_init

d%rho(:) = 0
!print *, d%scf_iter, d%ks_energies
do i = 1, size(d%no)
    n = d%no(i)
    l = d%lo(i)
    if (d%dirac) then
        if (d%so(i) == 1) then
            relat = 2
        else
            relat = 3
        end if
    else
        relat = 0
    end if
    Ein = d%ks_energies(i)
    Emax_init = d%Emax_init(i)
    Emin_init = d%Emin_init(i)

    call solve_radial_eigenproblem(n, l, Ein, d%reigen_eps, &
        d%reigen_max_iter, &
        d%R, d%Rp, d%V_tot, &
        d%Z, d%c, relat, d%perturb, Emin_init, Emax_init, &
        converged, d%ks_energies(i), P, Q)
    if (converged /= 0) then
        print *, "converged=", converged
        print *, d%scf_iter, n, l, relat
        !print *, "skipping the state"
        !Y = 0
        call stop_error("V2rho: Radial eigen problem didn't converge")
    end if
    if (relat == 0) then
        Y = P / d%R
    else
        Y = sqrt(P**2 + Q**2) / d%R
    end if
    d%rho = d%rho + d%fo(i) * Y**2
    d%orbitals(:, i) = Y
end do
d%rho = d%rho / (4*pi)
end subroutine

function KS_step(V, i, d) result(F)
! Calculates the residual R(V) = KS(V) - V
real(dp), intent(in) :: V(:)
integer, intent(in) :: i
type(dft_data_t), intent(inout) :: d
real(dp) :: F(size(V))

d%scf_iter = i
! We converge upon the non-Coulombic part of the potential:
d%V_tot = d%V_coulomb + V
call V2rho(d)
call rho2V(d)
F = d%V_xc + d%V_h - V
end function

subroutine total_energy(fo, ks_energies, V_in, V_h, V_coulomb, e_xc, R, Rp, n, &
        T_s, E_ee, E_en, EE_xc, Etot)
! This is a variational, quadratically convergent form of total energy
real(dp), intent(in) :: fo(:), ks_energies(:) ! occupations, energies
real(dp), intent(in) :: V_in(:) ! Total input effective potential
real(dp), intent(in) :: V_h(:) ! Hartree energy, solution of Poiss. eq.
real(dp), intent(in) :: V_coulomb(:) ! Coulomb inter. -Z/r  (negative)
real(dp), intent(in) :: e_xc(:) ! XC density
real(dp), intent(in) :: R(:), Rp(:), n(:) ! Radial grid, number density (positive)
real(dp), intent(out) :: Etot ! Total energy
real(dp), intent(out) :: T_s, E_ee, E_en, EE_xc ! Parts of the total energy

real(dp) :: rho(size(n))
real(dp) :: E_c, E_band!, Exc2
rho = -n

E_band = sum(fo * ks_energies)
T_s = E_band + 4*pi * integrate(Rp, V_in * rho * R**2)

E_ee = -2*pi * integrate(Rp, V_h * rho * R**2)
E_en =  4*pi * integrate(Rp, (-V_coulomb) * rho * R**2)
E_c = E_ee + E_en

EE_xc = -4*pi * integrate(Rp, e_xc * rho * R**2)

Etot = T_s + E_c + EE_xc
end subroutine

subroutine get_Vxc(R, rho, relat, c, exc, Vxc)
real(dp), intent(in) :: R(:) ! radial grid
real(dp), intent(in) :: rho(:) ! charge density
logical, intent(in) :: relat ! .true. return RLDA, otherwise LDA
real(dp), intent(in) :: c ! speed of light
real(dp), intent(out) :: Vxc(:), exc(:)

integer :: i
do i = 1, size(R)
    call getvxc_scalar(rho(i), relat, c, exc(i), Vxc(i))
end do
end subroutine

function get_Vh(R, Rp, rho) result(V)
real(dp), intent(in) :: R(:), Rp(:), rho(:)
real(dp) :: V(size(R))
V = rpoisson_outward_pc(R, Rp, rho)
end function

real(dp) function get_Y(y, b, c)
real(dp), intent(in) :: y, b, c
get_Y = y**2 + b*y + c
end function

subroutine getvxc_scalar(n, relat, c_light, exc, Vxc)
! Calculates XC LDA density and potential from the charge density "n".
real(dp), intent(in) :: n ! charge density (scalar)
real(dp), intent(in) :: c_light ! speed of light
logical, intent(in) :: relat ! if .true. returns RLDA, otherwise LDA
real(dp), intent(out) :: exc ! XC density
real(dp), intent(out) :: Vxc ! XC potential

real(dp), parameter :: y0 = -0.10498_dp
real(dp), parameter :: b = 3.72744_dp
real(dp), parameter :: c = 12.9352_dp
real(dp), parameter :: A = 0.0621814_dp

real(dp) :: Q, rs, y, ec, ex, Vc, Vx, beta, mu, R, S

if (n < epsilon(1._dp)) then
    exc = 0
    Vxc = 0
    return
end if

Q = sqrt(4*c - b**2)
rs = (3/(4*pi*n))**(1.0_dp/3)
y = sqrt(rs)
ec = A/2 * (log(y**2/get_Y(y, b, c)) + 2*b/Q * atan(Q/(2*y+b))  &
   - b*y0/get_Y(y0, b, c) * ( &
            log((y-y0)**2 / get_Y(y, b, c)) &
            + 2*(b+2*y0) / Q * atan(Q/(2*y+b)) &
          ) )
Vc = ec - A/6 * (c*(y-y0)-b*y0*y)/((y-y0)*get_Y(y, b, c))
ex = -3/(4*pi) * (3*pi**2*n)**(1.0_dp/3)
Vx = 4*ex/3

if (relat) then
    beta = -4 * pi * ex / (3 * c_light)
    mu = sqrt(1 + beta**2)
    R = 1 - 3 * ((beta * mu - log(beta + mu)) / (beta ** 2))**2 / 2
    S = 3 * log(beta + mu) / (2 * beta * mu) - 1.0_dp/2

    ex = ex * R
    Vx = Vx * S
end if
exc = ex + ec
Vxc = Vx + Vc
end subroutine

end module
