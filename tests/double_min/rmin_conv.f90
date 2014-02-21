program rmin_conv

! Calculates rmin convergence for the lowest 26 eigenvalues of the double
! minimum potential. Results are saved to 'energies_rmin.txt'.

use dftatom
use utils, only: savetxt
use interpolation, only: spline3, loadtxt
implicit none

! For non-singular potential, we need to set Z=0:
integer, parameter :: Z = 0
! Mesh parameters:
real(dp), parameter :: r_max = 12._dp, a = 2._dp
integer, parameter :: NN = 50000

real(dp), parameter :: c = 137.035999037_dp, eps = 1e-12_dp
integer :: n, l, relat, converged, i, maxiter, maxn
real(dp) :: r(NN+1), u(size(r)), Ein, E, P(size(r)), Q(size(r))
real(dp) :: Rp(NN+1)
real(dp) :: Emin_init, Emax_init
real(dp), allocatable :: data(:, :), Eall(:, :)
real(dp) :: E0, au, mu
real(dp) :: r_min

call loadtxt("potential_table.txt", data)

maxn = 26
maxiter = 50
allocate(Eall(0:maxn, maxiter))
do i = 1, maxiter
    r_min = 1.2_dp - (1.2_dp-0.47_dp) * (i-1)/(maxiter-1)
    Eall(0, i) = r_min

    R = mesh_exp(r_min, r_max, a, NN)
    Rp = mesh_exp_deriv(r_min, r_max, a, NN)
    u = spline3(data(:, 1), data(:, 2), R)

    au = 219474.62_dp
    E0 = -0.625_dp
    mu = 1836.12_dp / 2
    u = u * mu
    print *, "rmin = ", r_min
    print *, "Double minimum potential"
    print *, "Mesh parameters (r_min, r_max, a, N):"
    print "(ES10.2, F10.2, ES10.2, I10)", r_min, r_max, a, NN
    print *
    print "(A3, A18, A18)", "n", "E [a.u.]", "D [cm^-1]"
    print *
    do n = 1, maxn
        l = 0
        relat = 0
        Ein = -0.5
        Emin_init = -100*mu
        Emax_init = 100
        call solve_radial_eigenproblem(n, l, Ein, eps, 100, r, rp, u, &
            Z, c, relat, .false., Emin_init, Emax_init, converged, E, P, Q)
        if (converged /= 0) then
            print "(I3, A18, A18)", n-1, "0", "0"
            Eall(n, i) = 0
        else
            E = E / mu
            Eall(n, i) = E
            print "(I3, F18.6, F18.6)", n-1, E, (E0-E)*au
        end if
    end do
end do
call savetxt("energies_rmin.txt", Eall)
end program
