program N_conv

! Calculates N convergence for the lowest 26 eigenvalues of the double
! minimum potential. Results are saved to 'energies_N.txt'.

use dftatom
use utils, only: savetxt
use interpolation, only: spline3, loadtxt
implicit none

! For non-singular potential, we need to set Z=0:
integer, parameter :: Z = 0
! Mesh parameters:
real(dp), parameter :: r_max = 12._dp, a = 2._dp
integer :: NN = 50000

real(dp), parameter :: c = 137.035999037_dp, eps = 1e-12_dp
integer :: n, l, relat, converged, i, maxiter
real(dp) :: Ein, E
real(dp) :: Emin_init, Emax_init
real(dp), allocatable :: data(:, :), Eall(:, :), r(:), Rp(:), u(:), P(:), Q(:)
real(dp) :: E0, au, mu
integer, parameter :: maxn = 26
real(dp), parameter :: r_min(maxn) = [ &
    0.679_dp, &
    0.991_dp, &
    0.932_dp, &
    0.768_dp, &
    0.768_dp, &
    0.738_dp, &
    0.649_dp, &
    0.649_dp, &
    0.634_dp, &
    0.619_dp, &
    0.604_dp, &
    0.604_dp, &
    0.589_dp, &
    0.589_dp, &
    0.574_dp, &
    0.574_dp, &
    0.559_dp, &
    0.559_dp, &
    0.544_dp, &
    0.544_dp, &
    0.544_dp, &
    0.530_dp, &
    0.530_dp, &
    0.530_dp, &
    0.530_dp, &
    0.515_dp]

call loadtxt("potential_table.txt", data)

maxiter = 50
allocate(Eall(0:maxn, maxiter))
do i = 1, maxiter
    NN = 300 + (50000-300) * (i-1)/(maxiter-1)
    allocate(r(NN+1), Rp(NN+1), u(NN+1), P(NN+1), Q(NN+1))
    Eall(0, i) = NN
    au = 219474.62_dp
    E0 = -0.625_dp
    mu = 1836.12_dp / 2
    print *, "NN = ", NN
    print *, "Double minimum potential"
    print *, "Mesh parameters (r_max, a, N):"
    print "(F10.2, ES10.2, I10)", r_max, a, NN
    print *
    print "(A3, A18, A18)", "n", "E [a.u.]", "D [cm^-1]"
    print *
    do n = 1, maxn
        R = mesh_exp(r_min(n), r_max, a, NN)
        Rp = mesh_exp_deriv(r_min(n), r_max, a, NN)
        u = spline3(data(:, 1), data(:, 2), R)
        u = u * mu
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
    deallocate(r, Rp, u, P, Q)
end do
call savetxt("energies_N.txt", Eall)
end program
