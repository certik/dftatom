program uranium
! Calculates Hydrogen like relativistic energies for Z=92 (U)
!
! The purpose of this test is to check that:
! a) all U energies can converge to 1e-6
! b) that it converges when using very rough initial estimate
use dftatom
implicit none

! Atomic number:
integer, parameter :: Z = 92
! Mesh parameters:
real(dp), parameter :: r_min = 1e-8_dp, r_max = 50.0_dp, a = 6.2e7_dp
integer, parameter :: NN = 10000

real(dp), parameter :: c = 137.035999037_dp, eps = 1e-6_dp
integer :: n, l, relat, converged, relat_max
real(dp) :: r(NN+1), u(size(r)), Ein, E, E_exact, error, P(size(r)), Q(size(r))
real(dp) :: Rp(NN+1)


r = mesh_exp(r_min, r_max, a, NN)
Rp = mesh_exp_deriv(r_min, r_max, a, NN)
u(:) = -Z/r

print *, "Hydrogen like relativistic energies for Z=92 (U)"
print *, "Mesh parameters (r_min, r_max, a, N):"
print "(ES10.2, F10.2, ES10.2, I10)", r_min, r_max, a, NN
print *
print "(A3, A3, A3, A15, A15, A10)", "n", "l", "k", "E", "E_exact", "Error"
print *
do n = 1, 16
    do l = 0, n-1
        if (l == 0) then
            relat_max = 2
        else
            relat_max = 3
        end if
        do relat = 2, relat_max
            E_exact = E_nl(c, n, l, Z, relat)
            Ein = -100
            call solve_radial_eigenproblem(n, l, Ein, eps, 100, R, Rp, u, &
                Z, c, relat, .true., -10000._dp, 0._dp, converged, E, P, Q)
            error = abs(E - E_exact)
            if (converged /= 0) then
                print *, n, l, relat
                print *, "converged=", converged
                call stop_error("Not converged")
            end if
            print "(I3, I3, I3, F15.6, F15.6, ES10.2)", n, l, relat-2, &
                E, E_exact, error
            if (error > eps) call stop_error("Error is higher than 1e-6")
        end do
    end do
end do
end program
