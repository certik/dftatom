program oscillator
use dftatom
implicit none

! Atomic number:
integer, parameter :: Z = 0
! Omega:
real(dp), parameter :: omega = 1.0_dp
! Mesh parameters:
real(dp), parameter :: r_min = 1e-8_dp, r_max = 10._dp, a = 80._dp
integer, parameter :: NN = 5000

real(dp), parameter :: c = 137.03599907_dp, eps = 1e-10_dp
integer :: n, n2, l, relat, converged, relat_max, k, kappa
real(dp) :: r(NN+1), u(size(r)), Ein, E, E_exact_nonrel, error, P(size(r)), &
    Q(size(r))
real(dp) :: Rp(NN+1)
real(dp) :: Emin_init, Emax_init
real(dp) :: E_exact(7, -7:7) ! E_exact(n, kappa)



E_exact = 0
E_exact(7, -7) = 7.4996755269_dp
E_exact(6:7, -6) = [ &
 6.4997620499_dp, &
 8.4994625812_dp &
]
E_exact(5:7, -5) = [ &
     5.4998352632_dp, &
     7.4995757160_dp, &
     9.4992363376_dp &
]
E_exact(4:7, -4) = [ &
     4.4998951661_dp, &
     6.4996755428_dp, &
     8.4993760822_dp, &
    10.4989967965_dp &
]
E_exact(3:7, -3) = [ &
     3.4999417582_dp, &
     5.4997620612_dp, &
     7.4995025208_dp, &
     9.4991631492_dp, &
    11.4987439585_dp &
]
E_exact(2:7, -2) = [ &
     2.4999750389_dp, &
     4.4998352706_dp, &
     6.4996156527_dp, &
     8.4993161976_dp, &
    10.4989369173_dp, &
    12.4984778242_dp &
]
E_exact(1:7, -1) = [ &
    1.4999950078_dp, &
    3.4998951705_dp, &
    5.4997154776_dp, &
    7.4994559413_dp, &
    9.4991165738_dp, &
    11.4986973873_dp, &
    13.4981983940_dp &
]
E_exact(2:7, 1) = [ &
     2.4999351051_dp, &
     4.4997953424_dp, &
     6.4995757249_dp, &
     8.4992762722_dp, &
    10.4988969952_dp, &
    12.4984379048_dp &
]
E_exact(3:7, 2) = [ &
     3.4998752033_dp, &
     5.4996955116_dp, &
     7.4994359765_dp, &
     9.4990966102_dp, &
    11.4986774249_dp &
]
E_exact(4:7, 3) = [ &
     4.4998019930_dp, &
     6.4995823772_dp, &
     8.4992829240_dp, &
    10.4989036457_dp &
]
E_exact(5:7, 4) = [ &
 5.4997154739_dp, &
 7.4994559363_dp, &
 9.4991165675_dp &
]
E_exact(6:7, 5) = [ &
 6.4996156467_dp, &
 8.4993161897_dp &
]
E_exact(7, 6) = 7.4995025118_dp



r = mesh_exp(r_min, r_max, a, NN)
Rp = mesh_exp_deriv(r_min, r_max, a, NN)
u = omega**2 * r**2 / 2

print *, "Relativistic linear spherically symmetric harmonic oscillator"
print *, "Mesh parameters (r_min, r_max, a, N):"
print "(ES10.2, F10.2, ES10.2, I10)", r_min, r_max, a, NN
print *
print *, " n  l  k kappa        E             E_exact     Error      E_exact_nonrel"
print *
do n = 1, 7
    do l = 0, n-1
        if (l == 0) then
            relat_max = 2
        else
            relat_max = 3
        end if
        do relat = 2, relat_max
            k = relat - 2
            if (k == 0) then
                kappa = -l - 1
            else
                kappa = l
            end if
            E_exact_nonrel = omega * (2*n - l - 1._dp/2)
            Ein = -100
            Emin_init = 0
            Emax_init = 100
            if (kappa > 0) then
                n2 = n + 1
            else
                n2 = n
            end if
            call solve_radial_eigenproblem(n2, l, Ein, eps, 100, R, Rp, u, &
                Z, c, relat, .false., Emin_init, Emax_init, converged, E, P, Q)
            error = abs(E - E_exact(n, kappa))
            if (converged /= 0) call stop_error("Not converged")
            print "(I3, I3, I3, I3, F17.8, F17.8, ES10.2, F17.8)", n, l, k, &
                kappa, E, E_exact(n, kappa), error, E_exact_nonrel
            if (error > 1e-8_dp) call stop_error("Error is higher than 1e-8")
        end do
    end do
end do
end program
