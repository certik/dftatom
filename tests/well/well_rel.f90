program well
use special, only: spherical_bessel_jn_zeros
use dftatom
implicit none

! Atomic number:
integer, parameter :: Z = 0
! Mesh parameters:
real(dp), parameter :: r_min = 1e-8_dp, r_max = 10, a = 6.2e7_dp
integer, parameter :: NN = 6000

real(dp), parameter :: c = 137.03599907_dp, eps = 1e-10_dp
integer :: n, n2, l, relat, converged, relat_max, k, kappa
real(dp) :: r(NN+1), u(size(r)), Ein, E, E_exact_nonrel, error, P(size(r)), &
    Q(size(r))
real(dp) :: Rp(NN+1)
real(dp) :: Emin_init, Emax_init
real(dp), allocatable :: zeros(:, :)
real(dp) :: E_exact(7, -7:7) ! E_exact(n, kappa)



E_exact = 0
E_exact(7, -7) = 0.5525904113_dp
E_exact(6:7, -6) = [ &
 0.4376510015_dp, &
 0.8406357081_dp &
]
E_exact(5:7, -5) = [ &
 0.3347685757_dp, &
 0.6850117637_dp, &
 1.1309235187_dp &
]
E_exact(4:7, -4) = [ &
 0.2441543810_dp, &
 0.5425739559_dp, &
 0.9381557573_dp, &
 1.4319901885_dp &
]
E_exact(3:7, -3) = [ &
 0.1660865751_dp, &
 0.4135916009_dp, &
 0.7592590218_dp, &
 1.2034759695_dp, &
 1.7463192010_dp &
]
E_exact(2:7, -2) = [ &
 0.1009533714_dp, &
 0.2983952090_dp, &
 0.5944899358_dp, &
 0.9892629990_dp, &
 1.4827135257_dp, &
 2.0748352993_dp &
]
E_exact(1:7, -1) = [ &
 0.0493479580_dp, &
 0.1973910436_dp, &
 0.4441269495_dp, &
 0.7895517628_dp, &
 1.2336600278_dp, &
 1.7764447596_dp, &
 2.4178974112_dp &
]
E_exact(2:7, 1) = [ &
 0.1009533712_dp, &
 0.2983952108_dp, &
 0.5944899397_dp, &
 0.9892630035_dp, &
 1.4827135323_dp, &
 2.0748353048_dp &
]
E_exact(3:7, 2) = [ &
 0.1660865751_dp, &
 0.4135916010_dp, &
 0.7592590218_dp, &
 1.2034759695_dp, &
 1.7463192010_dp &
]
E_exact(4:7, 3) = [ &
 0.2441543811_dp, &
 0.5425739559_dp, &
 0.9381557573_dp, &
 1.4319901885_dp &
]
E_exact(5:7, 4) = [ &
 0.3347685757_dp, &
 0.6850117636_dp, &
 1.1309235187_dp &
]
E_exact(6:7, 5) = [ &
 0.4376510015_dp, &
 0.8406357081_dp &
]
E_exact(7, 6) = 0.5525904113_dp



r = mesh_exp(r_min, r_max, a, NN)
Rp = mesh_exp_deriv(r_min, r_max, a, NN)
u = 0

allocate(zeros(10, 0:7))
zeros = spherical_bessel_jn_zeros(7, 10, 1e-13_dp)
print *, "Relativistic infinite spherical potential well"
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
            E_exact_nonrel = (zeros(n-l, l) / r_max)**2 / 2
            k = relat - 2
            if (k == 0) then
                kappa = -l - 1
            else
                kappa = l
            end if
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
            if (error > 1e-6_dp) call stop_error("Error is higher than 1e-6")
        end do
    end do
end do
end program
