program bench_poisson
use types, only: dp
use mesh, only: mesh_exp, mesh_exp_deriv
use rpoisson, only: rpoisson_kernel1, rpoisson_kernel2
implicit none

real(dp), parameter :: r_min = 1e-7_dp, r_max = 10, a = 2.7e6_dp
real(dp), allocatable, dimension(:) :: R, Rp, u1, u2, u1p, u2p, rho
real(dp) :: t1, t2, dt1, dt2
integer :: N, i, iter

print *, "N iter dt1[s] dt2[s] total_time[s] total_mem[MB]"
do N = 5000, 50000, 5000
    allocate(R(N), Rp(N), u1(N), u2(N), u1p(N), u2p(N), rho(N))
    R = mesh_exp(r_min, r_max, a, N-1)
    Rp = mesh_exp_deriv(r_min, r_max, a, N-1)
    rho = exp(-R)
    u1 = 1
    u2 = 1
    u1p = 0
    u2p = 0
    iter = 4000 * 5000/N
    call cpu_time(t1)
    do i = 1, iter
        call rpoisson_kernel1(R, Rp, rho, u1, u2, u1p, u2p)
    end do
    call cpu_time(t2)
    dt1 = (t2-t1)/iter
    call cpu_time(t1)
    do i = 1, iter
        call rpoisson_kernel2(R, Rp, rho, u1, u2, u1p, u2p)
    end do
    call cpu_time(t2)
    dt2 = (t2-t1)/iter
    print "(i6, i6, es12.4, es12.4, f6.2, f6.2)", N, iter, dt1, dt2, dt1*iter, &
        N*8*7/1024._dp**2
    deallocate(R, Rp, u1, u2, u1p, u2p, rho)
end do
end program
