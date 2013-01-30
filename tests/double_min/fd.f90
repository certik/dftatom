program fd

! Calculates the the lowest 26 eigenvalues of the double minimum potential
! using the finite difference (FD) method as described in [1].
! The interpolated potential is saved to 'potential_grid.txt'.

! [1] Tobin, F. L., & Hinze, J. (1975). The eigenvalue problem for a double
! minimum potential. The Journal of Chemical Physics, 63(2), 1034.
! doi:10.1063/1.431399

use types, only: dp
use utils, only: stop_error, newunit, savetxt
use mesh, only: mesh_exp
use interpolation, only: spline3, lagrange3interp, hermite5interp, loadtxt
use lapack, only: dstevd
implicit none
real(dp) :: rmin, rmax
integer :: N
real(dp), allocatable :: xe(:), lam(:), U(:)
integer :: i
real(dp) :: mu, E0, au
real(dp), allocatable :: data(:, :)
real(dp) :: E, d
call loadtxt("potential_table.txt", data)

rmin = 0
rmax = 12
N = 100000
au = 219474.62_dp
E0 = -0.625_dp
mu = 1836.12_dp / 2
print *, "     R         D          E         dE/dR"
do i = 1, size(data, 1)
    E = data(i, 2)
    print "(f10.4, f10.1, f12.7, f12.7)", data(i, 1), (E0-E)*au, E, data(i, 3)
end do

d = (rmax-rmin)/(N-1)
allocate(xe(N), U(N), lam(N))
xe = mesh_exp(rmin, rmax, 1._dp, N-1)  ! uniform mesh
U = spline3(data(:, 1), data(:, 2), xe)
!U = lagrange3interp(data(:, 1), data(:, 2), xe)
!U = hermite5interp(data(:, 1), data(:, 2), data(:, 3), .false., xe)
call savetxt("potential_grid.txt", reshape([xe, U], [size(xe), 2]))


print *, "solving"
lam = radial_schroedinger_fd(d, U, mu)
print *, "done"
print *
print *, "N =", N
print "('d =', f10.4)", d
print "('d =', f14.8)", d
print "('d =', f14.8)", xe(2) - xe(1)
print *, "Eigenvalues (n, E [a.u.], D [cm^-1]):"
do i = 1, 26
    E = lam(i)
    write(*,'(1x,i5,":",f18.8, f18.2)') i-1, E, (E0-E)*au
end do

contains

function radial_schroedinger_fd(d, U, mu) result(lam)
! Solves the radial Schroedinger equation with l=0 and potential 'U' given on a
! uniform mesh with element size 'd' using finite differences as described in
! [1]. Returns all eigenvalues in ascending order.
real(dp), intent(in) :: d  ! The element size of a uniform mesh
real(dp), intent(in) :: U(:) ! The potential given on a uniform mesh
real(dp), intent(in) :: mu ! The particle mass in Hartree atomic units (use
                           ! mu=1 for electrons)
real(dp) :: lam(size(U)) ! The eigenvalues
real(dp), allocatable :: work(:), c(:, :)
integer, allocatable :: iwork(:)
integer :: i, info, N
lam = 2*(1 + mu*d**2*U) ! Diagonal
N = size(U)
allocate(work(1 + 4*N), iwork(3 + 5*N))
call dstevd("N", N, lam, [(-1._dp, i=1,N-1)], c, N, work, size(work), iwork, &
    size(iwork), info) ! "N"="Compute eigenvalues only" -> "c" is not referenced
if (info /= 0) then
    print *, "INFO =", info
    call stop_error("info /= 0")
end if
lam = lam / (2*mu*d**2)
end function

end program
