module constants

! Contains the mathematical constant "pi".

use types, only: dp
implicit none
private
public pi

! Contain more digits than double precision, so that
! it is rounded correctly:
real(dp), parameter :: pi   = 3.1415926535897932384626433832795_dp

end module
