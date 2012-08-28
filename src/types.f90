module types

! This module defines the 'dp' double precision type.

implicit none
private
public dp

integer, parameter :: dp=kind(0.d0)          ! double precision

end module
