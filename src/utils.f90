module utils

! Various utilities for general use in Fortran programs.

use types, only: dp
implicit none
private
public upcase, lowcase, white_char, blank, num_strings, get_string, &
    stop_error, savetxt, newunit, assert, str

interface str
    module procedure str_int, str_real, str_real_n
end interface

contains

function upcase(s) result(t)
! Returns string 's' in uppercase
character(*), intent(in) :: s
character(len(s)) :: t
integer :: i, diff
t = s; diff = ichar('A')-ichar('a')
do i = 1, len(t)
    if (ichar(t(i:i)) >= ichar('a') .and. ichar(t(i:i)) <= ichar('z')) then
        ! if lowercase, make uppercase
        t(i:i) = char(ichar(t(i:i)) + diff)
    end if
end do
end function

function lowcase(s) result(t)
! Returns string 's' in lowercase
character(*), intent(in) :: s
character(len(s)) :: t
integer :: i, diff
t = s; diff = ichar('A')-ichar('a')
do i = 1, len(t)
    if (ichar(t(i:i)) >= ichar('A') .and. ichar(t(i:i)) <= ichar('Z')) then
        ! if uppercase, make lowercase
        t(i:i) = char(ichar(t(i:i)) - diff)
    end if
end do
end function

logical function white_char(char) ! white character
! returns .true. if char is space (32) or tab (9), .false. otherwise
character, intent(in) :: char
if (iachar(char) == 32 .or. iachar(char) == 9) then
    white_char = .true.
else
    white_char = .false.
end if
end function

logical function blank(string)
! Returns true if string contains only white characters
character(*), intent(in) :: string
integer :: i
do i = 1, len(string)
    if (.not. white_char(string(i:i))) exit
end do
blank = (i>len(string))
end function

integer function num_strings(s) result(n)
! Returns number of substrings contained in input string 's' delimited
! by white space.
character(*), intent(in) :: s    ! input string
character(len(s)+2) :: t         ! temporary string to facilitate analysis
integer :: i
t = " " // s // " "
n = 0
do i = 1, len(t)-1
    if (white_char(t(i:i)) .and. .not. white_char(t(i+1:i+1))) n = n + 1
end do
end function

!--------------------------------------------------------------------------------------------------!

subroutine get_string(s,is,ss)
! Returns first substring ss in string s, delimited by white space, starting at
! index is in s. If ss is found, is is set to (index of last character of ss in
! s) + 1; else is is set to 0. If is is out of range on input, routine
! terminates with is = -1.
character(*), intent(in) :: s   ! input string
integer, intent(inout) :: is    ! on input: starting index for search for ss in
                                ! s on output: (index of last character of ss in
                                ! s) + 1
character(*), intent(out) :: ss ! first substring in s, starting from index is
character(len(s)+1) :: t        ! temporary string to facilitate search
integer i, i1, i2
logical prevwhite, curwhite
if (is <= 0 .or. is > len(s)) then
    ss = ""; is = -1; return
end if
t = s // " "
if (is == 1) then
    prevwhite = .true.
else
    prevwhite = white_char(t(is-1:is-1))
end if
i1 = 0; i2 = 0
do i = is, len(t)
    curwhite = white_char(t(i:i))
    if (prevwhite .and. .not. curwhite) i1 = i   ! beginning of substring
    if (i1>0 .and. curwhite) then                ! end of substring
        i2 = i-1; exit
    end if
    prevwhite=curwhite
end do
if (i2 > 0) then
    ss = t(i1:i2); is = i2+1
else
    ss = ""; is = 0
end if
end subroutine

integer function newunit(unit) result(n)
! Returns lowest i/o unit number not in use (to be used in older compilers).
!
! Starting at 10 to avoid lower numbers which are sometimes reserved.
! Note: largest valid unit number may be system-dependent.
!
! Arguments
! ---------
!
! If present, the new unit will be returned into it
integer, intent(out), optional :: unit
!
! Example
! -------
!
! integer :: u
! open(newunit(u), file="log.txt", status="old")
! read(u, *) a, b
! close(u)
!
! In new compilers, just use the "newunit" keyword argument:
!
! integer :: u
! open(newunit=u, file="log.txt", status="old")
! read(u, *) a, b
! close(u)

logical inuse
integer, parameter :: nmin=10   ! avoid lower numbers which are sometimes reserved
integer, parameter :: nmax=999  ! may be system-dependent
do n = nmin, nmax
    inquire(unit=n, opened=inuse)
    if (.not. inuse) then
        if (present(unit)) unit=n
        return
    end if
end do
call stop_error("newunit ERROR: available unit not found.")
end function

subroutine stop_error(msg)
! Aborts the program with nonzero exit code
!
! The statement "stop msg" will return 0 exit code when compiled using
! gfortran. stop_error() uses the statement "stop 1" which returns an exit code
! 1 and a print statement to print the message.
!
! Example
! -------
!
! call stop_error("Invalid argument")

character(len=*) :: msg ! Message to print on stdout
print *, msg
stop 1
end subroutine

subroutine savetxt(filename, d)
! Saves a 2D array into a textfile.
!
! Arguments
! ---------
!
character(len=*), intent(in) :: filename  ! File to save the array to
real(dp), intent(in) :: d(:, :)           ! The 2D array to save
!
! Example
! -------
!
! real(dp) :: data(3, 2)
! call savetxt("log.txt", data)

integer :: s, i
open(newunit(s), file=filename, status="replace")
do i = 1, size(d, 1)
    write(s, *) d(i, :)
end do
close(s)
end subroutine

subroutine assert(condition)
! If condition == .false., it aborts the program.
!
! Arguments
! ---------
!
logical, intent(in) :: condition
!
! Example
! -------
!
! call assert(a == 5)

if (.not. condition) call stop_error("Assert failed.")
end subroutine

pure integer function str_int_len(i) result(sz)
! Returns the length of the string representation of 'i'
integer, intent(in) :: i
integer, parameter :: MAX_STR = 100
character(MAX_STR) :: s
! If 's' is too short (MAX_STR too small), Fortan will abort with:
! "Fortran runtime error: End of record"
write(s, '(i0)') i
sz = len_trim(s)
end function

pure function str_int(i) result(s)
! Converts integer "i" to string
integer, intent(in) :: i
character(len=str_int_len(i)) :: s
write(s, '(i0)') i
end function

pure integer function str_real_len(r, fmt) result(sz)
! Returns the length of the string representation of 'i'
real(dp), intent(in) :: r
character(len=*), intent(in) :: fmt
integer, parameter :: MAX_STR = 100
character(MAX_STR) :: s
! If 's' is too short (MAX_STR too small), Fortan will abort with:
! "Fortran runtime error: End of record"
write(s, fmt) r
sz = len_trim(s)
end function

pure function str_real(r) result(s)
! Converts the real number "r" to string with 7 decimal digits.
real(dp), intent(in) :: r
character(len=*), parameter :: fmt="(f0.6)"
character(len=str_real_len(r, fmt)) :: s
write(s, fmt) r
end function

pure function str_real_n(r, n) result(s)
! Converts the real number "r" to string with 'n' decimal digits.
real(dp), intent(in) :: r
integer, intent(in) :: n
character(len=str_real_len(r, "(f0." // str_int(n) // ")")) :: s
write(s, "(f0." // str_int(n) // ")") r
end function

end module
