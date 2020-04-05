!******************************************************************************
! precision utility
!******************************************************************************
!
! Description:
!   This module comprises the implementation of a set of data which defines a
!   number of different symbolic constants used as kind constants for system-
!   independent precision specification. The kind constants are defined to
!   improve the portability between 32- and 64-bit platforms.
!
!------------------------------------------------------------------------------
module precision

implicit none

! Public Declarations:
public

! Public Parameters:

! Kinds for specified integer ranges.
integer, parameter :: i1 = selected_int_kind(2)  ! 1-byte
integer, parameter :: i2 = selected_int_kind(4)  ! 2-byte
integer, parameter :: i4 = selected_int_kind(9)  ! 4-byte

! Kinds for specified real precisions.
integer, parameter :: r4 = selected_real_kind(6,30)   ! 6 digits
integer, parameter :: r8 = selected_real_kind(12,30)  ! 12 digits

! System default real kinds.
integer, parameter :: rs = kind(0.)    ! real single precision
integer, parameter :: rd = kind(0.d0)  ! real double precision

! Kind for working real precision.
integer, parameter :: wp = rd

end module precision

