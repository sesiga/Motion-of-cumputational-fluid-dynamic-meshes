!******************************************************************************
! constants utility
!******************************************************************************
!
! Description:
!   This module comprises the implementation of a set of data which defines a
!   number of different symbolic constants used as constants in numerical
!   calculations. 
!
!------------------------------------------------------------------------------
module constants

    use precision

    implicit none

    ! Public Parameters:

    ! Big constants. These are less than huge so that we have some room for
    ! comparisons and other arithmetic without overflowing.
    integer, parameter  :: I_ZERO = 0
    integer, parameter  :: BIG_INTEGER = huge(I_ZERO)/16
    real(wp), parameter :: R_ZERO = 0.0
    real(wp), parameter :: BIG_REAL = huge(R_ZERO)/16.0

    ! Frequently used universal constants.
    real(wp), parameter :: PI = acos(-1.0_wp)

    ! Frequently used character length constant.
    integer, parameter :: LENGTH = 480

    ! Tecplot element types
    character(LENGTH), parameter :: lineseg = "line"
    character(LENGTH), parameter :: triangle = "triangle"
    character(LENGTH), parameter :: quadrilateral = "quadrilateral"
    character(LENGTH), parameter :: tetra = "tetraeder"
    character(LENGTH), parameter :: hexa = "hexaeder"
  
end module constants
