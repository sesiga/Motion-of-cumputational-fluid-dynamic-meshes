!******************************************************************************
! setup utility
!******************************************************************************
!
! Description:
!   This module comprises the implementation of a setup procedure.
!
!------------------------------------------------------------------------------
module setup

    use constants, only : LENGTH
    
    use precision, only : wp

    use utilities, only:               &
        reporterror,                   &
        freeunit
    
    use RadialBasisFunctions

    implicit none

    ! Local Parameters:
    character(*), parameter, private ::  MODname = 'setup'

    ! Global Procedures:
    public :: loadSetup

contains

!******************************************************************************
! Subprogram implementing a setup procedure
!******************************************************************************
!
!------------------------------------------------------------------------------
subroutine loadSetup(settings, deformationFileName, meshFileName, &
              & outputFileName, RBFfunction, RadiusCoeff)

    !*  Subroutine arguments
    character(*), intent(in)  :: settings
    character(*), intent(out) :: deformationFileName
    character(*), intent(out) :: meshFileName
    character(*), intent(out) :: outputFileName
    
    procedure(RBF), pointer :: RBFfunction
    real(wp), intent(out) :: RadiusCoeff
    !*  End Subroutine arguments

    !   Local scalars:
    integer           :: unit
    integer           :: ios
    character(LENGTH) :: folder
    character(LENGTH) :: input_deformationFileName
    character(LENGTH) :: input_meshFileName
    character(LENGTH) :: output_meshFileName
    character(LENGTH) :: input_FunctionName
    real(wp) :: input_RadiusCoefficient
    
    !   Namelist:
    namelist /INI/ folder
    namelist /INI/ input_deformationFileName
    namelist /INI/ input_meshFileName
    namelist /INI/ output_meshFileName
    namelist /INI/ input_FunctionName
    namelist /INI/ input_RadiusCoefficient

!- End of header --------------------------------------------------------------

    unit = freeunit()
    
    open(unit, file=settings, iostat=ios)
    if (ios /= 0) then
        call reporterror(MODname//': loadSetup reports: &
                         &"open" did not execute correctly.')
    end if
    read(unit, nml=INI, iostat=ios)
    if (ios /= 0) then
        call reporterror(MODname//': loadSetup reports: &
                         &"read" did not execute correctly.')
    end if
    close(unit, iostat=ios)
    if (ios /= 0) then
        call reporterror(MODname//': loadSetup reports: &
                         &"close" did not execute correctly.')
    end if

    deformationFileName = trim(folder)//trim(input_deformationFileName)
    meshFileName = trim(folder)//trim(input_meshFileName)
    outputFileName = trim(folder)//trim(output_meshFileName)
    RadiusCoeff = input_RadiusCoefficient
    
    !select function to use
    select case(input_FunctionName)
    case('Spline')
        RBFfunction => Spline
    case('Wendland C0')
        RBFfunction => WendlandC0
    case('Wendland C2')
        RBFfunction => WendlandC2
    case('Wendland C4')
        RBFfunction => WendlandC4
    case default
        call reporterror(MODname//': loadSetup reports: &
                         &"select case" did not execute correctly.')
    end select
    
end subroutine loadSetup

end module setup
