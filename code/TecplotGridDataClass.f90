!******************************************************************************
! TecplotGridDataClass
!******************************************************************************
!
! Description:
!   This module comprises the implementation of the features (attributes and
!   operations) that characterize the objects of this class. A data grid
!   encapsulates the concept of a geometric grid, that is, a set of points in a
!   point space (geometry) and, optionally, a set of connections between them
!   (topology). These are its attributes. A grid requires for its existence
!   points, but not connections, which is an optional attribute---i.e. we can
!   talk about a grid of scattered points.
!
!------------------------------------------------------------------------------
module TecplotGridDataClass

    use precision, only:               &
        wp

    use utilities, only:               &
        reporterror

    use constants, only:               &
        LENGTH

    implicit none

    ! Local Parameters:
    character(*), parameter, private :: MODname = 'TecplotGridDataClass' 

    !   Global (i.e. public) Declarations:
    public

    type TecplotZoneHeader
        character(LENGTH) :: header
        character(LENGTH) :: datapacking
        character(LENGTH) :: zonetype
		integer, pointer  :: IJK(:)
        integer           :: no_of_points
        integer           :: no_of_elements
    end type TecplotZoneHeader

    type TecplotDataField
		type(TecplotZoneHeader) :: fieldheader
        character(LENGTH)       :: elementType
        integer, pointer        :: connections(:,:) 
        integer, pointer        :: global_id(:)	        ! ipoint en el array points(:,:) de la malla volumetrica
    end type TecplotDataField

    type TecplotDataGrid
        real(wp), pointer               :: points(:,:) 
        real(wp), pointer               :: values(:,:) 
        type(TecplotDataField), pointer :: fields(:)  
    end type TecplotDataGrid

    interface finalize
        module procedure finalize_tecplotgrid,          &
                         finalize_tecplotfield
    end interface finalize

contains

!******************************************************************************
! finalize_tecplotgrid
!******************************************************************************
!
! Description:
!    finalize grid, freeing all resources.
!------------------------------------------------------------------------------
subroutine finalize_tecplotgrid(grid)

    !* Subroutine arguments:
    type(TecplotDataGrid), intent(inout) :: grid
    !* End Subroutine arguments:
    
	! local scalars:
	integer :: i

!- End of header --------------------------------------------------------------


    if (associated(grid%points)) deallocate(grid%points)

    if (associated(grid%values)) deallocate(grid%values)
            
	do i=1,size(grid%fields)
        call finalize(grid%fields(i))           
	end do

end subroutine finalize_tecplotgrid

!******************************************************************************
! finalize_tecplotfield
!******************************************************************************
!
! Description:
!    finalize field, freeing all resources.
!------------------------------------------------------------------------------
subroutine finalize_tecplotfield(field)

    !* Subroutine arguments:
    type(TecplotDataField), intent(inout) :: field
    !* End Subroutine arguments:
    
!- End of header --------------------------------------------------------------

    if (associated(field%connections)) deallocate(field%connections)

    if (associated(field%global_id)) deallocate(field%global_id)

end subroutine finalize_tecplotfield

end module TecplotGridDataClass
