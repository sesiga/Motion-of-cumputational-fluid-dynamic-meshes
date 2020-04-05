!******************************************************************************
! TecplotFileDataClass
!******************************************************************************
!
! Description:
!   This module comprises the implementation of the features (attributes and
!   operations) that characterize the objects of this class. A data field
!   encapsulates the concept of a mathematical field, that is, a mapping from
!   some domain (grid) to some data space (values) and, optionally, a name
!   (label) for the field. These are its attributes. A field requires for its
!   existence both a grid and values, but not a label, which is an optional
!   attribute.
!
!------------------------------------------------------------------------------
module TecplotFileDataClass

    use precision, only:               &
        wp

    use constants, only:               &
        LENGTH

    use utilities, only:               &
        reporterror

    ! Self consists of a DataGrid (component).
    use TecplotGridDataClass

    implicit none

    ! Local Parameters:
    character(*), parameter, private :: MODname = 'TecplotFileDataClass' 

    ! Global (i.e. public) Declarations:
    public

    ! Global Type Definitions:
    type TecplotDataset
        integer               :: id           ! id del dataset (unit file)
        character(len=LENGTH) :: name         ! nombre del dataset (title) 
        integer               :: nvars        ! numero de variables del dataset
        integer               :: dimspace
    end type TecplotDataset

    type TecplotVariables
        character(len=LENGTH) :: name         ! variable name
        integer               :: id           ! variable id: 1-3:dimension 4-6:desplazamiento 0:global_id -1:otro
    end type TecplotVariables

    type TecplotGriddata
        integer :: no_of_points                       ! # of points
        integer :: no_of_elements                     ! # of elements
        integer :: no_of_tetraeders                   ! # of tetrahedrons
        integer :: no_of_hexaeders                    ! # of hexahedrons
        integer :: no_of_triangles                    ! # of triangles
        integer :: no_of_quadrilaterals               ! # of quadrilaterlas
        integer :: no_of_linesegs                     ! # of linsegs
        integer :: points_per_tetraeder               ! # of nodes per tetraeder
        integer :: points_per_hexaeder                ! # of nodes per hexahedron
        integer :: points_per_triangle                ! # of nodes per triangle
        integer :: points_per_quadrilateral           ! # of nodes per quadrilateral
        integer :: points_per_lineseg                 ! # of nodes per lineseg
    end type TecplotGriddata

    type  TecplotFileClass
        character(LENGTH)               :: filename
        type(TecplotDataset)            :: dataset
        type(TecplotVariables), pointer :: variables(:)
        type(TecplotGriddata)           :: griddata
        type(TecplotDataGrid)           :: grid
    end type TecplotFileClass

    interface finalize
        module procedure finalize_tecplotfile
    end interface finalize

contains

!******************************************************************************
! finalize_tecplotfile
!******************************************************************************
!
! Description:
!    finalize tecplot file, freeing all resources.
!------------------------------------------------------------------------------
subroutine finalize_tecplotfile(self)

    !* Subroutine arguments:
    type(TecplotFileClass), intent(inout) :: self
    !* End Subroutine arguments:
    
!- End of header --------------------------------------------------------------

    if (associated(self%variables)) deallocate(self%variables)

    call finalize(self%grid)

end subroutine finalize_tecplotfile

!******************************************************************************
! get_dataset
!******************************************************************************
!
! Description:
!    finalize tecplot file, freeing all resources.
!------------------------------------------------------------------------------
function get_dataset(self) result(dataset)

    !* Subroutine arguments:
    type(TecplotFileClass), intent(in), target :: self

    type(TecplotDataset), pointer              :: dataset
    !* End Subroutine arguments:
    
!- End of header --------------------------------------------------------------

    dataset => self%dataset

end function get_dataset

end module TecplotFileDataClass
