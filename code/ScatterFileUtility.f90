!******************************************************************************
! ScatterFile utility
!******************************************************************************
!
! Description:
!   This module comprises the implementation of the features (attributes and
!   operations) that characterize the objects of this class. A Tecplot data
!   file encapsulates the concept of a Tecplot data file, its data features are
!   defined by the Tecplot data file format specification.
!
!------------------------------------------------------------------------------
module ScatterFileUtility

    !------modules used--------------------------

    use precision, only : wp

    use constants

    use utilities
	
	use ScatterFileDataClass

    ! Resume file
    use ResumeFileUtility

	!-----end modules used-------------------------

    implicit none

    ! Private Parameters:
    character(*), parameter, private ::  MODname = 'ScatterFileUtility'

    ! Local Procedures:
    private :: Scatter_Open,   &
		Scatter_Close, &
		Scatter_InitializeFile
    
    ! Global (i.e. public) Declarations:
	public :: Scatter_LoadFileData, &
        Scatter_Finalize

	contains

!******************************************************************************
! loadScatterDataFile
!******************************************************************************
!
! Description:
!   Loads the scatter file and its name
!------------------------------------------------------------------------------
	
subroutine Scatter_LoadFileData(self, fileName)

    !* Subroutine arguments
    type(ScatterFileClass), intent(out) :: self
    character(*), intent(in)            :: fileName
    !*  End Subroutine arguments

!- End of header --------------------------------------------------------------
    
    !Initialize self components
    call Scatter_InitializeFile(self)
	
	self%filename = fileName !copy the file's name
	
	!Open the Scatter File
	call Scatter_Open(self%FileName,self%DataSet%id)
    
    !Loads the scatter data into de arrays
    call Scatter_ReadLoad(self)
    
    !Closes the scatter file
    call Scatter_Close(self%DataSet%id)

end subroutine Scatter_LoadFileData
!------------------------------------------------------------------------------

!****************************************************************
!scatter_open
!****************************************************************
!
!Description
!Opens the scatter file and assigns a ID 
!------------------------------------------------------------------------

subroutine Scatter_Open(path, ncid)

    !* Subroutine arguments:
    character(*), intent(in) :: path
    integer, intent(out)     :: ncid    
    !* End Subroutine arguments:

    ! Local scalars:
    integer :: ierr
    
!- End of header --------------------------------------------------------------

    ncid = freeunit()
    open(unit=ncid, file=path, status='old', action='read', iostat=ierr)

    if (ierr > 0) then
        call reporterror(MODname//': scatter_open reports: &
                         &Unsuccessful openning Scatter file.')
	end if

end subroutine Scatter_Open
!-------------------------------------------------------------------------------

!************************************************************
!Scatter_Close
!************************************************************
!
!Description
!Closes the scatter file

subroutine Scatter_Close(ncid)

   !subroutine arguments
   integer, intent(in) :: ncid
   !end subroutine arguments
   
   !end of header-----------------------------------------
   
   close(unit=ncid)
   
end subroutine Scatter_Close
!-----------------------------------------------------------

!*****************************************************
!Scatter_initialize
!*****************************************************
!
!Description
!Initialize self components

subroutine Scatter_InitializeFile(self)

    !* Subroutine arguments:
    type(ScatterFileClass), intent(out) :: self
    !* End Subroutine arguments:

!- End of header --------------------------------------------------------------
    
    self%filename = ""
    
    self%dataset%id = 0    
    self%dataset%nvars = 0
	
	self%ScatterData%no_of_points = 0
   
    self%ScatterData%points => null()
    self%ScatterData%values => null()
   
    self%ScatterData%global_id => null()
    
end subroutine Scatter_InitializeFile
!-----------------------------------------------------------------------------------

!************************************************************
!Scatter_ReadLoad
!************************************************************
!
!Description
!Reads scatter data, allocatates arrays and loads scatter data into arrays

subroutine Scatter_ReadLoad(self)

   !subroutine arguments
   type(ScatterFileClass), intent(inout) :: self
   !end subroutine arguments
   
   !Local scalars
   integer :: k, i
   integer :: ierr
   character(len=LENGTH) :: string
   
   !end of header----------------------------------------
   
   !loads number of variables
   rewind(self%DataSet%id)
   
   do k=1,1
	   read(unit=self%DataSet%id,fmt='(a)') string
	   string = removeblanks(string)
   end do
   
   k=0
   do i=1,len_trim(string)
	   if(string(i:i) == '"') k=k+1
   end do
   
   self%DataSet%nvars = (k-2)/4
   
   !loads number of points
   rewind(self%DataSet%id)
   
   read(unit=self%DataSet%id, fmt=*, iostat=ierr)
   do
       
	   read(unit=self%DataSet%id, fmt=*, iostat=ierr)
	   if(ierr>0)then
           call reporterror(MODname//': scatter_open reports: &
                         &Unsuccessful openning Scatter file.')
       elseif(ierr<0)then
           exit
       end if 
       
	   self%ScatterData%no_of_points = self%ScatterData%no_of_points + 1
       
   end do
   
   !loads data into arrays
   rewind(self%DataSet%id)
   
   allocate(self%ScatterData%points(self%ScatterData%no_of_points,self%DataSet%nvars))
   allocate(self%ScatterData%values(self%ScatterData%no_of_points,self%DataSet%nvars))
   allocate(self%ScatterData%global_id(self%ScatterData%no_of_points))
   
   read(unit=self%DataSet%id, fmt=*)       
   do i=1,self%ScatterData%no_of_points
       read(unit=self%DataSet%id, fmt=*) self%ScatterData%points(i,1:self%DataSet%nvars), &
                                  self%ScatterData%values(i,1:self%DataSet%nvars), &
                                  self%ScatterData%global_id(i)      
   end do
   
end subroutine Scatter_ReadLoad
!----------------------------------------------------------------------------------

!*****************************************************
!Scatter_Finalize
!*****************************************************
!
!Description
!Deallocates scatter data arrays

subroutine Scatter_Finalize(self)

   !subroutine arguments
   type(ScatterfileClass), intent(in) :: self
   !end subroutine arguments
   
   !end of header------------------------------------
   
   deallocate(self%ScatterData%points)
   deallocate(self%ScatterData%values)
   deallocate(self%ScatterData%global_id)
   
end subroutine Scatter_Finalize
!---------------------------------------------------------

end module ScatterFileUtility
