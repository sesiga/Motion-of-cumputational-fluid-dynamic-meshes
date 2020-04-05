!******************************************************************************
! TecplotFile utility
!******************************************************************************
!
! Description:
!   This module comprises the implementation of the features (attributes and
!   operations) that characterize the objects of this class. A Tecplot data
!   file encapsulates the concept of a Tecplot data file, its data features are
!   defined by the Tecplot data file format specification.
!
!------------------------------------------------------------------------------
module TecplotFileUtility

    use precision, only:               &
        wp

    use constants

    use utilities

    ! File data
    use TecplotFileDataClass

    ! Grid data
    use TecplotGridDataClass

    ! Resume file
    use ResumeFileUtility
    
    !quality data class
    use QualityDataClass

    implicit none

    ! Private Parameters:
    character(*), parameter, private ::  MODname = 'TecplotDataFileClass'

    ! Local Procedures:
    private ::                              &
        tp_open,                            &
        tp_inquire_dataset,                 &
        tp_inquire_title,                   &
        tp_inquire_nvars,                   &
        tp_inquire_variables,               &
        tp_inquire_dimspace,                &
        tp_inquire_grid,                    &
        tp_inquire_nofields,                &
        tp_inquire_fields,                  &
        tp_inquire_header,                  &
        tp_inquire_connections,             &
        tp_inquire_datavalues,              &
        process_datavalues,                 &
        tp_inquire_griddata,                &
        record_parameter,                   &
        generate_connections,               &
        IndexPoint,                         &
        read_variables,                     &
        readBLOCK,                          &
        tp_close,                           &
        Tecplotfile_initialize,             &
        initialize_griddata,                &
        initialize_grid,                    &
        initialize_variables,               &
        initialize_fields,                  &
        loadTecplotDataFile,                &
        saveTecplotDataFile

    ! Procedure overloading:
    interface initialize
        module procedure Tecplotfile_initialize,    &
                         initialize_griddata,       &
                         initialize_grid,           &
                         initialize_variables,      &
                         initialize_fields      
    end interface initialize

    ! Global (i.e. public) Declarations:
    public

    ! Procedure overloading:
    interface loadDataFile
        module procedure loadTecplotDataFile
    end interface loadDataFile

    ! Procedure overloading:
    interface saveDataFile
        module procedure saveTecplotDataFile
    end interface saveDataFile

contains

!******************************************************************************
! loadTecplotDataFile
!******************************************************************************
!
! Description:
!   Carga el fichero Tecplot
!------------------------------------------------------------------------------
subroutine loadTecplotDataFile(self, fileName)

    !* Subroutine arguments
    type(TecplotFileClass), intent(out) :: self
    character(*), intent(in)            :: fileName
    !*  End Subroutine arguments

!- End of header --------------------------------------------------------------

    ! print*, '    inicializa self'
    ! inicializa todas las componentes de self
    call initialize(self)

    ! copia el nombre del fichero
    self%filename = filename  

    ! print*, '    tp_open'
    ! abre el fichero de datos y le asigna un id
    call tp_open(path=self%filename, ncid=self%dataset%id)

    ! print*, '    tp_inquire_dataset'
    ! Extrae dataset (excepto dimspace)
    call tp_inquire_dataset(ncid=self%dataset%id, dataset=self%dataset) 

    ! allocata e inicializa variables
    call initialize(self%variables, self%dataset%nvars)

    ! print*, '    tp_inquire_variables'
    ! Extrae el vector variables
    call tp_inquire_variables(ncid=self%dataset%id, variables=self%variables) 

    ! print*, '    tp_inquire_dimspace'
    ! Extrae dimspace
    call tp_inquire_dimspace(variables=self%variables, dimspace=self%dataset%dimspace) 

    ! print*, '    tp_inquire_grid'
    ! Extrae grid
    call tp_inquire_grid(dataset=self%dataset, variables=self%variables, grid=self%grid) 

    ! print*, '    tp_inquire_griddata'
	! Carga las componentes de griddata
	call tp_inquire_griddata(grid=self%grid, griddata=self%griddata)
    
    ! print*, '    tf_close'
    ! Close the file, freeing all resources.
    call tp_close(self=self)

    ! escribe resume file
    call write_resumefile(self=self)

end subroutine loadTecplotDataFile
!------------------------------------------------------------------------------

!******************************************************************************
! tp_open
!******************************************************************************
!
! Description:
!   abre el fichero de datos y le asigna un id
!------------------------------------------------------------------------------
subroutine tp_open(path, ncid)

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
        call reporterror(MODname//': tp_open reports: &
                         &Unsuccessful openning Tecplot file.')
    end if

end subroutine tp_open

!******************************************************************************
! tp_inquire_dataset
!******************************************************************************
!
! Description:
!   extrae dataset
!------------------------------------------------------------------------------
subroutine tp_inquire_dataset(ncid, dataset)

    !* Subroutine arguments:
    integer, intent(in)                 :: ncid    
    type(TecplotDataset), intent(inout) :: dataset
    !* End Subroutine arguments:
    
!- End of header --------------------------------------------------------------

    ! TITLE
    call tp_inquire_title(ncid, title=dataset%name)

    ! nvars = # VARIABLES
    call tp_inquire_nvars(ncid, nvars=dataset%nvars)

end subroutine tp_inquire_dataset

!******************************************************************************
! tp_inquire_title
!******************************************************************************
!
! Description:
!   extrae dataset%title (titulo del fichero)
!------------------------------------------------------------------------------
subroutine tp_inquire_title(ncid, title)

    !* Subroutine arguments:
    integer, intent(in)       :: ncid    
    character(*), intent(out) :: title
    !* End Subroutine arguments:

    ! Local scalars:
    integer :: i1
    integer :: i2
    logical :: istitle
    
    character(LENGTH) :: line
    character(LENGTH) :: string
    
!- End of header --------------------------------------------------------------

	rewind(ncid)
	title = ""

    ! busca TITLE
    istitle = .false.
    do
        read(unit=ncid, fmt='(a)') line
		string = uppercase(line)
		string = removeblanks(string)
        if (index(string,"TITLE=") == 1) then
            istitle = .true.
            exit
        end if

		if (index(string,"ZONE") == 1) exit
    end do

    if (istitle) then
        i1 = index(line,'"')
        i2 = i1 + index(line(i1+1:),'"')

        title = line(i1+1:i2-1)
    end if

end subroutine tp_inquire_title

!******************************************************************************
! tp_inquire_nvars
!******************************************************************************
!
! Description:
!   extrae dataset%nvars (numero de variables del fichero)
!------------------------------------------------------------------------------
subroutine tp_inquire_nvars(ncid, nvars)

    !* Subroutine arguments:
    integer, intent(in)  :: ncid    
    integer, intent(out) :: nvars
    !* End Subroutine arguments:

    ! Local scalars:
    integer :: i
    integer :: linelength
    
    character(LENGTH) :: line
    character(LENGTH) :: string
    
!- End of header --------------------------------------------------------------

	rewind(ncid)

	! busca VARIABLES
    do
        read(unit=ncid, fmt='(a)') line
		string = uppercase(line)
		string = removeblanks(string)
        if (index(string,"VARIABLES=") == 1) exit

		if (index(string,"ZONE") == 1) then
		    call reporterror(MODname//': tf90_inquire_nvars: &
		                     &VARIABLES not found.')
		end if
    end do

	! cuenta el numero de variables
	nvars = 0
	linelength = len_trim(line)
    do i=1,linelength
		if (line(i:i) == '"') nvars = nvars + 1
    end do
	nvars = nvars/2

end subroutine tp_inquire_nvars

!******************************************************************************
! tp_inquire_variables
!******************************************************************************
!
! Description:
!   extrae el vector variables
!------------------------------------------------------------------------------
subroutine tp_inquire_variables(ncid, variables)

    !* Subroutine arguments:
    integer, intent(in)             :: ncid    
    type(TecplotVariables), pointer :: variables(:)
    !* End Subroutine arguments:

    ! Local scalars:
    integer :: ivar
    integer :: i1
    integer :: i2

    character(LENGTH) :: line
    character(LENGTH) :: string    
    
!- End of header --------------------------------------------------------------

	rewind(ncid)

	! busca VARIABLES
    do
        read(unit=ncid, fmt='(a)') line
		string = uppercase(line)
		string = removeblanks(string)
        if (index(string,"VARIABLES=") == 1) exit

		if (index(string,"ZONE") == 1) then
		    call reporterror(MODname//': tf_inquire_variables: &
		                     &VARIABLES not found.')
		end if
    end do

	! carga VARIABLES
	ivar = 0	
    do
		i2 = 0
		do
			i1 = index(line(i2+1:), '"')
			if (i1 == 0) exit
			i1 = i1 + i2
			i2 = index(line(i1+1:), '"')
			ivar = ivar + 1
			i2 = i1 + i2
			variables(ivar)%name = line(i1+1:i2-1)

            variables(ivar)%id = -1
            select case (variables(ivar)%name)
            case ("x", "X")
                variables(ivar)%id = 1
            case ("y", "Y")
                variables(ivar)%id = 2
            case ("z", "Z")
                variables(ivar)%id = 3
            case ("dx", "dX")
                variables(ivar)%id = 4
            case ("dy", "dY")
                variables(ivar)%id = 5
            case ("dz", "dZ")
                variables(ivar)%id = 6
            case ("global_id")
                variables(ivar)%id = 0
            end select    
		end do

	    read(unit=ncid, fmt='(a)') line
		string = uppercase(line)
		string = adjustl(string)

		if (index(string,"ZONE") == 1) exit
    end do

end subroutine tp_inquire_variables

!******************************************************************************
! tp_inquire_dimspace
!******************************************************************************
!
! Description:
!   extrae el la dimension del problema
!------------------------------------------------------------------------------
subroutine tp_inquire_dimspace(variables, dimspace)

    !* Subroutine arguments:
    type(TecplotVariables), pointer :: variables(:)
    integer, intent(out)            :: dimspace
    !* End Subroutine arguments:

    ! Local scalars:
    integer :: maxid
    
!- End of header --------------------------------------------------------------

    maxid = maxval(variables(:)%id)

    dimspace = 0
    select case (maxid)
    case (1, 4)
        dimspace = 1
    case (2, 5)
        dimspace = 2
    case (3, 6)
        dimspace = 3
    end select

    if (dimspace == 0) then
        call reporterror(MODname//': tf_inquire_dimspace: &
                            &Error in dimspace.')
    end if

end subroutine tp_inquire_dimspace

!******************************************************************************
! tp_inquire_grid
!******************************************************************************
!
! Description:
!   extrae la componente grid
!------------------------------------------------------------------------------
subroutine tp_inquire_grid(dataset, variables, grid)

    !* Subroutine arguments:
    type(TecplotDataset), intent(in)   :: dataset
    type(TecplotVariables), intent(in) :: variables(:)
    type(TecplotDataGrid), intent(out) :: grid
    !* End Subroutine arguments:

    ! Local scalars:
    integer :: no_of_fields    
    integer :: no_of_datavalues

    ! Local arrays:
    real(wp), pointer :: datavalues(:,:)
    
!- End of header --------------------------------------------------------------

	! Extrae numero de fields (zones)
	call tp_inquire_nofields(dataset%id, no_of_fields)

	! Alocata en inicializa fields
	call initialize(no_of_fields, fields=grid%fields)

	! Carga componentes de fields: fieldheader + elementType 
	! Carga componentes de fields: connections + global_id
	call tp_inquire_fields(dataset=dataset, fields=grid%fields,         &
                           no_of_datavalues=no_of_datavalues)
    
    ! Carga el vector completo de datavalues 
	call tp_inquire_datavalues(dataset=dataset, fields=grid%fields,     &
                              no_of_datavalues=no_of_datavalues,        &
                              datavalues=datavalues)
    
    ! Procesa datavalues y extrae: points + values + global_id
    call process_datavalues(dataset, variables, datavalues, grid)

	! finaliza datavalues
	deallocate(datavalues)

end subroutine tp_inquire_grid

!******************************************************************************
! tp_inquire_nofields
!******************************************************************************
!
! Description:
!	Extrae el numero de fields (zones de tecplot)  
!------------------------------------------------------------------------------
subroutine tp_inquire_nofields(ncid, no_of_fields)

    !* Subroutine arguments:
    integer, intent(in)  :: ncid
	integer, intent(out) :: no_of_fields
    !* End Subroutine arguments:
    
	! Local scalars:
	integer :: iend

	! Local characters:
	character(LENGTH) :: line

!- End of header --------------------------------------------------------------

	rewind(ncid)

	no_of_fields = 0
	do
		read(ncid,'(a)',iostat=iend) line
		if (iend < 0) exit
		line = adjustl(line)
		if (index(line,'ZONE ') == 1) no_of_fields = no_of_fields + 1
	end do    

end subroutine tp_inquire_nofields

!******************************************************************************
! tp_inquire_fields
!******************************************************************************
!
! Description:
!	- Carga la informacion de fields
!   - Devuelve no_of_datavalues: no. total de points con repeticiones
!------------------------------------------------------------------------------
subroutine tp_inquire_fields(dataset, fields, no_of_datavalues)

    !* Subroutine arguments:
    type(TecplotDataset), intent(in)    :: dataset
    type(TecplotDataField), intent(out) :: fields(:)
    integer, intent(out)                :: no_of_datavalues
    !* End Subroutine arguments:
    
	! Local scalars:
	integer :: iend
    integer :: ifield

	! Local characters:
	character(LENGTH) :: line
	character(LENGTH) :: string
	character(LENGTH) :: header

!- End of header --------------------------------------------------------------

	rewind(dataset%id)

	! Primero extrae la informacion del header de cada field
	no_of_datavalues = 0
    ifield = 0
	do 
		read(unit=dataset%id,fmt='(a)',iostat=iend) line
		if (iend < 0) exit

		string = adjustl(line)

		if (index(string,'ZONE ') == 1) then

            ifield = ifield + 1

			header = " "			
			header = trim(header)//trim(line)
			do
				read(unit=dataset%id,fmt='(a)',iostat=iend) line
				string = adjustl(line)
				if (isnumbermark(string(1:1))) then
					backspace(unit=dataset%id)
					exit
				else
					header = trim(header)//' '//trim(line)
				end if					
			end do

			! copia todo el header de zone en una unica linea
			! carga las componentes no array de fileds(i)
			call tp_inquire_header(header, field=fields(ifield))

			! carga  global_id
			 call tp_inquire_global_id(dataset=dataset,                              &
			 		  	              fieldheader=fields(ifield)%fieldheader,        &
			 				          global_id=fields(ifield)%global_id)

			! carga  connections
			 call tp_inquire_connections(dataset=dataset,                               &
			 		  	                 fieldheader=fields(ifield)%fieldheader,        &
			 					         elementType=fields(ifield)%elementType,	    &
			 				             connections=fields(ifield)%connections)
            
			! recuenta para poder alocatar el # total de points
			no_of_datavalues = no_of_datavalues + fields(ifield)%fieldheader%no_of_points

		end if

        if (ifield == size(fields)) exit

	end do    

end subroutine tp_inquire_fields

!******************************************************************************
! tp_inquire_header
!******************************************************************************
!
! Description:
!   Extrae toda la informacion del header del field    
!------------------------------------------------------------------------------
subroutine tp_inquire_header(header, field)

    !* Subroutine arguments:
    character(*), intent(in)            :: header
	type(TecplotDataField), intent(out) :: field
    !* End Subroutine arguments:
    
	! Local scalars:
	integer :: i
	integer :: i1
	integer :: i2
    integer :: dimspace

	! Local characters:
	character(LENGTH) :: scratch
	character(LENGTH) :: line
	character(LENGTH) :: subscratch

    ! Local arrays:
    integer :: ijk(3)

!- End of header --------------------------------------------------------------

	! Copia la linea completa
	field%fieldheader%header = trim(header)	

	! Extrae la informacion
	line = uppercase(header)

    i1 = index(line, '"')
    i2 = 0
    if (i1 /= 0) then
        i2 = i1 + index(line(i1+1:), '"') 
    end if

    scratch = line(i2+1:)
	
	! quita los espacios en blanco antes y despues del =
	scratch = remove_selectedblanks(scratch, selectedchar="=")

	! cambia las comas (si las hay) por espacios en blanco
	scratch = exchange_chars(scratch, removedchar=",", newchar=" ")
    
	! Identify the zone format (if it exists, else POINT).
    call record_parameter(line=scratch, string='DATAPACKING', default='POINT',  &
                         actual=field%fieldheader%datapacking)

	!   Identify the zone element type (if it exists, else ORDERED).
    call record_parameter(line=scratch, string='ZONETYPE', default='ORDERED',   &
                          actual=field%fieldheader%zonetype)

    select case (field%fieldheader%zonetype)

    case ('ORDERED')

        ijk = 0
        dimspace = 1
        call record_parameter(line=scratch, string='I', default='1', actual=subscratch)
        read(subscratch, *) ijk(1)
        field%elementType = lineseg

        call record_parameter(line=scratch, string='J', default='1', actual=subscratch)
        read(subscratch, *) ijk(2)
        if (ijk(2) > 1) then
            field%elementType = quadrilateral
            dimspace = dimspace + 1
        end if 

        call record_parameter(line=scratch, string='K', default='1', actual=subscratch)
        read(subscratch, *) ijk(3)
        if (ijk(3) > 1)  then
            field%elementType = hexa
            dimspace = dimspace + 1
        end if 

        allocate(field%fieldheader%ijk(dimspace))
        field%fieldheader%IJK = ijk(1:dimspace)

        field%fieldheader%no_of_points = product(field%fieldheader%IJK)

        field%fieldheader%no_of_elements = 1
		do i=1,dimspace
	        field%fieldheader%no_of_elements =  &
                field%fieldheader%no_of_elements*(field%fieldheader%IJK(i)-1)
		end do
        
    case default
    
        call record_parameter(line=scratch, string='N', default='0', actual=subscratch)
        read(subscratch, *) field%fieldheader%no_of_points

        call record_parameter(line=scratch, string='E', default='0', actual=subscratch)
        read(subscratch, *) field%fieldheader%no_of_elements

        select case (field%fieldheader%zonetype)
        case ('FELINESEG')
            field%elementType = lineseg

        case ('FETRIANGLE')
            field%elementType = triangle
        
        case ('FEQUADRILATERAL')
            field%elementType = quadrilateral

        case ('FETETRAHEDRON')
            field%elementType = tetra

        case ('FEBRICK')
            field%elementType = hexa

        case default
            call reporterror(MODname//': tp_inquire_header reports: &
                            &Not valid element type.')
        end select

    end select

end subroutine tp_inquire_header

!******************************************************************************
! tp_inquire_global_id
!******************************************************************************
!
! Description:
!   Extrae global_id
!------------------------------------------------------------------------------
subroutine tp_inquire_global_id(dataset, fieldheader, global_id)

    !* Subroutine arguments:
	type(TecplotDataset), intent(in)    :: dataset
	type(TecplotZoneHeader), intent(in) :: fieldheader
	integer, pointer      			    :: global_id(:)
    !* End Subroutine arguments:
    
	! Local scalars:
	integer :: i
	integer :: ierr
    
	! local allocatable:
	real(wp), allocatable :: datos(:)

!- End of header --------------------------------------------------------------

	! Reserva memoria para global_id
	allocate(global_id(fieldheader%no_of_points), stat=ierr)
    if (ierr /= 0) then
        call reporterror(MODname//': tp_inquire_global_id reports: &
                         &Unsuccessful allocation of global_id.')
    end if

    allocate(datos(2*dataset%dimspace))

    do i=1,fieldheader%no_of_points
        read(dataset%id,*) datos(:), global_id(i)
    end do
		 
    deallocate(datos)

end subroutine tp_inquire_global_id

!******************************************************************************
! tp_inquire_connections
!******************************************************************************
!
! Description:
!   Extrae o genera las conectividades
!------------------------------------------------------------------------------
subroutine tp_inquire_connections(dataset, fieldheader, elementType, connections)

    !* Subroutine arguments:
	type(TecplotDataset), intent(in)    :: dataset
	type(TecplotZoneHeader), intent(in) :: fieldheader
	character(*), intent(in)            :: elementType
	integer, pointer      			    :: connections(:,:)
    !* End Subroutine arguments:
    
	! Local scalars:
	integer :: i
	integer :: ierr
	integer :: nvertex
	logical :: isconnectline

	! local characters:
	character(LENGTH) :: line

!- End of header --------------------------------------------------------------

	select case (elementType)
	case (lineseg)
		nvertex = 2
	case (triangle)
		nvertex = 3
	case (tetra, quadrilateral)
		nvertex = 4
	case (hexa)		
		nvertex = 8
	end select

	! Reserva memoria para connections
	allocate(connections(fieldheader%no_of_elements,nvertex), stat=ierr)
    if (ierr /= 0) then
        call reporterror(MODname//': tp_inquire_connections reports: &
                         &Unsuccessful allocation of connections.')
    end if
	
	select case (fieldheader%zonetype)
	case ('ORDERED') 

		call generate_connections(ijk=fieldheader%ijk, connections=connections)

	case default

		isconnectline = .false.
		do
			read(dataset%id,'(a)') line
		
			do i=1,len_trim(line)
				if (isdigit(line(i:i))) then
					isconnectline = .true.
					backspace(unit=dataset%id)
					exit
				end if
			end do
		
		    if (isconnectline) exit
		end do		

		do i=1,fieldheader%no_of_elements
			read(dataset%id,*) connections(i,:)            
		end do
	end select

end subroutine tp_inquire_connections

!******************************************************************************
! tp_inquire_datavalues
!******************************************************************************
!
! Description:
!	Extrae el array completo de datavalues
!------------------------------------------------------------------------------
subroutine tp_inquire_datavalues(dataset, fields, no_of_datavalues, datavalues)

    !* Subroutine arguments:
    type(TecplotDataset), intent(in)   :: dataset
    type(TecplotDataField), intent(in) :: fields(:)
    integer, intent(in)                :: no_of_datavalues
    real(wp), pointer                  :: datavalues(:,:)
    !* End Subroutine arguments:
    
	! Local scalars:
	integer :: ifield
    integer :: ierr
    integer :: i1
    integer :: i2
    integer :: iend

	! Local characters:
	character(LENGTH) :: line
	character(LENGTH) :: string

!- End of header --------------------------------------------------------------

	! Reserva memoria para cargar todos los datos de variables 
	allocate(datavalues(no_of_datavalues, 2*dataset%dimspace), stat=ierr)
    if (ierr /= 0) then
        call reporterror(MODname//': tp_inquire_points: &
                         &Unsuccessful allocation of rpoints.')
    end if
    
	rewind(dataset%id)

	ifield = 0
	i1 = 0
	do
		read(unit=dataset%id,fmt='(a)',iostat=iend) line
		if (iend < 0) exit

		string = adjustl(line)

		if (index(string,'ZONE ') == 1) then
			ifield = ifield + 1

			do
				read(unit=dataset%id,fmt='(a)',iostat=iend) line
				string = adjustl(line)
				if (isnumbermark(string(1:1))) then
					backspace(unit=dataset%id)
					exit
				end if					
			end do

			i2 = i1 + fields(ifield)%fieldheader%no_of_points
			call read_variables(dataset, fields(ifield)%fieldheader,        &
                                datavalues=datavalues(i1+1:i2,:))
			i1 = i2
		end if

		if (ifield == size(fields)) exit

	end do

end subroutine tp_inquire_datavalues

!******************************************************************************
! record_parameter
!******************************************************************************
!
! Description:
!   Busca en "line" el valor de "string" y lo asigna a "actual".
!   Si no tiene valor, le asigna el valor por defecto "default"
!------------------------------------------------------------------------------
subroutine record_parameter(line, string, default, actual)

	!*  Subroutine arguments
    character(*), intent(in)  :: line
    character(*), intent(in)  :: string
    character(*), intent(in)  :: default
    character(*), intent(out) :: actual
	!*  End Subroutine arguments

	!   Local scalars:
	integer :: i
	integer :: i1
	integer :: i2

    character(lENGTH) :: dummy

!- End of header --------------------------------------------------------------
    
    dummy = " "//trim(string)//'='
    i1 = index(line,trim(dummy))

    actual = ''
    if (i1 /= 0) then

	    i2 = index(line(i1+1:),"=")
		i2 = i1 + i2 + 1

		i = i2
		do
			if (.not.isalpha(line(i:i))) exit
			i = i + 1
		end do

		read(line(i2:i), '(a)') actual
    else
        actual = default
    end if
    
end subroutine record_parameter

!******************************************************************************
! generate_connections
!******************************************************************************
!
! Description:
!	genra las connections a partir de un ijk
!------------------------------------------------------------------------------
subroutine generate_connections(ijk, connections)

	!*  Subroutine arguments
	integer, intent(in)  :: ijk(:)
	integer, intent(out) :: connections(:,:)
	!*  End Subroutine arguments

	!   Local scalars:
	integer :: i
	integer :: j
	integer :: k
	integer :: ielement
	integer :: ndims

	! local arrays:
	integer :: ijkP(3)

!- End of header --------------------------------------------------------------

	ndims = count(ijk /= 1)

	select case (ndims)

	case (3)

		ielement = 0
		do k=1,ijk(3)-1
			do j=1,ijk(2)-1
				do i=1,ijk(1)-1
					ielement = ielement + 1

					ijkP = [i,j,k]
					connections(ielement,1) = indexPoint(ijkP, ijk)

					ijkP = [i+1,j,k]
					connections(ielement,2) = indexPoint(ijkP, ijk)

					ijkP = [i+1,j+1,k]
					connections(ielement,3) = indexPoint(ijkP, ijk)

					ijkP = [i,j+1,k]
					connections(ielement,4) = indexPoint(ijkP, ijk)

					ijkP = [i,j,k+1]
					connections(ielement,5) = indexPoint(ijkP, ijk)

					ijkP = [i+1,j,k+1]
					connections(ielement,6) = indexPoint(ijkP, ijk)

					ijkP = [i+1,j+1,k+1]
					connections(ielement,7) = indexPoint(ijkP, ijk)

					ijkP = [i,j+1,k+1]
					connections(ielement,8) = indexPoint(ijkP, ijk)
				end do
			end do
		end do

	case (2)

		ielement = 0
		do j=1,ijk(2)-1
			do i=1,ijk(1)-1
				ielement = ielement + 1

				ijkP = [i,j,1]
				connections(ielement,1) = indexPoint(ijkP, ijk)

				ijkP = [i+1,j,1]
				connections(ielement,2) = indexPoint(ijkP, ijk)

				ijkP = [i+1,j+1,1]
				connections(ielement,3) = indexPoint(ijkP, ijk)

				ijkP = [i,j+1,1]
				connections(ielement,4) = indexPoint(ijkP, ijk)
			end do
		end do

	case (1)

		ielement = 0
		do i=1,ijk(1)-1
			ielement = ielement + 1

			ijkP = [i,1,1]
			connections(ielement,1) = indexPoint(ijkP, ijk)

			ijkP = [i+1,1,1]
			connections(ielement,2) = indexPoint(ijkP, ijk)
		end do

	end select

end subroutine generate_connections

!******************************************************************************
! IndexPoint
!******************************************************************************
!
! Description:
!   Given (i,j,k) return the index of that point in points(:,:)
!------------------------------------------------------------------------------
function IndexPoint(ijkP, ijk)

	!*  Function arguments
    integer, intent(in) :: ijkP(:)
    integer, intent(in) :: ijk(:)
    integer  :: IndexPoint
	!*  Enf Function arguments

!- End of header --------------------------------------------------------------

    IndexPoint = ijkP(1) + (ijkP(2)-1)*ijk(1) + &
                (ijkP(3)-1)*ijk(1)*ijk(2)

end function IndexPoint

!******************************************************************************
! read_variables
!******************************************************************************
!
! Description:
!   - if (present(datavalues)) read datavalues	
!   - if (.not.present(datavalues)) salta las lineas con los datos de las variables
!------------------------------------------------------------------------------
subroutine read_variables(dataset, fieldheader, datavalues)

	!*  Subroutine arguments
	type(TecplotDataset), intent(in)     :: dataset
	type(TecplotZoneHeader), intent(in)  :: fieldheader
	real(wp), intent(out), optional      :: datavalues(:,:)
	!*  End Subroutine arguments

	!   Local scalars:
	integer :: i

!- End of header --------------------------------------------------------------

	select case (fieldheader%datapacking)
	case ("POINT")	

		if (present(datavalues)) then		
			do i=1,fieldheader%no_of_points
				read(dataset%id,*) datavalues(i,:)
			end do
		else
			do i=1,fieldheader%no_of_points
				read(dataset%id,*)
			end do
		end if

	case ("BLOCK")
		
		if (present(datavalues)) then		
			call readBLOCK(ncid=dataset%id, nvars=dataset%nvars,    &
						   no_of_points=fieldheader%no_of_points,   &
						   datavalues=datavalues)				
		else
			call readBLOCK(ncid=dataset%id, nvars=dataset%nvars,    &
						   no_of_points=fieldheader%no_of_points)				
		end if
	end select

end subroutine read_variables

!******************************************************************************
! readBLOCK
!******************************************************************************
!
! Description:
!   lee los datos de una zona BLOCK
!------------------------------------------------------------------------------
subroutine readBLOCK(ncid, nvars, no_of_points, datavalues)

    !*  Subroutine arguments
    integer, intent(in) 			:: ncid
    integer, intent(in) 			:: nvars
    integer, intent(in) 			:: no_of_points
	real(wp), intent(out), optional :: datavalues(:,:)
    !*  End Subroutine arguments

    ! Local scalars:
    integer :: i
	integer :: ii
	integer :: j
    integer :: k
    integer :: ncols
    integer :: nrows
    integer :: resto
	logical :: isvarline

    ! Local characters:
    character(LENGTH) :: line

!- End of header --------------------------------------------------------------

    do k=1,nvars
		isvarline = .false.
		do
			read(ncid,'(a)') line
		
			do i=1,len_trim(line)
				if (isdigit(line(i:i))) then
					isvarline = .true.
					backspace(unit=ncid)
					exit
				end if
			end do
		
			if (isvarline) exit
		end do		

		! cuenta columnas: busca digito seguido de blanco 
		ncols = 0
		do i=1,len_trim(line)
			if (isdigit(line(i:i))) then
				if (line(i+1:i+1) == " ") ncols = ncols + 1
			end if
		end do

        nrows = no_of_points/ncols
        resto = mod(no_of_points, ncols)

		if (present(datavalues)) then
			ii = 0
		    do i=1,nrows
		        read(ncid,*) (datavalues(ii+j,k), j=1,ncols)
				ii = ii + ncols	
		    end do	    
		    if (resto > 0) then
		        read(ncid,*) datavalues(ii+1:,k)
		    end if
        else
		    do i=1,nrows
		        read(ncid,*)
		    end do	    
		    if (resto > 0) then
		        read(ncid,*)
		    end if
		end if
    end do

end subroutine readBLOCK

!******************************************************************************
! process_datavalues
!******************************************************************************
!
! Description:
!	- Une todos los campos y crea una unica matriz de points y values, 
!		eliminando los puntos repetidos.    
!------------------------------------------------------------------------------
subroutine process_datavalues(dataset, variables, datavalues, grid)

    !* Subroutine arguments:
    type(TecplotDataset), intent(in)     :: dataset
    type(TecplotVariables), intent(in)   :: variables(:)
    real(wp), intent(in)                 :: datavalues(:,:)
    type(TecplotDataGrid), intent(inout) :: grid
    !* End Subroutine arguments:

	! local scalars:
    integer :: i
	integer :: ierr
    integer :: ipoint
    integer :: id
    integer :: ifield
	integer :: no_of_fields
	integer :: no_of_coords
    integer :: no_of_datavalues
	integer :: no_of_points

	! local alloctables:
    integer, allocatable  :: indexpoint(:)

!- End of header --------------------------------------------------------------

	no_of_coords = dataset%dimspace 
	no_of_datavalues = size(datavalues, dim=1)
	no_of_fields = size(grid%fields)

    allocate(indexpoint(no_of_datavalues), stat=ierr)
    indexpoint = 0

    ! Cuenta el numero de id's distintos
    no_of_points = 0
    do ifield=1,no_of_fields
        do i=1,grid%fields(ifield)%fieldheader%no_of_points

            id = grid%fields(ifield)%global_id(i)

            if (indexpoint(id) == 0) then
                no_of_points = no_of_points + 1
                indexpoint(id) = id
            end if
        end do
    end do
    deallocate(indexpoint)

    ! Reserva memoria para points
	allocate(grid%points(no_of_points, no_of_coords), stat=ierr)
    if (ierr /= 0) then
        call reporterror(MODname//': process_datavalues: &
                         &Unsuccessful allocation of points.')
    end if

    ! Reserva memoria para values
	allocate(grid%values(no_of_points, no_of_coords), stat=ierr)
    if (ierr /= 0) then
        call reporterror(MODname//': process_datavalues: &
                         &Unsuccessful allocation of values.')
    end if

    ! carga points y values
    ipoint = 0
    do ifield=1,no_of_fields
        do i=1,size(grid%fields(ifield)%global_id)
            ipoint = ipoint + 1            
            grid%points(grid%fields(ifield)%global_id(i),:) = &
                                    datavalues(ipoint,1:no_of_coords)
            grid%values(grid%fields(ifield)%global_id(i),:) = &
                                    datavalues(ipoint,no_of_coords+1:2*no_of_coords)
        end do
    end do

end subroutine process_datavalues

!******************************************************************************
! tp_inquire_griddata
!******************************************************************************
!
! Description:
!    Extrae griddata
!------------------------------------------------------------------------------
subroutine tp_inquire_griddata(grid, griddata)

    !* Subroutine arguments:
    type(TecplotDataGrid), intent(in)  :: grid
    type(TecplotGriddata), intent(out) :: griddata
    !* End Subroutine arguments:
    
	! Local scalars:
	integer :: i

!- End of header --------------------------------------------------------------

    griddata%no_of_points = size(grid%points, dim=1)

    griddata%no_of_elements = 0
    griddata%no_of_tetraeders = 0
    griddata%no_of_hexaeders = 0
    griddata%no_of_triangles = 0
    griddata%no_of_quadrilaterals = 0
    griddata%no_of_linesegs = 0

	do i=1,size(grid%fields)
		griddata%no_of_elements = griddata%no_of_elements +      			&
								  grid%fields(i)%fieldheader%no_of_elements

		select case (grid%fields(i)%elementType)
		case (tetra)
			griddata%no_of_tetraeders = griddata%no_of_tetraeders + 		&
								   grid%fields(i)%fieldheader%no_of_elements

		case (hexa)
			griddata%no_of_hexaeders = griddata%no_of_hexaeders + 		    &
								   grid%fields(i)%fieldheader%no_of_elements

		case (triangle)
			griddata%no_of_triangles = griddata%no_of_triangles +   		&
								   grid%fields(i)%fieldheader%no_of_elements

		case (quadrilateral)
			griddata%no_of_quadrilaterals = griddata%no_of_quadrilaterals +	&
								   grid%fields(i)%fieldheader%no_of_elements

		case (lineseg)
			griddata%no_of_linesegs = griddata%no_of_linesegs + 			&
								   grid%fields(i)%fieldheader%no_of_elements

		end select
	end do

    griddata%points_per_tetraeder = 4
    griddata%points_per_hexaeder = 8
    griddata%points_per_triangle = 3
    griddata%points_per_quadrilateral = 4
    griddata%points_per_lineseg = 2
 
end subroutine tp_inquire_griddata

!******************************************************************************
! tf90_close
!******************************************************************************
!
! Description:
!    Close the file
!------------------------------------------------------------------------------
subroutine tp_close(self)

    !* Subroutine arguments:
    type(TecplotFileClass), intent(in) :: self
    !* End Subroutine arguments:
    
!- End of header --------------------------------------------------------------
    
    close(unit=self%dataset%id)

end subroutine tp_close

!******************************************************************************
! Tecplotfile_initialize
!******************************************************************************
!
! Description: 
!   Inicializa las componentes de self
!------------------------------------------------------------------------------
subroutine Tecplotfile_initialize(self)

    !* Subroutine arguments:
    type(TecplotFileClass), intent(out) :: self
    !* End Subroutine arguments:

!- End of header --------------------------------------------------------------
    
    self%filename = ""
    
    self%dataset%id = 0    
    self%dataset%name = ""
    self%dataset%dimspace = 0
    self%dataset%nvars = 0

    self%variables => null()
    
    call initialize(grid=self%grid)
    
end subroutine Tecplotfile_initialize

!******************************************************************************
! initialize_griddata
!******************************************************************************
!
! Description:
!    initialize griddata.
!------------------------------------------------------------------------------
subroutine initialize_griddata(griddata)

    !* Subroutine arguments:
    type(TecplotGriddata), intent(out) :: griddata
    !* End Subroutine arguments:
    
!- End of header --------------------------------------------------------------

    griddata%no_of_points = 0
    griddata%no_of_elements = 0
    griddata%no_of_tetraeders = 0
    griddata%no_of_hexaeders = 0
    griddata%no_of_triangles = 0
    griddata%no_of_quadrilaterals = 0
    griddata%no_of_linesegs = 0
    griddata%points_per_tetraeder = 0
    griddata%points_per_hexaeder = 0
    griddata%points_per_triangle = 0
    griddata%points_per_quadrilateral = 0
    griddata%points_per_lineseg = 0

end subroutine initialize_griddata

!******************************************************************************
! initialize_grid
!******************************************************************************
!
! Description:
!    initialize grid.
!
!------------------------------------------------------------------------------
subroutine initialize_grid(grid)

    !* Subroutine arguments:
    type(TecplotDatagrid), intent(out) :: grid
    !* End Subroutine arguments:
    
!- End of header --------------------------------------------------------------

    grid%points  => null()
    grid%values  => null()
    grid%fields  => null()

end subroutine initialize_grid

!******************************************************************************
! initialize_variables
!******************************************************************************
!
! Description:
!    alocata e initialize variables
!------------------------------------------------------------------------------
subroutine initialize_variables(variables, nvars)

    !* Subroutine arguments:
    type(TecplotVariables), pointer :: variables(:)
    integer, intent(in)             :: nvars
    !* End Subroutine arguments:
    
    ! Local scalars:
    integer :: i
    integer :: ierr

!- End of header --------------------------------------------------------------

    allocate(variables(nvars), stat=ierr)
    if (ierr /= 0) then
        call reporterror(MODname//': initialize_variables reports: &
                            &Unsuccessful allocation of self%variables.')
    end if

    do i=1,nvars
        variables(i)%name = ""
        variables(i)%id = 0
    end do

end subroutine initialize_variables

!******************************************************************************
! initialize_fields
!******************************************************************************
!
! Description:
!    save memory for fields and initialize grid.
!------------------------------------------------------------------------------
subroutine initialize_fields(no_fields, fields)

    !* Subroutine arguments:
    integer, intent(in)             :: no_fields
    type(TecplotDataField), pointer :: fields(:)
    !* End Subroutine arguments:
    
    ! Local scalars:
    integer :: i
    integer :: ierr
    
!- End of header --------------------------------------------------------------

    ! reserva memoria para fields (todos campos)
    allocate(fields(no_fields), stat=ierr)
    if (ierr /= 0) then
      call reporterror(MODname//': initialize_fields reports: &
                        &Unsuccessful allocation of self%grid%fields.')
    end if

    do i=1,no_fields		
		fields(i)%fieldheader%header = ""
		fields(i)%fieldheader%IJK => null()
		fields(i)%fieldheader%dataPacking = ""
		fields(i)%fieldheader%zoneType = ""
		fields(i)%fieldheader%no_of_points = 0
		fields(i)%fieldheader%no_of_elements = 0
		fields(i)%elementType = ""
        fields(i)%connections => null()
        fields(i)%global_id => null()
    end do

end subroutine initialize_fields

!******************************************************************************
! saveTecplotDataFile
!******************************************************************************
!
! Description:
!   
!------------------------------------------------------------------------------
subroutine saveTecplotDataFile(self, fileName, calidad)

    !*  Subroutine arguments
    type(TecplotFileClass), intent(in) :: self
    character(*), intent(in)           :: fileName
    type(calidadfields), intent(in) :: calidad
    !*  End Subroutine arguments

    ! Local scalars:
    integer :: unit
    integer :: i

    character(LENGTH) :: variables

!- End of header --------------------------------------------------------------

    ! get file free unit
    unit = freeunit()

    ! open file
    open(unit, file=fileName, status='unknown', action='write')

    ! write TITLE
    write(unit,'(a)') 'TITLE = "'//trim(self%dataset%name)//'"'

    ! write VARIABLES
    variables = 'VARIABLES = '
    do i=1,self%dataset%nvars
        variables = trim(variables)//' "'//trim(self%variables(i)%name)//'"'    
    end do
    variables = trim(variables)//' "f_shape" "f_skew" "f_relative_size" '
    write(unit,'(a)') trim(variables)

    ! write fileds
    do i=1,size(self%grid%fields)
        !write(unit,'(a)') trim(self%grid%fields(i)%fieldheader%header)
        write(unit,'(5A15,2X,I5,2X,A5,2X,I5,4A35)') 'ZONE T=', '"',  trim(self%grid%fields(i)%elementType), '"', 'N=', &
            self%grid%fields(i)%fieldheader%no_of_points, 'E=', self%grid%fields(i)%fieldheader%no_of_elements, &
            'DATAPACKING = BLOCK', 'ZONETYPE = ', self%grid%fields(i)%fieldheader%zonetype, 'VARLOCATION=([6-8]=CELLCENTERED)'
        
        self%grid%fields(i)%fieldheader%datapacking = 'BLOCK'
        
        select case (self%grid%fields(i)%fieldheader%datapacking)
        case ('POINT')
            call write_point(unit, self%grid%points, self%grid%values,  &
                             self%grid%fields(i))        
        case ('BLOCK')
            call write_block(unit, self%grid%points, self%grid%values,  & 
                             self%grid%fields(i), calidad%fields(i))        
        end select

        select case (self%grid%fields(i)%fieldheader%zonetype)
        case ('ORDERED')
        case default
            call write_connections(unit, self%grid%fields(i))
        end select    
    end do

end subroutine saveTecplotDataFile

!******************************************************************************
! write_point
!******************************************************************************
!
! Description:
!   
!------------------------------------------------------------------------------
subroutine write_point(unit, points, values, field)

    !*  Subroutine arguments
    integer, intent(in)                :: unit
    real(wp), intent(in)               :: points(:,:)
    real(wp), intent(in)               :: values(:,:)
    type(TecplotDataField), intent(in) :: field
    !*  End Subroutine arguments

    !   Local scalars:
    integer  :: i
    character(LENGTH) :: formato

!- End of header --------------------------------------------------------------

    select case (size(points, dim=2))
    case (2)
        formato = '(4E25.16,i8)'
    case (3)
        formato = '(6E25.16,i8)'
    end select

    do i=1,field%fieldheader%no_of_points
        write(unit,fmt=formato) points(field%global_id(i),:),       &
                                values(field%global_id(i),:),       &
                                field%global_id(i)
    end do

end subroutine write_point

!******************************************************************************
! write_block
!******************************************************************************
!
! Description:
!   
!------------------------------------------------------------------------------
subroutine write_block(unit, points, values, field, calidad)

    !*  Subroutine arguments
    integer, intent(in)                :: unit
    real(wp), intent(in)               :: points(:,:)
    real(wp), intent(in)               :: values(:,:)
    type(TecplotDataField), intent(in) :: field
    type(calidadElementos), intent(in) :: calidad
    !*  End Subroutine arguments

    !   Local scalars:
    integer :: i
    integer :: ii    
    integer :: k
    integer :: ndim
    integer :: ncols
    integer :: nrows, nrows_calidad
    integer :: resto, resto_calidad

!- End of header --------------------------------------------------------------

    ncols = 5
    nrows = field%fieldheader%no_of_points/ncols
    resto = mod(field%fieldheader%no_of_points,ncols)
    
    nrows_calidad = field%fieldheader%no_of_elements/ncols
    resto_calidad = mod(field%fieldheader%no_of_elements,ncols)

    ndim = size(points,2)

    ! escribe points
    do k=1,ndim
        i = 1
        do ii=1,nrows
            write(unit,'(10E25.16)') points(field%global_id(i:i+4),k)
            i = i + 5
        end do
        if (resto > 0) write(unit,'(10E25.16)') points(field%global_id(i:),k)
    end do

    ! escribe values
    do k=1,ndim
        i = 1
        do ii=1,nrows
            write(unit,'(10E25.16)') values(field%global_id(i:i+4),k)
            i = i + 5
        end do
        if (resto > 0) write(unit,'(10E25.16)') values(field%global_id(i:),k)
    end do

    ! escribe global_id
    i = 1
    do ii=1,nrows
        write(unit,'(5i8)') field%global_id(i:i+4)
        i = i + 5
    end do
    if (resto > 0) write(unit,'(10i8)') field%global_id(i:)
    
    !escribe f_shape
    i = 1
    do ii=1,nrows_calidad
        write(unit,'(10E25.16)') calidad%f_shape(i:i+4)
        i = i + 5
    end do
    if (resto_calidad > 0) write(unit,'(10E25.16)') calidad%f_shape(i:)
    
    !escribe f_skew
    i = 1
    do ii=1,nrows_calidad
        write(unit,'(10E25.16)') calidad%f_skew(i:i+4)
        i = i + 5
    end do
    if (resto_calidad > 0) write(unit,'(10E25.16)') calidad%f_skew(i:)
    
    !escribe f_relative_size
    i = 1
    do ii=1,nrows_calidad
        write(unit,'(10E25.16)') calidad%f_relative_size(i:i+4)
        i = i + 5
    end do
    if (resto_calidad > 0) write(unit,'(10E25.16)') calidad%f_relative_size(i:)

end subroutine write_block

!******************************************************************************
! write_connections
!******************************************************************************
!
! Description:
!   
!------------------------------------------------------------------------------
subroutine write_connections(unit, field)

    !*  Subroutine arguments
    integer, intent(in)                :: unit
    type(TecplotDataField), intent(in) :: field
    !*  End Subroutine arguments

    !   Local scalars:
    integer  :: i

!- End of header --------------------------------------------------------------

    do i=1,field%fieldheader%no_of_elements
        write(unit,'(10i8)') field%connections(i,:)
    end do

end subroutine write_connections

!******************************************************************************
! write_resumefile
!******************************************************************************
!
! Description:
!   escribe los datos del fichero Tecplot
!------------------------------------------------------------------------------
subroutine write_resumefile(self)

    !*  Subroutine arguments
    type(TecplotFileClass), intent(in) :: self
    !*  End Subroutine arguments

    ! Local scalars:
    integer :: unit
    integer :: i

!- End of header --------------------------------------------------------------

    ! Pantalla
    unit = 6
    do i=1,2
    !---------------------------------------------------------
    ! escribe en pantalla toda la informacion recopilada
    !---------------------------------------------------------
        write(unit,*)
        write(unit,*) '-----------------------------------------------'
        write(unit,*) ' tecplot file information                      '
        write(unit,*) '-----------------------------------------------'
        write(unit,*)
        write(unit,'(2a)')   '    file name                : ', trim(self%filename)
        write(unit,*)
        write(unit,'(a,i4)') '    # dimensions             : ', self%dataset%dimspace
        write(unit,'(a,i4)') '    # variables              : ', self%dataset%nvars
        write(unit,*)
        write(unit,*) '-----------------------------------------------'
        write(unit,*) ' mesh information                              '
        write(unit,*) '-----------------------------------------------'
        write(unit,*)
        write(unit,'(a,i8)') '    # nodes                  : ', self%griddata%no_of_points
        write(unit,'(a,i8)') '    # total elements         : ', self%griddata%no_of_elements
        write(unit,*)
        write(unit,'(a,i8)') '    # volume elements        : ', self%griddata%no_of_tetraeders +     &
															    self%griddata%no_of_hexaeders
        write(unit,'(a,i8)') '       # tetrahedra          : ', self%griddata%no_of_tetraeders
        write(unit,'(a,i8)') '       # hexaheders          : ', self%griddata%no_of_hexaeders
        write(unit,*)
        write(unit,'(a,i8)') '    # surface elements       : ', self%griddata%no_of_triangles +      &
																self%griddata%no_of_quadrilaterals + &
																self%griddata%no_of_linesegs
        write(unit,'(a,i8)') '       # triangles           : ', self%griddata%no_of_triangles
        write(unit,'(a,i8)') '       # quadrilaterals      : ', self%griddata%no_of_quadrilaterals
        write(unit,'(a,i8)') '       # linesegs            : ', self%griddata%no_of_linesegs
        write(unit,*)

        ! Fichero
        unit = get_resumefile()
    end do

end subroutine write_resumefile

end module TecplotFileUtility
