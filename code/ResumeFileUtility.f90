!******************************************************************************
! ResumeFileUtility
!******************************************************************************
!
! Description:
!------------------------------------------------------------------------------
module ResumeFileUtility

    use utilities, only:    &
        freeunit
    
    implicit none

    integer, save, private :: unit

contains

!******************************************************************************
! open_resumefile
!******************************************************************************
!
! Description:
!       abre el fichero resume_file y escribe el preambulo
!------------------------------------------------------------------------------
function open_resumefile() result(ierr)

    !* Subroutine arguments:
    integer :: ierr
    !* End Subroutine arguments:
    
    !Local scalars:
    integer :: i
    integer :: funit

!--- End of header---------------------------------------------------------

    ! busca unidad
    unit = freeunit()

    ! abre fichero
    open(unit=unit, file='resume_file.dat', status='unknown', action='write', iostat=ierr)

    funit = 6
    do i=1,2
        write(funit,*)
        write(funit,'(a)') "------------------------------------------"
        write(funit,'(a)') "   Tecplot File Processing Application    "            
        write(funit,'(a)') "------------------------------------------"
        write(funit,*)
        
        ! Fichero
        funit = unit
    end do
    
end function open_resumefile

!******************************************************************************
! get_resumefile
!******************************************************************************
!
! Description:
!       devuelve la unidad del fichero resume_file
!------------------------------------------------------------------------------
function get_resumefile() result(iunit)

    !* Subroutine arguments:
    integer :: iunit
    !* End Subroutine arguments:

!--- End of header---------------------------------------------------------

    iunit = unit

end function get_resumefile

!******************************************************************************
! close_resumefile
!******************************************************************************
!
! Description:
!       cierra el fichero resume_file
!------------------------------------------------------------------------------
function close_resumefile() result(ierr)

    !* Subroutine arguments:
    integer :: ierr
    !* End Subroutine arguments:

    ! Local scalars:
    integer :: i
    integer :: iunit

!--- End of header---------------------------------------------------------

    ! Termination protocol.

    iunit = 6
    do i=1,2
        write(iunit,*)
        write(iunit,'(a)') "The program has finished successfully!"
        write(iunit,*)
        iunit = unit
    end do

    ! cierra fichero
    close(unit=unit, iostat=ierr)

end function close_resumefile

end module ResumeFileUtility
