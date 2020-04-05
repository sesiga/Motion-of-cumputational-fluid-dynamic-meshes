!***************************************************************************
!QualityProcess
!***************************************************************************
!
!This module allocates and deallocates quality arrays

module QualityProcess

   !modules used
   use precision, only : wp
   
   use utilities, only : reporterror
   
   use QualityDataClass
   !end modules used
   
   implicit none
   
   character(*), parameter, private :: MODname = 'QualityProcess'
   
    contains
    !end of header------------------------------------
    
    !*************************************************************************
    !QualityPreProcess
    !
    !allocates arrays
    subroutine QualityPreProcess(matrix, elementType, n)
    
       !subroutine arguments
       type(QualityMatrices), intent(inout) :: matrix
       
       character(*), intent(in) :: elementType !type of element
       
       integer, intent(in) :: n !number of points per element
       !end subroutine arguments
       
       !local scalars
       integer :: ierr
       integer :: i
       
       !end of header-----------------------------------
       
       select case(elementType)
           
       !case triangle
       case('triangle') 
           
           allocate(matrix%jacobian(3), stat=ierr)
           if(ierr/=0)then
              call reporterror(MODname// ': QualityPreProcess reports: &
                   &Unsuccessful allocating matrix%jacobian')
           end if
           
           allocate(matrix%decompose(3), stat=ierr)
           if(ierr/=0)then
              call reporterror(MODname// ': QualityPreProcess reports: &
                   &Unsuccessful allocating matrix%decompose')
           end if
           
           allocate(matrix%T(2,2), stat=ierr)
           if(ierr/=0)then
              call reporterror(MODname// ': QualityPreProcess reports: &
                   &Unsuccessful allocating matrix%T')
           end if
           
           allocate(matrix%Rw(2,2), stat=ierr)
           if(ierr/=0)then
              call reporterror(MODname// ': QualityPreProcess reports: &
                   &Unsuccessful allocating matrix%Rw')
           end if
           
           allocate(matrix%R(2,2), stat=ierr)
           if(ierr/=0)then
              call reporterror(MODname// ': QualityPreProcess reports: &
                   &Unsuccessful allocating matrix%R')
           end if
           
           do i=1,n
               allocate(matrix%jacobian(i)%A(2,2), stat=ierr)
               if(ierr/=0)then
                  call reporterror(MODname// ': QualityProcess reports: &
                       &Unsuccessful allocating matrix%jacobian%A')
               end if
               
               allocate(matrix%jacobian(i)%W(2,2), stat=ierr)
               if(ierr/=0)then
                  call reporterror(MODname// ': QualityProcess reports: &
                       &Unsuccessful allocating matrix%jacobian%W')
               end if
               
               allocate(matrix%jacobian(i)%lambda(2,2), stat=ierr)
               if(ierr/=0)then
                  call reporterror(MODname// ': QualityProcess reports: &
                       &Unsuccessful allocating matrix%jacobian%lambda')
               end if
               
               allocate(matrix%jacobian(i)%lambdaW(2,2), stat=ierr)
               if(ierr/=0)then
                  call reporterror(MODname// ': QualityProcess reports: &
                       &Unsuccessful allocating matrix%jacobian%lambdaW')
               end if
               
               allocate(matrix%decompose(i)%Q(2,2), stat=ierr)
               if(ierr/=0)then
                  call reporterror(MODname// ': QualityProcess reports: &
                       &Unsuccessful allocating matrix%decompose%Q')
               end if
               
               allocate(matrix%decompose(i)%Qw(2,2), stat=ierr)
               if(ierr/=0)then
                  call reporterror(MODname// ': QualityProcess reports: &
                       &Unsuccessful allocating matrix%decompose%Qw')
               end if
               
               allocate(matrix%decompose(i)%D(2,2), stat=ierr)
               if(ierr/=0)then
                  call reporterror(MODname// ': QualityProcess reports: &
                       &Unsuccessful allocating matrix%decompose%D')
               end if
               
               allocate(matrix%decompose(i)%Dw(2,2), stat=ierr)
               if(ierr/=0)then
                  call reporterror(MODname// ': QualityProcess reports: &
                       &Unsuccessful allocating matrix%decompose%Dw')
               end if
           end do
           
       end select
       
    end subroutine QualityPreProcess
    !*****************************************************************************
    
    !*****************************************************************************
    !QualityPostProcess
    !
    !deallocates arrays
    subroutine QualityPostProcess(matrix, elementType, n)
    
       !subroutine arguments
       type(QualityMatrices), intent(inout) :: matrix
       
       character(*), intent(in) :: elementType 
       
       integer, intent(in) :: n !number of points per element
       !end subroutine arguments
       
       !local  scalars
       integer :: i
       
       !end of header-------------------------------------
       
       select case(elementType)
           
       !case triangle
       case('triangle')
           
           do i=1,n
               deallocate(matrix%jacobian(i)%A)
               deallocate(matrix%jacobian(i)%W)
               deallocate(matrix%jacobian(i)%lambda)
               deallocate(matrix%jacobian(i)%lambdaW)
               deallocate(matrix%decompose(i)%Q)
               deallocate(matrix%decompose(i)%Qw)
               deallocate(matrix%decompose(i)%D)
               deallocate(matrix%decompose(i)%Dw)
           end do
           
           deallocate(matrix%jacobian)
           
           deallocate(matrix%decompose)
           
           deallocate(matrix%T)
           
           deallocate(matrix%R)
           
           deallocate(matrix%Rw)
           
       end select
            
    end subroutine QualityPostProcess
    !*****************************************************************************
    
end module QualityProcess