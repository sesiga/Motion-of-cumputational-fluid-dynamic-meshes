!*************************************************************
!RBFprocess
!*************************************************************
!
!this module loads scatter and tecplot data into 
!arrays and viceversa
    
module RBFprocess

   !modules used------------------------------------
   use precision, only : wp
   
   use utilities, only : reporterror
   
   use TecplotFileDataClass
   
   use ScatterFileDataClass
   
   use QualityDataClass
   !end moddules used-------------------------------
   
   implicit none
   
   !private parameters
   character(*), parameter, private :: MODname = 'RBFprocess'
   
   !public 
   public :: PreProcess, PostProcess
   
    contains
    !end of header---------------------------------
    
    !*************************************************************
    !PreProcess
    !*************************************************************
    !
    !Loads initial scatter data into arrays
    subroutine PreProcess(mesh, scatter, PointsOut, &
        ValuesIn, ValuesOut, PointsIn, PointsOut_id, pointsIn_id, calidad)
    
       !subroutine arguments
       type(TecplotFileClass), intent(in) :: mesh
       type(ScatterFileClass), intent(in) :: scatter
       
       real(wp), allocatable, intent(out) :: PointsOut(:,:), &
                                             ValuesIn(:,:), &
                                             ValuesOut(:,:), &
                                             PointsIn(:,:)
       
       integer, allocatable, intent(out) :: pointsOut_id(:)
       integer, allocatable, intent(out) :: pointsIn_id(:)
       type(calidadfields), intent(inout) :: calidad
       !end subroutine arguments
       
       !local scalars       
       integer :: global_id
       integer :: i, j, k, m, n, kmax, jmax, p
       integer :: dim1, dim2
       integer :: ierr
       
       !end of header--------------------------------------
       
       !allocata arrays de calidad
       allocate(calidad%fields( size(mesh%grid%fields) ))
       do i=1,size(calidad%fields)
           allocate(calidad%fields(i)%f_shape( mesh%grid%fields(i)%fieldheader%no_of_elements ))
           allocate(calidad%fields(i)%f_skew( mesh%grid%fields(i)%fieldheader%no_of_elements ))
           allocate(calidad%fields(i)%f_relative_size( mesh%grid%fields(i)%fieldheader%no_of_elements ))
       end do
       
       !allocates and loads mesh points into PointsOut
       dim1 = size(mesh%grid%points, dim=1)-size(scatter%ScatterData%Points, dim=1)
       dim2 = size(mesh%grid%points, dim=2)
       
       allocate(PointsOut(dim1, dim2), stat=ierr)
          if(ierr/=0)then
              call reporterror(MODname// ': PreProcess reports: &
                   &Unsuccessful allocating PointsOut')
          end if
          
      allocate(PointsOut_id(dim1), stat=ierr)
          if(ierr/=0)then
              call reporterror(MODname// ': PreProcess reports: &
                   &Unsuccessful allocating PointsOut_id')
          end if
          
          n = 0
          kmax = size(mesh%grid%fields)
          jmax = size(scatter%ScatterData%Points, dim=1)
          do k=1,kmax
              m = 0
              bucle_2: do i=1,mesh%grid%fields(k)%fieldheader%no_of_points
                  m = m+1
                  do j=1,jmax
                      
                      global_id = mesh%grid%fields(k)%global_id(m)
                  
                      if(global_id == scatter%ScatterData%global_id(j))then
                          exit
                      elseif(j == jmax)then
                          do p=1,n
                              if(pointsOut_id(p) == global_id) cycle bucle_2
                          end do
                          
                          n = n+1
                          !PointsOut(n,:) = mesh%grid%points(m,:)
                          PointsOut(n,:) = mesh%grid%points(global_id,:)
                          PointsOut_id(n) = global_id !global id del nodo en el vector PointsOut
                      end if
                      
                  end do
              end do bucle_2
          end do
          
       !allocates and loads scatter values into ValuesIn
       dim1 = size(scatter%ScatterData%Points, dim=1)
       dim2 = size(scatter%ScatterData%Points, dim=2)
       
       allocate(ValuesIn(dim1, dim2), stat=ierr)
          if(ierr/=0)then
              call reporterror(MODname// ': PreProcess reports: &
                   &Unsuccessful allocating ValuesIn')
          end if
          
          do i=1,dim1
              ValuesIn(i,:) = scatter%ScatterData%Values(i,:)
          end do
          
       !allocates ValuesOut
       dim1 = size(mesh%grid%points, dim=1)-size(scatter%ScatterData%Points, dim=1)
       dim2 = size(mesh%grid%points, dim=2)
       
       allocate(ValuesOut(dim1, dim2), stat=ierr)
          if(ierr/=0)then
              call reporterror(MODname// ': PreProcess reports: &
                   &Unsuccessful allocating ValuesOut')
          end if
          
       !allocates and loads PointsIn
       dim1 = size(scatter%ScatterData%Points, dim=1)
       dim2 = size(scatter%ScatterData%Points, dim=2)
       
       allocate(PointsIn(dim1, dim2), stat=ierr)
        if(ierr/=0)then
              call reporterror(MODname// ': PreProcess reports: &
                   &Unsuccessful allocating PointsIn')
        end if
        
      allocate(PointsIn_id(dim1), stat=ierr)
        if(ierr/=0)then
              call reporterror(MODname// ': PreProcess reports: &
                   &Unsuccessful allocating PointsIn_id')
        end if
        
        do i=1,dim1
            PointsIn(i,:) = scatter%ScatterData%Points(i,:)
            PointsIn_id(i) = scatter%ScatterData%global_id(i) !global id del nodo en el vector PointsIn
        end do
                           
    end subroutine PreProcess
    !******************************************************************************
    
    !******************************************************************************
    !PostProcess
    !******************************************************************************
    !
    !Loads processed data into arrays
    subroutine PostProcess(mesh, scatter, PointsIn, PointsOut, &
        ValuesIn, ValuesOut, PointsOut_id, PointsIn_id)
    
       !subroutine arguments
       type(TecplotFileClass), intent(inout) :: mesh
       type(ScatterFileClass), intent(in) :: scatter
       
       real(wp), allocatable, intent(inout) :: PointsOut(:,:), &
                                               ValuesOut(:,:), &
                                               ValuesIn(:,:), &
                                               PointsIn(:,:)
       
       integer, allocatable, intent(inout) :: pointsOut_id(:)
       integer, allocatable, intent(inout) :: pointsIn_id(:)
       !end subroutine arguments
       
       !local scalars
       integer ::  global_id
       logical :: value_state(size(mesh%grid%values)) ! .true. si ya se ha guardado la deformacion del nodo
       
       integer :: i, j, k, m, n1, n2, kmax, jmax
       
       !end of header---------------------------------------------
       
       !inicializa las variables
       value_state = .false.
       
       !loads deformated values into mesh
       n1 = 0; n2 = 0
       kmax = size(mesh%grid%fields)
       jmax = size(scatter%ScatterData%Points, dim=1)
       bucle_1: do k=1,kmax
          bucle_2: do i=1,mesh%grid%fields(k)%fieldheader%no_of_points
              bucle_3: do j=1,jmax
                  
                  global_id = mesh%grid%fields(k)%global_id(i)
                  
                  if(value_state(global_id)) cycle bucle_2
                  
                  if(global_id == scatter%ScatterData%global_id(j))then
                      n1 = n1+1
                      !mesh%grid%values(global_id,:) = valuesIn(n1,:) !scatter%ScatterData%values(j,:)
                      mesh%grid%values(global_id,:) = valuesIn(j,:)
                      value_state(global_id) = .true.
                      exit
                  elseif(j == jmax)then
                      n2 = n2+1
                      mesh%grid%values(global_id,:) = ValuesOut(n2,:)
                      value_state(global_id) = .true.
                  end if
                      
              end do bucle_3
          end do bucle_2
       end do bucle_1
       
       print *, 'centros', size(pointsin, dim=1)
       print *, 'evaluacion', size(pointsout, dim=1)
       
       !deallocates arrays
       deallocate(PointsOut, ValuesIn, ValuesOut, PointsIn, PointsOut_id, PointsIn_id)
       
        end subroutine PostProcess
    !*********************************************************************************    
        
        !******************************************************************************
    !dealocata calidad
    !
    subroutine qualityDeallocate(calidad)
    
       !subroutine arguments
       type(calidadfields), intent(inout) :: calidad
       !end subroutine arguments
       
       !local scalars
       integer :: i
       
       !end of header----------------------------------------
       do i=1,size(calidad%fields)
           deallocate(calidad%fields(i)%f_shape)
           deallocate(calidad%fields(i)%f_skew)
           deallocate(calidad%fields(i)%f_relative_size)
       end do
       
       deallocate(calidad%fields)
       
    end subroutine qualityDeallocate
                        
end module RBFprocess