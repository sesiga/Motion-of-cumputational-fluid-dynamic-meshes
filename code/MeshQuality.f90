!***************************************************************************
!MeshQuality
!***************************************************************************
!
!This module determinates the quality of the mesh 
    
module MeshQuality

   !modules used
   use precision, only : wp
   
   use utilities, only : reporterror, &
                         removeBlanks
   
   use TecplotFileDataClass
   
   use TecplotGridDataClass
   
   use QualityUtility
   
   use QualityProcess
   
   use QualityMetrics
   
   use QualityDataClass
   !end modules used
   
   implicit none 
   
   !private parameters
   character(*), parameter, private :: MODname = 'MeshQuality'
   
   !public parameters
   
    contains
    
    !end of header----------------------------------------
    
    !**********************************************************************
    !Quality
    !**********************************************************************
    !
    !Implements the procedure to obtain the mesh quality
    subroutine Quality(mesh, calidad)
    
       !subroutine argumentss
       type(TecplotFileClass), intent(in) :: mesh
       type(calidadfields), intent(inout) :: calidad
       !end subroutine arguments
       
       !local scalars
       type(QualityMatrices) :: matrix
       
       character(LENGTH) :: elementType
       
       integer :: n, i
             
       !end of header-----------------------------------------
       
       !calculates quality for each field
       !and its different figures
       do i=1,size(mesh%grid%fields)
           
           !type of figure
           elementType = removeblanks(mesh%grid%fields(i)%elementType)
           
           !number of points of the elements
           n = size(mesh%grid%fields(i)%connections, dim=2)
           
           !allocates matrices
           call QualityPreProcess(matrix, elementType, n)
           
           !differences among different number of points per figure
           select case(elementType)
           case('triangle')
               call QualityTriangle(mesh, matrix, i, calidad)
           case('quadrilateral')
               call QualityQuadrilateral(mesh, i, calidad)
           end select
                      
           !deallocate matrices
           call QualityPostProcess(matrix, elementType, n)
           
       end do
       
    end subroutine Quality
    !**********************************************************************
    
    !**********************************************************************
    !QualityTriangle
    !
    !calculates mesh quality for triangle elements
    subroutine QualityTriangle(mesh, matrix, field, calidad)
    
       !subroutine arguments
       type(TecplotFileClass), intent(in) :: mesh
       type(QualityMatrices), intent(inout) :: matrix
       integer, intent(in) :: field !field del que se calcula la calidad
       type(calidadfields), intent(inout) :: calidad
       !end subroutine arguments
       
       !local scalars
       integer, parameter :: n = 3 !number of points per figure
       integer, parameter :: dim = 2 !dimension
       
       type(FigurePoints) :: A_Points(3) 
       type(FigurePoints) :: W_Points(3)
       
       integer :: global_id(3)
       integer :: i, j, k
       
       !end of header-------------------------------------
       
       !loads points of each figure
           do j=1,size(mesh%grid%fields(field)%connections, dim=1)
               
               !global id of the points of the figure
               global_id(:) = mesh%grid%fields(field)%connections(j,:)
               
               !points to create matrix A, W
               do k=1,n  
                   A_Points(k)%dim2(:) = mesh%grid%points(global_id(k),:) + mesh%grid%values(global_id(k),:)
                   W_points(k)%dim2(:) = mesh%grid%points(global_id(k),:)
               end do
               
               !creates matrices A, W
               call JacobianMatrix2D(A_Points, W_points, matrix)
               
               !creates matrix lambda
               call LambdaMatrix(matrix)
               
               !creates matrix T
               call Tmatrix(matrix)
               
               !creates matrix R
               call Rmatrix2D(matrix)
               
               !creates matrix Q, D
               call DecomposeMatrix2D(matrix)
               
               !calcula la calidad
               !shape
               calidad%fields(field)%f_shape(j) = metricShape(matrix, dim)
               
               !skew
               calidad%fields(field)%f_skew(j) = metricSkew(matrix, dim)
               
               !relative size
               calidad%fields(field)%f_relative_size(j) = metricRelativeSize(matrix)
               
           end do
       
    end subroutine QualityTriangle
    !**********************************************************************
           
    !**********************************************************************
    !QualityQuadrilateral
    !
    !calcula la calidad para los cuadrilateros
    subroutine qualityQuadrilateral(mesh, field, calidad)
    
       !subroutine arguments 
       type(TecplotFileClass), intent(in) :: mesh
       integer, intent(in) :: field
       type(calidadfields), intent(inout) :: calidad
       !end subroutine arguments
       
       !local scalars
       integer :: i, j
       
       real(wp) :: nodesIn_A(4,2) !coordenadas de los nodos
       real(wp) :: nodesIn_W(4,2)
       
       integer :: global_id(4) !global id de los nodos
       
       type(matricesCuadrilateros) :: matrix
       
       !end of header-----------------------------------------
       
       open(1651,file='calidad.dat')
       
       !recorre todo el field
       do i=1,size(mesh%grid%fields(field)%connections, dim=1)
           
           !global id de los nodos
           global_id = mesh%grid%fields(field)%connections(i,:)
           
           !crea el vector con las coordenadas de los nodos del elemento
           do j=1,4
               nodesIn_A(j,:) = mesh%grid%points(global_id(j),:) + mesh%grid%values(global_id(j),:)
               nodesIn_W(j,:) = mesh%grid%points(global_id(j),:) 
           end do
           
           !crea las matrices
           call Matrix_cuadrilateros(nodesIn_A, nodesIn_W, matrix)
           
           !calcula la calidad
           !shape
           calidad%fields(field)%f_shape(i) = metricShape_cuadrilateros(matrix)
           
           !skew
           calidad%fields(field)%f_skew(i) = metricSkew_cuadrilateros(matrix)
           
           !relative size
           calidad%fields(field)%f_relative_size(i) = metricRelativeSize_cuadrilateros(matrix)
           
           
           write(1651,'(i8,3E25.16)') i, calidad%fields(field)%f_shape(i), calidad%fields(field)%f_skew(i), calidad%fields(field)%f_relative_size(i)
           
       end do
       
       close(1651)
       
    end subroutine qualityQuadrilateral
       
end module MeshQuality