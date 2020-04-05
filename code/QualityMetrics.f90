!QualityMetrics
!
!this module contains the procedures to calculate
!the quality of the mesh
!**********************************************************************************

module QualityMetrics

   !modules used
   use precision, only : wp
   
   use utilities, only : MatrixProduct, &
                         InvertMatrix, &
                         det2, &
                         det3
   
   use QualityDataClass
   !end modules used
   
   implicit none
   
   character(*), private, parameter :: MODname = 'QualityMetrics'
   
    contains
    !end of header--------------------------------------
    
    !*************************************************************************
    !MetricShape
    !
    function MetricShape(matrix, dim) result(shape)
    
       !function arguments
       type(QualityMatrices), intent(in) :: matrix
       integer, intent(in) :: dim !dimension
       
       real(wp) :: shape
       !end function arguments
       
       !end of header---------------------------------
       
       shape =  dim/kappa(matrix%T)
       
    end function MetricShape
    
    !************************************************************************
    !MetricSkew
    !
    function MetricSkew(matrix, dim) result(Skew)
    
       !function arguments
       type(QualityMatrices), intent(in) :: matrix
       integer, intent(in) :: dim
       
       real(wp) :: Skew
       !end function arguments
       
       !local scalars
       integer :: i
       
       !end of header----------------------------
       
       Skew = 1.0_wp
       do i=1,size(matrix%decompose)
           Skew = Skew*(dim/kappa(MatrixProduct(matrix%decompose(i)%Q, InvertMatrix(matrix%decompose(i)%Qw))))
       end do
       
    end function MetricSkew
    
    !**************************************************************************
    !MetricLengthRatio
    !
    function MetricLengthRatio(matrix, dim) result(LR)
    
       !function arguments
       type(QualityMatrices), intent(in) :: matrix
       integer, intent(in) :: dim
       
       real(wp) :: LR
       !end function arguments
       
       !local scalars
       integer :: i
       
       !end of header--------------------------------
       
       LR = 1.0_wp
       do i=1,size(matrix%decompose)
           LR = LR*(dim/kappa(MatrixProduct(matrix%decompose(i)%D, InvertMatrix(matrix%decompose(i)%Dw))))
       end do
       
    end function MetricLengthRatio
    
    !**************************************************************************
    !MetricOrientation
    !
    function MetricOrientation(matrix, dim) result(Ori)
    
       !function arguments
       type(QualityMatrices), intent(in) :: matrix
       integer, intent(in) :: dim
       
       real(wp) :: Ori
       !end function arguments
       
       !local scalars
       real(wp) :: Xo(size(matrix%T, dim=1), size(matrix%T, dim=2))
       
       !end of header---------------------------
       
       Xo = MatrixProduct(matrix%R, InvertMatrix(matrix%Rw))
       
       Ori = 1.0_wp+(trace(Xo)-dim)/4.0_wp
       
    end function MetricOrientation
    
    !**************************************************************************
    !MetricVolume
    !
    function MetricVolume(matrix, dim) result(Vol)
    
       !function arguments
       type(QualityMatrices), intent(in) :: matrix
       integer, intent(in) :: dim
       
       real(wp) :: Vol
       !end function arguments
       
       !end of header------------------------------
       
       Vol = det2(matrix%T(:,1), matrix%T(:,2))
       
    end function MetricVolume
    
    !**************************************************************************
    !Relative size
    !
    function MetricRelativeSize(matrix) result(RS)
    
       !function arguments
       type(QualityMatrices), intent(in) :: matrix
       real(wp) :: RS
       !end function arguments
       
       !local scalars
       real(wp) :: alpha !det(A)
       real(wp) :: omega !det(W)
       
       !end of header
       
       !det(A)
       alpha = det2(matrix%jacobian(1)%A(:,1), matrix%jacobian(1)%A(:,2))
       
       !det(W)
       omega = det2(matrix%jacobian(1)%W(:,1), matrix%jacobian(1)%W(:,2))
       
       !RS
       RS = min(alpha/omega, omega/alpha)
    
    end function MetricRelativeSize
    
    !**************************************************************************
    !Kappa
    !
    function kappa(A)
    
       !function arguments
       real(wp), intent(in) :: A(:,:)
       
       real(wp) :: kappa
       !end function arguments
       
       !local scalars
       real(wp) :: B(size(A, dim=1), size(A, dim=2))
       real(wp) :: C(size(A, dim=1), size(A, dim=2))
       real(wp) :: D(size(A, dim=1), size(A, dim=2))
       real(wp) :: E(size(A, dim=1), size(A, dim=2))
       
       !end of header----------------------------------
       
       B = MatrixProduct(transpose(A), A)
       C = InvertMatrix(A)
       D = InvertMatrix(transpose(A))
       E = MatrixProduct(C, D)
       
       kappa = sqrt(trace(B))*sqrt(trace(E))
       
    end function kappa
    
    !*************************************************************************
    !trace
    !
    function trace(A)
    
       !function arguments
       real(wp), intent(in) :: A(:,:)
       
       real(wp) :: trace
       !end function arguments
       
       !local scalars
       integer :: i
       
       !end of header------------------------------
       
       trace = 0.0_wp
       do i=1,size(A, dim=1)
           trace = trace + A(i,i)
       end do
       
    end function trace
    
    !******************************************************************************
    !metricShape_cuadrilateros
    !
    function metricShape_cuadrilateros(matrix) result(shape)
    
       !function arguments
       type(matricesCuadrilateros), intent(in) :: matrix
       real(wp) :: shape
       !end function arguments
       
       !local scalars
       integer :: i
       
       !end of header----------------------------------------
       
       shape = 0_wp
       do i=1,4
           shape = shape + 2_wp/kappa(matrix%T(i)%matriz)
       end do
       
       shape = shape/4_wp
              
    end function metricShape_cuadrilateros
    
    !*********************************************************************************
    !metricSkew_cuadrilateros
    !
    function metricSkew_cuadrilateros(matrix) result(skew)
    
       !function arguments
       type(matricesCuadrilateros), intent(in) :: matrix
       real(wp) :: skew
       !end function arguments
       
       !local scalars
       integer :: i
       
       !end of header----------------------------------------
       
       skew = 1_wp
       do i=1,4
           skew = skew*( 2_wp / kappa(matrixProduct(matrix%Q(i)%matriz, matrix%Qwinv(i)%matriz)) )
       end do      
       
    end function metricSkew_cuadrilateros
    
    !********************************************************************
    !metricRelativeSize_cuadrilateros
    !
    function metricRelativeSize_cuadrilateros(matrix) result(RS)
    
       !function arguments
       type(matricesCuadrilateros), intent(in) :: matrix
       real(wp) :: RS
       !end function arguments
       
       !local scalars
       integer :: i
       
       real(wp) :: tau
       
       !end of header----------------------------------------
       
       tau = 0_wp
       do i=1,4
           tau = tau + matrix%A(i)%det/matrix%W(i)%det
       end do
       tau = tau/4_wp
       
       if(tau < 1)then
           RS = tau
       else
           RS = 1_wp/tau
       end if
       
    end function metricRelativeSize_cuadrilateros
    
end module QualityMetrics