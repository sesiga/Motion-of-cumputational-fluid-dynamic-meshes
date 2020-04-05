!********************************************************************
!QualityUtility
!********************************************************************
!
!Performs operations to support MeshQuality module
    
module QualityUtility

   !modules used
   use precision, only : wp
   
   use utilities, only : reporterror, &
                         vector_product, &
                         det2, &
                         det3, &
                         MatrixProduct, &
                         InvertMatrix
   
   use QualityDataClass
   
   use QualityMetrics
   !end modules used
   
   implicit none
   
   !private parameters
   character(*), parameter, private :: MODname = 'QualityUtility'
   
    contains
    !end of header-------------------------------------
    
    !*****************************************************************
    !JacobianMatrix2D
    !*****************************************************************
    !
    !Creates matrices A, W: triangle
    subroutine JacobianMatrix2D(PA, PW, matrix)
    
       !subroutine arguments
       type(QualityMatrices), intent(out) :: matrix
       
       type(FigurePoints), intent(in) :: PA(:) !points to form matrix A
       type(FigurePoints), intent(in) :: PW(:) !points to form matrix W
       !end subroutine arguments
       
       !local scalars
       real(wp) :: M(2,2)
       
       integer :: i
       
       !end of header--------------------------------------------
       
           !defines matrix M
           M = reshape((/-1d0, 1d0, &
                         -1d0, 0d0/), (/2,2/))
           
           !creates matrices A
           matrix%jacobian(1)%A = reshape((/PA(2)%dim2(1)-PA(1)%dim2(1), PA(2)%dim2(2)-PA(1)%dim2(2), &
                                            PA(3)%dim2(1)-PA(1)%dim2(1), PA(3)%dim2(2)-PA(1)%dim2(2)/), (/2,2/))
                                            
           do i=2,3
               matrix%jacobian(i)%A = MatrixProduct(matrix%jacobian(i-1)%A, M)
           end do
                                                       
           !creates matrices W
           matrix%jacobian(1)%W = reshape((/PW(2)%dim2(1)-PW(1)%dim2(1), PW(2)%dim2(2)-PW(1)%dim2(2), &
                                            PW(3)%dim2(1)-PW(1)%dim2(1), PW(3)%dim2(2)-PW(1)%dim2(2)/), (/2,2/))
                                            
           do i=2,3
               matrix%jacobian(i)%W = MatrixProduct(matrix%jacobian(i-1)%W, M)
           end do
                  
    end subroutine JacobianMatrix2D
    !*****************************************************************************
    
    !*****************************************************************************
    !DecomposeMatrix2D
    !*****************************************************************************
    !
    !Decompose A, W into relevant matrices: triangle
    subroutine DecomposeMatrix2D(matrix)
    
       !subroutine arguments
       type(QualityMatrices), intent(inout) :: matrix
       !end subroutine arguments
       
       !local scalars
       real(wp) :: alpha !det(A)
       real(wp) :: omega !det(W)
       real(wp) :: s1, s2
       
       integer :: i
       
       !end of header----------------------------------
       
       !det(A)
       alpha = det2(matrix%jacobian(1)%A(:,1), matrix%jacobian(1)%A(:,2))
       
       !det(W)
       omega = det2(matrix%jacobian(1)%W(:,1), matrix%jacobian(1)%W(:,2))
       
       do i=1,3
           s1 = sqrt(matrix%jacobian(i)%lambda(1,1)*matrix%jacobian(i)%lambda(2,2))
           
           !matrix Q, column 1
           matrix%decompose(i)%Q(1,1) = 1.0_wp
           matrix%decompose(i)%Q(2,1) = 0.0_wp
           !column 2
           matrix%decompose(i)%Q(1,2) = matrix%jacobian(i)%lambda(1,2)/s1
           matrix%decompose(i)%Q(2,2) = alpha/s1
           
           !matrix D, column 1
           matrix%decompose(i)%D(1,1) = 1.0_wp
           matrix%decompose(i)%D(2,1) = 0.0_wp
           !column 2
           matrix%decompose(i)%D(1,2) = 0.0_wp
           matrix%decompose(i)%D(2,2) = sqrt(matrix%jacobian(i)%lambda(2,2)/matrix%jacobian(i)%lambda(1,1))
           
           s2 = sqrt(matrix%jacobian(i)%lambdaW(1,1)*matrix%jacobian(i)%lambdaW(2,2))
           
           !matrix Qw, column 1
           matrix%decompose(i)%Qw(1,1) = 1.0_wp
           matrix%decompose(i)%Qw(2,1) = 0.0_wp
           !column 2
           matrix%decompose(i)%Qw(1,2) = matrix%jacobian(i)%lambdaW(1,2)/s2
           matrix%decompose(i)%Qw(2,2) = omega/s2
           
           !matrix Dw, column 1
           matrix%decompose(i)%Dw(1,1) = 1.0_wp
           matrix%decompose(i)%Dw(2,1) = 0.0_wp
           !column 2
           matrix%decompose(i)%Dw(1,2) = 0.0_wp
           matrix%decompose(i)%Dw(2,2) = sqrt(matrix%jacobian(i)%lambdaW(2,2)/matrix%jacobian(i)%lambdaW(1,1))
       end do
           
    end subroutine DecomposeMatrix2D
    !******************************************************************************
    
    !******************************************************************************
    !Tmatrix
    !******************************************************************************
    !
    !calculates matrix T
    subroutine Tmatrix(matrix)
    
       !subroutine arguments
       type(QualityMatrices), intent(inout) :: matrix
       !end subroutine arguments
       
       !local scalars
       real(wp) :: A(size(matrix%jacobian(1)%A, dim=1), size(matrix%jacobian(1)%A, dim=2))
       
       !end of header-------------------------------------
       
       !invert matrix W
       A = InvertMatrix(matrix%jacobian(1)%W)
       
       !calculates matrix T=A*Winv
       matrix%T = MatrixProduct(matrix%jacobian(1)%A, A)
       
    end subroutine Tmatrix
    !******************************************************************************
    
    !******************************************************************************
    !LambdaMatrix
    !
    !creates lambda matrices
    subroutine LambdaMatrix(matrix)
    
       !subroutine arguments
       type(QualityMatrices), intent(inout) :: matrix
       !end subroutine arguments
       
       !local scalars
       integer :: i
       
       !end of header--------------------------------------
       
       !calculates matrix lambda
       do i=1,size(matrix%jacobian)
           matrix%jacobian(i)%lambda = MatrixProduct(transpose(matrix%jacobian(i)%A), matrix%jacobian(i)%A)
           matrix%jacobian(i)%lambdaW = MatrixProduct(transpose(matrix%jacobian(i)%W), matrix%jacobian(i)%W)
       end do
       
    end subroutine
    !*****************************************************************************
    
    !*****************************************************************************
    !Rmatrix2D
    !
    !calculates matrix R: triangle
    subroutine Rmatrix2D(matrix)
    
       !subroutine arguments
       type(QualityMatrices), intent(inout) :: matrix
       !end subroutine arguments
       
       !local scalars
       
       !end of header-----------------------------------
       
       !calculates matrix R
       matrix%R(1,1) = matrix%jacobian(1)%A(1,1)
       matrix%R(2,1) = matrix%jacobian(1)%A(2,1)
       matrix%R(1,2) = - matrix%jacobian(1)%A(2,1)
       matrix%R(2,2) = matrix%jacobian(1)%A(1,1)
       
       matrix%R = matrix%R/sqrt(matrix%jacobian(1)%A(1,1))
       
       !calculates matrix Rw
       matrix%Rw(1,1) = matrix%jacobian(1)%W(1,1)
       matrix%Rw(2,1) = matrix%jacobian(1)%W(2,1)
       matrix%Rw(1,2) = - matrix%jacobian(1)%W(2,1)
       matrix%Rw(2,2) = matrix%jacobian(1)%W(1,1)
       
       matrix%Rw = matrix%Rw/sqrt(matrix%jacobian(1)%W(1,1))
       
    end subroutine Rmatrix2D
    !*****************************************************************************
    
    !*****************************************************************************
    !Matrix_cuadrilateros
    !
    !crea la matriz jacobiana para el caso de cuadrilateros
    subroutine Matrix_cuadrilateros(nodesIn_A, nodesIn_W, matrix)
    
       !subroutine arguments 
       real(wp), intent(in) :: nodesIn_A(:,:)
       real(wp), intent(in) :: nodesIn_W(:,:)
       type(matricesCuadrilateros), intent(inout) :: matrix
       !end subroutine arguments
       
       !local scalars
       integer :: i
       
       real(wp) :: nodes_A(7,2) !matriz que contiene a los nodesIn y repite los tres primeros para facilitar el calculo
       real(wp) :: nodes_W(7,2)
       
       !end of header---------------------------------------
       
       !crea nodes
       do i=1,4
           nodes_A(i,:) =  nodesIn_A(i,:)
           nodes_W(i,:) =  nodesIn_W(i,:)
       end do
       nodes_A(5,:) = nodesIn_A(1,:) ; nodes_W(5,:) = nodesIn_W(1,:)
       nodes_A(6,:) = nodesIn_A(2,:) ; nodes_W(6,:) = nodesIn_W(2,:)
       nodes_A(7,:) = nodesIn_A(3,:) ; nodes_W(7,:) = nodesIn_W(3,:)
    
       !crea las matrices
       do i=1,4
           
           !matriz jacobiana: A
           matrix%A(i)%matriz(1,1) = nodes_A(i+1,1) - nodes_A(i,1)
           matrix%A(i)%matriz(2,1) = nodes_A(i+1,2) - nodes_A(i,2)
           matrix%A(i)%matriz(1,2) = nodes_A(i+3,1) - nodes_A(i,1)
           matrix%A(i)%matriz(2,2) = nodes_A(i+3,2) - nodes_A(i,2)
           
           !matriz jacobiana: W
           matrix%W(i)%matriz(1,1) = nodes_W(i+1,1) - nodes_W(i,1)
           matrix%W(i)%matriz(2,1) = nodes_W(i+1,2) - nodes_W(i,2)
           matrix%W(i)%matriz(1,2) = nodes_W(i+3,1) - nodes_W(i,1)
           matrix%W(i)%matriz(2,2) = nodes_W(i+3,2) - nodes_W(i,2)
           
           !tensor metrico: MT, MTw
           matrix%MT(i)%matriz = matrixProduct(transpose(matrix%A(i)%matriz), matrix%A(i)%matriz)
           
           matrix%MTw(i)%matriz = matrixProduct(transpose(matrix%W(i)%matriz), matrix%W(i)%matriz)
           
           !determinante de W
           !matrix%W(i)%det = det2(matrix%W(i)%matriz(:,1), matrix%W(i)%matriz(:,2))
           matrix%W(i)%det = sqrt(trace(matrix%MTw(i)%matriz))
           
           !determinante de A
           matrix%A(i)%det = sqrt(trace(matrix%MT(i)%matriz))
           
           !matriz Winv
           matrix%Winv(i)%matriz(1,1) = matrix%W(i)%matriz(2,2)/matrix%W(i)%det
           matrix%Winv(i)%matriz(2,1) = -matrix%W(i)%matriz(2,1)/matrix%W(i)%det
           matrix%Winv(i)%matriz(1,2) = -matrix%W(i)%matriz(1,2)/matrix%W(i)%det
           matrix%Winv(i)%matriz(2,2) = matrix%W(i)%matriz(1,1)/matrix%W(i)%det
                      
           !determinante de A
           !matrix%A(i)%det = det2(matrix%A(i)%matriz(:,1), matrix%A(i)%matriz(:,2))
           
           !matriz T
           matrix%T(i)%matriz = matrixProduct(matrix%A(i)%matriz, matrix%Winv(i)%matriz)
           
           !matriz Q
           matrix%Q(i)%matriz(1,1) = 1_wp
           matrix%Q(i)%matriz(2,1) = 0_wp
           matrix%Q(i)%matriz(1,2) = matrix%MT(i)%matriz(1,2)/sqrt(matrix%MT(i)%matriz(1,1)*matrix%MT(i)%matriz(2,2))
           matrix%Q(i)%matriz(2,2) = matrix%A(i)%det/sqrt(matrix%MT(i)%matriz(1,1)*matrix%MT(i)%matriz(2,2))
           
           !matriz Qw
           matrix%Qw(i)%matriz(1,1) = 1_wp
           matrix%Qw(i)%matriz(2,1) = 0_wp
           matrix%Qw(i)%matriz(1,2) = matrix%MTw(i)%matriz(1,2)/sqrt(matrix%MTw(i)%matriz(1,1)*matrix%MTw(i)%matriz(2,2))
           matrix%Qw(i)%matriz(2,2) = matrix%W(i)%det/sqrt(matrix%MTw(i)%matriz(1,1)*matrix%MTw(i)%matriz(2,2))
           
           !determinante Qw
           !matrix%Qw(i)%det = det2(matrix%Qw(i)%matriz(:,1), matrix%Qw(i)%matriz(:,2))
           matrix%Qw(i)%det = sqrt(trace(matrixProduct(transpose(matrix%Qw(i)%matriz),matrix%Qw(i)%matriz)))
           
           !matriz Qwinv
           matrix%Qwinv(i)%matriz(1,1) = matrix%Qw(i)%matriz(2,2)/matrix%Qw(i)%det
           matrix%Qwinv(i)%matriz(2,1) = -matrix%Qw(i)%matriz(2,1)/matrix%Qw(i)%det
           matrix%Qwinv(i)%matriz(1,2) = -matrix%Qw(i)%matriz(1,2)/matrix%Qw(i)%det
           matrix%Qwinv(i)%matriz(2,2) = matrix%Qw(i)%matriz(1,1)/matrix%Qw(i)%det
           
       end do 
    
    end subroutine Matrix_cuadrilateros

end module QualityUtility