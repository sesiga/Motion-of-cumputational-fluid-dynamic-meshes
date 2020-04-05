!******************************************************************************
!QualityDataClass
!
!This module contains different arrays used 
!to calculate quality mesh

module QualityDataClass

   !modules used
   use precision, only : wp
   !end modules used
   
   implicit none
   
   !jacobian matrices
   type JacobianMatrices
       real(wp), pointer :: A(:,:)
       real(wp), pointer :: W(:,:)
       real(wp), pointer :: lambda(:,:)
       real(wp), pointer :: lambdaW(:,:)
   end type JacobianMatrices
   
   !decompose matrices
   type DecomposeMatrices
       real(wp), pointer :: Q(:,:)
       real(wp), pointer :: Qw(:,:)
       real(wp), pointer :: D(:,:)
       real(wp), pointer :: Dw(:,:)
   end type
   
   !quality type
   type QualityMatrices
       type(JacobianMatrices), pointer :: Jacobian(:)
       type(DecomposeMatrices), pointer :: Decompose(:) 
       real(wp), pointer :: T(:,:)
       real(wp), pointer :: R(:,:)
       real(wp), pointer :: Rw(:,:)
   end type QualityMatrices
   
   !points type
   type FigurePoints
       real(wp) :: dim2(2) !case 2D
       real(wp) :: dim3(3) !case 3D
   end type FigurePoints
   
   !calidad elementos
   type calidadElementos
       real(wp), pointer :: f_shape(:)
       real(wp), pointer :: f_skew(:)
       real(wp), pointer :: f_relative_size(:)
   end type calidadElementos
   
   type calidadFields
       type(calidadElementos), pointer :: fields(:)
   end type calidadFields
   
   !cuadrilateros
   !matriz 2x2
   type matriz2x2
       real(wp) :: matriz(2,2)
       real(wp) :: det
   end type matriz2x2
         
   type matricesCuadrilateros
       type(matriz2x2) :: A(4) !matriz jacobiana
       type(matriz2x2) :: W(4), Winv(4)
       type(matriz2x2) :: MT(4), MTw(4) !metric tensor: Atras*A, Wtras*W
       type(matriz2x2) :: T(4)
       type(matriz2x2) :: Q(4), Qw(4)
       type(matriz2x2) :: Qinv(4), Qwinv(4)
   end type matricesCuadrilateros
   
end module QualityDataClass