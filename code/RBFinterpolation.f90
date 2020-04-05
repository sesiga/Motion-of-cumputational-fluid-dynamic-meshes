!****************************************************************************
!RBFinterpolation
!****************************************************************************
!
!description
!this module performs the interpolacion procedure. 
!it uses radial basis functions to get the new points position
!by knowing the deformation of some points
    
module RBFinterpolation

   !modules used----------------------------------------
   use precision, only : wp
   
   use RadialBasisFunctions
   
   use LinearAlgebra
   
   use utilities
   !end modules used------------------------------------

   implicit none
   
   public :: RBFinterpolator
   
   private :: CoeffComputation, &
       EvaluationMatrix  
   
    contains
    !end of header--------------------------------------
    
    !*******************************************************************
    !RBFinterpolator solve the system: h=C·w
    !output: w coeffients vector
    !*******************************************************************
    subroutine RBFinterpolator(RBFfunction, PointsIn, PointsOut, &
        ValuesIn, ValuesOut, rhoCoeff)
    
       !subroutine arguments
       procedure(RBF), pointer :: RBFfunction
       
       real(wp), intent(in) :: PointsOut(:,:), &
                               ValuesIn(:,:), &
                               PointsIn(:,:), &
                               rhoCoeff
       
       real(wp), intent(out) :: ValuesOut(:,:)
       !end subroutine arguments
       
       !local scalars
       integer :: ierr
       
       real(wp) :: h(size(ValuesIn,dim=1)+1, size(ValuesIn,dim=2))
       real(wp) :: C(size(ValuesIn,dim=1)+1, size(valuesIn,dim=1)+1)
       
       real(wp) :: A(size(ValuesOut,dim=1), size(ValuesIn,dim=1)+1)
       
       real(wp) :: rho
       !end local scalars
       
       !end of header------------------------------
       
       !calculates rho coefficient
       !call rhoComputation(PointsIn, rho, rhoCoeff)
       
       !creates h, C arrays
       call CoeffComputation(C, h, RBFfunction, PointsIn, ValuesIn, rhoCoeff)
       
       !solve system h=C·w
       call solve_symmetric_system(C, h, ierr)
       if(ierr<0) print *, 'error while allocatating lapack solver'
       
       !creates A array
       call EvaluationMatrix(A, RBFfunction, PointsOut, PointsIn, rhoCoeff)
       
       !solve system h=A·w
       ValuesOut = matmul(A, h)
       
    end subroutine RBFinterpolator
    !*****************************************************************************
        
    !*****************************************************************************
    !rhoComputation
    !calculates rho coefficient
    !*****************************************************************************
    subroutine rhoComputation(PointsOut, rho, rhoCoeff)
    
       !subroutine arguments
       real(wp), intent(in) :: PointsOut(:,:), &
                               rhoCoeff
       
       real(wp), intent(out) :: rho
       !end subroutine arguments
       
       !locals scalars
       real(wp) :: max2(2), min2(2) !2 dimensions case
       real(wp) :: max3(3), min3(3) !3 dimensions case
       
       !end of header----------------------------------
       
       !calculates extremes at each direction
       !and then calculates rho
       if(size(PointsOut, dim=2) == 2)then    
           
           max2 = maxval(PointsOut, dim=1)
           min2 = minval(PointsOut, dim=1)
           
           rho = norm2(max2-min2)*rhoCoeff
           
       elseif(size(PointsOut, dim=2) == 3)then
           
           max3 = maxval(PointsOut, dim=1)
           min3 = minval(PointsOut, dim=1)
           
           rho = norm2(max3-min3)*rhoCoeff
           
       end if
              
    end subroutine
    !*****************************************************************************
    
    !*****************************************************************************
    !CoeffComputation creates arrays h, C 
    !to perform coefficients computation
    !*****************************************************************************
    subroutine CoeffComputation(C, h, RBFfunction, &
                             PointsIn, ValuesIn, rho)
    
       !subroutine arguments
       procedure(RBF), pointer :: RBFfunction
    
       real(wp), intent(in) :: ValuesIn(:,:), &
                               PointsIn(:,:)
       real(wp), intent(in) :: rho
       
       real(wp), intent(out) :: C(:,:)
       real(wp), intent(out) :: h(:,:)
       !subroutine arguments
       
       !local scalars
       integer :: i, j
       
       !end of header-----------------------------------
       
       h(1,:) = 0d0
       do i=2,size(h,dim=1)
           h(i,:) = ValuesIn(i-1,:)
       end do
       
       C(1,:) = 1d0
       C(:,1) = 1d0
       C(1,1) = 0d0
       do i=2,size(C,dim=1)
           do j=2,size(C,dim=2)
               C(i,j) = RBFfunction(PointsIn(i-1,:), PointsIn(j-1,:), rho)
           end do
       end do
       
    end subroutine CoeffComputation
    !**************************************************************************
                             
    !**************************************************************************
    !EvaluationNodes creates array A
    !**************************************************************************
    subroutine EvaluationMatrix(A, RBFfunction, PointsOut, PointsIn, rho)
    
       !subroutine arguments
       procedure(RBF), pointer :: RBFfunction
       
       real(wp), intent(in) :: PointsOut(:,:), PointsIn(:,:)
       real(wp), intent(in) :: rho
       
       real(wp), intent(out) :: A(:,:)
       !end subroutine arguments
       
       !local scalars
       integer :: i, j
       
       !end of header---------------------------------------
       
       A(:,1) = 1d0
       do i=1,size(A,dim=1)
           do j=2,size(A,dim=2)
               
               A(i,j) = RBFfunction(PointsOut(i,:), PointsIn(j-1,:), rho)
               
           end do
       end do
       
    end subroutine EvaluationMatrix
    !**************************************************************************************

end module RBFinterpolation