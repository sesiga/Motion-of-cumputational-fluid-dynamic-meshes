!*********************************************************
!RBFfunctions
!*********************************************************
!
!This module contains radial basis funtions
!to interpolate

module RadialBasisFunctions

   !modules used
   use precision, only : wp
   
   use utilities, only : reporterror
   !end modules used 
   
   implicit none
   
   !private parameters
   character(*), parameter, private :: MODname = 'RadialBasisFunctions'
   
   !procedure to select function
   interface
       pure function RBF(x1, x2, rho) result(phi)
           use precision
           real(wp), intent(in) :: x1(:), x2(:)
           real(wp), intent(in) :: rho
           real(wp) :: phi
       end function
   end interface
   !end procedure 
   
    contains
    !end of header--------------------------
    
    !***************************************************************
    !Spline
    !***************************************************************
    pure function spline(x1, x2, rho) result(phi)
    
        !function arguments
        real(wp), intent(in) :: x1(:), x2(:)
        real(wp), intent(in) :: rho
        
        real(wp) :: phi
        !end function arguments
        
        phi = norm2(x1-x2)
        
    end function spline
    !**************************************************************
    
    !***************************************************************
    !Wendland C0
    !***************************************************************
    pure function WendlandC0(x1, x2, rho) result(phi)
    
        !function arguments
        real(wp), intent(in) :: x1(:), x2(:)
        real(wp), intent(in) :: rho
        
        real(wp) :: phi
        !end function arguments
        
        if( norm2(x1-x2)/rho < 1d0 )then
            phi = (1d0-norm2(x1-x2)/rho)**2
        else
            phi = 0d0
        end if
                
    end function WendlandC0
    !**************************************************************
    
    !***************************************************************
    !Wendland C2
    !***************************************************************
    pure function WendlandC2(x1, x2, rho) result(phi)
    
        !function arguments
        real(wp), intent(in) :: x1(:), x2(:)
        real(wp), intent(in) :: rho
        
        real(wp) :: phi
        !end function arguments
        
        if( norm2(x1-x2)/rho < 1d0 )then
            phi = (1d0-norm2(x1-x2)/rho)**4*(4d0*norm2(x1-x2)/rho+1d0)
        else
            phi = 0d0
        end if
                
    end function WendlandC2
    !**************************************************************
    
    !***************************************************************
    !Wendland C4
    !***************************************************************
    pure function WendlandC4(x1, x2, rho) result(phi)
    
        !function arguments
        real(wp), intent(in) :: x1(:), x2(:)
        real(wp), intent(in) :: rho
        
        real(wp) :: phi
        !end function arguments
        
        if( norm2(x1-x2)/rho < 1d0 )then
            phi = (1d0-norm2(x1-x2)/rho)**6*(35d0*(norm2(x1-x2))**2+18d0*norm2(x1-x2)+3d0)
        else
            phi = 0d0
        end if
                
    end function WendlandC4
    !**************************************************************
    
    !**************************************************************
    !hardy
    !
    pure function hardy(x1, x2, rho) result(phi)
    
       !function arguments
       real(wp), intent(in) :: x1(:), x2(:)
       real(wp), intent(in) :: rho
        
       real(wp) :: phi
       !end function arguments
       
       phi = sqrt(rho**2+(norm2(x1-x2))**2)
       
    end function hardy
   
end module RadialBasisFunctions