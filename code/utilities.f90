!******************************************************************************
! utilities utility
!******************************************************************************
!
! Description:
!   This module comprises the implementation of a group of commonly used
!   procedures.
!------------------------------------------------------------------------------
module utilities

    use precision, only: &
        wp

    implicit none

    ! Private Parameters:
    character(*), parameter, private :: MODname ='utilities'

    ! Public Declarations:
    public::              		    &                  
        reporterror,      		    & 
        freeunit,         		    &
        AreEqual,         		    &
        AreSimilar,                 &
        det2,             		    &
        det3,             		    &
        det,              		    &
        vector_product,   		    &
        norm,             		    &
        uppercase,        		    &
        lowercase,          	    &
        isdigit,          		    &
        isnumbermark,     		    &
        removeblanks,      	        &
        strcat,                     &
        numberofreals,              &
		remove_selectedblanks,      &
        exchange_chars,             &
		isalpha,                    &
        MatrixProduct,              &
        InvertMatrix

    interface det
        module procedure det23
    end interface

    private :: det23
    
contains

!******************************************************************************
! Report error
!******************************************************************************
!
! Description:
!   It reports an error string and stops the execution of the program.
!------------------------------------------------------------------------------
subroutine reporterror(error)

    !*Subroutine arguments
    character(len=*), intent(in) :: error
    !*End Subroutine arguments

!- End of header --------------------------------------------------------------

    write(*, *) 'ERROR: ', error
    stop 'Program terminated by reporterror.'

end subroutine reporterror

!******************************************************************************
! Free unit
!******************************************************************************
!
! Description:
!   It generates a free unit between 10 and 100.
!------------------------------------------------------------------------------
function freeunit() result(unit)

    !*Function arguments
    integer :: unit
    !*End Function arguments

    ! Local scalars:
    integer :: i
    logical :: logic
	integer,save :: init = 10

!- End of header --------------------------------------------------------------

    do i=init,100
        inquire(i, opened=logic)
        if (.not.logic) exit
    end do
    unit = i
	init = i + 1

end function freeunit

!******************************************************************************
! AreEqual
!******************************************************************************
!
! Description:
!   It returns true if a = b (taking into account roundoff errors).
!------------------------------------------------------------------------------
function AreEqual(a,b) result(equal)

    !*Function arguments     
    real(wp), intent(in) :: a
    real(wp), intent(in) :: b
    
    ! Result:
    logical :: equal
    !*End Function arguments

    ! Local scalars:
    real(wp) :: eps

!- End of header --------------------------------------------------------------

    eps = epsilon(1.0_wp)
    
    equal = .false.
    if ((a == 0.0_wp) .or. (b == 0.0_wp)) then
    
        if (abs(a-b) < eps) equal = .true.
    
    else
    
        if ((abs(a-b) < eps*abs(a)) .or. (abs(a-b) < eps*abs(b))) equal = .true.
        
    end if
        
end function AreEqual

!******************************************************************************
! AreSimilar
!******************************************************************************
!
! Description:
!   It returns true if a = b (taking into account roundoff errors).
!------------------------------------------------------------------------------
function AreSimilar(a,b) result(equal)

    !*Function arguments     
    real(wp), intent(in) :: a
    real(wp), intent(in) :: b
    
    ! Result:
    logical :: equal
    !*End Function arguments

    ! Local scalars:
    real(wp) :: eps

!- End of header --------------------------------------------------------------

    eps = 1E-6
    
    equal = .false.
    if (abs(a-b) < eps) equal = .true.
        
end function AreSimilar

!******************************************************************************
! det2
!******************************************************************************
!
! Description:
!   Calculates a 2x2d vectors determinant 
!------------------------------------------------------------------------------
function det2(u, v) result(det)

    !* Function arguments
    real(wp), intent(in) :: u(:)
    real(wp), intent(in) :: v(:)
    real(wp) :: det  
    !* End Function arguments

!- End of header ------------------------------------------------------

    det = u(1)*v(2) - u(2)*v(1)

end function det2

!******************************************************************************
! det3
!******************************************************************************
!
! Description:
!   Calculates a 3x3d vectors determinant 
!------------------------------------------------------------------------------
function det3(u, v, w) result(det)

    !* Function arguments
    real(wp), intent(in) :: u(:)
    real(wp), intent(in) :: v(:)
    real(wp), intent(in) :: w(:)
    real(wp) :: det   
    !* End Function arguments

!- End of header ------------------------------------------------------

    det = u(1)*det2(v(2:3), w(2:3))
    det = det - u(2)*det2(v(1:3:2), w(1:3:2))
    det = det + u(3)*det2(v(1:2), w(1:2))

end function det3

!******************************************************************************
! det23
!******************************************************************************
!
! Description:
!   Calculates a 2x2 or 3x3 matrix determinant 
!------------------------------------------------------------------------------
function det23(A) result(det)

    real(wp), intent(in) :: A(:,:)

    ! Result:
    real(wp) :: det

!- End of header --------------------------------------------------------------

    select case (size(A,1))    
    case(2)    
        det = det2(A(1,:), A(2,:)) 
    case(3)
        det = det3(A(1,:), A(2,:), A(3,:))
    end select
    
end function det23

!******************************************************************************
! vector_product
!******************************************************************************
!
! Description:
!   Calculate a 3D vector_product(u,v)
!------------------------------------------------------------------------------
pure funCTIOn vector_product(u, v) result(w)

    real(wp), intent(in) :: u(3)
    real(wp), intent(in) :: v(3)
    real(wp)             :: w(3)

!- End of header --------------------------------------------------------------

    w = [u(2)*v(3) - u(3)*v(2),   &
         u(3)*v(1) - u(1)*v(3),   &
         u(1)*v(2) - u(2)*v(1)]

end function vector_product

!******************************************************************************
! norm
!******************************************************************************
!
! Description:
!   calcula la norma de un vector
!------------------------------------------------------------------------------
function norm(v)

    !* Function arguments
    real(wp), intent(in) :: v(:)
    real(wp) :: norm
    !* End Function arguments

!- End of header --------------------------------------------------------------

    norm = sqrt(dot_product(v,v))

end function norm

!******************************************************************************
! Upper case
!******************************************************************************
!
! Description:
!   It returns the input string in upper case.
!------------------------------------------------------------------------------
pure function uppercase(string) result(newstring)

    !*Function arguments
    character(*), intent(in)   :: string

    ! Result:
    character(len=len(string)) :: newstring
    !*End Function arguments

    ! Local scalars:
    integer :: i

!- End of header --------------------------------------------------------------

    newstring = string

    do i = 1, len_trim(string)
        select case(string(i:i))
        case("a":"z")
            newstring(i:i) = char(ichar(string(i:i))-32)
        end select
    end do

end function uppercase

!******************************************************************************
! lowercase
!******************************************************************************
!
! Description:
!   It returns the input string in lower case.
!------------------------------------------------------------------------------
function lowercase(string) result(newstring)

    !*Function arguments     
    character(*), intent(in)   :: string

    !Result:
    character(len=len(string)) :: newstring
    !*End Function arguments     

    !Local scalars:
    integer :: i

!- End of header --------------------------------------------------------------

    newstring = string

    do i = 1, len_trim(string)
        select case(string(i:i))
        case("A":"Z")
            newstring(i:i) = char(ichar(string(i:i))+32)
        end select
    end do

end function lowercase

!******************************************************************************
! isdigit
!******************************************************************************
!
! Description:
!   It returns true if the character char is a digit and false if not.
!------------------------------------------------------------------------------
function isdigit(char) result(status)

    !*Function arguments     
    character(len=1), intent(in) :: char

    !Result:
    logical :: status
    !*End Function arguments     

!- End of header --------------------------------------------------------------

    status = .false.

    select case(char)
    case("0":"9")
        status = .true.
    end select

end function isdigit

!******************************************************************************
! isalpha
!******************************************************************************
!
! Description:
!   It returns true if the character char is an alpha-numeric character.
!------------------------------------------------------------------------------
function isalpha(char) result(status)

    !*Function arguments     
    character(len=1), intent(in) :: char

    !Result:
    logical :: status
    !*End Function arguments     

!- End of header --------------------------------------------------------------

    status = .false.

    select case(char)
    case("0":"9")
        status = .true.
	case("a":"z", "A":"Z")
		status = .true.
    end select

end function isalpha

!******************************************************************************
! isnumbermark
!******************************************************************************
!
! Description:
!   It returns true if the character char is a numeber mark and false if not.
!------------------------------------------------------------------------------
function isnumbermark(char) result(status)

    !*Function arguments     
    character(len=1), intent(in) :: char

    !Result:
    logical :: status
    !*End Function arguments     

!- End of header --------------------------------------------------------------

    status = .false.

    select case(char)
    case("0":"9", "+", "-", ".")
        status = .true.
    end select

end function isnumbermark

!******************************************************************************
! removeblanks
!******************************************************************************
!
! Description:
!   It returns the same string without blanks.
!------------------------------------------------------------------------------
function removeblanks(string) result(newstring)

    !*Function arguments     
    character(*), intent(in)   :: string

    !Result:
    character(len=len(string)) :: newstring
    !*End Function arguments     

    !Local scalars:
    integer :: i
    integer :: j

!- End of header --------------------------------------------------------------

    newstring = " "

    j = 0
    do i=1,len_trim(string)
        if (string(i:i) /= " ") then
            j = j+1
            newstring(j:j) = string(i:i)
        end if
    end do

end function removeblanks

!******************************************************************************
! remove_selectedblanks
!******************************************************************************
!
! Description:
!   It returns the same string without blanks before and after the 
!	selected character.
!------------------------------------------------------------------------------
function remove_selectedblanks(string, selectedchar) result(newstring)

    !*Function arguments     
    character(*), intent(in ) :: string
    character, intent(in)     :: selectedchar
    character(len=len(string)) :: newstring
    !*End Function arguments     

    !Local scalars:
    integer :: i
    integer :: i1
    integer :: i2
	integer :: lenstring
    character(len=len(string)) :: scratch

!- End of header --------------------------------------------------------------

    newstring = " "
	lenstring = len_trim(string)

	i2 = 0
	do
		i1 = index(string(i2+1:), selectedchar)
		if (i1 == 0) exit

		i1 = i1 + i2
		
		i = i1 + 1
		do
			if ((string(i:i) /= " ") .and. (string(i+1:i+1) == " ")) exit
			i = i + 1
		end do

		scratch = ""
		scratch(1:i-i2) = string(i2+1:i)

		newstring = trim(newstring)//"  "//removeblanks(scratch(1:i-i2))

		i2 = i
	end do

end function remove_selectedblanks

!******************************************************************************
! exchange_chars
!******************************************************************************
!
! Description:
!   It concatenates the stringleft and stringright.
!------------------------------------------------------------------------------
function exchange_chars(scratch, removedchar, newchar) result(newstring)

    !*Function arguments     
    character(*), intent(in) :: scratch
    character(*), intent(in) :: removedchar
    character(*), intent(in) :: newchar

    !Result:
    character(len=len(scratch)) :: newstring
    !*End Function arguments     

    !Local scalars:
    integer :: i

!- End of header --------------------------------------------------------------

    newstring = scratch
    do i=1,len(scratch)
        if (newstring(i:i) == removedchar) newstring(i:i)=newchar
    end do

end function exchange_chars

!******************************************************************************
! strcat
!******************************************************************************
!
! Description:
!   It concatenates the stringleft and stringright.
!------------------------------------------------------------------------------
function strcat(stringleft, stringright) result(newstring)

    !*Function arguments     
    character(*), intent(in) :: stringleft
    character(*), intent(in) :: stringright

    !Result:
    character(len=len_trim(stringleft)+len_trim(stringright)) :: newstring
    !*End Function arguments     

    !Local scalars:
    integer :: len1
    integer :: len2

!- End of header --------------------------------------------------------------

    newstring = " "
    
    len1 = len_trim(stringleft)
    len2 = len_trim(stringright)

    newstring(1:len1) = adjustl(stringleft)
    newstring(len1+1:len1+len2) = adjustl(stringright)

end function strcat

!******************************************************************************
! numberofreals
!******************************************************************************
!
! Description:
!   It returns the number of reals contained in string.
!------------------------------------------------------------------------------
function numberofreals(string) result(number)

    !*Function arguments     
    character(*), intent(in) :: string

    !Result:
    integer :: number
    !*End Function arguments     

    !Local scalars:
    integer :: i1
    integer :: i2

!- End of header --------------------------------------------------------------

    number = 0
    i1 = 1
    i2 = 0
    do
        i2 = index(string(i1:), '.', back=.false.)        
        if (i2 == 0) exit
        i1 = i1 + i2 + 1
        number = number + 1
    end do

end function numberofreals
!************************************************************************************

!************************************************************************************
!MatrixProduct
!
!performs the product of two matrices: A*B=C
function MatrixProduct(A, B) result(C)

   !subroutine arguments
   real(wp), intent(in) :: A(:,:)
   real(wp), intent(in) :: B(:,:)
   
   real(wp) :: C(size(A, dim=1), size(B, dim=2))
   !end subroutine arguments
   
   !local scalars
   integer :: i, j
   
   !end of header--------------------------------
   
   !check if matrices are compatible
   if(size(A, dim=2) /= size(B, dim=1))then
       call reporterror(MODname// ': MatrixProduct reports: &
                   &matrices index do not fit')
       stop
   end if
   
   do i=1,size(A, dim=1)
       do j=1, size(B, dim=2)
           
           C(i,j) = dot_product(A(i,:), B(:,j))
           
       end do
   end do
   
end function MatrixProduct
!************************************************************************************

!************************************************************************************
!InvertMatrix
!
!inverts the matrix A
function InvertMatrix(A) result(Ainv)

   !subroutine arguments
   real(wp), intent(in) :: A(:,:)
   
   real(wp) :: Ainv(size(A, dim=1), size(A, dim=2))
   !end subroutine arguments
   
   !local scalars
   real(wp) :: factor
   real(wp) :: B(size(A, dim=1), size(A, dim=2))
   
   integer :: j, k, m, n
   
   !end of header--------------------------------------
   
   !checks A is square matrix
   if(size(A, dim=1) /= size(A, dim=2))then
       call reporterror(MODname// ': InvertMatrix reports: &
                   &matrices index do not fit')
       stop
   end if
   
   B = A
      
   !identity matrix
   Ainv = 0.0_wp
   do j=1,size(Ainv, dim=1)
       Ainv(j,j) = 1.0_wp
   end do
   
   !calculates lower Ainv
   do j=1,size(A, dim=1)
       do k=j,size(A, dim=1)
           factor = B(k,k)
           Ainv(k,:) = Ainv(k,:)/factor
           B(k,:) = B(k,:)/factor
           do m=k+1,size(A, dim=1)
               factor = B(m,k)
               do n=1, size(A, dim=2)
                   Ainv(m,n) = Ainv(m,n)-factor*Ainv(k,n)
                   B(m,n) = B(m,n)-factor*B(k,n)
               end do
           end do
       end do
   end do
   
   !calculates upper Ainv
   do j=size(A, dim=1),1,-1
       do k=j,1,-1
           factor = B(k,k)
           Ainv(k,:) = Ainv(k,:)/factor
           B(k,:) = B(k,:)/factor
           do m=k-1,1,-1
               factor = B(m,k)
               do n=size(A, dim=2),1,-1
                   Ainv(m,n) = Ainv(m,n)-factor*Ainv(k,n)
                   B(m,n) = B(m,n)-factor*B(k,n)
               end do
           end do
       end do
   end do
   
end function InvertMatrix
!***********************************************************************************

end module utilities
