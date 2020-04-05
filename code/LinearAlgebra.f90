module LinearAlgebra

    use precision, only: &
    ! Imported Parameters:
    wp

    use utilities, only: &
    ! Imported Procedures:
    reporterror

    use lapack_solver

    implicit none

    ! Local Parameters:
    character(*), parameter, private :: &
    MODname = 'LinearAlgebra'          ! Local module label

    public :: solve_symmetric_system, invert_symmetric_matrix
    public :: condition_number_symmetric_matrix, get_conditionnumber

    private

    interface solve_symmetric_system
        module procedure simple_ssysv
        module procedure simple_dsysv
    end interface

    interface invert_symmetric_matrix
        module procedure ssysvInvert
        module procedure dsysvInvert
    end interface

    interface condition_number_symmetric_matrix
        module procedure simple_dsycon
    end interface

    private :: direct_method
    private :: inverse_method

contains

! lapack_sysv: SYmmetric SolVe system of linear equations A*X=B
subroutine simple_ssysv(A, B, INFO)
    use precision, only : wp => rs
    real(wp), intent(inout) :: A(:,:)
    real(wp), intent(inout) :: B(:,:)
    integer, intent(out)    :: INFO

    integer :: ipiv(size(A,1)), worksize
    real(wp), allocatable :: work(:)
    real(wp) :: worksize_arr(1)

    call lapack_sysv('U', size(A,1), size(B,2), A, max(1,size(A,1)), &
                    ipiv, B, max(1,size(B,1)), worksize_arr, -1, INFO)
    if (info == 0) then
        worksize = int(worksize_arr(1))
        allocate(work(worksize), stat=info)
        if (info == 0) then
            call lapack_sysv('U', size(A,1), size(B,2), A, max(1,size(A,1)), &
                            ipiv, B, max(1,size(B,1)), work, size(work), INFO)
            deallocate(work)
        else
            info = -100
        end if
    end if
end subroutine

subroutine simple_dsysv(A, B, INFO)
    use precision, only : wp => rd
    real(wp), intent(inout) :: A(:,:)
    real(wp), intent(inout) :: B(:,:)
    integer, intent(out)    :: INFO

    integer :: ipiv(size(A,1)), worksize
    real(wp), allocatable :: work(:)
    real(wp) :: worksize_arr(1)

    call lapack_sysv('U', size(A,1), size(B,2), A, max(1,size(A,1)), &
                    ipiv, B, max(1,size(B,1)), worksize_arr, -1, INFO)
    if (info == 0) then
        worksize = int(worksize_arr(1))
        allocate(work(worksize), stat=info)
        if (info == 0) then
            call lapack_sysv('U', size(A,1), size(B,2), A, max(1,size(A,1)), &
                            ipiv, B, max(1,size(B,1)), work, size(work), INFO)
            deallocate(work)
        else
            info = -100
        end if
    end if

end subroutine

subroutine ssysvInvert(A, INFO)
    use precision, only : wp => rs
    real(wp), intent(inout) :: A(:,:)
    integer, intent(out)    :: INFO

    real(wp) :: B(size(A,1),size(A,2))
    integer :: i

    B = A
    A = 0.0
    forall(i=1:size(A,1)) A(i,i) = 1.0
    call solve_symmetric_system(B, A, info)
end subroutine

subroutine dsysvInvert(A, INFO)
    use precision, only : wp => rd
    real(wp), intent(inout) :: A(:,:)
    integer, intent(out)    :: INFO

    real(wp) :: B(size(A,1),size(A,2))
    integer :: i

    B = A
    A = 0.0
    forall(i=1:size(A,1)) A(i,i) = 1.0
    call solve_symmetric_system(B, A, info)
end subroutine

subroutine simple_dsycon(A, RCOND, INFO)
    use precision, only : wp => rd
    real(wp), intent(inout) :: A(:,:)
    real(wp), intent(out)   :: RCOND
    integer, intent(out)    :: INFO

    integer  :: ipiv(size(A,1))
    integer  :: iwork(size(A,1))
    real(wp) :: WORK_size(size(A,1))
    real(wp) :: work_sycon(2*size(A,1))
    real(wp) :: nwork(1)
    real(wp) :: ANORM
    integer :: Lwork

    real(wp), allocatable :: work_sytrf(:)

    ! ANORM = || A ||_1
    ANORM = lapack_lange('1', size(A,1), size(A,2), A, max(1,size(A,1)), nwork)

    ! Calcula Lwork
    call lapack_dsytrf('U', size(A,1), A, max(1,size(A,1)), IPIV, WORK_size, LWORK=-1, INFO=INFO)

    Lwork = int(work_size(1))
    allocate(WORK_sytrf(Lwork))

    ! Factoriza
    call lapack_dsytrf('U', size(A,1), A, max(1,size(A,1)), IPIV, WORK_sytrf, LWORK, INFO)

    deallocate(WORK_sytrf)
	
    call lapack_sycon('U', size(A,1), A, max(1,size(A,1)), ipiv, &
                      ANORM, RCOND, work_sycon, iwork, INFO )

end subroutine

!******************************************************************************
! function condition_number
!******************************************************************************
!
! Description:
!   calcula el condition_number de una matriz A simetrica
!
! Method:
!   cn = lambda_n/lambda_1
!------------------------------------------------------------------------------
function get_conditionnumber(A) result(cn)

    !*Function arguments
    real(wp), intent(in) :: A(:,:)
    real(wp)             :: cn
    !*End Function arguments

    ! Local scacalrs:
    real(wp) :: p
    real(wp) :: lambda_max
    real(wp) :: lambda_min
    integer  :: stat

    ! Local arrays:
    real(wp) :: u(size(A,1))

!- End of header --------------------------------------------------------------

    u = 1.0_wp

    call direct_method(A, u, lambda_max, stat)

    u = 1.0_wp
    p = 0.000010*abs(lambda_max)
    call inverse_method(A, p, u, lambda_min, stat)

    cn = abs(lambda_max/lambda_min)

end function get_conditionnumber

!******************************************************************************
! direct_method
!******************************************************************************
!
! Description:
!   calcula el autovalor mas grande de A
!
!------------------------------------------------------------------------------
subroutine direct_method(A, u, sigma, stat)

    !*Subroutine arguments
    real(wp), intent(in)    :: A(:,:)
    real(wp), intent(inout) :: u(:)
    real(wp), intent(out)   :: sigma
    integer, intent(out)    :: stat
    !*End Subroutine arguments

    ! Local scalars
    real(wp) :: sigma_km1
    real(wp) :: sigma_km2
    real(wp) :: C
    real(wp) :: criterio
    real(wp) :: eps

    integer :: k
    integer :: nMaxIter

    ! Local arrays:
    real(wp) :: v(size(u))

!- End of header --------------------------------------------------------------

    nMaxIter = 1000
    eps = 1e-6_wp

    stat = 0
    k=0
        v = matmul(A,u)
        sigma = dot_product(v,u)/dot_product(u,u)
        u = v/maxval(abs(v))
        sigma_km2 = sigma

    k=1
        v = matmul(A,u)
        sigma = dot_product(v,u)/dot_product(u,u)
        u = v/maxval(abs(v))
        sigma_km1 = sigma

    do k=2,nMaxIter
        v = matmul(A,u)
        sigma = dot_product(v,u)/dot_product(u,u)
        u = v/maxval(abs(v))

        C = abs((sigma - sigma_km1)/(sigma_km1 - sigma_km2))
        criterio = abs((sigma - sigma_km1)/sigma)*C/abs(1d0-C)

        if (criterio <= eps) exit

        sigma_km2 = sigma_km1
        sigma_km1 = sigma

    end do

    if (k > nMaxIter) stat = 1

end subroutine direct_method

!******************************************************************************
! inverse_method
!******************************************************************************
!
! Description:
!   calcula el autovalor mas pequenho de A
!
!------------------------------------------------------------------------------
subroutine inverse_method(A, p, u, sigma, stat)

    !*Subroutine arguments
    real(wp), intent(in)    :: A(:,:)
    real(wp), intent(in)    :: p
    real(wp), intent(inout) :: u(:)
    real(wp), intent(out)   :: sigma
    integer, intent(out)    :: stat
    !*End Subroutine arguments

    ! Local scalars
    real(wp) :: sigma_km1
    real(wp) :: sigma_km2
    real(wp) :: C
    real(wp) :: criterio
    real(wp) :: eps

    integer :: n
    integer :: k
    integer :: nMaxIter
    integer :: info

    ! Local arrays:
    real(wp), allocatable :: Ap(:,:)
    real(wp), allocatable :: Id(:,:)
    real(wp), allocatable :: Av(:,:)
    real(wp), allocatable :: v(:)

!- End of header --------------------------------------------------------------

    nMaxIter = 1000
    eps = 1e-6_wp
    stat = 0

    n = size(A,1)

    ! alocata
    allocate(Ap(n,n))
    allocate(Id(n,n))
    allocate(Av(n,1))
    allocate(v(n))

    ! matriz p*I
    Id = 0.0_wp
    do k=1,n
        Id(k,k) = p
    end do

    ! matriz Ap
    Ap = A - Id

    k=0
        Av(:,1) = u
        call solve_symmetric_system(Ap, Av, info)      ! lapack casera
        if (info /= 0) then
            if (info < 0) then
                ! If info = -i, the ith argument had an illegal value
                call reporterror(MODname//': inverse_method reports on behalf &
                                &of LA: The ith argument had an illegal &
                                &value.')
            else if (info > 0) then
                ! If info = i, U(i,i) is exactly zero. The factorization
                ! has been completed, but the factor U is exactly singular,
                ! so the solution could not be computed.
               call reporterror(MODname//': inverse_method reports on behalf &
                                &of LA: The factor U(i,i) is exactly &
                                &zero.')
            end if
        end if

        v = Av(:,1)
        sigma = dot_product(v,u)/dot_product(u,u)
        u = v/maxval(abs(v))
        sigma_km2 = sigma

    ! matriz Ap
    Ap = A - Id

    k=1
        Av(:,1) = u
        call solve_symmetric_system(Ap, Av, info)      ! lapack casera
        if (info /= 0) then
            if (info < 0) then
                ! If info = -i, the ith argument had an illegal value
                call reporterror(MODname//': inverse_method reports on behalf &
                                &of LA: The ith argument had an illegal &
                                &value.')
            else if (info > 0) then
                ! If info = i, U(i,i) is exactly zero. The factorization
                ! has been completed, but the factor U is exactly singular,
                ! so the solution could not be computed.
               call reporterror(MODname//': inverse_method reports on behalf &
                                &of LA: The factor U(i,i) is exactly &
                                &zero.')
            end if
        end if

        v = Av(:,1)
        sigma = dot_product(v,u)/dot_product(u,u)
        u = v/maxval(abs(v))
        sigma_km1 = sigma

    do k=2,nMaxIter

        ! matriz Ap
        Ap = A - Id

        Av(:,1) = u
        call solve_symmetric_system(Ap, Av, info)      ! lapack casera
        if (info /= 0) then
            if (info < 0) then
                ! If info = -i, the ith argument had an illegal value
                call reporterror(MODname//': inverse_method reports on behalf &
                                &of LA: The ith argument had an illegal &
                                &value.')
            else if (info > 0) then
                ! If info = i, U(i,i) is exactly zero. The factorization
                ! has been completed, but the factor U is exactly singular,
                ! so the solution could not be computed.
               call reporterror(MODname//': inverse_method reports on behalf &
                                &of LA: The factor U(i,i) is exactly &
                                &zero.')
            end if
        end if

        v = Av(:,1)
        sigma = dot_product(v,u)/dot_product(u,u)
        u = v/maxval(abs(v))

        C = abs((sigma - sigma_km1)/(sigma_km1 - sigma_km2))
        criterio = abs((sigma - sigma_km1)/sigma)*C/abs(1d0-C)

        if (criterio <= eps) exit

        sigma_km2 = sigma_km1
        sigma_km1 = sigma

    end do

    sigma = p + 1.0_wp/sigma

    if (k > nMaxIter) stat = 1

end subroutine inverse_method

end module LinearAlgebra
