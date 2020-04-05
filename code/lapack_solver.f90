module lapack_solver
      
      public :: lapack_sysv, lapack_sycon, lapack_dsytrf, lapack_lange
           
      interface lapack_sysv 
          module procedure SSYSV
          module procedure DSYSV
      end interface lapack_sysv 

      interface lapack_dsytrf
        module procedure DSYTRF
      end interface lapack_dsytrf

      interface lapack_sycon
        module procedure DSYCON
      end interface lapack_sycon

      interface lapack_lange
        module procedure DLANGE
      end interface lapack_lange

private
      
contains
                             
!----------------------------------------------------------------------
!   lapack_dsysv_full
!----------------------------------------------------------------------
subroutine DISNANs( DIN,  DISNAN)
!  
!     .. Scalar Arguments ..
      real(8), intent(in)  :: DIN
      logical, intent(out) :: DISNAN
!     ..
!  
!  Purpose
!  =======
!  
!  DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
!  otherwise.  To be replaced by the Fortran 2003 intrinsic in the
!  future.
!  
!  .. Executable Statements ..
      DISNAN = DLAISNAN(DIN,DIN)

end subroutine DISNANs
 
logical function DLAISNAN( DIN1, DIN2 )
!  
!     .. Scalar Arguments ..
      real(8) :: DIN1, DIN2
!     ..
!  
!  
!  .. Executable Statements ..
      DLAISNAN = (DIN1.NE.DIN2)

end function DLAISNAN
      
subroutine DLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )
!  
!     .. Scalar Arguments ..
      character          UPLO
      integer :: INFO, KB, LDA, LDW, N, NB
!     ..
!     .. Array Arguments ..
      integer ::   IPIV( * )
      real(8) ::   A( LDA, * ), W( LDW, * )

!     .. Parameters ..
      real(8), parameter :: ZERO = 0.0D+0
      real(8), parameter :: ONE = 1.0D+0
      real(8), parameter :: EIGHT = 8.0D+0
      real(8), parameter :: SEVTEN = 17.0D+0
!     ..
!     .. Local Scalars ..
      integer :: IMAX, J, JB, JJ, JMAX, JP, K, KK, KKW, KP
      integer :: KSTEP, KW
      real(8) :: ABSAKK, ALPHA, COLMAX, D11, D21, D22, R1
      real(8) :: ROWMAX, T
!     ..
!     .. External Functions ..
!m      logical ::            LSAME
!m      integer ::            IDAMAX
!m      EXTERNAL           LSAME, IDAMAX
!     ..
!     .. External Subroutines ..
!m      EXTERNAL           DCOPY, DGEMM, DGEMV, DSCAL, DSWAP
!     ..
!     .. Intrinsic Functions ..
!m      INTRINSIC          abs, max, min, sqrt
!     ..
!     .. Executable Statements ..
!  
      INFO = 0
!  
!     Initialize ALPHA for use in choosing pivot block size.
!  
      ALPHA = (ONE + sqrt(SEVTEN))/EIGHT
!  
      if ( LSAME( UPLO, 'U' ) ) then
!  
!        Factorize the trailing columns of A using the upper triangle
!        of A and working backwards, and compute the matrix W = U12*D
!        for use in updating A11
!  
!        K is the main loop index, decreasing from N in steps of 1 or 2
!  
!        KW is the column of W which corresponds to column K of A
!  
         K = N
   10    CONTINUE
         KW = NB + K - N
!  
!        Exit from loop
!  
         if ((K <= N-NB+1 .and. NB < N) .or. (K < 1))  GO TO 30
!  
!        Copy column K of A to column KW of W and update it
!  
         call DCOPY( K, A( 1, K ), 1, W( 1, KW ), 1 ) 
         if (K < N)                                                     &
           call DGEMV( 'No transpose', K, N-K, -ONE, A( 1, K+1 ), LDA,  &
                        W( K, KW+1 ), LDW, ONE, W( 1, KW ), 1 )
!  
         KSTEP = 1
!  
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!  
         ABSAKK = abs(W(K,KW) )
!  
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!  
         if (K > 1) then
            IMAX = IDAMAX(K-1, W(1,KW), 1 )
            COLMAX = abs(W(IMAX,KW))
         else
            COLMAX = ZERO
         end if
!  
         if (max(ABSAKK,COLMAX ) == ZERO) then
!  
!           Column K is zero: set INFO and continue
!  
            if ( INFO == 0 ) INFO = K
            KP = K
         else
            if ( ABSAKK >= ALPHA*COLMAX ) then
!  
!              no interchange, use 1-by-1 pivot block
!  
               KP = K
            else
!  
!              Copy column IMAX to column KW-1 of W and update it
!  
               call DCOPY( IMAX, A( 1, IMAX ), 1, W( 1, KW-1 ), 1 )
               call DCOPY( K-IMAX, A( IMAX, IMAX+1 ), LDA,                &
                          W( IMAX+1, KW-1 ), 1 ) 
               if (K < N)                                               &
                  call DGEMV( 'No transpose', K, N-K, -ONE, A( 1, K+1 ),  &
                              LDA, W( IMAX, KW+1 ), LDW, ONE,             &
                              W( 1, KW-1 ), 1 )
!  
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!  
               JMAX = IMAX + IDAMAX( K-IMAX, W( IMAX+1, KW-1 ), 1 )
               ROWMAX = abs( W( JMAX, KW-1 ) )
               if (IMAX > 1) then
                  JMAX = IDAMAX( IMAX-1, W( 1, KW-1 ), 1 )
                  ROWMAX = max( ROWMAX, abs( W( JMAX, KW-1 ) ) )
               end if
!  
               if ( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) then
!  
!                 no interchange, use 1-by-1 pivot block
!  
                  KP = K
               else if ( abs( W( IMAX, KW-1 ) ) >= ALPHA*ROWMAX ) then
!  
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!  
                  KP = IMAX
!  
!                 copy column KW-1 of W to column KW
!  
                  call DCOPY( K, W( 1, KW-1 ), 1, W( 1, KW ), 1 )
               else
!  
!                 interchange rows and columns K-1 and IMAX, use 2-by-2
!                 pivot block
!  
                  KP = IMAX
                  KSTEP = 2
               end if
            end if
!  
            KK = K - KSTEP + 1
            KKW = NB + KK - N
!  
!           Updated column KP is already stored in column KKW of W
!  
            if ( KP.NE.KK ) then
!  
!              Copy non-updated column KK to column KP
!  
               A( KP, K ) = A( KK, K )
               call DCOPY( K-1-KP, A( KP+1, KK ), 1, A( KP, KP+1 ), LDA )
               call DCOPY( KP, A( 1, KK ), 1, A( 1, KP ), 1 )
!  
!              Interchange rows KK and KP in last KK columns of A and W
!  
               call DSWAP( N-KK+1, A( KK, KK ), LDA, A( KP, KK ), LDA )
               call DSWAP( N-KK+1, W( KK, KKW ), LDW, W( KP, KKW ), LDW )
            end if
!  
            if ( KSTEP == 1 ) then
!  
!              1-by-1 pivot block D(k): column KW of W now holds
!  
!              W(k) = U(k)*D(k)
!  
!              where U(k) is the k-th column of U
!  
!              Store U(k) in column k of A
!  
               call DCOPY( K, W( 1, KW ), 1, A( 1, K ), 1 )
               R1 = ONE / A( K, K )
               call DSCAL( K-1, R1, A( 1, K ), 1 )
            else
!  
!              2-by-2 pivot block D(k): columns KW and KW-1 of W now
!              hold
!  
!              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
!  
!              where U(k) and U(k-1) are the k-th and (k-1)-th columns
!              of U
!  
               if ( K > 2 ) then
!  
!                 Store U(k) and U(k-1) in columns k and k-1 of A
!  
                  D21 = W( K-1, KW )
                  D11 = W( K, KW ) / D21
                  D22 = W( K-1, KW-1 ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
                  DO 20 J = 1, K - 2
                     A( J, K-1 ) = D21*( D11*W( J, KW-1 )-W( J, KW ) )
                     A( J, K ) = D21*( D22*W( J, KW )-W( J, KW-1 ) )
   20             CONTINUE
               end if
!  
!              Copy D(k) to A
!  
               A( K-1, K-1 ) = W( K-1, KW-1 )
               A( K-1, K ) = W( K-1, KW )
               A( K, K ) = W( K, KW )
            end if
         end if
!  
!        Store details of the interchanges in IPIV
!  
         if ( KSTEP == 1 ) then
            IPIV( K ) = KP
         else
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         end if
!  
!        Decrease K and return to the start of the main loop
!  
         K = K - KSTEP
         GO TO 10
!  
   30    CONTINUE
!  
!        Update the upper triangle of A11 (= A(1:k,1:k)) as
!  
!        A11 := A11 - U12*D*U12' = A11 - U12*W'
!  
!        computing blocks of NB columns at a time
!  
         DO 50 J = ( ( K-1 ) / NB )*NB + 1, 1, -NB
            JB = min( NB, K-J+1 )
!  
!           Update the upper triangle of the diagonal block
!  
            DO 40 JJ = J, J + JB - 1
               call DGEMV( 'No transpose', JJ-J+1, N-K, -ONE,          & 
                           A( J, K+1 ), LDA, W( JJ, KW+1 ), LDW, ONE,  &
                           A( J, JJ ), 1 )
   40       CONTINUE
!  
!           Update the rectangular superdiagonal block
!  
            call DGEMM( 'No transpose', 'Transpose', J-1, JB, N-K, -ONE, &
                        A( 1, K+1 ), LDA, W( J, KW+1 ), LDW, ONE,        &
                        A( 1, J ), LDA )
   50    CONTINUE
!  
!        Put U12 in standard form by partially undoing the interchanges
!        in columns k+1:n
!  
         J = K + 1
   60    CONTINUE
         JJ = J
         JP = IPIV( J )
         if ( JP < 0 ) then
            JP = -JP
            J = J + 1
         end if
         J = J + 1
         if ( JP.NE.JJ .and. J <= N )  &
            call DSWAP( N-J+1, A( JP, J ), LDA, A( JJ, J ), LDA )
         if ( J <= N ) GO TO 60
!  
!        Set KB to the number of columns factorized
!  
         KB = N - K
!  
      else
!  
!        Factorize the leading columns of A using the lower triangle
!        of A and working forwards, and compute the matrix W = L21*D
!        for use in updating A22
!  
!        K is the main loop index, increasing from 1 in steps of 1 or 2
!  
         K = 1
   70    CONTINUE
!  
!        Exit from loop
!  
         if ( ( K >= NB .and. NB < N ) .or. K > N ) GO TO 90
!  
!        Copy column K of A to column K of W and update it
!  
         call DCOPY( N-K+1, A( K, K ), 1, W( K, K ), 1 )
         call DGEMV( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ), LDA,  &
                     W( K, 1 ), LDW, ONE, W( K, K ), 1 )
!  
         KSTEP = 1
!  
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!  
         ABSAKK = abs( W( K, K ) )
!  
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!  
         if ( K < N ) then
            IMAX = K + IDAMAX( N-K, W( K+1, K ), 1 )
            COLMAX = abs( W( IMAX, K ) )
         else
            COLMAX = ZERO
         end if
!  
         if ( max( ABSAKK, COLMAX ) == ZERO ) then
!  
!           Column K is zero: set INFO and continue
!  
            if ( INFO == 0 ) INFO = K
            KP = K
         else
            if ( ABSAKK >= ALPHA*COLMAX ) then
!  
!              no interchange, use 1-by-1 pivot block
!  
               KP = K
            else
!  
!              Copy column IMAX to column K+1 of W and update it
!  
               call DCOPY( IMAX-K, A( IMAX, K ), LDA, W( K, K+1 ), 1 )
               call DCOPY( N-IMAX+1, A( IMAX, IMAX ), 1, W( IMAX, K+1 ), 1 )
               call DGEMV( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ),      &
                           LDA, W( IMAX, 1 ), LDW, ONE, W( K, K+1 ), 1 )
!  
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!  
               JMAX = K - 1 + IDAMAX( IMAX-K, W( K, K+1 ), 1 )
               ROWMAX = abs( W( JMAX, K+1 ) )
               if ( IMAX < N ) then
                  JMAX = IMAX + IDAMAX( N-IMAX, W( IMAX+1, K+1 ), 1 )
                  ROWMAX = max( ROWMAX, abs( W( JMAX, K+1 ) ) )
               end if
!  
               if ( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) then
!  
!                 no interchange, use 1-by-1 pivot block
!  
                  KP = K
               else if ( abs( W( IMAX, K+1 ) ) >= ALPHA*ROWMAX ) then
!  
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!  
                  KP = IMAX
!  
!                 copy column K+1 of W to column K
!  
                  call DCOPY( N-K+1, W( K, K+1 ), 1, W( K, K ), 1 )
               else
!  
!                 interchange rows and columns K+1 and IMAX, use 2-by-2
!                 pivot block
!  
                  KP = IMAX
                  KSTEP = 2
               end if
            end if
!  
            KK = K + KSTEP - 1
!  
!           Updated column KP is already stored in column KK of W
!  
            if ( KP.NE.KK ) then
!  
!              Copy non-updated column KK to column KP
!  
               A( KP, K ) = A( KK, K )
               call DCOPY( KP-K-1, A( K+1, KK ), 1, A( KP, K+1 ), LDA )
               call DCOPY( N-KP+1, A( KP, KK ), 1, A( KP, KP ), 1 )
!  
!              Interchange rows KK and KP in first KK columns of A and W
!  
               call DSWAP( KK, A( KK, 1 ), LDA, A( KP, 1 ), LDA )
               call DSWAP( KK, W( KK, 1 ), LDW, W( KP, 1 ), LDW )
            end if
!  
            if ( KSTEP == 1 ) then
!  
!              1-by-1 pivot block D(k): column k of W now holds
!  
!              W(k) = L(k)*D(k)
!  
!              where L(k) is the k-th column of L
!  
!              Store L(k) in column k of A
!  
               call DCOPY( N-K+1, W( K, K ), 1, A( K, K ), 1 )
               if ( K < N ) then
                  R1 = ONE / A( K, K )
                  call DSCAL( N-K, R1, A( K+1, K ), 1 )
               end if
            else
!  
!              2-by-2 pivot block D(k): columns k and k+1 of W now hold
!  
!              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
!  
!              where L(k) and L(k+1) are the k-th and (k+1)-th columns
!              of L
!  
               if ( K < N-1 ) then
!  
!                 Store L(k) and L(k+1) in columns k and k+1 of A
!  
                  D21 = W( K+1, K )
                  D11 = W( K+1, K+1 ) / D21
                  D22 = W( K, K ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
                  DO 80 J = K + 2, N
                     A( J, K ) = D21*( D11*W( J, K )-W( J, K+1 ) )
                     A( J, K+1 ) = D21*( D22*W( J, K+1 )-W( J, K ) )
   80             CONTINUE
               end if
!  
!              Copy D(k) to A
!  
               A( K, K ) = W( K, K )
               A( K+1, K ) = W( K+1, K )
               A( K+1, K+1 ) = W( K+1, K+1 )
            end if
         end if
!  
!        Store details of the interchanges in IPIV
!  
         if ( KSTEP == 1 ) then
            IPIV( K ) = KP
         else
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         end if
!  
!        Increase K and return to the start of the main loop
!  
         K = K + KSTEP
         GO TO 70
!  
   90    CONTINUE
!  
!        Update the lower triangle of A22 (= A(k:n,k:n)) as
!  
!        A22 := A22 - L21*D*L21' = A22 - L21*W'
!  
!        computing blocks of NB columns at a time
!  
         DO 110 J = K, N, NB
            JB = min( NB, N-J+1 )
!  
!           Update the lower triangle of the diagonal block
!  
            DO 100 JJ = J, J + JB - 1
               call DGEMV( 'No transpose', J+JB-JJ, K-1, -ONE,     &
                           A( JJ, 1 ), LDA, W( JJ, 1 ), LDW, ONE,  &
                           A( JJ, JJ ), 1 )
  100       CONTINUE
!  
!           Update the rectangular subdiagonal block
!  
            if ( J+JB <= N ) &
               call DGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB,     &
                           K-1, -ONE, A( J+JB, 1 ), LDA, W( J, 1 ), LDW,  &
                           ONE, A( J+JB, J ), LDA )
  110    CONTINUE
!  
!        Put L21 in standard form by partially undoing the interchanges
!        in columns 1:k-1
!  
         J = K - 1
  120    CONTINUE
         JJ = J
         JP = IPIV( J )
         if ( JP < 0 ) then
            JP = -JP
            J = J - 1
         end if
         J = J - 1
         if ( JP.NE.JJ .and. J >= 1 )  &
            call DSWAP( J, A( JP, 1 ), LDA, A( JJ, 1 ), LDA )
         if ( J >= 1 ) GO TO 120
!  
!        Set KB to the number of columns factorized
!  
         KB = K - 1
!  
      end if
      RETURN
!  
!     End of DLASYF
!  
      end subroutine DLASYF
   
   
      subroutine DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,   &
                        LWORK, INFO )
!  
!  -- LAPACK driver routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!  
    
!     .. Scalar Arguments ..
      character          UPLO
      integer ::            INFO, LDA, LDB, LWORK, N, NRHS
!     ..
!     .. Array Arguments ..
      integer ::            IPIV( * )
      real(8) ::   A( LDA, * ), B( LDB, * ), WORK( * )
!     ..
!  
!  Purpose
!  =======
!  
!  DSYSV computes the solution to a real system of linear equations
!     A * X = B,
!  where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
!  matrices.
!  
!  The diagonal pivoting method is used to factor A as
!     A = U * D * U**T,  if UPLO = 'U', or
!     A = L * D * L**T,  if UPLO = 'L',
!  where U (or L) is a product of permutation and unit upper (lower)
!  triangular matrices, and D is symmetric and block diagonal with
!  1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
!  used to solve the system of equations A * X = B.
!  
!  Arguments
!  =========
!  
!  UPLO    (input) character*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!  
!  N       (input) integer ::
!          The number of linear equations, i.e., the order of the
!          matrix A.  N >= 0.
!  
!  NRHS    (input) integer ::
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!  
!  A       (input/output) real(8) :: array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!  
!          On exit, if INFO = 0, the block diagonal matrix D and the
!          multipliers used to obtain the factor U or L from the
!          factorization A = U*D*U**T or A = L*D*L**T as computed by
!          DSYTRF.
!  
!  LDA     (input) integer ::
!          The leading dimension of the array A.  LDA >= max(1,N).
!  
!  IPIV    (output) integer :: array, dimension (N)
!          Details of the interchanges and the block structure of D, as
!          determined by DSYTRF.  If IPIV(k) > 0, then rows and columns
!          k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1
!          diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0,
!          then rows and columns k-1 and -IPIV(k) were interchanged and
!          D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and
!          IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and
!          -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2
!          diagonal block.
!  
!  B       (input/output) real(8) :: array, dimension (LDB,NRHS)
!          On entry, the N-by-NRHS right hand side matrix B.
!          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!  
!  LDB     (input) integer ::
!          The leading dimension of the array B.  LDB >= max(1,N).
!  
!  WORK    (workspace/output) real(8) :: array, dimension (max(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!  
!  LWORK   (input) integer ::
!          The length of WORK.  LWORK >= 1, and for best performance
!          LWORK >= max(1,N*NB), where NB is the optimal blocksize for
!          DSYTRF.
!  
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!  
!  INFO    (output) integer ::
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
!               has been completed, but the block diagonal matrix D is
!               exactly singular, so the solution could not be computed.
!  
!  =====================================================================
!  
!     .. Local Scalars ..
      logical ::            LQUERY
      integer ::            LWKOPT, NB
!     ..
!     .. External Functions ..
!m      logical ::            LSAME
!m      integer ::            ILAENV
!m      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
!m      EXTERNAL           DSYTRF, DSYTRS, XERBLA
!     ..
!     .. Intrinsic Functions ..
!m     INTRINSIC          max
!     ..
!     .. Executable Statements ..
!  
!     Test the input parameters.
!  
      INFO = 0
      LQUERY = ( LWORK == -1 )
      if ( .NOT.LSAME( UPLO, 'U' ) .and. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( NRHS < 0 ) then
         INFO = -3
      else if ( LDA < max( 1, N ) ) then
         INFO = -5
      else if ( LDB < max( 1, N ) ) then
         INFO = -8
      else if ( LWORK < 1 .and. .NOT.LQUERY ) then
         INFO = -10
      end if
!  
      if ( INFO == 0 ) then
         if ( N == 0 ) then
            LWKOPT = 1
         else
            NB = ILAENV( 1, 'DSYTRF', UPLO, N, -1, -1, -1 )
            LWKOPT = N*NB
         end if
         WORK( 1 ) = LWKOPT
      end if
!  
      if ( INFO.NE.0 ) then
         call XERBLA( 'DSYSV ', -INFO )
         RETURN
      else if ( LQUERY ) then
         RETURN
      end if
!  
!     Compute the factorization A = U*D*U' or A = L*D*L'.
!  
      call DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
      if ( INFO == 0 ) then
!  
!        Solve the system A*X = B, overwriting B with X.
!  
         call DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!  
      end if
!  
      WORK( 1 ) = LWKOPT
!  
      RETURN
!  
!     End of DSYSV
!  
      end  subroutine DSYSV
      
      
      subroutine DSYTF2( UPLO, N, A, LDA, IPIV, INFO )
!  
!  -- LAPACK routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!  
!     .. Scalar Arguments ..
      character  :: UPLO
      integer    :: INFO, LDA, N
!     ..
!     .. Array Arguments ..
      integer ::   IPIV( * )
      real(8) ::   A( LDA, * )
!     ..
!  
!  Purpose
!  =======
!  
!  DSYTF2 computes the factorization of a real symmetric matrix A using
!  the Bunch-Kaufman diagonal pivoting method:
!  
!     A = U*D*U'  or  A = L*D*L'
!  
!  where U (or L) is a product of permutation and unit upper (lower)
!  triangular matrices, U' is the transpose of U, and D is symmetric and
!  block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!  
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!  
!  Arguments
!  =========
!  
!  UPLO    (input) character*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!  
!  N       (input) integer ::
!          The order of the matrix A.  N >= 0.
!  
!  A       (input/output) real(8) :: array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!  
!          On exit, the block diagonal matrix D and the multipliers used
!          to obtain the factor U or L (see below for further details).
!  
!  LDA     (input) integer ::
!          The leading dimension of the array A.  LDA >= max(1,N).
!  
!  IPIV    (output) integer :: array, dimension (N)
!          Details of the interchanges and the block structure of D.
!          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!          interchanged and D(k,k) is a 1-by-1 diagonal block.
!          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!  
!  INFO    (output) integer ::
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
!               has been completed, but the block diagonal matrix D is
!               exactly singular, and division by zero will occur if it
!               is used to solve a system of equations.
!  
!  Further Details
!  ===============
!  
!  09-29-06 - patch from
!    Bobby Cheng, MathWorks
!  
!    Replace l.204 and l.372
!         if ( max( ABSAKK, COLMAX ) == ZERO ) then
!    by
!         if ( (max( ABSAKK, COLMAX ) == ZERO) .or. DISNAN(ABSAKK) ) then
!  
!  01-01-96 - Based on modifications by
!    J. Lewis, Boeing Computer Services Company
!    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!  1-96 - Based on modifications by J. Lewis, Boeing Computer Services
!         Company
!  
!  If UPLO = 'U', then A = U*D*U', where
!     U = P(n)*U(n)* ... *P(k)U(k)* ...,
!  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!  
!             (   I    v    0   )   k-s
!     U(k) =  (   0    I    0   )   s
!             (   0    0    I   )   n-k
!                k-s   s   n-k
!  
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!  and A(k,k), and v overwrites A(1:k-2,k-1:k).
!  
!  If UPLO = 'L', then A = L*D*L', where
!     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
!  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
!  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!  
!             (   I    0     0   )  k-1
!     L(k) =  (   0    I     0   )  s
!             (   0    v     I   )  n-k-s+1
!                k-1   s  n-k-s+1
!  
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!  
!  =====================================================================
!  
!     .. Parameters ..
      real(8) ::   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      real(8) ::   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
!     ..
!     .. Local Scalars ..
      logical ::            UPPER
      integer ::            I, IMAX, J, JMAX, K, KK, KP, KSTEP
      real(8) ::   ABSAKK, ALPHA, COLMAX, D11, D12, D21, D22, R1
      real(8) ::   ROWMAX, T, WK, WKM1, WKP1
      
      logical :: nuestra
!     ..
!     .. External Functions ..
!m      logical ::            LSAME, DISNAN
!m      integer ::            IDAMAX
!m      EXTERNAL           LSAME, IDAMAX, DISNAN
!     ..
!     .. External Subroutines ..
!m      EXTERNAL           DSCAL, DSWAP, DSYR, XERBLA
!     ..
!     .. Intrinsic Functions ..
!m      INTRINSIC          abs, max, sqrt
!     ..
!     .. Executable Statements ..
!  
!     Test the input parameters.
!          
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .and. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( LDA < max( 1, N ) ) then
         INFO = -4
      end if
      if ( INFO.NE.0 ) then
         call XERBLA( 'DSYTF2', -INFO )
         RETURN
      end if
!  
!     Initialize ALPHA for use in choosing pivot block size.
!  
      ALPHA = ( ONE+sqrt( SEVTEN ) ) / EIGHT
!  
      if ( UPPER ) then
!  
!        Factorize A as U*D*U' using the upper triangle of A
!  
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2
!  
         K = N
   10    CONTINUE
!  
!        If K < 1, exit from loop
!  
         if ( K < 1 ) GO TO 70
         KSTEP = 1
!  
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!  
         ABSAKK = abs( A( K, K ) )
!  
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!  
         if ( K > 1 ) then
            IMAX = IDAMAX( K-1, A( 1, K ), 1 )
            COLMAX = abs( A( IMAX, K ) )
         else
            COLMAX = ZERO
         end if
!  
        call DISNANs(ABSAKK, nuestra)
         if ( (max( ABSAKK, COLMAX ) == ZERO) .or. nuestra ) then
!  
!           Column K is zero or contains a NaN: set INFO and continue
!  
            if ( INFO == 0 ) INFO = K
            KP = K
         else
            if ( ABSAKK >= ALPHA*COLMAX ) then
!  
!              no interchange, use 1-by-1 pivot block
!  
               KP = K
            else
!  
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!  
               JMAX = IMAX + IDAMAX( K-IMAX, A( IMAX, IMAX+1 ), LDA )
               ROWMAX = abs( A( IMAX, JMAX ) )
               if ( IMAX > 1 ) then
                  JMAX = IDAMAX( IMAX-1, A( 1, IMAX ), 1 )
                  ROWMAX = max( ROWMAX, abs( A( JMAX, IMAX ) ) )
               end if
!  
               if ( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) then
!  
!                 no interchange, use 1-by-1 pivot block
!  
                  KP = K
               else if ( abs( A( IMAX, IMAX ) ) >= ALPHA*ROWMAX ) then
!  
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!  
                  KP = IMAX
               else
!  
!                 interchange rows and columns K-1 and IMAX, use 2-by-2
!                 pivot block
!  
                  KP = IMAX
                  KSTEP = 2
               end if
            end if
!  
            KK = K - KSTEP + 1
            if ( KP.NE.KK ) then
!  
!              Interchange rows and columns KK and KP in the leading
!              submatrix A(1:k,1:k)
!  
               call DSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 )
               call DSWAP( KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ), LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               if ( KSTEP == 2 ) then
                  T = A( K-1, K )
                  A( K-1, K ) = A( KP, K )
                  A( KP, K ) = T
               end if
            end if
!  
!           Update the leading submatrix
!  
            if ( KSTEP == 1 ) then
!  
!              1-by-1 pivot block D(k): column k now holds
!  
!              W(k) = U(k)*D(k)
!  
!              where U(k) is the k-th column of U
!  
!              Perform a rank-1 update of A(1:k-1,1:k-1) as
!  
!              A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'
!  
               R1 = ONE / A( K, K )
               call DSYR( UPLO, K-1, -R1, A( 1, K ), 1, A, LDA )
!  
!              Store U(k) in column k
!  
               call DSCAL( K-1, R1, A( 1, K ), 1 )
            else
!  
!              2-by-2 pivot block D(k): columns k and k-1 now hold
!  
!              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
!  
!              where U(k) and U(k-1) are the k-th and (k-1)-th columns
!              of U
!  
!              Perform a rank-2 update of A(1:k-2,1:k-2) as
!  
!              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
!                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
!  
               if ( K > 2 ) then
!  
                  D12 = A( K-1, K )
                  D22 = A( K-1, K-1 ) / D12
                  D11 = A( K, K ) / D12
                  T = ONE / ( D11*D22-ONE )
                  D12 = T / D12
!  
                  DO 30 J = K - 2, 1, -1
                     WKM1 = D12*( D11*A( J, K-1 )-A( J, K ) )
                     WK = D12*( D22*A( J, K )-A( J, K-1 ) )
                     DO 20 I = J, 1, -1
                        A(I,J) = A(I,J) - A(I,K)*WK - A(I,K-1)*WKM1
   20                CONTINUE
                     A( J, K ) = WK
                     A( J, K-1 ) = WKM1
   30             CONTINUE
!  
               end if
!  
            end if
         end if
!  
!        Store details of the interchanges in IPIV
!  
         if ( KSTEP == 1 ) then
            IPIV( K ) = KP
         else
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         end if
!  
!        Decrease K and return to the start of the main loop
!  
         K = K - KSTEP
         GO TO 10
!  
      else
!  
!        Factorize A as L*D*L' using the lower triangle of A
!  
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2
!  
         K = 1
   40    CONTINUE
!  
!        If K > N, exit from loop
!  
         if ( K > N ) GO TO 70
         KSTEP = 1
!  
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!  
         ABSAKK = abs( A( K, K ) )
!  
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!  
         if ( K < N ) then
            IMAX = K + IDAMAX( N-K, A( K+1, K ), 1 )
            COLMAX = abs( A( IMAX, K ) )
         else
            COLMAX = ZERO
         end if
!  
        call DISNANs(ABSAKK, nuestra)
         if ( (max( ABSAKK, COLMAX ) == ZERO) .or.  nuestra) then
!  
!           Column K is zero or contains a NaN: set INFO and continue
!  
            if ( INFO == 0 ) INFO = K
            KP = K
         else
            if ( ABSAKK >= ALPHA*COLMAX ) then
!  
!              no interchange, use 1-by-1 pivot block
!  
               KP = K
            else
!  
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!  
               JMAX = K - 1 + IDAMAX( IMAX-K, A( IMAX, K ), LDA )
               ROWMAX = abs( A( IMAX, JMAX ) )
               if ( IMAX < N ) then
                  JMAX = IMAX + IDAMAX( N-IMAX, A( IMAX+1, IMAX ), 1 )
                  ROWMAX = max( ROWMAX, abs( A( JMAX, IMAX ) ) )
               end if
!  
               if ( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) then
!  
!                 no interchange, use 1-by-1 pivot block
!  
                  KP = K
               else if ( abs( A( IMAX, IMAX ) ) >= ALPHA*ROWMAX ) then
!  
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!  
                  KP = IMAX
               else
!  
!                 interchange rows and columns K+1 and IMAX, use 2-by-2
!                 pivot block
!  
                  KP = IMAX
                  KSTEP = 2
               end if
            end if
!  
            KK = K + KSTEP - 1
            if ( KP.NE.KK ) then
!  
!              Interchange rows and columns KK and KP in the trailing
!              submatrix A(k:n,k:n)
!  
               if ( KP < N )  &
                  call DSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )
               call DSWAP( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ), LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               if ( KSTEP == 2 ) then
                  T = A( K+1, K )
                  A( K+1, K ) = A( KP, K )
                  A( KP, K ) = T
               end if
            end if
!  
!           Update the trailing submatrix
!  
            if ( KSTEP == 1 ) then
!  
!              1-by-1 pivot block D(k): column k now holds
!  
!              W(k) = L(k)*D(k)
!  
!              where L(k) is the k-th column of L
!  
               if ( K < N ) then
!  
!                 Perform a rank-1 update of A(k+1:n,k+1:n) as
!  
!                 A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/D(k))*W(k)'
!  
                  D11 = ONE / A( K, K )
                  call DSYR( UPLO, N-K, -D11, A( K+1, K ), 1,   &
                             A( K+1, K+1 ), LDA )
!  
!                 Store L(k) in column K
!  
                  call DSCAL( N-K, D11, A( K+1, K ), 1 )
               end if
            else
!  
!              2-by-2 pivot block D(k)
!  
               if ( K < N-1 ) then
!  
!                 Perform a rank-2 update of A(k+2:n,k+2:n) as
!  
!                 A := A - ( (A(k) A(k+1))*D(k)**(-1) ) * (A(k) A(k+1))'
!  
!                 where L(k) and L(k+1) are the k-th and (k+1)-th
!                 columns of L
!  
                  D21 = A( K+1, K )
                  D11 = A( K+1, K+1 ) / D21
                  D22 = A( K, K ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
!  
                  DO 60 J = K + 2, N
!  
                     WK = D21*( D11*A( J, K )-A( J, K+1 ) )
                     WKP1 = D21*( D22*A( J, K+1 )-A( J, K ) )
!  
                     DO 50 I = J, N
                        A(I,J) = A(I,J) - A(I,K)*WK - A(I,K+1)*WKP1
   50                CONTINUE
!  
                     A( J, K ) = WK
                     A( J, K+1 ) = WKP1
!  
   60             CONTINUE
               end if
            end if
         end if
!  
!        Store details of the interchanges in IPIV
!  
         if ( KSTEP == 1 ) then
            IPIV( K ) = KP
         else
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         end if
!  
!        Increase K and return to the start of the main loop
!  
         K = K + KSTEP
         GO TO 40
!  
      end if
!  
   70 CONTINUE
!  
      RETURN
!  
!     End of DSYTF2
!  
      end subroutine DSYTF2
      
      
      subroutine DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
!  
!  -- LAPACK routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!  
!     .. Scalar Arguments ..
      character          UPLO
      integer ::            INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      integer ::            IPIV( * )
      real(8) ::   A( LDA, * ), WORK( * )
!     ..
!  
!  Purpose
!  =======
!  
!  DSYTRF computes the factorization of a real symmetric matrix A using
!  the Bunch-Kaufman diagonal pivoting method.  The form of the
!  factorization is
!  
!     A = U*D*U**T  or  A = L*D*L**T
!  
!  where U (or L) is a product of permutation and unit upper (lower)
!  triangular matrices, and D is symmetric and block diagonal with
!  1-by-1 and 2-by-2 diagonal blocks.
!  
!  This is the blocked version of the algorithm, calling Level 3 BLAS.
!  
!  Arguments
!  =========
!  
!  UPLO    (input) character*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!  
!  N       (input) integer ::
!          The order of the matrix A.  N >= 0.
!  
!  A       (input/output) real(8) :: array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!  
!          On exit, the block diagonal matrix D and the multipliers used
!          to obtain the factor U or L (see below for further details).
!  
!  LDA     (input) integer ::
!          The leading dimension of the array A.  LDA >= max(1,N).
!  
!  IPIV    (output) integer :: array, dimension (N)
!          Details of the interchanges and the block structure of D.
!          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!          interchanged and D(k,k) is a 1-by-1 diagonal block.
!          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!  
!  WORK    (workspace/output) real(8) :: array, dimension (max(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!  
!  LWORK   (input) integer ::
!          The length of WORK.  LWORK >=1.  For best performance
!          LWORK >= N*NB, where NB is the block size returned by ILAENV.
!  
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!  
!  INFO    (output) integer ::
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
!                has been completed, but the block diagonal matrix D is
!                exactly singular, and division by zero will occur if it
!                is used to solve a system of equations.
!  
!  Further Details
!  ===============
!  
!  If UPLO = 'U', then A = U*D*U', where
!     U = P(n)*U(n)* ... *P(k)U(k)* ...,
!  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!  
!             (   I    v    0   )   k-s
!     U(k) =  (   0    I    0   )   s
!             (   0    0    I   )   n-k
!                k-s   s   n-k
!  
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!  and A(k,k), and v overwrites A(1:k-2,k-1:k).
!  
!  If UPLO = 'L', then A = L*D*L', where
!     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
!  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
!  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!  
!             (   I    0     0   )  k-1
!     L(k) =  (   0    I     0   )  s
!             (   0    v     I   )  n-k-s+1
!                k-1   s  n-k-s+1
!  
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!  
!  =====================================================================
!  
!     .. Local Scalars ..
      logical ::            LQUERY, UPPER
      integer ::            IINFO, IWS, J, K, KB, LDWORK, LWKOPT, NB, NBMIN
!     ..
!     .. External Functions ..
!m      logical ::            LSAME
!m      integer ::            ILAENV
!m      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
!m      EXTERNAL           DLASYF, DSYTF2, XERBLA
!     ..
!     .. Intrinsic Functions ..
!m      INTRINSIC          max
!     ..
!     .. Executable Statements ..
!  
!     Test the input parameters.
!  
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK == -1 )
      if ( .NOT.UPPER .and. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( LDA < max( 1, N ) ) then
         INFO = -4
      else if ( LWORK < 1 .and. .NOT.LQUERY ) then
         INFO = -7
      end if
!  
      if ( INFO == 0 ) then
!  
!        Determine the block size
!  
         NB = ILAENV( 1, 'DSYTRF', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      end if
!  
      if ( INFO.NE.0 ) then
         call XERBLA( 'DSYTRF', -INFO )
         RETURN
      else if ( LQUERY ) then
         RETURN
      end if
!  
      NBMIN = 2
      LDWORK = N
      if ( NB > 1 .and. NB < N ) then
         IWS = LDWORK*NB
         if ( LWORK < IWS ) then
            NB = max( LWORK / LDWORK, 1 )
            NBMIN = max( 2, ILAENV( 2, 'DSYTRF', UPLO, N, -1, -1, -1 ) )
         end if
      else
         IWS = 1
      end if
      if ( NB < NBMIN ) NB = N
!  
      if ( UPPER ) then
!  
!        Factorize A as U*D*U' using the upper triangle of A
!  
!        K is the main loop index, decreasing from N to 1 in steps of
!        KB, where KB is the number of columns factorized by DLASYF;
!        KB is either NB or NB-1, or K for the last block
!  
         K = N
   10    CONTINUE
!  
!        If K < 1, exit from loop
!  
         if ( K < 1 ) GO TO 40
!  
         if ( K > NB ) then
!  
!           Factorize columns k-kb+1:k of A and use blocked code to
!           update columns 1:k-kb
!  
            call DLASYF( UPLO, K, NB, KB, A, LDA, IPIV, WORK, LDWORK,  &
                         IINFO )
         else
!  
!           Use unblocked code to factorize columns 1:k of A
!  
            call DSYTF2( UPLO, K, A, LDA, IPIV, IINFO )
            KB = K
         end if
!  
!        Set INFO on the first occurrence of a zero pivot
!  
         if ( INFO == 0 .and. IINFO > 0 ) INFO = IINFO
!  
!        Decrease K and return to the start of the main loop
!  
         K = K - KB
         GO TO 10
!  
      else
!  
!        Factorize A as L*D*L' using the lower triangle of A
!  
!        K is the main loop index, increasing from 1 to N in steps of
!        KB, where KB is the number of columns factorized by DLASYF;
!        KB is either NB or NB-1, or N-K+1 for the last block
!  
         K = 1
   20    CONTINUE
!  
!        If K > N, exit from loop
!  
         if ( K > N ) GO TO 40
!  
         if ( K <= N-NB ) then
!  
!           Factorize columns k:k+kb-1 of A and use blocked code to
!           update columns k+kb:n
!  
            call DLASYF( UPLO, N-K+1, NB, KB, A( K, K ), LDA, IPIV( K ),  &
                         WORK, LDWORK, IINFO )
         else
!  
!           Use unblocked code to factorize columns k:n of A
!  
            call DSYTF2( UPLO, N-K+1, A( K, K ), LDA, IPIV( K ), IINFO )
            KB = N - K + 1
         end if
!  
!        Set INFO on the first occurrence of a zero pivot
!  
         if ( INFO == 0 .and. IINFO > 0 )  INFO = IINFO + K - 1
!  
!        Adjust IPIV
!  
         DO 30 J = K, K + KB - 1
            if ( IPIV( J ) > 0 ) then
               IPIV( J ) = IPIV( J ) + K - 1
            else
               IPIV( J ) = IPIV( J ) - K + 1
            end if
   30    CONTINUE
!  
!        Increase K and return to the start of the main loop
!  
         K = K + KB
         GO TO 20
!  
      end if
!  
   40 CONTINUE
      WORK( 1 ) = LWKOPT
      RETURN
!  
!     End of DSYTRF
!  
      end subroutine DSYTRF
      subroutine DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!  
!  -- LAPACK routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!  
!     .. Scalar Arguments ..
      character          UPLO
      integer ::            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      integer ::            IPIV( * )
      real(8) ::   A( LDA, * ), B( LDB, * )
!     ..
!  
!  Purpose
!  =======
!  
!  DSYTRS solves a system of linear equations A*X = B with a real
!  symmetric matrix A using the factorization A = U*D*U**T or
!  A = L*D*L**T computed by DSYTRF.
!  
!  Arguments
!  =========
!  
!  UPLO    (input) character*1
!          Specifies whether the details of the factorization are stored
!          as an upper or lower triangular matrix.
!          = 'U':  Upper triangular, form is A = U*D*U**T;
!          = 'L':  Lower triangular, form is A = L*D*L**T.
!  
!  N       (input) integer ::
!          The order of the matrix A.  N >= 0.
!  
!  NRHS    (input) integer ::
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!  
!  A       (input) real(8) :: array, dimension (LDA,N)
!          The block diagonal matrix D and the multipliers used to
!          obtain the factor U or L as computed by DSYTRF.
!  
!  LDA     (input) integer ::
!          The leading dimension of the array A.  LDA >= max(1,N).
!  
!  IPIV    (input) integer :: array, dimension (N)
!          Details of the interchanges and the block structure of D
!          as determined by DSYTRF.
!  
!  B       (input/output) real(8) :: array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!  
!  LDB     (input) integer ::
!          The leading dimension of the array B.  LDB >= max(1,N).
!  
!  INFO    (output) integer ::
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!  
!  =====================================================================
!  
!     .. Parameters ..
      real(8) ::   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      logical ::            UPPER
      integer ::            J, K, KP
      real(8) ::   AK, AKM1, AKM1K, BK, BKM1, DENOM
!     ..
!     .. External Functions ..
!m      logical ::            LSAME
!m      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!m      EXTERNAL           DGEMV, DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
!m      INTRINSIC          max
!     ..
!     .. Executable Statements ..
!  
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .and. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( NRHS < 0 ) then
         INFO = -3
      else if ( LDA < max( 1, N ) ) then
         INFO = -5
      else if ( LDB < max( 1, N ) ) then
         INFO = -8
      end if
      if ( INFO.NE.0 ) then
         call XERBLA( 'DSYTRS', -INFO )
         RETURN
      end if
!  
!     Quick return if possible
!  
      if ( N == 0 .or. NRHS == 0 ) RETURN
!  
      if ( UPPER ) then
!  
!        Solve A*X = B, where A = U*D*U'.
!
!        First solve U*D*X = B, overwriting B with X.
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = N
   10    CONTINUE
!
!        If K < 1, exit from loop.
!
         if ( K < 1 ) GO TO 30
!
         if ( IPIV( K ) > 0 ) then
!
!           1 x 1 diagonal block
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            if ( KP.NE.K ) &
               call DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in column K of A.
!
            call DGER( K-1, NRHS, -ONE, A( 1, K ), 1, B( K, 1 ), LDB,  &
                       B( 1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
            call DSCAL( NRHS, ONE / A( K, K ), B( K, 1 ), LDB )
            K = K - 1
         else
!
!           2 x 2 diagonal block
!
!           Interchange rows K-1 and -IPIV(K).
!
            KP = -IPIV( K )
            if ( KP.NE.K-1 ) &
               call DSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in columns K-1 and K of A.
!
            call DGER( K-2, NRHS, -ONE, A( 1, K ), 1, B( K, 1 ), LDB, &
                       B( 1, 1 ), LDB )
            call DGER( K-2, NRHS, -ONE, A( 1, K-1 ), 1, B( K-1, 1 ),  &
                       LDB, B( 1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
            AKM1K = A( K-1, K )
            AKM1 = A( K-1, K-1 ) / AKM1K
            AK = A( K, K ) / AKM1K
            DENOM = AKM1*AK - ONE
            DO 20 J = 1, NRHS
               BKM1 = B( K-1, J ) / AKM1K
               BK = B( K, J ) / AKM1K
               B( K-1, J ) = ( AK*BKM1-BK ) / DENOM
               B( K, J ) = ( AKM1*BK-BKM1 ) / DENOM
   20       CONTINUE
            K = K - 2
         end if
!
         GO TO 10
   30    CONTINUE
!
!        Next solve U'*X = B, overwriting B with X.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = 1
   40    CONTINUE
!
!        If K > N, exit from loop.
!
         if ( K > N ) GO TO 50
!
         if ( IPIV( K ) > 0 ) then
!
!           1 x 1 diagonal block
!
!           Multiply by inv(U'(K)), where U(K) is the transformation
!           stored in column K of A.
!
            call DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB, A( 1, K ), &
                        1, ONE, B( K, 1 ), LDB )
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            if ( KP.NE.K ) &
               call DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K + 1
         else
!
!           2 x 2 diagonal block
!
!           Multiply by inv(U'(K+1)), where U(K+1) is the transformation
!           stored in columns K and K+1 of A.
!
            call DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB, A( 1, K ),  &
                       1, ONE, B( K, 1 ), LDB )
            call DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB,             &
                        A( 1, K+1 ), 1, ONE, B( K+1, 1 ), LDB )
!
!           Interchange rows K and -IPIV(K).
!
            KP = -IPIV( K )
            if ( KP.NE.K ) &
               call DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K + 2
         end if
!
         GO TO 40
   50    CONTINUE
!
      else
!
!        Solve A*X = B, where A = L*D*L'.
!
!        First solve L*D*X = B, overwriting B with X.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = 1
   60    CONTINUE
!
!        If K > N, exit from loop.
!
         if ( K > N ) GO TO 80
!
         if ( IPIV( K ) > 0 ) then
!
!           1 x 1 diagonal block
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            if ( KP.NE.K )  &
               call DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in column K of A.
!
            if ( K < N )  &
               call DGER( N-K, NRHS, -ONE, A( K+1, K ), 1, B( K, 1 ),  &
                          LDB, B( K+1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
            call DSCAL( NRHS, ONE / A( K, K ), B( K, 1 ), LDB )
            K = K + 1
         else
!
!           2 x 2 diagonal block
!
!           Interchange rows K+1 and -IPIV(K).
!
            KP = -IPIV( K )
            if ( KP.NE.K+1 ) &
               call DSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in columns K and K+1 of A.
!
            if ( K < N-1 ) then
               call DGER( N-K-1, NRHS, -ONE, A( K+2, K ), 1, B( K, 1 ),  &
                          LDB, B( K+2, 1 ), LDB )
               call DGER( N-K-1, NRHS, -ONE, A( K+2, K+1 ), 1,   &
                          B( K+1, 1 ), LDB, B( K+2, 1 ), LDB )
            end if
!
!           Multiply by the inverse of the diagonal block.
!
            AKM1K = A( K+1, K )
            AKM1 = A( K, K ) / AKM1K
            AK = A( K+1, K+1 ) / AKM1K
            DENOM = AKM1*AK - ONE
            DO 70 J = 1, NRHS
               BKM1 = B( K, J ) / AKM1K
               BK = B( K+1, J ) / AKM1K
               B( K, J ) = ( AK*BKM1-BK ) / DENOM
               B( K+1, J ) = ( AKM1*BK-BKM1 ) / DENOM
   70       CONTINUE
            K = K + 2
         end if
!
         GO TO 60
   80    CONTINUE
!
!        Next solve L'*X = B, overwriting B with X.
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = N
   90    CONTINUE
!
!        If K < 1, exit from loop.
!
         if ( K < 1 )  GO TO 100
!
         if ( IPIV( K ) > 0 ) then
!
!           1 x 1 diagonal block
!
!           Multiply by inv(L'(K)), where L(K) is the transformation
!           stored in column K of A.
!
            if ( K < N ) &
               call DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ), &
                           LDB, A( K+1, K ), 1, ONE, B( K, 1 ), LDB )
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            if ( KP.NE.K ) &
               call DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K - 1
         else
!
!           2 x 2 diagonal block
!
!           Multiply by inv(L'(K-1)), where L(K-1) is the transformation
!           stored in columns K-1 and K of A.
!
            if ( K < N ) then
               call DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ),   &
                           LDB, A( K+1, K ), 1, ONE, B( K, 1 ), LDB )
               call DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ),   &
                           LDB, A( K+1, K-1 ), 1, ONE, B( K-1, 1 ),     &
                           LDB )
            end if
!
!           Interchange rows K and -IPIV(K).
!
            KP = -IPIV( K )
            if ( KP.NE.K )   &
               call DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K - 2
         end if
!
         GO TO 90
  100    CONTINUE
      end if
!
      RETURN
!
!     End of DSYTRS
!
      end subroutine DSYTRS











!----------------------------------------------------------------------
!   lapack_sysv_full
!----------------------------------------------------------------------
 subroutine SISNANs( SIN, SISNAN )
!
!  -- LAPACK auxiliary routine (version 3.2.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2010
!
!     .. Scalar Arguments ..
      real, intent(in) :: SIN
      logical, intent(out) :: SISNAN
!     ..
!
!  Purpose
!  =======
!
!  SISNAN returns .TRUE. if its argument is NaN, and .FALSE.
!  otherwise.  To be replaced by the Fortran 2003 intrinsic in the
!  future.
!
!  Arguments
!  =========
!
!  SIN     (input) real ::
!          Input to test for NaN.
!
!  =====================================================================
!
!  .. External Functions ..
!m      logical :: SLAISNAN
!m      EXTERNAL SLAISNAN
!  ..
!  .. Executable Statements ..
      SISNAN = SLAISNAN(SIN,SIN)

end subroutine SISNANs

logical function SLAISNAN( SIN1, SIN2 )
!
!  -- LAPACK auxiliary routine (version 3.2.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2010
!
!     .. Scalar Arguments ..
      real ::               SIN1, SIN2
!     ..
!
!  Purpose
!  =======
!
!  This routine is not for general use.  It exists solely to avoid
!  over-optimization in SISNAN.
!
!  SLAISNAN checks for NaNs by comparing its two arguments for
!  inequality.  NaN is the only floating-point value where NaN != NaN
!  returns .TRUE.  To check for NaNs, pass the same variable as both
!  arguments.
!
!  A compiler must assume that the two arguments are
!  not the same variable, and the test will not be optimized away.
!  Interprocedural or whole-program optimization may delete this
!  test.  The ISNAN functions will be replaced by the correct
!  Fortran 03 intrinsic once the intrinsic is widely available.
!
!  Arguments
!  =========
!
!  SIN1     (input) real ::
!
!  SIN2     (input) real ::
!          Two numbers to compare for inequality.
!
!  =====================================================================
!
!  .. Executable Statements ..
      SLAISNAN = (SIN1.NE.SIN2)
      RETURN
      end function SLAISNAN
      subroutine SLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )
!
!  -- LAPACK routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!
!     .. Scalar Arguments ..
      character          UPLO
      integer ::            INFO, KB, LDA, LDW, N, NB
!     ..
!     .. Array Arguments ..
      integer ::            IPIV( * )
      real ::               A( LDA, * ), W( LDW, * )
!     ..
!
!  Purpose
!  =======
!
!  SLASYF computes a partial factorization of a real symmetric matrix A
!  using the Bunch-Kaufman diagonal pivoting method. The partial
!  factorization has the form:
!
!  A  =  ( I  U12 ) ( A11  0  ) (  I    0   )  if UPLO = 'U', or:
!        ( 0  U22 ) (  0   D  ) ( U12' U22' )
!
!  A  =  ( L11  0 ) (  D   0  ) ( L11' L21' )  if UPLO = 'L'
!        ( L21  I ) (  0  A22 ) (  0    I   )
!
!  where the order of D is at most NB. The actual order is returned in
!  the argument KB, and is either NB or NB-1, or N if N <= NB.
!
!  SLASYF is an auxiliary routine called by SSYTRF. It uses blocked code
!  (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
!  A22 (if UPLO = 'L').
!
!  Arguments
!  =========
!
!  UPLO    (input) character*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) integer ::
!          The order of the matrix A.  N >= 0.
!
!  NB      (input) integer ::
!          The maximum number of columns of the matrix A that should be
!          factored.  NB should be at least 2 to allow for 2-by-2 pivot
!          blocks.
!
!  KB      (output) integer ::
!          The number of columns of A that were actually factored.
!          KB is either NB-1 or NB, or N if N <= NB.
!
!  A       (input/output) real :: array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit, A contains details of the partial factorization.
!
!  LDA     (input) integer ::
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (output) integer :: array, dimension (N)
!          Details of the interchanges and the block structure of D.
!          If UPLO = 'U', only the last KB elements of IPIV are set;
!          if UPLO = 'L', only the first KB elements are set.
!
!          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!          interchanged and D(k,k) is a 1-by-1 diagonal block.
!          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!
!  W       (workspace) real :: array, dimension (LDW,NB)
!
!  LDW     (input) integer ::
!          The leading dimension of the array W.  LDW >= max(1,N).
!
!  INFO    (output) integer ::
!          = 0: successful exit
!          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
!               has been completed, but the block diagonal matrix D is
!               exactly singular.
!
!  =====================================================================
!
!     .. Parameters ..
      real ::               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      real ::               EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0E+0, SEVTEN = 17.0E+0 )
!     ..
!     .. Local Scalars ..
      integer ::            IMAX, J, JB, JJ, JMAX, JP, K, KK, KKW, KP
      integer ::            KSTEP, KW
      real ::               ABSAKK, ALPHA, COLMAX, D11, D21, D22, R1
      real ::               ROWMAX, T
!     ..
!     .. External Functions ..
!m      logical ::            LSAME
!m      integer ::            ISAMAX
!m      EXTERNAL           LSAME, ISAMAX
!     ..
!     .. External Subroutines ..
!m      EXTERNAL           SCOPY, SGEMM, SGEMV, SSCAL, SSWAP
!     ..
!     .. Intrinsic Functions ..
!m      INTRINSIC          abs, max, min, sqrt
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
!     Initialize ALPHA for use in choosing pivot block size.
!
      ALPHA = ( ONE+sqrt( SEVTEN ) ) / EIGHT
!
      if ( LSAME( UPLO, 'U' ) ) then
!
!        Factorize the trailing columns of A using the upper triangle
!        of A and working backwards, and compute the matrix W = U12*D
!        for use in updating A11
!
!        K is the main loop index, decreasing from N in steps of 1 or 2
!
!        KW is the column of W which corresponds to column K of A
!
         K = N
   10    CONTINUE
         KW = NB + K - N
!
!        Exit from loop
!
         if ( ( K <= N-NB+1 .and. NB < N ) .or. K < 1 ) GO TO 30
!
!        Copy column K of A to column KW of W and update it
!
         call SCOPY( K, A( 1, K ), 1, W( 1, KW ), 1 )
         if ( K < N ) &
            call SGEMV( 'No transpose', K, N-K, -ONE, A( 1, K+1 ), LDA, &
                        W( K, KW+1 ), LDW, ONE, W( 1, KW ), 1 )
!
         KSTEP = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
         ABSAKK = abs( W( K, KW ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
         if ( K > 1 ) then
            IMAX = ISAMAX( K-1, W( 1, KW ), 1 )
            COLMAX = abs( W( IMAX, KW ) )
         else
            COLMAX = ZERO
         end if
!
         if ( max( ABSAKK, COLMAX ) == ZERO ) then
!
!           Column K is zero: set INFO and continue
!
            if ( INFO == 0 ) INFO = K
            KP = K
         else
            if ( ABSAKK >= ALPHA*COLMAX ) then
!
!              no interchange, use 1-by-1 pivot block
!
               KP = K
            else
!
!              Copy column IMAX to column KW-1 of W and update it
!
               call SCOPY( IMAX, A( 1, IMAX ), 1, W( 1, KW-1 ), 1 )
               call SCOPY( K-IMAX, A( IMAX, IMAX+1 ), LDA,                &
                           W( IMAX+1, KW-1 ), 1 )
               if ( K < N )                                               &
                  call SGEMV( 'No transpose', K, N-K, -ONE, A( 1, K+1 ),  &
                              LDA, W( IMAX, KW+1 ), LDW, ONE,             &
                              W( 1, KW-1 ), 1 )
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
               JMAX = IMAX + ISAMAX( K-IMAX, W( IMAX+1, KW-1 ), 1 )
               ROWMAX = abs( W( JMAX, KW-1 ) )
               if ( IMAX > 1 ) then
                  JMAX = ISAMAX( IMAX-1, W( 1, KW-1 ), 1 )
                  ROWMAX = max( ROWMAX, abs( W( JMAX, KW-1 ) ) )
               end if
!
               if ( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) then
!
!                 no interchange, use 1-by-1 pivot block
!
                  KP = K
               else if ( abs( W( IMAX, KW-1 ) ) >= ALPHA*ROWMAX ) then
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                  KP = IMAX
!
!                 copy column KW-1 of W to column KW
!
                  call SCOPY( K, W( 1, KW-1 ), 1, W( 1, KW ), 1 )
               else
!
!                 interchange rows and columns K-1 and IMAX, use 2-by-2
!                 pivot block
!
                  KP = IMAX
                  KSTEP = 2
               end if
            end if
!
            KK = K - KSTEP + 1
            KKW = NB + KK - N
!
!           Updated column KP is already stored in column KKW of W
!
            if ( KP.NE.KK ) then
!
!              Copy non-updated column KK to column KP
!
               A( KP, K ) = A( KK, K )
               call SCOPY( K-1-KP, A( KP+1, KK ), 1, A( KP, KP+1 ), LDA )
               call SCOPY( KP, A( 1, KK ), 1, A( 1, KP ), 1 )
!
!              Interchange rows KK and KP in last KK columns of A and W
!
               call SSWAP( N-KK+1, A( KK, KK ), LDA, A( KP, KK ), LDA )
               call SSWAP( N-KK+1, W( KK, KKW ), LDW, W( KP, KKW ), LDW )
            end if
!
            if ( KSTEP == 1 ) then
!
!              1-by-1 pivot block D(k): column KW of W now holds
!
!              W(k) = U(k)*D(k)
!
!              where U(k) is the k-th column of U
!
!              Store U(k) in column k of A
!
               call SCOPY( K, W( 1, KW ), 1, A( 1, K ), 1 )
               R1 = ONE / A( K, K )
               call SSCAL( K-1, R1, A( 1, K ), 1 )
            else
!
!              2-by-2 pivot block D(k): columns KW and KW-1 of W now
!              hold
!
!              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
!
!              where U(k) and U(k-1) are the k-th and (k-1)-th columns
!              of U
!
               if ( K > 2 ) then
!
!                 Store U(k) and U(k-1) in columns k and k-1 of A
!
                  D21 = W( K-1, KW )
                  D11 = W( K, KW ) / D21
                  D22 = W( K-1, KW-1 ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
                  DO 20 J = 1, K - 2
                     A( J, K-1 ) = D21*( D11*W( J, KW-1 )-W( J, KW ) )
                     A( J, K ) = D21*( D22*W( J, KW )-W( J, KW-1 ) )
   20             CONTINUE
               end if
!
!              Copy D(k) to A
!
               A( K-1, K-1 ) = W( K-1, KW-1 )
               A( K-1, K ) = W( K-1, KW )
               A( K, K ) = W( K, KW )
            end if
         end if
!
!        Store details of the interchanges in IPIV
!
         if ( KSTEP == 1 ) then
            IPIV( K ) = KP
         else
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         end if
!
!        Decrease K and return to the start of the main loop
!
         K = K - KSTEP
         GO TO 10
!
   30    CONTINUE
!
!        Update the upper triangle of A11 (= A(1:k,1:k)) as
!
!        A11 := A11 - U12*D*U12' = A11 - U12*W'
!
!        computing blocks of NB columns at a time
!
         DO 50 J = ( ( K-1 ) / NB )*NB + 1, 1, -NB
            JB = min( NB, K-J+1 )
!
!           Update the upper triangle of the diagonal block
!
            DO 40 JJ = J, J + JB - 1
               call SGEMV( 'No transpose', JJ-J+1, N-K, -ONE,           &
                           A( J, K+1 ), LDA, W( JJ, KW+1 ), LDW, ONE,   &
                           A( J, JJ ), 1 )
   40       CONTINUE
!
!           Update the rectangular superdiagonal block
!
            call SGEMM( 'No transpose', 'Transpose', J-1, JB, N-K, -ONE,  &
                        A( 1, K+1 ), LDA, W( J, KW+1 ), LDW, ONE,         &
                        A( 1, J ), LDA )
   50    CONTINUE
!
!        Put U12 in standard form by partially undoing the interchanges
!        in columns k+1:n
!
         J = K + 1
   60    CONTINUE
         JJ = J
         JP = IPIV( J )
         if ( JP < 0 ) then
            JP = -JP
            J = J + 1
         end if
         J = J + 1
         if ( JP.NE.JJ .and. J <= N )  &
            call SSWAP( N-J+1, A( JP, J ), LDA, A( JJ, J ), LDA )
         if ( J <= N ) GO TO 60
!
!        Set KB to the number of columns factorized
!
         KB = N - K
!
      else
!
!        Factorize the leading columns of A using the lower triangle
!        of A and working forwards, and compute the matrix W = L21*D
!        for use in updating A22
!
!        K is the main loop index, increasing from 1 in steps of 1 or 2
!
         K = 1
   70    CONTINUE
!
!        Exit from loop
!
         if ( ( K >= NB .and. NB < N ) .or. K > N ) GO TO 90
!
!        Copy column K of A to column K of W and update it
!
         call SCOPY( N-K+1, A( K, K ), 1, W( K, K ), 1 )
         call SGEMV( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ), LDA,  &
                     W( K, 1 ), LDW, ONE, W( K, K ), 1 )
!
         KSTEP = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
         ABSAKK = abs( W( K, K ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
         if ( K < N ) then
            IMAX = K + ISAMAX( N-K, W( K+1, K ), 1 )
            COLMAX = abs( W( IMAX, K ) )
         else
            COLMAX = ZERO
         end if
!
         if ( max( ABSAKK, COLMAX ) == ZERO ) then
!
!           Column K is zero: set INFO and continue
!
            if ( INFO == 0 ) INFO = K
            KP = K
         else
            if ( ABSAKK >= ALPHA*COLMAX ) then
!
!              no interchange, use 1-by-1 pivot block
!
               KP = K
            else
!
!              Copy column IMAX to column K+1 of W and update it
!
               call SCOPY( IMAX-K, A( IMAX, K ), LDA, W( K, K+1 ), 1 )
               call SCOPY( N-IMAX+1, A( IMAX, IMAX ), 1, W( IMAX, K+1 ),  &
                           1 )
               call SGEMV( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ),   &
                           LDA, W( IMAX, 1 ), LDW, ONE, W( K, K+1 ), 1 )
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
               JMAX = K - 1 + ISAMAX( IMAX-K, W( K, K+1 ), 1 )
               ROWMAX = abs( W( JMAX, K+1 ) )
               if ( IMAX < N ) then
                  JMAX = IMAX + ISAMAX( N-IMAX, W( IMAX+1, K+1 ), 1 )
                  ROWMAX = max( ROWMAX, abs( W( JMAX, K+1 ) ) )
               end if
!
               if ( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) then
!
!                 no interchange, use 1-by-1 pivot block
!
                  KP = K
               else if ( abs( W( IMAX, K+1 ) ) >= ALPHA*ROWMAX ) then
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                  KP = IMAX
!
!                 copy column K+1 of W to column K
!
                  call SCOPY( N-K+1, W( K, K+1 ), 1, W( K, K ), 1 )
               else
!
!                 interchange rows and columns K+1 and IMAX, use 2-by-2
!                 pivot block
!
                  KP = IMAX
                  KSTEP = 2
               end if
            end if
!
            KK = K + KSTEP - 1
!
!           Updated column KP is already stored in column KK of W
!
            if ( KP.NE.KK ) then
!
!              Copy non-updated column KK to column KP
!
               A( KP, K ) = A( KK, K )
               call SCOPY( KP-K-1, A( K+1, KK ), 1, A( KP, K+1 ), LDA )
               call SCOPY( N-KP+1, A( KP, KK ), 1, A( KP, KP ), 1 )
!
!              Interchange rows KK and KP in first KK columns of A and W
!
               call SSWAP( KK, A( KK, 1 ), LDA, A( KP, 1 ), LDA )
               call SSWAP( KK, W( KK, 1 ), LDW, W( KP, 1 ), LDW )
            end if
!
            if ( KSTEP == 1 ) then
!
!              1-by-1 pivot block D(k): column k of W now holds
!
!              W(k) = L(k)*D(k)
!
!              where L(k) is the k-th column of L
!
!              Store L(k) in column k of A
!
               call SCOPY( N-K+1, W( K, K ), 1, A( K, K ), 1 )
               if ( K < N ) then
                  R1 = ONE / A( K, K )
                  call SSCAL( N-K, R1, A( K+1, K ), 1 )
               end if
            else
!
!              2-by-2 pivot block D(k): columns k and k+1 of W now hold
!
!              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
!
!              where L(k) and L(k+1) are the k-th and (k+1)-th columns
!              of L
!
               if ( K < N-1 ) then
!
!                 Store L(k) and L(k+1) in columns k and k+1 of A
!
                  D21 = W( K+1, K )
                  D11 = W( K+1, K+1 ) / D21
                  D22 = W( K, K ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
                  DO 80 J = K + 2, N
                     A( J, K ) = D21*( D11*W( J, K )-W( J, K+1 ) )
                     A( J, K+1 ) = D21*( D22*W( J, K+1 )-W( J, K ) )
   80             CONTINUE
               end if
!
!              Copy D(k) to A
!
               A( K, K ) = W( K, K )
               A( K+1, K ) = W( K+1, K )
               A( K+1, K+1 ) = W( K+1, K+1 )
            end if
         end if
!
!        Store details of the interchanges in IPIV
!
         if ( KSTEP == 1 ) then
            IPIV( K ) = KP
         else
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         end if
!
!        Increase K and return to the start of the main loop
!
         K = K + KSTEP
         GO TO 70
!
   90    CONTINUE
!
!        Update the lower triangle of A22 (= A(k:n,k:n)) as
!
!        A22 := A22 - L21*D*L21' = A22 - L21*W'
!
!        computing blocks of NB columns at a time
!
         DO 110 J = K, N, NB
            JB = min( NB, N-J+1 )
!
!           Update the lower triangle of the diagonal block
!
            DO 100 JJ = J, J + JB - 1
               call SGEMV( 'No transpose', J+JB-JJ, K-1, -ONE,       &
                           A( JJ, 1 ), LDA, W( JJ, 1 ), LDW, ONE,    &
                           A( JJ, JJ ), 1 )
  100       CONTINUE
!
!           Update the rectangular subdiagonal block
!
            if ( J+JB <= N )  &
               call SGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB,     &
                           K-1, -ONE, A( J+JB, 1 ), LDA, W( J, 1 ), LDW,  &
                           ONE, A( J+JB, J ), LDA )
  110    CONTINUE
!
!        Put L21 in standard form by partially undoing the interchanges
!        in columns 1:k-1
!
         J = K - 1
  120    CONTINUE
         JJ = J
         JP = IPIV( J )
         if ( JP < 0 ) then
            JP = -JP
            J = J - 1
         end if
         J = J - 1
         if ( JP.NE.JJ .and. J >= 1 ) &
            call SSWAP( J, A( JP, 1 ), LDA, A( JJ, 1 ), LDA )
         if ( J >= 1 ) GO TO 120
!
!        Set KB to the number of columns factorized
!
         KB = K - 1
!
      end if
      RETURN
!
!     End of SLASYF
!
      end subroutine SLASYF
      subroutine SSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,   &
                        LWORK, INFO )
!
!  -- LAPACK driver routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!
!     .. Scalar Arguments ..
      character          UPLO
      integer ::            INFO, LDA, LDB, LWORK, N, NRHS
!     ..
!     .. Array Arguments ..
      integer ::            IPIV( * )
      real ::               A( LDA, * ), B( LDB, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  SSYSV computes the solution to a real system of linear equations
!     A * X = B,
!  where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
!  matrices.
!
!  The diagonal pivoting method is used to factor A as
!     A = U * D * U**T,  if UPLO = 'U', or
!     A = L * D * L**T,  if UPLO = 'L',
!  where U (or L) is a product of permutation and unit upper (lower)
!  triangular matrices, and D is symmetric and block diagonal with 
!  1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
!  used to solve the system of equations A * X = B.
!
!  Arguments
!  =========
!
!  UPLO    (input) character*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) integer ::
!          The number of linear equations, i.e., the order of the
!          matrix A.  N >= 0.
!
!  NRHS    (input) integer ::
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input/output) real :: array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the block diagonal matrix D and the
!          multipliers used to obtain the factor U or L from the
!          factorization A = U*D*U**T or A = L*D*L**T as computed by
!          SSYTRF.
!
!  LDA     (input) integer ::
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (output) integer :: array, dimension (N)
!          Details of the interchanges and the block structure of D, as
!          determined by SSYTRF.  If IPIV(k) > 0, then rows and columns
!          k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1
!          diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0,
!          then rows and columns k-1 and -IPIV(k) were interchanged and
!          D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and
!          IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and
!          -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2
!          diagonal block.
!
!  B       (input/output) real :: array, dimension (LDB,NRHS)
!          On entry, the N-by-NRHS right hand side matrix B.
!          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!
!  LDB     (input) integer ::
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  WORK    (workspace/output) real :: array, dimension (max(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) integer ::
!          The length of WORK.  LWORK >= 1, and for best performance
!          LWORK >= max(1,N*NB), where NB is the optimal blocksize for
!          SSYTRF.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) integer ::
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
!               has been completed, but the block diagonal matrix D is
!               exactly singular, so the solution could not be computed.
!
!  =====================================================================
!
!     .. Local Scalars ..
      logical ::            LQUERY
      integer ::            LWKOPT, NB
!     ..
!     .. External Functions ..
!m      logical ::            LSAME
!m       integer ::            ILAENV
!m       EXTERNAL           ILAENV, LSAME
!     ..
!     .. External Subroutines ..
!m       EXTERNAL           SSYTRF, SSYTRS, XERBLA
!     ..
!     .. Intrinsic Functions ..
!m       INTRINSIC          max
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      LQUERY = ( LWORK == -1 )
      if ( .NOT.LSAME( UPLO, 'U' ) .and. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( NRHS < 0 ) then
         INFO = -3
      else if ( LDA < max( 1, N ) ) then
         INFO = -5
      else if ( LDB < max( 1, N ) ) then
         INFO = -8
      else if ( LWORK < 1 .and. .NOT.LQUERY ) then
         INFO = -10
      end if
!
      if ( INFO == 0 ) then
         if ( N == 0 ) then
            LWKOPT = 1
         else
            NB = ILAENV( 1, 'SSYTRF', UPLO, N, -1, -1, -1 )
            LWKOPT = N*NB
         end if
         WORK( 1 ) = LWKOPT
      end if
!
      if ( INFO.NE.0 ) then
         call XERBLA( 'SSYSV ', -INFO )
         RETURN
      else if ( LQUERY ) then
         RETURN
      end if
!
!     Compute the factorization A = U*D*U' or A = L*D*L'.
!
      call SSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
      if ( INFO == 0 ) then
!
!        Solve the system A*X = B, overwriting B with X.
!
         call SSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
      end if
!
      WORK( 1 ) = LWKOPT
!
      RETURN
!
!     End of SSYSV
!
      end subroutine SSYSV
      subroutine SSYTF2( UPLO, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!
!     .. Scalar Arguments ..
      character          UPLO
      integer ::            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      integer ::            IPIV( * )
      real ::               A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  SSYTF2 computes the factorization of a real symmetric matrix A using
!  the Bunch-Kaufman diagonal pivoting method:
!
!     A = U*D*U'  or  A = L*D*L'
!
!  where U (or L) is a product of permutation and unit upper (lower)
!  triangular matrices, U' is the transpose of U, and D is symmetric and
!  block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) character*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) integer ::
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) real :: array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, the block diagonal matrix D and the multipliers used
!          to obtain the factor U or L (see below for further details).
!
!  LDA     (input) integer ::
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (output) integer :: array, dimension (N)
!          Details of the interchanges and the block structure of D.
!          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!          interchanged and D(k,k) is a 1-by-1 diagonal block.
!          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!
!  INFO    (output) integer ::
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
!               has been completed, but the block diagonal matrix D is
!               exactly singular, and division by zero will occur if it
!               is used to solve a system of equations.
!
!  Further Details
!  ===============
!
!  09-29-06 - patch from
!    Bobby Cheng, MathWorks
!
!    Replace l.204 and l.372
!         if ( max( ABSAKK, COLMAX ) == ZERO ) then
!    by
!         if ( (max( ABSAKK, COLMAX ) == ZERO) .or. SISNAN(ABSAKK) ) then
!
!  01-01-96 - Based on modifications by
!    J. Lewis, Boeing Computer Services Company
!    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!  1-96 - Based on modifications by J. Lewis, Boeing Computer Services
!         Company
!
!  If UPLO = 'U', then A = U*D*U', where
!     U = P(n)*U(n)* ... *P(k)U(k)* ...,
!  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!
!             (   I    v    0   )   k-s
!     U(k) =  (   0    I    0   )   s
!             (   0    0    I   )   n-k
!                k-s   s   n-k
!
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!  and A(k,k), and v overwrites A(1:k-2,k-1:k).
!
!  If UPLO = 'L', then A = L*D*L', where
!     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
!  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
!  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!
!             (   I    0     0   )  k-1
!     L(k) =  (   0    I     0   )  s
!             (   0    v     I   )  n-k-s+1
!                k-1   s  n-k-s+1
!
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!
!  =====================================================================
!
!     .. Parameters ..
      real ::               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      real ::               EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0E+0, SEVTEN = 17.0E+0 )
!     ..
!     .. Local Scalars ..
      logical ::            UPPER
      integer ::            I, IMAX, J, JMAX, K, KK, KP, KSTEP
      real ::               ABSAKK, ALPHA, COLMAX, D11, D12, D21, D22, R1
      real ::               ROWMAX, T, WK, WKM1, WKP1
      logical :: nuestra
!     ..
!     .. External Functions ..
!m       logical ::            LSAME, SISNAN
!m       integer ::            ISAMAX
!m       EXTERNAL           LSAME, ISAMAX, SISNAN
!     ..
!     .. External Subroutines ..
!m       EXTERNAL           SSCAL, SSWAP, SSYR, XERBLA
!     ..
!     .. Intrinsic Functions ..
!m       INTRINSIC          abs, max, sqrt
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .and. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( LDA < max( 1, N ) ) then
         INFO = -4
      end if
      if ( INFO.NE.0 ) then
         call XERBLA( 'SSYTF2', -INFO )
         RETURN
      end if
!
!     Initialize ALPHA for use in choosing pivot block size.
!
      ALPHA = ( ONE+sqrt( SEVTEN ) ) / EIGHT
!
      if ( UPPER ) then
!
!        Factorize A as U*D*U' using the upper triangle of A
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2
!
         K = N
   10    CONTINUE
!
!        If K < 1, exit from loop
!
         if ( K < 1 ) GO TO 70
         KSTEP = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
         ABSAKK = abs( A( K, K ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
         if ( K > 1 ) then
            IMAX = ISAMAX( K-1, A( 1, K ), 1 )
            COLMAX = abs( A( IMAX, K ) )
         else
            COLMAX = ZERO
         end if
!
         call SISNANs(ABSAKK, nuestra)
         if ( (max( ABSAKK, COLMAX ) == ZERO) .or. nuestra ) then
!
!           Column K is zero or contains a NaN: set INFO and continue
!
            if ( INFO == 0 ) INFO = K
            KP = K
         else
            if ( ABSAKK >= ALPHA*COLMAX ) then
!
!              no interchange, use 1-by-1 pivot block
!
               KP = K
            else
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
               JMAX = IMAX + ISAMAX( K-IMAX, A( IMAX, IMAX+1 ), LDA )
               ROWMAX = abs( A( IMAX, JMAX ) )
               if ( IMAX > 1 ) then
                  JMAX = ISAMAX( IMAX-1, A( 1, IMAX ), 1 )
                  ROWMAX = max( ROWMAX, abs( A( JMAX, IMAX ) ) )
               end if
!
               if ( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) then
!
!                 no interchange, use 1-by-1 pivot block
!
                  KP = K
               else if ( abs( A( IMAX, IMAX ) ) >= ALPHA*ROWMAX ) then
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                  KP = IMAX
               else
!
!                 interchange rows and columns K-1 and IMAX, use 2-by-2
!                 pivot block
!
                  KP = IMAX
                  KSTEP = 2
               end if
            end if
!
            KK = K - KSTEP + 1
            if ( KP.NE.KK ) then
!
!              Interchange rows and columns KK and KP in the leading
!              submatrix A(1:k,1:k)
!
               call SSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 )
               call SSWAP( KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ), LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               if ( KSTEP == 2 ) then
                  T = A( K-1, K )
                  A( K-1, K ) = A( KP, K )
                  A( KP, K ) = T
               end if
            end if
!
!           Update the leading submatrix
!
            if ( KSTEP == 1 ) then
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = U(k)*D(k)
!
!              where U(k) is the k-th column of U
!
!              Perform a rank-1 update of A(1:k-1,1:k-1) as
!
!              A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'
!
               R1 = ONE / A( K, K )
               call SSYR( UPLO, K-1, -R1, A( 1, K ), 1, A, LDA )
!
!              Store U(k) in column k
!
               call SSCAL( K-1, R1, A( 1, K ), 1 )
            else
!
!              2-by-2 pivot block D(k): columns k and k-1 now hold
!
!              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
!
!              where U(k) and U(k-1) are the k-th and (k-1)-th columns
!              of U
!
!              Perform a rank-2 update of A(1:k-2,1:k-2) as
!
!              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
!                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
!
               if ( K > 2 ) then
!
                  D12 = A( K-1, K )
                  D22 = A( K-1, K-1 ) / D12
                  D11 = A( K, K ) / D12
                  T = ONE / ( D11*D22-ONE )
                  D12 = T / D12
!
                  DO 30 J = K - 2, 1, -1
                     WKM1 = D12*( D11*A( J, K-1 )-A( J, K ) )
                     WK = D12*( D22*A( J, K )-A( J, K-1 ) )
                     DO 20 I = J, 1, -1
                        A(I,J) = A(I,J) - A(I,K)*WK - A(I,K-1)*WKM1
   20                CONTINUE
                     A( J, K ) = WK
                     A( J, K-1 ) = WKM1
   30             CONTINUE
!
               end if
!
            end if
         end if
!
!        Store details of the interchanges in IPIV
!
         if ( KSTEP == 1 ) then
            IPIV( K ) = KP
         else
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         end if
!
!        Decrease K and return to the start of the main loop
!
         K = K - KSTEP
         GO TO 10
!
      else
!
!        Factorize A as L*D*L' using the lower triangle of A
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2
!
         K = 1
   40    CONTINUE
!
!        If K > N, exit from loop
!
         if ( K > N ) GO TO 70
         KSTEP = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
         ABSAKK = abs( A( K, K ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
         if ( K < N ) then
            IMAX = K + ISAMAX( N-K, A( K+1, K ), 1 )
            COLMAX = abs( A( IMAX, K ) )
         else
            COLMAX = ZERO
         end if
!
         call SISNANs(ABSAKK, nuestra)
         if ( (max( ABSAKK, COLMAX ) == ZERO) .or. nuestra ) then
!
!           Column K is zero or contains a NaN: set INFO and continue
!
            if ( INFO == 0 )  INFO = K
            KP = K
         else
            if ( ABSAKK >= ALPHA*COLMAX ) then
!
!              no interchange, use 1-by-1 pivot block
!
               KP = K
            else
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
               JMAX = K - 1 + ISAMAX( IMAX-K, A( IMAX, K ), LDA )
               ROWMAX = abs( A( IMAX, JMAX ) )
               if ( IMAX < N ) then
                  JMAX = IMAX + ISAMAX( N-IMAX, A( IMAX+1, IMAX ), 1 )
                  ROWMAX = max( ROWMAX, abs( A( JMAX, IMAX ) ) )
               end if
!
               if ( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) then
!
!                 no interchange, use 1-by-1 pivot block
!
                  KP = K
               else if ( abs( A( IMAX, IMAX ) ) >= ALPHA*ROWMAX ) then
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                  KP = IMAX
               else
!
!                 interchange rows and columns K+1 and IMAX, use 2-by-2
!                 pivot block
!
                  KP = IMAX
                  KSTEP = 2
               end if
            end if
!
            KK = K + KSTEP - 1
            if ( KP.NE.KK ) then
!
!              Interchange rows and columns KK and KP in the trailing
!              submatrix A(k:n,k:n)
!
               if ( KP < N )  &
                  call SSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )
               call SSWAP( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ), LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               if ( KSTEP == 2 ) then
                  T = A( K+1, K )
                  A( K+1, K ) = A( KP, K )
                  A( KP, K ) = T
               end if
            end if
!
!           Update the trailing submatrix
!
            if ( KSTEP == 1 ) then
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = L(k)*D(k)
!
!              where L(k) is the k-th column of L
!
               if ( K < N ) then
!
!                 Perform a rank-1 update of A(k+1:n,k+1:n) as
!
!                 A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/D(k))*W(k)'
!
                  D11 = ONE / A( K, K )
                  call SSYR(UPLO, N-K, -D11, A(K+1,K), 1, A(K+1,K+1), LDA)
!
!                 Store L(k) in column K
!
                  call SSCAL( N-K, D11, A( K+1, K ), 1 )
               end if
            else
!
!              2-by-2 pivot block D(k)
!
               if ( K < N-1 ) then
!
!                 Perform a rank-2 update of A(k+2:n,k+2:n) as
!
!                 A := A - ( (A(k) A(k+1))*D(k)**(-1) ) * (A(k) A(k+1))'
!
!                 where L(k) and L(k+1) are the k-th and (k+1)-th
!                 columns of L
!
                  D21 = A( K+1, K )
                  D11 = A( K+1, K+1 ) / D21
                  D22 = A( K, K ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
!
                  DO 60 J = K + 2, N
!
                     WK = D21*( D11*A( J, K )-A( J, K+1 ) )
                     WKP1 = D21*( D22*A( J, K+1 )-A( J, K ) )
!
                     DO 50 I = J, N
                        A(I,J) = A(I,J) - A(I,K)*WK - A(I,K+1)*WKP1
   50                CONTINUE
!
                     A( J, K ) = WK
                     A( J, K+1 ) = WKP1
!
   60             CONTINUE
               end if
            end if
         end if
!
!        Store details of the interchanges in IPIV
!
         if ( KSTEP == 1 ) then
            IPIV( K ) = KP
         else
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         end if
!
!        Increase K and return to the start of the main loop
!
         K = K + KSTEP
         GO TO 40
!
      end if
!
   70 CONTINUE
!
      RETURN
!
!     End of SSYTF2
!
      end subroutine SSYTF2
      subroutine SSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!
!     .. Scalar Arguments ..
      character          UPLO
      integer ::            INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      integer ::            IPIV( * )
      real ::               A( LDA, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  SSYTRF computes the factorization of a real symmetric matrix A using
!  the Bunch-Kaufman diagonal pivoting method.  The form of the
!  factorization is
!
!     A = U*D*U**T  or  A = L*D*L**T
!
!  where U (or L) is a product of permutation and unit upper (lower)
!  triangular matrices, and D is symmetric and block diagonal with 
!  1-by-1 and 2-by-2 diagonal blocks.
!
!  This is the blocked version of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) character*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) integer ::
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) real :: array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, the block diagonal matrix D and the multipliers used
!          to obtain the factor U or L (see below for further details).
!
!  LDA     (input) integer ::
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (output) integer :: array, dimension (N)
!          Details of the interchanges and the block structure of D.
!          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!          interchanged and D(k,k) is a 1-by-1 diagonal block.
!          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!
!  WORK    (workspace/output) real :: array, dimension (max(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) integer ::
!          The length of WORK.  LWORK >=1.  For best performance
!          LWORK >= N*NB, where NB is the block size returned by ILAENV.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) integer ::
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
!                has been completed, but the block diagonal matrix D is
!                exactly singular, and division by zero will occur if it
!                is used to solve a system of equations.
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', then A = U*D*U', where
!     U = P(n)*U(n)* ... *P(k)U(k)* ...,
!  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!
!             (   I    v    0   )   k-s
!     U(k) =  (   0    I    0   )   s
!             (   0    0    I   )   n-k
!                k-s   s   n-k
!
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!  and A(k,k), and v overwrites A(1:k-2,k-1:k).
!
!  If UPLO = 'L', then A = L*D*L', where
!     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
!  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
!  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!
!             (   I    0     0   )  k-1
!     L(k) =  (   0    I     0   )  s
!             (   0    v     I   )  n-k-s+1
!                k-1   s  n-k-s+1
!
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!
!  =====================================================================
!
!     .. Local Scalars ..
      logical ::            LQUERY, UPPER
      integer ::            IINFO, IWS, J, K, KB, LDWORK, LWKOPT, NB, NBMIN
!     ..
!     .. External Functions ..
!m       logical ::            LSAME
!m       integer ::            ILAENV
!m       EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
!m       EXTERNAL           SLASYF, SSYTF2, XERBLA
!     ..
!     .. Intrinsic Functions ..
!m       INTRINSIC          max
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK == -1 )
      if ( .NOT.UPPER .and. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( LDA < max( 1, N ) ) then
         INFO = -4
      else if ( LWORK < 1 .and. .NOT.LQUERY ) then
         INFO = -7
      end if
!
      if ( INFO == 0 ) then
!
!        Determine the block size
!
         NB = ILAENV( 1, 'SSYTRF', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      end if
!
      if ( INFO.NE.0 ) then
         call XERBLA( 'SSYTRF', -INFO )
         RETURN
      else if ( LQUERY ) then
         RETURN
      end if
!
      NBMIN = 2
      LDWORK = N
      if ( NB > 1 .and. NB < N ) then
         IWS = LDWORK*NB
         if ( LWORK < IWS ) then
            NB = max( LWORK / LDWORK, 1 )
            NBMIN = max( 2, ILAENV( 2, 'SSYTRF', UPLO, N, -1, -1, -1 ) )
         end if
      else
         IWS = 1
      end if
      if ( NB < NBMIN ) NB = N
!
      if ( UPPER ) then
!
!        Factorize A as U*D*U' using the upper triangle of A
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        KB, where KB is the number of columns factorized by SLASYF;
!        KB is either NB or NB-1, or K for the last block
!
         K = N
   10    CONTINUE
!
!        If K < 1, exit from loop
!
         if ( K < 1 ) GO TO 40
!
         if ( K > NB ) then
!
!           Factorize columns k-kb+1:k of A and use blocked code to
!           update columns 1:k-kb
!
            call SLASYF( UPLO, K, NB, KB, A, LDA, IPIV, WORK, LDWORK, &
                         IINFO )
         else
!
!           Use unblocked code to factorize columns 1:k of A
!
            call SSYTF2( UPLO, K, A, LDA, IPIV, IINFO )
            KB = K
         end if
!
!        Set INFO on the first occurrence of a zero pivot
!
         if ( INFO == 0 .and. IINFO > 0 ) INFO = IINFO
!
!        Decrease K and return to the start of the main loop
!
         K = K - KB
         GO TO 10
!
      else
!
!        Factorize A as L*D*L' using the lower triangle of A
!
!        K is the main loop index, increasing from 1 to N in steps of
!        KB, where KB is the number of columns factorized by SLASYF;
!        KB is either NB or NB-1, or N-K+1 for the last block
!
         K = 1
   20    CONTINUE
!
!        If K > N, exit from loop
!
         if ( K > N ) GO TO 40
!
         if ( K <= N-NB ) then
!
!           Factorize columns k:k+kb-1 of A and use blocked code to
!           update columns k+kb:n
!
            call SLASYF( UPLO, N-K+1, NB, KB, A( K, K ), LDA, IPIV( K ), &
                         WORK, LDWORK, IINFO )
         else
!
!           Use unblocked code to factorize columns k:n of A
!
            call SSYTF2( UPLO, N-K+1, A( K, K ), LDA, IPIV( K ), IINFO )
            KB = N - K + 1
         end if
!
!        Set INFO on the first occurrence of a zero pivot
!
         if ( INFO == 0 .and. IINFO > 0 ) INFO = IINFO + K - 1
!
!        Adjust IPIV
!
         DO 30 J = K, K + KB - 1
            if ( IPIV( J ) > 0 ) then
               IPIV( J ) = IPIV( J ) + K - 1
            else
               IPIV( J ) = IPIV( J ) - K + 1
            end if
   30    CONTINUE
!
!        Increase K and return to the start of the main loop
!
         K = K + KB
         GO TO 20
!
      end if
!
   40 CONTINUE
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of SSYTRF
!
      end subroutine SSYTRF
      subroutine SSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!
!     .. Scalar Arguments ..
      character          UPLO
      integer ::            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      integer ::            IPIV( * )
      real ::               A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  SSYTRS solves a system of linear equations A*X = B with a real
!  symmetric matrix A using the factorization A = U*D*U**T or
!  A = L*D*L**T computed by SSYTRF.
!
!  Arguments
!  =========
!
!  UPLO    (input) character*1
!          Specifies whether the details of the factorization are stored
!          as an upper or lower triangular matrix.
!          = 'U':  Upper triangular, form is A = U*D*U**T;
!          = 'L':  Lower triangular, form is A = L*D*L**T.
!
!  N       (input) integer ::
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) integer ::
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) real :: array, dimension (LDA,N)
!          The block diagonal matrix D and the multipliers used to
!          obtain the factor U or L as computed by SSYTRF.
!
!  LDA     (input) integer ::
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) integer :: array, dimension (N)
!          Details of the interchanges and the block structure of D
!          as determined by SSYTRF.
!
!  B       (input/output) real :: array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) integer ::
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) integer ::
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      real ::               ONE
      PARAMETER          ( ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      logical ::            UPPER
      integer ::            J, K, KP
      real ::               AK, AKM1, AKM1K, BK, BKM1, DENOM
!     ..
!     .. External Functions ..
!m       logical ::            LSAME
!m       EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!m       EXTERNAL           SGEMV, SGER, SSCAL, SSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
!m       INTRINSIC          max
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .and. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( NRHS < 0 ) then
         INFO = -3
      else if ( LDA < max( 1, N ) ) then
         INFO = -5
      else if ( LDB < max( 1, N ) ) then
         INFO = -8
      end if
      if ( INFO.NE.0 ) then
         call XERBLA( 'SSYTRS', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 .or. NRHS == 0 ) RETURN
!
      if ( UPPER ) then
!
!        Solve A*X = B, where A = U*D*U'.
!
!        First solve U*D*X = B, overwriting B with X.
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = N
   10    CONTINUE
!
!        If K < 1, exit from loop.
!
         if ( K < 1 ) GO TO 30
!
         if ( IPIV( K ) > 0 ) then
!
!           1 x 1 diagonal block
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            if ( KP.NE.K ) &
               call SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in column K of A.
!
            call SGER( K-1, NRHS, -ONE, A( 1, K ), 1, B( K, 1 ), LDB, &
                       B( 1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
            call SSCAL( NRHS, ONE / A( K, K ), B( K, 1 ), LDB )
            K = K - 1
         else
!
!           2 x 2 diagonal block
!
!           Interchange rows K-1 and -IPIV(K).
!
            KP = -IPIV( K )
            if ( KP.NE.K-1 ) &
               call SSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in columns K-1 and K of A.
!
            call SGER( K-2, NRHS, -ONE, A( 1, K ), 1, B( K, 1 ), LDB,  &
                       B( 1, 1 ), LDB )
            call SGER( K-2, NRHS, -ONE, A( 1, K-1 ), 1, B( K-1, 1 ),   &
                       LDB, B( 1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
            AKM1K = A( K-1, K )
            AKM1 = A( K-1, K-1 ) / AKM1K
            AK = A( K, K ) / AKM1K
            DENOM = AKM1*AK - ONE
            DO 20 J = 1, NRHS
               BKM1 = B( K-1, J ) / AKM1K
               BK = B( K, J ) / AKM1K
               B( K-1, J ) = ( AK*BKM1-BK ) / DENOM
               B( K, J ) = ( AKM1*BK-BKM1 ) / DENOM
   20       CONTINUE
            K = K - 2
         end if
!
         GO TO 10
   30    CONTINUE
!
!        Next solve U'*X = B, overwriting B with X.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = 1
   40    CONTINUE
!
!        If K > N, exit from loop.
!
         if ( K > N ) GO TO 50
!
         if ( IPIV( K ) > 0 ) then
!
!           1 x 1 diagonal block
!
!           Multiply by inv(U'(K)), where U(K) is the transformation
!           stored in column K of A.
!
            call SGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB, A( 1, K ),  &
                        1, ONE, B( K, 1 ), LDB )
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            if ( KP.NE.K ) &
               call SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K + 1
         else
!
!           2 x 2 diagonal block
!
!           Multiply by inv(U'(K+1)), where U(K+1) is the transformation
!           stored in columns K and K+1 of A.
!
            call SGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB, A( 1, K ),  &
                        1, ONE, B( K, 1 ), LDB )
            call SGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB,             &
                        A( 1, K+1 ), 1, ONE, B( K+1, 1 ), LDB )
!
!           Interchange rows K and -IPIV(K).
!
            KP = -IPIV( K )
            if ( KP.NE.K ) &
               call SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K + 2
         end if
!
         GO TO 40
   50    CONTINUE
!
      else
!
!        Solve A*X = B, where A = L*D*L'.
!
!        First solve L*D*X = B, overwriting B with X.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = 1
   60    CONTINUE
!
!        If K > N, exit from loop.
!
         if ( K > N ) GO TO 80
!
         if ( IPIV( K ) > 0 ) then
!
!           1 x 1 diagonal block
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            if ( KP.NE.K ) &
               call SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in column K of A.
!
            if ( K < N ) &
               call SGER( N-K, NRHS, -ONE, A( K+1, K ), 1, B( K, 1 ),  &
                          LDB, B( K+1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
            call SSCAL( NRHS, ONE / A( K, K ), B( K, 1 ), LDB )
            K = K + 1
         else
!
!           2 x 2 diagonal block
!
!           Interchange rows K+1 and -IPIV(K).
!
            KP = -IPIV( K )
            if ( KP.NE.K+1 )  &
               call SSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in columns K and K+1 of A.
!
            if ( K < N-1 ) then
               call SGER( N-K-1, NRHS, -ONE, A( K+2, K ), 1, B( K, 1 ),  &
                          LDB, B( K+2, 1 ), LDB )
               call SGER( N-K-1, NRHS, -ONE, A( K+2, K+1 ), 1,           &
                          B( K+1, 1 ), LDB, B( K+2, 1 ), LDB )
            end if
!
!           Multiply by the inverse of the diagonal block.
!
            AKM1K = A( K+1, K )
            AKM1 = A( K, K ) / AKM1K
            AK = A( K+1, K+1 ) / AKM1K
            DENOM = AKM1*AK - ONE
            DO 70 J = 1, NRHS
               BKM1 = B( K, J ) / AKM1K
               BK = B( K+1, J ) / AKM1K
               B( K, J ) = ( AK*BKM1-BK ) / DENOM
               B( K+1, J ) = ( AKM1*BK-BKM1 ) / DENOM
   70       CONTINUE
            K = K + 2
         end if
!
         GO TO 60
   80    CONTINUE
!
!        Next solve L'*X = B, overwriting B with X.
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = N
   90    CONTINUE
!
!        If K < 1, exit from loop.
!
         if ( K < 1 ) GO TO 100
!
         if ( IPIV( K ) > 0 ) then
!
!           1 x 1 diagonal block
!
!           Multiply by inv(L'(K)), where L(K) is the transformation
!           stored in column K of A.
!
            if ( K < N )  &
               call SGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ),  &
                           LDB, A( K+1, K ), 1, ONE, B( K, 1 ), LDB )
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            if ( KP.NE.K ) &
               call SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K - 1
         else
!
!           2 x 2 diagonal block
!
!           Multiply by inv(L'(K-1)), where L(K-1) is the transformation
!           stored in columns K-1 and K of A.
!
            if ( K < N ) then
               call SGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ),  &
                           LDB, A( K+1, K ), 1, ONE, B( K, 1 ), LDB )
               call SGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ),  &
                           LDB, A( K+1, K-1 ), 1, ONE, B( K-1, 1 ),    &
                           LDB )
            end if
!
!           Interchange rows K and -IPIV(K).
!
            KP = -IPIV( K )
            if ( KP.NE.K ) &
               call SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K - 2
         end if
!
         GO TO 90
  100    CONTINUE
      end if
!
      RETURN
!
!     End of SSYTRS
!
      end subroutine SSYTRS

























!---------------------------------------------------------------------
      logical function LSAME(CA,CB)
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      character CA,CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) character*1
!
!  CB      (input) character*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC ICHAR
!     ..
!     .. Local Scalars ..
      integer :: INTA,INTB,ZCODE
!     ..
!
!     Test if the characters are equal
!
      LSAME = CA  ==  CB
      if (LSAME) RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR('Z')
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR(CA)
      INTB = ICHAR(CB)
!
      if (ZCODE == 90 .or. ZCODE == 122) then
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
          if (INTA >= 97 .and. INTA <= 122) INTA = INTA - 32
          if (INTB >= 97 .and. INTB <= 122) INTB = INTB - 32
!
      else if (ZCODE == 233 .or. ZCODE == 169) then
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
          if (INTA >= 129 .and. INTA <= 137 .or.  &
              INTA >= 145 .and. INTA <= 153 .or.  &
              INTA >= 162 .and. INTA <= 169) INTA = INTA + 64
          if (INTB >= 129 .and. INTB <= 137 .or.  &
              INTB >= 145 .and. INTB <= 153 .or.  &
              INTB >= 162 .and. INTB <= 169) INTB = INTB + 64
!
      else if (ZCODE == 218 .or. ZCODE == 250) then
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
          if (INTA >= 225 .and. INTA <= 250) INTA = INTA - 32
          if (INTB >= 225 .and. INTB <= 250) INTB = INTB - 32
      end if
      LSAME = INTA  ==  INTB
!
!     RETURN
!
!     End of LSAME
!
      end function LSAME




!---------------------------------------------------------------------
      integer function IEEECK( ISPEC, ZERO, ONE )
!
!  -- LAPACK auxiliary routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!
!     .. Scalar Arguments ..
      integer ::            ISPEC
      real ::               ONE, ZERO
!     ..
!
!  Purpose
!  =======
!
!  IEEECK is called from the ILAENV to verify that Infinity and
!  possibly NaN arithmetic is safe (i.e. will not trap).
!
!  Arguments
!  =========
!
!  ISPEC   (input) integer ::
!          Specifies whether to test just for inifinity arithmetic
!          or whether to test for infinity and NaN arithmetic.
!          = 0: Verify infinity arithmetic only.
!          = 1: Verify infinity and NaN arithmetic.
!
!  ZERO    (input) real ::
!          Must contain the value 0.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  ONE     (input) real ::
!          Must contain the value 1.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  RETURN VALUE:  integer ::
!          = 0:  Arithmetic failed to produce the correct answers
!          = 1:  Arithmetic produced the correct answers
!
!     .. Local Scalars ..
      real ::   NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF
      real ::   NEGZRO, NEWZRO, POSINF
!     ..
!     .. Executable Statements ..
      IEEECK = 1
!
      POSINF = ONE / ZERO
      if ( POSINF <= ONE ) then
         IEEECK = 0
         RETURN
      end if
!
      NEGINF = -ONE / ZERO
      if ( NEGINF >= ZERO ) then
         IEEECK = 0
         RETURN
      end if
!
      NEGZRO = ONE / ( NEGINF+ONE )
      if ( NEGZRO.NE.ZERO ) then
         IEEECK = 0
         RETURN
      end if
!
      NEGINF = ONE / NEGZRO
      if ( NEGINF >= ZERO ) then
         IEEECK = 0
         RETURN
      end if
!
      NEWZRO = NEGZRO + ZERO
      if ( NEWZRO.NE.ZERO ) then
         IEEECK = 0
         RETURN
      end if
!
      POSINF = ONE / NEWZRO
      if ( POSINF <= ONE ) then
         IEEECK = 0
         RETURN
      end if
!
      NEGINF = NEGINF*POSINF
      if ( NEGINF >= ZERO ) then
         IEEECK = 0
         RETURN
      end if
!
      POSINF = POSINF*POSINF
      if ( POSINF <= ONE ) then
         IEEECK = 0
         RETURN
      end if
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
      if ( ISPEC == 0 )  RETURN
!
      NAN1 = POSINF + NEGINF
!
      NAN2 = POSINF / NEGINF
!
      NAN3 = POSINF / POSINF
!
      NAN4 = POSINF*ZERO
!
      NAN5 = NEGINF*NEGZRO
!
      NAN6 = NAN5*0.0
!
      if ( NAN1 == NAN1 ) then
         IEEECK = 0
         RETURN
      end if
!
      if ( NAN2 == NAN2 ) then
         IEEECK = 0
         RETURN
      end if
!
      if ( NAN3 == NAN3 ) then
         IEEECK = 0
         RETURN
      end if
!
      if ( NAN4 == NAN4 ) then
         IEEECK = 0
         RETURN
      end if
!
      if ( NAN5 == NAN5 ) then
         IEEECK = 0
         RETURN
      end if
!
      if ( NAN6 == NAN6 ) then
         IEEECK = 0
         RETURN
      end if
!
      RETURN
      end function IEEECK
 

!---------------------------------------------------------------------
      integer function ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!
!  -- LAPACK auxiliary routine (version 3.2.1)                        --
!
!  -- April 2009                                                      --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      character*( * )    NAME, OPTS
      integer ::         ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  ILAENV returns an integer ::
!  if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
!  if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) integer ::
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines (DEPRECATED)
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR method
!               for nonsymmetric eigenvalue problems (DEPRECATED)
!          = 9: maximum size of the subproblems at the bottom of the
!               computation tree in the divide-and-conquer algorithm
!               (used by xGELSD and xGESDD)
!          =10: ieee NaN arithmetic can be trusted not to trap
!          =11: infinity arithmetic can be trusted not to trap
!          12 <= ISPEC <= 16:
!               xHSEQR or one of its subroutines,
!               see IPARMQ for detailed explanation
!
!  NAME    (input) character*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) character*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) integer ::
!  N2      (input) integer ::
!  N3      (input) integer ::
!  N4      (input) integer ::
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      if ( NB <= 1 ) NB = max( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
      integer ::            I, IC, IZ, NB, NBMIN, NX
      logical ::            CNAME, SNAME
      character          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, min, real
!     ..
!     .. External Functions ..
!m       integer ::            IEEECK, IPARMQ
!m       EXTERNAL           IEEECK, IPARMQ
!     ..
!     .. Executable Statements ..
!
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120,                 &
              130, 140, 150, 160, 160, 160, 160, 160 ) ISPEC
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
   10 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      if ( IZ == 90 .or. IZ == 122 ) then
!
!        ASCII character set
!
         if ( IC >= 97 .and. IC <= 122 ) then
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               if ( IC >= 97 .and. IC <= 122 ) SUBNAM( I: I ) = CHAR( IC-32 )
   20       CONTINUE
         end if
!
      else if ( IZ == 233 .or. IZ == 169 ) then
!
!        EBCDIC character set
!
         if ( ( IC >= 129 .and. IC <= 137 ) .or.   &
             ( IC >= 145 .and. IC <= 153 ) .or.   &
             ( IC >= 162 .and. IC <= 169 ) ) then
            SUBNAM( 1: 1 ) = CHAR( IC+64 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               if ( ( IC >= 129 .and. IC <= 137 ) .or.  &
                   ( IC >= 145 .and. IC <= 153 ) .or.  &
                   ( IC >= 162 .and. IC <= 169 ) ) SUBNAM(I:I) = CHAR(IC+64)
   30       CONTINUE
         end if
!
      else if ( IZ == 218 .or. IZ == 250 ) then
!
!        Prime machines:  ASCII+128
!
         if ( IC >= 225 .and. IC <= 250 ) then
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               if ( IC >= 225 .and. IC <= 250 )  SUBNAM(I:I) = CHAR(IC-32)
   40       CONTINUE
         end if
      end if
!
      C1 = SUBNAM( 1: 1 )
      SNAME = C1 == 'S' .or. C1 == 'D'
      CNAME = C1 == 'C' .or. C1 == 'Z'
      if ( .NOT.( CNAME .or. SNAME ) ) RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
!
      GO TO ( 50, 60, 70 )ISPEC
!
   50 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      if ( C2 == 'GE' ) then
         if ( C3 == 'TRF' ) then
            if ( SNAME ) then
               NB = 64
            else
               NB = 64
            end if
         else if ( C3 == 'QRF' .or. C3 == 'RQF' .or. C3 == 'LQF' .or.  &
                  C3 == 'QLF' ) then
            if ( SNAME ) then
               NB = 32
            else
               NB = 32
            end if
         else if ( C3 == 'HRD' ) then
            if ( SNAME ) then
               NB = 32
            else
               NB = 32
            end if
         else if ( C3 == 'BRD' ) then
            if ( SNAME ) then
               NB = 32
            else
               NB = 32
            end if
         else if ( C3 == 'TRI' ) then
            if ( SNAME ) then
               NB = 64
            else
               NB = 64
            end if
         end if
      else if ( C2 == 'PO' ) then
         if ( C3 == 'TRF' ) then
            if ( SNAME ) then
               NB = 64
            else
               NB = 64
            end if
         end if
      else if ( C2 == 'SY' ) then
         if ( C3 == 'TRF' ) then
            if ( SNAME ) then
               NB = 64
            else
               NB = 64
            end if
         else if ( SNAME .and. C3 == 'TRD' ) then
            NB = 32
         else if ( SNAME .and. C3 == 'GST' ) then
            NB = 64
         end if
      else if ( CNAME .and. C2 == 'HE' ) then
         if ( C3 == 'TRF' ) then
            NB = 64
         else if ( C3 == 'TRD' ) then
            NB = 32
         else if ( C3 == 'GST' ) then
            NB = 64
         end if
      else if ( SNAME .and. C2 == 'OR' ) then
         if ( C3( 1: 1 ) == 'G' ) then
            if ( C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 ==    &
                'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR' )   &
                 then
               NB = 32
            end if
         else if ( C3( 1: 1 ) == 'M' ) then
            if ( C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 ==    &
                'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR' )   &
                 then
               NB = 32
            end if
         end if
      else if ( CNAME .and. C2 == 'UN' ) then
         if ( C3( 1: 1 ) == 'G' ) then
            if ( C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 ==   &
                'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR' )  &
                 then
               NB = 32
            end if
         else if ( C3( 1: 1 ) == 'M' ) then
            if ( C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 ==   &
                'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR' )  &
                 then
               NB = 32
            end if
         end if
      else if ( C2 == 'GB' ) then
         if ( C3 == 'TRF' ) then
            if ( SNAME ) then
               if ( N4 <= 64 ) then
                  NB = 1
               else
                  NB = 32
               end if
            else
               if ( N4 <= 64 ) then
                  NB = 1
               else
                  NB = 32
               end if
            end if
         end if
      else if ( C2 == 'PB' ) then
         if ( C3 == 'TRF' ) then
            if ( SNAME ) then
               if ( N2 <= 64 ) then
                  NB = 1
               else
                  NB = 32
               end if
            else
               if ( N2 <= 64 ) then
                  NB = 1
               else
                  NB = 32
               end if
            end if
         end if
      else if ( C2 == 'TR' ) then
         if ( C3 == 'TRI' ) then
            if ( SNAME ) then
               NB = 64
            else
               NB = 64
            end if
         end if
      else if ( C2 == 'LA' ) then
         if ( C3 == 'UUM' ) then
            if ( SNAME ) then
               NB = 64
            else
               NB = 64
            end if
         end if
      else if ( SNAME .and. C2 == 'ST' ) then
         if ( C3 == 'EBZ' ) then
            NB = 1
         end if
      end if
      ILAENV = NB
      RETURN
!
   60 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      if ( C2 == 'GE' ) then
         if ( C3 == 'QRF' .or. C3 == 'RQF' .or. C3 == 'LQF' .or. C3 == 'QLF' ) then
            if ( SNAME ) then
               NBMIN = 2
            else
               NBMIN = 2
            end if
         else if ( C3 == 'HRD' ) then
            if ( SNAME ) then
               NBMIN = 2
            else
               NBMIN = 2
            end if
         else if ( C3 == 'BRD' ) then
            if ( SNAME ) then
               NBMIN = 2
            else
               NBMIN = 2
            end if
         else if ( C3 == 'TRI' ) then
            if ( SNAME ) then
               NBMIN = 2
            else
               NBMIN = 2
            end if
         end if
      else if ( C2 == 'SY' ) then
         if ( C3 == 'TRF' ) then
            if ( SNAME ) then
               NBMIN = 8
            else
               NBMIN = 8
            end if
         else if ( SNAME .and. C3 == 'TRD' ) then
            NBMIN = 2
         end if
      else if ( CNAME .and. C2 == 'HE' ) then
         if ( C3 == 'TRD' ) then
            NBMIN = 2
         end if
      else if ( SNAME .and. C2 == 'OR' ) then
         if ( C3( 1: 1 ) == 'G' ) then
            if ( C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 ==   &  
                'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR' )  &
                 then
               NBMIN = 2
            end if
         else if ( C3( 1: 1 ) == 'M' ) then
            if ( C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 ==   &
                'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR' )  &
                 then
               NBMIN = 2
            end if
         end if
      else if ( CNAME .and. C2 == 'UN' ) then
         if ( C3( 1: 1 ) == 'G' ) then
            if ( C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 ==   &
                'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR' )  &
                 then
               NBMIN = 2
            end if
         else if ( C3( 1: 1 ) == 'M' ) then
            if ( C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 ==   &
                'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR' )  &
                 then
               NBMIN = 2
            end if
         end if
      end if
      ILAENV = NBMIN
      RETURN
!
   70 CONTINUE
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      if ( C2 == 'GE' ) then
         if ( C3 == 'QRF' .or. C3 == 'RQF' .or. C3 == 'LQF' .or. C3 ==   &
             'QLF' ) then
            if ( SNAME ) then
               NX = 128
            else
               NX = 128
            end if
         else if ( C3 == 'HRD' ) then
            if ( SNAME ) then
               NX = 128
            else
               NX = 128
            end if
         else if ( C3 == 'BRD' ) then
            if ( SNAME ) then
               NX = 128
            else
               NX = 128
            end if
         end if
      else if ( C2 == 'SY' ) then
         if ( SNAME .and. C3 == 'TRD' ) then
            NX = 32
         end if
      else if ( CNAME .and. C2 == 'HE' ) then
         if ( C3 == 'TRD' ) then
            NX = 32
         end if
      else if ( SNAME .and. C2 == 'OR' ) then
         if ( C3( 1: 1 ) == 'G' ) then
            if ( C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 ==   &
                'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR' )  &
                 then
               NX = 128
            end if
         end if
      else if ( CNAME .and. C2 == 'UN' ) then
         if ( C3( 1: 1 ) == 'G' ) then
            if ( C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 ==   &
                'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR' )  &
                 then
               NX = 128
            end if
         end if
      end if
      ILAENV = NX
      RETURN
!
   80 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
   90 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!
  100 CONTINUE
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( real( min( N1, N2 ) )*1.6E0 )
      RETURN
!
  110 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  120 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
  130 CONTINUE
!
!     ISPEC = 9:  maximum size of the subproblems at the bottom of the
!                 computation tree in the divide-and-conquer algorithm
!                 (used by xGELSD and xGESDD)
!
      ILAENV = 25
      RETURN
!
  140 CONTINUE
!
!     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
!
!     ILAENV = 0
      ILAENV = 1
      if ( ILAENV == 1 ) then
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      end if
      RETURN
!
  150 CONTINUE
!
!     ISPEC = 11: infinity arithmetic can be trusted not to trap
!
!     ILAENV = 0
      ILAENV = 1
      if ( ILAENV == 1 ) then
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      end if
      RETURN
!
  160 CONTINUE
!
!     12 <= ISPEC <= 16: xHSEQR or one of its subroutines. 
!
      ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      RETURN
!
!     End of ILAENV
!
      end function ILAENV


!---------------------------------------------------------------------
      integer function IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
!
!  -- LAPACK auxiliary routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!     
!     .. Scalar Arguments ..
      integer ::            IHI, ILO, ISPEC, LWORK, N
      character          NAME*( * ), OPTS*( * )
!
!  Purpose
!  =======
!
!       This program sets problem and machine dependent parameters
!       useful for xHSEQR and its subroutines. It is called whenever 
!       ILAENV is called with 12 <= ISPEC <= 16
!
!  Arguments
!  =========
!
!       ISPEC  (input) integer scalar
!              ISPEC specifies which tunable parameter IPARMQ should
!              return.
!
!              ISPEC=12: (INMIN)  Matrices of order nmin or less
!                        are sent directly to xLAHQR, the implicit
!                        double shift QR algorithm.  NMIN must be
!                        at least 11.
!
!              ISPEC=13: (INWIN)  Size of the deflation window.
!                        This is best set greater than or equal to
!                        the number of simultaneous shifts NS.
!                        Larger matrices benefit from larger deflation
!                        windows.
!
!              ISPEC=14: (INIBL) Determines when to stop nibbling and
!                        invest in an (expensive) multi-shift QR sweep.
!                        If the aggressive early deflation subroutine
!                        finds LD converged eigenvalues from an order
!                        NW deflation window and LD > (NW*NIBBLE)/100,
!                        then the next QR sweep is skipped and early
!                        deflation is applied immediately to the
!                        remaining active diagonal block.  Setting
!                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
!                        multi-shift QR sweep whenever early deflation
!                        finds a converged eigenvalue.  Setting
!                        IPARMQ(ISPEC=14) greater than or equal to 100
!                        prevents TTQRE from skipping a multi-shift
!                        QR sweep.
!
!              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
!                        a multi-shift QR iteration.
!
!              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
!                        following meanings.
!                        0:  During the multi-shift QR sweep,
!                            xLAQR5 does not accumulate reflections and
!                            does not use matrix-matrix multiply to
!                            update the far-from-diagonal matrix
!                            entries.
!                        1:  During the multi-shift QR sweep,
!                            xLAQR5 and/or xLAQRaccumulates reflections and uses
!                            matrix-matrix multiply to update the
!                            far-from-diagonal matrix entries.
!                        2:  During the multi-shift QR sweep.
!                            xLAQR5 accumulates reflections and takes
!                            advantage of 2-by-2 block structure during
!                            matrix-matrix multiplies.
!                        (If xTRMM is slower than xGEMM, then
!                        IPARMQ(ISPEC=16)=1 may be more efficient than
!                        IPARMQ(ISPEC=16)=2 despite the greater level of
!                        arithmetic work implied by the latter choice.)
!
!       NAME    (input) character string
!               Name of the calling subroutine
!
!       OPTS    (input) character string
!               This is a concatenation of the string arguments to
!               TTQRE.
!
!       N       (input) integer scalar
!               N is the order of the Hessenberg matrix H.
!
!       ILO     (input) integer ::
!       IHI     (input) integer ::
!               It is assumed that H is already upper triangular
!               in rows and columns 1:ILO-1 and IHI+1:N.
!
!       LWORK   (input) integer scalar
!               The amount of workspace available.
!
!  Further Details
!  ===============
!
!       Little is known about how best to choose these parameters.
!       It is possible to use different values of the parameters
!       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
!
!       It is probably best to choose different parameters for
!       different matrices and different parameters at different
!       times during the iteration, but this has not been
!       implemented --- yet.
!
!
!       The best choices of most of the parameters depend
!       in an ill-understood way on the relative execution
!       rate of xLAQR3 and xLAQR5 and on the nature of each
!       particular eigenvalue problem.  Experiment may be the
!       only practical way to determine which choices are most
!       effective.
!
!       Following is a list of default values supplied by IPARMQ.
!       These defaults may be adjusted in order to attain better
!       performance in any particular computational environment.
!
!       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
!                        Default: 75. (Must be at least 11.)
!
!       IPARMQ(ISPEC=13) Recommended deflation window size.
!                        This depends on ILO, IHI and NS, the
!                        number of simultaneous shifts returned
!                        by IPARMQ(ISPEC=15).  The default for
!                        (IHI-ILO+1) <= 500 is NS.  The default
!                        for (IHI-ILO+1) > 500 is 3*NS/2.
!
!       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
!
!       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
!                        a multi-shift QR iteration.
!
!                        If IHI-ILO+1 is ...
!
!                        greater than      ...but less    ... the
!                        or equal to ...      than        default is
!
!                                0               30       NS =   2+
!                               30               60       NS =   4+
!                               60              150       NS =  10
!                              150              590       NS =  **
!                              590             3000       NS =  64
!                             3000             6000       NS = 128
!                             6000             infinity   NS = 256
!
!                    (+)  By default matrices of this order are
!                         passed to the implicit double shift routine
!                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
!                         values of NS are used only in case of a rare
!                         xLAHQR failure.
!
!                    (**) The asterisks (**) indicate an ad-hoc
!                         function increasing from 10 to 64.
!
!       IPARMQ(ISPEC=16) Select structured matrix multiply.
!                        (See ISPEC=16 above for details.)
!                        Default: 3.
!
!     ================================================================
!     .. Parameters ..
      integer ::            INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14 )
      PARAMETER          ( ISHFTS = 15, IACC22 = 16 )
      integer ::            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14 )
      PARAMETER           ( NIBBLE = 14, KNWSWP = 500 )
      real ::               TWO
      PARAMETER          ( TWO = 2.0 )
!     ..
!     .. Local Scalars ..
      integer ::            NH, NS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          LOG, max, MOD, NINT, real 
!     ..
!     .. Executable Statements ..
      if ( ( ISPEC == ISHFTS ) .or. ( ISPEC == INWIN ) .or. ( ISPEC == IACC22 ) ) then
!
!        ==== Set the number simultaneous shifts ====
!
         NH = IHI - ILO + 1
         NS = 2
         if ( NH >= 30 ) NS = 4
         if ( NH >= 60 ) NS = 10
         if ( NH >= 150 ) &
            NS = max( 10, NH / NINT( LOG( real( NH ) ) / LOG( TWO ) ) )
         if ( NH >= 590 ) NS = 64
         if ( NH >= 3000 ) NS = 128
         if ( NH >= 6000 )  NS = 256
         NS = max( 2, NS-MOD( NS, 2 ) )
      end if
!
      if ( ISPEC == INMIN ) then
!
!
!        ===== Matrices of order smaller than NMIN get sent
!        .     to xLAHQR, the classic double shift algorithm.
!        .     This must be at least 11. ====
!
         IPARMQ = NMIN
!
      else if ( ISPEC == INIBL ) then
!
!        ==== INIBL: skip a multi-shift qr iteration and
!        .    whenever aggressive early deflation finds
!        .    at least (NIBBLE*(window size)/100) deflations. ====
!
         IPARMQ = NIBBLE
!
      else if ( ISPEC == ISHFTS ) then
!
!        ==== NSHFTS: The number of simultaneous shifts =====
!
         IPARMQ = NS
!
      else if ( ISPEC == INWIN ) then
!
!        ==== NW: deflation window size.  ====
!
         if ( NH <= KNWSWP ) then
            IPARMQ = NS
         else
            IPARMQ = 3*NS / 2
         end if
!
      else if ( ISPEC == IACC22 ) then
!
!        ==== IACC22: Whether to accumulate reflections
!        .     before updating the far-from-diagonal elements
!        .     and whether to use 2-by-2 block structure while
!        .     doing it.  A small amount of work could be saved
!        .     by making this choice dependent also upon the
!        .     NH=IHI-ILO+1.
!
         IPARMQ = 0
         if ( NS >= KACMIN ) IPARMQ = 1
         if ( NS >= K22MIN ) IPARMQ = 2
!
      else
!        ===== invalid value of ispec =====
         IPARMQ = -1
!
      end if
!
!     ==== End of IPARMQ ====
!
      end function IPARMQ






                             
!----------------------------------------------------------------------
!   lapack_auxiliares
!----------------------------------------------------------------------




!---------------------------------------------------------------------
      subroutine DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!     .. Scalar Arguments ..
      real(8) :: ALPHA,BETA
      integer :: K,LDA,LDB,LDC,M,N
      character TRANSA,TRANSB
!     ..
!     .. Array Arguments ..
      real(8) :: A(LDA,*),B(LDB,*),C(LDC,*)
!     ..
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Arguments
!  ==========
!
!  TRANSA - character*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - character*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - integer ::.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - integer ::.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - integer ::.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - real(8) ::.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - real(8) :: array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - integer ::.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - real(8) :: array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - integer ::.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - real(8) ::.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - real(8) :: array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - integer ::.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
!m       logical :: LSAME
!m      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!m       EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
!m       INTRINSIC max
!     ..
!     .. Local Scalars ..
      real(8) :: TEMP
      integer :: I,INFO,J,L,NCOLA,NROWA,NROWB
      logical :: NOTA,NOTB
!     ..
!     .. Parameters ..
      real(8) :: ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      if (NOTA) then
          NROWA = M
          NCOLA = K
      else
          NROWA = K
          NCOLA = M
      end if
      if (NOTB) then
          NROWB = K
      else
          NROWB = N
      end if
!
!     Test the input parameters.
!
      INFO = 0
      if ((.NOT.NOTA) .and. (.NOT.LSAME(TRANSA,'C')) .and.         &
          (.NOT.LSAME(TRANSA,'T'))) then
          INFO = 1
      else if ((.NOT.NOTB) .and. (.NOT.LSAME(TRANSB,'C')) .and.    &
               (.NOT.LSAME(TRANSB,'T'))) then
          INFO = 2
      else if (M < 0) then
          INFO = 3
      else if (N < 0) then
          INFO = 4
      else if (K < 0) then
          INFO = 5
      else if (LDA < max(1,NROWA)) then
          INFO = 8
      else if (LDB < max(1,NROWB)) then
          INFO = 10
      else if (LDC < max(1,M)) then
          INFO = 13
      end if
      if (INFO.NE.0) then
          call XERBLA('DGEMM ',INFO)
          RETURN
      end if
!
!     Quick return if possible.
!
      if ((M == 0) .or. (N == 0) .or.   &
          (((ALPHA == ZERO).or. (K == 0)).and. (BETA == ONE))) RETURN
!
!     And if  alpha.eq.zero.
!
      if (ALPHA == ZERO) then
          if (BETA == ZERO) then
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          else
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          end if
          RETURN
      end if
!
!     Start the operations.
!
      if (NOTB) then
          if (NOTA) then
!
!           Form  C := alpha*A*B + beta*C.
!
              DO 90 J = 1,N
                  if (BETA == ZERO) then
                      DO 50 I = 1,M
                          C(I,J) = ZERO
   50                 CONTINUE
                  else if (BETA.NE.ONE) then
                      DO 60 I = 1,M
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  end if
                  DO 80 L = 1,K
                      if (B(L,J).NE.ZERO) then
                          TEMP = ALPHA*B(L,J)
                          DO 70 I = 1,M
                              C(I,J) = C(I,J) + TEMP*A(I,L)
   70                     CONTINUE
                      end if
   80             CONTINUE
   90         CONTINUE
          else
!
!           Form  C := alpha*A'*B + beta*C
!
              DO 120 J = 1,N
                  DO 110 I = 1,M
                      TEMP = ZERO
                      DO 100 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  100                 CONTINUE
                      if (BETA == ZERO) then
                          C(I,J) = ALPHA*TEMP
                      else
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      end if
  110             CONTINUE
  120         CONTINUE
          end if
      else
          if (NOTA) then
!
!           Form  C := alpha*A*B' + beta*C
!
              DO 170 J = 1,N
                  if (BETA == ZERO) then
                      DO 130 I = 1,M
                          C(I,J) = ZERO
  130                 CONTINUE
                  else if (BETA.NE.ONE) then
                      DO 140 I = 1,M
                          C(I,J) = BETA*C(I,J)
  140                 CONTINUE
                  end if
                  DO 160 L = 1,K
                      if (B(J,L).NE.ZERO) then
                          TEMP = ALPHA*B(J,L)
                          DO 150 I = 1,M
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  150                     CONTINUE
                      end if
  160             CONTINUE
  170         CONTINUE
          else
!
!           Form  C := alpha*A'*B' + beta*C
!
              DO 200 J = 1,N
                  DO 190 I = 1,M
                      TEMP = ZERO
                      DO 180 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  180                 CONTINUE
                      if (BETA == ZERO) then
                          C(I,J) = ALPHA*TEMP
                      else
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      end if
  190             CONTINUE
  200         CONTINUE
          end if
      end if
!
      RETURN
!
!     End of DGEMM .
!
      end subroutine DGEMM




!---------------------------------------------------------------------
      subroutine DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!     .. Scalar Arguments ..
      real(8) :: ALPHA,BETA
      integer :: INCX,INCY,LDA,M,N
      character TRANS
!     ..
!     .. Array Arguments ..
      real(8) :: A(LDA,*),X(*),Y(*)
!     ..
!
!  Purpose
!  =======
!
!  DGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Arguments
!  ==========
!
!  TRANS  - character*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - integer ::.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - integer ::.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real(8) ::.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - real(8) :: array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - integer ::.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - real(8) :: array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - integer ::.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - real(8) ::.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - real(8) :: array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - integer ::.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      real(8) :: ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      real(8) :: TEMP
      integer :: I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
!     ..
!     .. External Functions ..
!m       logical :: LSAME
!m       EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!m       EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
!m       INTRINSIC max
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      if (.NOT.LSAME(TRANS,'N') .and. .NOT.LSAME(TRANS,'T') .and.  &
          .NOT.LSAME(TRANS,'C')) then
          INFO = 1
      else if (M < 0) then
          INFO = 2
      else if (N < 0) then
          INFO = 3
      else if (LDA < max(1,M)) then
          INFO = 6
      else if (INCX == 0) then
          INFO = 8
      else if (INCY == 0) then
          INFO = 11
      end if
      if (INFO.NE.0) then
          call XERBLA('DGEMV ',INFO)
          RETURN
      end if
!
!     Quick return if possible.
!
      if ((M == 0) .or. (N == 0) .or.  &
          ((ALPHA == ZERO).and. (BETA == ONE))) RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      if (LSAME(TRANS,'N')) then
          LENX = N
          LENY = M
      else
          LENX = M
          LENY = N
      end if
      if (INCX > 0) then
          KX = 1
      else
          KX = 1 - (LENX-1)*INCX
      end if
      if (INCY > 0) then
          KY = 1
      else
          KY = 1 - (LENY-1)*INCY
      end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      if (BETA.NE.ONE) then
          if (INCY == 1) then
              if (BETA == ZERO) then
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              else
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              end if
          else
              IY = KY
              if (BETA == ZERO) then
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              else
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              end if
          end if
      end if
      if (ALPHA == ZERO) RETURN
      if (LSAME(TRANS,'N')) then
!
!        Form  y := alpha*A*x + y.
!
          JX = KX
          if (INCY == 1) then
              DO 60 J = 1,N
                  if (X(JX).NE.ZERO) then
                      TEMP = ALPHA*X(JX)
                      DO 50 I = 1,M
                          Y(I) = Y(I) + TEMP*A(I,J)
   50                 CONTINUE
                  end if
                  JX = JX + INCX
   60         CONTINUE
          else
              DO 80 J = 1,N
                  if (X(JX).NE.ZERO) then
                      TEMP = ALPHA*X(JX)
                      IY = KY
                      DO 70 I = 1,M
                          Y(IY) = Y(IY) + TEMP*A(I,J)
                          IY = IY + INCY
   70                 CONTINUE
                  end if
                  JX = JX + INCX
   80         CONTINUE
          end if
      else
!
!        Form  y := alpha*A'*x + y.
!
          JY = KY
          if (INCX == 1) then
              DO 100 J = 1,N
                  TEMP = ZERO
                  DO 90 I = 1,M
                      TEMP = TEMP + A(I,J)*X(I)
   90             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  100         CONTINUE
          else
              DO 120 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  DO 110 I = 1,M
                      TEMP = TEMP + A(I,J)*X(IX)
                      IX = IX + INCX
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  120         CONTINUE
          end if
      end if
!
      RETURN
!
!     End of DGEMV .
!
      end subroutine DGEMV



!---------------------------------------------------------------------
      subroutine DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!     .. Scalar Arguments ..
      real(8) :: ALPHA
      integer :: INCX,INCY,LDA,M,N
!     ..
!     .. Array Arguments ..
      real(8) :: A(LDA,*),X(*),Y(*)
!     ..
!
!  Purpose
!  =======
!
!  DGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Arguments
!  ==========
!
!  M      - integer ::.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - integer ::.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real(8) ::.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - real(8) :: array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - integer ::.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - real(8) :: array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - integer ::.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - real(8) :: array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - integer ::.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      real(8) :: ZERO
      PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      real(8) :: TEMP
      integer :: I,INFO,IX,J,JY,KX
!     ..
!     .. External Subroutines ..
!m       EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
!m       INTRINSIC max
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      if (M < 0) then
          INFO = 1
      else if (N < 0) then
          INFO = 2
      else if (INCX == 0) then
          INFO = 5
      else if (INCY == 0) then
          INFO = 7
      else if (LDA < max(1,M)) then
          INFO = 9
      end if
      if (INFO.NE.0) then
          call XERBLA('DGER  ',INFO)
          RETURN
      end if
!
!     Quick return if possible.
!
      if ((M == 0) .or. (N == 0) .or. (ALPHA == ZERO)) RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      if (INCY > 0) then
          JY = 1
      else
          JY = 1 - (N-1)*INCY
      end if
      if (INCX == 1) then
          DO 20 J = 1,N
              if (Y(JY).NE.ZERO) then
                  TEMP = ALPHA*Y(JY)
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              end if
              JY = JY + INCY
   20     CONTINUE
      else
          if (INCX > 0) then
              KX = 1
          else
              KX = 1 - (M-1)*INCX
          end if
          DO 40 J = 1,N
              if (Y(JY).NE.ZERO) then
                  TEMP = ALPHA*Y(JY)
                  IX = KX
                  DO 30 I = 1,M
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              end if
              JY = JY + INCY
   40     CONTINUE
      end if
!
      RETURN
!
!     End of DGER  .
!
      end subroutine DGER


!---------------------------------------------------------------------
      subroutine DSCAL(N,DA,DX,INCX)
!     .. Scalar Arguments ..
      real(8) :: DA
      integer :: INCX,N
!     ..
!     .. Array Arguments ..
      real(8) :: DX(*)
!     ..
!
!  Purpose
!  =======
!*
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      integer :: I,M,MP1,NINCX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      if (N <= 0 .or. INCX <= 0) RETURN
      if (INCX == 1) GO TO 20
!
!        code for increment not equal to 1
!
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,5)
      if (M == 0) GO TO 40
      DO 30 I = 1,M
          DX(I) = DA*DX(I)
   30 CONTINUE
      if (N < 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          DX(I) = DA*DX(I)
          DX(I+1) = DA*DX(I+1)
          DX(I+2) = DA*DX(I+2)
          DX(I+3) = DA*DX(I+3)
          DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      end subroutine DSCAL


!---------------------------------------------------------------------
      subroutine DSWAP(N,DX,INCX,DY,INCY)
!     .. Scalar Arguments ..
      integer :: INCX,INCY,N
!     ..
!     .. Array Arguments ..
      real(8) :: DX(*),DY(*)
!     ..
!
!  Purpose
!  =======
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      real(8) :: DTEMP
      integer :: I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      if (N <= 0) RETURN
      if (INCX == 1 .and. INCY == 1) GO TO 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      IX = 1
      IY = 1
      if (INCX < 0) IX = (-N+1)*INCX + 1
      if (INCY < 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DTEMP = DX(IX)
          DX(IX) = DY(IY)
          DY(IY) = DTEMP
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
   20 M = MOD(N,3)
      if (M == 0) GO TO 40
      DO 30 I = 1,M
          DTEMP = DX(I)
          DX(I) = DY(I)
          DY(I) = DTEMP
   30 CONTINUE
      if (N < 3) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
          DTEMP = DX(I)
          DX(I) = DY(I)
          DY(I) = DTEMP
          DTEMP = DX(I+1)
          DX(I+1) = DY(I+1)
          DY(I+1) = DTEMP
          DTEMP = DX(I+2)
          DX(I+2) = DY(I+2)
          DY(I+2) = DTEMP
   50 CONTINUE
      RETURN
      end subroutine DSWAP


!---------------------------------------------------------------------
      subroutine DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
!     .. Scalar Arguments ..
      real(8) :: ALPHA
      integer :: INCX,LDA,N
      character UPLO
!     ..
!     .. Array Arguments ..
      real(8) :: A(LDA,*),X(*)
!     ..
!
!  Purpose
!  =======
!
!  DSYR   performs the symmetric rank 1 operation
!
!     A := alpha*x*x' + A,
!
!  where alpha is a real scalar, x is an n element vector and A is an
!  n by n symmetric matrix.
!
!  Arguments
!  ==========
!
!  UPLO   - character*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - integer ::.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real(8) ::.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - real(8) :: array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - integer ::.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  A      - real(8) :: array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced. On exit, the
!           upper triangular part of the array A is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced. On exit, the
!           lower triangular part of the array A is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDA    - integer ::.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      real(8) :: ZERO
      PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      real(8) :: TEMP
      integer :: I,INFO,IX,J,JX,KX
!     ..
!     .. External Functions ..
!m       logical :: LSAME
!m       EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!m       EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
!m       INTRINSIC max
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      if (.NOT.LSAME(UPLO,'U') .and. .NOT.LSAME(UPLO,'L')) then
          INFO = 1
      else if (N < 0) then
          INFO = 2
      else if (INCX == 0) then
          INFO = 5
      else if (LDA < max(1,N)) then
          INFO = 7
      end if
      if (INFO.NE.0) then
          call XERBLA('DSYR  ',INFO)
          RETURN
      end if
!
!     Quick return if possible.
!
      if ((N == 0) .or. (ALPHA == ZERO)) RETURN
!
!     Set the start point in X if the increment is not unity.
!
      if (INCX <= 0) then
          KX = 1 - (N-1)*INCX
      else if (INCX.NE.1) then
          KX = 1
      end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
      if (LSAME(UPLO,'U')) then
!
!        Form  A  when A is stored in upper triangle.
!
          if (INCX == 1) then
              DO 20 J = 1,N
                  if (X(J).NE.ZERO) then
                      TEMP = ALPHA*X(J)
                      DO 10 I = 1,J
                          A(I,J) = A(I,J) + X(I)*TEMP
   10                 CONTINUE
                  end if
   20         CONTINUE
          else
              JX = KX
              DO 40 J = 1,N
                  if (X(JX).NE.ZERO) then
                      TEMP = ALPHA*X(JX)
                      IX = KX
                      DO 30 I = 1,J
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
   30                 CONTINUE
                  end if
                  JX = JX + INCX
   40         CONTINUE
          end if
      else
!
!        Form  A  when A is stored in lower triangle.
!
          if (INCX == 1) then
              DO 60 J = 1,N
                  if (X(J).NE.ZERO) then
                      TEMP = ALPHA*X(J)
                      DO 50 I = J,N
                          A(I,J) = A(I,J) + X(I)*TEMP
   50                 CONTINUE
                  end if
   60         CONTINUE
          else
              JX = KX
              DO 80 J = 1,N
                  if (X(JX).NE.ZERO) then
                      TEMP = ALPHA*X(JX)
                      IX = JX
                      DO 70 I = J,N
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
   70                 CONTINUE
                  end if
                  JX = JX + INCX
   80         CONTINUE
          end if
      end if
!
      RETURN
!
!     End of DSYR  .
!
      end subroutine DSYR


!---------------------------------------------------------------------
      integer function IDAMAX(N,DX,INCX)
!     .. Scalar Arguments ..
      integer :: INCX,N
!     ..
!     .. Array Arguments ..
      real(8) :: DX(*)
!     ..
!
!  Purpose
!  =======
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      real(8) :: DMAX
      integer :: I,IX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DABS
!     ..
      IDAMAX = 0
      if (N < 1 .or. INCX <= 0) RETURN
      IDAMAX = 1
      if (N == 1) RETURN
      if (INCX == 1) GO TO 20
!
!        code for increment not equal to 1
!
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
          if (DABS(DX(IX)) <= DMAX) GO TO 5
          IDAMAX = I
          DMAX = DABS(DX(IX))
    5     IX = IX + INCX
   10 CONTINUE
      RETURN
!
!        code for increment equal to 1
!
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
          if (DABS(DX(I)) <= DMAX) GO TO 30
          IDAMAX = I
          DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
      end function IDAMAX


!---------------------------------------------------------------------
      integer function ISAMAX(N,SX,INCX)
!     .. Scalar Arguments ..
      integer :: INCX,N
!     ..
!     .. Array Arguments ..
      real :: SX(*)
!     ..
!
!  Purpose
!  =======
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      real :: SMAX
      integer :: I,IX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC abs
!     ..
      ISAMAX = 0
      if (N < 1 .or. INCX <= 0) RETURN
      ISAMAX = 1
      if (N == 1) RETURN
      if (INCX == 1) GO TO 20
!
!        code for increment not equal to 1
!
      IX = 1
      SMAX = abs(SX(1))
      IX = IX + INCX
      DO 10 I = 2,N
          if (abs(SX(IX)) <= SMAX) GO TO 5
          ISAMAX = I
          SMAX = abs(SX(IX))
    5     IX = IX + INCX
   10 CONTINUE
      RETURN
!
!        code for increment equal to 1
!
   20 SMAX = abs(SX(1))
      DO 30 I = 2,N
          if (abs(SX(I)) <= SMAX) GO TO 30
          ISAMAX = I
          SMAX = abs(SX(I))
   30 CONTINUE
      RETURN
      end function ISAMAX



!---------------------------------------------------------------------
      subroutine SCOPY(N,SX,INCX,SY,INCY)
!     .. Scalar Arguments ..
      integer :: INCX,INCY,N
!     ..
!     .. Array Arguments ..
      real :: SX(*),SY(*)
!     ..
!
!  Purpose
!  =======
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to 1.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      integer :: I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      if (N <= 0) RETURN
      if (INCX == 1 .and. INCY == 1) GO TO 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      IX = 1
      IY = 1
      if (INCX < 0) IX = (-N+1)*INCX + 1
      if (INCY < 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          SY(IY) = SX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,7)
      if (M == 0) GO TO 40
      DO 30 I = 1,M
          SY(I) = SX(I)
   30 CONTINUE
      if (N < 7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
          SY(I) = SX(I)
          SY(I+1) = SX(I+1)
          SY(I+2) = SX(I+2)
          SY(I+3) = SX(I+3)
          SY(I+4) = SX(I+4)
          SY(I+5) = SX(I+5)
          SY(I+6) = SX(I+6)
   50 CONTINUE
      RETURN
      end subroutine SCOPY


!---------------------------------------------------------------------
      subroutine SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!     .. Scalar Arguments ..
      real :: ALPHA,BETA
      integer :: K,LDA,LDB,LDC,M,N
      character TRANSA,TRANSB
!     ..
!     .. Array Arguments ..
      real :: A(LDA,*),B(LDB,*),C(LDC,*)
!     ..
!
!  Purpose
!  =======
!
!  SGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Arguments
!  ==========
!
!  TRANSA - character*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - character*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - integer ::.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - integer ::.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - integer ::.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - real ::            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - real ::             array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - integer ::.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - real ::             array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - integer ::.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - real ::            .
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - real ::             array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - integer ::.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
!m       logical :: LSAME
!m       EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!m       EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
!m       INTRINSIC max
!     ..
!     .. Local Scalars ..
      real :: TEMP
      integer :: I,INFO,J,L,NCOLA,NROWA,NROWB
      logical :: NOTA,NOTB
!     ..
!     .. Parameters ..
      real :: ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      if (NOTA) then
          NROWA = M
          NCOLA = K
      else
          NROWA = K
          NCOLA = M
      end if
      if (NOTB) then
          NROWB = K
      else
          NROWB = N
      end if
!
!     Test the input parameters.
!
      INFO = 0
      if ((.NOT.NOTA) .and. (.NOT.LSAME(TRANSA,'C')) .and.       &
          (.NOT.LSAME(TRANSA,'T'))) then
          INFO = 1
      else if ((.NOT.NOTB) .and. (.NOT.LSAME(TRANSB,'C')) .and.  &
               (.NOT.LSAME(TRANSB,'T'))) then
          INFO = 2
      else if (M < 0) then
          INFO = 3
      else if (N < 0) then
          INFO = 4
      else if (K < 0) then
          INFO = 5
      else if (LDA < max(1,NROWA)) then
          INFO = 8
      else if (LDB < max(1,NROWB)) then
          INFO = 10
      else if (LDC < max(1,M)) then
          INFO = 13
      end if
      if (INFO.NE.0) then
          call XERBLA('SGEMM ',INFO)
          RETURN
      end if
!
!     Quick return if possible.
!
      if ((M == 0) .or. (N == 0) .or.  &
          (((ALPHA == ZERO).or. (K == 0)).and. (BETA == ONE))) RETURN
!
!     And if  alpha.eq.zero.
!
      if (ALPHA == ZERO) then
          if (BETA == ZERO) then
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          else
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          end if
          RETURN
      end if
!
!     Start the operations.
!
      if (NOTB) then
          if (NOTA) then
!
!           Form  C := alpha*A*B + beta*C.
!
              DO 90 J = 1,N
                  if (BETA == ZERO) then
                      DO 50 I = 1,M
                          C(I,J) = ZERO
   50                 CONTINUE
                  else if (BETA.NE.ONE) then
                      DO 60 I = 1,M
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  end if
                  DO 80 L = 1,K
                      if (B(L,J).NE.ZERO) then
                          TEMP = ALPHA*B(L,J)
                          DO 70 I = 1,M
                              C(I,J) = C(I,J) + TEMP*A(I,L)
   70                     CONTINUE
                      end if
   80             CONTINUE
   90         CONTINUE
          else
!
!           Form  C := alpha*A'*B + beta*C
!
              DO 120 J = 1,N
                  DO 110 I = 1,M
                      TEMP = ZERO
                      DO 100 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  100                 CONTINUE
                      if (BETA == ZERO) then
                          C(I,J) = ALPHA*TEMP
                      else
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      end if
  110             CONTINUE
  120         CONTINUE
          end if
      else
          if (NOTA) then
!
!           Form  C := alpha*A*B' + beta*C
!
              DO 170 J = 1,N
                  if (BETA == ZERO) then
                      DO 130 I = 1,M
                          C(I,J) = ZERO
  130                 CONTINUE
                  else if (BETA.NE.ONE) then
                      DO 140 I = 1,M
                          C(I,J) = BETA*C(I,J)
  140                 CONTINUE
                  end if
                  DO 160 L = 1,K
                      if (B(J,L).NE.ZERO) then
                          TEMP = ALPHA*B(J,L)
                          DO 150 I = 1,M
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  150                     CONTINUE
                      end if
  160             CONTINUE
  170         CONTINUE
          else
!
!           Form  C := alpha*A'*B' + beta*C
!
              DO 200 J = 1,N
                  DO 190 I = 1,M
                      TEMP = ZERO
                      DO 180 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  180                 CONTINUE
                      if (BETA == ZERO) then
                          C(I,J) = ALPHA*TEMP
                      else
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      end if
  190             CONTINUE
  200         CONTINUE
          end if
      end if
!
      RETURN
!
!     End of SGEMM .
!
      end subroutine SGEMM


!---------------------------------------------------------------------
      subroutine SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!     .. Scalar Arguments ..
      real :: ALPHA,BETA
      integer :: INCX,INCY,LDA,M,N
      character TRANS
!     ..
!     .. Array Arguments ..
      real :: A(LDA,*),X(*),Y(*)
!     ..
!
!  Purpose
!  =======
!
!  SGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Arguments
!  ==========
!
!  TRANS  - character*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - integer ::.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - integer ::.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real ::            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - real ::             array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - integer ::.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - real ::             array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - integer ::.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - real ::            .
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - real ::             array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - integer ::.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      real :: ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      real :: TEMP
      integer :: I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
!     ..
!     .. External Functions ..
!m       logical :: LSAME
!m       EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!m       EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
!m       INTRINSIC max
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      if (.NOT.LSAME(TRANS,'N') .and. .NOT.LSAME(TRANS,'T') .and.   &
          .NOT.LSAME(TRANS,'C')) then
          INFO = 1
      else if (M < 0) then
          INFO = 2
      else if (N < 0) then
          INFO = 3
      else if (LDA < max(1,M)) then
          INFO = 6
      else if (INCX == 0) then
          INFO = 8
      else if (INCY == 0) then
          INFO = 11
      end if
      if (INFO.NE.0) then
          call XERBLA('SGEMV ',INFO)
          RETURN
      end if
!
!     Quick return if possible.
!
      if ((M == 0) .or. (N == 0) .or.  &
          ((ALPHA == ZERO).and. (BETA == ONE))) RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      if (LSAME(TRANS,'N')) then
          LENX = N
          LENY = M
      else
          LENX = M
          LENY = N
      end if
      if (INCX > 0) then
          KX = 1
      else
          KX = 1 - (LENX-1)*INCX
      end if
      if (INCY > 0) then
          KY = 1
      else
          KY = 1 - (LENY-1)*INCY
      end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      if (BETA.NE.ONE) then
          if (INCY == 1) then
              if (BETA == ZERO) then
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              else
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              end if
          else
              IY = KY
              if (BETA == ZERO) then
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              else
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              end if
          end if
      end if
      if (ALPHA == ZERO) RETURN
      if (LSAME(TRANS,'N')) then
!
!        Form  y := alpha*A*x + y.
!
          JX = KX
          if (INCY == 1) then
              DO 60 J = 1,N
                  if (X(JX).NE.ZERO) then
                      TEMP = ALPHA*X(JX)
                      DO 50 I = 1,M
                          Y(I) = Y(I) + TEMP*A(I,J)
   50                 CONTINUE
                  end if
                  JX = JX + INCX
   60         CONTINUE
          else
              DO 80 J = 1,N
                  if (X(JX).NE.ZERO) then
                      TEMP = ALPHA*X(JX)
                      IY = KY
                      DO 70 I = 1,M
                          Y(IY) = Y(IY) + TEMP*A(I,J)
                          IY = IY + INCY
   70                 CONTINUE
                  end if
                  JX = JX + INCX
   80         CONTINUE
          end if
      else
!
!        Form  y := alpha*A'*x + y.
!
          JY = KY
          if (INCX == 1) then
              DO 100 J = 1,N
                  TEMP = ZERO
                  DO 90 I = 1,M
                      TEMP = TEMP + A(I,J)*X(I)
   90             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  100         CONTINUE
          else
              DO 120 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  DO 110 I = 1,M
                      TEMP = TEMP + A(I,J)*X(IX)
                      IX = IX + INCX
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  120         CONTINUE
          end if
      end if
!
      RETURN
!
!     End of SGEMV .
!
      end  subroutine SGEMV


!---------------------------------------------------------------------
      subroutine SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!     .. Scalar Arguments ..
      real :: ALPHA
      integer :: INCX,INCY,LDA,M,N
!     ..
!     .. Array Arguments ..
      real :: A(LDA,*),X(*),Y(*)
!     ..
!
!  Purpose
!  =======
!
!  SGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Arguments
!  ==========
!
!  M      - integer ::.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - integer ::.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real ::            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - real ::             array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - integer ::.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - real ::             array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - integer ::.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - real ::             array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - integer ::.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      real :: ZERO
      PARAMETER (ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      real :: TEMP
      integer :: I,INFO,IX,J,JY,KX
!     ..
!     .. External Subroutines ..
!m       EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
!m       INTRINSIC max
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      if (M < 0) then
          INFO = 1
      else if (N < 0) then
          INFO = 2
      else if (INCX == 0) then
          INFO = 5
      else if (INCY == 0) then
          INFO = 7
      else if (LDA < max(1,M)) then
          INFO = 9
      end if
      if (INFO.NE.0) then
          call XERBLA('SGER  ',INFO)
          RETURN
      end if
!
!     Quick return if possible.
!
      if ((M == 0) .or. (N == 0) .or. (ALPHA == ZERO)) RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      if (INCY > 0) then
          JY = 1
      else
          JY = 1 - (N-1)*INCY
      end if
      if (INCX == 1) then
          DO 20 J = 1,N
              if (Y(JY).NE.ZERO) then
                  TEMP = ALPHA*Y(JY)
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              end if
              JY = JY + INCY
   20     CONTINUE
      else
          if (INCX > 0) then
              KX = 1
          else
              KX = 1 - (M-1)*INCX
          end if
          DO 40 J = 1,N
              if (Y(JY).NE.ZERO) then
                  TEMP = ALPHA*Y(JY)
                  IX = KX
                  DO 30 I = 1,M
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              end if
              JY = JY + INCY
   40     CONTINUE
      end if
!
      RETURN
!
!     End of SGER  .
!
      end subroutine SGER


!---------------------------------------------------------------------
      subroutine SSCAL(N,SA,SX,INCX)
!     .. Scalar Arguments ..
      real :: SA
      integer :: INCX,N
!     ..
!     .. Array Arguments ..
      real :: SX(*)
!     ..
!
!  Purpose
!  =======
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to 1.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      integer :: I,M,MP1,NINCX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      if (N <= 0 .or. INCX <= 0) RETURN
      if (INCX == 1) GO TO 20
!
!        code for increment not equal to 1
!
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          SX(I) = SA*SX(I)
   10 CONTINUE
      RETURN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,5)
      if (M == 0) GO TO 40
      DO 30 I = 1,M
          SX(I) = SA*SX(I)
   30 CONTINUE
      if (N < 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          SX(I) = SA*SX(I)
          SX(I+1) = SA*SX(I+1)
          SX(I+2) = SA*SX(I+2)
          SX(I+3) = SA*SX(I+3)
          SX(I+4) = SA*SX(I+4)
   50 CONTINUE
      RETURN
      end subroutine SSCAL


!---------------------------------------------------------------------
      subroutine SSWAP(N,SX,INCX,SY,INCY)
!     .. Scalar Arguments ..
      integer :: INCX,INCY,N
!     ..
!     .. Array Arguments ..
      real :: SX(*),SY(*)
!     ..
!
!  Purpose
!  =======
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal to 1.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      real :: STEMP
      integer :: I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      if (N <= 0) RETURN
      if (INCX == 1 .and. INCY == 1) GO TO 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      IX = 1
      IY = 1
      if (INCX < 0) IX = (-N+1)*INCX + 1
      if (INCY < 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          STEMP = SX(IX)
          SX(IX) = SY(IY)
          SY(IY) = STEMP
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
   20 M = MOD(N,3)
      if (M == 0) GO TO 40
      DO 30 I = 1,M
          STEMP = SX(I)
          SX(I) = SY(I)
          SY(I) = STEMP
   30 CONTINUE
      if (N < 3) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
          STEMP = SX(I)
          SX(I) = SY(I)
          SY(I) = STEMP
          STEMP = SX(I+1)
          SX(I+1) = SY(I+1)
          SY(I+1) = STEMP
          STEMP = SX(I+2)
          SX(I+2) = SY(I+2)
          SY(I+2) = STEMP
   50 CONTINUE
      RETURN
      end subroutine SSWAP


!---------------------------------------------------------------------
      subroutine SSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
!     .. Scalar Arguments ..
      real :: ALPHA
      integer :: INCX,LDA,N
      character UPLO
!     ..
!     .. Array Arguments ..
      real :: A(LDA,*),X(*)
!     ..
!
!  Purpose
!  =======
!
!  SSYR   performs the symmetric rank 1 operation
!
!     A := alpha*x*x' + A,
!
!  where alpha is a real scalar, x is an n element vector and A is an
!  n by n symmetric matrix.
!
!  Arguments
!  ==========
!
!  UPLO   - character*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - integer ::.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real ::            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - real ::             array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - integer ::.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  A      - real ::             array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced. On exit, the
!           upper triangular part of the array A is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced. On exit, the
!           lower triangular part of the array A is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDA    - integer ::.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      real :: ZERO
      PARAMETER (ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      real :: TEMP
      integer :: I,INFO,IX,J,JX,KX
!     ..
!     .. External Functions ..
!m       logical :: LSAME
!m       EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!m       EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
!m       INTRINSIC max
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      if (.NOT.LSAME(UPLO,'U') .and. .NOT.LSAME(UPLO,'L')) then
          INFO = 1
      else if (N < 0) then
          INFO = 2
      else if (INCX == 0) then
          INFO = 5
      else if (LDA < max(1,N)) then
          INFO = 7
      end if
      if (INFO.NE.0) then
          call XERBLA('SSYR  ',INFO)
          RETURN
      end if
!
!     Quick return if possible.
!
      if ((N == 0) .or. (ALPHA == ZERO)) RETURN
!
!     Set the start point in X if the increment is not unity.
!
      if (INCX <= 0) then
          KX = 1 - (N-1)*INCX
      else if (INCX.NE.1) then
          KX = 1
      end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
      if (LSAME(UPLO,'U')) then
!
!        Form  A  when A is stored in upper triangle.
!
          if (INCX == 1) then
              DO 20 J = 1,N
                  if (X(J).NE.ZERO) then
                      TEMP = ALPHA*X(J)
                      DO 10 I = 1,J
                          A(I,J) = A(I,J) + X(I)*TEMP
   10                 CONTINUE
                  end if
   20         CONTINUE
          else
              JX = KX
              DO 40 J = 1,N
                  if (X(JX).NE.ZERO) then
                      TEMP = ALPHA*X(JX)
                      IX = KX
                      DO 30 I = 1,J
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
   30                 CONTINUE
                  end if
                  JX = JX + INCX
   40         CONTINUE
          end if
      else
!
!        Form  A  when A is stored in lower triangle.
!
          if (INCX == 1) then
              DO 60 J = 1,N
                  if (X(J).NE.ZERO) then
                      TEMP = ALPHA*X(J)
                      DO 50 I = J,N
                          A(I,J) = A(I,J) + X(I)*TEMP
   50                 CONTINUE
                  end if
   60         CONTINUE
          else
              JX = KX
              DO 80 J = 1,N
                  if (X(JX).NE.ZERO) then
                      TEMP = ALPHA*X(JX)
                      IX = JX
                      DO 70 I = J,N
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
   70                 CONTINUE
                  end if
                  JX = JX + INCX
   80         CONTINUE
          end if
      end if
!
      RETURN
!
!     End of SSYR  .
!
      end subroutine SSYR


!---------------------------------------------------------------------
      subroutine XERBLA(SRNAME,INFO)
!
!  -- LAPACK auxiliary routine (preliminary version) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      integer :: INFO
      character*6 SRNAME
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) character*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) integer ::
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
!
      WRITE (*,FMT=9999) SRNAME,INFO
!
      STOP
!
 9999 FORMAT (' ** On entry to ',A6,' parameter number ',I2,' had ',  &
             'an illegal value')
!
!     End of XERBLA
!
      end subroutine XERBLA





!----------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------
!-------------------------------   NUEVO: CALCULO DE CONDITION NUMBER -------------------------------------------
!----------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------

! =========== DOCUMENTATION ===========
!
! Online html documentation available at
! http://www.netlib.org/lapack/explore-html/
!
! \htmlonly
! Download DSYCON + dependencies
!
! Definition:
! ===========
!
! SUBROUTINE DSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK,
! IWORK, INFO )
!
! .. Scalar Arguments ..
! CHARACTER UPLO
! INTEGER INFO, LDA, N
! DOUBLE PRECISION ANORM, RCOND
! ..
! .. Array Arguments ..
! INTEGER IPIV( * ), IWORK( * )
! DOUBLE PRECISION A( LDA, * ), WORK( * )
! ..
!
!
!> \par Purpose:
! =============
!>
!> \verbatim
!>
!> DSYCON estimates the reciprocal of the condition number (in the
!> 1-norm) of a real symmetric matrix A using the factorization
!> A = U*D*U**T or A = L*D*L**T computed by DSYTRF.
!>
!> An estimate is obtained for norm(inv(A)), and the reciprocal of the
!> condition number is computed as RCOND = 1 / (ANORM! norm(inv(A))).
!> \endverbatim
!
! Arguments:
! ==========
!
!> \param[in] UPLO
!> \verbatim
!> UPLO is CHARACTER*1
!> Specifies whether the details of the factorization are stored
!> as an upper or lower triangular matrix.
!> = 'U': Upper triangular, form is A = U*D*U**T;
!> = 'L': Lower triangular, form is A = L*D*L**T.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!> N is INTEGER
!> The order of the matrix A. N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!> A is DOUBLE PRECISION array, dimension (LDA,N)
!> The block diagonal matrix D and the multipliers used to
!> obtain the factor U or L as computed by DSYTRF.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!> LDA is INTEGER
!> The leading dimension of the array A. LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!> IPIV is INTEGER array, dimension (N)
!> Details of the interchanges and the block structure of D
!> as determined by DSYTRF.
!> \endverbatim
!>
!> \param[in] ANORM
!> \verbatim
!> ANORM is DOUBLE PRECISION
!> The 1-norm of the original matrix A.
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!> RCOND is DOUBLE PRECISION
!> The reciprocal of the condition number of the matrix A,
!> computed as RCOND = 1/(ANORM! AINVNM), where AINVNM is an
!> estimate of the 1-norm of inv(A) computed in this routine.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!> WORK is DOUBLE PRECISION array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!> IWORK is INTEGER array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!> INFO is INTEGER
!> = 0: successful exit
!> < 0: if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
! Authors:
! ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2011
!
!> \ingroup doubleSYcomputational
!
! =====================================================================
!  =====================================================================
      SUBROUTINE DSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK, &
                        IWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   ANORM, RCOND
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * ), IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, KASE
      DOUBLE PRECISION   AINVNM
!     ..
!     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
!     ..
!     .. External Functions ..
!m      LOGICAL            LSAME
!m      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!m      EXTERNAL           DLACN2, DSYTRS, XERBLA
!     ..
!     .. Intrinsic Functions ..
!m      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYCON', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.LE.ZERO ) THEN
         RETURN
      END IF

!
!     Check that the diagonal matrix D is nonsingular.
!
      IF( UPPER ) THEN
!
!        Upper triangular storage: examine D from bottom to top
!
         DO 10 I = N, 1, -1
            IF( IPIV( I ).GT.0 .AND. A( I, I ).EQ.ZERO ) RETURN
   10    CONTINUE
      ELSE
!
!        Lower triangular storage: examine D from top to bottom.
!
         DO 20 I = 1, N
            IF( IPIV( I ).GT.0 .AND. A( I, I ).EQ.ZERO ) RETURN
   20    CONTINUE
      END IF
!
!     Estimate the 1-norm of the inverse.
!
      KASE = 0
   30 CONTINUE
      CALL DLACN2( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
!
!        Multiply by inv(L*D*L**T) or inv(U*D*U**T).
!
         CALL DSYTRS( UPLO, N, 1, A, LDA, IPIV, WORK, N, INFO )
         GO TO 30
      END IF
!
!     Compute the estimate of the reciprocal condition number.
!
      IF( AINVNM.NE.ZERO ) &
        RCOND = ( ONE / AINVNM ) / ANORM
!
      RETURN
!
!     End of DSYCON
!
      END subroutine DSYCON


!  Definition:
!  ===========
!
!       SUBROUTINE DLACN2( N, V, X, ISGN, EST, KASE, ISAVE )
!
!       .. Scalar Arguments ..
!       INTEGER            KASE, N
!       DOUBLE PRECISION   EST
!       ..
!       .. Array Arguments ..
!       INTEGER            ISGN( * ), ISAVE( 3 )
!       DOUBLE PRECISION   V( ! ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLACN2 estimates the 1-norm of a square, real matrix A.
!> Reverse communication is used for evaluating matrix-vector products.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The order of the matrix.  N >= 1.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension (N)
!>         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
!>         (W is not returned).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (N)
!>         On an intermediate return, X should be overwritten by
!>               A * X,   if KASE=1,
!>               A**T * X,  if KASE=2,
!>         and DLACN2 must be re-called with all the other parameters
!>         unchanged.
!> \endverbatim
!>
!> \param[out] ISGN
!> \verbatim
!>          ISGN is INTEGER array, dimension (N)
!> \endverbatim
!>
!> \param[in,out] EST
!> \verbatim
!>          EST is DOUBLE PRECISION
!>         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be
!>         unchanged from the previous call to DLACN2.
!>         On exit, EST is an estimate (a lower bound) for norm(A).
!> \endverbatim
!>
!> \param[in,out] KASE
!> \verbatim
!>          KASE is INTEGER
!>         On the initial call to DLACN2, KASE should be 0.
!>         On an intermediate return, KASE will be 1 or 2, indicating
!>         whether X should be overwritten by A * X  or A**T * X.
!>         On the final return from DLACN2, KASE will again be 0.
!> \endverbatim
!>
!> \param[in,out] ISAVE
!> \verbatim
!>          ISAVE is INTEGER array, dimension (3)
!>         ISAVE is used to save variables between calls to DLACN2
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date September 2012
!
!> \ingroup doubleOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Originally named SONEST, dated March 16, 1988.
!>
!>  This is a thread safe version of DLACON, which uses the array ISAVE
!>  in place of a SAVE statement, as follows:
!>
!>     DLACON     DLACN2
!>      JUMP     ISAVE(1)
!>      J        ISAVE(2)
!>      ITER     ISAVE(3)
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!>     Nick Higham, University of Manchester
!
!> \par References:
!  ================
!>
!>  N.J. Higham, "FORTRAN codes for estimating the one-norm of
!>  a real or complex matrix, with applications to condition estimation",
!>  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
!>
!  =====================================================================
      SUBROUTINE DLACN2( N, V, X, ISGN, EST, KASE, ISAVE )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            KASE, N
      DOUBLE PRECISION   EST
!     ..
!     .. Array Arguments ..
      INTEGER            ISGN( * ), ISAVE( 3 )
      DOUBLE PRECISION   V( * ), X( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, JLAST
      DOUBLE PRECISION   ALTSGN, ESTOLD, TEMP
!     ..
!     .. External Functions ..
!m      INTEGER            IDAMAX
!m      DOUBLE PRECISION   DASUM
!m      EXTERNAL           IDAMAX, DASUM
!     ..
!     .. External Subroutines ..
!m      EXTERNAL           DCOPY
!     ..
!     .. Intrinsic Functions ..
!m      INTRINSIC          ABS, DBLE, NINT, SIGN
!     ..
!     .. Executable Statements ..
!
      IF( KASE.EQ.0 ) THEN
         DO 10 I = 1, N
            X( I ) = ONE / DBLE( N )
   10    CONTINUE
         KASE = 1
         ISAVE( 1 ) = 1
         RETURN
      END IF
!
      GO TO ( 20, 40, 70, 110, 140 )ISAVE( 1 )
!
!     ................ ENTRY   (ISAVE( 1 ) = 1)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A!X.
!
   20 CONTINUE
      IF( N.EQ.1 ) THEN
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
!        ... QUIT
         GO TO 150
      END IF
      EST = DASUM( N, X, 1 )
!
      DO 30 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
   30 CONTINUE
      KASE = 2
      ISAVE( 1 ) = 2
      RETURN
!
!     ................ ENTRY   (ISAVE( 1 ) = 2)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
!
   40 CONTINUE
      ISAVE( 2 ) = IDAMAX( N, X, 1 )
      ISAVE( 3 ) = 2
!
!     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
!
   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = ZERO
   60 CONTINUE
      X( ISAVE( 2 ) ) = ONE
      KASE = 1
      ISAVE( 1 ) = 3
      RETURN
!
!     ................ ENTRY   (ISAVE( 1 ) = 3)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
   70 CONTINUE
      CALL DCOPY( N, X, 1, V, 1 )
      ESTOLD = EST
      EST = DASUM( N, V, 1 )
      DO 80 I = 1, N
         IF( NINT( SIGN( ONE, X( I ) ) ).NE.ISGN( I ) ) GO TO 90
   80 CONTINUE
!     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 120
!
   90 CONTINUE
!     TEST FOR CYCLING.
      IF( EST.LE.ESTOLD )  GO TO 120
!
      DO 100 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
  100 CONTINUE
      KASE = 2
      ISAVE( 1 ) = 4
      RETURN
!
!     ................ ENTRY   (ISAVE( 1 ) = 4)
!     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
!
  110 CONTINUE
      JLAST = ISAVE( 2 )
      ISAVE( 2 ) = IDAMAX( N, X, 1 )
      IF( ( X( JLAST ).NE.ABS( X( ISAVE( 2 ) ) ) ) .AND.   &
         ( ISAVE( 3 ).LT.ITMAX ) ) THEN
         ISAVE( 3 ) = ISAVE( 3 ) + 1
         GO TO 50
      END IF
!
!     ITERATION COMPLETE.  FINAL STAGE.
!
  120 CONTINUE
      ALTSGN = ONE
      DO 130 I = 1, N
         X( I ) = ALTSGN*( ONE+DBLE( I-1 ) / DBLE( N-1 ) )
         ALTSGN = -ALTSGN
  130 CONTINUE
      KASE = 1
      ISAVE( 1 ) = 5
      RETURN
!
!     ................ ENTRY   (ISAVE( 1 ) = 5)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
  140 CONTINUE
      TEMP = TWO*( DASUM( N, X, 1 ) / DBLE( 3*N ) )
      IF( TEMP.GT.EST ) THEN
         CALL DCOPY( N, X, 1, V, 1 )
         EST = TEMP
      END IF
!
  150 CONTINUE
      KASE = 0
      RETURN
!
!     End of DLACN2
!
      END subroutine DLACN2

!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DASUM takes the sum of the absolute values.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2011
!
!> \ingroup double_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
!
!  -- Reference BLAS level1 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,M,MP1,NINCX
!     ..
!     .. Intrinsic Functions ..
!m      INTRINSIC DABS,MOD
!     ..
      DASUM = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
!        code for increment equal to 1
!
!
!        clean-up loop
!
         M = MOD(N,6)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DTEMP + DABS(DX(I))
            END DO
            IF (N.LT.6) THEN
               DASUM = DTEMP
               RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,6
            DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I+1)) +  &
                    DABS(DX(I+2)) + DABS(DX(I+3)) +        &
                    DABS(DX(I+4)) + DABS(DX(I+5))
         END DO
      ELSE
!
!        code for increment not equal to 1
!
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            DTEMP = DTEMP + DABS(DX(I))
         END DO
      END IF
      DASUM = DTEMP
      RETURN
      END function DASUM

!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*),DY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DCOPY copies a vector, x, to a vector, y.
!>    uses unrolled loops for increments equal to one.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2011
!
!> \ingroup double_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
!
!  -- Reference BLAS level1 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
!m      INTRINSIC MOD
!     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
         M = MOD(N,7)
         IF (M.NE.0) THEN
            DO I = 1,M
               DY(I) = DX(I)
            END DO
            IF (N.LT.7) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,7
            DY(I) = DX(I)
            DY(I+1) = DX(I+1)
            DY(I+2) = DX(I+2)
            DY(I+3) = DX(I+3)
            DY(I+4) = DX(I+4)
            DY(I+5) = DX(I+5)
            DY(I+6) = DX(I+6)
         END DO
      ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DY(IY) = DX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END subroutine DCOPY

!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            LDA, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLANGE  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> real matrix A.
!> \endverbatim
!>
!> \return DLANGE
!> \verbatim
!>
!>    DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!>             (
!>             ( norm1(A),         NORM = '1', 'O' or 'o'
!>             (
!>             ( normI(A),         NORM = 'I' or 'i'
!>             (
!>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!>
!> where  norm1  denotes the  one norm of a matrix (maximum column sum),
!> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!> normF  denotes the  Frobenius norm of a matrix (square root of sum of
!> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NORM
!> \verbatim
!>          NORM is CHARACTER*1
!>          Specifies the value to be returned in DLANGE as described
!>          above.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.  When M = 0,
!>          DLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.  When N = 0,
!>          DLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(M,1).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
!>          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
!>          referenced.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date September 2012
!
!> \ingroup doubleGEauxiliary
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   SCALE, SUM, VALUE, TEMP
!     ..
!     .. External Subroutines ..
!m      EXTERNAL           DLASSQ
!     ..
!     .. External Functions ..
!m     LOGICAL            LSAME, DISNAN
!m     EXTERNAL           LSAME, DISNAN
!     ..
!     .. Intrinsic Functions ..
!m      INTRINSIC          ABS, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( MIN( M, N ).EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, M
               TEMP = ABS( A( I, J ) )
               IF( VALUE.LT.TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
!
!        Find norm1(A).
!
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, M
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            IF( VALUE.LT.SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
         DO 50 I = 1, M
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, M
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, M
            TEMP = WORK( I )
            IF( VALUE.LT.TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL DLASSQ( M, A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE!SQRT( SUM )
      END IF
!
      DLANGE = VALUE
      RETURN
!
!     End of DLANGE
!
      END function DLANGE

!  =====================================================================
      LOGICAL FUNCTION DISNAN( DIN )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   DIN
!     ..
!
!  =====================================================================
!
!  .. External Functions ..
!m      LOGICAL DLAISNAN
!m      EXTERNAL DLAISNAN
!  ..
!  .. Executable Statements ..
      DISNAN = DLAISNAN(DIN,DIN)
      RETURN
      END function DISNAN


!  =====================================================================
      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
!     ..
!     .. External Functions ..
!m      LOGICAL            DISNAN
!m      EXTERNAL           DISNAN
!     ..
!     .. Intrinsic Functions ..
!m      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            ABSXI = ABS( X( IX ) )
            IF( ABSXI.GT.ZERO.OR.DISNAN( ABSXI ) ) THEN
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
!
!     End of DLASSQ
!
      END subroutine DLASSQ

end module lapack_solver
