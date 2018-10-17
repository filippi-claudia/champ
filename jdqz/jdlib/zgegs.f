*> \brief <b> ZGEEVX computes the eigenvalues and, optionally, the left
*and/or right eigenvectors for GE matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZGEGS + dependencies
*> <a
*href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgegs.f">
*> [TGZ]</a>
*> <a
*href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgegs.f">
*> [ZIP]</a>
*> <a
*href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgegs.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHA,
*       BETA,
*                         VSL, LDVSL, VSR, LDVSR, WORK, LWORK, RWORK,
*                         INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBVSL, JOBVSR
*       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
*      $                   BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ),
*      $                   WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> This routine is deprecated and has been replaced by routine ZGGES.
*>
*> ZGEGS computes the eigenvalues, Schur form, and, optionally, the
*> left and or/right Schur vectors of a complex matrix pair (A,B).
*> Given two square matrices A and B, the generalized Schur
*> factorization has the form
*>
*>    A = Q*S*Z**H,  B = Q*T*Z**H
*>
*> where Q and Z are unitary matrices and S and T are upper triangular.
*> The columns of Q are the left Schur vectors
*> and the columns of Z are the right Schur vectors.
*>
*> If only the eigenvalues of (A,B) are needed, the driver routine
*> ZGEGV should be used instead.  See ZGEGV for a description of the
*> eigenvalues of the generalized nonsymmetric eigenvalue problem
*> (GNEP).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBVSL
*> \verbatim
*>          JOBVSL is CHARACTER*1
*>          = 'N':  do not compute the left Schur vectors;
*>          = 'V':  compute the left Schur vectors (returned in VSL).
*> \endverbatim
*>
*> \param[in] JOBVSR
*> \verbatim
*>          JOBVSR is CHARACTER*1
*>          = 'N':  do not compute the right Schur vectors;
*>          = 'V':  compute the right Schur vectors (returned in VSR).
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A, B, VSL, and VSR.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA, N)
*>          On entry, the matrix A.
*>          On exit, the upper triangular matrix S from the generalized
*>          Schur factorization.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (LDB, N)
*>          On entry, the matrix B.
*>          On exit, the upper triangular matrix T from the generalized
*>          Schur factorization.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16 array, dimension (N)
*>          The complex scalars alpha that define the eigenvalues of
*>          GNEP.  ALPHA(j) = S(j,j), the diagonal element of the Schur
*>          form of A.
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is COMPLEX*16 array, dimension (N)
*>          The non-negative real scalars beta that define the
*>          eigenvalues of GNEP.  BETA(j) = T(j,j), the diagonal
*element
*>          of the triangular factor T.
*>
*>          Together, the quantities alpha = ALPHA(j) and beta =
*BETA(j)
*>          represent the j-th eigenvalue of the matrix pair (A,B), in
*>          one of the forms lambda = alpha/beta or mu = beta/alpha.
*>          Since either lambda or mu may overflow, they should not,
*>          in general, be computed.
*> \endverbatim
*>
*> \param[out] VSL
*> \verbatim
*>          VSL is COMPLEX*16 array, dimension (LDVSL,N)
*>          If JOBVSL = 'V', the matrix of left Schur vectors Q.
*>          Not referenced if JOBVSL = 'N'.
*> \endverbatim
*>
*> \param[in] LDVSL
*> \verbatim
*>          LDVSL is INTEGER
*>          The leading dimension of the matrix VSL. LDVSL >= 1, and
*>          if JOBVSL = 'V', LDVSL >= N.
*> \endverbatim
*>
*> \param[out] VSR
*> \verbatim
*>          VSR is COMPLEX*16 array, dimension (LDVSR,N)
*>          If JOBVSR = 'V', the matrix of right Schur vectors Z.
*>          Not referenced if JOBVSR = 'N'.
*> \endverbatim
*>
*> \param[in] LDVSR
*> \verbatim
*>          LDVSR is INTEGER
*>          The leading dimension of the matrix VSR. LDVSR >= 1, and
*>          if JOBVSR = 'V', LDVSR >= N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,2*N).
*>          For good performance, LWORK must generally be larger.
*>          To compute the optimal value of LWORK, call ILAENV to get
*>          blocksizes (for ZGEQRF, ZUNMQR, and CUNGQR.)  Then compute:
*>          NB  -- MAX of the blocksizes for ZGEQRF, ZUNMQR, and
*CUNGQR;
*>          the optimal LWORK is N*(NB+1).
*>
*>          If LWORK = -1, then a workspace query is assumed; the
*routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no
*error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (3*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          =1,...,N:
*>                The QZ iteration failed.  (A,B) are not in Schur
*>                form, but ALPHA(j) and BETA(j) should be correct for
*>                j=INFO+1,...,N.
*>          > N:  errors that usually indicate LAPACK problems:
*>                =N+1: error return from ZGGBAL
*>                =N+2: error return from ZGEQRF
*>                =N+3: error return from ZUNMQR
*>                =N+4: error return from ZUNGQR
*>                =N+5: error return from ZGGHRD
*>                =N+6: error return from ZHGEQZ (other than failed
*>                                               iteration)
*>                =N+7: error return from ZGGBAK (computing VSL)
*>                =N+8: error return from ZGGBAK (computing VSR)
*>                =N+9: error return from ZLASCL (various places)
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16GEeigen
*
*  =====================================================================
      SUBROUTINE zgegs( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHA, BETA,
     $                  VSL, LDVSL, VSR, LDVSR, WORK, LWORK, RWORK,
     $                  INFO )
*
*  -- LAPACK driver routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,
*  --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG
*  Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVSL, JOBVSR
      INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( lda, * ), ALPHA( * ), B( ldb, * ),
     $                   beta( * ), vsl( ldvsl, * ), vsr( ldvsr, * ),
     $                   work( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
      COMPLEX*16         CZERO, CONE
      parameter( czero = ( 0.0d0, 0.0d0 ),
     $                   cone = ( 1.0d0, 0.0d0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILASCL, ILBSCL, ILVSL, ILVSR, LQUERY
      INTEGER            ICOLS, IHI, IINFO, IJOBVL, IJOBVR, ILEFT, ILO,
     $                   iright, irows, irwork, itau, iwork, lopt,
     $                   lwkmin, lwkopt, nb, nb1, nb2, nb3
      DOUBLE PRECISION   ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS,
     $                   safmin, smlnum
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zgeqrf, zggbak, zggbal, zgghrd, zhgeqz,
     $                   zlacpy, zlascl, zlaset, zungqr, zunmqr
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           lsame, ilaenv, dlamch, zlange
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          int, max
*     ..
*     .. Executable Statements ..
*
*     Decode the input arguments
*
      IF( lsame( jobvsl, 'N' ) ) THEN
         ijobvl = 1
         ilvsl = .false.
      ELSE IF( lsame( jobvsl, 'V' ) ) THEN
         ijobvl = 2
         ilvsl = .true.
      ELSE
         ijobvl = -1
         ilvsl = .false.
      END IF
*
      IF( lsame( jobvsr, 'N' ) ) THEN
         ijobvr = 1
         ilvsr = .false.
      ELSE IF( lsame( jobvsr, 'V' ) ) THEN
         ijobvr = 2
         ilvsr = .true.
      ELSE
         ijobvr = -1
         ilvsr = .false.
      END IF
*
*     Test the input arguments
*
      lwkmin = max( 2*n, 1 )
      lwkopt = lwkmin
      work( 1 ) = lwkopt
      lquery = ( lwork.EQ.-1 )
      info = 0
      IF( ijobvl.LE.0 ) THEN
         info = -1
      ELSE IF( ijobvr.LE.0 ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( ldb.LT.max( 1, n ) ) THEN
         info = -7
      ELSE IF( ldvsl.LT.1 .OR. ( ilvsl .AND. ldvsl.LT.n ) ) THEN
         info = -11
      ELSE IF( ldvsr.LT.1 .OR. ( ilvsr .AND. ldvsr.LT.n ) ) THEN
         info = -13
      ELSE IF( lwork.LT.lwkmin .AND. .NOT.lquery ) THEN
         info = -15
      END IF
*
      IF( info.EQ.0 ) THEN
         nb1 = ilaenv( 1, 'ZGEQRF', ' ', n, n, -1, -1 )
         nb2 = ilaenv( 1, 'ZUNMQR', ' ', n, n, n, -1 )
         nb3 = ilaenv( 1, 'ZUNGQR', ' ', n, n, n, -1 )
         nb = max( nb1, nb2, nb3 )
         lopt = n*( nb+1 )
         work( 1 ) = lopt
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZGEGS ', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      eps = dlamch( 'E' )*dlamch( 'B' )
      safmin = dlamch( 'S' )
      smlnum = n*safmin / eps
      bignum = one / smlnum
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      anrm = zlange( 'M', n, n, a, lda, rwork )
      ilascl = .false.
      IF( anrm.GT.zero .AND. anrm.LT.smlnum ) THEN
         anrmto = smlnum
         ilascl = .true.
      ELSE IF( anrm.GT.bignum ) THEN
         anrmto = bignum
         ilascl = .true.
      END IF
*
      IF( ilascl ) THEN
         CALL zlascl( 'G', -1, -1, anrm, anrmto, n, n, a, lda, iinfo )
         IF( iinfo.NE.0 ) THEN
            info = n + 9
            RETURN
         END IF
      END IF
*
*     Scale B if max element outside range [SMLNUM,BIGNUM]
*
      bnrm = zlange( 'M', n, n, b, ldb, rwork )
      ilbscl = .false.
      IF( bnrm.GT.zero .AND. bnrm.LT.smlnum ) THEN
         bnrmto = smlnum
         ilbscl = .true.
      ELSE IF( bnrm.GT.bignum ) THEN
         bnrmto = bignum
         ilbscl = .true.
      END IF
*
      IF( ilbscl ) THEN
         CALL zlascl( 'G', -1, -1, bnrm, bnrmto, n, n, b, ldb, iinfo )
         IF( iinfo.NE.0 ) THEN
            info = n + 9
            RETURN
         END IF
      END IF
*
*     Permute the matrix to make it more nearly triangular
*
      ileft = 1
      iright = n + 1
      irwork = iright + n
      iwork = 1
      CALL zggbal( 'P', n, a, lda, b, ldb, ilo, ihi, rwork( ileft ),
     $             rwork( iright ), rwork( irwork ), iinfo )
      IF( iinfo.NE.0 ) THEN
         info = n + 1
         GO TO 10
      END IF
*
*     Reduce B to triangular form, and initialize VSL and/or VSR
*
      irows = ihi + 1 - ilo
      icols = n + 1 - ilo
      itau = iwork
      iwork = itau + irows
      CALL zgeqrf( irows, icols, b( ilo, ilo ), ldb, work( itau ),
     $             work( iwork ), lwork+1-iwork, iinfo )
      IF( iinfo.GE.0 )
     $   lwkopt = max( lwkopt, int( work( iwork ) )+iwork-1 )
      IF( iinfo.NE.0 ) THEN
         info = n + 2
         GO TO 10
      END IF
*
      CALL zunmqr( 'L', 'C', irows, icols, irows, b( ilo, ilo ), ldb,
     $             work( itau ), a( ilo, ilo ), lda, work( iwork ),
     $             lwork+1-iwork, iinfo )
      IF( iinfo.GE.0 )
     $   lwkopt = max( lwkopt, int( work( iwork ) )+iwork-1 )
      IF( iinfo.NE.0 ) THEN
         info = n + 3
         GO TO 10
      END IF
*
      IF( ilvsl ) THEN
         CALL zlaset( 'Full', n, n, czero, cone, vsl, ldvsl )
         CALL zlacpy( 'L', irows-1, irows-1, b( ilo+1, ilo ), ldb,
     $                vsl( ilo+1, ilo ), ldvsl )
         CALL zungqr( irows, irows, irows, vsl( ilo, ilo ), ldvsl,
     $                work( itau ), work( iwork ), lwork+1-iwork,
     $                iinfo )
         IF( iinfo.GE.0 )
     $      lwkopt = max( lwkopt, int( work( iwork ) )+iwork-1 )
         IF( iinfo.NE.0 ) THEN
            info = n + 4
            GO TO 10
         END IF
      END IF
*
      IF( ilvsr )
     $   CALL zlaset( 'Full', n, n, czero, cone, vsr, ldvsr )
*
*     Reduce to generalized Hessenberg form
*
      CALL zgghrd( jobvsl, jobvsr, n, ilo, ihi, a, lda, b, ldb, vsl,
     $             ldvsl, vsr, ldvsr, iinfo )
      IF( iinfo.NE.0 ) THEN
         info = n + 5
         GO TO 10
      END IF
*
*     Perform QZ algorithm, computing Schur vectors if desired
*
      iwork = itau
      CALL zhgeqz( 'S', jobvsl, jobvsr, n, ilo, ihi, a, lda, b, ldb,
     $             alpha, beta, vsl, ldvsl, vsr, ldvsr, work( iwork ),
     $             lwork+1-iwork, rwork( irwork ), iinfo )
      IF( iinfo.GE.0 )
     $   lwkopt = max( lwkopt, int( work( iwork ) )+iwork-1 )
      IF( iinfo.NE.0 ) THEN
         IF( iinfo.GT.0 .AND. iinfo.LE.n ) THEN
            info = iinfo
         ELSE IF( iinfo.GT.n .AND. iinfo.LE.2*n ) THEN
            info = iinfo - n
         ELSE
            info = n + 6
         END IF
         GO TO 10
      END IF
*
*     Apply permutation to VSL and VSR
*
      IF( ilvsl ) THEN
         CALL zggbak( 'P', 'L', n, ilo, ihi, rwork( ileft ),
     $                rwork( iright ), n, vsl, ldvsl, iinfo )
         IF( iinfo.NE.0 ) THEN
            info = n + 7
            GO TO 10
         END IF
      END IF
      IF( ilvsr ) THEN
         CALL zggbak( 'P', 'R', n, ilo, ihi, rwork( ileft ),
     $                rwork( iright ), n, vsr, ldvsr, iinfo )
         IF( iinfo.NE.0 ) THEN
            info = n + 8
            GO TO 10
         END IF
      END IF
*
*     Undo scaling
*
      IF( ilascl ) THEN
         CALL zlascl( 'U', -1, -1, anrmto, anrm, n, n, a, lda, iinfo )
         IF( iinfo.NE.0 ) THEN
            info = n + 9
            RETURN
         END IF
         CALL zlascl( 'G', -1, -1, anrmto, anrm, n, 1, alpha, n, iinfo)
         IF( iinfo.NE.0 ) THEN
            info = n + 9
            RETURN
         END IF
      END IF
*
      IF( ilbscl ) THEN
         CALL zlascl( 'U', -1, -1, bnrmto, bnrm, n, n, b, ldb, iinfo )
         IF( iinfo.NE.0 ) THEN
            info = n + 9
            RETURN
         END IF
         CALL zlascl( 'G', -1, -1, bnrmto, bnrm, n, 1, beta, n, iinfo )
         IF( iinfo.NE.0 ) THEN
            info = n + 9
            RETURN
         END IF
      END IF
*
   10 CONTINUE
      work( 1 ) = lwkopt
*
      RETURN
*
*     End of ZGEGS
*
      END
