*> \brief <b> ZGEEVX computes the eigenvalues and, optionally, the left
*and/or right eigenvectors for GE matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZGEGV + dependencies
*> <a
*href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgegv.f">
*> [TGZ]</a>
*> <a
*href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgegv.f">
*> [ZIP]</a>
*> <a
*href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgegv.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,
*                         VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO
*                         )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBVL, JOBVR
*       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
*      $                   BETA( * ), VL( LDVL, * ), VR( LDVR, * ),
*      $                   WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> This routine is deprecated and has been replaced by routine ZGGEV.
*>
*> ZGEGV computes the eigenvalues and, optionally, the left and/or
*right
*> eigenvectors of a complex matrix pair (A,B).
*> Given two square matrices A and B,
*> the generalized nonsymmetric eigenvalue problem (GNEP) is to find
*the
*> eigenvalues lambda and corresponding (non-zero) eigenvectors x such
*> that
*>    A*x = lambda*B*x.
*>
*> An alternate form is to find the eigenvalues mu and corresponding
*> eigenvectors y such that
*>    mu*A*y = B*y.
*>
*> These two forms are equivalent with mu = 1/lambda and x = y if
*> neither lambda nor mu is zero.  In order to deal with the case that
*> lambda or mu is zero or small, two values alpha and beta are
*returned
*> for each eigenvalue, such that lambda = alpha/beta and
*> mu = beta/alpha.
*>
*> The vectors x and y in the above equations are right eigenvectors of
*> the matrix pair (A,B).  Vectors u and v satisfying
*>    u**H*A = lambda*u**H*B  or  mu*v**H*A = v**H*B
*> are left eigenvectors of (A,B).
*>
*> Note: this routine performs "full balancing" on A and B
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBVL
*> \verbatim
*>          JOBVL is CHARACTER*1
*>          = 'N':  do not compute the left generalized eigenvectors;
*>          = 'V':  compute the left generalized eigenvectors (returned
*>                  in VL).
*> \endverbatim
*>
*> \param[in] JOBVR
*> \verbatim
*>          JOBVR is CHARACTER*1
*>          = 'N':  do not compute the right generalized eigenvectors;
*>          = 'V':  compute the right generalized eigenvectors
*(returned
*>                  in VR).
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A, B, VL, and VR.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA, N)
*>          On entry, the matrix A.
*>          If JOBVL = 'V' or JOBVR = 'V', then on exit A
*>          contains the Schur form of A from the generalized Schur
*>          factorization of the pair (A,B) after balancing.  If no
*>          eigenvectors were computed, then only the diagonal elements
*>          of the Schur form will be correct.  See ZGGHRD and ZHGEQZ
*>          for details.
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
*>          If JOBVL = 'V' or JOBVR = 'V', then on exit B contains the
*>          upper triangular matrix obtained from B in the generalized
*>          Schur factorization of the pair (A,B) after balancing.
*>          If no eigenvectors were computed, then only the diagonal
*>          elements of B will be correct.  See ZGGHRD and ZHGEQZ for
*>          details.
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
*>          GNEP.
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is COMPLEX*16 array, dimension (N)
*>          The complex scalars beta that define the eigenvalues of
*GNEP.
*>
*>          Together, the quantities alpha = ALPHA(j) and beta =
*BETA(j)
*>          represent the j-th eigenvalue of the matrix pair (A,B), in
*>          one of the forms lambda = alpha/beta or mu = beta/alpha.
*>          Since either lambda or mu may overflow, they should not,
*>          in general, be computed.
*> \endverbatim
*>
*> \param[out] VL
*> \verbatim
*>          VL is COMPLEX*16 array, dimension (LDVL,N)
*>          If JOBVL = 'V', the left eigenvectors u(j) are stored
*>          in the columns of VL, in the same order as their
*eigenvalues.
*>          Each eigenvector is scaled so that its largest component
*has
*>          abs(real part) + abs(imag. part) = 1, except for
*eigenvectors
*>          corresponding to an eigenvalue with alpha = beta = 0, which
*>          are set to zero.
*>          Not referenced if JOBVL = 'N'.
*> \endverbatim
*>
*> \param[in] LDVL
*> \verbatim
*>          LDVL is INTEGER
*>          The leading dimension of the matrix VL. LDVL >= 1, and
*>          if JOBVL = 'V', LDVL >= N.
*> \endverbatim
*>
*> \param[out] VR
*> \verbatim
*>          VR is COMPLEX*16 array, dimension (LDVR,N)
*>          If JOBVR = 'V', the right eigenvectors x(j) are stored
*>          in the columns of VR, in the same order as their
*eigenvalues.
*>          Each eigenvector is scaled so that its largest component
*has
*>          abs(real part) + abs(imag. part) = 1, except for
*eigenvectors
*>          corresponding to an eigenvalue with alpha = beta = 0, which
*>          are set to zero.
*>          Not referenced if JOBVR = 'N'.
*> \endverbatim
*>
*> \param[in] LDVR
*> \verbatim
*>          LDVR is INTEGER
*>          The leading dimension of the matrix VR. LDVR >= 1, and
*>          if JOBVR = 'V', LDVR >= N.
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
*>          blocksizes (for ZGEQRF, ZUNMQR, and ZUNGQR.)  Then compute:
*>          NB  -- MAX of the blocksizes for ZGEQRF, ZUNMQR, and
*ZUNGQR;
*>          The optimal LWORK is  MAX( 2*N, N*(NB+1) ).
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
*>          RWORK is DOUBLE PRECISION array, dimension (8*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          =1,...,N:
*>                The QZ iteration failed.  No eigenvectors have been
*>                calculated, but ALPHA(j) and BETA(j) should be
*>                correct for j=INFO+1,...,N.
*>          > N:  errors that usually indicate LAPACK problems:
*>                =N+1: error return from ZGGBAL
*>                =N+2: error return from ZGEQRF
*>                =N+3: error return from ZUNMQR
*>                =N+4: error return from ZUNGQR
*>                =N+5: error return from ZGGHRD
*>                =N+6: error return from ZHGEQZ (other than failed
*>                                               iteration)
*>                =N+7: error return from ZTGEVC
*>                =N+8: error return from ZGGBAK (computing VL)
*>                =N+9: error return from ZGGBAK (computing VR)
*>                =N+10: error return from ZLASCL (various calls)
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
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Balancing
*>  ---------
*>
*>  This driver calls ZGGBAL to both permute and scale rows and columns
*>  of A and B.  The permutations PL and PR are chosen so that PL*A*PR
*>  and PL*B*R will be upper triangular except for the diagonal blocks
*>  A(i:j,i:j) and B(i:j,i:j), with i and j as close together as
*>  possible.  The diagonal scaling matrices DL and DR are chosen so
*>  that the pair  DL*PL*A*PR*DR, DL*PL*B*PR*DR have elements close to
*>  one (except for the elements that start out zero.)
*>
*>  After the eigenvalues and eigenvectors of the balanced matrices
*>  have been computed, ZGGBAK transforms the eigenvectors back to what
*>  they would have been (in perfect arithmetic) if they had not been
*>  balanced.
*>
*>  Contents of A and B on Exit
*>  -------- -- - --- - -- ----
*>
*>  If any eigenvectors are computed (either JOBVL='V' or JOBVR='V' or
*>  both), then on exit the arrays A and B will contain the complex
*Schur
*>  form[*] of the "balanced" versions of A and B.  If no eigenvectors
*>  are computed, then only the diagonal blocks will be correct.
*>
*>  [*] In other words, upper triangular form.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE zgegv( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,
     $                  VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
*
*  -- LAPACK driver routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,
*  --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG
*  Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( lda, * ), ALPHA( * ), B( ldb, * ),
     $                   beta( * ), vl( ldvl, * ), vr( ldvr, * ),
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
      LOGICAL            ILIMIT, ILV, ILVL, ILVR, LQUERY
      CHARACTER          CHTEMP
      INTEGER            ICOLS, IHI, IINFO, IJOBVL, IJOBVR, ILEFT, ILO,
     $                   in, iright, irows, irwork, itau, iwork, jc, jr,
     $                   lopt, lwkmin, lwkopt, nb, nb1, nb2, nb3
      DOUBLE PRECISION   ABSAI, ABSAR, ABSB, ANRM, ANRM1, ANRM2, BNRM,
     $                   bnrm1, bnrm2, eps, safmax, safmin, salfai,
     $                   salfar, sbeta, scale, temp
      COMPLEX*16         X
*     ..
*     .. Local Arrays ..
      LOGICAL            LDUMMA( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zgeqrf, zggbak, zggbal, zgghrd, zhgeqz,
     $                   zlacpy, zlascl, zlaset, ztgevc, zungqr, zunmqr
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           lsame, ilaenv, dlamch, zlange
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, dimag, int, max
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
*     ..
*     .. Statement Function definitions ..
      abs1( x ) = abs( dble( x ) ) + abs( dimag( x ) )
*     ..
*     .. Executable Statements ..
*
*     Decode the input arguments
*
      IF( lsame( jobvl, 'N' ) ) THEN
         ijobvl = 1
         ilvl = .false.
      ELSE IF( lsame( jobvl, 'V' ) ) THEN
         ijobvl = 2
         ilvl = .true.
      ELSE
         ijobvl = -1
         ilvl = .false.
      END IF
*
      IF( lsame( jobvr, 'N' ) ) THEN
         ijobvr = 1
         ilvr = .false.
      ELSE IF( lsame( jobvr, 'V' ) ) THEN
         ijobvr = 2
         ilvr = .true.
      ELSE
         ijobvr = -1
         ilvr = .false.
      END IF
      ilv = ilvl .OR. ilvr
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
      ELSE IF( ldvl.LT.1 .OR. ( ilvl .AND. ldvl.LT.n ) ) THEN
         info = -11
      ELSE IF( ldvr.LT.1 .OR. ( ilvr .AND. ldvr.LT.n ) ) THEN
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
         lopt = max( 2*n, n*( nb+1 ) )
         work( 1 ) = lopt
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZGEGV ', -info )
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
      safmin = safmin + safmin
      safmax = one / safmin
*
*     Scale A
*
      anrm = zlange( 'M', n, n, a, lda, rwork )
      anrm1 = anrm
      anrm2 = one
      IF( anrm.LT.one ) THEN
         IF( safmax*anrm.LT.one ) THEN
            anrm1 = safmin
            anrm2 = safmax*anrm
         END IF
      END IF
*
      IF( anrm.GT.zero ) THEN
         CALL zlascl( 'G', -1, -1, anrm, one, n, n, a, lda, iinfo )
         IF( iinfo.NE.0 ) THEN
            info = n + 10
            RETURN
         END IF
      END IF
*
*     Scale B
*
      bnrm = zlange( 'M', n, n, b, ldb, rwork )
      bnrm1 = bnrm
      bnrm2 = one
      IF( bnrm.LT.one ) THEN
         IF( safmax*bnrm.LT.one ) THEN
            bnrm1 = safmin
            bnrm2 = safmax*bnrm
         END IF
      END IF
*
      IF( bnrm.GT.zero ) THEN
         CALL zlascl( 'G', -1, -1, bnrm, one, n, n, b, ldb, iinfo )
         IF( iinfo.NE.0 ) THEN
            info = n + 10
            RETURN
         END IF
      END IF
*
*     Permute the matrix to make it more nearly triangular
*     Also "balance" the matrix.
*
      ileft = 1
      iright = n + 1
      irwork = iright + n
      CALL zggbal( 'P', n, a, lda, b, ldb, ilo, ihi, rwork( ileft ),
     $             rwork( iright ), rwork( irwork ), iinfo )
      IF( iinfo.NE.0 ) THEN
         info = n + 1
         GO TO 80
      END IF
*
*     Reduce B to triangular form, and initialize VL and/or VR
*
      irows = ihi + 1 - ilo
      IF( ilv ) THEN
         icols = n + 1 - ilo
      ELSE
         icols = irows
      END IF
      itau = 1
      iwork = itau + irows
      CALL zgeqrf( irows, icols, b( ilo, ilo ), ldb, work( itau ),
     $             work( iwork ), lwork+1-iwork, iinfo )
      IF( iinfo.GE.0 )
     $   lwkopt = max( lwkopt, int( work( iwork ) )+iwork-1 )
      IF( iinfo.NE.0 ) THEN
         info = n + 2
         GO TO 80
      END IF
*
      CALL zunmqr( 'L', 'C', irows, icols, irows, b( ilo, ilo ), ldb,
     $             work( itau ), a( ilo, ilo ), lda, work( iwork ),
     $             lwork+1-iwork, iinfo )
      IF( iinfo.GE.0 )
     $   lwkopt = max( lwkopt, int( work( iwork ) )+iwork-1 )
      IF( iinfo.NE.0 ) THEN
         info = n + 3
         GO TO 80
      END IF
*
      IF( ilvl ) THEN
         CALL zlaset( 'Full', n, n, czero, cone, vl, ldvl )
         CALL zlacpy( 'L', irows-1, irows-1, b( ilo+1, ilo ), ldb,
     $                vl( ilo+1, ilo ), ldvl )
         CALL zungqr( irows, irows, irows, vl( ilo, ilo ), ldvl,
     $                work( itau ), work( iwork ), lwork+1-iwork,
     $                iinfo )
         IF( iinfo.GE.0 )
     $      lwkopt = max( lwkopt, int( work( iwork ) )+iwork-1 )
         IF( iinfo.NE.0 ) THEN
            info = n + 4
            GO TO 80
         END IF
      END IF
*
      IF( ilvr )
     $   CALL zlaset( 'Full', n, n, czero, cone, vr, ldvr )
*
*     Reduce to generalized Hessenberg form
*
      IF( ilv ) THEN
*
*        Eigenvectors requested -- work on whole matrix.
*
         CALL zgghrd( jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, vl,
     $                ldvl, vr, ldvr, iinfo )
      ELSE
         CALL zgghrd( 'N', 'N', irows, 1, irows, a( ilo, ilo ), lda,
     $                b( ilo, ilo ), ldb, vl, ldvl, vr, ldvr, iinfo )
      END IF
      IF( iinfo.NE.0 ) THEN
         info = n + 5
         GO TO 80
      END IF
*
*     Perform QZ algorithm
*
      iwork = itau
      IF( ilv ) THEN
         chtemp = 'S'
      ELSE
         chtemp = 'E'
      END IF
      CALL zhgeqz( chtemp, jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb,
     $             alpha, beta, vl, ldvl, vr, ldvr, work( iwork ),
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
         GO TO 80
      END IF
*
      IF( ilv ) THEN
*
*        Compute Eigenvectors
*
         IF( ilvl ) THEN
            IF( ilvr ) THEN
               chtemp = 'B'
            ELSE
               chtemp = 'L'
            END IF
         ELSE
            chtemp = 'R'
         END IF
*
         CALL ztgevc( chtemp, 'B', ldumma, n, a, lda, b, ldb, vl, ldvl,
     $                vr, ldvr, n, in, work( iwork ), rwork( irwork ),
     $                iinfo )
         IF( iinfo.NE.0 ) THEN
            info = n + 7
            GO TO 80
         END IF
*
*        Undo balancing on VL and VR, rescale
*
         IF( ilvl ) THEN
            CALL zggbak( 'P', 'L', n, ilo, ihi, rwork( ileft ),
     $                   rwork( iright ), n, vl, ldvl, iinfo )
            IF( iinfo.NE.0 ) THEN
               info = n + 8
               GO TO 80
            END IF
            DO 30 jc = 1, n
               temp = zero
               DO 10 jr = 1, n
                  temp = max( temp, abs1( vl( jr, jc ) ) )
   10          CONTINUE
               IF( temp.LT.safmin )
     $            GO TO 30
               temp = one / temp
               DO 20 jr = 1, n
                  vl( jr, jc ) = vl( jr, jc )*temp
   20          CONTINUE
   30       CONTINUE
         END IF
         IF( ilvr ) THEN
            CALL zggbak( 'P', 'R', n, ilo, ihi, rwork( ileft ),
     $                   rwork( iright ), n, vr, ldvr, iinfo )
            IF( iinfo.NE.0 ) THEN
               info = n + 9
               GO TO 80
            END IF
            DO 60 jc = 1, n
               temp = zero
               DO 40 jr = 1, n
                  temp = max( temp, abs1( vr( jr, jc ) ) )
   40          CONTINUE
               IF( temp.LT.safmin )
     $            GO TO 60
               temp = one / temp
               DO 50 jr = 1, n
                  vr( jr, jc ) = vr( jr, jc )*temp
   50          CONTINUE
   60       CONTINUE
         END IF
*
*        End of eigenvector calculation
*
      END IF
*
*     Undo scaling in alpha, beta
*
*     Note: this does not give the alpha and beta for the unscaled
*     problem.
*
*     Un-scaling is limited to avoid underflow in alpha and beta
*     if they are significant.
*
      DO 70 jc = 1, n
         absar = abs( dble( alpha( jc ) ) )
         absai = abs( dimag( alpha( jc ) ) )
         absb = abs( dble( beta( jc ) ) )
         salfar = anrm*dble( alpha( jc ) )
         salfai = anrm*dimag( alpha( jc ) )
         sbeta = bnrm*dble( beta( jc ) )
         ilimit = .false.
         scale = one
*
*        Check for significant underflow in imaginary part of ALPHA
*
         IF( abs( salfai ).LT.safmin .AND. absai.GE.
     $       max( safmin, eps*absar, eps*absb ) ) THEN
            ilimit = .true.
            scale = ( safmin / anrm1 ) / max( safmin, anrm2*absai )
         END IF
*
*        Check for significant underflow in real part of ALPHA
*
         IF( abs( salfar ).LT.safmin .AND. absar.GE.
     $       max( safmin, eps*absai, eps*absb ) ) THEN
            ilimit = .true.
            scale = max( scale, ( safmin / anrm1 ) /
     $              max( safmin, anrm2*absar ) )
         END IF
*
*        Check for significant underflow in BETA
*
         IF( abs( sbeta ).LT.safmin .AND. absb.GE.
     $       max( safmin, eps*absar, eps*absai ) ) THEN
            ilimit = .true.
            scale = max( scale, ( safmin / bnrm1 ) /
     $              max( safmin, bnrm2*absb ) )
         END IF
*
*        Check for possible overflow when limiting scaling
*
         IF( ilimit ) THEN
            temp = ( scale*safmin )*max( abs( salfar ), abs( salfai ),
     $             abs( sbeta ) )
            IF( temp.GT.one )
     $         scale = scale / temp
            IF( scale.LT.one )
     $         ilimit = .false.
         END IF
*
*        Recompute un-scaled ALPHA, BETA if necessary.
*
         IF( ilimit ) THEN
            salfar = ( scale*dble( alpha( jc ) ) )*anrm
            salfai = ( scale*dimag( alpha( jc ) ) )*anrm
            sbeta = ( scale*beta( jc ) )*bnrm
         END IF
         alpha( jc ) = dcmplx( salfar, salfai )
         beta( jc ) = sbeta
   70 CONTINUE
*
   80 CONTINUE
      work( 1 ) = lwkopt
*
      RETURN
*
*     End of ZGEGV
*
      END
