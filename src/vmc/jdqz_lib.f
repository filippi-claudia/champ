      subroutine error (m)
c
c     Coded by Diederik R. Fokkema
c
c     $Id: error.f,v 1.4 1995/07/26 09:26:26 fokkema Exp $
c
c     .. Parameters ..
c
      implicit none
      character m*(*)

      write(6,*) 'HAHA'
ctex@ \begin{manpage}{ERROR} 
ctex@ \subtitle{Name}
ctex@    ERROR --- Type an error message and stop
ctex@
ctex@ \subtitle{Declaration}
ctex@    %declaration
ctex@ \subtitle{Parameters}
ctex@    \variable{m}
ctex@       character string. On entry m must contain the message string.
ctex@
ctex@ \subtitle{Description}
ctex@    This subroutine types an error message and stops.
ctex@
ctex@ \end{manpage}
ctex@ \begin{verbatim}
ctex@    % actual code
ctex@ \end{verbatim}
c
c     .. Local ..
c
c     None.
c
c     .. Called subroutines
c
c     None.
c
c     .. Executable statements
c
      print 10, m
 10   format (/,1x,'Error: ',a,/)
c
c     --- Stop
c
      stop
      end
      subroutine makemm (n, k, w, v, m, zm, ldm)
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/08/03 23:33:20 caveman>
c
c     .. Parameters ..
c
      implicit none
      integer n, k, ldm
      double complex w(n,*), v(n,*), m(ldm,*), zm(ldm,*)
c
c     .. Local ..
c
      integer i, j
      double complex zdotc
c
c     .. Executable statements ..
c
      do i=1,k
         do j=1,k
            if (i.eq.k.or.j.eq.k)
     $           m(i,j) = zdotc (n, w(1,i), 1, v(1,j), 1)
            zm(i,j) = m(i,j)
         enddo
      enddo
c
c     --- Return
c
      end
      subroutine mkqkz (n, k, q, kq, qkz, invqkz, ldqkz, ipiv)
c
c     Coded by Diederik Fokkema
c
c     $Id: mkqkz.f,v 1.1 1995/08/05 09:08:21 caveman Exp $
c
c     Time-stamp: <95/08/03 23:52:52 caveman>
c
c     .. Parameters ..
c
      implicit none
      integer n, k, ldqkz, ipiv(*)
      double complex q(n,*), kq(n,*), qkz(ldqkz,*), invqkz(ldqkz,*)
c
c     .. local ..
c
      integer i, j, info
      double complex zdotc
c
c     .. Executable statements ..
c
      do i=1,k
         do j=1,k
            if (i.eq.k.or.j.eq.k) qkz(i,j) =
     $           zdotc (n, q(1,i), 1, kq(1,j), 1)
            invqkz(i,j) = qkz(i,j)
         enddo
      enddo
      call zgetrf (k, k, invqkz, ldqkz, ipiv, info)
c
c     --- Return
c
      end

      subroutine myexc (n, s, t, z, q, ldz, ifst, ilst)
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/10/31 23:51:12 caveman>
c
c     .. Parameters ..
c
      implicit none
      integer n, ldz, ifst, ilst
      double complex s(ldz,*), t(ldz,*), z(ldz,*), q(ldz,*)
c
c     .. Local ..
c
      logical tlts
      integer k, m1, m2, m3
      double complex f, f1, f2, c1, c2, r, sn
      real*8 cs
c
c     .. Executable statements ..
c
      if (n.eq.1 .or. ifst.eq.ilst) return
      if (ifst.lt.ilst) then
         m1 = 0
         m2 = -1
         m3 = 1
      else
         m1 = -1
         m2 = 0
         m3 = -1
      endif
      do k = ifst+m1, ilst+m2, m3
         f = max(abs(t(k+1,k+1)),abs(s(k+1,k+1)))
         f1 = t(k+1,k+1)/f
         f2 = s(k+1,k+1)/f
         tlts = .true.
         if (abs(f1).gt.abs(f2)) tlts = .false.
         c1 = f1*s(k,k) - f2*t(k,k)
         c2 = f1*s(k,k+1) - f2*t(k,k+1)
         call zlartg (c2, c1, cs, sn, r)
         call zrot (k+1, s(1,k+1), 1, s(1,k), 1, cs, sn)
         call zrot (k+1, t(1,k+1), 1, t(1,k), 1, cs, sn)
         call zrot (n, q(1,k+1), 1, q(1,k), 1, cs, sn)
         if (tlts) then
            c1 = s(k,k)
            c2 = s(k+1,k) 
         else
            c1 = t(k,k)
            c2 = t(k+1,k) 
         endif
         call zlartg (c1, c2, cs, sn, r)
         call zrot (n-k+1, s(k,k), ldz, s(k+1,k), ldz, cs, sn)
         call zrot (n-k+1, t(k,k), ldz, t(k+1,k), ldz, cs, sn)
         call zrot (n, z(1,k), 1, z(1,k+1), 1, cs, dconjg(sn))
      enddo
c
c     --- Return
c
      end
      subroutine psolve (n, x, nq, q, kz, invqkz, ldqkz, ipiv, f)
c
c     Coded by Diederik Fokkema
c
c     .. Parameters ..
c
      implicit none
      integer n, nq, ldqkz, ipiv(*)
      double complex x(*), q(n,*), kz(n,*),
     $     invqkz(ldqkz,*), f(*)
c
c     .. local .. 
c
      integer info
      double complex zero, one
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))
c
c     .. Executable Statements ..
c
      call precon( n, x )
      call zgemv ('c', n, nq, one, q, n, x, 1, zero, f, 1)
      call zgetrs('n', nq, 1, invqkz, ldqkz, ipiv, f, ldqkz, info)
      call zgemv ('n', n, nq, -one, kz, n, f, 1, one, x, 1)
c
c     --- Return
c
      end


      subroutine qzsort (ta, tb, k, s, t, z, q, ldz, alpha, beta,
     $     order)
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/08/03 23:34:03 caveman>
c
c     .. Parameters ..
c
      implicit none
      integer k, ldz, order
      double complex ta, tb, s(ldz,*), t(ldz,*), z(ldz,*),
     $     q(ldz,*), alpha(*), beta(*)
c
c     .. Local ..
c
      integer i, j, select
c
c     .. Executable statements ..
c
      do i = 1,k
         do j = 1,k
            alpha(j) = s(j,j)
            beta(j) = t(j,j)
         enddo
         j = select (k-i+1, ta, tb, alpha(i), beta(i), order) + i-1
         call myexc (k, s, t, z, q, ldz, j, i)
      enddo
c
c     --- Return
c
      end
      integer function select (n, sa, sb, a, b, order)
c
c     Coded by Diederik Fokkema
c     Modified Martin van Gijzen, test on division by zero
c
c     .. Parameters ..
c
      implicit none
      integer n, order
      double complex sa, sb, a(*), b(*)
c
c     .. Local ..
c
      integer i, j
      real*8 dtmp, optval
c
c     .. Executable statements ..
c
      j = 1
      if ( order .le. 0 ) then
         optval =  1.d99
      else
         optval = -1.d99
      end if
      if (order.eq.0) then
c
c...        Nearest to target
c
         do i=1,n
            if ( b(i) .ne. 0.d0 ) then
               dtmp = abs(a(i)/b(i)-sa/sb)
               if (dtmp.lt.optval) then
                  j=i
                  optval = dtmp
               end if
            end if
         enddo
      elseif (order.eq.-1) then
c
c...        Smallest real part
c
         do i=1,n
            if ( b(i) .eq. 0.d0 ) then
               if ( dble(a(i)) .lt. 0.d0 ) then
                  j=i
                  optval = -1.d99
               end if
            else
               dtmp = dble(a(i)/b(i))
               if (dtmp.lt.optval) then
                  j=i
                  optval = dtmp
               end if
            end if
         enddo
      elseif (order.eq.1) then
c
c...        Largest real part
c
         do i=1,n
            if ( b(i) .eq. 0.d0 ) then
               if ( dble(a(i)) .gt. 0.d0 ) then
                  j=i
                  optval = 1.d99
               end if
            else
               dtmp = dble(a(i)/b(i))
               if (dtmp.gt.optval) then
                  j=i
                  optval = dtmp
               end if
            end if
         enddo
      elseif (order.eq.-2) then
c
c...        Smallest imaginari part
c
         do i=1,n
            if ( b(i) .eq. 0.d0 ) then
               if ( imag(a(i)) .lt. 0.d0 ) then
                  j=i
                  optval = -1.d99
               end if
            else
               dtmp = imag(a(i)/b(i))
               if (dtmp.lt.optval) then
                  j=i
                  optval = dtmp
               end if
            end if
         enddo
      elseif (order.eq.2) then
c
c...        Largest imaginari part
c
         do i=1,n
            if ( b(i) .eq. 0.d0 ) then
               if ( imag(a(i)) .gt. 0.d0 ) then
                  j=i
                  optval = 1.d99
               end if
            else
               dtmp = imag(a(i)/b(i))
               if (dtmp.gt.optval) then
                  j=i
                  optval = dtmp
               end if
            end if
         enddo
      else
         call error ('unknown order in select')
      endif
      select = j
c
c     --- Return
c
      end
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
      subroutine zmgs (n, k, v, w, job )
c
c     Coded by Diederik Fokkema
c     Modified 05-23-96: M. Kooper: job =1 corrected, array YWORK added
c
c     .. Parameters ..
c
      implicit none
      integer n, k, job
      double complex v(n,*), w(*)
c
c     .. Local ..
c
      integer i
      real*8 s0, s1, dznrm2
      double complex znrm
c
c     .. Executable statements ..
c
      s1 = dznrm2 (n, w, 1)
      do i=1, k
	 s0 = s1
	 call zortho (n, v(1,i), w, s0, s1, znrm)
      enddo
      if (job.eq.0) then
	 return
      else
         znrm  = 1.d0/s1
         call zscal (n, znrm, w, 1)
      endif
c
c     --- Return
c
      return
      end

      subroutine zortho (n, v, w, s0, s1, znrm)
c
c     Coded by Diederik Fokkema
c
c     .. Parameters ..
c
      implicit none
      integer n
      real*8 s0, s1
      double complex v(*), w(*), znrm
c
c     .. Local ..
c
      real*8 kappa, dznrm2
      double complex ztmp, zdotc
      parameter (kappa=1d2)
c
c     .. Executable statements ..
c
      znrm = zdotc (n, v, 1, w, 1)
      call zaxpy (n, (-znrm), v, 1, w, 1)
      s1 = dznrm2 (n, w, 1)
      if (s1.gt.s0/kappa) then
	 goto 100
      else
	 s0 = s1
	 ztmp = zdotc (n, v, 1, w, 1)
         znrm = znrm + ztmp
	 call zaxpy (n, (-ztmp), v, 1, w, 1)
	 s1 = dznrm2 (n, w, 1)
	 if (s1.gt.s0/kappa) then
	    goto 100
	 else
	    call error ('zero vector in zmgs')
	 endif
      endif
 100  continue
c
c     --- Return
c
      return
      end

      subroutine zones (n, x)
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/07/30 21:45:46 caveman>
c
c     .. Parameters ..
c
      implicit none
      integer n
      double complex x(*)
c
c     .. Local ..
c
      integer i
c
c     .. Executable statements ..
c
      do i=1,n
         x(i) = (1.0d0,0.0d0)
      enddo
c
c     --- Return
c
      end

      subroutine zxpay(n,dx,incx,da,dy,incy)
c
c     modified by:  D.R. Fokkema
c     01/06/94
c
c     a vector plus constant times a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none
      double complex dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        dy(iy) = da*dy(iy) + dx(ix)
        ix = ix + incx
        iy = iy + incy
      enddo
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dy(i) = da*dy(i) + dx(i)
      enddo
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do i = mp1,n,4
        dy(i) = da*dy(i) + dx(i)
        dy(i + 1) = da*dy(i + 1) + dx(i + 1)
        dy(i + 2) = da*dy(i + 2) + dx(i + 2)
        dy(i + 3) = da*dy(i + 3) + dx(i + 3)
      enddo
      return
      end
      subroutine zzeros (n, x)
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/07/30 21:48:00 caveman>
c
c     .. Parameters ..
c
      implicit none
      integer n
      double complex x(*)
c
c     .. Local ..
c
      integer i
c
c     .. Executable statements ..
c
      do i=1,n
         x(i) = (0.0d0,0.0d0)
      enddo
c
c     --- Return
c
      end

