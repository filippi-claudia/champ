*DECK ICAMAX
      INTEGER FUNCTION ICAMAX (N, CX, INCX)
C***BEGIN PROLOGUE  ICAMAX
C***PURPOSE  Find the smallest index of the component of a complex
C            vector having the maximum sum of magnitudes of real
C            and imaginary parts.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A2
C***TYPE      COMPLEX (ISAMAX-S, IDAMAX-D, ICAMAX-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       CX  complex vector with N elements
C     INCX  storage spacing between elements of CX
C
C     --Output--
C   ICAMAX  smallest index (zero if N .LE. 0)
C
C     Returns the smallest index of the component of CX having the
C     largest sum of magnitudes of real and imaginary parts.
C     ICAMAX = first I, I = 1 to N, to maximize
C     ABS(REAL(CX(IX+(I-1)*INCX))) + ABS(IMAG(CX(IX+(I-1)*INCX))),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  ICAMAX
      COMPLEX CX(*)
      REAL SMAX, XMAG
      INTEGER I, INCX, IX, N
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C***FIRST EXECUTABLE STATEMENT  ICAMAX
      ICAMAX = 0
      IF (N .LE. 0) RETURN
      ICAMAX = 1
      IF (N .EQ. 1) RETURN
C
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      SMAX = CABS1(CX(IX))
      IX = IX + INCX
      DO 10 I = 2,N
        XMAG = CABS1(CX(IX))
        IF (XMAG .GT. SMAX) THEN
          ICAMAX = I
          SMAX = XMAG
        ENDIF
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
   20 SMAX = CABS1(CX(1))
      DO 30 I = 2,N
        XMAG = CABS1(CX(I))
        IF (XMAG .GT. SMAX) THEN
          ICAMAX = I
          SMAX = XMAG
        ENDIF
   30 CONTINUE
      RETURN
      END
