      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      use precision_kinds, only: dp
!  From Numerical Recipes
      implicit none

, n, nmax
      real(dp) :: I, K, N, P, QN
      real(dp) :: SIG, SPLINE, SUBROUTINE, U
      real(dp) :: UN, X, Y, Y2
      real(dp) :: YP1, YPN
      real(dp), dimension(N) :: X
      real(dp), dimension(N) :: Y
      real(dp), dimension(N) :: Y2
      real(dp), dimension(NMAX) :: U
      real(dp), parameter :: NMAX = 5001

!  Compilers that do not allow automatic arrays:
!  Compilers that do allow automatic arrays:
!      DIMENSION X(N),Y(N),Y2(N),U(N)

!  Set lower boundary cond to be natural or to have specified 1st deriv
      IF (YP1.GT..99D30) THEN
        Y2(1)=0
        U(1)=0
       ELSE
        Y2(1)=-0.5D0
        U(1)=(3.D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF

!  Decomposition loop of tridiagonal algorithm
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2
        Y2(I)=(SIG-1)/P
        U(I)=(6*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     &  /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
   11 CONTINUE

!  Set upper boundary cond to be natural or to have specified 1st deriv
      IF (YPN.GT..99D30) THEN
        QN=0
        UN=0
      ELSE
        QN=0.5D0
        UN=(3/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF

!  Backsubstitution loop of tridiagonal algorithm
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
   12 CONTINUE

      RETURN
      END

c-----------------------------------------------------------------------

      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
!  From Numerical Recipes
      use precision_kinds, only: dp
      implicit none


      real(dp) :: A, B, H, K, KHI
      real(dp) :: KLO, N, SPLINT, SUBROUTINE
      real(dp) :: X, XA, Y, Y2A
      real(dp) :: YA
      DIMENSION XA(N),YA(N),Y2A(N)

!  Find right place in table by binary search
      KLO=1
      KHI=N
    1 IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF

      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) call fatal_error ('The x values in splint are not distinct')
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     &      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6

      RETURN
      END


