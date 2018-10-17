      SUBROUTINE V_SPLINE(k_bc1,k_bcn,n,x,f,wk)
!***********************************************************************
!V_SPLINE evaluates the coefficients for a 1d cubic interpolating spline
!References:
!  Forsythe, Malcolm, Moler, Computer Methods for Mathematical
!    Computations, Prentice-Hall, 1977, p.76
!  Engeln-Muellges, Uhlig, Numerical Algorithms with Fortran, Springer,
!    1996, p.251
!  W.A.Houlberg, D.McCune 3/2000
!Input:
!  k_bc1-option for BC at x(1)
!       =-1 periodic, ignore k_bcn
!       =0 not-a-knot
!       =1 s'(x1) = input value of f(2,1)
!       =2 s''(x1) = input value of f(3,1)
!       =3 s'(x1) = 0.0
!       =4 s''(x1) = 0.0
!       =5 match first derivative to first 2 points
!       =6 match second derivative to first 3 points
!       =7 match third derivative to first 4 points
!       =else use not-a-knot
!  k_bcn-option for boundary condition at x(n)
!       =0 not-a-knot
!       =1 s'(x1) = input value of f(2,1)
!       =2 s''(x1) = input value of f(3,1)
!       =3 s'(x1) = 0.0
!       =4 s''(x1) = 0.0
!       =5 match first derivative to first 2 points
!       =6 match second derivative to first 3 points
!       =7 match third derivative to first 4 points
!       =else use knot-a-knot
!  n-number of data points or knots-(n.ge.2)
!  x(n)-abscissas of the knots in strictly increasing order
!  f(1,i)-ordinates of the knots
!  f(2,1)-input value of s'(x1) for k_bc1=1
!  f(2,n)-input value of s'(xn) for k_bcn=1
!  f(3,1)-input value of s''(x1) for k_bc1=2
!  f(3,n)-input value of s''(xn) for k_bcn=2
!  wk(n)-scratch work area for periodic BC
!Output:
!  f(2,i)=s'(x(i))
!  f(3,i)=s''(x(i))
!  f(4,i)=s'''(x(i))
!Comments:
!  s(x)=f(1,i)+f(2,i)*(x-x(i))+f(3,i)*(x-x(i))**2/2!
!       +f(4,i)*(x-x(i))**3/3! for x(i).le.x.le.x(i+1)
!  W_SPLINE can be used to evaluate the spline and its derivatives
!  The cubic spline is twice differentiable (C2)
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        k_bc1,                   k_bcn,
     &               n
      REAL           x(*),                    wk(*),
     &               f(4,*)
!Declaration in local variables
      INTEGER        i,                       ib,
     &               imax,                    imin
      REAL           a1,                      an,
     &               b1,                      bn,
     &               q,                       t,
     &               hn
!Set default range
      imin=1
      imax=n
!Set first and second BC values
      a1=0.0
      b1=0.0
      an=0.0
      bn=0.0
      IF(k_bc1.eq.1) THEN
        a1=f(2,1)
      ELSEIF(k_bc1.eq.2) THEN
        b1=f(3,1)
      ELSEIF(k_bc1.eq.5) THEN
        a1=(f(1,2)-f(1,1))/(x(2)-x(1))
      ELSEIF(k_bc1.eq.6) THEN
        b1=2.0*((f(1,3)-f(1,2))/(x(3)-x(2))
     &         -(f(1,2)-f(1,1))/(x(2)-x(1)))/(x(3)-x(1))
      ENDIF
      IF(k_bcn.eq.1) THEN
        an=f(2,n)
      ELSEIF(k_bcn.eq.2) THEN
        bn=f(3,n)
      ELSEIF(k_bcn.eq.5) THEN
        an=(f(1,n)-f(1,n-1))/(x(n)-x(n-1))
      ELSEIF(k_bcn.eq.6) THEN
        bn=2.0*((f(1,n)-f(1,n-1))/(x(n)-x(n-1))
     &         -(f(1,n-1)-f(1,n-2))/(x(n-1)-x(n-2)))/(x(n)-x(n-2))
      ENDIF
!Clear f(2:4,n)
      f(2,n)=0.0
      f(3,n)=0.0
      f(4,n)=0.0
      IF(n.eq.2) THEN
!Coefficients for n=2
        f(2,1)=(f(1,2)-f(1,1))/(x(2)-x(1))
        f(3,1)=0.0
        f(4,1)=0.0
        f(2,2)=f(2,1)
        f(3,2)=0.0
        f(4,2)=0.0
      ELSEIF(n.gt.2) THEN
!Set up tridiagonal system for A*y=B where y(i) are the second
!  derivatives at the knots
!  f(2,i) are the diagonal elements of A
!  f(4,i) are the off-diagonal elements of A
!  f(3,i) are the B elements/3, and will become c/3 upon solution
        f(4,1)=x(2)-x(1)
        f(3,2)=(f(1,2)-f(1,1))/f(4,1)
        DO i=2,n-1
          f(4,i)=x(i+1)-x(i)
          f(2,i)=2.0*(f(4,i-1)+f(4,i))
          f(3,i+1)=(f(1,i+1)-f(1,i))/f(4,i)
          f(3,i)=f(3,i+1)-f(3,i)
        ENDDO
!  BC's
!    Left
        IF(k_bc1.eq.-1) THEN
          f(2,1)=2.0*(f(4,1)+f(4,n-1))
          f(3,1)=(f(1,2)-f(1,1))/f(4,1)-(f(1,n)-f(1,n-1))/f(4,n-1)
          wk(1)=f(4,n-1)
          DO i=2,n-3
            wk(i)=0.0
          ENDDO
          wk(n-2)=f(4,n-2)
          wk(n-1)=f(4,n-1)
        ELSEIF(k_bc1.eq.1.or.k_bc1.eq.3.or.k_bc1.eq.5) THEN
          f(2,1)=2.0*f(4,1)
          f(3,1)=(f(1,2)-f(1,1))/f(4,1)-a1
        ELSEIF(k_bc1.eq.2.or.k_bc1.eq.4.or.k_bc1.eq.6) THEN
          f(2,1)=2.0*f(4,1)
          f(3,1)=f(4,1)*b1/3.0
          f(4,1)=0.0
        ELSEIF(k_bc1.eq.7) THEN
          f(2,1)=-f(4,1)
          f(3,1)=f(3,3)/(x(4)-x(2))-f(3,2)/(x(3)-x(1))
          f(3,1)=f(3,1)*f(4,1)**2/(x(4)-x(1))
        ELSE                             ! not a knot:
          imin=2
          f(2,2)=f(4,1)+2.0*f(4,2)
          f(3,2)=f(3,2)*f(4,2)/(f(4,1)+f(4,2))
        ENDIF
!    Right
        IF(k_bcn.eq.1.or.k_bcn.eq.3.or.k_bcn.eq.5) THEN
          f(2,n)=2.0*f(4,n-1)
          f(3,n)=-(f(1,n)-f(1,n-1))/f(4,n-1)+an
        ELSEIF(k_bcn.eq.2.or.k_bcn.eq.4.or.k_bcn.eq.6) THEN
          f(2,n)=2.0*f(4,n-1)
          f(3,n)=f(4,n-1)*bn/3.0
          f(4,n-1)=0.0
        ELSEIF(k_bcn.eq.7) THEN
          f(2,n)=-f(4,n-1)
          f(3,n)=f(3,n-1)/(x(n)-x(n-2))-f(3,n-2)/(x(n-1)-x(n-3))
          f(3,n)=-f(3,n)*f(4,n-1)**2/(x(n)-x(n-3))
        ELSEIF(k_bc1.ne.-1) THEN         ! not a knot:
          imax=n-1
          f(2,n-1)=2.0*f(4,n-2)+f(4,n-1)
          f(3,n-1)=f(3,n-1)*f(4,n-2)/(f(4,n-1)+f(4,n-2))
        ENDIF
!  Limit solution for only three points in domain
        IF(n.eq.3) THEN
          f(3,1)=0.0
          f(3,n)=0.0
        ENDIF
        IF(k_bc1.eq.-1) THEN
!Solve system of equations for second derivatives at the knots
!  Periodic BC
!    Forward elimination
          DO i=2,n-2
            t=f(4,i-1)/f(2,i-1)
            f(2,i)=f(2,i)-t*f(4,i-1)
            f(3,i)=f(3,i)-t*f(3,i-1)
            wk(i)=wk(i)-t*wk(i-1)
            q=wk(n-1)/f(2,i-1)
            wk(n-1)=-q*f(4,i-1)
            f(2,n-1)=f(2,n-1)-q*wk(i-1)
            f(3,n-1)=f(3,n-1)-q*f(3,i-1)
          ENDDO
!    Correct the n-1 element
          wk(n-1)=wk(n-1)+f(4,n-2)
!    Complete the forward elimination
!    wk(n-1) and wk(n-2) are the off-diag elements of the lower corner
          t=wk(n-1)/f(2,n-2)
          f(2,n-1)=f(2,n-1)-t*wk(n-2)
          f(3,n-1)=f(3,n-1)-t*f(3,n-2)
!    Back substitution
          f(3,n-1)=f(3,n-1)/f(2,n-1)
          f(3,n-2)=(f(3,n-2)-wk(n-2)*f(3,n-1))/f(2,n-2)
          DO ib=3,n-1
            i=n-ib
            f(3,i)=(f(3,i)-f(4,i)*f(3,i+1)-wk(i)*f(3,n-1))/f(2,i)
          ENDDO
          f(3,n)=f(3,1)
        ELSE
!  Non-periodic BC
!    Forward elimination
!    For Not-A-Knot BC the off-diagonal end elements are not equal
          DO i=imin+1,imax
            IF((i.eq.n-1).and.(imax.eq.n-1)) THEN
              t=(f(4,i-1)-f(4,i))/f(2,i-1)
            ELSE
              t=f(4,i-1)/f(2,i-1)
            ENDIF
            IF((i.eq.imin+1).and.(imin.eq.2)) THEN
              f(2,i)=f(2,i)-t*(f(4,i-1)-f(4,i-2))
            ELSE
              f(2,i)=f(2,i)-t*f(4,i-1)
            ENDIF
            f(3,i)=f(3,i)-t*f(3,i-1)
          ENDDO
!    Back substitution
          f(3,imax)=f(3,imax)/f(2,imax)
          DO ib=1,imax-imin
            i=imax-ib
            IF((i.eq.2).and.(imin.eq.2)) THEN
              f(3,i)=(f(3,i)-(f(4,i)-f(4,i-1))*f(3,i+1))/f(2,i)
            ELSE
              f(3,i)=(f(3,i)-f(4,i)*f(3,i+1))/f(2,i)
            ENDIF
          ENDDO
!    Reset d array to step size
          f(4,1)=x(2)-x(1)
          f(4,n-1)=x(n)-x(n-1)
!    Set f(3,1) for not-a-knot
          IF(k_bc1.le.0.or.k_bc1.gt.5) THEN
            f(3,1)=(f(3,2)*(f(4,1)+f(4,2))-f(3,3)*f(4,1))/f(4,2)
          ENDIF
!    Set f(3,n) for not-a-knot
          IF(k_bcn.le.0.or.k_bcn.gt.5) THEN
            f(3,n)=f(3,n-1)+(f(3,n-1)-f(3,n-2))*f(4,n-1)/f(4,n-2)
          ENDIF
        ENDIF
!f(3,i) is now the sigma(i) of the text and f(4,i) is the step size
!Compute polynomial coefficients
        DO i=1,n-1
          f(2,i)=
     >        (f(1,i+1)-f(1,i))/f(4,i)-f(4,i)*(f(3,i+1)+2.0*f(3,i))
          f(4,i)=(f(3,i+1)-f(3,i))/f(4,i)
          f(3,i)=6.0*f(3,i)
          f(4,i)=6.0*f(4,i)
        ENDDO
        IF(k_bc1.eq.-1) THEN
          f(2,n)=f(2,1)
          f(3,n)=f(3,1)
          f(4,n)=f(4,1)
        ELSE
           hn=x(n)-x(n-1)
           f(2,n)=f(2,n-1)+hn*(f(3,n-1)+0.5*hn*f(4,n-1))
           f(3,n)=f(3,n-1)+hn*f(4,n-1)
           f(4,n)=f(4,n-1)
           IF(k_bcn.eq.1.or.k_bcn.eq.3.or.k_bcn.eq.5) THEN
              f(2,n)=an
           ELSE IF(k_bcn.eq.2.or.k_bcn.eq.4.or.k_bcn.eq.6) THEN
              f(3,n)=bn
           ENDIF
        ENDIF
      ENDIF
      RETURN
      END
