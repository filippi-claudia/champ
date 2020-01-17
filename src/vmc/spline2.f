      subroutine spline2(x, y, n, yp1, ypn, y2, work)
C Cubic spline routine based on p. 88 of Numerical Recipes.
C Given arrays x and y of length n containing y=f(x) with x's in
C ascending order and given yp1 and ypn for first derivative of interpolating
C function at the endpoints, returns array y2 of length n which contains
C the second derivatives of the interpolating function at the tabulated
C points x.  If yp1 and/or ypn are 1.e30 or larger, routine sets corresponding
C boundary condition for a natural spline, with zero second derivative on
C that boundary.
C The cubic spline fit to the function is then given by
C
C  y = A y  + B y    + C y'' + D y''
C	  j	 j+1	  j	  j+1
C
C with A=(x(j+1)-x)/(x(j+1)-x(j)), B=1-A=(x-x(j))/(x(j+1)-x(j)),
C C=(A^3-A)(x(j+1)-x(j))^2/6, and D=(B^3-B)(x(j+1)-x(j))^2/6.
C
C The first derivative is therefore (with dx = x(j+1)-x(j))
C
C  y' = (y(j+1)-y(j))/dx + (3A^2-1)dx y''(j)/6 + (3B^2-1)dx y''(j+1)/6
C
C and the second derivative is
C
C  y'' = A y''(j) + B y''(j+1)
C
C Input:
c  x(n)=x values in ascending order.
c  y(n)=y values at x points.
c  n=number of incoming data points.
c  yp1=y' at x(1) or else > 1e30 (latter uses natural spline).
c  ypn=y' at x(n) or else > 1e30 (as above).
c  Note that use of a "natural" spline has little to recommend it.
c  work(n)=work space.
c Output:
c  y2(n)=spline fit array of y'' values.

      implicit none
      integer n
      real*8 x(n),y(n),y2(n),work(n),yp1,ypn
      integer i,k
      real*8 sig,p,qn,workn

      if(yp1 .gt. 1.d+30) then
c       lower boundary condition is either natural ...
        y2(1) = 0.d0
        work(1) = 0.d0
      else
c       or else to have a specified first derivative.
        y2(1) =  - 0.5d0
        work(1) = (3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

c     Decomposition loop of tridiagonal algorithm:
      do 100 i=2,n-1
        sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p = sig*y2(i-1) + 2.d0
        y2(i) = (sig-1.d0)/p
        work(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &   /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*work(i-1))/p
  100 continue

      if(ypn.gt.1.d+30) then
c       Set upper boundary condition to be natural ...
        qn = 0.d0
        workn = 0.d0
      else
c       Or else to have a specified first derivative
        qn = 0.5d0
        workn=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif

      y2(n) = (workn-qn*work(n-1))/(qn*y2(n-1)+1.d0)
c     Backsubstitution loop of tridiagonal algorithm:
      do 200 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+work(k)
  200 continue

      return
      end
