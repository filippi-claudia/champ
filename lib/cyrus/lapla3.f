      subroutine lapla3(y,del2,eps,h,n)
      implicit real*8(a-h,o-z)
c Calculates the Laplacian of y using 3 pt 1st and 2nd derivatives.
c Uses forward and backward formulae at end points and central
c differences for all the rest.
c                                             Cyrus Feb 1992

c  y     = function whose Laplacian is to be calculated
c  del2  = Laplacian of function
c  eps   = initial radial coordinate
c  h     = step size in radial coordinate
c  n     = number of points
c The radial mesh is uniform i.e.  r(i) = eps+(i-1)*h
      dimension y(*),del2(*)

      do 10 i=1,n
        r=eps+(i-1)*h
        if(i.eq.1) then
          del2(i)=((2*y(i)-5*y(i+1)+4*y(i+2)-y(i+3))/h
     1    + (-3*y(i)+4*y(i+1)-y(i+2))/r) /h
         elseif(i.eq.n) then
          del2(i)=((2*y(i)-5*y(i-1)+4*y(i-2)-y(i-3))/h
     1    + ( 3*y(i)-4*y(i-1)+y(i-2))/r) /h
         else
          del2(i)=((y(i-1)-2*y(i)+y(i+1))/h + (y(i+1)-y(i-1))/r) /h
        endif
   10   continue
      return
      end
