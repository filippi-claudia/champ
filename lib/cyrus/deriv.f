      subroutine deriv(y,der,eps,h,n)
      implicit real*8(a-h,o-z)
c Calculates the derivative of y using 5 pt on uniform mesh.
c The formulae can be found in Abramowitz and Stegun pg 914
c                                             Cyrus Feb 1992

c  y     = function whose derivative is to be calculated
c  der   = derivative of function
c  eps   = initial radial coordinate
c  h     = step size in radial coordinate
c  n     = number of points
c The radial mesh is uniform i.e.  r(i) = eps+(i-1)*h
      parameter(
     1 c00=1.d0,  c01=-8.d0,   c02=0.d0,   c03=8.d0,   c04=-1.d0,
     1 c10=-3.d0, c11=-10.d0,  c12=18.d0,  c13=-6.d0,  c14=1.d0,
     1 c20=-25.d0,c21=48.d0,   c22=-36.d0, c23=16.d0,  c24=-3.d0)
      dimension y(*),der(*)

      hinv12=1.d0/(h*12.d0)

      do 10 i=1,n
        r=eps+(i-1)*h
        if(i.eq.1) then
          der(i)=  c20*y(i)+c21*y(i+1)+c22*y(i+2)+c23*y(i+3)+c24*y(i+4)
         elseif(i.eq.2) then
          der(i)=  c10*y(i-1)+c11*y(i)+c12*y(i+1)+c13*y(i+2)+c14*y(i+3)
         elseif(i.eq.n-1) then
          der(i)=-(c10*y(i+1)+c11*y(i)+c12*y(i-1)+c13*y(i-2)+c14*y(i-3))
         elseif(i.eq.n) then
          der(i)=-(c20*y(i)+c21*y(i-1)+c22*y(i-2)+c23*y(i-3)+c24*y(i-4))
         else
          der(i)=  c00*y(i-2)+c01*y(i-1)+c02*y(i)+c03*y(i+1)+c04*y(i+2)
        endif
   10   der(i)=hinv12*der(i)
      return
      end
