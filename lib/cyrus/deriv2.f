      subroutine deriv2(y,d2,eps,h,n)
      implicit real*8(a-h,o-z)
c Calculates the 2nd deriv. of y using 5 pt formulae.
c The formulae can be found in Abramowitz and Stegun pg 914
c                                             Cyrus Feb 1992

c  y     = function whose 2nd deriv. is to be calculated
c  d2    = 2nd derivative of function
c  eps   = initial radial coordinate
c  h     = step size in radial coordinate
c  n     = number of points

c The radial mesh is uniform i.e.  r(i) = eps+(i-1)*h

c The interior pts have errors that converge as h^4 whereas
c the 1st, 2nd, nth and (n-1)th pts have errors that converge as h^3.
c The errors of the 1st and nth pts are an order of mag. larger than
c those of the 2nd and (n-1)th pts respectively.

      parameter (one=1.d0,twelve=12.d0,
     1 d00=-1.d0, d01=16.d0,   d02=-30.d0, d03=16.d0,  d04=-1.d0,
     1 d10=11.d0, d11=-20.d0,  d12=6.d0,   d13=4.d0,   d14=-1.d0,
     1 d20=35.d0, d21=-104.d0, d22=114.d0, d23=-56.d0, d24=11.d0)
      dimension y(*),d2(*)

      oneby12hsq=one/(twelve*h*h)

      do 10 i=1,n
        if(i.eq.1) then
          d2(i)= (d20*y(i)+d21*y(i+1)+d22*y(i+2)+d23*y(i+3)+d24*y(i+4))
         elseif(i.eq.2) then
          d2(i)= (d10*y(i-1)+d11*y(i)+d12*y(i+1)+d13*y(i+2)+d14*y(i+3))
         elseif(i.eq.n-1) then
          d2(i)= (d10*y(i+1)+d11*y(i)+d12*y(i-1)+d13*y(i-2)+d14*y(i-3))
         elseif(i.eq.n) then
          d2(i)= (d20*y(i)+d21*y(i-1)+d22*y(i-2)+d23*y(i-3)+d24*y(i-4))
         else
          d2(i)= (d00*y(i-2)+d01*y(i-1)+d02*y(i)+d03*y(i+1)+d04*y(i+2))
        endif
   10   d2(i)=d2(i)*oneby12hsq
      return
      end
