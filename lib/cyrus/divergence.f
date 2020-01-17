      subroutine divergence(y,div,eps,h,n)
      implicit real*8(a-h,o-z)
c Calculates the divergence of a purely radial vector field
c using 5 pt formulae on a uniform radial mesh.
c Div(y)=dy/dr+2*y/r
c                                             Cyrus Feb 1992

c  y     = radial part of field whose divergence is to be calculated
c  div   = divergence of function
c  eps   = initial radial coordinate
c  h     = step size in radial coordinate
c  n     = number of points
c The radial mesh is uniform i.e.  r(i) = eps+(i-1)*h
      parameter(two=2.d0)
      dimension y(*),div(*)

      call deriv(y,div,eps,h,n)

      do 10 i=1,n
        r=eps+(i-1)*h
   10   div(i)=div(i)+two*y(i)/r

      return
      end
