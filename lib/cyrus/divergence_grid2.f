      subroutine divergence_grid2(y,div,r,r0,h,n)
      implicit real*8(a-h,o-z)
c Calculates the divergence of a purely radial vector field
c using 5 pt formulae on a uniform radial mesh.
c Div(y)=dy/dr+2*y/r
c                                             Cyrus Feb 1992

c  y     = radial part of field whose divergence is to be calculated
c  div   = divergence of function
c  r     = radial mesh
c  r0    = parameter of radial mesh
c  h     = exponential step size in radial coordinate
c  n     = number of points
c The radial mesh is exponential i.e. r_i=r_0*(exp((i-1)*h)-1)
c                           where     r_0=r_N/(exp((N-1)*h)-1)

      parameter(two=2.d0)
      dimension y(*),div(*),r(*)

      call deriv_grid2(y,div,r,r0,h,n)

      do 10 i=1,n
   10   div(i)=div(i)+two*y(i)/r(i)

      return
      end
