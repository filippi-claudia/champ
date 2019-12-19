      subroutine laplac_grid2(y,del2,r,r0,h,n)
      implicit real*8(a-h,o-z)
c Calculates the Laplacian of y using 5 pt 1st and 2nd derivatives.
c The formulae can be found in Abramowitz and Stegun pg 914
c                                             Cyrus Feb 1992

c  y     = function whose Laplacian is to be calculated
c  del2  = Laplacian of function
c  r     = radial mesh
c  r0    = parameter of radial mesh
c  h     = exponential step size in radial coordinate
c  n     = number of points
c The radial mesh is exponential i.e. r_i=r_0*(exp((i-1)*h)-1)
c                           where     r_0=r_N/(exp((N-1)*h)-1)
c Laplac(y) = d2y_dr2 + (2/r)*dy_dr
c           = d2y_di2*di_dr**2 + dy_di*(d2i_dr2+(2/r)*di_dr)

      parameter(one=1.d0,two=2.d0,d1b12=1.d0/12.d0,
     1 c00=1.d0,  c01=-8.d0,   c02=0.d0,   c03=8.d0,   c04=-1.d0,
     1 c10=-3.d0, c11=-10.d0,  c12=18.d0,  c13=-6.d0,  c14=1.d0,
     1 c20=-25.d0,c21=48.d0,   c22=-36.d0, c23=16.d0,  c24=-3.d0,
     1 d00=-1.d0, d01=16.d0,   d02=-30.d0, d03=16.d0,  d04=-1.d0,
     1 d10=11.d0, d11=-20.d0,  d12=6.d0,   d13=4.d0,   d14=-1.d0,
     1 d20=35.d0, d21=-104.d0, d22=114.d0, d23=-56.d0, d24=11.d0)
      dimension y(*),del2(*),r(*)

      do 10 i=1,n
        if(i.eq.1) then
          d2= (d20*y(i)+d21*y(i+1)+d22*y(i+2)+d23*y(i+3)+d24*y(i+4))
          d1= (c20*y(i)+c21*y(i+1)+c22*y(i+2)+c23*y(i+3)+c24*y(i+4))
         elseif(i.eq.2) then
          d2= (d10*y(i-1)+d11*y(i)+d12*y(i+1)+d13*y(i+2)+d14*y(i+3))
          d1= (c10*y(i-1)+c11*y(i)+c12*y(i+1)+c13*y(i+2)+c14*y(i+3))
         elseif(i.eq.n-1) then
          d2= (d10*y(i+1)+d11*y(i)+d12*y(i-1)+d13*y(i-2)+d14*y(i-3))
          d1=-(c10*y(i+1)+c11*y(i)+c12*y(i-1)+c13*y(i-2)+c14*y(i-3))
         elseif(i.eq.n) then
          d2= (d20*y(i)+d21*y(i-1)+d22*y(i-2)+d23*y(i-3)+d24*y(i-4))
          d1=-(c20*y(i)+c21*y(i-1)+c22*y(i-2)+c23*y(i-3)+c24*y(i-4))
         else
          d2= (d00*y(i-2)+d01*y(i-1)+d02*y(i)+d03*y(i+1)+d04*y(i+2))
          d1= (c00*y(i-2)+c01*y(i-1)+c02*y(i)+c03*y(i+1)+c04*y(i+2))
        endif
        twobyr=two/r(i)
        di_dr=one/(h*(r0+r(i)))
        d2i_dr2=-di_dr/(r0+r(i))
   10   del2(i)=d1b12*(d2*di_dr**2+d1*(d2i_dr2+twobyr*di_dr))
      return
      end
