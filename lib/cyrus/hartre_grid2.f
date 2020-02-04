      subroutine hartre_grid2(y,vh,r,r0,h,n)
      implicit real*8(a-h,o-z)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c:: Calculates the Hartree potential for density rho.
c::  y     = function whose Laplacian is to be calculated
c::  vh    = Hartree potential
c::  r     = radial mesh
c::  r0    = parameter of radial mesh
c::  h     = exponential step size in radial coordinate
c::  n     = number of points
c:: The radial mesh is exponential i.e. r_i=r_0*(exp((i-1)*h)-1)
c::                           where     r_0=r_N/(exp((N-1)*h)-1)
c::
c:: Simpson integration routine wts 1/3, 4/3, 2/3, 4/3 ... 4/3, 1/3
c:: Uses Newton's  4 pt integ. for the end pts if n is even
c:: Wts for this are 3/8, 9/8, 9/8, 3/8 
c:: Hence if n is even correction for wts is:                         
c:: c3=(1/3+3/8-2/3), c2=(9/8-4/3), c1=(9/8-2/3), c0=(3/8-4/3)        
c:: ( =1/8=.1125,       =-5/8=-.625   =11/8=1.375   =-23/8=-2.875) /3 
c:: The integral from pt.1 to pt. 2 obtained by fitting a quadratic to
c:: the 1st 3 pts is (5*y(1)+8*y(2)-y(3))*h/12
c:: Error is (y'''')*(h**4)/90                        if n is odd     
c::   "   "          "        +(y'''')*(h**5)*(3/80)  if n is even    
c::                                            Cyrus  Feb92           
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      parameter(c3=1.d0/8.d0,c2=-5.d0/8.d0,c1=11.d0/8.d0,c0=-23.d0/8.d0
     &,d1=1.25d0,d2=2.d0,d3=-.25d0
     &,zero=0.d0,two=2.d0,thrd=1.d0/3.d0,half=.5d0
     &,fourpi=12.5663706143592d0)

      dimension y(*),vh(*),r(*)

      y1(i)=y(i)*r(i)   *(r0+r(i))
      y2(i)=y(i)*r(i)**2*(r0+r(i))

      h3=h*thrd
      n1=n-1

c Calculate integral from 0 to r, (4*pi/r) Integ[dr' r'^2 rho(r')]
      vh(1)=zero
      odd=-half*y2(1)
      eve=zero
      do 10 i=1,n1,2
        odd=odd+y2(i)
        eve=eve+y2(i+1)
        s=two*(odd+two*eve)
        if(i+2.le.n) vh(i+2)=(s+y2(i+2)) / r(i+2)
        if(i.eq.1) then
          vh(2)=(d1*y2(1)+d2*y2(2)+d3*y2(3)) / r(i+1)
         else
          vh(i+1)=(s+c3*y2(i-2)+c2*y2(i-1)+c1*y2(i)+c0*y2(i+1)) / r(i+1)
        endif
   10   continue

c Add in contribution from r to final pt., (4*pi) Integ[dr' r' rho(r')]
c Do this by subtracting
c contrib. from 0 to r, and adding contrib. from 0 to final pt.
      odd=-half*y1(1)
      eve=zero
      do 20 i=1,n1,2
        odd=odd+y1(i)
        eve=eve+y1(i+1)
        s=two*(odd+two*eve)
        if(i+2.le.n) vh(i+2)=vh(i+2)-(s+y1(i+2))
        if(i.eq.1) then
          vh(2)=vh(2)-(d1*y1(1)+d2*y1(2)+d3*y1(3))
         else
          vh(i+1)=vh(i+1)-(s+c3*y1(i-2)+c2*y1(i-1)+c1*y1(i)+c0*y1(i+1))
        endif
   20   continue

c Calcul integ. from 0 to final pt.
      if(mod(n,2).eq.1) then
        vh1=s+y1(n)
       else
        vh1=s+c3*y1(n-3)+c2*y1(n-2)+c1*y1(n-1)+c0*y1(n)
      endif

      do 30 i=1,n
   30   vh(i)=fourpi*(vh(i)+vh1)*h3
      return
      end
