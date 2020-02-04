      subroutine intpol(x,y,n,xi,yi,ni,m)
c Does a Lagrange interpolation.
c x and y are arrays of length n with the original data.
c xi is the point at which the interpolated value yi is needed.
c ni is the closest mesh point below xi.
c m = order of polynomial using m+1 original grid pts
c If m < 1, then m is reset to 3
c Probably a good idea to have m odd, so there are an equal number of
c original mesh pts below and above xi.
c If m=3 the 4 points used to perform the interpolation are
c ni-1,ni,ni+1,ni+2 (2 below, 2 above) except at the end points.
c If m=4 the 5 points used to perform the interpolation are
c ni-2,ni-1,ni,ni+1,ni+2 (3 below, 2 above) except at the end points.

      implicit real*8(a-h,o-z)
      dimension x(n),y(n)
      data zero,one/0.d0,1.d0/

      if(m.lt.1) m=3
      n1=min0(n-m,max0(1,ni-m/2))
      n2=n1+m
c     write(6,'(i5,9d12.4)') n1,xi,(x(ii),ii=n1,n2)
      yi=zero
      do 20 i=n1,n2
        prod=one
        do 10 j=n1,n2
          if(j.ne.i) prod=prod*(xi-x(j))/(x(i)-x(j))
   10     continue
   20   yi=yi+prod*y(i)
      return
      end
