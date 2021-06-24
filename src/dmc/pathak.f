      subroutine pathak(distance,f,epsilon)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'

      parameter(one=1.d0)

      x1=distance/epsilon

      if(x1.gt.1.d0) then
        f=1.d0
      else
        x2=x1*x1
        x4=x2*x2
        x6=x2*x4
        f=7.d0*x6-15.d0*x4+9.d0*x2
      endif

      return
      end