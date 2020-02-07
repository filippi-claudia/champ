      subroutine deriv_jastrow(x,v,d2,div_vj,value)
c Written by Claudia Filippi

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      implicit real*8(a-h,o-z)

      parameter (zero=0.d0)

      include 'vmc.h'
      include 'force.h'
      include 'pseudo.h'

      common /contr2/ ijas,icusp,icusp2,isc,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch


      dimension x(3,*),v(3,*),div_vj(MELEC)

      do 10 i=1,nelec
        v(1,i)=zero
        v(2,i)=zero
   10   v(3,i)=zero
      d2=zero

      call deriv_jastrow4(x,v,d2,value)

      return
      end
