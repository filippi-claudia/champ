      subroutine jastrowe(iel,x,v,d2,value,iflag)
c Written by Claudia Filippi by modifying jastrow

      implicit real*8(a-h,o-z)

      include 'vmc.h'

      parameter (zero=0.d0)

      common /contr2/ ijas,icusp,icusp2,isc,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch

      include 'pseudo.h'

      include 'force.h'

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr

      dimension x(3,*),v(3,*)

      do 10 i=1,nelec
        v(1,i)=zero
        v(2,i)=zero
   10   v(3,i)=zero

      call jastrow4e(iel,x,v,d2,value,iflag)

      return
      end
