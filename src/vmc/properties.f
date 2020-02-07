cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     additional properties 
c     Friedemann Schautz
c     
c
c     properties so far:
c     1   2   3   4      5      6  
c     <x> <y> <z> <x**2> <y**2> <z**2> 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prop_compute(coord)
      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      implicit real*8(a-h,o-z)

      include 'properties.h'
c     electron coordinates
      dimension coord(3,*)

      if(iprop.eq.0) return
      do 5 i=1,nprop
       vprop(i)=0.d0
 5    enddo

      do 10 i=1,nelec
       do 20 m=1,3      
        vprop(m)  = vprop(m)+coord(m,i)
        vprop(3+m)= vprop(3+m) + coord(m,i)**2
 20    enddo
 10   enddo
      end

c-----------------------------------------------------------------------
      subroutine prop_init(iflg)
      implicit real*8(a-h,o-z)
      include 'properties.h'

      if(iprop.eq.0) return

      do i=1,nprop
       vprop_sum(i)=0.d0
      enddo

C$ iflg = 0: init *cum, *cm2 as well
      if(iflg.gt.0) return
      do i=1,nprop
       vprop_cum(i)=0.d0
       vprop_cm2(i)=0.d0
      enddo
      end

c-----------------------------------------------------------------------
      subroutine prop_cum(w)
      implicit real*8(a-h,o-z)
      include 'properties.h'

      if(iprop.eq.0) return
      do i=1,nprop
       vprop_now = vprop_sum(i)/w
       vprop_cm2(i)=vprop_cm2(i)+ vprop_sum(i)*vprop_now
       vprop_cum(i)=vprop_cum(i)+ vprop_sum(i)
      enddo
      end

c-----------------------------------------------------------------------
      subroutine prop_avrg(wcum,iblk,pav,perr)
      implicit real*8(a-h,o-z)
      include 'properties.h'
      dimension pav(MAXPROP),perr(MAXPROP)

      err(x,x2)=dsqrt(abs(x2/wcum-(x/wcum)**2)/iblk)

      if(iprop.eq.0) return
      do i=1,nprop
       perr(i)=err(vprop_cum(i),vprop_cm2(i))
       pav(i)=vprop_cum(i)/wcum
      enddo
      end
c-----------------------------------------------------------------------
      subroutine prop_dump(iu)
      implicit real*8(a-h,o-z)
      include 'properties.h'
      if(iprop.eq.0) return
      write(iu) nprop
      write(iu) (vprop_cum(i),vprop_cm2(i),i=1,nprop)
      end
      subroutine prop_rstrt(iu)
      implicit real*8(a-h,o-z)
      include 'properties.h'
      if(iprop.eq.0) return
      read(iu) nprop
      read(iu) (vprop_cum(i),vprop_cm2(i),i=1,nprop)
      end
c-----------------------------------------------------------------------
      subroutine prop_fin(passes,iblk,efin,eerr)
      implicit real*8(a-h,o-z)
      include 'properties.h'

      if(iprop.eq.0) return
      write(6,'(''--- additional properties ---'')')
      ipropprt_sav=ipropprt
      ipropprt=-1
      call prop_prt(passes,iblk,6)
      ipropprt=ipropprt_sav

      end
c-----------------------------------------------------------------------
      subroutine prop_prt(w,iblk,iu)
      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      implicit real*8(a-h,o-z)

c compute averages and print then out
      include 'properties.h'
      dimension pav(MAXPROP),perr(MAXPROP)

      common /icount_prop/ icount

c ipropprt 0 no printout
c          1 each iteration full printout
c         >1 after ipropprt iterations reduced printout
c         -1 force printout

      if(iprop.eq.0.or.ipropprt.eq.0) return

      if(ipropprt.gt.0.and.icount.ne.ipropprt) then
        icount=icount+1
        return
      endif

      icount=1

      call prop_avrg(w,iblk,pav(1),perr(1))
      

      write(iu,10)
      write(iu,20) 'X  ',pav(1),perr(1),pav(1)/dble(nelec),perr(1)
     $     /dble(nelec)
      write(iu,20) 'Y  ',pav(2),perr(2),pav(2)/dble(nelec),perr(2)
     $     /dble(nelec)
      write(iu,20) 'Z  ',pav(3),perr(3),pav(3)/dble(nelec),perr(3)
     $     /dble(nelec)
      write(iu,20) 'XX ',pav(4),perr(4),pav(4)/dble(nelec),perr(4)
     $     /dble(nelec)
      write(iu,20) 'YY ',pav(5),perr(5),pav(5)/dble(nelec),perr(5)
     $     /dble(nelec)
      write(iu,20) 'ZZ ',pav(6),perr(6),pav(6)/dble(nelec),perr(6)
     $     /dble(nelec)

c....dipole
      write(iu,50) 'center of nuclear charge: ',cc_nuc
      write(iu,30)
      dipx=cc_nuc(1)*nelec*2.5417 - pav(1) *2.5417
      dipy=cc_nuc(2)*nelec*2.5417 - pav(2) *2.5417
      dipz=cc_nuc(3)*nelec*2.5417 - pav(3) *2.5417
      dip=dsqrt(dipx**2+dipy**2+dipz**2)
      diperr=dabs (perr(1)*2.5417 * dipx / dip) +
     $       dabs (perr(2)*2.5417 * dipy / dip) + 
     $       dabs (perr(3)*2.5417 * dipz / dip) 
      write(iu,40) 'Dip X ',dipx,perr(1)*2.5417
      write(iu,40) 'Dip Y ',dipy,perr(2)*2.5417
      write(iu,40) 'Dip Z ',dipz,perr(3)*2.5417
      write(iu,40) 'Dip   ',dip,diperr



 10   format('-------- property operator averages  ----------')
 20   format(a3,'   ',f16.8,' +- ',f16.8,' ( ',f16.8,' +- ',f16.8,' )')
 30   format('-------- dipole operator averages  ----------')
 40   format(a6,'   ',f16.8,' +- ',f16.8)
 50   format(a25,' ',3f12.8)

      end

!*********************************************************************
        subroutine prop_cc_nuc(znuc,cent,iwctype,mctype,mcent, 
     &  ncent,cc_nuc)
!*********************************************************************

        implicit none

        integer  mctype,mcent,ncent
        integer  iwctype(mcent)
        real*8 znuc(mctype),cent(3,mcent)
        real*8 cc_nuc(3), tmp
        integer i,id

        cc_nuc(1)=0.d0
        cc_nuc(2)=0.d0
        cc_nuc(3)=0.d0
        tmp=0.d0
        do i=1,ncent
          id=znuc(iwctype(i))
          if(znuc(iwctype(i)).gt.2) id=znuc(iwctype(i))+2
          cc_nuc(1)=cc_nuc(1)+znuc(iwctype(i))*cent(1,i)
          cc_nuc(2)=cc_nuc(2)+znuc(iwctype(i))*cent(2,i)
          cc_nuc(3)=cc_nuc(3)+znuc(iwctype(i))*cent(3,i)
          tmp=tmp+znuc(iwctype(i))
        enddo
        cc_nuc(1)=cc_nuc(1)/tmp
        cc_nuc(2)=cc_nuc(2)/tmp
        cc_nuc(3)=cc_nuc(3)/tmp
!        write (*,*) 'Center of nuclear charge:', cc_nuc(:),tmp


        return
        end
