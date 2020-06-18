      subroutine mmpol_prt(iblk,wgcum,wgcm2)
      implicit real*8(a-h,o-z)
 
      include 'dmc.h'
      include 'force.h'
      include 'mmpol.h'
      data hatokc/627.509541d0/

      dimension wgcum(MFORCE),wgcm2(MFORCE)

      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar

      rn_eff(w,w2)=w**2/w2
      error(x,x2,w,w2)=dsqrt(max((x2/w-(x/w)**2)/(rn_eff(w,w2)-1),0.d0))
      errg(x,x2,i)=error(x,x2,wgcum(i),wgcm2(i))

      if(immpol.eq.0.or.immpolprt.eq.0) return

      dmmpol_ave=dmmpol_cum/wgcum(1)
      cmmpol_ave=cmmpol_cum/wgcum(1)

      if(iblk.eq.1) then
        dmmpol_err=0
        cmmpol_err=0
       else
        dmmpol_err=errg(dmmpol_cum,dmmpol_cm2,1)
        cmmpol_err=errg(cmmpol_cum,cmmpol_cm2,1)

        idmmpol_err=nint(100000*dmmpol_err)
        icmmpol_err=nint(100000*cmmpol_err)
      endif

      evalg_eff=nconf*nstep*rn_eff(wgcum(1),wgcm2(1))
      rtevalg_eff1=dsqrt(evalg_eff-1)

      dmmpol_kcal=dmmpol_ave*hatokc
      cmmpol_kcal=cmmpol_ave*hatokc
      dekcal=dmmpol_err*hatokc
      cekcal=cmmpol_err*hatokc
      dckcal=dmmpol_kcal+cmmpol_kcal

      write(6,*)'    <H(QM/MM)/dipoles/charges/tot +- err (kcal/mol) '
      write(6,1000) dmmpol_kcal,dekcal,cmmpol_kcal,cekcal,dckcal
      write(6,*)'    <H(QM/MM)/dipoles/charges/tot +- err (hartree) '
      write(6,1000) dmmpol_ave,dmmpol_err,cmmpol_ave,cmmpol_err,dckcal/hatokc
      write(6,*)
 1000 format (f15.7,' +- ',2f15.7,' +- ',f15.7,' --> ',f15.7)

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_fin(iblk,wgcum,wgcm2)

      implicit real*8(a-h,o-z)

      include 'dmc.h'
      include 'force.h'
      include 'mmpol.h'

      dimension wgcum(MFORCE),wgcm2(MFORCE)

      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar

      if(immpol.eq.0) return
    
      immpolprt_sav=immpolprt
      immpolprt=1
      call mmpol_prt(iblk,wgcum(1),wgcm2(1))
      immpolprt=immpolprt_sav

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_save(iw)
      implicit real*8(a-h,o-z)
 
      include 'dmc.h'
      include 'mmpol.h'
      common /mmpol_hpsi/QMdp,QMq,eek_pol(3,MCHMM)
      common /mmpolo/ dmmpolo(MWALK),cmmpolo(MWALK),
     &         eeko1(MWALK,MCHMM),eeko2(MWALK,MCHMM),eeko3(MWALK,MCHMM)

      if(immpol.eq.0) return

      dmmpolo(iw)=QMdp
      cmmpolo(iw)=QMq

      do i=1,nchmm
        eeko1(iw,i)=eek_pol(1,i)
        eeko2(iw,i)=eek_pol(2,i)
        eeko3(iw,i)=eek_pol(3,i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_sum(p,q,iw)

      use mmpol_averages, only: cmmpol_cum, cmmpol_cm2, eek2_cum, dmmpol_sum, eek1_cm2, eek_sum, eek2_cm2
      use mmpol_averages, only: cmmpol_sum, dmmpol_cum, dmmpol_cm2, eek3_cum, eek1_cum, eek3_cm2

      implicit real*8(a-h,o-z)
 
      include 'dmc.h'
      include 'mmpol.h'
      common /mmpol_hpsi/QMdp,QMq,eek_pol(3,MCHMM)
      common /mmpolo/ dmmpolo(MWALK),cmmpolo(MWALK),
     &         eeko1(MWALK,MCHMM),eeko2(MWALK,MCHMM),eeko3(MWALK,MCHMM)

      if(immpol.eq.0) return

      dmmpol_sum=dmmpol_sum+p*QMdp+q*dmmpolo(iw)
      cmmpol_sum=cmmpol_sum+p*QMq+q*cmmpolo(iw)

      do i=1,nchmm
        eek_sum(1,i)= eek_sum(1,i)+p*eek_pol(1,i)+q*eeko1(iw,i)
        eek_sum(2,i)= eek_sum(2,i)+p*eek_pol(2,i)+q*eeko2(iw,i)
        eek_sum(3,i)= eek_sum(3,i)+p*eek_pol(3,i)+q*eeko3(iw,i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_cum(wsum)

      use mmpol_averages, only: cmmpol_cum, cmmpol_cm2, eek2_cum, dmmpol_sum, eek1_cm2, eek_sum, eek2_cm2
      use mmpol_averages, only: cmmpol_sum, dmmpol_cum, dmmpol_cm2, eek3_cum, eek1_cum, eek3_cm2

      implicit real*8(a-h,o-z)
 
      include 'dmc.h'
      include 'mmpol.h'

      if(immpol.eq.0) return

      dmmpolnow=dmmpol_sum/wsum
      cmmpolnow=cmmpol_sum/wsum

      dmmpol_cm2=dmmpol_cm2+dmmpol_sum*dmmpolnow
      cmmpol_cm2=cmmpol_cm2+cmmpol_sum*cmmpolnow

      dmmpol_cum=dmmpol_cum+dmmpol_sum
      cmmpol_cum=cmmpol_cum+cmmpol_sum

      do i=1,nchmm
        eek_now1=eek_sum(1,i)/wsum
        eek1_cm2(i)=eek1_cm2(i)+eek_sum(1,i)*eek_now1
        eek1_cum(i)=eek1_cum(i)+eek_sum(1,i)

        eek_now2=eek_sum(2,i)/wsum
        eek2_cm2(i)=eek2_cm2(i)+eek_sum(2,i)*eek_now2
        eek2_cum(i)=eek2_cum(i)+eek_sum(2,i)

        eek_now3=eek_sum(3,i)/wsum
        eek3_cm2(i)=eek3_cm2(i)+eek_sum(3,i)*eek_now3
        eek3_cum(i)=eek3_cum(i)+eek_sum(3,i)

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_splitj(iw,iw2)
      implicit real*8(a-h,o-z)

      include 'dmc.h'
      include 'mmpol.h'
      common /mmpolo/ dmmpolo(MWALK),cmmpolo(MWALK),
     &         eeko1(MWALK,MCHMM),eeko2(MWALK,MCHMM),eeko3(MWALK,MCHMM)


      dmmpolo(iw2)=dmmpolo(iw)
      cmmpolo(iw2)=cmmpolo(iw)

      do i=1,nchmm
        eeko1(iw2,i)=eeko1(iw,i)
        eeko2(iw2,i)=eeko2(iw,i)
        eeko3(iw2,i)=eeko3(iw,i)
      enddo

      return
      end
