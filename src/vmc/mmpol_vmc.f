      subroutine mmpol_prt(wcum,iblk)

      use mmpol_cntrl, only: icall_mm, ich_mmpol, immpol, immpolprt, isites_mmpol
      use mmpol_averages, only: cmmpol_cm2, cmmpol_cum, cmmpol_sum, dmmpol_cm2, dmmpol_cum, dmmpol_sum
      use mmpol_averages, only: eek1_cm2, eek1_cum, eek2_cm2, eek2_cum, eek3_cm2, eek3_cum, eek_sum

      implicit real*8(a-h,o-z)

      include 'mmpol.h'

      err(x,x2)=dsqrt(abs(x2/wcum-(x/wcum)**2)/iblk)

      if(immpol.eq.0.or.immpolprt.eq.0) return

      dmmpol_ave=dmmpol_cum/wcum
      cmmpol_ave=cmmpol_cum/wcum

      if(iblk.eq.1) then
        dmmpol_err=0
        cmmpol_err=0

       else
        dmmpol_err=err(dmmpol_cum,dmmpol_cm2)
        cmmpol_err=err(cmmpol_cum,cmmpol_cm2)

        idmmpol_err=nint(100000*dmmpol_err)
        icmmpol_err=nint(100000*cmmpol_err)
      endif

      rtpass=dsqrt(wcum)
      if(immpol.eq.2)then
        write(6,'(''mmpol dE (induced dipoles) ='',t17,f12.7,'' +-'',f11.7,f9.5)') dmmpol_ave,dmmpol_err,dmmpol_err*rtpass
        write(6,'(''mmpol dE  (fixed charges)  ='',t17,f12.7,'' +-'',f11.7,f9.5)') cmmpol_ave,cmmpol_err,cmmpol_merr*rtpass
        write(6,'(''<H(QM/MM)>'',t17,f12.7,'' +-'',f11.7,f9.5)') dmmpol_ave+cmmpol_ave
      endif
      if(immpol.eq.3)then
        qmmpol_ave=dmmpol_ave+cmmpol_ave
        qmmpol_err=dsqrt(dmmpol_err**2.0d0+cmmpol_err**2.0d0)
        write(6,'(''mmpol <H(QM/MM)> ='',t17,f12.7,'' +-'',f11.7,f9.5)') qmmpol_ave,qmmpol_err,qmmpol_err*rtpass
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_fin(wcum,iblk)

      use contrl, only: idump, irstar, isite, n_conf, nblk, nblkeq, nconf_new, nstep
      use mmpol_cntrl, only: icall_mm, ich_mmpol, immpol, immpolprt, isites_mmpol
      use mmpol_parms, only: chmm, nchmm, rqq, x_mmpol
      use mmpol_averages, only: cmmpol_cm2, cmmpol_cum, cmmpol_sum, dmmpol_cm2, dmmpol_cum, dmmpol_sum
      use mmpol_averages, only: eek1_cm2, eek1_cum, eek2_cm2, eek2_cum, eek3_cm2, eek3_cum, eek_sum

      implicit real*8(a-h,o-z)

      include 'mmpol.h'
    

      data hatokc/627.509541d0/

      dimension eek1_ave(MCHMM),eek1_err(MCHMM)
      dimension eek2_ave(MCHMM),eek2_err(MCHMM)
      dimension eek3_ave(MCHMM),eek3_err(MCHMM)
      dimension eek_ave(3,MCHMM),eek_err(3,MCHMM)

      err(x,x2)=dsqrt(abs(x2/wcum-(x/wcum)**2)/iblk)

      if(immpol.eq.0) return

      rtpass=dsqrt(wcum)

      dmmpol_ave=dmmpol_cum/wcum
      cmmpol_ave=cmmpol_cum/wcum

      dmmpol_err=err(dmmpol_cum,dmmpol_cm2)
      cmmpol_err=err(cmmpol_cum,cmmpol_cm2)

      idmmpol_err=nint(100000*dmmpol_err)
      icmmpol_err=nint(100000*cmmpol_err)

      if(immpol.eq.2)then
        write(6,*)
        write(6,'(''mmpol dE (induced dipoles) ='',t17,f12.7,'' +-'',f11.7,f9.5)') dmmpol_ave,dmmpol_err,dmmpol_err*rtpass
        write(6,'(''mmpol dE  (fixed charges)  ='',t17,f12.7,'' +-'',f11.7,f9.5)') cmmpol_ave,cmmpol_err,cmmpol_merr*rtpass
      endif

      if(immpol.eq.3)then
        qmmpol_ave=dmmpol_ave+cmmpol_ave
        qmmpol_err=dsqrt(dmmpol_err**2.0d0+cmmpol_err**2.0d0)
        write(6,'(''mmpol <H(QM/MM)> ='',t17,f12.7,'' +-'',f11.7,f9.5)') qmmpol_ave,qmmpol_err,qmmpol_err*rtpass
      endif

      if(immpol.ne.3)then
        do  i=1,nchmm
          eek1_ave(i)=eek1_cum(i)/wcum
          eek1_err(i)=err(eek1_cum(i),eek1_cm2(i))
          eek_ave(1,i)=eek1_cum(i)/wcum
          eek_err(1,i)=err(eek1_cum(i),eek1_cm2(i))

          eek2_ave(i)=eek2_cum(i)/wcum
          eek2_err(i)=err(eek2_cum(i),eek2_cm2(i))
          eek_ave(2,i)=eek2_cum(i)/wcum
          eek_err(2,i)=err(eek2_cum(i),eek2_cm2(i))

          eek3_ave(i)=eek3_cum(i)/wcum
          eek3_err(i)=err(eek3_cum(i),eek3_cm2(i))
          eek_ave(3,i)=eek3_cum(i)/wcum
          eek_err(3,i)=err(eek3_cum(i),eek3_cm2(i))
        enddo
 	call mmpol_dipoles(eek_ave,eek_err)

	dmmpol_kcal=dmmpol_ave*hatokc
	cmmpol_kcal=cmmpol_ave*hatokc
	dekcal=dmmpol_err*hatokc
	cekcal=cmmpol_err*hatokc
	dckcal=dmmpol_kcal+cmmpol_kcal
	write(6,*)'    <H(QM/MM)/dipoles/charges/tot +- err (kcal/mol) '
	write(6,1000)dmmpol_kcal,dekcal,cmmpol_kcal,cekcal,dckcal
	write(6,*)'    <H(QM/MM)/dipoles/charges/tot +- err (hartree) '
	write(6,1000)dmmpol_ave,dmmpol_err,cmmpol_ave,cmmpol_err,dckcal/hatokc
        write(6,*)
      endif


 1000 format(9F12.5)

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_save
      use mmpol_hpsi, only: eek_pol, peQMdp, peQMq
      use mmpolo, only: cmmpolo, dmmpolo, eeko
      use mmpol_cntrl, only: icall_mm, ich_mmpol, immpol, immpolprt, isites_mmpol
      use mmpol_parms, only: chmm, nchmm, rqq, x_mmpol
      implicit real*8(a-h,o-z)





      include 'mmpol.h'

      if(immpol.eq.0) return
      dmmpolo=QMdp
      cmmpolo=QMq

      do i=1,nchmm
        eeko(1,i)=eek_pol(1,i)
        eeko(2,i)=eek_pol(2,i)
        eeko(3,i)=eek_pol(3,i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_sum(p,q)

      use mmpol_hpsi, only: eek_pol, peQMdp, peQMq
      use mmpolo, only: cmmpolo, dmmpolo, eeko
      use mmpol_cntrl, only: icall_mm, ich_mmpol, immpol, immpolprt, isites_mmpol
      use mmpol_parms, only: chmm, nchmm, rqq, x_mmpol
      use mmpol_averages, only: cmmpol_cm2, cmmpol_cum, cmmpol_sum, dmmpol_cm2, dmmpol_cum, dmmpol_sum
      use mmpol_averages, only: eek1_cm2, eek1_cum, eek2_cm2, eek2_cum, eek3_cm2, eek3_cum, eek_sum

      implicit real*8(a-h,o-z)

      include 'mmpol.h'

      if(immpol.eq.0) return
      dmmpol_sum=dmmpol_sum+p*QMdp+q*dmmpolo
      cmmpol_sum=cmmpol_sum+p*QMq+q*cmmpolo

      do i=1,nchmm
        eek_sum(1,i)= eek_sum(1,i)+p*eek_pol(1,i)+q*eeko(1,i)
        eek_sum(2,i)= eek_sum(2,i)+p*eek_pol(2,i)+q*eeko(2,i)
        eek_sum(3,i)= eek_sum(3,i)+p*eek_pol(3,i)+q*eeko(3,i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_cum(wsum)

      use mmpol_cntrl, only: icall_mm, ich_mmpol, immpol, immpolprt, isites_mmpol
      use mmpol_parms, only: chmm, nchmm, rqq, x_mmpol
      use mmpol_averages, only: cmmpol_cm2, cmmpol_cum, cmmpol_sum, dmmpol_cm2, dmmpol_cum, dmmpol_sum
      use mmpol_averages, only: eek1_cm2, eek1_cum, eek2_cm2, eek2_cum, eek3_cm2, eek3_cum, eek_sum

      implicit real*8(a-h,o-z)

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
