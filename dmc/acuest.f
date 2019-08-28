      subroutine acuest
c Written by Cyrus Umrigar, modified by Claudia Filippi
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'
      parameter (zero=0.d0,one=1.d0)

c routine to accumulate estimators for energy etc.
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /contrldmc/ tau,rttau,taueff(MFORCE),tautot,nfprod,idmc,ipq
     &,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
      common /iterat/ ipass,iblk
      common /estsum/ wsum,w_acc_sum,wfsum,wgsum(MFORCE),wg_acc_sum,wdsum,
     &wgdsum, wsum1(MFORCE),w_acc_sum1,wfsum1,wgsum1(MFORCE),wg_acc_sum1,
     &wdsum1, esum,efsum,egsum(MFORCE),esum1(MFORCE),efsum1,egsum1(MFORCE),
     &ei1sum,ei2sum,ei3sum, pesum(MFORCE),tpbsum(MFORCE),tjfsum(MFORCE),r2sum,
     &risum,tausum(MFORCE)
      common /estcum/ wcum,w_acc_cum,wfcum,wgcum(MFORCE),wg_acc_cum,wdcum,
     &wgdcum, wcum1,w_acc_cum1,wfcum1,wgcum1(MFORCE),wg_acc_cum1,
     &wdcum1, ecum,efcum,egcum(MFORCE),ecum1,efcum1,egcum1(MFORCE),
     &ei1cum,ei2cum,ei3cum, pecum(MFORCE),tpbcum(MFORCE),tjfcum(MFORCE),r2cum,
     &ricum,taucum(MFORCE)
      common /estcm2/ wcm2,wfcm2,wgcm2(MFORCE),wdcm2,wgdcm2, wcm21,
     &wfcm21,wgcm21(MFORCE),wdcm21, ecm2,efcm2,egcm2(MFORCE), ecm21,
     &efcm21,egcm21(MFORCE),ei1cm2,ei2cm2,ei3cm2, pecm2(MFORCE),tpbcm2(MFORCE),
     &tjfcm2(MFORCE),r2cm2,ricm2
      common /derivest/ derivsum(10,MFORCE),derivcum(10,MFORCE),derivcm2(MFORCE),
     &derivtotave_num_old(MFORCE)
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk
      common /age/ iage(MWALK),ioldest,ioldestmx
      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /forcest/ fgcum(MFORCE),fgcm2(MFORCE)
      common /force_dmc/ itausec,nwprod

c statement function for error calculation
      rn_eff(w,w2)=w**2/w2
      error(x,x2,w,w2)=dsqrt(max((x2/w-(x/w)**2)/(rn_eff(w,w2)-1),0.d0))
      errg(x,x2,i)=error(x,x2,wgcum(i),wgcm2(i))

c wt   = weight of configurations
c xsum = sum of values of x from dmc
c xnow = average of values of x from dmc
c xcum = accumulated sums of xnow
c xcm2 = accumulated sums of xnow**2
c xave = current average value of x
c xerr = current error of x

      iblk=iblk+1
      npass=iblk*nstep

      wnow=wsum/nstep
      wfnow=wfsum/nstep
      enow=esum/wsum
      efnow=efsum/wfsum
      ei1now=wfsum/wdsum
      ei2now=wgsum(1)/wgdsum
      rinow=risum/wgsum(1)
      r2now=r2sum/wgsum(1)

      wcm2=wcm2+wsum**2
      wfcm2=wfcm2+wfsum**2
      ecm2=ecm2+esum*enow
      efcm2=efcm2+efsum*efnow
      ei1cm2=ei1cm2+ei1now**2
      ei2cm2=ei2cm2+ei2now**2
      r2cm2=r2cm2+r2sum*r2now
      ricm2=ricm2+risum*rinow

      wcum=wcum+wsum
      wfcum=wfcum+wfsum
      wdcum=wdcum+wdsum
      wgdcum=wgdcum+wgdsum
      ecum=ecum+esum
      efcum=efcum+efsum
      ei1cum=ei1cum+ei1now
      ei2cum=ei2cum+ei2now
      r2cum=r2cum+r2sum
      ricum=ricum+risum

      do 10 ifr=1,nforce

        wgnow=wgsum(ifr)/nstep
        egnow=egsum(ifr)/wgsum(ifr)
        penow=pesum(ifr)/wgsum(ifr)
        tpbnow=tpbsum(ifr)/wgsum(ifr)
        tjfnow=tjfsum(ifr)/wgsum(ifr)

        wgcm2(ifr)=wgcm2(ifr)+wgsum(ifr)**2
        egcm2(ifr)=egcm2(ifr)+egsum(ifr)*egnow
        pecm2(ifr)=pecm2(ifr)+pesum(ifr)*penow
        tpbcm2(ifr)=tpbcm2(ifr)+tpbsum(ifr)*tpbnow
        tjfcm2(ifr)=tjfcm2(ifr)+tjfsum(ifr)*tjfnow

        wgcum(ifr)=wgcum(ifr)+wgsum(ifr)
        egcum(ifr)=egcum(ifr)+egsum(ifr)
        pecum(ifr)=pecum(ifr)+pesum(ifr)
        tpbcum(ifr)=tpbcum(ifr)+tpbsum(ifr)
        tjfcum(ifr)=tjfcum(ifr)+tjfsum(ifr)

        do 5 k=1,10
  5        derivcum(k,ifr)=derivcum(k,ifr)+derivsum(k,ifr)

        if(iblk.eq.1) then
          egerr=0
          peerr=0
          tpberr=0
          tjferr=0
         else
          egerr=errg(egcum(ifr),egcm2(ifr),ifr)
          peerr=errg(pecum(ifr),pecm2(ifr),ifr)
          tpberr=errg(tpbcum(ifr),tpbcm2(ifr),ifr)
          tjferr=errg(tjfcum(ifr),tjfcm2(ifr),ifr)
        endif

        egave=egcum(ifr)/wgcum(ifr)
        peave=pecum(ifr)/wgcum(ifr)
        tpbave=tpbcum(ifr)/wgcum(ifr)
        tjfave=tjfcum(ifr)/wgcum(ifr)

        if(ifr.gt.1) then
          fgcum(ifr)=fgcum(ifr)+wgsum(1)*(egnow-egsum(1)/wgsum(1))
          fgcm2(ifr)=fgcm2(ifr)+wgsum(1)*(egnow-egsum(1)/wgsum(1))**2
          fgave=egcum(1)/wgcum(1)-egcum(ifr)/wgcum(ifr)
          fgave=fgave/deltot(ifr)
          if(iblk.eq.1) then
            fgerr=0
            ifgerr=0
           else
            fgerr=errg(fgcum(ifr),fgcm2(ifr),1)
            fgerr=fgerr/abs(deltot(ifr))
            ifgerr=nint(1.0d9* fgerr)
          endif
          egave1=egcum(1)/wgcum(1)
          if(iblk.eq.1) derivtotave_num_old(ifr)=0.d0
          derivtotave_num=-(derivcum(1,ifr)-derivcum(1,1)+derivcum(2,ifr)-derivcum(2,1)-egave1*(derivcum(3,ifr)-derivcum(3,1)))
          derivtotave=derivtotave_num/wgcum(1)
          delta_derivtotave_num=derivtotave_num-derivtotave_num_old(ifr)
          derivtotave_num_old(ifr)=derivtotave_num
          derivcm2(ifr)=derivcm2(ifr)+delta_derivtotave_num**2/wgsum(1)
          if(iblk.eq.1) then
            derivgerr=0
            iderivgerr=0
           else
            derivgerr=errg(derivtotave,derivcm2(ifr),1)
            derivgerr=derivgerr/abs(deltot(ifr))
            iderivgerr=nint(1.0d9* derivgerr)
          endif
         else
          call optjas_cum(wgsum(1),egnow)
          call optorb_cum(wgsum(1),egsum(1))
          call optci_cum(wgsum(1))

          call prop_cum(wgsum(ifr))
          call pcm_cum(wgsum(ifr))
          call mmpol_cum(wgsum(ifr))
        endif

        taucum(ifr)=taucum(ifr)+tausum(ifr)

c write out header first time

        if (iblk.eq.1.and.ifr.eq.1)
     &  write(6,'(t5,''egnow'',t15,''egave'',t21,''(egerr)'' ,t32
     &  ,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr)'',t66
     &  ,''tjfave'',t72,''(tjferr)'',t83,''fgave'',t97,''(fgerr)'',
     &  t109,''fgave_n'',''(fgerr_n)'',t130,''npass'',t140,''wgsum'',t150,''ioldest'')')

c write out current values of averages etc.

        iegerr=nint(100000* egerr)
        ipeerr=nint(100000* peerr)
        itpber=nint(100000*tpberr)
        itjfer=nint(100000*tjferr)

        if(ifr.eq.1) then
          write(6,'(f10.5,4(f10.5,''('',i5,'')''),46x,3i10)') egnow,
     &    egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,npass,
     &    nint(wgsum(ifr)),ioldest

          call prop_prt_dmc(iblk,0,wgcum,wgcm2)
          call pcm_prt(iblk,wgcum,wgcm2)
          call mmpol_prt(iblk,wgcum,wgcm2)
         else
          write(6,'(f10.5,4(f10.5,''('',i5,'')''),f14.9,
     &    ''('',i9,'')'',f14.9,''('',i9,'')'',6x,i10)') egnow,
     &    egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,
     &    fgave,ifgerr,derivtotave,iderivgerr,nint(wgsum(ifr))
        endif
   10 continue

c zero out xsum variables for metrop

      wsum=zero
      wfsum=zero
      wdsum=zero
      wgdsum=zero
      esum=zero
      efsum=zero
      ei1sum=zero
      ei2sum=zero
      r2sum=zero
      risum=zero

      do 20 ifr=1,nforce
        egsum(ifr)=zero
        wgsum(ifr)=zero
        pesum(ifr)=zero
        tpbsum(ifr)=zero
        tjfsum(ifr)=zero
        tausum(ifr)=zero
        do 20 k=1,10
   20     derivsum(k,ifr)=zero

      call optorb_init(1)
      call optci_init(1)

      call prop_init(1)
      call pcm_init(1)
      call mmpol_init(1)

      return
      end
