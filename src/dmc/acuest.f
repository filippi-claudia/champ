      subroutine acuest
c MPI version created by Claudia Filippi starting from serial version
c routine to accumulate estimators for energy etc.

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'
      include 'mpif.h'
      parameter (zero=0.d0,one=1.d0)

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
      common /derivest/ derivsum(10,MFORCE),derivcum(10,MFORCE)
     &,derivcm2(MFORCE),derivtotave_num_old(MFORCE),derivcm(MFORCE)
      common /age/ iage(MWALK),ioldest,ioldestmx
      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /forcest/ fgcum(MFORCE),fgcm2(MFORCE)
      common /sigma_branch/ sigma

      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /derivanaly/ deriv_esum(10,3,MCENT,PTH),deriv_ecum(10,3,MCENT,PTH),
     &esnake(3,MCENT,MWALK,PTH),ehist(3,MCENT,MWALK,0:MFORCE_WT_PRD,PTH),
     &deriv_eold(3,MCENT,MWALK),pold(MWALK,PTH),deriv_cm(3,MCENT,PTH),
     &deriv_cm2(3,MCENT,PTH),eps_pathak(PTH),ipathak
      common /force_analy/ iforce_analy
      
      character*12 mode
      common /contr3/ mode

      common /mpiblk/ iblk_proc
      logical wid
      common /mpiconf/ idtask,nproc,wid

      dimension egcollect(MFORCE),wgcollect(MFORCE),pecollect(MFORCE),
     &tpbcollect(MFORCE),tjfcollect(MFORCE),eg2collect(MFORCE),wg2collect(MFORCE),
     &pe2collect(MFORCE),tpb2collect(MFORCE),tjf2collect(MFORCE),fsum(MFORCE),
     &f2sum(MFORCE),eg2sum(MFORCE),wg2sum(MFORCE),pe2sum(MFORCE),tpb2sum(MFORCE),
     &tjf2sum(MFORCE),taucollect(MFORCE),fcollect(MFORCE),f2collect(MFORCE),
     &derivcollect(10,MFORCE),derivsumave(MFORCE),derivsumave2(MFORCE),
     &derivsumcollect(MFORCE),derivsum2collect(MFORCE),deriv_ecol(10,3,MCENT,PTH),
     &deriv_ave(10,3,MCENT,PTH),deriv_avetot(3,MCENT,PTH),deriv_sumcol2(3,MCENT,PTH),
     &deriv_sumave(3,MCENT,PTH),deriv_sumave2(3,MCENT,PTH),deriv_sumcol(3,MCENT,PTH),
     &deriv_err(3,MCENT,PTH)

c statement function for error calculation
      rn_eff(w,w2)=w**2/w2
      error(x,x2,w,w2)=dsqrt(max((x2/w-(x/w)**2)/(rn_eff(w,w2)-1),0.d0))
      errg(x,x2,i)=error(x,x2,wgcum(i),wgcm2(i))

      if(mode.eq.'dmc_one_mpi2') then
        call acuest_gpop
        return
      endif

c wt   = weight of configurations
c xsum = sum of values of x from dmc
c xnow = average of values of x from dmc
c xcum = accumulated sums of xnow
c xcm2 = accumulated sums of xnow**2
c xave = current average value of x
c xerr = current error of x

      iblk=iblk+1
      iblk_proc=iblk_proc+nproc

      npass=iblk_proc*nstep

      wnow=wsum/nstep
      wfnow=wfsum/nstep
      enow=esum/wsum
      efnow=efsum/wfsum
      ei1now=wfsum/wdsum
      ei2now=wgsum(1)/wgdsum
      rinow=risum/wgsum(1)
      r2now=r2sum/wgsum(1)

      ei1cm2=ei1cm2+ei1now**2
      ei2cm2=ei2cm2+ei2now**2
      r2cm2=r2cm2+r2sum*r2now
      ricm2=ricm2+risum*rinow

      wdcum=wdcum+wdsum
      wgdcum=wgdcum+wgdsum
      ei1cum=ei1cum+ei1now
      ei2cum=ei2cum+ei2now
      r2cum=r2cum+r2sum
      ricum=ricum+risum

      w2sum=wsum**2
      wf2sum=wfsum**2
      e2sum=esum*enow
      ef2sum=efsum*efnow

      do 10 ifr=1,nforce
        wgnow=wgsum(ifr)/nstep
        egnow=egsum(ifr)/wgsum(ifr)
        penow=pesum(ifr)/wgsum(ifr)
        tpbnow=tpbsum(ifr)/wgsum(ifr)
        tjfnow=tjfsum(ifr)/wgsum(ifr)

        wg2sum(ifr)=wgsum(ifr)**2
        eg2sum(ifr)=egsum(ifr)*egnow
        pe2sum(ifr)=pesum(ifr)*penow
        tpb2sum(ifr)=tpbsum(ifr)*tpbnow
        tjf2sum(ifr)=tjfsum(ifr)*tjfnow
        if(ifr.gt.1) then
          fsum(ifr)=wgsum(1)*(egnow-egsum(1)/wgsum(1))
          f2sum(ifr)=wgsum(1)*(egnow-egsum(1)/wgsum(1))**2
          derivsumave(ifr)=-(derivsum(1,ifr)-derivsum(1,1)+derivsum(2,ifr)-derivsum(2,1)-enow*(derivsum(3,ifr)-derivsum(3,1)))
          derivsumave2(ifr)=(derivsumave(ifr))**2/wgsum(1)
        endif
  10  continue

      if(iforce_analy.eq.1) then
        if(ipathak.gt.1) then
c        open(120,file='force_sum2',form='formatted',status='unknown')
          do 101 iph=1,ipathak
            do 102 ic=1,ncent
              do 102 k=1,3
                deriv_sumave(k,ic,iph)=-(deriv_esum(1,k,ic,iph)+deriv_esum(2,k,ic,iph)
     &                                -egsum(1)/wgsum(1)*deriv_esum(3,k,ic,iph))
  102           deriv_sumave2(k,ic,iph)=(deriv_sumave(k,ic,iph))**2/wgsum(1)
  101       write(6,'(i,f10.4,4f18.10)') iblk,eps_pathak(iph),-deriv_esum(1,3,2,iph)/wgsum(1),-deriv_esum(2,3,2,iph)/wgsum(1),
     &             egsum(1)/wgsum(1)*deriv_esum(3,3,2,iph)/wgsum(1),deriv_sumave(3,2,iph)/wgsum(1)
        else
          do 103 ic=1,ncent
            do 103 k=1,3
              deriv_sumave(k,ic,1)=-(deriv_esum(1,k,ic,1)+deriv_esum(2,k,ic,1)
     &                             -egsum(1)/wgsum(1)*deriv_esum(3,k,ic,1))/wgsum(1)
  103         deriv_sumave2(k,ic,1)=(deriv_sumave(k,ic,1))**2/wgsum(1)
        endif
      endif


      call mpi_allreduce(wgsum,wgcollect,MFORCE
     &,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(egsum,egcollect,MFORCE
     &,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(tausum,taucollect,MFORCE
     &,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 12 ifr=1,nforce
        wgcum(ifr)=wgcum(ifr)+wgcollect(ifr)
        egcum(ifr)=egcum(ifr)+egcollect(ifr)
        taucum(ifr)=taucum(ifr)+taucollect(ifr)
  12  continue

      call mpi_reduce(pesum,pecollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tpbsum,tpbcollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tjfsum,tjfcollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(wg2sum,wg2collect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(eg2sum,eg2collect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(pe2sum,pe2collect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tpb2sum,tpb2collect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tjf2sum,tjf2collect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(fsum,fcollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(f2sum,f2collect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(derivsumave,derivsumcollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(derivsumave2,derivsum2collect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(derivsum,derivcollect,10*MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(esum,ecollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wsum,wcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(efsum,efcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wfsum,wfcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(e2sum,e2collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(w2sum,w2collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(ef2sum,ef2collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wf2sum,wf2collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      if(iforce_analy.eq.1) then
        call mpi_reduce(deriv_esum,deriv_ecol,10*3*MCENT*int(PTH)
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        call mpi_reduce(deriv_sumave,deriv_sumcol,3*MCENT*int(PTH)
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        call mpi_reduce(deriv_sumave2,deriv_sumcol2,3*MCENT*int(PTH)
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      endif 

      call optjas_cum(wgsum(1),egnow)
      call optorb_cum(wgsum(1),egsum(1))
      call optci_cum(wgsum(1))

      call prop_reduce(wgsum(1))
      call pcm_reduce(wgsum(1))
      call mmpol_reduce(wgsum(1))

      call max_sigma

      if(.not.wid) goto 17

      wcm2=wcm2+w2collect
      wfcm2=wfcm2+wf2collect
      ecm2=ecm2+e2collect
      efcm2=efcm2+ef2collect

      wcum=wcum+wcollect
      wfcum=wfcum+wfcollect
      ecum=ecum+ecollect
      efcum=efcum+efcollect

      if(iforce_analy.eq.1) then
        write(6,*) 'ANALYTICAL FORCE'
c        open(123,file='force_acum99',form='formatted',status='unknown')
        if(ipathak.gt.0) then
          do 202 iph=1,ipathak
            do 201 k=1,3
              do 201 ic=1,ncent
                deriv_cm(k,ic,iph)=deriv_cm(k,ic,iph)+deriv_sumcol(k,ic,iph)
                deriv_cm2(k,ic,iph)=deriv_cm2(k,ic,iph)+deriv_sumcol2(k,ic,iph)
                deriv_err(k,ic,iph)=errg(deriv_cm(k,ic,iph),deriv_cm2(k,ic,iph),1)
                do 201 j=1,3
                  deriv_ecum(j,k,ic,iph)=deriv_ecum(j,k,ic,iph)+deriv_ecol(j,k,ic,iph)
                  deriv_ave(j,k,ic,iph)=-(deriv_ecum(j,k,ic,iph)/wgcum(1))
                  if(j.eq.3) deriv_ave(j,k,ic,iph)=(egcum(1)/wgcum(1)*deriv_ecum(j,k,ic,iph)/wgcum(1))
  201             deriv_avetot(k,ic,iph)=-(deriv_ecum(1,k,ic,iph)+deriv_ecum(2,k,ic,iph)
     &                                     -egcum(1)/wgcum(1)*deriv_ecum(3,k,ic,iph))/wgcum(1)
c          write(123,'(i,f10.4,5f18.10)') iblk,eps_pathak(iph),deriv_ave(1,3,2,iph),deriv_ave(2,3,2,iph),
c     &                              deriv_ave(3,3,2,iph),deriv_avetot(3,2,iph),deriv_err(3,2,iph)
  202     write(6,'(f6.4,5f18.10)') eps_pathak(iph),deriv_ave(1,3,2,iph),deriv_ave(2,3,2,iph),
     &                              deriv_ave(3,3,2,iph),deriv_avetot(3,2,iph),deriv_err(3,2,iph)
        else
          do 203 k=1,3
            do 203 ic=1,ncent
              deriv_cm(k,ic,1)=deriv_cm(k,ic,1)+deriv_sumcol(k,ic,1)
              deriv_cm2(k,ic,1)=deriv_cm2(k,ic,1)+deriv_sumcol2(k,ic,1)
              deriv_err(k,ic,1)=errg(deriv_cm(k,ic,1),deriv_cm2(k,ic,1),1)
              do 203 j=1,3
                deriv_ecum(j,k,ic,1)=deriv_ecum(j,k,ic,1)+deriv_ecol(j,k,ic,1)
                deriv_ave(j,k,ic,1)=-(deriv_ecum(j,k,ic,1)/wgcum(1))
                if(j.eq.3) deriv_ave(j,k,ic,1)=(egcum(1)/wgcum(1)*deriv_ecum(j,k,ic,1)/wgcum(1))
  203           deriv_avetot(k,ic,1)=-(deriv_ecum(1,k,ic,1)+deriv_ecum(2,k,ic,1)
     &                                  -egcum(1)/wgcum(1)*deriv_ecum(3,k,ic,1))/wgcum(1)
        endif
      endif

      do 15 ifr=1,nforce
        wgcm2(ifr)=wgcm2(ifr)+wg2collect(ifr)
        egcm2(ifr)=egcm2(ifr)+eg2collect(ifr)
        pecm2(ifr)=pecm2(ifr)+pe2collect(ifr)
        tpbcm2(ifr)=tpbcm2(ifr)+tpb2collect(ifr)
        tjfcm2(ifr)=tjfcm2(ifr)+tjf2collect(ifr)

        pecum(ifr)=pecum(ifr)+pecollect(ifr)
        tpbcum(ifr)=tpbcum(ifr)+tpbcollect(ifr)
        tjfcum(ifr)=tjfcum(ifr)+tjfcollect(ifr)
        do 13 k=1,3
  13      derivcum(k,ifr)=derivcum(k,ifr)+derivcollect(k,ifr)

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
          fgcum(ifr)=fgcum(ifr)+fcollect(ifr)
          fgcm2(ifr)=fgcm2(ifr)+f2collect(ifr)
          fgave=egcum(1)/wgcum(1)-egcum(ifr)/wgcum(ifr)
          fgave=fgave/deltot(ifr)
          if(iblk.eq.1) then
            fgerr=0
            ifgerr=0
           else
            fgerr=errg(fgcum(ifr),fgcm2(ifr),1)
            fgerr=fgerr/abs(deltot(ifr))
            ifgerr=nint(1e12* fgerr)
          endif
          egave1=egcum(1)/wgcum(1)
          derivtotave_num=-(derivcum(1,ifr)-derivcum(1,1)+derivcum(2,ifr)-derivcum(2,1)-egave1*(derivcum(3,ifr)-derivcum(3,1)))
          derivtotave=derivtotave_num/wgcum(1)
c          open(90,file='force_numeric',form='formatted',status='unknown')
c          write(90,'(i5,1p6e20.10)') iblk,derivtotave*1e5
          write(6,*) 'NUMERICAL FORCE'
          write(6,*) 'ACCUM',0,(derivcum(2,ifr)-derivcum(2,1))*1e5,(derivcum(3,ifr)-derivcum(3,1))* 1e5
          write(6,*) 'FORCE',1,-(derivcum(1,ifr)-derivcum(1,1))/wgcum(1)*1e5
          write(6,*) 'FORCE',2,-(derivcum(2,ifr)-derivcum(2,1))/wgcum(1)*1e5
          write(6,*) 'FORCE',3,egave1*(derivcum(3,ifr)-derivcum(3,1))/wgcum(1)*1e5
          write(6,*) 'TOTAL',0,derivtotave*1e5
          derivcm(ifr)=derivcm(ifr)+derivsumcollect(ifr)
          derivcm2(ifr)=derivcm2(ifr)+derivsum2collect(ifr)
          if(iblk.eq.1) then
            derivgerr=0
            iderivgerr=0
           else
            derivgerr=errg(derivcm(ifr),derivcm2(ifr),1)
            derivgerr=derivgerr/abs(deltot(ifr))
            iderivgerr=nint(1e12* derivgerr)
          endif
         else
          call prop_prt_dmc(iblk,0,wgcum,wgcm2)
          call pcm_prt(iblk,wgcum,wgcm2)
          call mmpol_prt(iblk,wgcum,wgcm2)
        endif


c write out header first time

        if (iblk.eq.1.and.ifr.eq.1)
     &  write(6,'(t5,''egnow'',t15,''egave'',t21,''(egerr)'' ,t32
     &  ,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr)'',t66
     &  ,''tjfave'',t72,''(tjferr)'',t83,''fgave'',t96,''(fgerr)''
     &  ,t114,''fgave_n'',t131,''(fgerr_n)'',t146,''npass'',t156,''wgsum'',t164,''ioldest''
     &  ,t177,''max_sigma'')')

c write out current values of averages etc.

        iegerr=nint(100000* egerr)
        ipeerr=nint(100000* peerr)
        itpber=nint(100000*tpberr)
        itjfer=nint(100000*tjferr)

        if(ifr.eq.1) then
          write(6,'(f10.5,4(f10.5,''('',i5,'')''),62x,3i10,5x,f10.5)')
     &    egcollect(ifr)/wgcollect(ifr),
     &    egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,npass,
     &    nint(wgcollect(ifr)/nproc),ioldest,sigma
         else
          write(6,'(f10.5,4(f10.5,''('',i5,'')''),f17.12,
     &    ''('',i12,'')'',f17.12,''('',i12,'')'',10x,i10)')
     &    egcollect(ifr)/wgcollect(ifr),
     &    egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,
     &    fgave,ifgerr,derivtotave,iderivgerr,nint(wgcollect(ifr)/nproc)
        endif
   15 continue

c zero out xsum variables for metrop

   17 wsum=zero
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

      if(iforce_analy.eq.1) then
        if(ipathak.gt.0) then
          do 301 iph=1,ipathak
            do 301 j=1,3
              do 301 ic=1,ncent
                do 301 k=1,3
301               deriv_esum(j,k,ic,iph)=zero
        else
        do 302 j=1,3
          do 302 ic=1,ncent
            do 302 k=1,3
302           deriv_esum(j,k,ic,1)=zero
        endif
      endif 

      call optorb_init(1)
      call optci_init(1)

      call prop_init(1)
      call pcm_init(1)
      call mmpol_init(1)

      return
      end

c-----------------------------------------------------------------------
      subroutine max_sigma

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'
      include 'mpif.h'

      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
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
      logical wid
      common /mpiconf/ idtask,nproc,wid
      common /sigma_branch/ sigma
      dimension wg2collect(MFORCE), wg2sum(MFORCE), wgcm2temp(MFORCE)

c statement function for error calculation
      rn_eff(w,w2)=w**2/w2
      error(x,x2,w,w2)=dsqrt(max((x2/w-(x/w)**2)/(rn_eff(w,w2)-1),0.d0))
      errg1(x,x2)=error(x,x2,wgcum1temp,wgcm21temp)
      ifr=1

      call mpi_reduce(egcum1(ifr),egcum1temp,1,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(egcm21(ifr),egcm21temp,1,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wgcum1(ifr),wgcum1temp,1,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wgcm21(ifr),wgcm21temp,1,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)

      wg2sum(ifr)=wgsum(ifr)**2
      call mpi_reduce(wg2sum,wg2collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      if(wid) then
        wgcm2temp(ifr)=wgcm2(ifr)+wg2collect(ifr)
        egerr1=errg1(egcum1temp,egcm21temp)
c determining maximum sigma for local energy
        evalg_eff=nconf*nstep*rn_eff(wgcum(1),wgcm2temp(1))
        rtevalg_eff1=dsqrt(evalg_eff-1)
        sigma=egerr1*rtevalg_eff1
        do 150 id=1,nproc-1
  150     call mpi_isend(sigma,1,mpi_double_precision,id
     &    ,1,MPI_COMM_WORLD,irequest,ierr)
        else
          call mpi_recv(sigma,1,mpi_double_precision,0
     &  ,1,MPI_COMM_WORLD,istatus,ierr)
      endif

      end subroutine
