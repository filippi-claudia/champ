      subroutine acuest_gpop
c MPI version created by Claudia Filippi starting from serial version
c routine to accumulate estimators for energy etc.

      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X,
     &NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20,
     &radmax, delri, NEQSX, MTERMS, MCENT3, NCOEF, MEXCIT
      use dmc_mod, only: MWALK, MFPROD, MFPRD1, MPATH
      use const, only: delta, deltai, etrial, fbias, hb, imetro, ipr, nelec, pi
      use forcest, only: fgcm2, fgcum
      use forcepar, only: deltot, istrech, nforce
      use age, only: iage, ioldest, ioldestmx
      use contrldmc, only: iacc_rej, icross, icuspg, icut_br, icut_e, idiv_v, idmc, ipq,
     &itau_eff, nfprod, rttau, tau, taueff, tautot
      use iterat, only: iblk, ipass
      use estsum, only: efsum, efsum1, egsum, egsum1, ei1sum, ei2sum, ei3sum, esum1_dmc, esum_dmc,
     &pesum_dmc, r2sum, risum, tausum, tjfsum_dmc, tpbsum_dmc, w_acc_sum, w_acc_sum1, wdsum,
     &wdsum1, wfsum, wfsum1, wg_acc_sum, wg_acc_sum1, wgdsum, wgsum, wgsum1, wsum1, wsum_dmc
      use estcum, only: ecum1_dmc, ecum_dmc, efcum, efcum1, egcum, egcum1, ei1cum, ei2cum,
     &ei3cum, pecum_dmc, r2cum_dmc, ricum, taucum, tjfcum_dmc, tpbcum_dmc, w_acc_cum, w_acc_cum1,
     &wcum1, wcum_dmc, wdcum, wdcum1, wfcum, wfcum1, wg_acc_cum, wg_acc_cum1, wgcum, wgcum1,
     &wgdcum
      use force_dmc, only: itausec, nwprod
      use est2cm, only: ecm21_dmc, ecm2_dmc, efcm2, efcm21, egcm2, egcm21, ei1cm2, ei2cm2,
     &ei3cm2, pecm2_dmc, r2cm2_dmc, ricm2, tjfcm_dmc, tpbcm2_dmc, wcm2, wcm21, wdcm2, wdcm21,
     &wfcm2, wfcm21, wgcm2, wgcm21, wgdcm2
      use derivest, only: derivcm2, derivcum, derivsum, derivtotave_num_old
      use mpiconf, only: idtask, nproc, wid, NPROCX
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF

      use branch, only: eest, eigv, eold, ff, fprod, nwalk, pwt, wdsumo, wgdsumo, wt, wtgen,
     &wthist
      use contrl, only: idump, irstar, isite, nblk, nblkeq, nconf, nconf_new, nstep
      implicit real*8(a-h,o-z)



      include 'mpif.h'
      parameter (zero=0.d0,one=1.d0)



      dimension egcollect(MFORCE),wgcollect(MFORCE),pecollect(MFORCE),
     +tpbcollect(MFORCE),tjfcollect(MFORCE),taucollect(MFORCE),
     +derivcollect(10,MFORCE)

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

      call mpi_reduce(pesum_dmc,pecollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tpbsum_dmc,tpbcollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tjfsum_dmc,tjfcollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tausum,taucollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(derivsum,derivcollect,10*MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_allreduce(ioldest,ioldest_collect,1
     &,mpi_integer,mpi_max,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(ioldestmx,ioldestmx_collect,1
     &,mpi_integer,mpi_max,MPI_COMM_WORLD,ierr)

      ioldest=ioldest_collect
      ioldestmx=ioldestmx_collect

      call prop_reduce(dum)
      call pcm_reduce(dum)
      call mmpol_reduce(dum)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(.not.wid) goto 17

      wnow=wsum_dmc/nstep
      wfnow=wfsum/nstep
      enow=esum_dmc/wsum_dmc
      efnow=efsum/wfsum
      ei1now=wfsum/wdsum
      ei2now=wgsum(1)/wgdsum
      rinow=risum/wgsum(1)
      r2now=r2sum/wgsum(1)

      ei1cm2=ei1cm2+ei1now**2
      ei2cm2=ei2cm2+ei2now**2
      r2cm2_dmc=r2cm2_dmc+r2sum*r2now
      ricm2=ricm2+risum*rinow

      wdcum=wdcum+wdsum
      wgdcum=wgdcum+wgdsum
      ei1cum=ei1cum+ei1now
      ei2cum=ei2cum+ei2now
      r2cum_dmc=r2cum_dmc+r2sum
      ricum=ricum+risum

      wcm2=wcm2+wsum_dmc**2
      wfcm2=wfcm2+wfsum**2
      ecm2_dmc=ecm2_dmc+esum_dmc*enow
      efcm2=efcm2+efsum*efnow

      wcum_dmc=wcum_dmc+wsum_dmc
      wfcum=wfcum+wfsum
      ecum_dmc=ecum_dmc+esum_dmc
      efcum=efcum+efsum

      do 15 ifr=1,nforce

        pesum_dmc(ifr)=pecollect(ifr)
        tpbsum_dmc(ifr)=tpbcollect(ifr)
        tjfsum_dmc(ifr)=tjfcollect(ifr)
        tausum(ifr)=taucollect(ifr)
        do 13 k=1,3
  13      derivsum(k,ifr)=derivcollect(k,ifr)

        wgnow=wgsum(ifr)/nstep
        egnow=egsum(ifr)/wgsum(ifr)
        penow=pesum_dmc(ifr)/wgsum(ifr)
        tpbnow=tpbsum_dmc(ifr)/wgsum(ifr)
        tjfnow=tjfsum_dmc(ifr)/wgsum(ifr)

        wgcm2(ifr)=wgcm2(ifr)+wgsum(ifr)**2
        egcm2(ifr)=egcm2(ifr)+egsum(ifr)*egnow
        pecm2_dmc(ifr)=pecm2_dmc(ifr)+pesum_dmc(ifr)*penow
        tpbcm2_dmc(ifr)=tpbcm2_dmc(ifr)+tpbsum_dmc(ifr)*tpbnow
        tjfcm_dmc(ifr)=tjfcm_dmc(ifr)+tjfsum_dmc(ifr)*tjfnow

        wgcum(ifr)=wgcum(ifr)+wgsum(ifr)
        egcum(ifr)=egcum(ifr)+egsum(ifr)
        pecum_dmc(ifr)=pecum_dmc(ifr)+pesum_dmc(ifr)
        tpbcum_dmc(ifr)=tpbcum_dmc(ifr)+tpbsum_dmc(ifr)
        tjfcum_dmc(ifr)=tjfcum_dmc(ifr)+tjfsum_dmc(ifr)
        taucum(ifr)=taucum(ifr)+tausum(ifr)
        do 14 k=1,3
  14      derivcum(k,ifr)=derivcum(k,ifr)+derivsum(k,ifr)

        if(iblk.eq.1) then
          egerr=0
          peerr=0
          tpberr=0
          tjferr=0
         else
          egerr=errg(egcum(ifr),egcm2(ifr),ifr)
          peerr=errg(pecum_dmc(ifr),pecm2_dmc(ifr),ifr)
          tpberr=errg(tpbcum_dmc(ifr),tpbcm2_dmc(ifr),ifr)
          tjferr=errg(tjfcum_dmc(ifr),tjfcm_dmc(ifr),ifr)
        endif

        egave=egcum(ifr)/wgcum(ifr)
        peave=pecum_dmc(ifr)/wgcum(ifr)
        tpbave=tpbcum_dmc(ifr)/wgcum(ifr)
        tjfave=tjfcum_dmc(ifr)/wgcum(ifr)

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
            ifgerr=nint(100000* fgerr)
          endif
          egave1=egcum(1)/wgcum(1)
          derivtotave=-(derivcum(1,ifr)-derivcum(1,1)+derivcum(2,ifr)-derivcum(2,1)-egave1*(derivcum(3,ifr)-derivcum(3,1)))/wgcum(1)
         else
          call prop_cum(wgsum(ifr))
          call pcm_cum(wgsum(ifr))
          call mmpol_cum(wgsum(ifr))
        endif

c write out header first time

        if (iblk.eq.1.and.ifr.eq.1)
     &  write(6,'(t5,''egnow'',t15,''egave'',t21,''(egerr)'' ,t32
     &  ,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr)'',t66
     &  ,''tjfave'',t72,''(tjferr)'',t83,''fgave'',t89,''(fgerr)'',
     &  t101,''npass'',t111,''wgsum'',t121,''ioldest'')')

c write out current values of averages etc.

        iegerr=nint(100000* egerr)
        ipeerr=nint(100000* peerr)
        itpber=nint(100000*tpberr)
        itjfer=nint(100000*tjferr)

        if(ifr.eq.1) then
          write(6,'(f10.5,4(f10.5,''('',i5,'')''),17x,3i10)')
     &    egsum(ifr)/wgsum(ifr),
     &    egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,npass,
     &    nint(wgsum(ifr)),ioldest

          call prop_prt_dmc(iblk,0,wgcum,wgcm2)
          call pcm_prt(iblk,wgcum,wgcm2)
          call mmpol_prt(iblk,wgcum,wgcm2)
         else
          write(6,'(f10.5,5(f10.5,''('',i5,'')''),f10.5,10x,i10)')
     &    egsum(ifr)/wgsum(ifr),
     &    egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,
     &    fgave,ifgerr,derivtotave,nint(wgsum(ifr))
        endif
   15 continue

c     call flush(6)

c zero out xsum variables for metrop

   17 wsum_dmc=zero
      wfsum=zero
      wdsum=zero
      wgdsum=zero
      esum_dmc=zero
      efsum=zero
      ei1sum=zero
      ei2sum=zero
      r2sum=zero
      risum=zero

      do 18 ifr=1,nforce
        egsum(ifr)=zero
        wgsum(ifr)=zero
        pesum_dmc(ifr)=zero
        tpbsum_dmc(ifr)=zero
        tjfsum_dmc(ifr)=zero
        tausum(ifr)=zero
        do 18 k=1,10
   18     derivsum(k,ifr)=zero

      call prop_init(1)
      call pcm_init(1)
      call mmpol_init(1)
      
      return
      end
