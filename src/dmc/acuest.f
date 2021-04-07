      subroutine acuest
c MPI version created by Claudia Filippi starting from serial version
c routine to accumulate estimators for energy etc.

      use precision_kinds, only: dp
      use forcest, only: fgcm2, fgcum
      use forcepar, only: deltot, nforce
      use age, only: ioldest
      use estcum, only: iblk
      use estsum, only: efsum, egsum, ei1sum, ei2sum, esum_dmc
      use estsum, only: pesum_dmc, r2sum, risum, tausum, tjfsum_dmc, tpbsum_dmc, wdsum
      use estsum, only: wfsum, wgdsum, wgsum, wsum_dmc
      use estcum, only: ecum_dmc, efcum, egcum, ei1cum, ei2cum
      use estcum, only: pecum_dmc, r2cum_dmc, ricum, taucum, tjfcum_dmc, tpbcum_dmc
      use estcum, only: wcum_dmc, wdcum, wfcum, wgcum
      use estcum, only: wgdcum
      use est2cm, only: ecm2_dmc, efcm2, egcm2, ei1cm2, ei2cm2
      use est2cm, only: pecm2_dmc, r2cm2_dmc, ricm2, tjfcm_dmc, tpbcm2_dmc, wcm2
      use est2cm, only: wfcm2, wgcm2
      use derivest, only: derivcm2, derivcum, derivsum, derivtotave_num_old
      use mpiconf, only: nproc, wid
      use contr3, only: mode
      use mpiblk, only: iblk_proc
      use force_mod, only: MFORCE
      use contrl, only: nstep
      use mpi

      implicit none

      integer :: i, npass, ifr, ierr, k, ifgerr, iderivgerr
      integer :: iegerr, ipeerr, itpber, itjfer
      real(dp) :: egcollect, wgcollect, pecollect, tjfcollect, eg2collect
      real(dp) :: wg2collect, tpb2collect, tjf2collect, fsum, eg2sum, wg2sum, pe2sum
      real(dp) :: tpb2sum, taucollect, fcollect, f2collect, rn_eff, w
      real(dp) :: w2, error, x, x2, errg, wnow, wfnow
      real(dp) :: enow, efnow, ei1now, ei2now, rinow, r2now
      real(dp) :: w2sum, wf2sum, e2sum, ef2sum, wgnow, egnow, penow
      real(dp) :: tpbnow, tjfnow, tjf2sum, f2sum, tpbcollect, pe2collect
      real(dp) :: derivcollect, ecollect, wcollect, efcollect, wfcollect, e2collect, w2collect
      real(dp) :: ef2collect, wf2collect, egerr, peerr, tpberr, tjferr
      real(dp) :: egave, peave, tpbave, tjfave, fgave, fgerr
      real(dp) :: egave1, derivtotave_num, derivtotave, delta_derivtotave_num, derivgerr
      real(dp) :: fgerr_n
      real(dp), parameter :: zero = 0.0d0 
      real(dp), parameter :: one = 1.0d0 

      dimension egcollect(MFORCE),wgcollect(MFORCE),pecollect(MFORCE),
     &tpbcollect(MFORCE),tjfcollect(MFORCE),eg2collect(MFORCE),wg2collect(MFORCE),
     &pe2collect(MFORCE),tpb2collect(MFORCE),tjf2collect(MFORCE),fsum(MFORCE),
     &f2sum(MFORCE),eg2sum(MFORCE),wg2sum(MFORCE),pe2sum(MFORCE),tpb2sum(MFORCE),
     &tjf2sum(MFORCE),taucollect(MFORCE),fcollect(MFORCE),f2collect(MFORCE),
     &derivcollect(10,MFORCE)

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

      w2sum=wsum_dmc**2
      wf2sum=wfsum**2
      e2sum=esum_dmc*enow
      ef2sum=efsum*efnow

      do 10 ifr=1,nforce
        wgnow=wgsum(ifr)/nstep
        egnow=egsum(ifr)/wgsum(ifr)
        penow=pesum_dmc(ifr)/wgsum(ifr)
        tpbnow=tpbsum_dmc(ifr)/wgsum(ifr)
        tjfnow=tjfsum_dmc(ifr)/wgsum(ifr)

        wg2sum(ifr)=wgsum(ifr)**2
        eg2sum(ifr)=egsum(ifr)*egnow
        pe2sum(ifr)=pesum_dmc(ifr)*penow
        tpb2sum(ifr)=tpbsum_dmc(ifr)*tpbnow
        tjf2sum(ifr)=tjfsum_dmc(ifr)*tjfnow
        if(ifr.gt.1) then
          fsum(ifr)=wgsum(1)*(egnow-egsum(1)/wgsum(1))
          f2sum(ifr)=wgsum(1)*(egnow-egsum(1)/wgsum(1))**2
        endif
  10  continue

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

      call mpi_reduce(pesum_dmc,pecollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tpbsum_dmc,tpbcollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tjfsum_dmc,tjfcollect,MFORCE
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

      call mpi_reduce(derivsum,derivcollect,10*MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(esum_dmc,ecollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wsum_dmc,wcollect,1
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


      call optjas_cum(wgsum(1),egnow)
      call optorb_cum(wgsum(1),egsum(1))
      call optci_cum(wgsum(1))

      call prop_reduce(wgsum(1))
      call pcm_reduce(wgsum(1))
      call mmpol_reduce(wgsum(1))

      if(.not.wid) goto 17

      wcm2=wcm2+w2collect
      wfcm2=wfcm2+wf2collect
      ecm2_dmc=ecm2_dmc+e2collect
      efcm2=efcm2+ef2collect

      wcum_dmc=wcum_dmc+wcollect
      wfcum=wfcum+wfcollect
      ecum_dmc=ecum_dmc+ecollect
      efcum=efcum+efcollect

      do 15 ifr=1,nforce
        wgcm2(ifr)=wgcm2(ifr)+wg2collect(ifr)
        egcm2(ifr)=egcm2(ifr)+eg2collect(ifr)
        pecm2_dmc(ifr)=pecm2_dmc(ifr)+pe2collect(ifr)
        tpbcm2_dmc(ifr)=tpbcm2_dmc(ifr)+tpb2collect(ifr)
        tjfcm_dmc(ifr)=tjfcm_dmc(ifr)+tjf2collect(ifr)

        pecum_dmc(ifr)=pecum_dmc(ifr)+pecollect(ifr)
        tpbcum_dmc(ifr)=tpbcum_dmc(ifr)+tpbcollect(ifr)
        tjfcum_dmc(ifr)=tjfcum_dmc(ifr)+tjfcollect(ifr)
        do 13 k=1,3
  13      derivcum(k,ifr)=derivcum(k,ifr)+derivcollect(k,ifr)

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
     &  ,''tjfave'',t72,''(tjferr)'',t83,''fgave'',t96,''(fgerr)'',
     &  t114,''fgave_n'',t131,''(fgerr_n)'',t146,''npass'',t156,''wgsum'',t164,''ioldest'')')

c write out current values of averages etc.

        iegerr=nint(100000* egerr)
        ipeerr=nint(100000* peerr)
        itpber=nint(100000*tpberr)
        itjfer=nint(100000*tjferr)

        if(ifr.eq.1) then
          write(6,'(f10.5,4(f10.5,''('',i5,'')''),62x,3i10)')
     &    egcollect(ifr)/wgcollect(ifr),
     &    egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,npass,
     &    nint(wgcollect(ifr)/nproc),ioldest
         else
          write(6,'(f10.5,4(f10.5,''('',i5,'')''),f17.12,
     &    ''('',i12,'')'',f17.12,''('',i12,'')'',10x,i10)')
     &    egcollect(ifr)/wgcollect(ifr),
     &    egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,
     &    fgave,ifgerr,derivtotave,iderivgerr,nint(wgcollect(ifr)/nproc)
        endif
   15 continue

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

      do 20 ifr=1,nforce
        egsum(ifr)=zero
        wgsum(ifr)=zero
        pesum_dmc(ifr)=zero
        tpbsum_dmc(ifr)=zero
        tjfsum_dmc(ifr)=zero
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
