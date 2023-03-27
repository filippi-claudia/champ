      module acuest_mod
      contains
      subroutine acuest
c MPI version created by Claudia Filippi starting from serial version
c routine to accumulate estimators for energy etc.

      use acuest_gpop_mod, only: acuest_gpop
      use age,     only: ioldest
      use contrl_file, only: ounit
      use control, only: mode
      use control_dmc, only: dmc_nstep
      use derivest, only: derivcm2,derivcum,derivsum,derivtotave_num_old
      use est2cm,  only: ecm2_dmc,efcm2,egcm2,ei1cm2,ei2cm2,pecm2_dmc
      use est2cm,  only: r2cm2_dmc,ricm2,tpbcm2_dmc,wcm2,wfcm2
      use est2cm,  only: wgcm2
      use estcum,  only: ecum_dmc,efcum,egcum,ei1cum,ei2cum,iblk
      use estcum,  only: pecum_dmc,r2cum_dmc,ricum,taucum
      use estcum,  only: tpbcum_dmc,wcum_dmc,wdcum,wfcum,wgcum,wgdcum
      use estsum,  only: efsum,egsum,ei1sum,ei2sum,esum_dmc,pesum_dmc
      use estsum,  only: r2sum,risum,tausum,tpbsum_dmc,wdsum
      use estsum,  only: wfsum,wgdsum,wgsum,wsum_dmc
      use mmpol,   only: mmpol_init
      use mmpol_dmc, only: mmpol_prt
      use mmpol_reduce_mod, only: mmpol_reduce
      use mpi
      use mpiblk,  only: iblk_proc
      use mpiconf, only: nproc,wid
      use multiple_geo, only: MFORCE,fgcm2,fgcum,nforce
      use optci_mod, only: optci_cum,optci_init
      use optjas_mod, only: optjas_cum
      use optorb_f_mod, only: optorb_cum,optorb_init
      use pcm_dmc, only: pcm_prt
      use pcm_mod, only: pcm_init
      use pcm_reduce_mod, only: pcm_reduce
      use precision_kinds, only: dp
      use prop_dmc, only: prop_prt_dmc
      use prop_reduce_mod, only: prop_reduce
      use properties_mod, only: prop_init
      use contrldmc, only: idmc

      implicit none

      integer :: i, iderivgerr, iegerr, ierr, ifgerr
      integer :: ifr, ipeerr, itpber
      integer :: k, npass
      real(dp) :: delta_derivtotave_num, derivgerr, derivtotave, derivtotave_num, e2collect
      real(dp) :: e2sum, ecollect, ef2collect, ef2sum
      real(dp) :: efcollect, efnow, egave, egave1
      real(dp) :: egerr, egnow, ei1now, ei2now
      real(dp) :: enow, errg, error, fgave
      real(dp) :: fgerr, peave, peerr, penow
      real(dp) :: r2now, rinow, rn_eff
      real(dp) :: tpbave, tpberr
      real(dp) :: tpbnow, w, w2, w2collect
      real(dp) :: w2sum, wcollect, wf2collect, wf2sum
      real(dp) :: wfcollect, wfnow, wgnow, wnow
      real(dp) :: x, x2
      real(dp), dimension(MFORCE) :: egcollect
      real(dp), dimension(MFORCE) :: wgcollect
      real(dp), dimension(MFORCE) :: pecollect
      real(dp), dimension(MFORCE) :: tpbcollect
      real(dp), dimension(MFORCE) :: eg2collect
      real(dp), dimension(MFORCE) :: wg2collect
      real(dp), dimension(MFORCE) :: pe2collect
      real(dp), dimension(MFORCE) :: tpb2collect
      real(dp), dimension(MFORCE) :: fsum
      real(dp), dimension(MFORCE) :: f2sum
      real(dp), dimension(MFORCE) :: eg2sum
      real(dp), dimension(MFORCE) :: wg2sum
      real(dp), dimension(MFORCE) :: pe2sum
      real(dp), dimension(MFORCE) :: tpb2sum
      real(dp), dimension(MFORCE) :: taucollect
      real(dp), dimension(MFORCE) :: fcollect
      real(dp), dimension(MFORCE) :: f2collect
      real(dp), dimension(10, MFORCE) :: derivcollect
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0

c Statement function for error calculation, it might be reaplaced in the near future:
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

      npass=iblk_proc*dmc_nstep

      if(idmc.gt.0) then
         wnow=wsum_dmc/dmc_nstep
         wfnow=wfsum/dmc_nstep
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
      endif
      
      do ifr=1,nforce
        wgnow=wgsum(ifr)/dmc_nstep
        egnow=egsum(ifr)/wgsum(ifr)
        penow=pesum_dmc(ifr)/wgsum(ifr)
        tpbnow=tpbsum_dmc(ifr)/wgsum(ifr)

        wg2sum(ifr)=wgsum(ifr)**2
        eg2sum(ifr)=egsum(ifr)*egnow
        pe2sum(ifr)=pesum_dmc(ifr)*penow
        tpb2sum(ifr)=tpbsum_dmc(ifr)*tpbnow
        if(ifr.gt.1) then
          fsum(ifr)=wgsum(1)*(egnow-egsum(1)/wgsum(1))
          f2sum(ifr)=wgsum(1)*(egnow-egsum(1)/wgsum(1))**2
        endif
      enddo

      call mpi_allreduce(wgsum,wgcollect,MFORCE
     &,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(egsum,egcollect,MFORCE
     &,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(tausum,taucollect,MFORCE
     &,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do ifr=1,nforce
        wgcum(ifr)=wgcum(ifr)+wgcollect(ifr)
        egcum(ifr)=egcum(ifr)+egcollect(ifr)
        taucum(ifr)=taucum(ifr)+taucollect(ifr)
      enddo

      call mpi_reduce(pesum_dmc,pecollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tpbsum_dmc,tpbcollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(wg2sum,wg2collect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(eg2sum,eg2collect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(pe2sum,pe2collect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tpb2sum,tpb2collect,MFORCE
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

      do ifr=1,nforce
        wgcm2(ifr)=wgcm2(ifr)+wg2collect(ifr)
        egcm2(ifr)=egcm2(ifr)+eg2collect(ifr)
        pecm2_dmc(ifr)=pecm2_dmc(ifr)+pe2collect(ifr)
        tpbcm2_dmc(ifr)=tpbcm2_dmc(ifr)+tpb2collect(ifr)

        pecum_dmc(ifr)=pecum_dmc(ifr)+pecollect(ifr)
        tpbcum_dmc(ifr)=tpbcum_dmc(ifr)+tpbcollect(ifr)
        do k=1,3
          derivcum(k,ifr)=derivcum(k,ifr)+derivcollect(k,ifr)
        enddo

        if(iblk.eq.1) then
          egerr=0
          peerr=0
          tpberr=0
         else
          egerr=errg(egcum(ifr),egcm2(ifr),ifr)
          peerr=errg(pecum_dmc(ifr),pecm2_dmc(ifr),ifr)
          tpberr=errg(tpbcum_dmc(ifr),tpbcm2_dmc(ifr),ifr)
        endif

        egave=egcum(ifr)/wgcum(ifr)
        peave=pecum_dmc(ifr)/wgcum(ifr)
        tpbave=tpbcum_dmc(ifr)/wgcum(ifr)

        if(ifr.gt.1) then
          fgcum(ifr)=fgcum(ifr)+fcollect(ifr)
          fgcm2(ifr)=fgcm2(ifr)+f2collect(ifr)
          fgave=egcum(1)/wgcum(1)-egcum(ifr)/wgcum(ifr)
          if(iblk.eq.1) then
            fgerr=0
            ifgerr=0
           else
            fgerr=errg(fgcum(ifr),fgcm2(ifr),1)
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
            iderivgerr=nint(1e12* derivgerr)
          endif
         else
          call prop_prt_dmc(iblk,0,wgcum,wgcm2)
          call pcm_prt(iblk,wgcum,wgcm2)
          call mmpol_prt(iblk,wgcum,wgcm2)
        endif

c write out header first time

        if (iblk.eq.1.and.ifr.eq.1) then
          if(nforce.gt.1) then
            write(ounit,'(t5,''egnow'',t15,''egave'',t21,''(egerr)'' ,t32
     &      ,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr)'',t66
     &      ,''fgave'',t79,''(fgerr)'',t93,''fgave_n'',t106,''(fgerr_n)'',t123
     &      ,''npass'',t132,''wgsum'',t142,''ioldest'')')
          else
            write(ounit,'(t5,''egnow'',t15,''egave'',t21,''(egerr)'' ,t32
     &      ,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr)'',t67
     &      ,''npass'',t77,''wgsum'',t85,''ioldest'')')
          endif
        endif

c write out current values of averages etc.

        iegerr=nint(100000* egerr)
        ipeerr=nint(100000* peerr)
        itpber=nint(100000*tpberr)

        if(ifr.eq.1) then
          if(nforce.gt.1) then
            write(ounit,'(f10.5,3(f10.5,''('',i5,'')''),62x,3i10)')
     &      egcollect(ifr)/wgcollect(ifr),
     &      egave,iegerr,peave,ipeerr,tpbave,itpber,npass,
     &      nint(wgcollect(ifr)/nproc),ioldest
           else
            write(ounit,'(f10.5,3(f10.5,''('',i5,'')''),3i10)')
     &      egcollect(ifr)/wgcollect(ifr),
     &      egave,iegerr,peave,ipeerr,tpbave,itpber,npass,
     &      nint(wgcollect(ifr)/nproc),ioldest
          endif
         else
          write(ounit,'(f10.5,3(f10.5,''('',i5,'')''),f17.12,
     &    ''('',i12,'')'',f17.12,''('',i12,'')'',10x,i10)')
     &    egcollect(ifr)/wgcollect(ifr),
     &    egave,iegerr,peave,ipeerr,tpbave,itpber,
     &    fgave,ifgerr,derivtotave,iderivgerr,nint(wgcollect(ifr)/nproc)
        endif
      enddo

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

      do ifr=1,nforce
        egsum(ifr)=zero
        wgsum(ifr)=zero
        pesum_dmc(ifr)=zero
        tpbsum_dmc(ifr)=zero
        tausum(ifr)=zero
        do k=1,10
          derivsum(k,ifr)=zero
        enddo
      enddo

      call optorb_init(1)
      call optci_init(1)

      call prop_init(1)
      call pcm_init(1)
      call mmpol_init(1)

      return
      end
      end module
