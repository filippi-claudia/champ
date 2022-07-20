      module finwrt_mod
      contains
      subroutine finwrt
c MPI version created by Claudia Filippi starting from serial version
c routine to print out final results

      use vmc_mod, only: nrad
      use vmc_mod, only: delri
      use const, only: etrial
      use multiple_geo, only: fgcm2, fgcum
      use age, only: iage, ioldest, ioldestmx
      use contrl_per, only: iperiodic
      use contrldmc, only: idmc, nfprod, tau
      use system, only: ncent
      use estcum, only: iblk
      use config, only: vold_dmc, xold_dmc
      use stats, only: acc, nacc, nodecr, trymove
      use estcum, only: ecum1_dmc, ecum_dmc, efcum, efcum1, egcum, egcum1
      use estcum, only: pecum_dmc, taucum, tjfcum_dmc, tpbcum_dmc
      use estcum, only: wcum1, wcum_dmc, wfcum, wfcum1, wgcum, wgcum1
      use est2cm, only: ecm21_dmc, ecm2_dmc, efcm2, efcm21, egcm2, egcm21
      use est2cm, only: pecm2_dmc, tjfcm_dmc, tpbcm2_dmc, wcm2, wcm21
      use est2cm, only: wfcm2, wfcm21, wgcm2, wgcm21
      use step, only: rprob
      use mpiconf, only: nproc, wid
      use denupdn, only: rprobdn, rprobup
      use header, only: title
      use grdntspar, only: igrdtype, ngradnts
      use mpiblk, only: iblk_proc
      use multiple_geo, only: MFORCE
      use branch, only: eold, nwalk
      use optwf_corsam, only: energy, energy_err, force, force_err
!      use contrl, only: nblkeq, nconf, nstep
      use control_dmc, only: dmc_nblkeq, dmc_nconf, dmc_nstep
      use mpi
      use contrl_file, only: ounit
      use precision_kinds, only: dp
      use prop_dmc, only: prop_prt_dmc
      use pcm_dmc, only: pcm_fin
      use mmpol_dmc, only: mmpol_fin
      use finwrt_more_mod, only: finwrt_more
      use misc_grdnts, only: finwrt_grdnts_cart, finwrt_grdnts_zmat
      use misc_grdnts, only: finwrt_diaghess_zmat
      use control, only: ipr
      use control, only: mode
      use system, only: nelec
      use multiple_geo, only: nforce
      implicit none

      integer :: i, ierr, ifr, j, k
      integer :: nacc_collect, nodecr_collect
      real(dp) :: acc_collect, accav, accavn, delr, eave
      real(dp) :: eerr, eerr1, efave, eferr
      real(dp) :: eferr1, egave, egerr, egerr1
      real(dp) :: errc, errc1, errf, errf1
      real(dp) :: errg, errg1, error, errorn
      real(dp) :: errw, errw1, eval, eval_eff
      real(dp) :: eval_proc, evalf_eff, evalg_eff, fgave
      real(dp) :: fgerr, pass_proc, passes, peave
      real(dp) :: peerr, rn, rn_eff, rteval_eff1
      real(dp) :: rtevalf_eff1, rtevalg_eff1, rtpass1, term
      real(dp) :: tjfave, tjferr, tpbave, tpberr
      real(dp) :: trymove_collect, w, w2, wave
      real(dp) :: werr, werr1, wfave, wferr
      real(dp) :: wferr1, wgave, wgerr, wgerr1
      real(dp) :: x, x2
      real(dp), dimension(MFORCE) :: ffin_grdnts
      real(dp), dimension(MFORCE) :: ferr_grdnts
      real(dp), dimension(nrad) :: rprobcollect
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: half = .5d0

c Statement functions for error calculation, it might be reaplaced in the near future:
      rn_eff(w,w2)=w**2/w2
      error(x,x2,w,w2)=dsqrt(max((x2/w-(x/w)**2)/(rn_eff(w,w2)-1),0.d0))
      errorn(x,x2,rn)=dsqrt(max((x2/rn-(x/rn)**2)/(rn-1),0.d0))
      errc(x,x2)=error(x,x2,wcum_dmc,wcm2)
      errf(x,x2)=error(x,x2,wfcum,wfcm2)
      errg(x,x2,i)=error(x,x2,wgcum(i),wgcm2(i))
      errc1(x,x2)=error(x,x2,wcum1,wcm21)
      errf1(x,x2)=error(x,x2,wfcum1,wfcm21)
      errg1(x,x2,i)=error(x,x2,wgcum1(i),wgcm21(i))
      errw(x,x2)=errorn(x,x2,dfloat(iblk_proc))/dmc_nstep
      errw1(x,x2)=errorn(x,x2,pass_proc)

      do ifr=1,nforce
        energy(ifr)=0
        energy_err(ifr)=0
        ffin_grdnts(ifr)=0
        ferr_grdnts(ifr)=0
        force(ifr)=0
        force_err(ifr)=0
      enddo

      passes=dfloat(iblk*dmc_nstep)
      eval=dmc_nconf*passes
      if(mode.eq.'dmc_one_mpi1') then
        pass_proc=dfloat(iblk_proc*dmc_nstep)
        eval_proc=dmc_nconf*pass_proc
       else
        iblk_proc=iblk
        pass_proc=passes
        eval_proc=eval
      endif
c Either the next 3 lines or the 3 lines following them could be used.
c They should give nearly (but not exactly) the same result.
c Strictly the 1st 3 are for step-by-step quantities and the last 3 for blk-by-blk
c     eval_eff=nconf*rn_eff(wcum1,wcm21)
c     evalf_eff=nconf*rn_eff(wfcum1,wfcm21)
c     evalg_eff=nconf*rn_eff(wgcum1(1),wgcm21(1))
      eval_eff=dmc_nconf*dmc_nstep*rn_eff(wcum_dmc,wcm2)
      evalf_eff=dmc_nconf*dmc_nstep*rn_eff(wfcum,wfcm2)
      evalg_eff=dmc_nconf*dmc_nstep*rn_eff(wgcum(1),wgcm2(1))
      rtpass1=dsqrt(pass_proc-1)
      rteval_eff1=dsqrt(eval_eff-1)
      rtevalf_eff1=dsqrt(evalf_eff-1)
      rtevalg_eff1=dsqrt(evalg_eff-1)

      write(ounit,'(''passes,eval,pass_proc,eval_proc,eval_eff,
     &evalf_eff,evalg_eff'',19f9.0)')
     & passes,eval,pass_proc,eval_proc,eval_eff,evalf_eff,
     & evalg_eff

      if(mode.eq.'mpi_one_mpi2') then

c Collect radial charge density for atoms
      if(iperiodic.eq.0) then
        call mpi_reduce(rprob,rprobcollect,nrad,mpi_double_precision
     &  ,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do i=1,nrad
          rprob(i)=rprobcollect(i)
        enddo
      endif

      call mpi_reduce(nodecr,nodecr_collect,1,mpi_integer,mpi_sum,0,
     &MPI_COMM_WORLD,ierr)
      call mpi_reduce(trymove,trymove_collect,1,mpi_double_precision,mpi_sum,0,
     &MPI_COMM_WORLD,ierr)
      call mpi_reduce(acc,acc_collect,1,mpi_double_precision,mpi_sum,0,
     &MPI_COMM_WORLD,ierr)
      call mpi_reduce(nacc,nacc_collect,1,mpi_integer,mpi_sum,0,
     &MPI_COMM_WORLD,ierr)
      nodecr=nodecr_collect
      trymove=trymove_collect
      acc=acc_collect
      nacc=nacc_collect

      endif

      if (ipr.gt.-2)
     &  write(11,'(3i5,f11.5,f7.4,f10.7,
     &  '' nstep,nblk,nconf,etrial,tau,taueff'')')
     &  dmc_nstep,iblk,dmc_nconf,etrial,tau,taucum(1)/wgcum(1)

      if(.not.wid.and.mode.eq.'dmc_one_mpi2') return

      if(iperiodic.eq.0 .and. ncent.eq.1) then
        write(45,'(''  r   rprob'')')
        delr=one/delri
        term=one/(wgcum(1)*delr)
        do i=1,nrad
          write(45,'(f5.3,3f9.5)') delr*(i-half),rprob(i)*term,rprobup(i)*term,rprobdn(i)*term
        enddo
      endif

      if(idmc.ge.0) then
        write(ounit,'(10i6)') (iage(i),i=1,nwalk)
        do i=1,nwalk
          if(iage(i).gt.50) then
            write(ounit,'(i4,i6,f10.4,99f8.4)') i,iage(i),eold(i,1),
     &      ((xold_dmc(k,j,i,1),k=1,3),j=1,nelec)
            write(ounit,'(99f8.4)') ((vold_dmc(k,j,i,1),k=1,3),j=1,nelec)
          endif
        enddo

        write(ounit,'(''age of oldest walker (this generation, any gen)='',
     &   3i9)') ioldest,ioldestmx
      endif

c     write(ounit,'(''average of the squares of the accepted step-size='',
c    & f10.5)') dr2ac/trymove

      write(ounit,'(''taueff'',20f7.4)') (taucum(ifr)/wgcum(ifr),
     & ifr=1,nforce)

      if (wid) then
        accav=acc/trymove
        accavn=dfloat(nacc)/trymove
        wave=wcum_dmc/pass_proc
        wfave=wfcum/pass_proc
        eave=ecum_dmc/wcum_dmc
        efave=efcum/wfcum

        werr=errw(wcum_dmc,wcm2)
        wferr=errw(wfcum,wfcm2)
        werr1=errw1(wcum1,wcm21)
        wferr1=errw1(wfcum1,wfcm21)
        eerr=errc(ecum_dmc,ecm2_dmc)
        eferr=errf(efcum,efcm2)
        eerr1=errc1(ecum1_dmc,ecm21_dmc)
        eferr1=errf1(efcum1,efcm21)
      endif

      if(mode.eq.'dmc_one_mpi1') then
        write(ounit,'(''dmc_mov1_mpi '',2a10)') title
       else
        write(ounit,'(''dmc_mov1_mpi_globalpop '',2a10)') title
      endif
      write(ounit,'(''No/frac. of node crossings,acceptance='',i9,3f10.6)')
     &nodecr,dfloat(nodecr)/trymove,accav,accavn
      if(idmc.ge.0) then
        write(ounit,'(''No. of walkers at end of run='',i5)') nwalk

        write(ounit,'(''nwalk_eff/nwalk         ='',2f6.3)')
     &   rn_eff(wcum1,wcm21)/pass_proc,rn_eff(wcum_dmc,wcm2)/iblk_proc
        write(ounit,'(''nwalk_eff/nwalk with f  ='',2f6.3)')
     &   rn_eff(wfcum1,wfcm21)/pass_proc,rn_eff(wfcum,wfcm2)/iblk_proc
        write(ounit,'(''nwalk_eff/nwalk with fs ='',2f6.3)')
     &   rn_eff(wgcum1(1),wgcm21(1))/pass_proc,rn_eff(wgcum(1),
     &   wgcm2(1))/iblk_proc
      endif

      write(ounit,'(''nconf*passes'',t19,''passes  nconf nstep  nblk nblkeq  nproc  tau    taueff'',
     &/,2f12.0,5i6,2f9.5)')
     &eval,passes,dmc_nconf,dmc_nstep,iblk,dmc_nblkeq,nproc,tau,taucum(1)/wgcum(1)
      write(ounit,'(''physical variable         average     rms error   sigma*T_cor  sigma   T_cor'')')
      if(idmc.ge.0) then
        write(ounit,'(''weights ='',t22,f14.7,'' +-'',f11.7,2f9.5,f8.2)')
     &  wave,werr,werr*rtpass1,werr1*rtpass1,(werr/werr1)**2
        write(ounit,'(''wts with f ='',t22,f14.7,'' +-'',f11.7,2f9.5,f8.2)')
     &  wfave,wferr,wferr*rtpass1,wferr1*rtpass1,(wferr/wferr1)**2
        do ifr=1,nforce
          wgave=wgcum(ifr)/pass_proc
          wgerr=errw(wgcum(ifr),wgcm2(ifr))
          wgerr1=errw1(wgcum1(ifr),wgcm21(ifr))
          write(ounit,'(''wts with fs ='',t22,f14.7,'' +-'',f11.7,2f9.5,f8.2)')
     &    wgave,wgerr,wgerr*rtpass1,wgerr1*rtpass1,(wgerr/wgerr1)**2
        enddo
        write(ounit,'(''total energy (   0) ='',t24,f12.7,'' +-'',f11.7,
     &  2f9.5,f8.2)') eave,eerr,eerr*rteval_eff1,eerr1*rteval_eff1,
     &  (eerr/eerr1)**2
        write(ounit,'(''total energy (   1) ='',t24,f12.7,'' +-'',f11.7,
     &  2f9.5,f8.2)') efave,eferr,eferr*rtevalf_eff1,eferr1*rtevalf_eff1,
     &  (eferr/eferr1)**2
        endif

      do ifr=1,nforce
        egave=egcum(ifr)/wgcum(ifr)
        egerr=errg(egcum(ifr),egcm2(ifr),ifr)
        egerr1=errg1(egcum1(ifr),egcm21(ifr),ifr)
        write(ounit,'(''total energy ('',i4,'') ='',t24,f12.7,'' +-'',
     &  f11.7,2f9.5,f8.2)') nfprod,egave,egerr,egerr*rtevalg_eff1,
     &  egerr1*rtevalg_eff1,(egerr/egerr1)**2
        energy(ifr)=egave
        energy_err(ifr)=egerr
      enddo
      do ifr=1,nforce
        peave=pecum_dmc(ifr)/wgcum(ifr)
        tpbave=tpbcum_dmc(ifr)/wgcum(ifr)
        tjfave=tjfcum_dmc(ifr)/wgcum(ifr)

        peerr=errg(pecum_dmc(ifr),pecm2_dmc(ifr),ifr)
        tpberr=errg(tpbcum_dmc(ifr),tpbcm2_dmc(ifr),ifr)
        tjferr=errg(tjfcum_dmc(ifr),tjfcm_dmc(ifr),ifr)
        write(ounit,'(''potential energy ='',t24,f12.7,'' +-''
     &  ,f11.7,f9.5)') peave,peerr,peerr*rtevalg_eff1
        write(ounit,'(''jf kinetic energy ='',t24,f12.7,'' +-''
     &  ,f11.7,f9.5)') tjfave,tjferr,tjferr*rtevalg_eff1
        write(ounit,'(''pb kinetic energy ='',t24,f12.7,'' +-''
     &  ,f11.7,f9.5)') tpbave,tpberr,tpberr*rtevalg_eff1
      enddo
      do ifr=2,nforce
        fgave=egcum(1)/wgcum(1)-egcum(ifr)/wgcum(ifr)
        fgerr=errg(fgcum(ifr),fgcm2(ifr),1)
c save forces and forces errors for calculations of energy gradients.
c Done by Omar Valsson 2008-12-01
        if(ngradnts.gt.0) then
          ffin_grdnts(ifr-1)=fgave
          ferr_grdnts(ifr-1)=fgerr
        endif
        write(ounit,'(''force config'',i2,t24,e19.10
     &  ,'' +-'',e16.8,f9.5)') ifr,fgave,fgerr,fgerr*rtevalg_eff1
        force(ifr)=fgave
        force_err(ifr)=fgerr
      enddo

      call prop_prt_dmc(iblk,1,wgcum,wgcm2)
      call pcm_fin(iblk,wgcum,wgcm2)
      call mmpol_fin(iblk,wgcum,wgcm2)

      call finwrt_more

      if(ngradnts.gt.0 .and. igrdtype.eq.1) call finwrt_grdnts_cart(ffin_grdnts,ferr_grdnts)
      if(ngradnts.gt.0 .and. igrdtype.eq.2) call finwrt_grdnts_zmat(ffin_grdnts,ferr_grdnts)
      if(ngradnts.gt.0 .and. igrdtype.eq.2) call finwrt_diaghess_zmat(ffin_grdnts,ferr_grdnts)

      return
      end
      end module
