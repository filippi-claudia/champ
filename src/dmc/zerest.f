      subroutine zerest
c Written by Cyrus Umrigar, modified by Claudia Filippi

      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X,
     &NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20,
     &radmax, delri, NEQSX, MTERMS, MCENT3, NCOEF, MEXCIT
      use dmc_mod, only: MWALK, MFPROD, MFPRD1, MPATH
      use basis, only: zex, betaq, n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz,
     & n4s, n4p, n4fxxx, n4fyyy, n4fzzz, n4fxxy, n4fxxz, n4fyyx, n4fyyz,
     & n4fzzx, n4fzzy, n4fxyz, nsa, npa, ndzra, ndz2a, ndxya, ndxza, ndyza, ndx2a
      use forcest, only: fgcm2, fgcum
      use forcepar, only: deltot, istrech, nforce
      use age, only: iage, ioldest, ioldestmx
      use iterat, only: iblk, ipass
      use stats, only: acc, dfus2ac, dfus2un, dr2ac, dr2un, nacc, nbrnch, nodecr, trymove
      use estsum, only: efsum, efsum1, egsum, egsum1, ei1sum, ei2sum, ei3sum, esum1_dmc, esum_dmc,
     &pesum_dmc, r2sum, risum, tausum, tjfsum_dmc, tpbsum_dmc, w_acc_sum, w_acc_sum1, wdsum,
     &wdsum1, wfsum, wfsum1, wg_acc_sum, wg_acc_sum1, wgdsum, wgsum, wgsum1, wsum1, wsum_dmc
      use estcum, only: ecum1_dmc, ecum_dmc, efcum, efcum1, egcum, egcum1, ei1cum, ei2cum,
     &ei3cum, pecum_dmc, r2cum_dmc, ricum, taucum, tjfcum_dmc, tpbcum_dmc, w_acc_cum, w_acc_cum1,
     &wcum1, wcum_dmc, wdcum, wdcum1, wfcum, wfcum1, wg_acc_cum, wg_acc_cum1, wgcum, wgcum1,
     &wgdcum
      use est2cm, only: ecm21_dmc, ecm2_dmc, efcm2, efcm21, egcm2, egcm21, ei1cm2, ei2cm2,
     &ei3cm2, pecm2_dmc, r2cm2_dmc, ricm2, tjfcm_dmc, tpbcm2_dmc, wcm2, wcm21, wdcm2, wdcm21,
     &wfcm2, wfcm21, wgcm2, wgcm21, wgdcm2
      use derivest, only: derivcm2, derivcum, derivsum, derivtotave_num_old
      use step, only: ekin, ekin2, rprob, suc, trunfb, try
      use denupdn, only: rprobdn, rprobup
      use mpiblk, only: iblk_proc
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use pseudo_mod, only: MPS_L, MPS_QUAD, MPS_GRID, MGAUSS

      use branch, only: eest, eigv, eold, ff, fprod, nwalk, pwt, wdsumo, wgdsumo, wt, wtgen,
     &wthist
      implicit real*8(a-h,o-z)



      parameter (zero=0.d0,one=1.d0)

c routine to accumulate estimators for energy etc.



      iblk=0
      iblk_proc=0

c zero out estimators

      wcum1=zero
      wfcum1=zero
      wcum_dmc=zero
      wfcum=zero
      wdcum=zero
      wgdcum=zero
      ecum1_dmc=zero
      efcum1=zero
      ecum_dmc=zero
      efcum=zero
      ei1cum=zero
      ei2cum=zero
      ei3cum=zero
      r2cum_dmc=zero
      ricum=zero

      wcm21=zero
      wfcm21=zero
      wcm2=zero
      wfcm2=zero
      wdcm2=zero
      wgdcm2=zero
      ecm21_dmc=zero
      efcm21=zero
      ecm2_dmc=zero
      efcm2=zero
      ei1cm2=zero
      ei2cm2=zero
      ei3cm2=zero
      r2cm2_dmc=zero
      ricm2=zero

      wfsum1=zero
      wsum_dmc=zero
      wfsum=zero
      wdsum=zero
      wgdsum=zero
      efsum1=zero
      esum_dmc=zero
      efsum=zero
      ei1sum=zero
      ei2sum=zero
      ei3sum=zero
      r2sum=zero
      risum=zero

      do 85 ifr=1,nforce
        tausum(ifr)=zero
        taucum(ifr)=zero
        wgcum1(ifr)=zero
        wgcum(ifr)=zero
        egcum1(ifr)=zero
        egcum(ifr)=zero
        wgcm21(ifr)=zero
        wgcm2(ifr)=zero
        egcm21(ifr)=zero
        egcm2(ifr)=zero
        wsum1(ifr)=zero
        wgsum1(ifr)=zero
        wgsum(ifr)=zero
        esum1_dmc(ifr)=zero
        egsum1(ifr)=zero
        egsum(ifr)=zero
        pecum_dmc(ifr)=zero
        tpbcum_dmc(ifr)=zero
        tjfcum_dmc(ifr)=zero
        pecm2_dmc(ifr)=zero
        tpbcm2_dmc(ifr)=zero
        tjfcm_dmc(ifr)=zero
        pesum_dmc(ifr)=zero
        tpbsum_dmc(ifr)=zero
        tjfsum_dmc(ifr)=zero
        fgcum(ifr)=zero
        fgcm2(ifr)=zero
        derivcm2(ifr)=zero
        do 85 k=1,10
          derivsum(k,ifr)=zero
   85     derivcum(k,ifr)=zero

      nbrnch=0
      trymove=0
      acc=0
      nacc=0
      nodecr=0

c Zero out estimators for charge density of atom.
      do 90 i=1,nrad
        rprobup(i)=zero
        rprobdn(i)=zero
   90   rprob(i)=zero

      call optjas_init
      call optci_init(0)
      call optorb_init(0)
      call optx_jas_orb_init
      call optx_jas_ci_init
      call optx_orb_ci_init

      call prop_init(0)
      call pcm_init(0)
      call mmpol_init(0)

      return
      end
