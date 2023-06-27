module zerest_mod
contains
      subroutine zerest
! Written by Cyrus Umrigar, modified by Claudia Filippi

      use age,     only: iage,ioldest,ioldestmx
      use denupdn, only: rprobdn,rprobup
      use derivest, only: derivcm2,derivcum,derivsum
      use est2cm,  only: ecm21_dmc,ecm2_dmc,efcm2,efcm21,egcm2,egcm21
      use est2cm,  only: ei1cm2,ei2cm2,ei3cm2,pecm2_dmc,r2cm2_dmc,ricm2
      use est2cm,  only: tpbcm2_dmc,wcm2,wcm21,wdcm2,wfcm2
      use est2cm,  only: wfcm21,wgcm2,wgcm21,wgdcm2
      use estcum,  only: ecum1_dmc,ecum_dmc,efcum,efcum1,egcum,egcum1
      use estcum,  only: ei1cum,ei2cum,ei3cum,iblk,pecum_dmc,r2cum_dmc
      use estcum,  only: ricum,taucum,tpbcum_dmc,wcum1
      use estcum,  only: wcum_dmc,wdcum,wfcum,wfcum1,wgcum,wgcum1,wgdcum
      use estsum,  only: efsum,efsum1,egsum,egsum1,ei1sum,ei2sum,ei3sum
      use estsum,  only: esum1_dmc,esum_dmc,pesum_dmc,r2sum,risum,tausum
      use estsum,  only: tpbsum_dmc,wdsum,wfsum,wfsum1,wgdsum
      use estsum,  only: wgsum,wgsum1,wsum1,wsum_dmc
      use mmpol,   only: mmpol_init
      use mpiblk,  only: iblk_proc
      use multiple_geo, only: fgcm2,fgcum,nforce
      use optci_mod, only: optci_init
      use optjas_mod, only: optjas_init
      use optorb_f_mod, only: optorb_init
      use optx_jas_ci, only: optx_jas_ci_init
      use optx_jas_orb, only: optx_jas_orb_init
      use optx_orb_ci, only: optx_orb_ci_init
      use pcm_mod, only: pcm_init
      use precision_kinds, only: dp
      use properties_mod, only: prop_init
      use stats,   only: acc,nacc,nbrnch,nodecr,trymove
      use step,    only: rprob
      use vmc_mod, only: nrad


      implicit none

      integer :: i, ifr, k

      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0


! routine to accumulate estimators for energy etc.

      iblk=0
      iblk_proc=0

! zero out estimators

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

! debug
      iage=0
      ioldest=0
      ioldestmx=0

      do ifr=1,nforce
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
        pecm2_dmc(ifr)=zero
        tpbcm2_dmc(ifr)=zero
        pesum_dmc(ifr)=zero
        tpbsum_dmc(ifr)=zero
        fgcum(ifr)=zero
        fgcm2(ifr)=zero
        derivcm2(ifr)=zero
        do k=1,10
          derivsum(k,ifr)=zero
          derivcum(k,ifr)=zero
        enddo
      enddo

      nbrnch=0
      trymove=0
      acc=0
      nacc=0
      nodecr=0

! Zero out estimators for charge density of atom.
      do i=1,nrad
        rprobup(i)=zero
        rprobdn(i)=zero
        rprob(i)=zero
      enddo

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
end module
