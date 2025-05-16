module zerest_mod
contains
      subroutine zerest
! Written by Cyrus Umrigar, modified by Claudia Filippi

      use age,     only: iage,ioldest,ioldestmx
      use contrldmc, only: icut_e
      use derivest, only: derivcm2,derivcum,derivsum
      use est2cm,  only: ecm21_dmc,ecm2_dmc,efcm2,efcm21,egcm2,egcm21
      use est2cm,  only: pecm2_dmc
      use est2cm,  only: tpbcm2_dmc,wcm2,wcm21,wfcm2
      use est2cm,  only: wfcm21,wgcm2,wgcm21
      use estcum,  only: ecum1_dmc,ecum_dmc,efcum,efcum1,egcum,egcum1
      use estcum,  only: iblk,pecum_dmc
      use estcum,  only: taucum,tpbcum_dmc,wcum1
      use estcum,  only: wcum_dmc,wfcum,wfcum1,wgcum,wgcum1
      use estsum,  only: efsum,efsum1,egsum,egsum1
      use estsum,  only: esum1_dmc,esum_dmc,pesum_dmc,tausum
      use estsum,  only: tpbsum_dmc,wfsum,wfsum1
      use estsum,  only: wgsum,wgsum1,wsum1,wsum_dmc
      use fragments,  only: esum_i, ecum_i, esumfrag, ecumfrag
      use fragments,  only: egsum1frag, egsumfrag, egcum1frag, egcumfrag, egcm21frag, egcm2frag, nfrag
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
      use force_analytic, only: force_analy_init
      use system,    only: ncent
      use force_pth, only: PTH

      implicit none

      integer :: i, ifr, k, ic, iph

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
      ecum1_dmc=zero
      efcum1=zero
      ecum_dmc=zero
      efcum=zero

      wcm21=zero
      wfcm21=zero
      wcm2=zero
      wfcm2=zero
      ecm21_dmc=zero
      efcm21=zero
      ecm2_dmc=zero
      efcm2=zero

      wfsum1=zero
      wsum_dmc=zero
      wfsum=zero
      efsum1=zero
      esum_dmc=zero
      efsum=zero
      
      if (icut_e.lt.0) then
        ecum_i = zero
        esum_i = zero
      endif
      
      if (nfrag.gt.1) then
        ecumfrag = zero
        esumfrag = zero
        egcumfrag = zero
        egcum1frag = zero
        egsumfrag = zero
        egsum1frag = zero
        egcm2frag = zero
        egcm21frag = zero
      endif
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
      enddo

      do iph=1,PTH
        do k=1,3
          do ic=1,ncent
            do ifr=1,3
              derivsum(ifr,k,ic,iph)=zero
              derivcum(ifr,k,ic,iph)=zero
            enddo
            derivcm2(k,ic,iph)=zero
          enddo
        enddo
      enddo

      nbrnch=0
      trymove=0
      acc=0
      nacc=0
      nodecr=0

      call optjas_init
      call optci_init(0)
      call optorb_init(0)
      call optx_jas_orb_init
      call optx_jas_ci_init
      call optx_orb_ci_init

      call prop_init(0)
      call pcm_init(0)
      call mmpol_init(0)

      call force_analy_init(0)

      return
      end
end module
