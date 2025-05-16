module acues1_mod
contains
      subroutine acues1
! MPI version created by Claudia Filippi starting from serial version
! routine to accumulate estimators for energy etc.

      use fragments, only: esum_i, ecum_i, eest_i
      use fragments, only: esumfrag, ecumfrag, eestfrag
      use fragments, only: egcum1frag, egsum1frag, egsumfrag, egcm21frag, nfrag
      use precision_kinds, only: dp
      use const, only: etrial
      use contrldmc, only: icut_e
      use control, only: ipr, mode
      use multiple_geo, only: nforce
      use const,   only: etrial
      use contrl_file, only: ounit
      use contrldmc, only: idmc,nfprod
      use control, only: ipr,mode
      use estcum, only: ipass
      use estsum, only: efsum, efsum1, egsum, egsum1, esum1_dmc, esum_dmc
      use estsum, only: tausum, wfsum, wfsum1, wgsum, wgsum1, wsum1, wsum_dmc
      use estcum, only: ecum1_dmc, efcum1, egcum, egcum1
      use estcum, only: taucum, wcum1, wfcum1, wgcum, wgcum1
      use est2cm, only: ecm21_dmc, efcm21, egcm21, wcm21
      use est2cm, only: wfcm21, wgcm21
      use branch, only: eest, eigv, ff, fprod, wdsumo, wgdsumo, wtgen
      use contrl_file,    only: ounit

      use acues1_gpop_mod, only: acues1_gpop
      use branch,  only: eest,eigv,ff,fprod,wdsumo,wgdsumo,wtgen
      use multiple_geo, only: nforce
      use precision_kinds, only: dp

      implicit none

      integer :: ifr, ipmod, nfpro
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0

      if(mode.eq.'dmc_one_mpi2') then
        call acues1_gpop
        return
      endif
!     statistical fluctuations without blocking
      if(idmc.gt.0) then
         wcum1=wcum1+wsum1(1)
         wfcum1=wfcum1+wfsum1
         ecum1_dmc=ecum1_dmc+esum1_dmc(1)
         efcum1=efcum1+efsum1

         wcm21=wcm21+wsum1(1)**2
         wfcm21=wfcm21+wfsum1**2
         ecm21_dmc=ecm21_dmc+esum1_dmc(1)**2/wsum1(1)
         efcm21=efcm21+efsum1**2/wfsum1
      endif

      do ifr=1,nforce !TODO: Fragment this block
        wgcum1(ifr)=wgcum1(ifr)+wgsum1(ifr)
        egcum1(ifr)=egcum1(ifr)+egsum1(ifr)
        wgcm21(ifr)=wgcm21(ifr)+wgsum1(ifr)**2
        egcm21(ifr)=egcm21(ifr)+egsum1(ifr)**2/wgsum1(ifr)
      enddo
      
! collect block averages
      wsum_dmc=wsum_dmc+wsum1(1)
      wfsum=wfsum+wfsum1
      esum_dmc=esum_dmc+esum1_dmc(1)
      efsum=efsum+efsum1

      do ifr=1,nforce
        wgsum(ifr)=wgsum(ifr)+wgsum1(ifr)
        egsum(ifr)=egsum(ifr)+egsum1(ifr)
      enddo

      if (nfrag.gt.1) then
        egcum1frag(:)=egcum1frag(:)+egsum1frag(:)
        egcm21frag(:)=egcm21frag(:)+egsum1frag(:)**2/wgsum1(1)
        egsumfrag(:)=egsumfrag(:)+egsum1frag(:)
      endif

! Estimate eigenvalue of G from the energy
      ipmod=mod(ipass,nfprod)
      if(iabs(idmc).eq.1) then
        nfpro=min(nfprod,ipass)
        eigv=(wgsum1(1)/wtgen(ipmod))**(one/nfpro)
       else
        eest=(egcum(1)+egsum(1))/(wgcum(1)+wgsum(1))
        if (icut_e.lt.0) then
          eest_i(:) = (ecum_i(:) + esum_i(:))/(wgcum(1)+wgsum(1))
        endif
        if (nfrag.gt.1) then
          eestfrag(:) = (ecumfrag(:) + esumfrag(:)) / (wgcum(1) + wgsum(1))
        endif
        !print *, 'acues1', egsum(1), sum(esum_i(:)), egsum(1) - sum(esum_i(:))
        !print *, 'eest (acues1)', eest, sum(eest_i(:)), eest - sum(eest_i(:))
        eigv=dexp((etrial-eest)*(taucum(1)+tausum(1))/ &
                                (wgcum(1)+wgsum(1)))
        if(ipr.ge.1) write(ounit,'(''eigv'',9f14.6)') eigv,eest, &
        egcum(1),egsum(1),wgcum(1),wgsum(1),fprod
      endif

      wdsumo=wsum1(1)
      wgdsumo=wsum1(1)*fprod/ff(mod(ipass+1,nfprod))
      wtgen(ipmod)=wsum1(1)

! zero out step averages
      wfsum1=zero
      efsum1=zero
      do ifr=1,nforce
        wsum1(ifr)=zero
        wgsum1(ifr)=zero
        esum1_dmc(ifr)=zero
        egsum1(ifr)=zero
      enddo
      
      if (nfrag.gt.1) then
        egsum1frag(:)=zero
      endif

      return
      end
end module
