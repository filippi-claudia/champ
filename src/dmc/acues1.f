      subroutine acues1
c MPI version created by Claudia Filippi starting from serial version
c routine to accumulate estimators for energy etc.

      use precision_kinds, only: dp
      use const, only: etrial, ipr
      use forcepar, only: nforce
      use contrldmc, only: idmc
      use contrldmc, only: nfprod
      use estcum, only: ipass
      use estsum, only: efsum, efsum1, egsum, egsum1, esum1_dmc, esum_dmc
      use estsum, only: tausum, wdsum
      use estsum, only: wdsum1, wfsum, wfsum1, wgdsum, wgsum, wgsum1, wsum1, wsum_dmc
      use estcum, only: ecum1_dmc, efcum1, egcum, egcum1
      use estcum, only: ei3cum, taucum
      use estcum, only: wcum1, wfcum1, wgcum, wgcum1
      use est2cm, only: ecm21_dmc, efcm21, egcm21
      use est2cm, only: ei3cm2, wcm21
      use est2cm, only: wfcm21, wgcm21
      use contr3, only: mode
      use branch, only: eest, eigv, ff, fprod, wdsumo, wgdsumo, wtgen
      use contrl_file,    only: ounit
      implicit none

      integer :: ifr, ipmod, nfpro
      real(dp) :: wgdsum1
      real(dp) :: eisum = 0.0_dp  !eisum was uninit causing issue on line 64 (check with claudia)
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0

      if(mode.eq.'dmc_one_mpi2') then
        call acues1_gpop
        return
      endif
c statistical fluctuations without blocking
      wdsum1=wdsumo
      wgdsum1=wgdsumo

      wcum1=wcum1+wsum1(1)
      wfcum1=wfcum1+wfsum1
      ecum1_dmc=ecum1_dmc+esum1_dmc(1)
      efcum1=efcum1+efsum1
      ei3cum=ei3cum+wfsum1/wdsum1

      wcm21=wcm21+wsum1(1)**2
      wfcm21=wfcm21+wfsum1**2
      ecm21_dmc=ecm21_dmc+esum1_dmc(1)**2/wsum1(1)
      efcm21=efcm21+efsum1**2/wfsum1
      ei3cm2=ei3cm2+(wfsum1/wdsum1)**2
      do ifr=1,nforce
        wgcum1(ifr)=wgcum1(ifr)+wgsum1(ifr)
        egcum1(ifr)=egcum1(ifr)+egsum1(ifr)
        wgcm21(ifr)=wgcm21(ifr)+wgsum1(ifr)**2
        egcm21(ifr)=egcm21(ifr)+egsum1(ifr)**2/wgsum1(ifr)
      enddo

c collect block averages
      wsum_dmc=wsum_dmc+wsum1(1)
      wfsum=wfsum+wfsum1
      wdsum=wdsum+wdsumo
      wgdsum=wgdsum+wgdsum1
      esum_dmc=esum_dmc+esum1_dmc(1)
      efsum=efsum+efsum1
      eisum=eisum+wfsum1/wdsum1
      do ifr=1,nforce
        wgsum(ifr)=wgsum(ifr)+wgsum1(ifr)
        egsum(ifr)=egsum(ifr)+egsum1(ifr)
      enddo

c Estimate eigenvalue of G from the energy
      ipmod=mod(ipass,nfprod)
      if(iabs(idmc).eq.1) then
        nfpro=min(nfprod,ipass)
        eigv=(wgsum1(1)/wtgen(ipmod))**(one/nfpro)
       else
        eest=(egcum(1)+egsum(1))/(wgcum(1)+wgsum(1))
        eigv=dexp((etrial-eest)*(taucum(1)+tausum(1))/
     &                          (wgcum(1)+wgsum(1)))
        if(ipr.ge.1) write(ounit,'(''eigv'',9f14.6)') eigv,eest,
     &  egcum(1),egsum(1),wgcum(1),wgsum(1),fprod
      endif

      wdsumo=wsum1(1)
      wgdsumo=wsum1(1)*fprod/ff(mod(ipass+1,nfprod))
      wtgen(ipmod)=wsum1(1)

c zero out step averages
      wfsum1=zero
      wdsum1=zero
      efsum1=zero
      do ifr=1,nforce
        wsum1(ifr)=zero
        wgsum1(ifr)=zero
        esum1_dmc(ifr)=zero
        egsum1(ifr)=zero
      enddo

      return
      end
