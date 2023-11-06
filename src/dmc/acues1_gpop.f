      module acues1_gpop_mod
      contains
      subroutine acues1_gpop
c MPI version created by Claudia Filippi starting from serial version
c routine to accumulate estimators for energy etc.

      use const, only: etrial
      use contrl_file, only: errunit,ounit
      use control, only: ipr
      use multiple_geo, only: nforce, MFORCE
      use contrldmc, only: idmc
      use contrldmc, only: nfprod
      use estcum, only: ipass
      use estsum, only: efsum, efsum1, egsum, egsum1, esum1_dmc, esum_dmc
      use estsum, only: tausum
      use estsum, only: wfsum, wfsum1, wgsum, wgsum1, wsum1, wsum_dmc
      use estcum, only: ecum1_dmc, efcum1, egcum, egcum1
      use estcum, only: taucum, wcum1, wfcum1, wgcum, wgcum1
      use est2cm, only: ecm21_dmc, efcm21, egcm21, wcm21
      use est2cm, only: wfcm21, wgcm21
      use mpiconf, only: wid
      use branch, only: eest, eigv, ff, fprod, wdsumo, wgdsumo, wtgen
      use mpi
      use mpiconf, only: wid
      use multiple_geo, only: MFORCE,nforce
      use precision_kinds, only: dp
      use redistribute_mod, only: redistribute
      implicit none

      integer :: ierr, ifr, ipmod, mod, iabs, nfpro
      real(dp) :: ecollect
      real(dp), dimension(MFORCE) :: egcollect
      real(dp), dimension(MFORCE) :: wgcollect
      real(dp) :: wcollect, efcollect, wfcollect, taublock, accavn

      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0


c statistical fluctuations without blocking

      call mpi_reduce(esum1_dmc(1),ecollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wsum1(1),wcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(efsum1,efcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wfsum1,wfcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(wgsum1,wgcollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(egsum1,egcollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tausum(1),taublock,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)


      if(.not.wid) goto 23

      wsum1(1)=wcollect
      esum1_dmc(1)=ecollect
      efsum1=efcollect
      wfsum1=wfcollect
      do ifr=1,nforce
        wgsum1(ifr)=wgcollect(ifr)
        egsum1(ifr)=egcollect(ifr)
      enddo

      wcum1=wcum1+wsum1(1)
      wfcum1=wfcum1+wfsum1
      ecum1_dmc=ecum1_dmc+esum1_dmc(1)
      efcum1=efcum1+efsum1

      wcm21=wcm21+wsum1(1)**2
      wfcm21=wfcm21+wfsum1**2
      ecm21_dmc=ecm21_dmc+esum1_dmc(1)**2/wsum1(1)
      efcm21=efcm21+efsum1**2/wfsum1

      do ifr=1,nforce
        wgcum1(ifr)=wgcum1(ifr)+wgsum1(ifr)
        egcum1(ifr)=egcum1(ifr)+egsum1(ifr)
        wgcm21(ifr)=wgcm21(ifr)+wgsum1(ifr)**2
        egcm21(ifr)=egcm21(ifr)+egsum1(ifr)**2/wgsum1(ifr)
      enddo

c sum1 block averages
      wsum_dmc=wsum_dmc+wsum1(1)
      wfsum=wfsum+wfsum1
      esum_dmc=esum_dmc+esum1_dmc(1)
      efsum=efsum+efsum1

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
        eigv=dexp((etrial-eest)*(taucum(1)+taublock)/
     &                          (wgcum(1)+wgsum(1)))
        if(ipr.ge.1) write(ounit,'(''eigv'',9f14.6)') eigv,eest,accavn,
     &  egcum(1),egsum(1),wgcum(1),wgsum(1),fprod
      endif

      wdsumo=wsum1(1)
      wgdsumo=wsum1(1)*fprod/ff(mod(ipass+1,nfprod))
      wtgen(ipmod)=wsum1(1)

   23 call mpi_bcast(eigv,1,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(eest,1,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(wdsumo,1,mpi_double_precision,0,MPI_COMM_WORLD
     &,ierr)


      call redistribute


c zero out step averages
      wfsum1=zero
      efsum1=zero
      do ifr=1,nforce
        wsum1(ifr)=zero
        wgsum1(ifr)=zero
        esum1_dmc(ifr)=zero
        egsum1(ifr)=zero
      enddo

      return
      end
      end module
