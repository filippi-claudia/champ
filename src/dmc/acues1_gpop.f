      subroutine acues1_gpop
c MPI version created by Claudia Filippi starting from serial version
c routine to accumulate estimators for energy etc.

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
      implicit real*8(a-h,o-z)









      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'
      include 'mpif.h'
      parameter (zero=0.d0,one=1.d0)

      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /estcm2/ wcm2,wfcm2,wgcm2(MFORCE),wdcm2,wgdcm2, wcm21,
     &wfcm21,wgcm21(MFORCE),wdcm21, ecm2,efcm2,egcm2(MFORCE), ecm21,
     &efcm21,egcm21(MFORCE),ei1cm2,ei2cm2,ei3cm2, pecm2(MFORCE),tpbcm2(MFORCE),
     &tjfcm2(MFORCE),r2cm2,ricm2
      common /derivest/ derivsum(10,MFORCE),derivcum(10,MFORCE)
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk

      logical wid
      common /mpiconf/ idtask,nproc,wid

      dimension egcollect(MFORCE),wgcollect(MFORCE),pecollect(MFORCE),
     +tpbcollect(MFORCE),tjfcollect(MFORCE),taucollect(MFORCE),
     +derivcollect(10,MFORCE)

c statistical fluctuations without blocking
      wdsum1=wdsumo
      wgdsum1=wgdsumo

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

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(.not.wid) goto 23

      wsum1(1)=wcollect
      esum1_dmc(1)=ecollect
      efsum1=efcollect
      wfsum1=wfcollect
      do 20 ifr=1,nforce
        wgsum1(ifr)=wgcollect(ifr)
   20   egsum1(ifr)=egcollect(ifr)

      wcum1=wcum1+wsum1(1)
      wfcum1=wfcum1+wfsum1
      ecum1_dmc=ecum1_dmc+esum1_dmc(1)
      efcum1=efcum1+efsum1
      ei3cum=ei3cum+wfsum1/wdsum1

      wcm21=wcm21+wsum1(1)**2
      wfcm21=wfcm21+wfsum1**2
      ecm21=ecm21+esum1_dmc(1)**2/wsum1(1)
      efcm21=efcm21+efsum1**2/wfsum1
      ei3cm2=ei3cm2+(wfsum1/wdsum1)**2

      do 21 ifr=1,nforce
        wgcum1(ifr)=wgcum1(ifr)+wgsum1(ifr)
        egcum1(ifr)=egcum1(ifr)+egsum1(ifr)
        wgcm21(ifr)=wgcm21(ifr)+wgsum1(ifr)**2
   21   egcm21(ifr)=egcm21(ifr)+egsum1(ifr)**2/wgsum1(ifr)

c sum1 block averages
      wsum_dmc=wsum_dmc+wsum1(1)
      wfsum=wfsum+wfsum1
      wdsum=wdsum+wdsum1
      wgdsum=wgdsum+wgdsum1
      esum_dmc=esum_dmc+esum1_dmc(1)
      efsum=efsum+efsum1
      eisum=eisum+wfsum1/wdsum1
      do 22 ifr=1,nforce
        wgsum(ifr)=wgsum(ifr)+wgsum1(ifr)
   22   egsum(ifr)=egsum(ifr)+egsum1(ifr)

c Estimate eigenvalue of G from the energy
      ipmod=mod(ipass,nfprod)
      if(iabs(idmc).eq.1) then
        nfpro=min(nfprod,ipass)
        eigv=(wgsum1(1)/wtgen(ipmod))**(one/nfpro)
       else
        eest=(egcum(1)+egsum(1))/(wgcum(1)+wgsum(1))
        eigv=dexp((etrial-eest)*(taucum(1)+taublock)/
     &                          (wgcum(1)+wgsum(1)))
        if(ipr.ge.1) write(6,'(''eigv'',9f14.6)') eigv,eest,accavn,
     &  egcum(1),egsum(1),wgcum(1),wgsum(1),fprod
      endif

      wdsumo=wsum1(1)
      wgdsumo=wsum1(1)*fprod/ff(mod(ipass+1,nfprod))
      wtgen(ipmod)=wsum1(1)

   23 call mpi_bcast(eigv,1,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(eest,1,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(wdsumo,1,mpi_double_precision,0,MPI_COMM_WORLD
     &,ierr)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      call redistribute

      call mpi_barrier(MPI_COMM_WORLD,ierr)

c zero out step averages
      wfsum1=zero
      wdsum1=zero
      efsum1=zero
      do 24 ifr=1,nforce
        wsum1(ifr)=zero
        wgsum1(ifr)=zero
        esum1_dmc(ifr)=zero
   24   egsum1(ifr)=zero

      return
      end
