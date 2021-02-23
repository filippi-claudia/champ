      subroutine dumper
c MPI version created by Claudia Filippi starting from serial version
c routine to pick up and dump everything needed to restart
c job where it left off
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X,
     &NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20,
     &radmax, delri, NEQSX, MTERMS, MCENT3, NCOEF, MEXCIT
      use dmc_mod, only: MWALK, MFPROD, MFPRD1, MPATH
      use basis, only: zex, betaq, n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz,
     & n4s, n4p, n4fxxx, n4fyyy, n4fzzz, n4fxxy, n4fxxz, n4fyyx, n4fyyz,
     & n4fzzx, n4fzzy, n4fxyz, nsa, npa, ndzra, ndz2a, ndxya, ndxza, ndyza, ndx2a
      use const, only: delta, deltai, etrial, fbias, hb, imetro, ipr, nelec, pi
      use forcest, only: fgcm2, fgcum
      use forcepar, only: deltot, istrech, nforce
      use age, only: iage, ioldest, ioldestmx
      use contrldmc, only: iacc_rej, icross, icuspg, icut_br, icut_e, idiv_v, idmc, ipq,
     &itau_eff, nfprod, rttau, tau, taueff, tautot
      use atom, only: cent, iwctype, ncent, nctype, pecent, znuc

      use iterat, only: iblk, ipass
      use config, only: d2o, peo_dmc, psido_dmc, psijo_dmc, vold_dmc, xold_dmc

      use stats, only: acc, dfus2ac, dfus2un, dr2ac, dr2un, nacc, nbrnch, nodecr, trymove

      use estsum, only: efsum, efsum1, egsum, egsum1, ei1sum, ei2sum, ei3sum, esum1_dmc, esum_dmc,
     &pesum_dmc, r2sum, risum, tausum, tjfsum_dmc, tpbsum_dmc, w_acc_sum, w_acc_sum1, wdsum,
     &wdsum1, wfsum, wfsum1, wg_acc_sum, wg_acc_sum1, wgdsum, wgsum, wgsum1, wsum1, wsum_dmc
      use estcum, only: ecum1_dmc, ecum_dmc, efcum, efcum1, egcum, egcum1, ei1cum, ei2cum,
     &ei3cum, pecum_dmc, r2cum_dmc, ricum, taucum, tjfcum_dmc, tpbcum_dmc, w_acc_cum, w_acc_cum1,
     &wcum1, wcum_dmc, wdcum, wdcum1, wfcum, wfcum1, wg_acc_cum, wg_acc_cum1, wgcum, wgcum1,
     &wgdcum
      use force_dmc, only: itausec, nwprod
      use est2cm, only: ecm21_dmc, ecm2_dmc, efcm2, efcm21, egcm2, egcm21, ei1cm2, ei2cm2,
     &ei3cm2, pecm2_dmc, r2cm2_dmc, ricm2, tjfcm_dmc, tpbcm2_dmc, wcm2, wcm21, wdcm2, wdcm21,
     &wfcm2, wfcm21, wgcm2, wgcm21, wgdcm2
      use derivest, only: derivcm2, derivcum, derivsum, derivtotave_num_old

      use step, only: ekin, ekin2, rprob, suc, trunfb, try
      use mpiconf, only: idtask, nproc, wid, NPROCX

      use denupdn, only: rprobdn, rprobup
      use contr3, only: mode
      use mpiblk, only: iblk_proc
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use mstates_mod, only: MSTATES, MDETCSFX
      use pseudo_mod, only: MPS_L, MPS_QUAD, MPS_GRID, MGAUSS

      use qua, only: nquad, wq, xq, xq0, yq, yq0, zq, zq0

      use branch, only: eest, eigv, eold, ff, fprod, nwalk, pwt, wdsumo, wgdsumo, wt, wtgen,
     &wthist
      use casula, only: i_vpsp, icasula, t_vpsp
      use jacobsave, only: ajacob, ajacold
      use pseudo, only: lpot, nloc, vps, vpso

      use dets, only: cdet, ndet
      use elec, only: ndn, nup
      use coefs, only: coef, nbasis, norb
      use ghostatom, only: newghostype, nghostcent
      use jaspar, only: is, nspin1, nspin2, sspin, sspinn
      implicit real*8(a-h,o-z)











      include 'mpif.h'
C      include 'mpi_qmc.h'
      parameter (zero=0.d0,one=1.d0)
      parameter (small=1.e-6)

      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /velratio/ fratio(MWALK,MFORCE)
      common /jaspar1/ cjas1(MWF),cjas2(MWF)



      dimension irn(4,0:NPROCX),istatus(MPI_STATUS_SIZE)
      dimension irn_tmp(4,0:NPROCX)

      if(mode.eq.'dmc_one_mpi2') then
        call dumper_gpop
      call mpi_barrier(MPI_COMM_WORLD,ierr)
        return
      endif

      if(nforce.gt.1) call strech(xold_dmc,xold_dmc,ajacob,1,0)

      call savern(irn(1,idtask))

      nscounts=4
      call mpi_gather(irn(1,idtask),nscounts,mpi_integer
     &,irn_tmp,nscounts,mpi_integer,0,MPI_COMM_WORLD,ierr)

      if(.not.wid) then
        call mpi_isend(nwalk,1,mpi_integer,0
     &  ,1,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(xold_dmc,3*MELEC*nwalk,mpi_double_precision,0
     &  ,2,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(wt,nwalk,mpi_double_precision,0
     &  ,3,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(ff(0),nfprod,mpi_double_precision,0
     &  ,4,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(fprod,1,mpi_double_precision,0
     &  ,5,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(fratio,MWALK*nforce,mpi_double_precision,0
     &  ,6,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(eigv,1,mpi_double_precision,0
     &  ,7,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(eest,1,mpi_double_precision,0
     &  ,8,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(wdsumo,1,mpi_double_precision,0
     &  ,9,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(iage,nwalk,mpi_integer,0
     &  ,10,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(ioldest,1,mpi_integer,0
     &  ,11,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(ioldestmx,1,mpi_integer,0
     &  ,12,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(xq,nquad,mpi_double_precision,0
     &  ,13,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(yq,nquad,mpi_double_precision,0
     &  ,14,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(zq,nquad,mpi_double_precision,0
     &  ,15,MPI_COMM_WORLD,irequest,ierr)
       else
        open(unit=10,status='unknown',form='unformatted',file='restart_dmc')
        write(10) nproc
        write(10) nwalk
        write(10) (((xold_dmc(ic,i,iw,1),ic=1,3),i=1,nelec),iw=1,nwalk)
        write(10) nfprod,(ff(i),i=0,nfprod),(wt(i),i=1,nwalk),fprod
     &  ,eigv,eest,wdsumo
        write(10) (iage(i),i=1,nwalk),ioldest,ioldestmx
        write(10) nforce,((fratio(iw,ifr),iw=1,nwalk),ifr=1,nforce)
        if(nloc.gt.0)
     &  write(10) nquad,(xq(i),yq(i),zq(i),wq(i),i=1,nquad)
c       if(nforce.gt.1) write(10) nwprod
c    &  ,((pwt(i,j),i=1,nwalk),j=1,nforce)
c    &  ,(((wthist(i,l,j),i=1,nwalk),l=0,nwprod-1),j=1,nforce)
        do 450 id=1,nproc-1
          call mpi_recv(nwalk,1,mpi_integer,id
     &    ,1,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(xold_dmc,3*MELEC*nwalk,mpi_double_precision,id
     &    ,2,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(wt,nwalk,mpi_double_precision,id
     &    ,3,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(ff(0),nfprod,mpi_double_precision,id
     &    ,4,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(fprod,1,mpi_double_precision,id
     &    ,5,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(fratio,MWALK*nforce,mpi_double_precision,id
     &    ,6,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(eigv,1,mpi_double_precision,id
     &    ,7,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(eest,1,mpi_double_precision,id
     &    ,8,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(wdsumo,1,mpi_double_precision,id
     &    ,9,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(iage,nwalk,mpi_integer,id
     &    ,10,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(ioldest,1,mpi_integer,id
     &    ,11,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(ioldestmx,1,mpi_integer,id
     &    ,12,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(xq,nquad,mpi_double_precision,id
     &    ,13,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(yq,nquad,mpi_double_precision,id
     &    ,14,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(zq,nquad,mpi_double_precision,id
     &    ,15,MPI_COMM_WORLD,istatus,ierr)
          write(10) nwalk
          write(10) (((xold_dmc(ic,i,iw,1),ic=1,3),i=1,nelec),iw=1,nwalk)
          write(10) nfprod,(ff(i),i=0,nfprod),(wt(i),i=1,nwalk),fprod
     &    ,eigv,eest,wdsumo
          write(10) (iage(i),i=1,nwalk),ioldest,ioldestmx
          write(10) nforce,((fratio(iw,ifr),iw=1,nwalk),ifr=1,nforce)
          if(nloc.gt.0)
     &    write(10) nquad,(xq(i),yq(i),zq(i),wq(i),i=1,nquad)
c         if(nforce.gt.1) write(10) nwprod
c    &    ,((pwt(i,j),i=1,nwalk),j=1,nforce)
c    &    ,(((wthist(i,l,j),i=1,nwalk),l=0,nwprod-1),j=1,nforce)
  450   continue
      endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(.not.wid) return

      write(10) (wgcum(i),egcum(i),pecum_dmc(i),tpbcum_dmc(i),tjfcum_dmc(i)
     &,wgcm2(i),egcm2(i),pecm2_dmc(i),tpbcm2_dmc(i),tjfcm_dmc(i),taucum(i)
     &,i=1,nforce)
      write(10) ((irn_tmp(i,j),i=1,4),j=0,nproc-1)
      write(10) hb
      write(10) tau,rttau,idmc
      write(10) nelec,nconf,nforce
      write(10) (wtgen(i),i=0,nfprod),wgdsumo
      write(10) wcum_dmc,wfcum,wdcum,wgdcum
     &,wcum1/nproc,wfcum1/nproc,(wgcum1(i)/nproc,i=1,nforce),wdcum1
     &,ecum_dmc,efcum,ecum1_dmc/nproc,efcum1/nproc,(egcum1(i)/nproc,i=1,nforce)
     &,ei1cum,ei2cum,ei3cum,r2cum_dmc,ricum
      write(10) ipass,iblk,iblk_proc
      write(10) wcm2,wfcm2,wdcm2,wgdcm2,wcm21/nproc
     &,wfcm21/nproc,(wgcm21(i)/nproc,i=1,nforce),wdcm21, ecm2_dmc,efcm2
     &,ecm21_dmc/nproc,efcm21/nproc,(egcm21(i)/nproc,i=1,nforce)
     &,ei1cm2,ei2cm2,ei3cm2,r2cm2_dmc,ricm2
      write(10) (fgcum(i),i=1,nforce),(fgcm2(i),i=1,nforce)
     &,((derivcum(k,i),k=1,3),i=1,nforce),(derivcm2(i),i=1,nforce)
     &,(derivtotave_num_old(i),i=1,nforce)
      write(10) (rprob(i)/nproc,rprobup(i),rprobdn(i),i=1,nrad)
      write(10) dfus2ac,dfus2un,dr2ac,dr2un,acc
     &,trymove,nacc,nbrnch,nodecr
      call prop_dump(10)
      call pcm_dump(10)
      call mmpol_dump(10)
      write(10) ((coef(ib,i,1),ib=1,nbasis),i=1,norb)
      write(10) nbasis
      write(10) (zex(ib,1),ib=1,nbasis)
      write(10) nctype,ncent,newghostype,nghostcent,(iwctype(i),i=1,ncent+nghostcent)
      write(10) ((cent(k,ic),k=1,3),ic=1,ncent+nghostcent)
      write(10) pecent
      write(10) (znuc(i),i=1,nctype)
      write(10) (n1s(i),i=1,nctype)
      write(10) (n2s(i),i=1,nctype)
      write(10) ((n2p(ic,i),ic=1,3),i=1,nctype)
      write(10) (n3s(i),i=1,nctype)
      write(10) ((n3p(ic,i),ic=1,3),i=1,nctype)
      write(10) (n3dzr(i),i=1,nctype)
      write(10) (n3dx2(i),i=1,nctype)
      write(10) (n3dxy(i),i=1,nctype)
      write(10) (n3dxz(i),i=1,nctype)
      write(10) (n3dyz(i),i=1,nctype)
      write(10) (n4s(i),i=1,nctype)
      write(10) ((n4p(ic,i),ic=1,3),i=1,nctype)
      write(10) (nsa(i),i=1,nctype)
      write(10) ((npa(ic,i),ic=1,3),i=1,nctype)
      write(10) (ndzra(i),i=1,nctype)
      write(10) (ndx2a(i),i=1,nctype)
      write(10) (ndxya(i),i=1,nctype)
      write(10) (ndxza(i),i=1,nctype)
      write(10) (ndyza(i),i=1,nctype)
      write(10) (cdet(i,1,1),i=1,ndet)
      write(10) ndet,nup,ndn
      write(10) cjas1(1),cjas2(1)
      close (unit=10)
      write(6,'(1x,''successful dump to unit 10'')')

      return
      end
