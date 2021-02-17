      subroutine dumper
c MPI version created by Claudia Filippi starting from serial version
c routine to pick up and dump everything needed to restart
c job where it left off
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
      implicit real*8(a-h,o-z)










      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'pseudo.h'
      include 'mpif.h'
      include 'mpi_qmc.h'
      parameter (zero=0.d0,one=1.d0)
      parameter (small=1.e-6)

      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /velratio/ fratio(MWALK,MFORCE)
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /ghostatom/ newghostype,nghostcent
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,lpot(MCTYPE),nloc
      common /qua/ xq0(MPS_QUAD),yq0(MPS_QUAD),zq0(MPS_QUAD)
     &,xq(MPS_QUAD),yq(MPS_QUAD),zq(MPS_QUAD),wq(MPS_QUAD),nquad
      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /elec/ nup,ndn
      common /jaspar1/ cjas1(MWF),cjas2(MWF)
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /estcum/ wcum,w_acc_cum,wfcum,wgcum(MFORCE),wg_acc_cum,wdcum,
     &wgdcum, wcum1,w_acc_cum1,wfcum1,wgcum1(MFORCE),wg_acc_cum1,
     &wdcum1, ecum,efcum,egcum(MFORCE),ecum1,efcum1,egcum1(MFORCE),
     &ei1cum,ei2cum,ei3cum, pecum(MFORCE),tpbcum(MFORCE),tjfcum(MFORCE),r2cum,
     &ricum,taucum(MFORCE)
      common /estcm2/ wcm2,wfcm2,wgcm2(MFORCE),wdcm2,wgdcm2, wcm21,
     &wfcm21,wgcm21(MFORCE),wdcm21, ecm2,efcm2,egcm2(MFORCE), ecm21,
     &efcm21,egcm21(MFORCE),ei1cm2,ei2cm2,ei3cm2, pecm2(MFORCE),tpbcm2(MFORCE),
     &tjfcm2(MFORCE),r2cm2,ricm2
      common /derivest/ derivsum(10,MFORCE),derivcum(10,MFORCE)
     &,derivcm2(MFORCE),derivtotave_num_old(MFORCE)
      common /step/try(nrad),suc(nrad),trunfb(nrad),rprob(nrad),
     &ekin(nrad),ekin2(nrad)
      common /denupdn/ rprobup(nrad),rprobdn(nrad)
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk
      common /jacobsave/ ajacob,ajacold(MWALK,MFORCE)
      common /casula/ t_vpsp(MCENT,MPS_QUAD,MELEC),icasula,i_vpsp
      common /force_dmc/ itausec,nwprod

      character*12 mode
      common /contr3/ mode

      common /mpiblk/ iblk_proc
      logical wid
      common /mpiconf/ idtask,nproc,wid

      dimension irn(4,0:nprocx),istatus(MPI_STATUS_SIZE)
      dimension irn_tmp(4,0:nprocx)

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

      write(10) (wgcum(i),egcum(i),pecum(i),tpbcum(i),tjfcum(i)
     &,wgcm2(i),egcm2(i),pecm2(i),tpbcm2(i),tjfcm2(i),taucum(i)
     &,i=1,nforce)
      write(10) ((irn_tmp(i,j),i=1,4),j=0,nproc-1)
      write(10) hb
      write(10) tau,rttau,idmc
      write(10) nelec,nconf,nforce
      write(10) (wtgen(i),i=0,nfprod),wgdsumo
      write(10) wcum,wfcum,wdcum,wgdcum
     &,wcum1/nproc,wfcum1/nproc,(wgcum1(i)/nproc,i=1,nforce),wdcum1
     &,ecum,efcum,ecum1/nproc,efcum1/nproc,(egcum1(i)/nproc,i=1,nforce)
     &,ei1cum,ei2cum,ei3cum,r2cum,ricum
      write(10) ipass,iblk,iblk_proc
      write(10) wcm2,wfcm2,wdcm2,wgdcm2,wcm21/nproc
     &,wfcm21/nproc,(wgcm21(i)/nproc,i=1,nforce),wdcm21, ecm2,efcm2
     &,ecm21/nproc,efcm21/nproc,(egcm21(i)/nproc,i=1,nforce)
     &,ei1cm2,ei2cm2,ei3cm2,r2cm2,ricm2
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
