      module dumper_mod
      contains
      subroutine dumper
c MPI version created by Claudia Filippi starting from serial version
c routine to pick up and dump everything needed to restart
c job where it left off

      use age,     only: iage,ioldest,ioldestmx
      use basis,   only: ns, np, nd, nf, ng, zex
      use branch,  only: eest,eigv,ff,fprod,nwalk,wdsumo,wgdsumo,wt
      use branch,  only: wtgen
      use coefs,   only: nbasis
      use config,  only: xold_dmc
      use constants, only: hb
      use contrl_file, only: ounit
      use contrldmc, only: idmc,nfprod,rttau,tau
      use control, only: mode
      use control_dmc, only: dmc_nconf
      use dmc_mod, only: MWALK
      use dumper_gpop_mod, only: dumper_gpop
      use est2cm,  only: ecm21_dmc,ecm2_dmc,efcm2,efcm21,egcm2,egcm21
      use est2cm,  only: pecm2_dmc
      use est2cm,  only: tpbcm2_dmc,wcm2,wcm21
      use est2cm,  only: wfcm2,wfcm21,wgcm2,wgcm21
      use estcum,  only: ecum1_dmc,ecum_dmc,efcum,efcum1,egcum,egcum1
      use estcum,  only: iblk,ipass,pecum_dmc
      use estcum,  only: taucum,tpbcum_dmc
      use estcum,  only: wcum1,wcum_dmc,wfcum,wfcum1,wgcum,wgcum1
      use force_analytic, only: force_analy_dump
      use jacobsave, only: ajacob
      use mmpol,   only: mmpol_dump
      use mpi
      use mpiblk,  only: iblk_proc
      use mpiconf, only: idtask,nproc,wid
      use multiple_geo, only: fgcm2,fgcum,nforce,pecent
      use pcm_mod, only: pcm_dump
      use precision_kinds, only: dp
      use properties_mod, only: prop_dump
      use pseudo,  only: nloc
      use qua,     only: nquad,wq,xq,yq,zq
      use random_mod, only: savern
      use slater,  only: cdet,coef,ndet,norb
      use stats,   only: acc,dfus2ac,dfus2un,nacc,nbrnch
      use stats,   only: nodecr,trymove
      use strech_mod, only: strech
      use system,  only: cent,iwctype,ncent,nctype,ndn,nelec,newghostype
      use system,  only: nghostcent,nup,znuc
      use velratio, only: fratio
      use vmc_mod, only: nrad

      implicit none

      integer :: i, ib, ic, id, ierr
      integer :: ifr, irequest, iw, j
      integer :: k
      integer, dimension(8, 0:nproc) :: irn
      integer, dimension(MPI_STATUS_SIZE) :: istatus

      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: small = 1.e-6

      if(mode.eq.'dmc_one_mpi2') then
        call dumper_gpop
c       call mpi_barrier(MPI_COMM_WORLD,ierr)
        return
      endif

      if(nforce.gt.1) call strech(xold_dmc,xold_dmc,ajacob,1,0)

      call savern(irn(1,idtask))

      if(.not.wid) then
        call mpi_send(nwalk,1,mpi_integer,0
     &  ,1,MPI_COMM_WORLD,irequest,ierr)
        call mpi_send(xold_dmc,3*nelec*nwalk,mpi_double_precision,0
     &  ,2,MPI_COMM_WORLD,irequest,ierr)
        call mpi_send(wt,nwalk,mpi_double_precision,0
     &  ,3,MPI_COMM_WORLD,irequest,ierr)
        call mpi_send(ff(0),nfprod,mpi_double_precision,0
     &  ,4,MPI_COMM_WORLD,irequest,ierr)
        call mpi_send(fprod,1,mpi_double_precision,0
     &  ,5,MPI_COMM_WORLD,irequest,ierr)
        call mpi_send(fratio,MWALK*nforce,mpi_double_precision,0
     &  ,6,MPI_COMM_WORLD,irequest,ierr)
        call mpi_send(eigv,1,mpi_double_precision,0
     &  ,7,MPI_COMM_WORLD,irequest,ierr)
        call mpi_send(eest,1,mpi_double_precision,0
     &  ,8,MPI_COMM_WORLD,irequest,ierr)
        call mpi_send(wdsumo,1,mpi_double_precision,0
     &  ,9,MPI_COMM_WORLD,irequest,ierr)
        call mpi_send(iage,nwalk,mpi_integer,0
     &  ,10,MPI_COMM_WORLD,irequest,ierr)
        call mpi_send(ioldest,1,mpi_integer,0
     &  ,11,MPI_COMM_WORLD,irequest,ierr)
        call mpi_send(ioldestmx,1,mpi_integer,0
     &  ,12,MPI_COMM_WORLD,irequest,ierr)
        call mpi_send(xq,nquad,mpi_double_precision,0
     &  ,13,MPI_COMM_WORLD,irequest,ierr)
        call mpi_send(yq,nquad,mpi_double_precision,0
     &  ,14,MPI_COMM_WORLD,irequest,ierr)
        call mpi_send(zq,nquad,mpi_double_precision,0
     &  ,15,MPI_COMM_WORLD,irequest,ierr)
        call mpi_send(irn(:,idtask), 8, mpi_integer, 0, 16, MPI_COMM_WORLD, ierr)
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
        do id=1,nproc-1
          call mpi_recv(nwalk,1,mpi_integer,id
     &    ,1,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(xold_dmc,3*nelec*nwalk,mpi_double_precision,id
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
          call mpi_recv(irn(:, id), 8, mpi_integer, id,16,MPI_COMM_WORLD,istatus,ierr)
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
        enddo
      endif
c       call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(.not.wid) return

      write(10) (wgcum(i),egcum(i),pecum_dmc(i),tpbcum_dmc(i)
     &,wgcm2(i),egcm2(i),pecm2_dmc(i),tpbcm2_dmc(i),taucum(i)
     &,i=1,nforce)
      write(10) ((irn(i,j),i=1,8),j=0,nproc-1)
      write(10) hb
      write(10) tau,rttau,idmc
      write(10) nelec,dmc_nconf,nforce
      write(10) (wtgen(i),i=0,nfprod),wgdsumo
      write(10) wcum_dmc,wfcum
     &,wcum1/nproc,wfcum1/nproc,(wgcum1(i)/nproc,i=1,nforce)
     &,ecum_dmc,efcum,ecum1_dmc/nproc,efcum1/nproc,(egcum1(i)/nproc,i=1,nforce)
      write(10) ipass,iblk,iblk_proc
      write(10) wcm2,wfcm2,wcm21/nproc
     &,wfcm21/nproc,(wgcm21(i)/nproc,i=1,nforce),ecm2_dmc,efcm2
     &,ecm21_dmc/nproc,efcm21/nproc,(egcm21(i)/nproc,i=1,nforce)
      write(10) (fgcum(i),i=1,nforce),(fgcm2(i),i=1,nforce)
      write(10) dfus2ac,dfus2un,acc
     &,trymove,nacc,nbrnch,nodecr
      call prop_dump(10)
      call pcm_dump(10)
      call mmpol_dump(10)
      call force_analy_dump(10)
      write(10) ((coef(ib,i,1),ib=1,nbasis),i=1,norb)
      write(10) nbasis
      write(10) (zex(ib,1),ib=1,nbasis)
      write(10) nctype,ncent,newghostype,nghostcent,(iwctype(i),i=1,ncent+nghostcent)
      write(10) ((cent(k,ic),k=1,3),ic=1,ncent+nghostcent)
      write(10) pecent
      write(10) (znuc(i),i=1,nctype)
      write(10) (ns(i),i=1,nctype)
      write(10) (np(i),i=1,nctype)
      write(10) (nd(i),i=1,nctype)
      write(10) (nf(i),i=1,nctype)
      write(10) (ng(i),i=1,nctype)
      write(10) (cdet(i,1,1),i=1,ndet)
      write(10) ndet,nup,ndn
      close (unit=10)
      write(ounit,'(1x,''successful dump to unit 10'')')

      return
      end
      end module
