      subroutine acues1_reduce

      use forcepar, only: deltot, istrech, nforce
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'
      include 'mpif.h'

      common /iterat/ ipass,iblk
      common /stats/ dfus2ac,dfus2un,dr2ac,dr2un,acc,trymove,nacc,
     &nbrnch,nodecr
      common /estcum/ wcum,w_acc_cum,wfcum,wgcum(MFORCE),wg_acc_cum,wdcum,
     &wgdcum, wcum1,w_acc_cum1,wfcum1,wgcum1(MFORCE),wg_acc_cum1,
     &wdcum1, ecum,efcum,egcum(MFORCE),ecum1,efcum1,egcum1(MFORCE),
     &ei1cum,ei2cum,ei3cum, pecum(MFORCE),tpbcum(MFORCE),tjfcum(MFORCE),r2cum,
     &ricum,taucum(MFORCE)
      common /estcm2/ wcm2,wfcm2,wgcm2(MFORCE),wdcm2,wgdcm2, wcm21,
     &wfcm21,wgcm21(MFORCE),wdcm21, ecm2,efcm2,egcm2(MFORCE), ecm21,
     &efcm21,egcm21(MFORCE),ei1cm2,ei2cm2,ei3cm2, pecm2(MFORCE),tpbcm2(MFORCE),
     &tjfcm2(MFORCE),r2cm2,ricm2
      common /step/try(nrad),suc(nrad),trunfb(nrad),rprob(nrad),
     &ekin(nrad),ekin2(nrad)

      character*12 mode
      common /contr3/ mode

      logical wid
      common /mpiconf/ idtask,nproc,wid

      dimension eg1collect(MFORCE),eg21collect(MFORCE),wg1collect(MFORCE)
     &,wg21collect(MFORCE),taucollect(MFORCE),rprobcollect(nrad)

      if(mode.eq.'dmc_one_mpi2') return

      call mpi_reduce(ecum1,e1collect,1,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(ecm21,e21collect,1,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wcum1,w1collect,1,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wcm21,w21collect,1,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)

      ecum1=e1collect
      wcum1=w1collect
      ecm21=e21collect
      wcm21=w21collect

      call mpi_reduce(efcum1,ef1collect,1,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(efcm21,ef21collect,1,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wfcum1,wf1collect,1,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wfcm21,wf21collect,1,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)

      efcum1=ef1collect
      wfcum1=wf1collect
      efcm21=ef21collect
      wfcm21=wf21collect

      call mpi_reduce(egcum1,eg1collect,MFORCE,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(egcm21,eg21collect,MFORCE,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wgcum1,wg1collect,MFORCE,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wgcm21,wg21collect,MFORCE,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)

      do 1 ifr=1,nforce
        egcum1(ifr)=eg1collect(ifr)
        wgcum1(ifr)=wg1collect(ifr)
        egcm21(ifr)=eg21collect(ifr)
    1   wgcm21(ifr)=wg21collect(ifr)

c Collect radial charge density for atoms
      if(iperiodic.eq.0) then
        call mpi_reduce(rprob,rprobcollect,nrad,mpi_double_precision
     &  ,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do 2 i=1,nrad
    2     rprob(i)=rprobcollect(i)
      endif

      call mpi_reduce(nodecr,nodecr_collect,1,mpi_integer,mpi_sum,0,
     &MPI_COMM_WORLD,ierr)
      call mpi_reduce(trymove,trymove_collect,1,mpi_double_precision,mpi_sum,0,
     &MPI_COMM_WORLD,ierr)
      call mpi_reduce(acc,acc_collect,1,mpi_double_precision,mpi_sum,0,
     &MPI_COMM_WORLD,ierr)
      call mpi_reduce(nacc,nacc_collect,1,mpi_integer,mpi_sum,0,
     &MPI_COMM_WORLD,ierr)
      nodecr=nodecr_collect
      trymove=trymove_collect
      acc=acc_collect
      nacc=nacc_collect

      call optjas_reduce
      call optorb_reduce
      call optci_reduce
      call optx_jas_orb_reduce
      call optx_jas_ci_reduce
      call optx_orb_ci_reduce

      if(wid) then
        do 60 id=1,nproc-1
          call mpi_isend(egcum1,1,mpi_double_precision,id
     &    ,1,MPI_COMM_WORLD,irequest,ierr)
   60     call mpi_isend(wgcum1,1,mpi_double_precision,id
     &    ,2,MPI_COMM_WORLD,irequest,ierr)
       else
        call mpi_recv(egcum1,1,mpi_double_precision,0
     &  ,1,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(wgcum1,1,mpi_double_precision,0
     &  ,2,MPI_COMM_WORLD,istatus,ierr)
      endif

c     efin=egcum1(1)/wgcum1(1)
      efin=egcum(1)/wgcum(1)

      call optjas_fin(wgcum1(1),egcum1(1))
      call optci_fin(iblk,wgcum1(1),efin)
      call optorb_fin(wgcum1(1),egcum1(1))
      call optx_jas_ci_fin(wgcum1(1),efin)
      call optx_jas_orb_fin(wgcum1(1),egcum1(1))
      call optx_orb_ci_fin(wgcum1(1),efin)

      return
      end
