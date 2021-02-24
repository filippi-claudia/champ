      subroutine acues1_reduce

      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X,
     &NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20,
     &radmax, delri, NEQSX, MTERMS, MCENT3, NCOEF, MEXCIT
      use dmc_mod, only: MWALK, MFPROD, MFPRD1, MPATH
      use forcepar, only: deltot, istrech, nforce
      use estcum, only: iblk, ipass
      use stats, only: acc, dfus2ac, dfus2un, dr2ac, dr2un, nacc, nbrnch, nodecr, trymove
      use estcum, only: ecum1_dmc, ecum_dmc, efcum, efcum1, egcum, egcum1, ei1cum, ei2cum,
     &ei3cum, pecum_dmc, r2cum_dmc, ricum, taucum, tjfcum_dmc, tpbcum_dmc, w_acc_cum, w_acc_cum1,
     &wcum1, wcum_dmc, wdcum, wdcum1, wfcum, wfcum1, wg_acc_cum, wg_acc_cum1, wgcum, wgcum1,
     &wgdcum
      use est2cm, only: ecm21_dmc, ecm2_dmc, efcm2, efcm21, egcm2, egcm21, ei1cm2, ei2cm2,
     &ei3cm2, pecm2_dmc, r2cm2_dmc, ricm2, tjfcm_dmc, tpbcm2_dmc, wcm2, wcm21, wdcm2, wdcm21,
     &wfcm2, wfcm21, wgcm2, wgcm21, wgdcm2
      use step, only: ekin, ekin2, rprob, suc, trunfb, try
      use mpiconf, only: idtask, nproc, wid, NPROCX
      use contr3, only: mode
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF

      implicit real*8(a-h,o-z)


      include 'mpif.h'

      dimension eg1collect(MFORCE),eg21collect(MFORCE),wg1collect(MFORCE)
     &,wg21collect(MFORCE),taucollect(MFORCE),rprobcollect(nrad)

      if(mode.eq.'dmc_one_mpi2') return

      call mpi_reduce(ecum1_dmc,e1collect,1,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(ecm21_dmc,e21collect,1,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wcum1,w1collect,1,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wcm21,w21collect,1,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)

      ecum1_dmc=e1collect
      wcum1=w1collect
      ecm21_dmc=e21collect
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
