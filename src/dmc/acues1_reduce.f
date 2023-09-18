      module acues1_reduce_mod
      contains
      subroutine acues1_reduce

      use precision_kinds, only: dp
      use vmc_mod, only: nrad
      use multiple_geo, only: nforce, MFORCE
      use estcum, only: iblk
      use stats, only: acc, nacc, nodecr, trymove
      use estcum, only: ecum1_dmc, efcum1, egcum, egcum1
      use estcum, only: wcum1, wfcum1, wgcum, wgcum1
      use est2cm, only: ecm21_dmc, efcm21, egcm21
      use est2cm, only: wcm21
      use est2cm, only: wfcm21, wgcm21, wgcm2
      use mpiconf, only: nproc, wid
      use control, only: mode
      use contrl_per, only: iperiodic
      use control, only: mode
      use est2cm,  only: ecm21_dmc,efcm21,egcm21,wcm21,wfcm21,wgcm21
      use estcum,  only: ecum1_dmc,efcum1,egcum,egcum1,iblk,wcum1,wfcum1
      use estcum,  only: wgcum,wgcum1
      use mpi
      use mpiconf, only: nproc,wid
      use multiple_geo, only: MFORCE,nforce
      use optci_mod, only: optci_fin
      use optci_reduce_mod, only: optci_reduce
      use optjas_mod, only: optjas_fin
      use optjas_reduce_mod, only: optjas_reduce
      use optorb_f_mod, only: optorb_fin
      use optorb_reduce_mod, only: optorb_reduce
      use optx_jas_ci, only: optx_jas_ci_fin
      use optx_jas_ci_reduce_mod, only: optx_jas_ci_reduce
      use optx_jas_orb, only: optx_jas_orb_fin
      use optx_jas_orb_reduce_mod, only: optx_jas_orb_reduce
      use optx_orb_ci, only: optx_orb_ci_fin
      use optx_orb_ci_reduce_mod, only: optx_orb_ci_reduce
      use precision_kinds, only: dp
      use stats,   only: acc,nacc,nodecr,trymove
      use vmc_mod, only: nrad
      use force_analytic, only: force_analy_fin
      ! use force_analy_reduce_mod, only: force_analy_reduce


      implicit none

      integer :: ierr, ifr, i, nodecr_collect, nacc_collect, id, irequest
      integer, dimension(MPI_STATUS_SIZE) :: istatus
      real(dp) :: e1collect, e21collect, w1collect, w21collect, ef1collect
      real(dp) :: ef21collect, wf1collect, wf21collect, trymove_collect
      real(dp) :: acc_collect, efin
      real(dp) :: rn_eff

      real(dp), dimension(MFORCE) :: eg1collect
      real(dp), dimension(MFORCE) :: eg21collect
      real(dp), dimension(MFORCE) :: wg1collect
      real(dp), dimension(MFORCE) :: wg21collect


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

      do ifr=1,nforce
        egcum1(ifr)=eg1collect(ifr)
        wgcum1(ifr)=wg1collect(ifr)
        egcm21(ifr)=eg21collect(ifr)
        wgcm21(ifr)=wg21collect(ifr)
      enddo

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
      ! call force_analy_reduce

      if(wid) then
        do id=1,nproc-1
          call mpi_isend(egcum1,1,mpi_double_precision,id
     &    ,1,MPI_COMM_WORLD,irequest,ierr)
          call mpi_isend(wgcum1,1,mpi_double_precision,id
     &    ,2,MPI_COMM_WORLD,irequest,ierr)
        enddo
       else
        call mpi_recv(egcum1,1,mpi_double_precision,0
     &  ,1,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(wgcum1,1,mpi_double_precision,0
     &  ,2,MPI_COMM_WORLD,istatus,ierr)
      endif

      efin=egcum(1)/wgcum(1)
      rn_eff=(wgcum(1)**2/wgcm2(1))-1.d0

      call optjas_fin(wgcum1(1),egcum1(1))
      call optci_fin(iblk,wgcum1(1),efin)
      call optorb_fin(wgcum1(1),egcum1(1))
      call optx_jas_ci_fin(wgcum1(1),efin)
      call optx_jas_orb_fin(wgcum1(1),egcum1(1))
      call optx_orb_ci_fin(wgcum1(1),efin)
      if(wid) call force_analy_fin(wgcum(1),rn_eff,efin)

      return
      end
      end module
