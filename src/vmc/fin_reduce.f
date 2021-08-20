      subroutine fin_reduce
c MPI version written by Claudia Filippi

      use vmc_mod, only: nrad
      use csfs, only: nstates
      use mstates_mod, only: MSTATES
      use est2cm, only: ecm21
      use estcum, only: ecum1, iblk
      use estsig, only: ecm21s, ecum1s
      use forcepar, only: nforce
      use forcewt, only: wcum
      use mpiconf, only: nproc, wid
      use step, only: rprob, suc, try
      !use contrl, only: nstep
      use control_vmc, only: vmc_nstep
      use method_opt, only: method
      use mpi
      use custom_broadcast,   only: bcast

      use precision_kinds, only: dp
      implicit none

      integer :: i, id, ierr, istate
      integer, dimension(MPI_STATUS_SIZE) :: istatus
      real(dp) :: dble, efin, passes
      real(dp), dimension(nrad) :: rprobt
      real(dp), dimension(nrad) :: tryt
      real(dp), dimension(nrad) :: suct
      real(dp), dimension(MSTATES) :: collect




      call mpi_reduce(ecum1,collect,nstates,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      do 10 istate=1,nstates
  10    ecum1(istate)=collect(istate)

      call mpi_reduce(ecm21,collect,nstates,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      do 20 istate=1,nstates
  20    ecm21(istate)=collect(istate)

      call mpi_reduce(ecum1s,collect,nstates,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      do 30 istate=1,nstates
  30    ecum1s(istate)=collect(istate)

      call mpi_reduce(ecm21s,collect,nstates,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      do 40 istate=1,nstates
  40    ecm21s(istate)=collect(istate)

      call mpi_reduce(rprob,rprobt,nrad,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(suc,suct,nrad,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(try,tryt,nrad,mpi_double_precision
     &,mpi_sum,0,MPI_COMM_WORLD,ierr)

      do 50 i=1,nrad
        rprob(i)=rprobt(i)
        suc(i)=suct(i)
   50   try(i)=tryt(i)

      call optjas_reduce

      call optorb_reduce

      call optci_reduce

      if(method.ne.'sr_n'.and.method.ne.'lin_d') then
        call optx_jas_orb_reduce

        call optx_jas_ci_reduce

        call optx_orb_ci_reduce
      endif

      ! broadcast the ecum1 and wcum data using custom-made broadcaster
      call bcast(ecum1)
      call bcast(wcum)


      passes=dble(iblk)*dble(vmc_nstep)
      efin=ecum1(1)/wcum(1,1)

      call optjas_fin(wcum(1,1),ecum1)

      call optci_fin(iblk,wcum(1,1),efin)

      call optorb_fin(wcum(1,1),ecum1)

      if(method.ne.'sr_n'.and.method.ne.'lin_d') then
        call optx_jas_ci_fin(wcum(1,1),efin)

        call optx_jas_orb_fin(wcum(1,1),ecum1)

        call optx_orb_ci_fin(wcum(1,1),efin)
      endif

c     call efficiency_prt(passes)

c reduce pcm properties
      call pcm_reduce

      if(wid) call pcm_fin(wcum(1,1),iblk)

c reduce mmpol properties
      call mmpol_reduce

      if(wid) call mmpol_fin(wcum(1,1),iblk)

c reduce analytical forces
      call force_analy_reduce

      if(wid) call force_analy_fin(wcum(1,1),iblk,efin)

      return
      end
