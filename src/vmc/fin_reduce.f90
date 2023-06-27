module fin_reduce_mod
contains
      subroutine fin_reduce
! MPI version written by Claudia Filippi

      use control_vmc, only: vmc_nstep
      use csfs,    only: nstates
      use custom_broadcast, only: bcast
      use est2cm,  only: ecm21
      use estcum,  only: ecum1,iblk
      use estsig,  only: ecm21s,ecum1s
      use force_analy_reduce_mod, only: force_analy_reduce
      use force_analytic, only: force_analy_fin
      use forcewt, only: wcum
      use mmpol_reduce_mod, only: mmpol_reduce
      use mmpol_vmc, only: mmpol_fin
      use mpi
      use mpiconf, only: nproc,wid
      use mstates_mod, only: MSTATES
      use multiple_geo, only: nforce
      use optci_mod, only: optci_fin
      use optci_reduce_mod, only: optci_reduce
      use optjas_mod, only: optjas_fin
      use optjas_reduce_mod, only: optjas_reduce
      use optorb_f_mod, only: optorb_fin
      use optorb_reduce_mod, only: optorb_reduce
      use optwf_control, only: method
      use optx_jas_ci, only: optx_jas_ci_fin
      use optx_jas_ci_reduce_mod, only: optx_jas_ci_reduce
      use optx_jas_orb, only: optx_jas_orb_fin
      use optx_jas_orb_reduce_mod, only: optx_jas_orb_reduce
      use optx_orb_ci, only: optx_orb_ci_fin
      use optx_orb_ci_reduce_mod, only: optx_orb_ci_reduce
      use pcm_reduce_mod, only: pcm_reduce
      use pcm_vmc, only: pcm_fin
      use precision_kinds, only: dp
      use step,    only: rprob,suc,try
      use vmc_mod, only: nrad
!use contrl, only: nstep

      implicit none

      integer :: i, id, ierr, istate
      integer, dimension(MPI_STATUS_SIZE) :: istatus
      real(dp) :: dble, efin, passes
      real(dp), dimension(nrad) :: rprobt
      real(dp), dimension(nrad) :: tryt
      real(dp), dimension(nrad) :: suct
      real(dp), dimension(MSTATES) :: collect




      call mpi_reduce(ecum1,collect,nstates,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)
      do istate=1,nstates
        ecum1(istate)=collect(istate)
      enddo

      call mpi_reduce(ecm21,collect,nstates,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)
      do istate=1,nstates
        ecm21(istate)=collect(istate)
      enddo

      call mpi_reduce(ecum1s,collect,nstates,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)
      do istate=1,nstates
        ecum1s(istate)=collect(istate)
      enddo

      call mpi_reduce(ecm21s,collect,nstates,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)
      do istate=1,nstates
        ecm21s(istate)=collect(istate)
      enddo

      call mpi_reduce(rprob,rprobt,nrad,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(suc,suct,nrad,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(try,tryt,nrad,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)

      do i=1,nrad
        rprob(i)=rprobt(i)
        suc(i)=suct(i)
        try(i)=tryt(i)
      enddo

      call optjas_reduce

      call optorb_reduce

      call optci_reduce

      if(method.ne.'sr_n'.and.method.ne.'lin_d') then
        call optx_jas_orb_reduce

        call optx_jas_ci_reduce

        call optx_orb_ci_reduce
      endif

      call bcast(ecum1)
      call bcast(wcum)
!       if(wid) then
!         do 60 id=1,nproc-1
!           call mpi_send(ecum1,nstates,mpi_double_precision,id
!      &    ,1,MPI_COMM_WORLD,ierr)
! c    &    ,1,MPI_COMM_WORLD,irequest,ierr)
!    60     call mpi_send(wcum,nstates*nforce,mpi_double_precision,id
!      &    ,2,MPI_COMM_WORLD,ierr)
! c    &    ,2,MPI_COMM_WORLD,irequest,ierr)
!        else
!         call mpi_recv(ecum1,nstates,mpi_double_precision,0
!      &  ,1,MPI_COMM_WORLD,istatus,ierr)
!         call mpi_recv(wcum,nstates*nforce,mpi_double_precision,0
!      &  ,2,MPI_COMM_WORLD,istatus,ierr)
!       endif

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

!     call efficiency_prt(passes)

! reduce pcm properties
      call pcm_reduce

      if(wid) call pcm_fin(wcum(1,1),iblk)

! reduce mmpol properties
      call mmpol_reduce

      if(wid) call mmpol_fin(wcum(1,1),iblk)

! reduce analytical forces
      call force_analy_reduce

      if(wid) call force_analy_fin(wcum(1,1),iblk,efin)

      return
      end
end module
