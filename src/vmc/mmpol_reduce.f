      module mmpol_reduce_mod
      contains
      subroutine mmpol_reduce

      use mmpol_averages, only: cmmpol_cm2,cmmpol_cum,dmmpol_cm2
      use mmpol_averages, only: dmmpol_cum,eek1_cm2,eek1_cum,eek2_cm2
      use mmpol_averages, only: eek2_cum,eek3_cm2,eek3_cum
      use mmpol_cntrl, only: immpol
      use mmpol_mod, only: MCHMM
      use mmpol_parms, only: nchmm
      use mpi
      use mpiconf, only: wid
      use precision_kinds, only: dp

      implicit none

      integer :: i, ierr
      real(dp) :: cmmpol_collect, cpcm2_collect, dmmpol_collect, dpcm2_collect
      real(dp), dimension(MCHMM) :: eek1_collect
      real(dp), dimension(MCHMM) :: eek1cm2_collect
      real(dp), dimension(MCHMM) :: eek2_collect
      real(dp), dimension(MCHMM) :: eek2cm2_collect
      real(dp), dimension(MCHMM) :: eek3_collect
      real(dp), dimension(MCHMM) :: eek3cm2_collect

 
      if(immpol.eq.0) return

      call mpi_reduce(dmmpol_cum,dmmpol_collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(dmmpol_cm2,dpcm2_collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(cmmpol_cum,cmmpol_collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(cmmpol_cm2,cpcm2_collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(eek1_cum,eek1_collect,nchmm
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(eek2_cum,eek2_collect,nchmm
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(eek3_cum,eek3_collect,nchmm
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)


      if(wid) then
        do i=1,nchmm
        eek1_cum(i)=eek1_collect(i)
        eek2_cum(i)=eek2_collect(i)
        eek3_cum(i)=eek3_collect(i)
        enddo
       else
        do  i=1,nchmm
        eek1_cum(i)=0
        eek2_cum(i)=0
        eek3_cum(i)=0
        enddo
      endif

      call mpi_reduce(eek1_cm2,eek1cm2_collect,nchmm
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(eek2_cm2,eek2cm2_collect,nchmm
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(eek3_cm2,eek3cm2_collect,nchmm
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      if(wid) then
        do  i=1,nchmm
        eek1_cm2(i)=eek1cm2_collect(i)
        eek2_cm2(i)=eek2cm2_collect(i)
        eek3_cm2(i)=eek3cm2_collect(i)
        enddo
       else
        do i=1,nchmm
        eek1_cm2(i)=0.0d0
        eek2_cm2(i)=0.0d0
        eek3_cm2(i)=0.0d0
        enddo
      endif

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(wid) then
        dmmpol_cum=dmmpol_collect
        dmmpol_cm2=dpcm2_collect
        cmmpol_cum=cmmpol_collect
        cmmpol_cm2=cpcm2_collect
       else
        dmmpol_cum=0.0d0         
        dmmpol_cm2=0.0d0         
        cmmpol_cum=0.0d0         
        cmmpol_cm2=0.0d0         
      endif

      return
      end

      end module
