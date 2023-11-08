module pcm_reduce_mod
      use error,   only: fatal_error
contains
      subroutine pcm_reduce

      use contrl_file, only: ounit
      use mpi
      use mpiconf, only: wid
      use pcm,     only: MCHS
      use pcm_averages, only: enfpcm_cm2,enfpcm_cum,qopcm_cm2,qopcm_cum
      use pcm_averages, only: spcmcm2,spcmcum,vpcmcm2,vpcmcum
      use pcm_cntrl, only: ipcm
      use pcm_parms, only: nchs
      use precision_kinds, only: dp
      implicit none

      integer :: i, ierr
      real(dp) :: qopcm2collect, qopcmcollect, spcm2collect, spcmcollect, vpcm2collect
      real(dp) :: vpcmcollect
      real(dp), dimension(MCHS) :: collect

      if(ipcm.eq.0) return

      call mpi_reduce(spcmcum,spcmcollect,1 &
      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(spcmcm2,spcm2collect,1 &
      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(vpcmcum,vpcmcollect,1 &
      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(vpcmcm2,vpcm2collect,1 &
      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(qopcm_cum,qopcmcollect,1 &
      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(qopcm_cm2,qopcm2collect,1 &
      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(enfpcm_cum,collect,nchs &
      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      if(wid) then
        do i=1,nchs
          enfpcm_cum(i)=collect(i)
        enddo
       else
        do i=1,nchs
          enfpcm_cum(i)=0
        enddo
      endif

      call mpi_reduce(enfpcm_cm2,collect,nchs &
      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      if(wid) then
        do i=1,nchs
          enfpcm_cm2(i)=collect(i)
        enddo
       else
        do i=1,nchs
          enfpcm_cm2(i)=0
        enddo
      endif


      if(wid) then
        spcmcum=spcmcollect
        spcmcm2=spcm2collect
        vpcmcum=vpcmcollect
        vpcmcm2=vpcm2collect
        qopcm_cum=qopcmcollect
        qopcm_cm2=qopcm2collect
       else
        spcmcum=0
        spcmcm2=0
        vpcmcum=0
        vpcmcm2=0
        qopcm_cum=0
        qopcm_cm2=0
      endif

      return
      end

      subroutine pcm_reduce_chvol

      use contrl_file, only: ounit
      use mpi
      use mpiconf, only: nproc
      use pcm,     only: MCHV
      use pcm_cntrl, only: ipcm
      use pcm_fdc, only: qvol
      use pcm_parms, only: ch,nch,nchs,nchv,xpol
      use pcm_xv_new, only: xv_new
      use precision_kinds, only: dp
      implicit none

      integer :: i, ierr, nchv3
      integer, dimension(0:nproc) :: nchv_proc
      integer, dimension(0:nproc) :: icount
      integer, dimension(0:nproc) :: idispl
      real(dp), dimension(nproc) :: charge

      if(ipcm.eq.0) return

      call mpi_allgather(nchv,1,mpi_integer,nchv_proc,1,mpi_integer, &
        MPI_COMM_WORLD,ierr)

      nch=nchs
      do i=0,nproc-1
        nch=nch+nchv_proc(i)
      enddo
      if(nch.gt.MCHV) call fatal_error('PCM_REDUCE: increase MCHV')

!     call mpi_allgather(qvol,1,mpi_double_precision,charge,1,mpi_double_precision,
!    &  MPI_COMM_WORLD,ierr)
!     qvol=0
!     do 15 i=0,nproc-1
! 15    qvol=qvol+charge(i)

!     tmp=(1.0d0-fs)/nscv
!     qvol=tmp/nproc
      do i=nchs+1,nch
        ch(i)=qvol
      enddo

      write(ounit,*) ' *** pcm update *** ',nch-nchs,qvol

      idispl(0)=3*nchs
      icount(0)=3*nchv_proc(0)
      do i=1,nproc-1
        icount(i)=3*nchv_proc(i)
        idispl(i)=idispl(i-1)+3*nchv_proc(i-1)
      enddo
!     write(ounit,*) ' *** pcm update *** ',(icount(i),i=0,nproc-1)
!     write(ounit,*) ' *** pcm update *** ',(idispl(i),i=0,nproc-1)

      nchv3=3*nchv
      call mpi_allgatherv(xv_new,nchv3,mpi_double_precision, &
        xpol,icount,idispl,mpi_double_precision,MPI_COMM_WORLD,ierr)

      return
      end
end module
