module walksav_jas_mod

      use branch,  only: nwalk
      use dmc_mod, only: MWALK
      use jastrow_update, only: fijo,fjo,fso,fsumo
      use mpi
      use precision_kinds, only: dp
      use system,  only: nelec
      use velocity_jastrow, only: vj
      implicit none

      real(dp), allocatable, save :: fsow(:, :, :)
      real(dp), allocatable, save :: fijow(:, :, :, :)
      real(dp), allocatable, save :: fsumow(:)
      real(dp), allocatable, save :: fjow(:, :, :)
      real(dp), allocatable, save :: vjw(:, :, :)

contains
      subroutine walksav_jas(iw)
! Written by Claudia Filippi
      implicit none

      integer :: i, iw

      if(.not.allocated(fsow)) allocate(fsow(nelec, nelec, MWALK))
      if(.not.allocated(fijow)) allocate(fijow(3, nelec, nelec, MWALK))
      if(.not.allocated(fsumow)) allocate(fsumow(MWALK))
      if(.not.allocated(fjow)) allocate(fjow(3, nelec, MWALK))
      if(.not.allocated(vjw)) allocate(vjw(3, nelec, MWALK))

! real(dp), dimension(nelec, nelec, MWALK) :: fsow
! real(dp), dimension(3, nelec, nelec, MWALK) :: fijow
! real(dp), dimension(MWALK) :: fsumow
! real(dp), dimension(3, nelec, MWALK) :: fjow
! real(dp), dimension(3, nelec, MWALK) :: vjw

! save fsow,fijow,fsumow,fjow
! save vjw

      fsumow(iw)=fsumo(1)

      fjow(1:3,1:nelec,iw)=fjo(1:3,1:nelec,1)

      do i=2,nelec
        fsow(i,1:i-1,iw)=fso(i,1:i-1,1)
        fijow(1:3,i,1:i-1,iw)=fijo(1:3,i,1:i-1,1)
        fijow(1:3,1:i-1,i,iw)=fijo(1:3,1:i-1,i,1)
      enddo

      do i=1,nelec
        fsow(i,i,iw)=fso(i,i,1)
        fijow(1:3,i,i,iw)=fijo(1:3,i,i,1)
      enddo

      vjw(1:3,1:nelec,iw)=vj(1:3,1:nelec,1)

      end subroutine walksav_jas

      subroutine walkstrjas(iw)
      implicit none

      integer :: i, iw

      fsumo(1)=fsumow(iw)

      fjo(1:3,1:nelec,1)=fjow(1:3,1:nelec,iw)

      do i=2,nelec
        fso(i,1:i-1,1)=fsow(i,1:i-1,iw)
        fijo(1:3,i,1:i-1,1)=fijow(1:3,i,1:i-1,iw)
        fijo(1:3,1:i-1,i,1)=fijow(1:3,1:i-1,i,iw)
      enddo

      do i=1,nelec
        fso(i,i,1)=fsow(i,i,iw)
        fijo(1:3,i,i,1)=fijow(1:3,i,i,iw)
      enddo

      vj(1:3,1:nelec,1)=vjw(1:3,1:nelec,iw)

      end subroutine walkstrjas

      subroutine splitjjas(iw,iw2)
      implicit none

      integer :: i, iw, iw2

      fsumow(iw2)=fsumow(iw)

      fjow(1:3,1:nelec,iw2)=fjow(1:3,1:nelec,iw)

      do i=2,nelec
        fsow(i,1:i-1,iw2)=fsow(i,1:i-1,iw)
        fijow(1:3,i,1:i-1,iw2)=fijow(1:3,i,1:i-1,iw)
        fijow(1:3,1:i-1,i,iw2)=fijow(1:3,1:i-1,i,iw)
      enddo

      do i=1,nelec
        fsow(i,i,iw2)=fsow(i,i,iw)
        fijow(1:3,i,i,iw2)=fijow(1:3,i,i,iw)
      enddo

      vjw(1:3,1:nelec,iw2)=vjw(1:3,1:nelec,iw)

      end subroutine splitjjas

      subroutine send_jas(irecv)
      implicit none

      integer :: ierr, irecv, irequest
      integer :: itag

      itag=0
      call mpi_isend(fsumow(nwalk),1,mpi_double_precision,irecv &
      ,itag+1,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(fjow(1,1,nwalk),3*nelec,mpi_double_precision,irecv &
      ,itag+2,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(fsow(1,1,nwalk),nelec*nelec,mpi_double_precision &
      ,irecv,itag+3,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(fijow(1,1,1,nwalk),3*nelec*nelec &
      ,mpi_double_precision,irecv,itag+4,MPI_COMM_WORLD,irequest,ierr)

      call mpi_isend(vjw(1,1,nwalk),3*nelec,mpi_double_precision,irecv &
      ,itag+5,MPI_COMM_WORLD,irequest,ierr)

      end subroutine send_jas

      subroutine recv_jas(isend)
      implicit none

      integer :: ierr, isend
      integer :: itag
      integer, dimension(MPI_STATUS_SIZE) :: istatus

      itag=0
      call mpi_recv(fsumow(nwalk),1,mpi_double_precision,isend &
      ,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fjow(1,1,nwalk),3*nelec,mpi_double_precision,isend &
      ,itag+2,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fsow(1,1,nwalk),nelec*nelec,mpi_double_precision &
      ,isend,itag+3,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fijow(1,1,1,nwalk),3*nelec*nelec &
      ,mpi_double_precision,isend,itag+4,MPI_COMM_WORLD,istatus,ierr)

      call mpi_recv(vjw(1,1,nwalk),3*nelec,mpi_double_precision,isend &
      ,itag+5,MPI_COMM_WORLD,istatus,ierr)

      end subroutine recv_jas
end module walksav_jas_mod
