      module walksav_jas_mod
      contains
      subroutine walksav_jas(iw)
c Written by Claudia Filippi

      use dmc_mod, only: MWALK
      use system, only: nelec
      use branch, only: nwalk
      use jastrow_update, only: fijo, fjo, fso, fsumo
      use velocity_jastrow, only: vj
      use mpi

      use precision_kinds, only: dp
      implicit none

      integer :: i, ierr, irecv, irequest, isend
      integer :: itag, iw, iw2, j
      integer :: kk
      integer, dimension(MPI_STATUS_SIZE) :: istatus


      real(dp), allocatable, save :: fsow(:, :, :)
      real(dp), allocatable, save :: fijow(:, :, :, :)
      real(dp), allocatable, save :: fsumow(:)
      real(dp), allocatable, save :: fjow(:, :, :)
      real(dp), allocatable, save :: vjw(:, :, :)

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

      do i=1,nelec
        fjow(1,i,iw)=fjo(1,i,1)
        fjow(2,i,iw)=fjo(2,i,1)
        fjow(3,i,iw)=fjo(3,i,1)
      enddo

      do i=2,nelec
        do j=1,i-1
        fsow(i,j,iw)=fso(i,j,1)
        fijow(1,i,j,iw)=fijo(1,i,j,1)
        fijow(2,i,j,iw)=fijo(2,i,j,1)
        fijow(3,i,j,iw)=fijo(3,i,j,1)
        fijow(1,j,i,iw)=fijo(1,j,i,1)
        fijow(2,j,i,iw)=fijo(2,j,i,1)
        fijow(3,j,i,iw)=fijo(3,j,i,1)
        enddo
      enddo

      do i=1,nelec
        fsow(i,i,iw)=fso(i,i,1)
        fijow(1,i,i,iw)=fijo(1,i,i,1)
        fijow(2,i,i,iw)=fijo(2,i,i,1)
        fijow(3,i,i,iw)=fijo(3,i,i,1)
      enddo

      do i=1,nelec
        do kk=1,3
          vjw(kk,i,iw)=vj(kk,i,1)
        enddo
      enddo

      return

      entry walkstrjas(iw)

      fsumo(1)=fsumow(iw)

      do i=1,nelec
        fjo(1,i,1)=fjow(1,i,iw)
        fjo(2,i,1)=fjow(2,i,iw)
        fjo(3,i,1)=fjow(3,i,iw)
      enddo

      do i=2,nelec
        do j=1,i-1
        fso(i,j,1)=fsow(i,j,iw)
        fijo(1,i,j,1)=fijow(1,i,j,iw)
        fijo(2,i,j,1)=fijow(2,i,j,iw)
        fijo(3,i,j,1)=fijow(3,i,j,iw)
        fijo(1,j,i,1)=fijow(1,j,i,iw)
        fijo(2,j,i,1)=fijow(2,j,i,iw)
        fijo(3,j,i,1)=fijow(3,j,i,iw)
        enddo
      enddo

      do i=1,nelec
        fso(i,i,1)=fsow(i,i,iw)
        fijo(1,i,i,1)=fijow(1,i,i,iw)
        fijo(2,i,i,1)=fijow(2,i,i,iw)
        fijo(3,i,i,1)=fijow(3,i,i,iw)
      enddo

      do i=1,nelec
        do kk=1,3
          vj(kk,i,1)=vjw(kk,i,iw)
        enddo
      enddo

      return

      entry splitjjas(iw,iw2)

      fsumow(iw2)=fsumow(iw)

      do i=1,nelec
        fjow(1,i,iw2)=fjow(1,i,iw)
        fjow(2,i,iw2)=fjow(2,i,iw)
        fjow(3,i,iw2)=fjow(3,i,iw)
      enddo

      do i=2,nelec
        do j=1,i-1
        fsow(i,j,iw2)=fsow(i,j,iw)
        fijow(1,i,j,iw2)=fijow(1,i,j,iw)
        fijow(2,i,j,iw2)=fijow(2,i,j,iw)
        fijow(3,i,j,iw2)=fijow(3,i,j,iw)
        fijow(1,j,i,iw2)=fijow(1,j,i,iw)
        fijow(2,j,i,iw2)=fijow(2,j,i,iw)
        fijow(3,j,i,iw2)=fijow(3,j,i,iw)
        enddo
      enddo

      do i=1,nelec
        fsow(i,i,iw2)=fsow(i,i,iw)
        fijow(1,i,i,iw2)=fijow(1,i,i,iw)
        fijow(2,i,i,iw2)=fijow(2,i,i,iw)
        fijow(3,i,i,iw2)=fijow(3,i,i,iw)
      enddo

      do i=1,nelec
        do kk=1,3
          vjw(kk,i,iw2)=vjw(kk,i,iw)
        enddo
      enddo

      return

      entry send_jas(irecv)

      itag=0
      call mpi_isend(fsumow(nwalk),1,mpi_double_precision,irecv
     &,itag+1,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(fjow(1,1,nwalk),3*nelec,mpi_double_precision,irecv
     &,itag+2,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(fsow(1,1,nwalk),nelec*nelec,mpi_double_precision
     &,irecv,itag+3,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(fijow(1,1,1,nwalk),3*nelec*nelec
     &,mpi_double_precision,irecv,itag+4,MPI_COMM_WORLD,irequest,ierr)

      call mpi_isend(vjw(1,1,nwalk),3*nelec,mpi_double_precision,irecv
     &,itag+5,MPI_COMM_WORLD,irequest,ierr)

      return

      entry recv_jas(isend)

      itag=0
      call mpi_recv(fsumow(nwalk),1,mpi_double_precision,isend
     &,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fjow(1,1,nwalk),3*nelec,mpi_double_precision,isend
     &,itag+2,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fsow(1,1,nwalk),nelec*nelec,mpi_double_precision
     &,isend,itag+3,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fijow(1,1,1,nwalk),3*nelec*nelec
     &,mpi_double_precision,isend,itag+4,MPI_COMM_WORLD,istatus,ierr)

      call mpi_recv(vjw(1,1,nwalk),3*nelec,mpi_double_precision,isend
     &,itag+5,MPI_COMM_WORLD,istatus,ierr)

      return
      end
      end module
