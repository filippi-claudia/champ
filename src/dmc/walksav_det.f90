module walksav_det_mod
      use branch,  only: nwalk
      use csfs,    only: nstates
      use dmc_mod, only: mwalk
      use mpi
      use mstates_mod, only: MSTATES
      use multidet, only: ivirt, ndetiab, numrep_det, ndetsingle, ndetdouble
      use multimat, only: aa,wfmat
      use multislater, only: detiab
      use orbval,  only: dorb,orb
      use precision_kinds, only: dp
      use slater,  only: ddx,kref,ndet,norb,slmi
      use system,  only: ndn,nelec,nup
      use vmc_mod, only: MEXCIT,nmat_dim,norb_tot
      use ycompact, only: ymat
      implicit none

      integer, allocatable, save :: krefw(:)
      real(dp), allocatable, save :: slmuiw(:, :)
      real(dp), allocatable, save :: slmdiw(:, :)
      real(dp), allocatable, save :: ddxw(:, :, :)
      real(dp), allocatable, save :: detuw(:, :)
      real(dp), allocatable, save :: detdw(:, :)
      real(dp), allocatable, save :: aaw(:,:,:,:)
      real(dp), allocatable, save :: wfmatw(:,:,:,:)
      real(dp), allocatable, save :: ymatw(:,:,:,:,:)
      real(dp), allocatable, save :: orbw(:,:,:)
      real(dp), allocatable, save :: dorbw(:,:,:,:)

contains
      subroutine walksav_det(iw)
! Written by Claudia Filippi
      implicit none

      integer :: iab
      integer :: iw
      integer :: k, kcum
      integer :: ndim, nel, ndim2

      if(.not.allocated(aaw)) allocate(aaw(nelec,norb_tot,mwalk,2))
      if(.not.allocated(wfmatw)) allocate(wfmatw(ndet,MEXCIT**2,mwalk,2))
      if(.not.allocated(ymatw)) allocate(ymatw(norb_tot,nelec,mwalk,2,MSTATES))
      if(.not.allocated(orbw)) allocate(orbw(nelec,norb_tot,mwalk))
      if(.not.allocated(dorbw)) allocate(dorbw(norb,nelec,3,mwalk))

      if(.not.allocated(krefw)) allocate(krefw(mwalk), source=0)
      if(.not.allocated(slmuiw)) allocate(slmuiw(nmat_dim,mwalk))
      if(.not.allocated(slmdiw)) allocate(slmdiw(nmat_dim,mwalk))
      if(.not.allocated(ddxw)) allocate(ddxw(3, nelec,mwalk))
      if(.not.allocated(detuw)) allocate(detuw(ndet,mwalk))
      if(.not.allocated(detdw)) allocate(detdw(ndet,mwalk))

       detuw(1:ndet,iw)=detiab(1:ndet,1,1)
       detdw(1:ndet,iw)=detiab(1:ndet,2,1)

       krefw(iw)=kref
       slmuiw(1:nup*nup,iw)=slmi(1:nup*nup,1,1)
       slmdiw(1:ndn*ndn,iw)=slmi(1:ndn*ndn,2,1)
       ddxw(1:3,1:nelec,iw)=ddx(1:3,1:nelec,1)

       do iab=1,2
         nel=nup
         if(iab.eq.2) nel=ndn
         ymatw(ivirt(iab):norb,1:nel,iw,iab,1:nstates)= &
             ymat(ivirt(iab):norb,1:nel,iab,1:nstates)
         aaw(1:nel,ivirt(iab):norb,iw,iab)=aa(1:nel,ivirt(iab):norb,iab,1)

! loop over unique or unequivalent determinants
! single excitations
        if(ndetsingle(iab).ge.1) &
          wfmatw(1:ndetsingle(iab),1,iw,iab)=wfmat(1:ndetsingle(iab),1,iab,1)

! double excitations  
        kcum=ndetsingle(iab)+ndetdouble(iab)
        if(ndetdouble(iab).ge.1) &
          wfmatw(ndetsingle(iab)+1:kcum,1:4,iw,iab)= &
              wfmat(ndetsingle(iab)+1:kcum,1:4,iab,1)

! multiple excitations
        do k=kcum+1,ndetiab(iab)
          ndim=numrep_det(k,iab)
          ndim2=ndim*ndim
          wfmatw(k,1:ndim2,iw,iab)=wfmat(k,1:ndim2,iab,1)
        enddo

       enddo

       orbw(1:nelec,1:norb,iw)=orb(1:nelec,1:norb,1)
       dorbw(1:norb,1:nelec,1:3,iw)=dorb(1:norb,1:nelec,1:3,1)

      end subroutine walksav_det

      subroutine walkstrdet(iw)
      implicit none

      integer :: iab
      integer :: iw
      integer :: k, kcum
      integer :: ndim, nel, ndim2

      detiab(1:ndet,1,1)=detuw(1:ndet,iw)
      detiab(1:ndet,2,1)=detdw(1:ndet,iw)

      kref=krefw(iw)
      slmi(1:nup*nup,1,1)=slmuiw(1:nup*nup,iw)
      slmi(1:ndn*ndn,2,1)=slmdiw(1:ndn*ndn,iw)
      ddx(1:3,1:nelec,1)=ddxw(1:3,1:nelec,iw)

       do iab=1,2
         nel=nup
         if(iab.eq.2) nel=ndn
         ymat(ivirt(iab):norb,1:nel,iab,1:nstates)= &
             ymatw(ivirt(iab):norb,1:nel,iw,iab,1:nstates)
         aa(1:nel,ivirt(iab):norb,iab,1)=aaw(1:nel,ivirt(iab):norb,iw,iab)

! loop over unique or unequivalent determinants
! single excitations
        if(ndetsingle(iab).ge.1) &
          wfmat(1:ndetsingle(iab),1,iab,1)=wfmatw(1:ndetsingle(iab),1,iw,iab)

! double excitations
        kcum=ndetsingle(iab)+ndetdouble(iab)
        if(ndetdouble(iab).ge.1) &
          wfmat(ndetsingle(iab)+1:kcum,1:4,iab,1)= &
              wfmatw(ndetsingle(iab)+1:kcum,1:4,iw,iab)

! multiple excitations
        do k=kcum+1,ndetiab(iab)
          ndim=numrep_det(k,iab)
          ndim2=ndim*ndim
          wfmat(k,1:ndim2,iab,1)=wfmatw(k,1:ndim2,iw,iab)
        enddo

       enddo

       orb(1:nelec,1:norb,1)=orbw(1:nelec,1:norb,iw)
       dorb(1:norb,1:nelec,1:3,1)=dorbw(1:norb,1:nelec,1:3,iw)

      end subroutine walkstrdet

      subroutine splitjdet(iw,iw2)
      implicit none

      integer :: iab
      integer :: iw, iw2
      integer :: k, kcum
      integer :: ndim, nel, ndim2

      detuw(1:ndet,iw2)=detuw(1:ndet,iw)
      detdw(1:ndet,iw2)=detdw(1:ndet,iw)

      krefw(iw2)=krefw(iw)
      slmuiw(1:nup*nup,iw2)=slmuiw(1:nup*nup,iw)
      slmdiw(1:ndn*ndn,iw2)=slmdiw(1:ndn*ndn,iw)
      ddxw(1:3,1:nelec,iw2)=ddxw(1:3,1:nelec,iw)

       do iab=1,2
         nel=nup
         if(iab.eq.2) nel=ndn
         ymatw(ivirt(iab):norb,1:nel,iw2,iab,1:nstates)= &
             ymatw(ivirt(iab):norb,1:nel,iw,iab,1:nstates)
         aaw(1:nel,ivirt(iab):norb,iw2,iab)=aaw(1:nel,ivirt(iab):norb,iw,iab)

! loop over unique or unequivalent determinants
! single excitations
        if(ndetsingle(iab).ge.1) &
          wfmatw(1:ndetsingle(iab),1,iw2,iab)=wfmatw(1:ndetsingle(iab),1,iw,iab)

! double excitations
        kcum=ndetsingle(iab)+ndetdouble(iab)
        if(ndetdouble(iab).ge.1) &
          wfmatw(ndetsingle(iab)+1:kcum,1:4,iw2,iab)= &
              wfmatw(ndetsingle(iab)+1:kcum,1:4,iw,iab)

! multiple excitations
        do k=kcum+1,ndetiab(iab)
          ndim=numrep_det(k,iab)
          ndim2=ndim*ndim
          wfmatw(k,1:ndim2,iw2,iab)=wfmatw(k,1:ndim2,iw,iab)
        enddo

       enddo

       orbw(1:nelec,1:norb,iw2)=orbw(1:nelec,1:norb,iw)
       dorbw(1:norb,1:nelec,1:3,iw2)=dorbw(1:norb,1:nelec,1:3,iw)

      end subroutine splitjdet

      subroutine send_det(irecv)
      implicit none

      integer :: iab, ierr, irecv
      integer :: irequest, istate
      integer :: itag
      integer :: k
      integer :: ndim

      itag=0
      call mpi_isend(detuw(1,nwalk),ndet,mpi_double_precision,irecv &
      ,itag+1,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(detdw(1,nwalk),ndet,mpi_double_precision,irecv &
      ,itag+2,MPI_COMM_WORLD,irequest,ierr)
      itag=itag+2

      call mpi_isend(slmuiw(1,nwalk),nup*nup,mpi_double_precision &
        ,irecv,itag+1,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(slmdiw(1,nwalk),ndn*ndn,mpi_double_precision &
        ,irecv,itag+2,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(ddxw(1,1,nwalk),3*nelec,mpi_double_precision &
        ,irecv,itag+3,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(krefw(nwalk),1,mpi_integer &
        ,irecv,itag+4,MPI_COMM_WORLD,irequest,ierr)
      itag=itag+4

      call mpi_isend(aaw(1,1,nwalk,1),nelec*norb,mpi_double_precision &
       ,irecv,itag+1,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(aaw(1,1,nwalk,2),nelec*norb,mpi_double_precision &
       ,irecv,itag+2,MPI_COMM_WORLD,irequest,ierr)
      itag=itag+2

      do istate=1,nstates
        do iab=1,2
        itag=itag+1
        call mpi_isend(ymatw(1,1,nwalk,iab,istate),norb_tot*nelec,mpi_double_precision &
         ,irecv,itag,MPI_COMM_WORLD,irequest,ierr)
        enddo
      enddo

      do iab=1,2
        do k=1,ndetiab(iab)
          ndim=numrep_det(k,iab)
          itag=itag+1
          call mpi_isend(wfmatw(k,:,nwalk,iab),ndim*ndim,mpi_double_precision &
           ,irecv,itag,MPI_COMM_WORLD,irequest,ierr)
        enddo
      enddo

      call mpi_isend(orbw(1,1,nwalk),nelec*norb,mpi_double_precision &
        ,irecv,itag+1,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(dorbw(1,1,1,nwalk),3*nelec*norb,mpi_double_precision &
        ,irecv,itag+2,MPI_COMM_WORLD,irequest,ierr)

      end subroutine send_det

      subroutine recv_det(isend)
      implicit none

      integer :: iab, ierr, isend, istate
      integer :: istatus, itag
      integer :: k
      integer :: ndim

      dimension istatus(MPI_STATUS_SIZE)

      itag=0
      call mpi_recv(detuw(1,nwalk),ndet,mpi_double_precision,isend &
      ,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(detdw(1,nwalk),ndet,mpi_double_precision,isend &
      ,itag+2,MPI_COMM_WORLD,istatus,ierr)
      itag=itag+2

      call mpi_recv(slmuiw(1,nwalk),nup*nup,mpi_double_precision &
        ,isend,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(slmdiw(1,nwalk),ndn*ndn,mpi_double_precision &
        ,isend,itag+2,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(ddxw(1,1,nwalk),3*nelec,mpi_double_precision &
        ,isend,itag+3,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(krefw(nwalk),1,mpi_integer &
        ,isend,itag+4,MPI_COMM_WORLD,istatus,ierr)
      itag=itag+4

      call mpi_recv(aaw(1,1,nwalk,1),nelec*norb,mpi_double_precision &
       ,isend,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(aaw(1,1,nwalk,2),nelec*norb,mpi_double_precision &
       ,isend,itag+2,MPI_COMM_WORLD,istatus,ierr)
      itag=itag+2

      do istate=1,nstates
        do iab=1,2
        itag=itag+1
        call mpi_recv(ymatw(1,1,nwalk,iab,istate),norb_tot*nelec,mpi_double_precision &
         ,isend,itag,MPI_COMM_WORLD,istatus,ierr)
        enddo
      enddo

      do iab=1,2
        do k=1,ndetiab(iab)
          ndim=numrep_det(k,iab)
          itag=itag+1
          call mpi_recv(wfmatw(k,:,nwalk,iab),ndim*ndim,mpi_double_precision &
           ,isend,itag,MPI_COMM_WORLD,istatus,ierr)
        enddo
      enddo

      call mpi_recv(orbw(1,1,nwalk),nelec*norb,mpi_double_precision &
        ,isend,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(dorbw(1,1,1,nwalk),3*nelec*norb,mpi_double_precision &
        ,isend,itag+2,MPI_COMM_WORLD,istatus,ierr)

      end subroutine recv_det
end module walksav_det_mod
