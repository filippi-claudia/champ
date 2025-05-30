module walksav_det_mod
      use branch,  only: nwalk
      use csfs,    only: nstates
      use dmc_mod, only: mwalk
      use mpi
      use mstates_mod, only: MSTATES
      use multidet, only: ivirt,numrep_det
      use multimat, only: aa,wfmat
      use multislater, only: detiab
      use orbval,  only: dorb,orb
      use precision_kinds, only: dp
      use slater,  only: ddx,fp,kref,ndet,norb,slmi
      use system,  only: ndn,nelec,nup
      use vmc_mod, only: MEXCIT,nmat_dim,norb_tot
      use ycompact, only: ymat
      implicit none

      integer, allocatable, save :: krefw(:)
      real(dp), allocatable, save :: slmuiw(:, :)
      real(dp), allocatable, save :: slmdiw(:, :)
      real(dp), allocatable, save :: fpuw(:, :, :)
      real(dp), allocatable, save :: fpdw(:, :, :)
      real(dp), allocatable, save :: fppuw(:, :)
      real(dp), allocatable, save :: fppdw(:, :)
      real(dp), allocatable, save :: ddxw(:, :, :)
      real(dp), allocatable, save :: d2dx2w(:, :)
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

      integer :: i, iab, ierr, iorb, irecv
      integer :: irequest, irequest_array, isend, istate
      integer :: istatus, itag, iw, iw2
      integer :: j, k, kk
      integer :: ndim, nel, ndim2

      dimension istatus(MPI_STATUS_SIZE)
      dimension irequest_array(MPI_STATUS_SIZE)

      if(.not.allocated(aaw)) allocate(aaw(nelec,norb_tot,mwalk,2))
      if(.not.allocated(wfmatw)) allocate(wfmatw(ndet,MEXCIT**2,mwalk,2))
      if(.not.allocated(ymatw)) allocate(ymatw(norb_tot,nelec,mwalk,2,MSTATES))
      if(.not.allocated(orbw)) allocate(orbw(nelec,norb_tot,mwalk))
      if(.not.allocated(dorbw)) allocate(dorbw(3,nelec,norb_tot,mwalk))

      if(.not.allocated(krefw)) allocate(krefw(mwalk), source=0)
      if(.not.allocated(slmuiw)) allocate(slmuiw(nmat_dim,mwalk))
      if(.not.allocated(slmdiw)) allocate(slmdiw(nmat_dim,mwalk))
      if(.not.allocated(fpuw)) allocate(fpuw(3, nmat_dim,mwalk))
      if(.not.allocated(fpdw)) allocate(fpdw(3, nmat_dim,mwalk))
      if(.not.allocated(fppuw)) allocate(fppuw(nmat_dim,mwalk))
      if(.not.allocated(fppdw)) allocate(fppdw(nmat_dim,mwalk))
      if(.not.allocated(ddxw)) allocate(ddxw(3, nelec,mwalk))
      if(.not.allocated(d2dx2w)) allocate(d2dx2w(nelec,mwalk))
      if(.not.allocated(detuw)) allocate(detuw(ndet,mwalk))
      if(.not.allocated(detdw)) allocate(detdw(ndet,mwalk))



       do k=1,ndet
         detuw(k,iw)=detiab(k,1,1)
         detdw(k,iw)=detiab(k,2,1)
       enddo

       krefw(iw)=kref
       do j=1,nup*nup
         slmuiw(j,iw)=slmi(j,1,1)
         fpuw(1,j,iw)=fp(1,j,1,1)
         fpuw(2,j,iw)=fp(2,j,1,1)
         fpuw(3,j,iw)=fp(3,j,1,1)
       enddo
       do j=1,ndn*ndn
         slmdiw(j,iw)=slmi(j,2,1)
         fpdw(1,j,iw)=fp(1,j,2,1)
         fpdw(2,j,iw)=fp(2,j,2,1)
         fpdw(3,j,iw)=fp(3,j,2,1)
       enddo
       do i=1,nelec
         ddxw(1,i,iw)=ddx(1,i,1)
         ddxw(2,i,iw)=ddx(2,i,1)
         ddxw(3,i,iw)=ddx(3,i,1)
       enddo

       do iab=1,2
         nel=nup
         if(iab.eq.2) nel=ndn
         do j=ivirt(iab),norb
          do i=1,nel
            do istate=1,nstates
              ymatw(j,i,iw,iab,istate)=ymat(j,i,iab,istate)
            enddo
            aaw(i,j,iw,iab)=aa(i,j,iab,1)
          enddo
         enddo
         do k=1,kref-1
              ndim=numrep_det(k,iab)
              ndim2=ndim*ndim
              wfmatw(k,1:ndim2,iw,iab)=wfmat(k,1:ndim2,iab,1)
         enddo
         do k=kref+1,ndet
            ndim=numrep_det(k,iab)
            ndim2=ndim*ndim
            wfmatw(k,1:ndim2,iw,iab)=wfmat(k,1:ndim2,iab,1)
          enddo
       enddo

       do i=1,nelec
         do iorb=1,norb
           orbw(i,iorb,iw)=orb(i,iorb,1)
           do kk=1,3
             dorbw(kk,i,iorb,iw)=dorb(iorb,i,kk,1)
           enddo
         enddo
       enddo

      end subroutine

      subroutine walkstrdet(iw)
      implicit none

      integer :: i, iab, ierr, iorb, irecv
      integer :: irequest, irequest_array, isend, istate
      integer :: istatus, itag, iw, iw2
      integer :: j, k, kk
      integer :: ndim, nel, ndim2

      dimension istatus(MPI_STATUS_SIZE)
      dimension irequest_array(MPI_STATUS_SIZE)

      do k=1,ndet
        detiab(k,1,1)=detuw(k,iw)
        detiab(k,2,1)=detdw(k,iw)
      enddo

      kref=krefw(iw)
      do j=1,nup*nup
        slmi(j,1,1)=slmuiw(j,iw)
        fp(1,j,1,1)=fpuw(1,j,iw)
        fp(2,j,1,1)=fpuw(2,j,iw)
        fp(3,j,1,1)=fpuw(3,j,iw)
      enddo
      do j=1,ndn*ndn
        slmi(j,2,1)=slmdiw(j,iw)
        fp(1,j,2,1)=fpdw(1,j,iw)
        fp(2,j,2,1)=fpdw(2,j,iw)
        fp(3,j,2,1)=fpdw(3,j,iw)
      enddo
      do i=1,nelec
        ddx(1,i,1)=ddxw(1,i,iw)
        ddx(2,i,1)=ddxw(2,i,iw)
        ddx(3,i,1)=ddxw(3,i,iw)
      enddo

       do iab=1,2
         nel=nup
         if(iab.eq.2) nel=ndn
         do j=ivirt(iab),norb
          do i=1,nel
            do istate=1,nstates
              ymat(j,i,iab,istate)=ymatw(j,i,iw,iab,istate)
            enddo
            aa(i,j,iab,1)=aaw(i,j,iw,iab)
          enddo
         enddo
          do k=1,kref-1
             ndim=numrep_det(k,iab)
             ndim2=ndim*ndim
             wfmat(k,1:ndim2,iab,1)=wfmatw(k,1:ndim2,iw,iab)
          enddo
          do k=kref+1,ndet
              ndim=numrep_det(k,iab)
              ndim2=ndim*ndim
              wfmat(k,1:ndim2,iab,1)=wfmatw(k,1:ndim2,iw,iab)
          enddo
       enddo

       do i=1,nelec
         do iorb=1,norb
           orb(i,iorb,1)=orbw(i,iorb,iw)
           do kk=1,3
             dorb(iorb,i,kk,1)=dorbw(kk,i,iorb,iw)
           enddo
         enddo
       enddo

      end subroutine

      subroutine splitjdet(iw,iw2)
      implicit none

      integer :: i, iab, ierr, iorb, irecv
      integer :: irequest, irequest_array, isend, istate
      integer :: istatus, itag, iw, iw2
      integer :: j, k, kk
      integer :: ndim, nel, ndim2

      dimension istatus(MPI_STATUS_SIZE)
      dimension irequest_array(MPI_STATUS_SIZE)

      do k=1,ndet
        detuw(k,iw2)=detuw(k,iw)
        detdw(k,iw2)=detdw(k,iw)
      enddo

      krefw(iw2)=krefw(iw)
      do j=1,nup*nup
        slmuiw(j,iw2)=slmuiw(j,iw)
        fpuw(1,j,iw2)=fpuw(1,j,iw)
        fpuw(2,j,iw2)=fpuw(2,j,iw)
        fpuw(3,j,iw2)=fpuw(3,j,iw)
      enddo
      do j=1,ndn*ndn
        slmdiw(j,iw2)=slmdiw(j,iw)
        fpdw(1,j,iw2)=fpdw(1,j,iw)
        fpdw(2,j,iw2)=fpdw(2,j,iw)
        fpdw(3,j,iw2)=fpdw(3,j,iw)
      enddo
      do i=1,nelec
        ddxw(1,i,iw2)=ddxw(1,i,iw)
        ddxw(2,i,iw2)=ddxw(2,i,iw)
        ddxw(3,i,iw2)=ddxw(3,i,iw)
      enddo

       do iab=1,2
         nel=nup
         if(iab.eq.2) nel=ndn
         do j=ivirt(iab),norb
          do i=1,nel
            do istate=1,nstates
              ymatw(j,i,iw2,iab,istate)=ymatw(j,i,iw,iab,istate)
            enddo
            aaw(i,j,iw2,iab)=aaw(i,j,iw,iab)
          enddo
         enddo
          do k=1,krefw(iw)-1
              ndim=numrep_det(k,iab)
              ndim2=ndim*ndim
              wfmatw(k,1:ndim2,iw2,iab)=wfmatw(k,1:ndim2,iw,iab)
          enddo
          do k=krefw(iw)+1,ndet
              ndim=numrep_det(k,iab)
              ndim2=ndim*ndim
              wfmatw(k,1:ndim2,iw2,iab)=wfmatw(k,1:ndim2,iw,iab)
          enddo
       enddo

       do i=1,nelec
         do iorb=1,norb
           orbw(i,iorb,iw2)=orbw(i,iorb,iw)
           do kk=1,3
             dorbw(kk,i,iorb,iw2)=dorbw(kk,i,iorb,iw)
           enddo
         enddo
       enddo

      end subroutine

      subroutine send_det(irecv)
      implicit none

      integer :: i, iab, ierr, iorb, irecv
      integer :: irequest, irequest_array, isend, istate
      integer :: istatus, itag, iw, iw2
      integer :: j, k, kk
      integer :: ndim, nel, ndim2

      dimension istatus(MPI_STATUS_SIZE)
      dimension irequest_array(MPI_STATUS_SIZE)

      itag=0
      call mpi_isend(detuw(1,nwalk),ndet,mpi_double_precision,irecv &
      ,itag+1,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(detdw(1,nwalk),ndet,mpi_double_precision,irecv &
      ,itag+2,MPI_COMM_WORLD,irequest,ierr)
      itag=itag+2

      call mpi_isend(slmuiw(1,nwalk),nup*nup,mpi_double_precision &
        ,irecv,itag+1,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(fpuw(1,1,nwalk),3*nup*nup,mpi_double_precision &
        ,irecv,itag+2,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(slmdiw(1,nwalk),ndn*ndn,mpi_double_precision &
        ,irecv,itag+3,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(fpdw(1,1,nwalk),3*ndn*ndn,mpi_double_precision &
        ,irecv,itag+4,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(ddxw(1,1,nwalk),3*nelec,mpi_double_precision &
        ,irecv,itag+5,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(krefw(nwalk),1,mpi_integer &
        ,irecv,itag+6,MPI_COMM_WORLD,irequest,ierr)
      itag=itag+6

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
        do k=1,ndet
          ndim=numrep_det(k,iab)
          if(k.ne.krefw(nwalk).and.ndim.gt.0) then
            itag=itag+1
            call mpi_isend(wfmatw(k,:,nwalk,iab),ndim*ndim,mpi_double_precision &
           ,irecv,itag,MPI_COMM_WORLD,irequest,ierr)
          endif
        enddo
      enddo

      call mpi_isend(orbw(1,1,nwalk),nelec*norb,mpi_double_precision &
        ,irecv,itag+1,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(dorbw(1,1,1,nwalk),3*nelec*norb,mpi_double_precision &
        ,irecv,itag+2,MPI_COMM_WORLD,irequest,ierr)
      itag=itag+2

      end subroutine

      subroutine recv_det(isend)
      implicit none

      integer :: i, iab, ierr, iorb, irecv
      integer :: irequest, irequest_array, isend, istate
      integer :: istatus, itag, iw, iw2
      integer :: j, k, kk
      integer :: ndim, nel, ndim2

      dimension istatus(MPI_STATUS_SIZE)
      dimension irequest_array(MPI_STATUS_SIZE)

      itag=0
      call mpi_recv(detuw(1,nwalk),ndet,mpi_double_precision,isend &
      ,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(detdw(1,nwalk),ndet,mpi_double_precision,isend &
      ,itag+2,MPI_COMM_WORLD,istatus,ierr)
      itag=itag+2

      call mpi_recv(slmuiw(1,nwalk),nup*nup,mpi_double_precision &
        ,isend,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fpuw(1,1,nwalk),3*nup*nup,mpi_double_precision &
        ,isend,itag+2,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(slmdiw(1,nwalk),ndn*ndn,mpi_double_precision &
        ,isend,itag+3,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fpdw(1,1,nwalk),3*ndn*ndn,mpi_double_precision &
        ,isend,itag+4,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(ddxw(1,1,nwalk),3*nelec,mpi_double_precision &
        ,isend,itag+5,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(krefw(nwalk),1,mpi_integer &
        ,isend,itag+6,MPI_COMM_WORLD,irequest_array,ierr)
      itag=itag+6

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
        do k=1,ndet
          ndim=numrep_det(k,iab)
          if(k.ne.krefw(nwalk).and.ndim.gt.0) then
            itag=itag+1
            call mpi_recv(wfmatw(k,:,nwalk,iab),ndim*ndim,mpi_double_precision &
           ,isend,itag,MPI_COMM_WORLD,istatus,ierr)
        endif
        enddo
      enddo

      call mpi_recv(orbw(1,1,nwalk),nelec*norb_tot,mpi_double_precision &
        ,isend,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(dorbw(1,1,1,nwalk),3*nelec*norb_tot,mpi_double_precision &
        ,isend,itag+2,MPI_COMM_WORLD,istatus,ierr)
      itag=itag+2

      return
      end
end module
