      module redistribute_mod
      contains
      subroutine redistribute
c Figure out who communicates a walker to whom to achieve load balance.
c It communicates only 1 walker at a time between processors, but it does
c this as many times as needed to achieve load balance.  Most of the time
c only one go around will be needed.
c nwalk_stack > 0 there are  nwalk_stack processors with extra walkers on stack
c             < 0 there are -nwalk_stack processors with fewer walkers on stack
c iwalk_stack(i)  which processors are on the stack
c icommunicate(i) = 0 processor i does not communicate with others
c                 > 0 processor i sends walker to processor icommunicate(i)-1
c                 < 0 processor i receives walker from processor abs(icommunicate(i))-1
c nlo                 # of processors that should have nwalk_av_int walkers
c nhi                 # of processors that should have nwalk_av_int+1 walkers
c Written by Cyrus Umrigar and Claudia Filippi, Oct. 2001.

      use branch,  only: nwalk
      use contrl_file, only: ounit
      use control, only: ipr
      use move_walker, only: recv_walker,send_walker
      use mpi
      use mpiconf, only: idtask,nproc
      use walksav_det_mod, only: recv_det,send_det
      use walksav_jas_mod, only: recv_jas,send_jas
      implicit none

      integer :: i, icomm, ido_again, ierr, ihi
      integer :: ilo, nhi, nlo, nwalk_av_int
      integer :: nwalk_stack, nwalk_sum
      integer, dimension(0:nproc) :: nwalk_all
      integer, dimension(0:nproc) :: icommunicate_all
      integer, dimension(nproc) :: iwalk_stack



c     call mpi_allgather(nwalk,1,mpi_integer,nwalk_all,1,mpi_integer,
c    &MPI_COMM_WORLD,ierr)

      call mpi_gather(nwalk,1,mpi_integer,nwalk_all,1,mpi_integer,0,
     &MPI_COMM_WORLD,ierr)
      call mpi_bcast(nwalk_all(0),nproc,mpi_integer,0,MPI_COMM_WORLD,ierr)

      if(ipr.ge.1) write(ounit,'(''nwalk_all='',(10i4))') (nwalk_all(i),i=0,nproc-1)

      nwalk_sum=0
      do i=0,nproc-1
        nwalk_sum=nwalk_sum+nwalk_all(i)
      enddo

      nwalk_av_int=nwalk_sum/nproc
      nhi=nwalk_sum-nwalk_av_int*nproc
      nlo=(nwalk_sum-nhi*(nwalk_av_int+1))/nwalk_av_int
      if(ipr.ge.1) write(ounit,'(''nwalk_sum,nwalk_av_int,nlo,nhi='',9i4)')
     &nwalk_sum,nwalk_av_int,nlo,nhi

   20 ilo=0
      ihi=0
      nwalk_stack=0
      do i=0,nproc-1
        if(nwalk_all(i).ge.nwalk_av_int+1) then
          ihi=ihi+1
         elseif(nwalk_all(i).le.nwalk_av_int) then
          ilo=ilo+1
        endif

        icommunicate_all(i)=0
        if(nwalk_all(i).gt.nwalk_av_int+1 .or.
     &    (nwalk_all(i).eq.nwalk_av_int+1 .and. ihi.gt.nhi)) then
          if(nwalk_stack.lt.0) then
            icommunicate_all(i)=iwalk_stack(-nwalk_stack)+1
            icommunicate_all(iwalk_stack(-nwalk_stack))=-(i+1)
            nwalk_all(i)=nwalk_all(i)-1
            nwalk_all(iwalk_stack(-nwalk_stack))=
     &      nwalk_all(iwalk_stack(-nwalk_stack))+1
           else
            iwalk_stack(nwalk_stack+1)=i
          endif
          nwalk_stack=nwalk_stack+1
         elseif(nwalk_all(i).lt.nwalk_av_int .or.
     &    (nwalk_all(i).eq.nwalk_av_int .and. ilo.gt.nlo)) then
          if(nwalk_stack.gt.0) then
            icommunicate_all(i)=-(iwalk_stack(nwalk_stack)+1)
            icommunicate_all(iwalk_stack(nwalk_stack))=i+1
            nwalk_all(i)=nwalk_all(i)+1
            nwalk_all(iwalk_stack(nwalk_stack))=
     &      nwalk_all(iwalk_stack(nwalk_stack))-1
           else
            iwalk_stack(-nwalk_stack+1)=i
          endif
          nwalk_stack=nwalk_stack-1
        endif

      enddo
      if(ipr.ge.1) then
        write(ounit,'(''icommunicate_all='',(10i4))')
     &   (icommunicate_all(i),i=0,nproc-1)
        write(ounit,'(''iwalk_stack='',(10i4))')
     &   (iwalk_stack(i),i=1,nproc)
        write(ounit,'(''nwalk_all='',(10i4))') (nwalk_all(i),i=0,nproc-1)
        write(ounit,'(''nwalk_stack='',i4)') nwalk_stack
        write(ounit,*)
      endif

c     call flush(6)
c call routine to move walkers as specified in icommunicate_all
c and to update the values of nwalk on each processor
      icomm=icommunicate_all(idtask)
      if(icomm.lt.0) nwalk=nwalk+1

      if(icomm.gt.0) call send_walker(icomm-1)
      if(icomm.lt.0) call recv_walker(-icomm-1)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(icomm.gt.0) call send_det(icomm-1)
      if(icomm.lt.0) call recv_det(-icomm-1)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(icomm.gt.0) call send_jas(icomm-1)
      if(icomm.lt.0) call recv_jas(-icomm-1)

      if(icomm.gt.0) nwalk=nwalk-1

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      ido_again=0
      do i=0,nproc-1
      if(nwalk_all(i).gt.nwalk_av_int+1 .or.
     &  nwalk_all(i).lt.nwalk_av_int) ido_again=1
      enddo

      if(ido_again.eq.1) goto 20

      return
      end
      end module
