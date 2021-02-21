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

      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X,
     &NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20,
     &radmax, delri, NEQSX, MTERMS, MCENT3, NCOEF, MEXCIT
      use dmc_mod, only: MWALK, MFPROD, MFPRD1, MPATH
      use const, only: delta, deltai, etrial, fbias, hb, imetro, ipr, nelec, pi
      use forcest, only: fgcm2, fgcum
      use forcepar, only: deltot, istrech, nforce
      use force_dmc, only: itausec, nwprod
      use mpiconf, only: idtask, nproc, wid, NPROCX

      include 'force.h'
      include 'mpif.h'


      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk

      dimension nwalk_all(0:NPROCX),icommunicate_all(0:NPROCX),
     &iwalk_stack(NPROCX)

c     call mpi_allgather(nwalk,1,mpi_integer,nwalk_all,1,mpi_integer,
c    &MPI_COMM_WORLD,ierr)

      call mpi_gather(nwalk,1,mpi_integer,nwalk_all,1,mpi_integer,0,
     &MPI_COMM_WORLD,ierr)
      call mpi_bcast(nwalk_all(0),nproc,mpi_integer,0,MPI_COMM_WORLD,ierr)

      if(ipr.ge.1) write(6,'(''nwalk_all='',(10i4))') (nwalk_all(i),i=0,nproc-1)

      nwalk_sum=0
      do 10 i=0,nproc-1
   10   nwalk_sum=nwalk_sum+nwalk_all(i)

      nwalk_av_int=nwalk_sum/nproc
      nhi=nwalk_sum-nwalk_av_int*nproc
      nlo=(nwalk_sum-nhi*(nwalk_av_int+1))/nwalk_av_int
      if(ipr.ge.1) write(6,'(''nwalk_sum,nwalk_av_int,nlo,nhi='',9i4)')
     &nwalk_sum,nwalk_av_int,nlo,nhi

   20 ilo=0
      ihi=0
      nwalk_stack=0
      do 30 i=0,nproc-1
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

   30 continue
      if(ipr.ge.1) then
        write(6,'(''icommunicate_all='',(10i4))')
     &   (icommunicate_all(i),i=0,nproc-1)
        write(6,'(''iwalk_stack='',(10i4))')
     &   (iwalk_stack(i),i=1,nproc)
        write(6,'(''nwalk_all='',(10i4))') (nwalk_all(i),i=0,nproc-1)
        write(6,'(''nwalk_stack='',i4)') nwalk_stack
        write(6,*)
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
      do 40 i=0,nproc-1
   40 if(nwalk_all(i).gt.nwalk_av_int+1 .or.
     &  nwalk_all(i).lt.nwalk_av_int) ido_again=1

      if(ido_again.eq.1) goto 20

      return
      end
