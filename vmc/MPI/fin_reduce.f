      subroutine fin_reduce
c MPI version written by Claudia Filippi

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'mpif.h'

      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /estcum/ ecum1(MSTATES),ecum(MSTATES,MFORCE),pecum(MSTATES),tpbcum(MSTATES),tjfcum(MSTATES),r2cum,iblk
      common /est2cm/ ecm21(MSTATES),ecm2(MSTATES,MFORCE),pecm2(MSTATES),tpbcm2(MSTATES),tjfcm2(MSTATES),r2cm2
      common /estsig/ ecum1s(MSTATES),ecm21s(MSTATES)
      common /forcewt/ wsum(MSTATES,MFORCE),wcum(MSTATES,MFORCE)

      common /step/try(nrad),suc(nrad),trunfb(nrad),rprob(nrad),
     &ekin(nrad),ekin2(nrad)

      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      logical wid
      common /mpiconf/ idtask,nproc,wid

      dimension rprobt(nrad),tryt(nrad),suct(nrad)
      dimension collect(MSTATES)
      dimension istatus(MPI_STATUS_SIZE)

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

      if(wid) then
        do 60 id=1,nproc-1
          call mpi_send(ecum1,nstates,mpi_double_precision,id
     &    ,1,MPI_COMM_WORLD,ierr)
c    &    ,1,MPI_COMM_WORLD,irequest,ierr)
   60     call mpi_send(wcum,MSTATES*nforce,mpi_double_precision,id
     &    ,2,MPI_COMM_WORLD,ierr)
c    &    ,2,MPI_COMM_WORLD,irequest,ierr)
       else
        call mpi_recv(ecum1,nstates,mpi_double_precision,0
     &  ,1,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(wcum,MSTATES*nforce,mpi_double_precision,0
     &  ,2,MPI_COMM_WORLD,istatus,ierr)
      endif

      passes=dble(iblk)*dble(nstep)
      efin=ecum1(1)/wcum(1,1)

      call optjas_fin(wcum,ecum1)

      call optci_fin(iblk,wcum(1,1),efin)

      call optorb_fin(wcum(1,1),ecum1)

      if(method.ne.'sr_n'.and.method.ne.'lin_d') then
        call optx_jas_ci_fin(wcum(1,1),efin)

        call optx_jas_orb_fin(wcum(1,1),ecum1)

        call optx_orb_ci_fin(wcum(1,1),efin)
      endif

      call efficiency_prt(passes)

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
