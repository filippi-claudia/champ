      subroutine finwrt_more
c written by Claudia Filippi

      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use est2cm, only: ecm2, ecm21, pecm2, r2cm2, tjfcm2, tpbcm2
      use estcum, only: ecum, ecum1, iblk, pecum, r2cum, tjfcum, tpbcum
      use estpsi, only: apsi, aref, detref
      use estsig, only: ecm21s, ecum1s
      implicit real*8(a-h,o-z)






      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'mpif.h'

      logical wid
      common /mpiconf/ idtask,nproc,wid


      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /estsum/ esum1(MSTATES),esum(MSTATES,MFORCE),pesum(MSTATES),tpbsum(MSTATES),tjfsum(MSTATES),r2sum,acc

      common /optwf_corsam/ add_diag(MFORCE),energy(MFORCE),energy_err(MFORCE),force(MFORCE),force_err(MFORCE),sigma

      dimension istatus(MPI_STATUS_SIZE)

      passes=dfloat(iblk*nstep)
      write(6,'(''average psid, det_ref '',2d12.5)') (apsi(istate)*nproc/passes,istate=1,nstates),aref*nproc/passes
      write(6,'(''log detref '',2d12.5)') (detref(iab)*nproc/passes,iab=1,2)

      if(wid) then
        do 20 id=1,nproc-1
          call mpi_send(energy,3,mpi_double_precision,id
     &    ,1,MPI_COMM_WORLD,ierr)
          call mpi_send(energy_err,3,mpi_double_precision,id
     &    ,2,MPI_COMM_WORLD,ierr)
          call mpi_send(force,3,mpi_double_precision,id
     &    ,3,MPI_COMM_WORLD,ierr)
          call mpi_send(force_err,3,mpi_double_precision,id
     &    ,4,MPI_COMM_WORLD,ierr)
  20      call mpi_send(sigma,1,mpi_double_precision,id
     &    ,5,MPI_COMM_WORLD,ierr)
       else
        call mpi_recv(energy,3,mpi_double_precision,0
     &  ,1,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(energy_err,3,mpi_double_precision,0
     &  ,2,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(force,3,mpi_double_precision,0
     &  ,3,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(force_err,3,mpi_double_precision,0
     &  ,4,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(sigma,1,mpi_double_precision,0
     &  ,5,MPI_COMM_WORLD,istatus,ierr)
      endif

      return
      end
