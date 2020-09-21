      subroutine finwrt_more
c written by Claudia Filippi

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use csfs, only: nstates

      use estcum, only: iblk
      use estpsi, only: apsi, aref, detref
      use mpiconf, only: nproc, wid
      use optwf_corsam, only: energy, energy_err, force, force_err
      use contrl, only: nstep
      implicit real*8(a-h,o-z)



      include 'mpif.h'




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
