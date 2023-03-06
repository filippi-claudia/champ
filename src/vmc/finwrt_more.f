      module finwrt_more_mod
      contains
      subroutine finwrt_more
c written by Claudia Filippi

      use contrl_file, only: errunit,ounit
      use control_vmc, only: vmc_nstep
      use csfs,    only: nstates
      use estcum,  only: iblk
      use estpsi,  only: apsi,aref,detref
      use mpi
      use mpiconf, only: nproc
      use optwf_corsam, only: energy,energy_err,force,force_err,sigma
      use precision_kinds, only: dp
      use sa_check, only: energy_all,energy_err_all
      implicit none

      integer :: iab, ierr, istate
      real(dp) :: passes


      passes=dfloat(iblk*vmc_nstep)
      write(ounit,'(''average psid, det_ref '',2d12.5)') (apsi(istate)*nproc/passes,istate=1,nstates),aref*nproc/passes
      write(ounit,'(''log detref '',2d12.5)') (detref(iab)*nproc/passes,iab=1,2)

      call mpi_bcast(energy,size(energy),mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(energy_err,size(energy_err),mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(force,size(force),mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(force_err,size(force_err),mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(sigma,1,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(energy_all,nstates,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(energy_err_all,nstates,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      return
      end

      end module
