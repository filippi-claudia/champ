module finwrt_more_mod
contains
      subroutine finwrt_more
! written by Claudia Filippi

      use contrl_file, only: ounit, errunit
      use control_vmc, only: vmc_nstep
      use csfs, only: nstates
      use estcum, only: iblk
      use estpsi, only: apsi, aref, detref
      use mpi
      use mpiconf, only: nproc
      use optwf_corsam, only: energy, energy_err, force, force_err, sigma
      use precision_kinds, only: dp
      use sa_check, only: energy_all, energy_err_all
      use vmc_mod, only: nwftypeorb, stoo
      implicit none

      integer :: iab, ierr, istate, o
      real(dp) :: passes

      passes=dble(iblk*vmc_nstep)

      do istate=1, nstates
        o=stoo(istate)
        write(ounit,'(''average psid, det_ref '',2d12.5)') apsi(istate)*nproc/passes, aref(o)*nproc/passes
        write(ounit,'(''orb set,log detref '',i4,2d13.5)') istate, (detref(iab,o)*nproc/passes,iab=1,2)
      enddo

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
