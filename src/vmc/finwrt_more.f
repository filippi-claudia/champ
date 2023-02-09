      module finwrt_more_mod
      contains
      subroutine finwrt_more
c written by Claudia Filippi

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

      integer :: iab, ierr, istate
      real(dp) :: passes

      !STU use mapping for orbs in aref, detref
      ! in first write statement, apsi differs, but aref only differs if
      ! the orbitals differ
      ! istate and detref in second write statement only differ if orbs
      ! are different
      passes=dfloat(iblk*vmc_nstep)
      
      do istate=1, nstates
        write(ounit,'(''average psid, det_ref '',2d12.5)') apsi(istate)*nproc/passes, aref(stoo(istate))*nproc/passes
        write(ounit,'(''orb set,log detref '',i4,2d13.5)') istate, (detref(iab,stoo(istate))*nproc/passes,iab=1,2)
      enddo

      call mpi_bcast(energy,3,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(energy_err,3,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(force,3,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(force_err,3,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(sigma,1,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(energy_all,nstates,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(energy_err_all,nstates,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      return
      end

      end module
