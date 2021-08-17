      subroutine fatal_error(msg)
      use contrl_file,    only: ounit, errunit
      use mpi

      character msg*(*)

      write(ounit,'(''Fatal error: '',a)') msg
      write(errunit,'(''Fatal error: '',a)') msg
      call mpi_abort(MPI_COMM_WORLD,0,ierr)

      end
