      subroutine fatal_error(msg)
      use mpi

      character msg*(*)

      write(6,'(''Fatal error: '',a)') msg
      call mpi_abort(MPI_COMM_WORLD,0,ierr)

      end
