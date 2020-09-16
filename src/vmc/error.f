      subroutine fatal_error(msg)

      character msg*(*)

      include 'mpif.h'

      write(6,'(''Fatal error: '',a)') msg
      call mpi_abort(MPI_COMM_WORLD,0,ierr)

      end
