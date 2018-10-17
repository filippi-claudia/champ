      subroutine fatal_error(msg)

      implicit double precision (a-h,o-z)
      character msg*(*)

      include 'mpi_qmc.h'
      include 'mpif.h'

      write(6,'(''Fatal error: '',a)') msg
      call mpi_abort(MPI_COMM_WORLD,0,ierr)

      return
      end
