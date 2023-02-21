module error
      private
      public :: fatal_error, trexio_error
      contains
      subroutine fatal_error(msg)
      use contrl_file, only: errunit,ounit
      use mpi
            implicit none
            integer  :: ierr

            character msg*(*)

            write(ounit,'(''Fatal error: '',a)') msg
            write(errunit,'(''Fatal error: '',a)') msg
            error stop "Stopping with error stop"
      end

      subroutine trexio_error(trexio_rc, check_rc, message, filename, line)
            !> This subroutine handles the error in reading/writing with trexio data
            !> @author Ravindra Shinde (r.l.shinde@utwente.nl)
            !> @date 01 June 2022
            !> \param[in] trexio_rc : the return code from the trexio library
            !> \param[in] check_rc  : the return code to compare against trexio_rc
            !> \param[in] message   : the error message for printing
            !> \param[in] filename  : the name of the file where the error occurred
            !> \param[in] line      : the line number where the error occurred


      use contrl_file, only: errunit,ounit
      use mpi,     only: MPI_COMM_WORLD,mpi_abort
            implicit none

            integer, intent(in), value :: trexio_rc
            integer, intent(in), value :: check_rc
            integer, intent(in), value :: line
            character(len=*), intent(in), optional  :: message
            character(len=*), intent(in), optional  :: filename
            integer :: ierr

            if (trexio_rc /= check_rc) then
                  write(ounit,'(a)') "Error reading/writing data from trexio file :: ", trim(message)
                  write(errunit,'(a)') "Error reading/writing data from trexio file :: ", trim(message)
                  write(errunit,'(3a,i6)') "Debug source file :: ", trim(filename), " at line " , line
                  call mpi_abort(MPI_COMM_WORLD,-100,ierr)
            endif

      end subroutine trexio_error

end module
