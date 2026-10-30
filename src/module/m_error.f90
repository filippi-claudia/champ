!> @brief Module for error handling and program termination.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module provides centralized error handling routines for fatal
!> errors and TREXIO library errors. It ensures proper error reporting to output
!> files and graceful parallel program termination using MPI abort.
!>
!> Key subroutines:
!> - fatal_error(): Handles general fatal errors with immediate MPI abort
!> - trexio_error(): Handles TREXIO library-specific errors with detailed diagnostics
!>
!> @note All errors are logged to both standard output and error unit before
!> program termination. Uses MPI_Abort to ensure all parallel processes terminate.
module error
      private
      public :: fatal_error, trexio_error
      contains

      !> Handles fatal errors and terminates the program immediately.
      !>
      !> @details This subroutine prints a fatal error message to both the output
      !> and error units, then terminates all MPI processes using MPI_Abort. This
      !> ensures that all parallel processes are stopped immediately when a critical
      !> error is encountered.
      !>
      !> @param[in] msg Error message to be displayed before termination
      !>
      !> @note This subroutine does not return - it terminates the entire MPI job.
      !> The MPI abort code is set to 0.
      subroutine fatal_error(msg)
      use contrl_file,    only: ounit, errunit
      use mpi
            implicit none
            integer  :: ierr

            character msg*(*)

            write(ounit,'(''Fatal error: '',a)') msg
            write(errunit,'(''Fatal error: '',a)') msg
            call mpi_abort(MPI_COMM_WORLD,0,ierr)

      end

      !> Handles TREXIO library errors with detailed diagnostic information.
      !>
      !> @details This subroutine checks if a TREXIO return code matches an expected
      !> value. If they don't match, it prints detailed error information including
      !> the error message, source filename, and line number where the error occurred,
      !> then terminates all MPI processes. This provides comprehensive debugging
      !> information for TREXIO I/O failures.
      !>
      !> The error information is written to both output and error units. The source
      !> file and line number help developers quickly locate the problematic TREXIO
      !> operation in the code.
      !>
      !> @param[in] trexio_rc Return code from TREXIO library function call
      !> @param[in] check_rc Expected return code to compare against trexio_rc
      !> @param[in] message Optional error message describing the failed operation
      !> @param[in] filename Optional source filename where the error occurred
      !> @param[in] line Line number in the source file where the error occurred
      !>
      !> @note If trexio_rc equals check_rc, the subroutine returns normally.
      !> Otherwise, it terminates the program with MPI abort code -100.
      !>
      !> @author Ravindra Shinde 
      !> @email (r.l.shinde@utwente.nl)
      !> @date 01 June 2022
      subroutine trexio_error(trexio_rc, check_rc, message, filename, line)

      use contrl_file,   only: ounit, errunit
      use mpi,           only: mpi_abort, MPI_COMM_WORLD
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
