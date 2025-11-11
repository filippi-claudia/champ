!> @brief Module for high-precision wall-clock timing using MPI.
!> @author Ravindra Shinde
!> @date June 23, 2021
!>
!> @details This module provides wall-clock timing functionality using MPI_Wtime
!> for better precision than standard Fortran timing routines. It tracks execution
!> time at various checkpoints and can log elapsed time with optional iteration numbers.
!>
!> The module maintains multiple time checkpoints:
!> - Starting time reference
!> - Multiple checkpoint times for interval measurements
!> - Final time for total elapsed time calculation
module mpitimer

    use precision_kinds, only: dp

    implicit none

    !> Starting wall-clock time reference point.
    real(dp), save  :: time_start

    !> First checkpoint time for interval measurements.
    real(dp), save  :: time_check1

    !> Second checkpoint time for interval measurements.
    real(dp)        :: time_check2

    !> Final wall-clock time logged.
    real(dp), save  :: time_final

    private
    public :: time_start
    public :: time_check1
    public :: time_check2
    public :: time_final
    public :: time, elapsed_time

contains

    !> Returns current wall-clock time using MPI_Wtime.
    !> @return Current wall-clock time in seconds.
    double precision function time()
        use mpi
        implicit None
        time = MPI_Wtime()
    end function time

    !> Logs elapsed time between checkpoints with optional message and iteration number.
    !> @param[in] message Descriptive message to be printed with timing information.
    !> @param[in] iter Optional iteration number to include in output.
    subroutine elapsed_time(message, iter)
        use contrl_file,    only: ounit
        implicit None

        character(len=*), intent(in)    :: message
        integer, intent(in), optional   :: iter

        if (present(iter)) then
            time_check2 = time()
            write(ounit, '(a,i4,a,t60,f12.3,a,f12.3,a)') "REAL TIME ELAPSED [",iter,"] :: " // trim(message), time_check2 - time_start, " (sec)", time_check2 - time_check1, " (sec)"
            time_check1 = time_check2
        else
            time_check2 = time()
            write(ounit, '(a,t60, f12.3, a, f12.3, a)') "REAL TIME ELAPSED :: " // trim(message), time_check2 - time_start, " (sec)", time_check2 - time_check1, " (sec)"
            time_check1 = time_check2
        end if

        return
    end subroutine elapsed_time


end module mpitimer
