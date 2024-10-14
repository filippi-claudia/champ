!> This module provides the walltime (obtained from MPI_Wtime for better precision)
!> @author Ravindra Shinde
!> @date June 23 2021
module mpitimer

    use precision_kinds, only: dp

    !> the starting time logged
    real(dp), save  :: time_start

    !> the time at the last check1
    real(dp), save  :: time_check1

    !> the time at the last check2
    real(dp)        :: time_check2

    !> the final time logged
    real(dp), save  :: time_final

    private
    public :: time_start
    public :: time_check1
    public :: time_check2
    public :: time_final
    public :: time, elapsed_time

contains
    !> Subroutine that uses MPI_Wtime to get the current time
    double precision function time()
        use mpi
        implicit None
        time = MPI_Wtime()
    end function time

    !> Subroutine that logs the elapsed time between two points
    subroutine elapsed_time(message, iter)
        use contrl_file,    only: ounit
        implicit None

        !> @param message :: the message to be printed
        character(len=*), intent(in)    :: message

        !> @param iter :: the iteration number
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
