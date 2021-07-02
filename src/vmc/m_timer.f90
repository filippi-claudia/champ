module mpitimer
    !> This module provides the walltime (obtained from MPI_Wtime for better precision)
    !! @author Ravindra Shinde
    !! @date June 23 2021

    use precision_kinds, only: dp
    use mpi


    real(dp), save  :: time_start
    real(dp), save  :: time_check1
    real(dp)        :: time_check2
    real(dp), save  :: time_final

    private
    public :: time_start
    public :: time_check1
    public :: time_check2
    public :: time_final
    public :: time

contains
    double precision function time()
        implicit None
        time = MPI_Wtime()
    end function time
end module mpitimer