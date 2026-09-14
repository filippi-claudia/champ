!> @brief Names of the optional output files written by CHAMP.
!> @details The names collected here are read from the `%module outputs`
!> section of the input file. Giving each run its own set of output file
!> names makes it possible to run, for instance, a VMC and a DMC
!> calculation in the same directory without one overwriting the results
!> of the other.
module outputs

    implicit none

    !> Maximum length of a user-specified output file name
    integer, parameter :: MAX_FILENAME_LENGTH = 80

    !> Name of the file with the analytic forces written by force_analy_fin.
    !> Set with the `file_force_analytic` keyword of the `%module outputs`
    !> input section.
    character(len=MAX_FILENAME_LENGTH) :: file_force_analytic = 'force_analytic'

    private
    public :: MAX_FILENAME_LENGTH
    public :: file_force_analytic
    save

end module outputs
