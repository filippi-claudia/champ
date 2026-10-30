!> @brief Module for MPI configuration and process identification.
!> @author CHAMP developers
!> @date 2020
!>
!> @details This module manages MPI configuration including process IDs,
!> total number of processes, and master process identification.
module mpiconf

    implicit none

    !> Task ID (rank) of current MPI process.
    integer  :: idtask

    !> Total number of MPI processes.
    integer  :: nproc

    !> MPI error code.
    integer  :: ierr

    !> Flag indicating if current process is master (rank 0).
    logical  :: wid

    private
    public :: idtask, nproc, wid
    public :: mpiconf_init
    save
contains

    !> Initializes MPI configuration flags.
    subroutine mpiconf_init()
        wid = (idtask .eq. 0)
    end subroutine mpiconf_init
end module mpiconf

!> @brief Module for MPI block-level parallelization.
!> @author CHAMP developers
!> @date 2020
!>
!> @details This module manages the assignment of Monte Carlo blocks
!> to different MPI processes for parallel execution.
module mpiblk

    implicit none

    !> Index of block assigned to current processor.
    integer  :: iblk_proc

    private
    public :: iblk_proc
    save
end module mpiblk

!> @brief Module for generic MPI broadcast operations.
!> @author Ravindra Shinde
!> @date June 16, 2021
!>
!> @details This module provides a generic interface for broadcasting data
!> from the root MPI process to all other processes. It handles the MPI
!> broadcasting of data read from external files using normal Fortran read
!> commands. The generic procedure interface simplifies broadcasting by
!> automatically selecting the appropriate type-specific subroutine.
!>
!> Supported data types:
!> - Scalars: double, real, integer, int64, logical, character
!> - 1D arrays: double, real, integer, int64, character
!> - 2D arrays: double, real, integer
!> - 3D arrays: double, real, integer
!> - 4D arrays: double, real, integer
module custom_broadcast
      use mpi
      use precision_kinds, only: dp

    implicit none

    !> MPI error code for broadcast operations.
    integer     :: MPIerror

    private
    private     :: MPIerror
    save

    public       :: bcast

    !> Generic interface for broadcasting variables of any supported type.
    interface bcast
        module procedure bcast_double
        module procedure bcast_real
        module procedure bcast_integer
        module procedure bcast_int64
        module procedure bcast_logical
        module procedure bcast_character
        module procedure bcast_double_1d
        module procedure bcast_real_1d
        module procedure bcast_integer_1d
        module procedure bcast_int64_1d
        module procedure bcast_character_1d
        module procedure bcast_double_2d
        module procedure bcast_real_2d
        module procedure bcast_integer_2d
        module procedure bcast_double_3d
        module procedure bcast_real_3d
        module procedure bcast_integer_3d
        module procedure bcast_double_4d
        module procedure bcast_real_4d
        module procedure bcast_integer_4d
    end interface bcast

    contains

    !> Broadcasts a scalar double precision variable from root to all processes.
    !> @param[in] scalar Double precision scalar to broadcast.
    subroutine bcast_double(scalar)
        implicit none
        real(dp), intent(in)    :: scalar

        call MPI_BCAST(scalar, 1, MPI_Double_Precision, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_double

    !> Broadcasts a scalar single precision variable from root to all processes.
    !> @param[in] scalar Single precision scalar to broadcast.
    subroutine bcast_real(scalar)
        implicit none
        real, intent(in)         :: scalar

        call MPI_BCAST(scalar, 1, MPI_Real, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_real

    !> Broadcasts a scalar integer variable from root to all processes.
    !> @param[in] scalar Integer scalar to broadcast.
    subroutine bcast_integer(scalar)
        implicit none
        integer, intent(in)       :: scalar

        call MPI_BCAST(scalar, 1, MPI_INTEGER, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_integer

    !> Broadcasts a scalar 64-bit integer variable from root to all processes.
    !> @param[in] scalar 64-bit integer scalar to broadcast.
    subroutine bcast_int64(scalar)
        use iso_fortran_env, only: int64
        implicit none
        integer(int64), intent(in)       :: scalar

        call MPI_BCAST(scalar, 1, MPI_INTEGER8, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_int64

    !> Broadcasts a scalar logical variable from root to all processes.
    !> @param[in] scalar Logical scalar to broadcast.
    subroutine bcast_logical(scalar)
        implicit none
        logical, intent(in)       :: scalar

        call MPI_BCAST(scalar, 1, MPI_LOGICAL, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_logical

    !> Broadcasts a scalar character variable from root to all processes.
    !> @param[in] scalar Character string to broadcast.
    subroutine bcast_character(scalar)
        implicit none
        character (*), intent(in)     :: scalar
        integer                       :: clength

        clength = len(scalar)

        call MPI_BCAST(scalar, clength, MPI_CHARACTER, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_character

    !> Broadcasts a 1D double precision array from root to all processes.
    !> @param[in] array 1D double precision array to broadcast.
    subroutine bcast_double_1d(array)
        implicit none
        real(dp), dimension(:), intent(in) :: array
        integer                            :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Double_Precision, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_double_1d

    !> Broadcasts a 1D single precision array from root to all processes.
    !> @param[in] array 1D single precision array to broadcast.
    subroutine bcast_real_1d(array)
        implicit none
        real, dimension(:), intent(in)      :: array
        integer                             :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Real, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_real_1d

    !> Broadcasts a 1D integer array from root to all processes.
    !> @param[in] array 1D integer array to broadcast.
    subroutine bcast_integer_1d(array)
        implicit none
        integer, dimension(:), intent(in)   :: array
        integer                             :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Integer, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_integer_1d

    !> Broadcasts a 1D 64-bit integer array from root to all processes.
    !> @param[in] array 1D 64-bit integer array to broadcast.
    subroutine bcast_int64_1d(array)
        use iso_fortran_env, only: int64
        implicit none
        integer(int64), dimension(:), intent(in)   :: array
        integer                                    :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_INTEGER8, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_int64_1d

    !> Broadcasts a 1D character array from root to all processes.
    !> @param[in] array 1D character array to broadcast.
    subroutine bcast_character_1d(array)
        implicit none
        character (*), dimension(:), intent(in)     :: array
        integer                                     :: clength
        integer                                     :: nelements

        clength     = len(array(1))
        nelements   = size(array)

        call MPI_BCAST(array, clength*nelements, MPI_CHARACTER, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_character_1d

    !> Broadcasts a 2D double precision array from root to all processes.
    !> @param[in] array 2D double precision array to broadcast.
    subroutine bcast_double_2d(array)
        implicit none
        real(dp), dimension(:,:), intent(in)    :: array
        integer                                 :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Double_Precision, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_double_2d

    !> Broadcasts a 2D single precision array from root to all processes.
    !> @param[in] array 2D single precision array to broadcast.
    subroutine bcast_real_2d(array)
        implicit none
        real, dimension(:,:), intent(in)        :: array
        integer                                 :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Real, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_real_2d

    !> Broadcasts a 2D integer array from root to all processes.
    !> @param[in] array 2D integer array to broadcast.
    subroutine bcast_integer_2d(array)
        implicit none
        integer, dimension(:,:), intent(in)     :: array
        integer                                 :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Integer, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_integer_2d

    !> Broadcasts a 3D double precision array from root to all processes.
    !> @param[in] array 3D double precision array to broadcast.
    subroutine bcast_double_3d(array)
        implicit none
        real(dp), dimension(:,:,:), intent(in)    :: array
        integer                                   :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Double_Precision, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_double_3d

    !> Broadcasts a 3D single precision array from root to all processes.
    !> @param[in] array 3D single precision array to broadcast.
    subroutine bcast_real_3d(array)
        implicit none
        real, dimension(:,:,:), intent(in)        :: array
        integer                                   :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Real, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_real_3d

    !> Broadcasts a 3D integer array from root to all processes.
    !> @param[in] array 3D integer array to broadcast.
    subroutine bcast_integer_3d(array)
        implicit none
        integer, dimension(:,:,:), intent(in)     :: array
        integer                                 :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Integer, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_integer_3d

    !> Broadcasts a 4D double precision array from root to all processes.
    !> @param[in] array 4D double precision array to broadcast.
    subroutine bcast_double_4d(array)
        implicit none
        real(dp), dimension(:,:,:,:), intent(in)    :: array
        integer                                     :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Double_Precision, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_double_4d

    !> Broadcasts a 4D single precision array from root to all processes.
    !> @param[in] array 4D single precision array to broadcast.
    subroutine bcast_real_4d(array)
        implicit none
        real, dimension(:,:,:,:), intent(in)      :: array
        integer                                   :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Real, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_real_4d

    !> Broadcasts a 4D integer array from root to all processes.
    !> @param[in] array 4D integer array to broadcast.
    subroutine bcast_integer_4d(array)
        implicit none
        integer, dimension(:,:,:,:), intent(in) :: array
        integer                                 :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Integer, 0, MPI_Comm_World, MPIerror)
    end subroutine bcast_integer_4d

  end module custom_broadcast

