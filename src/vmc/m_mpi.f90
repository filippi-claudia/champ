module mpiconf
    !> Arguments: idtask, nproc, wid, NPROCX
    integer, parameter :: NPROCX = 1524
    integer  :: idtask
    integer  :: nproc
    logical  :: wid

    private
    public :: NPROCX
    public :: idtask, nproc, wid
    public :: mpiconf_init
    save
contains
    subroutine mpiconf_init()
        if(nproc.gt.NPROCX) call fatal_error('MAIN: nproc > NPROCX')
        wid = (idtask .eq. 0)
    end subroutine mpiconf_init
end module mpiconf

module mpiblk
    !> Arguments: iblk_proc

    integer  :: iblk_proc
    private

    public :: iblk_proc
    save
end module mpiblk



module custom_broadcast
!> This module is used to handle the MPI broadcasting of the data
!! read from the external data files using the normal Fortran read
!! commands. A generic procedure interface is used to simplify the
!! broadcasting.
!! @author Ravindra Shinde
!! @date 16 June 2021

    use precision_kinds, only: dp
    use mpi
    implicit none
    integer     :: MPIerror

    private
    private     :: MPIerror
    save

    public       :: bcast

!   Generic interface for all types of variables

    interface bcast
        module procedure bcast_double
        module procedure bcast_real
        module procedure bcast_integer
        module procedure bcast_logical
        module procedure bcast_character
        module procedure bcast_double_1d
        module procedure bcast_real_1d
        module procedure bcast_integer_1d
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

    subroutine bcast_double(scalar)
    !>  Broadcasts a scalar double precision variable from root processor
    !!  to all other processors.
        use mpi
        implicit none
        real(dp), intent(in)    :: scalar               ! scalar to be broadcast

        call MPI_BCAST(scalar, 1, MPI_Double_Precision, 0, MPI_Comm_World, MPIerror)
        call MPI_BARRIER(MPI_Comm_World, MPIerror)
    end subroutine bcast_double



    subroutine bcast_real(scalar)
    !>  Broadcasts a scalar real variable from root processor
    !!  to all other processors.
        use mpi
        implicit none
        real, intent(in)         :: scalar               ! scalar to be broadcast

        call MPI_BCAST(scalar, 1, MPI_Real, 0, MPI_Comm_World, MPIerror)
        call MPI_BARRIER(MPI_Comm_World, MPIerror)
    end subroutine bcast_real



    subroutine bcast_integer(scalar)
    !>  Broadcasts a scalar integer variable from root processor
    !!  to all other processors.
        use mpi
        implicit none
        integer, intent(in)       :: scalar               ! scalar to be broadcast

        call MPI_BCAST(scalar, 1, MPI_INTEGER, 0, MPI_Comm_World, MPIerror)
        call MPI_BARRIER(MPI_Comm_World, MPIerror)
    end subroutine bcast_integer



    subroutine bcast_logical(scalar)
    !>  Broadcasts a scalar logical variable from root processor
    !!  to all other processors.
        use mpi
        implicit none
        logical, intent(in)       :: scalar        ! scalar to be broadcast

        call MPI_BCAST(scalar, 1, MPI_LOGICAL, 0, MPI_Comm_World, MPIerror)
        call MPI_BARRIER(MPI_Comm_World, MPIerror)
    end subroutine bcast_logical



    subroutine bcast_character(scalar)
    !>  Broadcasts a scalar character variable from root processor
    !!  to all other processors.
        use mpi
        implicit none
        character (*), intent(in)     :: scalar    ! scalar to be broadcast
        integer                       :: clength   ! length of character

        clength = len(scalar)

        call MPI_BCAST(scalar, clength, MPI_CHARACTER, 0, MPI_Comm_World, MPIerror)
        call MPI_BARRIER(MPI_Comm_World, MPIerror)
    end subroutine bcast_character



    subroutine bcast_double_1d(array)
    !>  Broadcasts a vector double variable from root processor
    !!  to all other processors.
        use mpi
        implicit none
        real(dp), dimension(:), intent(in) :: array             ! array to be broadcast
        integer                            :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Double_Precision, 0, MPI_Comm_World, MPIerror)
        call MPI_BARRIER(MPI_Comm_World, MPIerror)
    end subroutine bcast_double_1d

    subroutine bcast_real_1d(array)
    !>  Broadcasts a vector single precision variable from root processor
    !!  to all other processors.
        use mpi
        implicit none
        real, dimension(:), intent(in)      :: array             ! array to be broadcast
        integer                             :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Real, 0, MPI_Comm_World, MPIerror)
        call MPI_BARRIER(MPI_Comm_World, MPIerror)
    end subroutine bcast_real_1d

    subroutine bcast_integer_1d(array)
    !>  Broadcasts a vector integer variable from root processor
    !!  to all other processors.
        use mpi
        implicit none
        integer, dimension(:), intent(in)   :: array             ! array to be broadcast
        integer                             :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Integer, 0, MPI_Comm_World, MPIerror)
        call MPI_BARRIER(MPI_Comm_World, MPIerror)
    end subroutine bcast_integer_1d

    subroutine bcast_character_1d(array)
    !>  Broadcasts a vector character variable from root processor
    !!  to all other processors.
        use mpi
        implicit none
        character (*), dimension(:), intent(in)     :: array    ! scalar to be broadcast
        integer                                     :: clength  ! length of character
        integer                                     :: nelements! number of elements

        clength     = len(array(1))
        nelements   = size(array)

        call MPI_BCAST(array, clength*nelements, MPI_CHARACTER, 0, MPI_Comm_World, MPIerror)
        call MPI_BARRIER(MPI_Comm_World, MPIerror)
    end subroutine bcast_character_1d

    subroutine bcast_double_2d(array)
    !>  Broadcasts 2D double precision array from root processor
    !!  to all other processors.
        use mpi
        implicit none
        real(dp), dimension(:,:), intent(in)    :: array             ! array to be broadcast
        integer                                 :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Double_Precision, 0, MPI_Comm_World, MPIerror)
        call MPI_BARRIER(MPI_Comm_World, MPIerror)
    end subroutine bcast_double_2d

    subroutine bcast_real_2d(array)
    !>  Broadcasts 2D single precision array from root processor
    !!  to all other processors.
        use mpi
        implicit none
        real, dimension(:,:), intent(in)        :: array             ! array to be broadcast
        integer                                 :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Real, 0, MPI_Comm_World, MPIerror)
        call MPI_BARRIER(MPI_Comm_World, MPIerror)
    end subroutine bcast_real_2d

    subroutine bcast_integer_2d(array)
    !>  Broadcasts 2D integer array from root processor
    !!  to all other processors.
        use mpi
        implicit none
        integer, dimension(:,:), intent(in)     :: array             ! array to be broadcast
        integer                                 :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Integer, 0, MPI_Comm_World, MPIerror)
        call MPI_BARRIER(MPI_Comm_World, MPIerror)
    end subroutine bcast_integer_2d


    subroutine bcast_double_3d(array)
    !>  Broadcasts 3D double precision array from root processor
    !!  to all other processors.
        use mpi
        implicit none
        real(dp), dimension(:,:,:), intent(in)    :: array             ! array to be broadcast
        integer                                   :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Double_Precision, 0, MPI_Comm_World, MPIerror)
        call MPI_BARRIER(MPI_Comm_World, MPIerror)
    end subroutine bcast_double_3d

    subroutine bcast_real_3d(array)
    !>  Broadcasts 3D single precision array from root processor
    !!  to all other processors.
        use mpi
        implicit none
        real, dimension(:,:,:), intent(in)        :: array             ! array to be broadcast
        integer                                   :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Real, 0, MPI_Comm_World, MPIerror)
        call MPI_BARRIER(MPI_Comm_World, MPIerror)
    end subroutine bcast_real_3d

    subroutine bcast_integer_3d(array)
    !>  Broadcasts 3D integer array from root processor
    !!  to all other processors.
        use mpi
        implicit none
        integer, dimension(:,:,:), intent(in)     :: array             ! array to be broadcast
        integer                                 :: nelements

        nelements = size(array)

        call MPI_BCAST(array, nelements, MPI_Integer, 0, MPI_Comm_World, MPIerror)
        call MPI_BARRIER(MPI_Comm_World, MPIerror)
    end subroutine bcast_integer_3d

    subroutine bcast_double_4d(array)
        !>  Broadcasts 4D double precision array from root processor
        !!  to all other processors.
            use mpi
            implicit none
            real(dp), dimension(:,:,:,:), intent(in)    :: array             ! array to be broadcast
            integer                                     :: nelements

            nelements = size(array)

            call MPI_BCAST(array, nelements, MPI_Double_Precision, 0, MPI_Comm_World, MPIerror)
            call MPI_BARRIER(MPI_Comm_World, MPIerror)
        end subroutine bcast_double_4d

        subroutine bcast_real_4d(array)
        !>  Broadcasts 4D single precision array from root processor
        !!  to all other processors.
            use mpi
            implicit none
            real, dimension(:,:,:,:), intent(in)      :: array             ! array to be broadcast
            integer                                   :: nelements

            nelements = size(array)

            call MPI_BCAST(array, nelements, MPI_Real, 0, MPI_Comm_World, MPIerror)
            call MPI_BARRIER(MPI_Comm_World, MPIerror)
        end subroutine bcast_real_4d

        subroutine bcast_integer_4d(array)
        !>  Broadcasts 4D integer array from root processor
        !!  to all other processors.
            use mpi
            implicit none
            integer, dimension(:,:,:,:), intent(in) :: array             ! array to be broadcast
            integer                                 :: nelements

            nelements = size(array)

            call MPI_BCAST(array, nelements, MPI_Integer, 0, MPI_Comm_World, MPIerror)
            call MPI_BARRIER(MPI_Comm_World, MPIerror)
        end subroutine bcast_integer_4d

  end module custom_broadcast

