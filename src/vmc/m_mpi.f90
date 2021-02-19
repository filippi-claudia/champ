module mpiconf
    !> Arguments: idtask, nproc, wid, NPROCX
    integer, parameter :: NPROCX = 1524
    integer  :: idtask
    integer  :: nproc
    logical  :: wid ! true if only one mpi task

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
