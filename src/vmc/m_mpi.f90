module mpiconf
    !> Arguments: idtask, nproc, wid
    integer, parameter :: NPROCX = 1024
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
        if (nproc .gt. NPROCX) stop 'nproc>NPROCX in main'
        wid = (idtask .eq. 0)
    end subroutine mpiconf_init
end module mpiconf
