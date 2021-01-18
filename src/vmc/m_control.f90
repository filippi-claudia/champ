module contrl
    !> Arguments: idump, irstar, isite, nconf, nblk, nblk_max, nblkeq, nconf_new, nstep,
    !> icharged_atom

    integer :: idump
    integer :: irstar
    integer :: isite
    integer :: nconf
    integer :: nblk
    integer :: nblkeq
    integer :: nblk_max
    integer :: nconf_new
    integer :: nstep
    integer :: icharged_atom 

    private
    public :: idump, irstar, isite, nconf, nblk, nblkeq, nconf_new, nstep
    public :: icharged_atom
    public :: nblk_max
    save
end module contrl

module contr2
    !> Arguments: i3body, ianalyt_lap, iaver, icusp, icusp2, ifock, ijas, irewgt, isc, istrch

    integer :: i3body
    integer :: ianalyt_lap
    integer :: iaver
    integer :: icusp
    integer :: icusp2
    integer :: ifock
    integer :: ijas
    integer :: irewgt
    integer :: isc
    integer :: istrch

    private
    public :: i3body, ianalyt_lap, iaver, icusp, icusp2, ifock, ijas, irewgt, isc, istrch
    save
end module contr2

module contr3
    !> Arguments: mode

    character*12 :: mode

    private
    public :: mode
    public :: init_control_mode
    save
contains
    subroutine init_control_mode(str_mode)
        implicit None
        character(12), intent(IN) :: str_mode
        mode = str_mode
    end subroutine init_control_mode

end module contr3

module contrl_per
    !> Arguments: iperiodic, ibasis

    integer :: iperiodic, ibasis

    private
    public   :: iperiodic, ibasis
    save
end module contrl_per

module contrldmc
    !> Arguments: iacc_rej, icross, icuspg, icut_br, icut_e, idiv_v, idmc, ipq, itau_eff, nfprod, rttau, tau, taueff, tautot
    use force_mod, only: MFORCE
    use precision_kinds, only: dp

    integer :: iacc_rej
    integer :: icross
    integer :: icuspg
    integer :: icut_br
    integer :: icut_e
    integer :: idiv_v
    integer :: idmc
    integer :: ipq
    integer :: itau_eff
    integer :: nfprod
    real(dp) :: rttau
    real(dp) :: tau
    real(dp), dimension(:), allocatable :: taueff !(MFORCE)
    real(dp) :: tautot

    private
    public :: iacc_rej, icross, icuspg, icut_br, icut_e, idiv_v, idmc, ipq, itau_eff, nfprod, rttau, tau, taueff, tautot
    public :: allocate_contrldmc, deallocate_contrldmc
    save
contains
    subroutine allocate_contrldmc()
        use force_mod, only: MFORCE
        use precision_kinds, only: dp
        if (.not. allocated(taueff)) allocate (taueff(MFORCE))
    end subroutine allocate_contrldmc

    subroutine deallocate_contrldmc()
        if (allocated(taueff)) deallocate (taueff)
    end subroutine deallocate_contrldmc

end module contrldmc

module contrl_file

    implicit None
    character(20) :: log_filename
    character(20) :: proc_filename

    private
    public :: log_filename, proc_filename
    public :: close_files
    public :: init_procfile, init_logfile
    save
contains

! Open all log/output files at once does not work because init_procfile
! needs to know the value of ipr flag from read_input.f.
!    subroutine init_files()
!        call init_logfile()
!        call init_procfile()
!    end subroutine init_files

    subroutine close_files()
        close (5)
        close (6)
        close (45)
    end subroutine close_files

    subroutine init_logfile()
        use mpiconf, only: wid

        !> Open the standard output and the log file only on the master
        if (wid) then
            log_filename = 'output.log'
        else
            log_filename = '/dev/null'
            close (6)
            open (6, file='/dev/null')
        endif
        open (45, file=log_filename, status='unknown')
    end subroutine init_logfile

    subroutine init_procfile()
        use mpiconf, only: idtask
        use const, only: ipr

        if (ipr .gt. 1) then
            if (idtask .lt. 10) then
                write (proc_filename, '(i1)') idtask
            elseif (idtask .lt. 100) then
                write (proc_filename, '(i2)') idtask
            elseif (idtask .lt. 1000) then
                write (proc_filename, '(i3)') idtask
            else
                write (proc_filename, '(i4)') idtask
            endif
            proc_filename = 'check.'//proc_filename(1:index(proc_filename, ' ') - 1)
            open (unit=88, form='formatted', file=proc_filename)
        endif
    end subroutine init_procfile

end module contrl_file

subroutine allocate_m_control()
    use contrldmc, only: allocate_contrldmc
    call allocate_contrldmc()
end subroutine allocate_m_control
