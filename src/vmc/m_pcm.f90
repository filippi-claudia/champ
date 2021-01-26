module pcm
    integer, parameter :: MCHS = 1
    integer, parameter :: MCHV = 1
    integer, parameter :: MSPHERE = 30
    private
    public :: MCHS, MCHV, MSPHERE
    save
end module pcm

module pcm_3dgrid
    !     flags and dimensions for the 3d grid objects
    use precision_kinds, only: dp
    integer, parameter :: MGRID_PCM = 1
    integer, parameter :: IUNDEFINED = -1234567890
    integer, parameter :: MGRID_PCM2 = MGRID_PCM*MGRID_PCM
    integer, parameter :: MGRID_PCM3 = MGRID_PCM2*MGRID_PCM
    real(dp), parameter :: UNDEFINED = -1234567890.d0
    real(dp) :: PCM_SHIFT 

    private
    public :: MGRID_PCM, MGRID_PCM2, MGRID_PCM3
    public :: UNDEFINED, IUNDEFINED, PCM_SHIFT
    save
end module pcm_3dgrid

module pcm_ah
    !> Arguments: ahca, bh
    use pcm, only: MCHS
    use precision_kinds, only: dp

    real(dp), dimension(:, :), allocatable :: ahca !(MCHS,MCHS)
    real(dp), dimension(:), allocatable :: bh !(MCHS)

    private
    public :: ahca, bh
    public :: allocate_pcm_ah, deallocate_pcm_ah
    save
contains
    subroutine allocate_pcm_ah()
        use pcm, only: MCHS
        use precision_kinds, only: dp
        if (.not. allocated(ahca)) allocate (ahca(MCHS, MCHS))
        if (.not. allocated(bh)) allocate (bh(MCHS))
    end subroutine allocate_pcm_ah

    subroutine deallocate_pcm_ah()
        if (allocated(bh)) deallocate (bh)
        if (allocated(ahca)) deallocate (ahca)
    end subroutine deallocate_pcm_ah

end module pcm_ah

module pcm_ameta
    !> Arguments: amdlg, eta
    use pcm, only: MCHS
    use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: amdlg !(MCHS)
    real(dp), dimension(:, :), allocatable :: eta !(3,MCHS)

    private
    public :: amdlg, eta
    public :: allocate_pcm_ameta, deallocate_pcm_ameta
    save
contains
    subroutine allocate_pcm_ameta()
        use pcm, only: MCHS
        use precision_kinds, only: dp
        if (.not. allocated(amdlg)) allocate (amdlg(MCHS))
        if (.not. allocated(eta)) allocate (eta(3, MCHS))
    end subroutine allocate_pcm_ameta

    subroutine deallocate_pcm_ameta()
        if (allocated(eta)) deallocate (eta)
        if (allocated(amdlg)) deallocate (amdlg)
    end subroutine deallocate_pcm_ameta

end module pcm_ameta

module pcm_averages
    !> Arguments: spcmsum, spcmcum, spcmcm2, vpcmsum, vpcmcum, vpcmcm2
    ! qopcm_sum, qopcm_cum, qopcm_cm2,
    ! enfpcm_sum(MCHS), enfpcm_cum(MCHS), enfpcm_cm2(MCHS)
    use pcm, only: MCHS
    use precision_kinds, only: dp

    real(dp) :: spcmsum
    real(dp) :: spcmcum
    real(dp) :: spcmcm2
    real(dp) :: vpcmsum
    real(dp) :: vpcmcum
    real(dp) :: vpcmcm2
    real(dp) :: qopcm_sum
    real(dp) :: qopcm_cum
    real(dp) :: qopcm_cm2
    real(dp), dimension(:), allocatable :: enfpcm_sum !(MCHS)
    real(dp), dimension(:), allocatable :: enfpcm_cum !(MCHS)
    real(dp), dimension(:), allocatable :: enfpcm_cm2 !(MCHS)

    private
    public :: spcmsum, spcmcum, spcmcm2, vpcmsum, vpcmcum, vpcmcm2
    public :: qopcm_sum, qopcm_cum, qopcm_cm2
    public :: enfpcm_sum, enfpcm_cum, enfpcm_cm2
    public :: allocate_pcm_averages, deallocate_pcm_averages
    save
contains
    subroutine allocate_pcm_averages()
        use pcm, only: MCHS
        use precision_kinds, only: dp
        if (.not. allocated(enfpcm_sum)) allocate (enfpcm_sum(MCHS))
        if (.not. allocated(enfpcm_cum)) allocate (enfpcm_cum(MCHS))
        if (.not. allocated(enfpcm_cm2)) allocate (enfpcm_cm2(MCHS))
    end subroutine allocate_pcm_averages

    subroutine deallocate_pcm_averages()
        if (allocated(enfpcm_cm2)) deallocate (enfpcm_cm2)
        if (allocated(enfpcm_cum)) deallocate (enfpcm_cum)
        if (allocated(enfpcm_sum)) deallocate (enfpcm_sum)
    end subroutine deallocate_pcm_averages

end module pcm_averages

module pcm_cntrl
    !> Arguments: ichpol, ipcm, ipcmprt, icall, isurf

    integer :: icall
    integer :: ichpol
    integer :: ipcm
    integer :: ipcmprt
    integer :: isurf

    private
    public :: ichpol, ipcm, ipcmprt, icall, isurf
    save
end module pcm_cntrl

module pcm_fdc
    !> Arguments: fs, rcol, feps, rcolt, rcolv, qfree, qvol
    use precision_kinds, only: dp

    real(dp) :: feps
    real(dp) :: fs
    real(dp) :: qfree
    real(dp) :: qvol
    real(dp) :: rcol
    real(dp) :: rcolt
    real(dp) :: rcolv

    private
    public :: fs, rcol, feps, rcolt, rcolv, qfree, qvol
    save
end module pcm_fdc

module pcm_force
    !> Arguments: sch_s
    use pcm, only: MCHS
    use force_mod, only: MFORCE
    use precision_kinds, only: dp

    real(dp), dimension(:, :), allocatable :: sch_s !(MCHS,MFORCE)

    private
    public :: sch_s
    public :: allocate_pcm_force, deallocate_pcm_force
    save
contains
    subroutine allocate_pcm_force()
        use pcm, only: MCHS
        use force_mod, only: MFORCE
        use precision_kinds, only: dp
        if (.not. allocated(sch_s)) allocate (sch_s(MCHS, MFORCE))
    end subroutine allocate_pcm_force

    subroutine deallocate_pcm_force()
        if (allocated(sch_s)) deallocate (sch_s)
    end subroutine deallocate_pcm_force

end module pcm_force

module pcm_grid3d_contrl
    !> Arguments: ipcm_3dgrid

    integer :: ipcm_3dgrid

    private
    public :: ipcm_3dgrid
    save
end module pcm_grid3d_contrl

module pcm_grid3d_array
    !> Arguments: pcm_cart_from_int
    use pcm_3dgrid, only: MGRID_PCM
    use precision_kinds, only: dp

    real(dp), dimension(:, :), allocatable :: pcm_cart_from_int !(MGRID_PCM,3)

    private
    public :: pcm_cart_from_int
    public :: allocate_pcm_grid3d_array, deallocate_pcm_grid3d_array
    save
contains
    subroutine allocate_pcm_grid3d_array()
        use pcm_3dgrid, only: MGRID_PCM
        use precision_kinds, only: dp
        if (.not. allocated(pcm_cart_from_int)) allocate (pcm_cart_from_int(MGRID_PCM, 3))
    end subroutine allocate_pcm_grid3d_array

    subroutine deallocate_pcm_grid3d_array()
        if (allocated(pcm_cart_from_int)) deallocate (pcm_cart_from_int)
    end subroutine deallocate_pcm_grid3d_array

end module pcm_grid3d_array

module pcm_grid3d_param
    !> Arguments: pcm_endpt, pcm_origin, ipcm_nstep3d, pcm_step3d
    use precision_kinds, only: dp

    integer, dimension(:), allocatable :: ipcm_nstep3d !(3)
    real(dp), dimension(:), allocatable :: pcm_endpt !(3)
    real(dp), dimension(:), allocatable :: pcm_origin !(3)
    real(dp), dimension(:), allocatable :: pcm_step3d !(3)

    private
    public :: pcm_endpt, pcm_origin, ipcm_nstep3d, pcm_step3d
    public :: allocate_pcm_grid3d_param, deallocate_pcm_grid3d_param
    save
contains
    subroutine allocate_pcm_grid3d_param()
        use precision_kinds, only: dp
        integer, parameter :: size = 3
        if (.not. allocated(ipcm_nstep3d)) allocate (ipcm_nstep3d(size))
        if (.not. allocated(pcm_endpt)) allocate (pcm_endpt(size))
        if (.not. allocated(pcm_origin)) allocate (pcm_origin(size))
        if (.not. allocated(pcm_step3d)) allocate (pcm_step3d(size))
    end subroutine allocate_pcm_grid3d_param

    subroutine deallocate_pcm_grid3d_param()
        if (allocated(pcm_step3d)) deallocate (pcm_step3d)
        if (allocated(pcm_origin)) deallocate (pcm_origin)
        if (allocated(pcm_endpt)) deallocate (pcm_endpt)
        if (allocated(ipcm_nstep3d)) deallocate (ipcm_nstep3d)
    end subroutine deallocate_pcm_grid3d_param

end module pcm_grid3d_param

module pcm_hpsi
    !> Arguments: enfpcm, pepcms, pepcmv, qopcm
    use pcm, only: MCHS
    use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: enfpcm !(MCHS)
    real(dp) :: pepcms
    real(dp) :: pepcmv
    real(dp) :: qopcm

    private
    public :: enfpcm, pepcms, pepcmv, qopcm
    public :: allocate_pcm_hpsi, deallocate_pcm_hpsi
    save
contains
    subroutine allocate_pcm_hpsi()
        use pcm, only: MCHS
        use precision_kinds, only: dp
        if (.not. allocated(enfpcm)) allocate (enfpcm(MCHS))
    end subroutine allocate_pcm_hpsi

    subroutine deallocate_pcm_hpsi()
        if (allocated(enfpcm)) deallocate (enfpcm)
    end subroutine deallocate_pcm_hpsi

end module pcm_hpsi

module pcm_inda
    use pcm, only: MCHS
    !> Arguments: inda

    integer, dimension(:), allocatable :: inda !(MCHS)

    private
    public :: inda
    public :: allocate_pcm_inda, deallocate_pcm_inda
    save
contains
    subroutine allocate_pcm_inda()
        use pcm, only: MCHS
        if (.not. allocated(inda)) allocate (inda(MCHS))
    end subroutine allocate_pcm_inda

    subroutine deallocate_pcm_inda()
        if (allocated(inda)) deallocate (inda)
    end subroutine deallocate_pcm_inda

end module pcm_inda

module m_pcm_num_spl
    !> Arguments: pcm_num_spl
    use pcm_3dgrid, only: MGRID_PCM
    use precision_kinds, only: dp

    real(dp), dimension(:, :, :, :), allocatable :: pcm_num_spl !(8,MGRID_PCM,MGRID_PCM,MGRID_PCM)

    private
    public :: pcm_num_spl
    public :: allocate_m_pcm_num_spl, deallocate_m_pcm_num_spl
    save
contains
    subroutine allocate_m_pcm_num_spl()
        use pcm_3dgrid, only: MGRID_PCM
        use precision_kinds, only: dp
        if (.not. allocated(pcm_num_spl)) allocate (pcm_num_spl(8, MGRID_PCM, MGRID_PCM, MGRID_PCM))
    end subroutine allocate_m_pcm_num_spl

    subroutine deallocate_m_pcm_num_spl()
        if (allocated(pcm_num_spl)) deallocate (pcm_num_spl)
    end subroutine deallocate_m_pcm_num_spl

end module m_pcm_num_spl

module pcm_num_spl2
    !> Arguments: bc, wk
    use precision_kinds, only: dp

    real(dp) :: bc
    real(dp) :: wk

    private
    public :: bc, wk
    save
end module pcm_num_spl2

module pcm_parms
    !> Arguments: re, nchv, nesph, ze, iscov, eps_solv, xpol,
    !             retk, ch, xe, nvopcm, nch, re2, ncopcm, surk, nscv, nchs, ye, nchs2, nchs1
    use pcm, only: MCHV, MSPHERE
    use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: ch !(MCHV)
    real(dp) :: eps_solv
    integer :: iscov
    integer :: nch
    integer :: nchs
    integer :: nchs1
    integer :: nchs2
    integer :: nchv
    integer :: ncopcm
    integer :: nesph
    integer :: nscv
    integer :: nvopcm
    real(dp), dimension(:), allocatable :: re !(MSPHERE)
    real(dp), dimension(:), allocatable :: re2 !(MSPHERE)
    real(dp) :: retk
    real(dp) :: surk
    real(dp), dimension(:), allocatable :: xe !(MSPHERE)
    real(dp), dimension(:, :), allocatable :: xpol !(3,MCHV)
    real(dp), dimension(:), allocatable :: ye !(MSPHERE)
    real(dp), dimension(:), allocatable :: ze !(MSPHERE)

    private
    public :: re, nchv, nesph, ze, iscov, eps_solv, xpol
    public :: ch, xe, nvopcm, nch, re2, ncopcm, surk
    public :: retk, nscv, nchs, ye, nchs2, nchs1
    public :: allocate_pcm_parms, deallocate_pcm_parms
    save
contains
    subroutine allocate_pcm_parms()
        use pcm, only: MCHV, MSPHERE
        use precision_kinds, only: dp
        if (.not. allocated(ch)) allocate (ch(MCHV))
        if (.not. allocated(re)) allocate (re(MSPHERE))
        if (.not. allocated(re2)) allocate (re2(MSPHERE))
        if (.not. allocated(xe)) allocate (xe(MSPHERE))
        if (.not. allocated(xpol)) allocate (xpol(3, MCHV))
        if (.not. allocated(ye)) allocate (ye(MSPHERE))
        if (.not. allocated(ze)) allocate (ze(MSPHERE))
    end subroutine allocate_pcm_parms

    subroutine deallocate_pcm_parms()
        if (allocated(ze)) deallocate (ze)
        if (allocated(ye)) deallocate (ye)
        if (allocated(xpol)) deallocate (xpol)
        if (allocated(xe)) deallocate (xe)
        if (allocated(re2)) deallocate (re2)
        if (allocated(re)) deallocate (re)
        if (allocated(ch)) deallocate (ch)
    end subroutine deallocate_pcm_parms

end module pcm_parms

module pcm_pot
    !> Arguments: penupol, penups, penupv
    use precision_kinds, only: dp

    real(dp) :: penupol
    real(dp) :: penups
    real(dp) :: penupv
    private

    public :: penupol, penups, penupv
    save
end module pcm_pot

module pcm_xv_new
    !> Arguments: xv_new
    use pcm, only: MCHV
    use precision_kinds, only: dp

    real(dp), dimension(:, :), allocatable :: xv_new !(3,MCHV)

    private
    public :: xv_new
    public :: allocate_pcm_xv_new, deallocate_pcm_xv_new
    save
contains
    subroutine allocate_pcm_xv_new()
        use pcm, only: MCHV
        use precision_kinds, only: dp
        if (.not. allocated(xv_new)) allocate (xv_new(3, MCHV))
    end subroutine allocate_pcm_xv_new

    subroutine deallocate_pcm_xv_new()
        if (allocated(xv_new)) deallocate (xv_new)
    end subroutine deallocate_pcm_xv_new

end module pcm_xv_new

module pcm_unit
    !> Arguments: pcmfile_cavity, pcmfile_chv, pcmfile_chs

    character*80 :: pcmfile_cavity
    character*80 :: pcmfile_chs
    character*80 :: pcmfile_chv

    private
    public :: pcmfile_cavity, pcmfile_chv, pcmfile_chs
    save
end module pcm_unit

module pcmo
    !> Arguments: enfpcmo, qopcmo, spcmo, vpcmo
    use pcm, only: MCHS
    use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: enfpcmo !(MCHS)
    real(dp) :: qopcmo
    real(dp) :: spcmo
    real(dp) :: vpcmo

    private
    public :: enfpcmo, qopcmo, spcmo, vpcmo
    public :: allocate_pcmo, deallocate_pcmo
    save
contains
    subroutine allocate_pcmo()
        use pcm, only: MCHS
        use precision_kinds, only: dp
        if (.not. allocated(enfpcmo)) allocate (enfpcmo(MCHS))
    end subroutine allocate_pcmo

    subroutine deallocate_pcmo()
        if (allocated(enfpcmo)) deallocate (enfpcmo)
    end subroutine deallocate_pcmo

end module pcmo

module spc
    !> Arguments: nsf, num

    integer :: nsf
    integer, dimension(:), allocatable :: num !(50)

    private
    public :: nsf, num
    public :: allocate_spc, deallocate_spc
    save
contains
    subroutine allocate_spc()
        if (.not. allocated(num)) allocate (num(50))
    end subroutine allocate_spc

    subroutine deallocate_spc()
        if (allocated(num)) deallocate (num)
    end subroutine deallocate_spc

end module spc

module spc1
    !> Arguments: csf, qsf, rsf
    use precision_kinds, only: dp

    real(dp), dimension(:, :, :), allocatable :: csf !(750,4,50)
    real(dp), dimension(:, :), allocatable :: qsf !(50,3)
    real(dp), dimension(:), allocatable :: rsf !(50)

    private
    public :: csf, qsf, rsf
    public :: allocate_spc1, deallocate_spc1
    save
contains
    subroutine allocate_spc1()
        use precision_kinds, only: dp
        if (.not. allocated(csf)) allocate (csf(750, 4, 50))
        if (.not. allocated(qsf)) allocate (qsf(50, 3))
        if (.not. allocated(rsf)) allocate (rsf(50))
    end subroutine allocate_spc1

    subroutine deallocate_spc1()
        if (allocated(rsf)) deallocate (rsf)
        if (allocated(qsf)) deallocate (qsf)
        if (allocated(csf)) deallocate (csf)
    end subroutine deallocate_spc1

end module spc1

module spc2
    !> Arguments: nxyz, sfxyz, usf
    use precision_kinds, only: dp

    integer :: nxyz
    real(dp), dimension(:, :), allocatable :: sfxyz !(5000,4)
    real(dp), dimension(:, :), allocatable :: usf !(5000,3)

    private
    public :: nxyz, sfxyz, usf
    public :: allocate_spc2, deallocate_spc2
    save
contains
    subroutine allocate_spc2()
        use precision_kinds, only: dp
        if (.not. allocated(sfxyz)) allocate (sfxyz(5000, 4))
        if (.not. allocated(usf)) allocate (usf(5000, 3))
    end subroutine allocate_spc2

    subroutine deallocate_spc2()
        if (allocated(usf)) deallocate (usf)
        if (allocated(sfxyz)) deallocate (sfxyz)
    end subroutine deallocate_spc2

end module spc2

subroutine allocate_m_pcm()
    use pcm_ah, only: allocate_pcm_ah
    use pcm_ameta, only: allocate_pcm_ameta
    use pcm_averages, only: allocate_pcm_averages
    use pcm_force, only: allocate_pcm_force
    use pcm_grid3d_array, only: allocate_pcm_grid3d_array
    use pcm_grid3d_param, only: allocate_pcm_grid3d_param
    use pcm_hpsi, only: allocate_pcm_hpsi
    use pcm_inda, only: allocate_pcm_inda
    use m_pcm_num_spl, only: allocate_m_pcm_num_spl
    use pcm_parms, only: allocate_pcm_parms
    use pcm_xv_new, only: allocate_pcm_xv_new
    use pcmo, only: allocate_pcmo
    use spc, only: allocate_spc
    use spc1, only: allocate_spc1
    use spc2, only: allocate_spc2

    call allocate_pcm_ah()
    call allocate_pcm_ameta()
    call allocate_pcm_averages()
    call allocate_pcm_force()
    call allocate_pcm_grid3d_array()
    call allocate_pcm_grid3d_param()
    call allocate_pcm_hpsi()
    call allocate_pcm_inda()
    call allocate_m_pcm_num_spl()
    call allocate_pcm_parms()
    call allocate_pcm_xv_new()
    call allocate_pcmo()
    call allocate_spc()
    call allocate_spc1()
    call allocate_spc2()
end subroutine allocate_m_pcm
