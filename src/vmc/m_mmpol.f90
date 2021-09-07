module mmpol_mod
    !> Arguments MCHMM, mmpolfile_sites, mmpolfile_chmm

    integer, parameter :: MCHMM = 1
    character*80 mmpolfile_sites
    character*80 mmpolfile_chmm
    private
    public :: MCHMM, mmpolfile_sites, mmpolfile_chmm
    save
end module mmpol_mod

module mmpol_cntrl
    !> Arguments: isites_mmpol, immpolprt, icall_mm, ich_mmpol, immpol

    integer :: icall_mm
    integer :: ich_mmpol
    integer :: immpol
    integer :: immpolprt
    integer :: isites_mmpol

    private
    public :: isites_mmpol, immpolprt, icall_mm, ich_mmpol, immpol
    save
end module mmpol_cntrl

module mmpol_dipol
    !> Arguments: dipo, alfa
    use mmpol_mod, only: MCHMM
    use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: alfa !(MCHMM)
    real(dp), dimension(:, :), allocatable :: dipo !(3,MCHMM)

    private
    public :: dipo, alfa
    public :: allocate_mmpol_dipol, deallocate_mmpol_dipol
    save
contains
    subroutine allocate_mmpol_dipol()
        use mmpol_mod, only: MCHMM
        if (.not. allocated(alfa)) allocate (alfa(MCHMM))
        if (.not. allocated(dipo)) allocate (dipo(3, MCHMM))
    end subroutine allocate_mmpol_dipol

    subroutine deallocate_mmpol_dipol()
        if (allocated(dipo)) deallocate (dipo)
        if (allocated(alfa)) deallocate (alfa)
    end subroutine deallocate_mmpol_dipol

end module mmpol_dipol

module mmpol_hpsi
    !> Arguments: eek_pol, peQMdp, peQMq
    use mmpol_mod, only: MCHMM
    use precision_kinds, only: dp

    real(dp), dimension(:, :), allocatable :: eek_pol !(3,MCHMM)
    real(dp) :: peQMdp
    real(dp) :: peQMq

    private
    public :: eek_pol, peQMdp, peQMq
    public :: allocate_mmpol_hpsi, deallocate_mmpol_hpsi
    save
contains
    subroutine allocate_mmpol_hpsi()
        use mmpol_mod, only: MCHMM
        if (.not. allocated(eek_pol)) allocate (eek_pol(3, MCHMM))
    end subroutine allocate_mmpol_hpsi

    subroutine deallocate_mmpol_hpsi()
        if (allocated(eek_pol)) deallocate (eek_pol)
    end subroutine deallocate_mmpol_hpsi

end module mmpol_hpsi

module mmpol_parms
    !> Arguments: x_mmpol, nchmm, chmm, rqq
    use mmpol_mod, only: MCHMM
    use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: chmm !(MCHMM)
    integer :: nchmm
    real(dp), dimension(:, :), allocatable :: rqq !(MCHMM,MCHMM)
    real(dp), dimension(:, :), allocatable :: x_mmpol !(3,MCHMM)

    private
    public :: x_mmpol, nchmm, chmm, rqq
    public :: allocate_mmpol_parms, deallocate_mmpol_parms
    save
contains
    subroutine allocate_mmpol_parms()
        use mmpol_mod, only: MCHMM
        if (.not. allocated(chmm)) allocate (chmm(MCHMM))
        if (.not. allocated(rqq)) allocate (rqq(MCHMM, MCHMM))
        if (.not. allocated(x_mmpol)) allocate (x_mmpol(3, MCHMM))
    end subroutine allocate_mmpol_parms

    subroutine deallocate_mmpol_parms()
        if (allocated(x_mmpol)) deallocate (x_mmpol)
        if (allocated(rqq)) deallocate (rqq)
        if (allocated(chmm)) deallocate (chmm)
    end subroutine deallocate_mmpol_parms

end module mmpol_parms

module mmpolo
    !> Arguments: cmmpolo, dmmpolo, eeko
    use mmpol_mod, only: MCHMM
    use precision_kinds, only: dp
    use dmc_mod, only: mwalk

    real(dp) :: cmmpolo
    real(dp) :: dmmpolo
    real(dp), dimension(:, :), allocatable :: eeko !(3,MCHMM)

    !> DMC variables:
    real(dp), dimension(:), allocatable :: cmmpolo_dmc !(mwalk)
    real(dp), dimension(:), allocatable :: dmmpolo_dmc !(mwalk)
    real(dp), dimension(:,:), allocatable :: eeko1 !(mwalk,MCHMM)
    real(dp), dimension(:,:), allocatable :: eeko2 !(mwalk,MCHMM)
    real(dp), dimension(:,:), allocatable :: eeko3 !(mwalk,MCHMM)

    private
    public :: cmmpolo, dmmpolo, eeko
    public :: cmmpolo_dmc, dmmpolo_dmc, eeko1, eeko2, eeko3
    public :: allocate_mmpolo, deallocate_mmpolo
    save
contains
    subroutine allocate_mmpolo()
        use mmpol_mod, only: MCHMM
        if (.not. allocated(eeko)) allocate(eeko(3, MCHMM))
        if (.not. allocated(cmmpolo_dmc)) allocate(cmmpolo_dmc(mwalk))
        if (.not. allocated(dmmpolo_dmc)) allocate(dmmpolo_dmc(mwalk))
        if (.not. allocated(eeko1)) allocate(eeko1(mwalk,MCHMM))
        if (.not. allocated(eeko2)) allocate(eeko2(mwalk,MCHMM))
        if (.not. allocated(eeko3)) allocate(eeko3(mwalk,MCHMM))
    end subroutine allocate_mmpolo

    subroutine deallocate_mmpolo()
        if (allocated(eeko)) deallocate(eeko)
        if (allocated(cmmpolo_dmc)) deallocate(cmmpolo_dmc)
        if (allocated(dmmpolo_dmc)) deallocate(dmmpolo_dmc)
        if (allocated(eeko1)) deallocate(eeko1)
        if (allocated(eeko2)) deallocate(eeko2)
        if (allocated(eeko3)) deallocate(eeko3)
    end subroutine deallocate_mmpolo

end module mmpolo

module mmpol_ahpol
    !> Arguments: ah_pol, bh_pol
    use mmpol_mod, only: MCHMM
    use precision_kinds, only: dp

    real(dp), dimension(:, :), allocatable :: ah_pol !(3*MCHMM,3*MCHMM)
    real(dp), dimension(:), allocatable :: bh_pol !(3*MCHMM)

    private
    public :: ah_pol, bh_pol
    public :: allocate_mmpol_ahpol, deallocate_mmpol_ahpol
    save
contains
    subroutine allocate_mmpol_ahpol()
        use mmpol_mod, only: MCHMM
        if (.not. allocated(ah_pol)) allocate (ah_pol(3*MCHMM, 3*MCHMM))
        if (.not. allocated(bh_pol)) allocate (bh_pol(3*MCHMM))
    end subroutine allocate_mmpol_ahpol

    subroutine deallocate_mmpol_ahpol()
        if (allocated(bh_pol)) deallocate (bh_pol)
        if (allocated(ah_pol)) deallocate (ah_pol)
    end subroutine deallocate_mmpol_ahpol

end module mmpol_ahpol

module mmpol_averages
    !> Arguments: cmmpol_cum, cmmpol_cm2, eek2_cum, dmmpol_sum, eek1_cm2, eek_sum, eek2_cm2, cmmpol_sum, dmmpol_cum, dmmpol_cm2, eek3_cum, eek1_cum, eek3_cm2
    use mmpol_mod, only: MCHMM
    use precision_kinds, only: dp

    real(dp) :: cmmpol_cm2
    real(dp) :: cmmpol_cum
    real(dp) :: cmmpol_sum
    real(dp) :: dmmpol_cm2
    real(dp) :: dmmpol_cum
    real(dp) :: dmmpol_sum
    real(dp), dimension(:), allocatable :: eek1_cm2 !(MCHMM)
    real(dp), dimension(:), allocatable :: eek1_cum !(MCHMM)
    real(dp), dimension(:), allocatable :: eek2_cm2 !(MCHMM)
    real(dp), dimension(:), allocatable :: eek2_cum !(MCHMM)
    real(dp), dimension(:), allocatable :: eek3_cm2 !(MCHMM)
    real(dp), dimension(:), allocatable :: eek3_cum !(MCHMM)
    real(dp), dimension(:, :), allocatable :: eek_sum !(3,MCHMM)

    private
    public :: cmmpol_cum, cmmpol_cm2, eek2_cum, dmmpol_sum, eek1_cm2, eek_sum
    public :: eek2_cm2, cmmpol_sum, dmmpol_cum, dmmpol_cm2, eek3_cum, eek1_cum, eek3_cm2
    public :: allocate_mmpol_averages, deallocate_mmpol_averages
    save
contains
    subroutine allocate_mmpol_averages()
        use mmpol_mod, only: MCHMM
        if (.not. allocated(eek1_cm2)) allocate (eek1_cm2(MCHMM))
        if (.not. allocated(eek1_cum)) allocate (eek1_cum(MCHMM))
        if (.not. allocated(eek2_cm2)) allocate (eek2_cm2(MCHMM))
        if (.not. allocated(eek2_cum)) allocate (eek2_cum(MCHMM))
        if (.not. allocated(eek3_cm2)) allocate (eek3_cm2(MCHMM))
        if (.not. allocated(eek3_cum)) allocate (eek3_cum(MCHMM))
        if (.not. allocated(eek_sum)) allocate (eek_sum(3, MCHMM))
    end subroutine allocate_mmpol_averages

    subroutine deallocate_mmpol_averages()
        if (allocated(eek_sum)) deallocate (eek_sum)
        if (allocated(eek3_cum)) deallocate (eek3_cum)
        if (allocated(eek3_cm2)) deallocate (eek3_cm2)
        if (allocated(eek2_cum)) deallocate (eek2_cum)
        if (allocated(eek2_cm2)) deallocate (eek2_cm2)
        if (allocated(eek1_cum)) deallocate (eek1_cum)
        if (allocated(eek1_cm2)) deallocate (eek1_cm2)
    end subroutine deallocate_mmpol_averages

end module mmpol_averages

module mmpol_fdc
    !> Arguments: a_cutoff, rcolm, screen1, screen2
    use mmpol_mod, only: MCHMM
    use precision_kinds, only: dp

    real(dp) :: a_cutoff
    real(dp) :: rcolm
    real(dp), dimension(:, :), allocatable :: screen1 !(MCHMM,MCHMM)
    real(dp), dimension(:, :), allocatable :: screen2 !(MCHMM,MCHMM)

    private
    public :: a_cutoff, rcolm, screen1, screen2
    public :: allocate_mmpol_fdc, deallocate_mmpol_fdc
    save
contains
    subroutine allocate_mmpol_fdc()
        use mmpol_mod, only: MCHMM
        if (.not. allocated(screen1)) allocate (screen1(MCHMM, MCHMM))
        if (.not. allocated(screen2)) allocate (screen2(MCHMM, MCHMM))
    end subroutine allocate_mmpol_fdc

    subroutine deallocate_mmpol_fdc()
        if (allocated(screen2)) deallocate (screen2)
        if (allocated(screen1)) deallocate (screen1)
    end subroutine deallocate_mmpol_fdc

end module mmpol_fdc

module mmpol_field
    !> Arguments: eqk_pol, enk_pol
    use mmpol_mod, only: MCHMM
    use precision_kinds, only: dp

    real(dp), dimension(:, :), allocatable :: enk_pol !(3,MCHMM)
    real(dp), dimension(:, :), allocatable :: eqk_pol !(3,MCHMM)

    private
    public :: eqk_pol, enk_pol
    public :: allocate_mmpol_field, deallocate_mmpol_field
    save
contains
    subroutine allocate_mmpol_field()
        use mmpol_mod, only: MCHMM
        if (.not. allocated(enk_pol)) allocate (enk_pol(3, MCHMM))
        if (.not. allocated(eqk_pol)) allocate (eqk_pol(3, MCHMM))
    end subroutine allocate_mmpol_field

    subroutine deallocate_mmpol_field()
        if (allocated(eqk_pol)) deallocate (eqk_pol)
        if (allocated(enk_pol)) deallocate (enk_pol)
    end subroutine deallocate_mmpol_field

end module mmpol_field

module mmpol_inds
    use mmpol_mod, only: MCHMM
    !> Arguments: inds_pol

    integer, dimension(:), allocatable :: inds_pol !(MCHMM)

    private
    public :: inds_pol
    public :: allocate_mmpol_inds, deallocate_mmpol_inds
    save
contains
    subroutine allocate_mmpol_inds()
        use mmpol_mod, only: MCHMM
        if (.not. allocated(inds_pol)) allocate (inds_pol(MCHMM))
    end subroutine allocate_mmpol_inds

    subroutine deallocate_mmpol_inds()
        if (allocated(inds_pol)) deallocate (inds_pol)
    end subroutine deallocate_mmpol_inds

end module mmpol_inds

module mmpol_pot
    !> Arguments: peqq, pepol_dp, pepol_q, penu_q, peq_dp, penu_dp, u_dd, u_self
    use precision_kinds, only: dp

    real(dp) :: penu_dp
    real(dp) :: penu_q
    real(dp) :: pepol_dp
    real(dp) :: pepol_q
    real(dp) :: peq_dp
    real(dp) :: peqq
    real(dp) :: u_dd
    real(dp) :: u_self

    private
    public :: peqq, pepol_dp, pepol_q, penu_q, peq_dp, penu_dp, u_dd, u_self
    save
end module mmpol_pot

subroutine allocate_m_mmpol()
    use mmpol_dipol, only: allocate_mmpol_dipol
    use mmpol_hpsi, only: allocate_mmpol_hpsi
    use mmpol_parms, only: allocate_mmpol_parms
    use mmpolo, only: allocate_mmpolo
    use mmpol_ahpol, only: allocate_mmpol_ahpol
    use mmpol_averages, only: allocate_mmpol_averages
    use mmpol_fdc, only: allocate_mmpol_fdc
    use mmpol_field, only: allocate_mmpol_field
    use mmpol_inds, only: allocate_mmpol_inds

    call allocate_mmpol_dipol()
    call allocate_mmpol_hpsi()
    call allocate_mmpol_parms()
    call allocate_mmpolo()
    call allocate_mmpol_ahpol()
    call allocate_mmpol_averages()
    call allocate_mmpol_fdc()
    call allocate_mmpol_field()
    call allocate_mmpol_inds()
end subroutine allocate_m_mmpol

subroutine deallocate_m_mmpol()
    use mmpol_dipol, only: deallocate_mmpol_dipol
    use mmpol_hpsi, only: deallocate_mmpol_hpsi
    use mmpol_parms, only: deallocate_mmpol_parms
    use mmpolo, only: deallocate_mmpolo
    use mmpol_ahpol, only: deallocate_mmpol_ahpol
    use mmpol_averages, only: deallocate_mmpol_averages
    use mmpol_fdc, only: deallocate_mmpol_fdc
    use mmpol_field, only: deallocate_mmpol_field
    use mmpol_inds, only: deallocate_mmpol_inds

    call deallocate_mmpol_dipol()
    call deallocate_mmpol_hpsi()
    call deallocate_mmpol_parms()
    call deallocate_mmpolo()
    call deallocate_mmpol_ahpol()
    call deallocate_mmpol_averages()
    call deallocate_mmpol_fdc()
    call deallocate_mmpol_field()
    call deallocate_mmpol_inds()
end subroutine deallocate_m_mmpol
