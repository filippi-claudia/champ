module optorb_mod
    ! flags and dimensions for orbital optimization
    ! maximal number of terms, max dim of reduced matrices

    ! integer, parameter :: MXREDUCED = 1
    ! integer, parameter :: MXMATDIM = MXREDUCED*(MXREDUCED + 1)
    ! integer, parameter :: MXMATDIM2 = MXMATDIM/2

    integer, parameter :: MXORBOP = 8000
    integer :: MXREDUCED
    integer :: MXMATDIM
    integer :: MXMATDIM2

    integer, parameter :: MXREP = 10
    private
    public :: MXORBOP, MXREDUCED, MXMATDIM, MXMATDIM2, MXREP
    public :: set_optorb_size
    save

contains
    subroutine set_optorb_size()

        use method_opt, only: method
        if (method .eq. 'linear') then
            MXREDUCED = MXORBOP
        else
            MXREDUCED = 1
        end if
        MXMATDIM = MXREDUCED*(MXREDUCED + 1)
        MXMATDIM2 = MXMATDIM/2

    end subroutine set_optorb_size
end module optorb_mod

module orb_mat_001
    !> Arguments: orb_o, orb_oe, orb_ho
    use optorb_mod, only: MXORBOP
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES

    real(dp), dimension(:, :), allocatable :: orb_ho !(MXORBOP,MSTATES)
    real(dp), dimension(:, :), allocatable :: orb_o !(MXORBOP,MSTATES)
    real(dp), dimension(:, :), allocatable :: orb_oe !(MXORBOP,MSTATES)

    private
    public :: orb_o, orb_oe, orb_ho
    public :: allocate_orb_mat_001, deallocate_orb_mat_001
    save
contains
    subroutine allocate_orb_mat_001()
        use optorb_mod, only: MXORBOP
        use precision_kinds, only: dp
        use mstates_mod, only: MSTATES
        if (.not. allocated(orb_ho)) allocate (orb_ho(MXORBOP, MSTATES))
        if (.not. allocated(orb_o)) allocate (orb_o(MXORBOP, MSTATES))
        if (.not. allocated(orb_oe)) allocate (orb_oe(MXORBOP, MSTATES))
    end subroutine allocate_orb_mat_001

    subroutine deallocate_orb_mat_001()
        if (allocated(orb_oe)) deallocate (orb_oe)
        if (allocated(orb_o)) deallocate (orb_o)
        if (allocated(orb_ho)) deallocate (orb_ho)
    end subroutine deallocate_orb_mat_001

end module orb_mat_001

module orb_mat_002
    !> Arguments: orb_ho_old, orb_o_old, orb_oe_old
    use optorb_mod, only: MXORBOP
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES

    real(dp), dimension(:, :), allocatable :: orb_ho_old !(MXORBOP,MSTATES)
    real(dp), dimension(:, :), allocatable :: orb_o_old !(MXORBOP,MSTATES)
    real(dp), dimension(:, :), allocatable :: orb_oe_old !(MXORBOP,MSTATES)

    private
    public :: orb_ho_old, orb_o_old, orb_oe_old
    public :: allocate_orb_mat_002, deallocate_orb_mat_002
    save
contains
    subroutine allocate_orb_mat_002()
        use optorb_mod, only: MXORBOP
        use precision_kinds, only: dp
        use mstates_mod, only: MSTATES
        if (.not. allocated(orb_ho_old)) allocate (orb_ho_old(MXORBOP, MSTATES))
        if (.not. allocated(orb_o_old)) allocate (orb_o_old(MXORBOP, MSTATES))
        if (.not. allocated(orb_oe_old)) allocate (orb_oe_old(MXORBOP, MSTATES))
    end subroutine allocate_orb_mat_002

    subroutine deallocate_orb_mat_002()
        if (allocated(orb_oe_old)) deallocate (orb_oe_old)
        if (allocated(orb_o_old)) deallocate (orb_o_old)
        if (allocated(orb_ho_old)) deallocate (orb_ho_old)
    end subroutine deallocate_orb_mat_002

end module orb_mat_002

module orb_mat_003
    !> Arguments: orb_o_sum, orb_o_cum
    use optorb_mod, only: MXORBOP
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES

    real(dp), dimension(:, :), allocatable :: orb_o_cum !(MXORBOP,MSTATES)
    real(dp), dimension(:, :), allocatable :: orb_o_sum !(MXORBOP,MSTATES)

    private
    public :: orb_o_sum, orb_o_cum
    public :: allocate_orb_mat_003, deallocate_orb_mat_003
    save
contains
    subroutine allocate_orb_mat_003()
        use optorb_mod, only: MXORBOP
        use precision_kinds, only: dp
        use mstates_mod, only: MSTATES
        if (.not. allocated(orb_o_cum)) allocate (orb_o_cum(MXORBOP, MSTATES))
        if (.not. allocated(orb_o_sum)) allocate (orb_o_sum(MXORBOP, MSTATES))
    end subroutine allocate_orb_mat_003

    subroutine deallocate_orb_mat_003()
        if (allocated(orb_o_sum)) deallocate (orb_o_sum)
        if (allocated(orb_o_cum)) deallocate (orb_o_cum)
    end subroutine deallocate_orb_mat_003

end module orb_mat_003

module orb_mat_004
    !> Arguments: orb_oe_sum, orb_oe_cum
    use optorb_mod, only: MXORBOP
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES

    real(dp), dimension(:, :), allocatable :: orb_oe_cum !(MXORBOP,MSTATES)
    real(dp), dimension(:, :), allocatable :: orb_oe_sum !(MXORBOP,MSTATES)

    private
    public :: orb_oe_sum, orb_oe_cum
    public :: allocate_orb_mat_004, deallocate_orb_mat_004
    save
contains
    subroutine allocate_orb_mat_004()
        use optorb_mod, only: MXORBOP
        use precision_kinds, only: dp
        use mstates_mod, only: MSTATES
        if (.not. allocated(orb_oe_cum)) allocate (orb_oe_cum(MXORBOP, MSTATES))
        if (.not. allocated(orb_oe_sum)) allocate (orb_oe_sum(MXORBOP, MSTATES))
    end subroutine allocate_orb_mat_004

    subroutine deallocate_orb_mat_004()
        if (allocated(orb_oe_sum)) deallocate (orb_oe_sum)
        if (allocated(orb_oe_cum)) deallocate (orb_oe_cum)
    end subroutine deallocate_orb_mat_004

end module orb_mat_004

module orb_mat_005
    !> Arguments: orb_ho_cum
    use optorb_mod, only: MXORBOP
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES

    real(dp), dimension(:, :), allocatable :: orb_ho_cum !(MXORBOP,MSTATES)

    private
    public :: orb_ho_cum
    public :: allocate_orb_mat_005, deallocate_orb_mat_005
    save
contains
    subroutine allocate_orb_mat_005()
        use optorb_mod, only: MXORBOP
        use precision_kinds, only: dp
        use mstates_mod, only: MSTATES
        if (.not. allocated(orb_ho_cum)) allocate (orb_ho_cum(MXORBOP, MSTATES))
    end subroutine allocate_orb_mat_005

    subroutine deallocate_orb_mat_005()
        if (allocated(orb_ho_cum)) deallocate (orb_ho_cum)
    end subroutine deallocate_orb_mat_005

end module orb_mat_005

module orb_mat_006
    !> Arguments: orb_oo_cum
    use optorb_mod, only: MXMATDIM2
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES

    real(dp), dimension(:, :), allocatable :: orb_oo_cum !(MXMATDIM2,MSTATES)

    private
    public :: orb_oo_cum
    public :: allocate_orb_mat_006, deallocate_orb_mat_006
    save
contains
    subroutine allocate_orb_mat_006()
        use optorb_mod, only: MXMATDIM2
        use precision_kinds, only: dp
        use mstates_mod, only: MSTATES
        if (.not. allocated(orb_oo_cum)) allocate (orb_oo_cum(MXMATDIM2, MSTATES))
    end subroutine allocate_orb_mat_006

    subroutine deallocate_orb_mat_006()
        if (allocated(orb_oo_cum)) deallocate (orb_oo_cum)
    end subroutine deallocate_orb_mat_006

end module orb_mat_006

module orb_mat_007
    !> Arguments: orb_oho_cum
    use optorb_mod, only: MXMATDIM
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES

    real(dp), dimension(:, :), allocatable :: orb_oho_cum !(MXMATDIM,MSTATES)

    private
    public :: orb_oho_cum
    public :: allocate_orb_mat_007, deallocate_orb_mat_007
    save
contains
    subroutine allocate_orb_mat_007()
        use optorb_mod, only: MXMATDIM
        use precision_kinds, only: dp
        use mstates_mod, only: MSTATES
        if (.not. allocated(orb_oho_cum)) allocate (orb_oho_cum(MXMATDIM, MSTATES))
    end subroutine allocate_orb_mat_007

    subroutine deallocate_orb_mat_007()
        if (allocated(orb_oho_cum)) deallocate (orb_oho_cum)
    end subroutine deallocate_orb_mat_007

end module orb_mat_007

module orb_mat_022
    use optorb_mod, only: MXORBOP
    !> Arguments: ideriv

    integer, dimension(:, :), allocatable :: ideriv !(2,MXORBOP)

    private
    public :: ideriv
    public :: allocate_orb_mat_022, deallocate_orb_mat_022
    save
contains
    subroutine allocate_orb_mat_022()
        use optorb_mod, only: MXORBOP
        if (.not. allocated(ideriv)) allocate (ideriv(2, MXORBOP))
    end subroutine allocate_orb_mat_022

    subroutine deallocate_orb_mat_022()
        if (allocated(ideriv)) deallocate (ideriv)
    end subroutine deallocate_orb_mat_022

end module orb_mat_022

module orb_mat_024
    !> Arguments: orb_oe_bsum, orb_f_bcum, orb_e_bsum, orb_w_bsum, orb_o_bsum, orb_f_bcm2
    use optorb_mod, only: MXORBOP
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES

    real(dp), dimension(:), allocatable :: orb_e_bsum !(MSTATES)
    real(dp), dimension(:, :), allocatable :: orb_f_bcm2 !(MXORBOP,MSTATES)
    real(dp), dimension(:, :), allocatable :: orb_f_bcum !(MXORBOP,MSTATES)
    real(dp), dimension(:, :), allocatable :: orb_o_bsum !(MXORBOP,MSTATES)
    real(dp), dimension(:, :), allocatable :: orb_oe_bsum !(MXORBOP,MSTATES)
    real(dp), dimension(:), allocatable :: orb_w_bsum !(MSTATES)

    private
    public :: orb_oe_bsum, orb_f_bcum, orb_e_bsum, orb_w_bsum, orb_o_bsum, orb_f_bcm2
    public :: allocate_orb_mat_024, deallocate_orb_mat_024
    save
contains
    subroutine allocate_orb_mat_024()
        use optorb_mod, only: MXORBOP
        use precision_kinds, only: dp
        use mstates_mod, only: MSTATES
        if (.not. allocated(orb_e_bsum)) allocate (orb_e_bsum(MSTATES))
        if (.not. allocated(orb_f_bcm2)) allocate (orb_f_bcm2(MXORBOP, MSTATES))
        if (.not. allocated(orb_f_bcum)) allocate (orb_f_bcum(MXORBOP, MSTATES))
        if (.not. allocated(orb_o_bsum)) allocate (orb_o_bsum(MXORBOP, MSTATES))
        if (.not. allocated(orb_oe_bsum)) allocate (orb_oe_bsum(MXORBOP, MSTATES))
        if (.not. allocated(orb_w_bsum)) allocate (orb_w_bsum(MSTATES))
    end subroutine allocate_orb_mat_024

    subroutine deallocate_orb_mat_024()
        if (allocated(orb_w_bsum)) deallocate (orb_w_bsum)
        if (allocated(orb_oe_bsum)) deallocate (orb_oe_bsum)
        if (allocated(orb_o_bsum)) deallocate (orb_o_bsum)
        if (allocated(orb_f_bcum)) deallocate (orb_f_bcum)
        if (allocated(orb_f_bcm2)) deallocate (orb_f_bcm2)
        if (allocated(orb_e_bsum)) deallocate (orb_e_bsum)
    end subroutine deallocate_orb_mat_024

end module orb_mat_024

module orb_mat_030
    !> Arguments: orb_wcum, orb_ecum
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES

    real(dp), dimension(:), allocatable :: orb_ecum !(MSTATES)
    real(dp), dimension(:), allocatable :: orb_wcum !(MSTATES)

    private
    public :: orb_wcum, orb_ecum
    public :: allocate_orb_mat_030, deallocate_orb_mat_030
    save
contains
    subroutine allocate_orb_mat_030()
        use precision_kinds, only: dp
        use mstates_mod, only: MSTATES
        if (.not. allocated(orb_ecum)) allocate (orb_ecum(MSTATES))
        if (.not. allocated(orb_wcum)) allocate (orb_wcum(MSTATES))
    end subroutine allocate_orb_mat_030

    subroutine deallocate_orb_mat_030()
        if (allocated(orb_wcum)) deallocate (orb_wcum)
        if (allocated(orb_ecum)) deallocate (orb_ecum)
    end subroutine deallocate_orb_mat_030

end module orb_mat_030

module orb_mat_033
    use optorb_mod, only: MXORBOP
    !> Arguments: irepcol_ref, ideriv_ref, ideriv_iab

    integer, dimension(:), allocatable :: ideriv_iab !(MXORBOP)
    integer, dimension(:, :), allocatable :: ideriv_ref !(MXORBOP,2)
    integer, dimension(:, :), allocatable :: irepcol_ref !(MXORBOP,2)

    private
    public :: irepcol_ref, ideriv_ref, ideriv_iab
    public :: allocate_orb_mat_033, deallocate_orb_mat_033
    save
contains
    subroutine allocate_orb_mat_033()
        use optorb_mod, only: MXORBOP
        if (.not. allocated(ideriv_iab)) allocate (ideriv_iab(MXORBOP))
        if (.not. allocated(ideriv_ref)) allocate (ideriv_ref(MXORBOP, 2))
        if (.not. allocated(irepcol_ref)) allocate (irepcol_ref(MXORBOP, 2))
    end subroutine allocate_orb_mat_033

    subroutine deallocate_orb_mat_033()
        if (allocated(irepcol_ref)) deallocate (irepcol_ref)
        if (allocated(ideriv_ref)) deallocate (ideriv_ref)
        if (allocated(ideriv_iab)) deallocate (ideriv_iab)
    end subroutine deallocate_orb_mat_033

end module orb_mat_033

module optorb
    !> Arguments: dmat_diag, irrep, orb_energy
    use precision_kinds, only: dp
    use vmc_mod, only: MORB

    real(dp), dimension(:), allocatable :: dmat_diag !(MORB)
    integer, dimension(:), allocatable :: irrep !(MORB)
    real(dp), dimension(:), allocatable :: orb_energy !(MORB)

    private
    public :: dmat_diag, irrep, orb_energy
    public :: allocate_optorb, deallocate_optorb
    save
contains
    subroutine allocate_optorb()
        use precision_kinds, only: dp
        use vmc_mod, only: MORB
        if (.not. allocated(dmat_diag)) allocate (dmat_diag(MORB))
        if (.not. allocated(irrep)) allocate (irrep(MORB))
        if (.not. allocated(orb_energy)) allocate (orb_energy(MORB))
    end subroutine allocate_optorb

    subroutine deallocate_optorb()
        if (allocated(orb_energy)) deallocate (orb_energy)
        if (allocated(irrep)) deallocate (irrep)
        if (allocated(dmat_diag)) deallocate (dmat_diag)
    end subroutine deallocate_optorb

end module optorb

module optorb_cblock   ! from optorb.h
    ! norbterm: number of terms (possibly after a transformation)
    ! norbprim: number of primitive terms (determinant ratios)
    integer :: norbterm
    integer :: norbprim

    ! PLT: From old common block /orb004/
    integer :: nefp_blocks
    integer :: nb_current
    integer :: norb_f_bcum

    ! reduced correlation matrix pointers
    ! threshold in terms of std dev. , limit for keeping operators
    ! if iuse_trafo: linearly transformed operators sampled instead of primitive
    !     replacement operators
    ! PLT: From old common blocks /orb006/ and /orb008/.
    integer :: isample_cmat
    integer :: nreduced
    integer :: iuse_trafiuse_trafoo

    ! Dumping block averages for error analysis.
    ! PLT: From old common blocks /orb009/ and /orb010/.
    integer :: idump_blockav
    integer :: iorbsample
    integer :: ns_current

    ! Printing flags:
    integer :: iorbprt
    integer :: iorbprt_sav

    private
    public :: norbterm, norbprim
    public :: nefp_blocks, nb_current, norb_f_bcum
    public :: isample_cmat, nreduced, iuse_trafiuse_trafoo
    public :: idump_blockav, iorbsample, ns_current
    public :: iorbprt, iorbprt_sav
    save
end module optorb_cblock

module optorb_mix
    !> Arguments: iwmix_virt, norbopt, norbvirt
    use vmc_mod, only: MORB

    integer, dimension(:, :), allocatable :: iwmix_virt !(MORB,MORB)
    integer :: norbopt
    integer :: norbvirt

    private
    public :: iwmix_virt, norbopt, norbvirt
    public :: allocate_optorb_mix, deallocate_optorb_mix
    save
contains
    subroutine allocate_optorb_mix()
        use vmc_mod, only: MORB
        if (.not. allocated(iwmix_virt)) allocate (iwmix_virt(MORB, MORB))
    end subroutine allocate_optorb_mix

    subroutine deallocate_optorb_mix()
        if (allocated(iwmix_virt)) deallocate (iwmix_virt)
    end subroutine deallocate_optorb_mix

end module optorb_mix

subroutine allocate_m_optorb()
    use orb_mat_001, only: allocate_orb_mat_001
    use orb_mat_002, only: allocate_orb_mat_002
    use orb_mat_003, only: allocate_orb_mat_003
    use orb_mat_004, only: allocate_orb_mat_004
    use orb_mat_005, only: allocate_orb_mat_005
    use orb_mat_006, only: allocate_orb_mat_006
    use orb_mat_007, only: allocate_orb_mat_007
    use orb_mat_022, only: allocate_orb_mat_022
    use orb_mat_024, only: allocate_orb_mat_024
    use orb_mat_030, only: allocate_orb_mat_030
    use orb_mat_033, only: allocate_orb_mat_033
    use optorb, only: allocate_optorb
    use optorb_mix, only: allocate_optorb_mix

    call allocate_orb_mat_001()
    call allocate_orb_mat_002()
    call allocate_orb_mat_003()
    call allocate_orb_mat_004()
    call allocate_orb_mat_005()
    call allocate_orb_mat_006()
    call allocate_orb_mat_007()
    call allocate_orb_mat_022()
    call allocate_orb_mat_024()
    call allocate_orb_mat_030()
    call allocate_orb_mat_033()
    call allocate_optorb()
    call allocate_optorb_mix()
end subroutine allocate_m_optorb

subroutine deallocate_m_optorb()
    use orb_mat_001, only: deallocate_orb_mat_001
    use orb_mat_002, only: deallocate_orb_mat_002
    use orb_mat_003, only: deallocate_orb_mat_003
    use orb_mat_004, only: deallocate_orb_mat_004
    use orb_mat_005, only: deallocate_orb_mat_005
    use orb_mat_006, only: deallocate_orb_mat_006
    use orb_mat_007, only: deallocate_orb_mat_007
    use orb_mat_022, only: deallocate_orb_mat_022
    use orb_mat_024, only: deallocate_orb_mat_024
    use orb_mat_030, only: deallocate_orb_mat_030
    use orb_mat_033, only: deallocate_orb_mat_033
    use optorb, only: deallocate_optorb
    use optorb_mix, only: deallocate_optorb_mix

    call deallocate_orb_mat_001()
    call deallocate_orb_mat_002()
    call deallocate_orb_mat_003()
    call deallocate_orb_mat_004()
    call deallocate_orb_mat_005()
    call deallocate_orb_mat_006()
    call deallocate_orb_mat_007()
    call deallocate_orb_mat_022()
    call deallocate_orb_mat_024()
    call deallocate_orb_mat_030()
    call deallocate_orb_mat_033()
    call deallocate_optorb()
    call deallocate_optorb_mix()
end subroutine deallocate_m_optorb