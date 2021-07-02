module gradhess_all
    !> Arguments: MPARMALL, grad, h, s
    ! use optorb_mod, only: MXREDUCED
    ! use optjas, only: MPARMJ
    ! use optci, only: MXCIREDUCED
    use precision_kinds, only: dp

    implicit none

    ! integer, parameter :: MPARMALL = MPARMJ + MXCIREDUCED + MXREDUCED
    integer :: MPARMALL
    real(dp), dimension(:), allocatable :: grad !(MPARMALL)
    real(dp), dimension(:, :), allocatable :: h !(MPARMALL,MPARMALL)
    real(dp), dimension(:, :), allocatable :: s !(MPARMALL,MPARMALL)

    private
    public :: MPARMALL, grad, h, s
    public :: allocate_gradhess_all, deallocate_gradhess_all, set_gradhess_all_size
    save
contains

    subroutine set_gradhess_all_size()
        use optci, only: MXCIREDUCED
        use optjas, only: MPARMJ
        use optorb_mod, only: MXREDUCED
        MPARMALL = MPARMJ + MXCIREDUCED + MXREDUCED
    end subroutine set_gradhess_all_size

    subroutine allocate_gradhess_all()
        if (.not. allocated(grad)) allocate (grad(MPARMALL))
        if (.not. allocated(h)) allocate (h(MPARMALL, MPARMALL))
        if (.not. allocated(s)) allocate (s(MPARMALL, MPARMALL))
    end subroutine allocate_gradhess_all

    subroutine deallocate_gradhess_all()
        if (allocated(s)) deallocate (s)
        if (allocated(h)) deallocate (h)
        if (allocated(grad)) deallocate (grad)
    end subroutine deallocate_gradhess_all

end module gradhess_all

module gradhessj
    !> Arguments: d2j, d2j_e, de, de_de, de_e, dj, dj_de, dj_dj, dj_dj_e, dj_e, dj_e2, e2
    use optjas, only: MPARMJ
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES

    implicit none

    real(dp), dimension(:, :, :), allocatable :: d2j !(MPARMJ,MPARMJ,MSTATES)
    real(dp), dimension(:, :, :), allocatable :: d2j_e !(MPARMJ,MPARMJ,MSTATES)
    real(dp), dimension(:, :), allocatable :: de !(MPARMJ,MSTATES)
    real(dp), dimension(:, :, :), allocatable :: de_de !(MPARMJ,MPARMJ,MSTATES)
    real(dp), dimension(:, :), allocatable :: de_e !(MPARMJ,MSTATES)
    real(dp), dimension(:, :), allocatable :: dj !(MPARMJ,MSTATES)
    real(dp), dimension(:, :, :), allocatable :: dj_de !(MPARMJ,MPARMJ,MSTATES)
    real(dp), dimension(:, :, :), allocatable :: dj_dj !(MPARMJ,MPARMJ,MSTATES)
    real(dp), dimension(:, :, :), allocatable :: dj_dj_e !(MPARMJ,MPARMJ,MSTATES)
    real(dp), dimension(:, :), allocatable :: dj_e !(MPARMJ,MSTATES)
    real(dp), dimension(:, :), allocatable :: dj_e2 !(MPARMJ,MSTATES)
    real(dp), dimension(:, :), allocatable :: e2 !(MPARMJ,MSTATES)

    private
    public :: d2j, d2j_e, de, de_de, de_e, dj, dj_de, dj_dj, dj_dj_e, dj_e, dj_e2, e2
    public :: allocate_gradhessj, deallocate_gradhessj
    save
contains
    subroutine allocate_gradhessj()
        use optjas, only: MPARMJ
        use mstates_mod, only: MSTATES
        if (.not. allocated(d2j)) allocate (d2j(MPARMJ, MPARMJ, MSTATES))
        if (.not. allocated(d2j_e)) allocate (d2j_e(MPARMJ, MPARMJ, MSTATES))
        if (.not. allocated(de)) allocate (de(MPARMJ, MSTATES))
        if (.not. allocated(de_de)) allocate (de_de(MPARMJ, MPARMJ, MSTATES))
        if (.not. allocated(de_e)) allocate (de_e(MPARMJ, MSTATES))
        if (.not. allocated(dj)) allocate (dj(MPARMJ, MSTATES))
        if (.not. allocated(dj_de)) allocate (dj_de(MPARMJ, MPARMJ, MSTATES))
        if (.not. allocated(dj_dj)) allocate (dj_dj(MPARMJ, MPARMJ, MSTATES))
        if (.not. allocated(dj_dj_e)) allocate (dj_dj_e(MPARMJ, MPARMJ, MSTATES))
        if (.not. allocated(dj_e)) allocate (dj_e(MPARMJ, MSTATES))
        if (.not. allocated(dj_e2)) allocate (dj_e2(MPARMJ, MSTATES))
        if (.not. allocated(e2)) allocate (e2(MPARMJ, MSTATES))
    end subroutine allocate_gradhessj

    subroutine deallocate_gradhessj()
        if (allocated(e2)) deallocate (e2)
        if (allocated(dj_e2)) deallocate (dj_e2)
        if (allocated(dj_e)) deallocate (dj_e)
        if (allocated(dj_dj_e)) deallocate (dj_dj_e)
        if (allocated(dj_dj)) deallocate (dj_dj)
        if (allocated(dj_de)) deallocate (dj_de)
        if (allocated(dj)) deallocate (dj)
        if (allocated(de_e)) deallocate (de_e)
        if (allocated(de_de)) deallocate (de_de)
        if (allocated(de)) deallocate (de)
        if (allocated(d2j_e)) deallocate (d2j_e)
        if (allocated(d2j)) deallocate (d2j)
    end subroutine deallocate_gradhessj

end module gradhessj

module gradhessjo
    !> Arguments: d1d2a_old, d1d2b_old, d2d2a_old, d2d2b_old, denergy_old, gvalue_old
    use optjas, only: MPARMJ
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES

    implicit none

    real(dp), dimension(:), allocatable :: d1d2a_old !(MCTYPE)
    real(dp), dimension(:), allocatable :: d1d2b_old !(2)
    real(dp), dimension(:), allocatable :: d2d2a_old !(MCTYPE)
    real(dp), dimension(:), allocatable :: d2d2b_old !(2)
    real(dp), dimension(:, :), allocatable :: denergy_old !(MPARMJ,MSTATES)
    real(dp), dimension(:), allocatable :: gvalue_old !(MPARMJ)

    private
    public   ::  d1d2a_old, d1d2b_old, d2d2a_old, d2d2b_old, denergy_old, gvalue_old
    public :: allocate_gradhessjo, deallocate_gradhessjo
    save
contains
    subroutine allocate_gradhessjo()
        use atom, only: nctype_tot
        use optjas, only: MPARMJ
        use mstates_mod, only: MSTATES
        if (.not. allocated(d1d2a_old)) allocate (d1d2a_old(nctype_tot))
        if (.not. allocated(d1d2b_old)) allocate (d1d2b_old(2))
        if (.not. allocated(d2d2a_old)) allocate (d2d2a_old(nctype_tot))
        if (.not. allocated(d2d2b_old)) allocate (d2d2b_old(2))
        if (.not. allocated(denergy_old)) allocate (denergy_old(MPARMJ, MSTATES))
        if (.not. allocated(gvalue_old)) allocate (gvalue_old(MPARMJ))
    end subroutine allocate_gradhessjo

    subroutine deallocate_gradhessjo()
        if (allocated(gvalue_old)) deallocate (gvalue_old)
        if (allocated(denergy_old)) deallocate (denergy_old)
        if (allocated(d2d2b_old)) deallocate (d2d2b_old)
        if (allocated(d2d2a_old)) deallocate (d2d2a_old)
        if (allocated(d1d2b_old)) deallocate (d1d2b_old)
        if (allocated(d1d2a_old)) deallocate (d1d2a_old)
    end subroutine deallocate_gradhessjo

end module gradhessjo

module gradhess_ci
    !> Arguments: grad_ci, h_ci, s_ci
    use optci, only: MXCITERM, MXCIREDUCED
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:), allocatable :: grad_ci !(MXCITERM)
    real(dp), dimension(:, :), allocatable :: h_ci !(MXCITERM,MXCIREDUCED)
    real(dp), dimension(:, :), allocatable :: s_ci !(MXCITERM,MXCIREDUCED)

    private
    public   ::  grad_ci, h_ci, s_ci
    public :: allocate_gradhess_ci, deallocate_gradhess_ci
    save
contains
    subroutine allocate_gradhess_ci()
        use optci, only: MXCITERM, MXCIREDUCED
        if (.not. allocated(grad_ci)) allocate (grad_ci(MXCITERM))
        if (.not. allocated(h_ci)) allocate (h_ci(MXCITERM, MXCIREDUCED))
        if (.not. allocated(s_ci)) allocate (s_ci(MXCITERM, MXCIREDUCED))
    end subroutine allocate_gradhess_ci

    subroutine deallocate_gradhess_ci()
        if (allocated(s_ci)) deallocate (s_ci)
        if (allocated(h_ci)) deallocate (h_ci)
        if (allocated(grad_ci)) deallocate (grad_ci)
    end subroutine deallocate_gradhess_ci

end module gradhess_ci

module gradhess_jas
    !> Arguments: grad_jas, h_jas, s_jas
    use optjas, only: MPARMJ
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:), allocatable :: grad_jas !(MPARMJ)
    real(dp), dimension(:, :), allocatable :: h_jas !(MPARMJ,MPARMJ)
    real(dp), dimension(:, :), allocatable :: s_jas !(MPARMJ,MPARMJ)

    private
    public   ::  grad_jas, h_jas, s_jas
    public :: allocate_gradhess_jas, deallocate_gradhess_jas
    save
contains
    subroutine allocate_gradhess_jas()
        use optjas, only: MPARMJ
        if (.not. allocated(grad_jas)) allocate (grad_jas(MPARMJ))
        if (.not. allocated(h_jas)) allocate (h_jas(MPARMJ, MPARMJ))
        if (.not. allocated(s_jas)) allocate (s_jas(MPARMJ, MPARMJ))
    end subroutine allocate_gradhess_jas

    subroutine deallocate_gradhess_jas()
        if (allocated(s_jas)) deallocate (s_jas)
        if (allocated(h_jas)) deallocate (h_jas)
        if (allocated(grad_jas)) deallocate (grad_jas)
    end subroutine deallocate_gradhess_jas

end module gradhess_jas

module gradhess_mix_jas_ci
    !> Arguments: h_mix_jas_ci, s_mix_jas_ci
    use optjas, only: MPARMJ
    use optci, only: MXCITERM
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: h_mix_jas_ci !(2*MPARMJ,MXCITERM)
    real(dp), dimension(:, :), allocatable :: s_mix_jas_ci !(MPARMJ,MXCITERM)

    private
    public   ::  h_mix_jas_ci, s_mix_jas_ci
    public :: allocate_gradhess_mix_jas_ci, deallocate_gradhess_mix_jas_ci
    save
contains
    subroutine allocate_gradhess_mix_jas_ci()
        use optjas, only: MPARMJ
        use optci, only: MXCITERM
        if (.not. allocated(h_mix_jas_ci)) allocate (h_mix_jas_ci(2*MPARMJ, MXCITERM))
        if (.not. allocated(s_mix_jas_ci)) allocate (s_mix_jas_ci(MPARMJ, MXCITERM))
    end subroutine allocate_gradhess_mix_jas_ci

    subroutine deallocate_gradhess_mix_jas_ci()
        if (allocated(s_mix_jas_ci)) deallocate (s_mix_jas_ci)
        if (allocated(h_mix_jas_ci)) deallocate (h_mix_jas_ci)
    end subroutine deallocate_gradhess_mix_jas_ci

end module gradhess_mix_jas_ci

module gradhess_mix_jas_orb
    !> Arguments: h_mix_jas_orb, s_mix_jas_orb
    use optorb_mod, only: MXREDUCED
    use optjas, only: MPARMJ
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: h_mix_jas_orb !(2*MPARMJ,MXREDUCED)
    real(dp), dimension(:, :), allocatable :: s_mix_jas_orb !(MPARMJ,MXREDUCED)

    private
    public   ::  h_mix_jas_orb, s_mix_jas_orb
    public :: allocate_gradhess_mix_jas_orb, deallocate_gradhess_mix_jas_orb
    save
contains
    subroutine allocate_gradhess_mix_jas_orb()
        use optorb_mod, only: MXREDUCED
        use optjas, only: MPARMJ
        if (.not. allocated(h_mix_jas_orb)) allocate (h_mix_jas_orb(2*MPARMJ, MXREDUCED))
        if (.not. allocated(s_mix_jas_orb)) allocate (s_mix_jas_orb(MPARMJ, MXREDUCED))
    end subroutine allocate_gradhess_mix_jas_orb

    subroutine deallocate_gradhess_mix_jas_orb()
        if (allocated(s_mix_jas_orb)) deallocate (s_mix_jas_orb)
        if (allocated(h_mix_jas_orb)) deallocate (h_mix_jas_orb)
    end subroutine deallocate_gradhess_mix_jas_orb

end module gradhess_mix_jas_orb

module gradhess_mix_orb_ci
    !> Arguments: h_mix_ci_orb, s_mix_ci_orb
    use optorb_mod, only: MXREDUCED
    use optci, only: MXCITERM
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: h_mix_ci_orb !(2*MXCITERM,MXREDUCED)
    real(dp), dimension(:, :), allocatable :: s_mix_ci_orb !(MXCITERM,MXREDUCED)

    private
    public   ::  h_mix_ci_orb, s_mix_ci_orb
    public :: allocate_gradhess_mix_orb_ci, deallocate_gradhess_mix_orb_ci
    save
contains
    subroutine allocate_gradhess_mix_orb_ci()
        use optorb_mod, only: MXREDUCED
        use optci, only: MXCITERM
        if (.not. allocated(h_mix_ci_orb)) allocate (h_mix_ci_orb(2*MXCITERM, MXREDUCED))
        if (.not. allocated(s_mix_ci_orb)) allocate (s_mix_ci_orb(MXCITERM, MXREDUCED))
    end subroutine allocate_gradhess_mix_orb_ci

    subroutine deallocate_gradhess_mix_orb_ci()
        if (allocated(s_mix_ci_orb)) deallocate (s_mix_ci_orb)
        if (allocated(h_mix_ci_orb)) deallocate (h_mix_ci_orb)
    end subroutine deallocate_gradhess_mix_orb_ci

end module gradhess_mix_orb_ci

module gradjerr
    !> Arguments: dj_bsum, dj_e_bsum, dj_e_save, dj_save, e_bsum, grad_jas_bcm2, grad_jas_bcum
    use optjas, only: MPARMJ
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES

    implicit none

    real(dp), dimension(:, :), allocatable :: dj_bsum !(MPARMJ,MSTATES)
    real(dp), dimension(:, :), allocatable :: dj_e_bsum !(MPARMJ,MSTATES)
    real(dp), dimension(:, :), allocatable :: dj_e_save !(MPARMJ,MSTATES)
    real(dp), dimension(:, :), allocatable :: dj_save !(MPARMJ,MSTATES)
    real(dp), dimension(:), allocatable :: e_bsum !(MSTATES)
    real(dp), dimension(:, :), allocatable :: grad_jas_bcm2 !(MPARMJ,MSTATES)
    real(dp), dimension(:, :), allocatable :: grad_jas_bcum !(MPARMJ,MSTATES)

    private
    public   ::  dj_bsum, dj_e_bsum, dj_e_save, dj_save, e_bsum, grad_jas_bcm2, grad_jas_bcum
    public :: allocate_gradjerr, deallocate_gradjerr
    save
contains
    subroutine allocate_gradjerr()
        use optjas, only: MPARMJ
        use mstates_mod, only: MSTATES
        if (.not. allocated(dj_bsum)) allocate (dj_bsum(MPARMJ, MSTATES))
        if (.not. allocated(dj_e_bsum)) allocate (dj_e_bsum(MPARMJ, MSTATES))
        if (.not. allocated(dj_e_save)) allocate (dj_e_save(MPARMJ, MSTATES))
        if (.not. allocated(dj_save)) allocate (dj_save(MPARMJ, MSTATES))
        if (.not. allocated(e_bsum)) allocate (e_bsum(MSTATES))
        if (.not. allocated(grad_jas_bcm2)) allocate (grad_jas_bcm2(MPARMJ, MSTATES))
        if (.not. allocated(grad_jas_bcum)) allocate (grad_jas_bcum(MPARMJ, MSTATES))
    end subroutine allocate_gradjerr

    subroutine deallocate_gradjerr()
        if (allocated(grad_jas_bcum)) deallocate (grad_jas_bcum)
        if (allocated(grad_jas_bcm2)) deallocate (grad_jas_bcm2)
        if (allocated(e_bsum)) deallocate (e_bsum)
        if (allocated(dj_save)) deallocate (dj_save)
        if (allocated(dj_e_save)) deallocate (dj_e_save)
        if (allocated(dj_e_bsum)) deallocate (dj_e_bsum)
        if (allocated(dj_bsum)) deallocate (dj_bsum)
    end subroutine deallocate_gradjerr

end module gradjerr

subroutine allocate_m_gradhess()
    use gradhess_all, only: allocate_gradhess_all
    use gradhessj, only: allocate_gradhessj
    use gradhessjo, only: allocate_gradhessjo
    use gradhess_ci, only: allocate_gradhess_ci
    use gradhess_jas, only: allocate_gradhess_jas
    use gradhess_mix_jas_ci, only: allocate_gradhess_mix_jas_ci
    use gradhess_mix_jas_orb, only: allocate_gradhess_mix_jas_orb
    use gradhess_mix_orb_ci, only: allocate_gradhess_mix_orb_ci
    use gradjerr, only: allocate_gradjerr

    implicit none

    call allocate_gradhess_all()
    call allocate_gradhessj()
    call allocate_gradhessjo()
    call allocate_gradhess_ci()
    call allocate_gradhess_jas()
    call allocate_gradhess_mix_jas_ci()
    call allocate_gradhess_mix_jas_orb()
    call allocate_gradhess_mix_orb_ci()
    call allocate_gradjerr()
end subroutine allocate_m_gradhess

subroutine deallocate_m_gradhess()
    use gradhess_all, only: deallocate_gradhess_all
    use gradhessj, only: deallocate_gradhessj
    use gradhessjo, only: deallocate_gradhessjo
    use gradhess_ci, only: deallocate_gradhess_ci
    use gradhess_jas, only: deallocate_gradhess_jas
    use gradhess_mix_jas_ci, only: deallocate_gradhess_mix_jas_ci
    use gradhess_mix_jas_orb, only: deallocate_gradhess_mix_jas_orb
    use gradhess_mix_orb_ci, only: deallocate_gradhess_mix_orb_ci
    use gradjerr, only: deallocate_gradjerr

    implicit none

    call deallocate_gradhess_all()
    call deallocate_gradhessj()
    call deallocate_gradhessjo()
    call deallocate_gradhess_ci()
    call deallocate_gradhess_jas()
    call deallocate_gradhess_mix_jas_ci()
    call deallocate_gradhess_mix_jas_orb()
    call deallocate_gradhess_mix_orb_ci()
    call deallocate_gradjerr()
end subroutine deallocate_m_gradhess
