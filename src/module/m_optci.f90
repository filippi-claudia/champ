module optci
    !> Arguments:
    !     flags and dimensions for generalized CI expectation values
    !     maximal number of terms, max dim of reduced matrices
    integer :: mxciterm

    integer :: mxcireduced
    integer :: ncimatdim

    private
    public :: mxciterm, mxcireduced, ncimatdim
    public :: set_optci_size
    save
contains
    subroutine set_optci_size()
      use optwf_control, only: method
        if (method .eq. 'linear') then
            mxcireduced = mxciterm
        else
            mxcireduced = 1
        end if
        ncimatdim = mxciterm*(mxcireduced + 1)/2
    end subroutine set_optci_size
end module optci

module ci000
    !> Arguments: iciprt, nciprim, nciterm
    !     iciprt : print flag
    !     nciterm: number of terms (possibly after a transformation)
    !     nciprim: number of primitive terms (determinant ratios)

    integer :: iciprt
    integer :: nciprim
    integer :: nciterm

    private
    public :: iciprt, nciprim, nciterm
    save
end module ci000

module ci001_blk
    !> Arguments: ci_oe, ci_o
      use optci,   only: mxcireduced,mxciterm
      use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: ci_o !(mxciterm)
    real(dp), dimension(:, :), allocatable :: ci_oe !(mxciterm,mxcireduced)

    private
    public :: ci_oe, ci_o
    public :: allocate_ci001_blk, deallocate_ci001_blk
    save
contains
    subroutine allocate_ci001_blk()
      use optci,   only: mxcireduced,mxciterm
        if (.not. allocated(ci_o)) allocate (ci_o(mxciterm))
        if (.not. allocated(ci_oe)) allocate (ci_oe(mxciterm, mxcireduced))
    end subroutine allocate_ci001_blk

    subroutine deallocate_ci001_blk()
        if (allocated(ci_oe)) deallocate (ci_oe)
        if (allocated(ci_o)) deallocate (ci_o)
    end subroutine deallocate_ci001_blk

end module ci001_blk

module ci002_blk
    !> Arguments: ci_o_old, ci_oe_old
      use optci,   only: mxcireduced,mxciterm
      use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: ci_o_old !(mxciterm)
    real(dp), dimension(:, :), allocatable :: ci_oe_old !(mxciterm,mxcireduced)

    private
    public :: ci_o_old, ci_oe_old
    public :: allocate_ci002_blk, deallocate_ci002_blk
    save
contains
    subroutine allocate_ci002_blk()
      use optci,   only: mxcireduced,mxciterm
        if (.not. allocated(ci_o_old)) allocate (ci_o_old(mxciterm))
        if (.not. allocated(ci_oe_old)) allocate (ci_oe_old(mxciterm, mxcireduced))
    end subroutine allocate_ci002_blk

    subroutine deallocate_ci002_blk()
        if (allocated(ci_oe_old)) deallocate (ci_oe_old)
        if (allocated(ci_o_old)) deallocate (ci_o_old)
    end subroutine deallocate_ci002_blk

end module ci002_blk

module ci003_blk
    !> Arguments: ci_e_old, ci_e
      use optci,   only: mxciterm
      use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: ci_e !(mxciterm)
    real(dp), dimension(:), allocatable :: ci_e_old !(mxciterm)

    private
    public :: ci_e_old, ci_e
    public :: allocate_ci003_blk, deallocate_ci003_blk
    save
contains
    subroutine allocate_ci003_blk()
      use optci,   only: mxciterm
        if (.not. allocated(ci_e)) allocate (ci_e(mxciterm))
        if (.not. allocated(ci_e_old)) allocate (ci_e_old(mxciterm))
    end subroutine allocate_ci003_blk

    subroutine deallocate_ci003_blk()
        if (allocated(ci_e_old)) deallocate (ci_e_old)
        if (allocated(ci_e)) deallocate (ci_e)
    end subroutine deallocate_ci003_blk

end module ci003_blk

module ci004_blk
    !> Arguments: ci_de, ci_de_old
      use optci,   only: mxciterm
      use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: ci_de !(mxciterm)
    real(dp), dimension(:), allocatable :: ci_de_old !(mxciterm)

    private
    public :: ci_de, ci_de_old
    public :: allocate_ci004_blk, deallocate_ci004_blk
    save
contains
    subroutine allocate_ci004_blk()
      use optci,   only: mxciterm
        if (.not. allocated(ci_de)) allocate (ci_de(mxciterm))
        if (.not. allocated(ci_de_old)) allocate (ci_de_old(mxciterm))
    end subroutine allocate_ci004_blk

    subroutine deallocate_ci004_blk()
        if (allocated(ci_de_old)) deallocate (ci_de_old)
        if (allocated(ci_de)) deallocate (ci_de)
    end subroutine deallocate_ci004_blk

end module ci004_blk

module ci005_blk
    !> Arguments: ci_o_cum, ci_o_sum
      use optci,   only: mxciterm
      use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: ci_o_cum !(mxciterm)
    real(dp), dimension(:), allocatable :: ci_o_sum !(mxciterm)

    private
    public :: ci_o_cum, ci_o_sum
    public :: allocate_ci005_blk, deallocate_ci005_blk
    save
contains
    subroutine allocate_ci005_blk()
      use optci,   only: mxciterm
        if (.not. allocated(ci_o_cum)) allocate (ci_o_cum(mxciterm))
        if (.not. allocated(ci_o_sum)) allocate (ci_o_sum(mxciterm))
    end subroutine allocate_ci005_blk

    subroutine deallocate_ci005_blk()
        if (allocated(ci_o_sum)) deallocate (ci_o_sum)
        if (allocated(ci_o_cum)) deallocate (ci_o_cum)
    end subroutine deallocate_ci005_blk

end module ci005_blk

module ci006_blk
    !> Arguments: ci_de_cum, ci_de_sum
      use optci,   only: mxciterm
      use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: ci_de_cum !(mxciterm)
    real(dp), dimension(:), allocatable :: ci_de_sum !(mxciterm)

    private
    public :: ci_de_cum, ci_de_sum
    public :: allocate_ci006_blk, deallocate_ci006_blk
    save
contains
    subroutine allocate_ci006_blk()
      use optci,   only: mxciterm
        if (.not. allocated(ci_de_cum)) allocate (ci_de_cum(mxciterm))
        if (.not. allocated(ci_de_sum)) allocate (ci_de_sum(mxciterm))
    end subroutine allocate_ci006_blk

    subroutine deallocate_ci006_blk()
        if (allocated(ci_de_sum)) deallocate (ci_de_sum)
        if (allocated(ci_de_cum)) deallocate (ci_de_cum)
    end subroutine deallocate_ci006_blk

end module ci006_blk

module ci008_blk
    !> Arguments: ci_oe_cm2, ci_oe_sum, ci_oe_cum
      use optci,   only: mxcireduced,mxciterm
      use precision_kinds, only: dp

    real(dp), dimension(:, :), allocatable :: ci_oe_cm2 !(mxciterm,mxcireduced)
    real(dp), dimension(:, :), allocatable :: ci_oe_cum !(mxciterm,mxcireduced)
    real(dp), dimension(:, :), allocatable :: ci_oe_sum !(mxciterm,mxcireduced)

    private
    public :: ci_oe_cm2, ci_oe_sum, ci_oe_cum
    public :: allocate_ci008_blk, deallocate_ci008_blk
    save
contains
    subroutine allocate_ci008_blk()
      use optci,   only: mxcireduced,mxciterm
        if (.not. allocated(ci_oe_cm2)) allocate (ci_oe_cm2(mxciterm, mxcireduced))
        if (.not. allocated(ci_oe_cum)) allocate (ci_oe_cum(mxciterm, mxcireduced))
        if (.not. allocated(ci_oe_sum)) allocate (ci_oe_sum(mxciterm, mxcireduced))
    end subroutine allocate_ci008_blk

    subroutine deallocate_ci008_blk()
        if (allocated(ci_oe_sum)) deallocate (ci_oe_sum)
        if (allocated(ci_oe_cum)) deallocate (ci_oe_cum)
        if (allocated(ci_oe_cm2)) deallocate (ci_oe_cm2)
    end subroutine deallocate_ci008_blk

end module ci008_blk

module ci009_blk
    !> Arguments: ci_oo_sum, ci_oo_cm2, ci_oo_cum
      use optci,   only: ncimatdim
      use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: ci_oo_cm2 !(ncimatdim)
    real(dp), dimension(:), allocatable :: ci_oo_cum !(ncimatdim)
    real(dp), dimension(:), allocatable :: ci_oo_sum !(ncimatdim)

    private
    public :: ci_oo_sum, ci_oo_cm2, ci_oo_cum
    public :: allocate_ci009_blk, deallocate_ci009_blk
    save
contains
    subroutine allocate_ci009_blk()
      use optci,   only: ncimatdim
        if (.not. allocated(ci_oo_cm2)) allocate (ci_oo_cm2(ncimatdim))
        if (.not. allocated(ci_oo_cum)) allocate (ci_oo_cum(ncimatdim))
        if (.not. allocated(ci_oo_sum)) allocate (ci_oo_sum(ncimatdim))
    end subroutine allocate_ci009_blk

    subroutine deallocate_ci009_blk()
        if (allocated(ci_oo_sum)) deallocate (ci_oo_sum)
        if (allocated(ci_oo_cum)) deallocate (ci_oo_cum)
        if (allocated(ci_oo_cm2)) deallocate (ci_oo_cm2)
    end subroutine deallocate_ci009_blk

end module ci009_blk

module ci010_blk
    !> Arguments: ci_ooe_cum, ci_ooe_sum
      use optci,   only: ncimatdim
      use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: ci_ooe_cum !(ncimatdim)
    real(dp), dimension(:), allocatable :: ci_ooe_sum !(ncimatdim)
    private

    public :: ci_ooe_cum, ci_ooe_sum
    public :: allocate_ci010_blk, deallocate_ci010_blk
    save
contains
    subroutine allocate_ci010_blk()
      use optci,   only: ncimatdim
        if (.not. allocated(ci_ooe_cum)) allocate (ci_ooe_cum(ncimatdim))
        if (.not. allocated(ci_ooe_sum)) allocate (ci_ooe_sum(ncimatdim))
    end subroutine allocate_ci010_blk

    subroutine deallocate_ci010_blk()
        if (allocated(ci_ooe_sum)) deallocate (ci_ooe_sum)
        if (allocated(ci_ooe_cum)) deallocate (ci_ooe_cum)
    end subroutine deallocate_ci010_blk

end module ci010_blk

module m_optci
contains
subroutine allocate_m_optci()
      use ci001_blk, only: allocate_ci001_blk
      use ci002_blk, only: allocate_ci002_blk
      use ci003_blk, only: allocate_ci003_blk
      use ci004_blk, only: allocate_ci004_blk
      use ci005_blk, only: allocate_ci005_blk
      use ci006_blk, only: allocate_ci006_blk
      use ci008_blk, only: allocate_ci008_blk
      use ci009_blk, only: allocate_ci009_blk
      use ci010_blk, only: allocate_ci010_blk

    call allocate_ci001_blk()
    call allocate_ci002_blk()
    call allocate_ci003_blk()
    call allocate_ci004_blk()
    call allocate_ci005_blk()
    call allocate_ci006_blk()
    call allocate_ci008_blk()
    call allocate_ci009_blk()
    call allocate_ci010_blk()
end subroutine allocate_m_optci

subroutine deallocate_m_optci()
      use ci001_blk, only: deallocate_ci001_blk
      use ci002_blk, only: deallocate_ci002_blk
      use ci003_blk, only: deallocate_ci003_blk
      use ci004_blk, only: deallocate_ci004_blk
      use ci005_blk, only: deallocate_ci005_blk
      use ci006_blk, only: deallocate_ci006_blk
      use ci008_blk, only: deallocate_ci008_blk
      use ci009_blk, only: deallocate_ci009_blk
      use ci010_blk, only: deallocate_ci010_blk

    call deallocate_ci001_blk()
    call deallocate_ci002_blk()
    call deallocate_ci003_blk()
    call deallocate_ci004_blk()
    call deallocate_ci005_blk()
    call deallocate_ci006_blk()
    call deallocate_ci008_blk()
    call deallocate_ci009_blk()
    call deallocate_ci010_blk()
end subroutine deallocate_m_optci
end module
