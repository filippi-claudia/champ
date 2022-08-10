module system
    !> Arguments: znuc, cent, pecent, iwctype, nctype, ncent
      use precision_kinds, only: dp

    implicit none

    integer  :: nelec
    integer :: ndn
    integer :: nup

    real(dp), dimension(:, :), allocatable :: cent
    real(dp), dimension(:), allocatable :: znuc
    
    integer, dimension(:), allocatable :: iwctype
    integer :: nctype, ncent
    integer :: nctype_tot, ncent_tot
    character(len=2), dimension(:), allocatable :: symbol
    character(len=2), dimension(:), allocatable :: atomtyp

    integer :: newghostype
    integer :: nghostcent

    private
    public :: nelec
    public   :: ndn, nup
    public   :: znuc, cent, iwctype, nctype, ncent, ncent_tot, nctype_tot, symbol, atomtyp
    public   :: allocate_atom, deallocate_atom
    public   :: newghostype, nghostcent
    save

contains
    subroutine allocate_atom()

        if (.not. allocated(cent)) allocate (cent(3, ncent_tot))
        if (.not. allocated(znuc)) allocate (znuc(nctype_tot))
        if (.not. allocated(iwctype)) allocate (iwctype(nctype_tot), source=0)
        if (.not. allocated(symbol)) allocate (symbol(ncent_tot))
    end subroutine allocate_atom

    subroutine deallocate_atom()
        if (allocated(iwctype)) deallocate (iwctype)
        if (allocated(znuc)) deallocate (znuc)
        if (allocated(cent)) deallocate (cent)
    end subroutine deallocate_atom

end module system


