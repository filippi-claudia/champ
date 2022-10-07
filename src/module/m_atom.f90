module system
    !> Arguments: znuc, cent, iwctype, nctype, ncent
    use precision_kinds, only: dp

    implicit none

    integer     :: nelec
    integer     :: ndn
    integer     :: nup
    integer     :: nctype
    integer     :: nctype_tot
    integer     :: ncent
    integer     :: ncent_tot
    integer     :: newghostype
    integer     :: nghostcent

    real(dp)    :: pecent

    real(dp), dimension(:, :), allocatable :: cent
    real(dp), dimension(:), allocatable :: znuc
    integer, dimension(:), allocatable :: iwctype
    character(len=2), dimension(:), allocatable :: symbol
    character(len=2), dimension(:), allocatable :: atomtyp

    private
    public   :: pecent
    public   :: znuc
    public   :: cent
    public   :: iwctype
    public   :: nctype
    public   :: ncent
    public   :: ncent_tot
    public   :: nctype_tot
    public   :: symbol
    public   :: atomtyp
    public   :: newghostype
    public   :: nghostcent
    public   :: allocate_atom, deallocate_atom
    public   :: nelec
    public   :: ndn
    public   :: nup
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


