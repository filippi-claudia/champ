module atom
    !> Arguments: znuc, cent, pecent, iwctype, nctype, ncent
    use precision_kinds, only: dp


    real(dp), dimension(:, :), allocatable :: cent
    real(dp), dimension(:), allocatable :: znuc
    real(dp) :: pecent
    integer, dimension(:), allocatable :: iwctype
    integer :: nctype, ncent
    integer :: nctype_tot, ncent_tot

    private
    public   :: znuc, cent, pecent, iwctype, nctype, ncent, ncent_tot, nctype_tot
    public :: allocate_atom, deallocate_atom
    save
contains
    subroutine allocate_atom()
        use precision_kinds, only: dp

        if (.not. allocated(cent)) allocate (cent(3, ncent_tot))
        if (.not. allocated(znuc)) allocate (znuc(nctype_tot))
        if (.not. allocated(iwctype)) allocate (iwctype(nctype_tot))
    end subroutine allocate_atom

    subroutine deallocate_atom()
        if (allocated(iwctype)) deallocate (iwctype)
        if (allocated(znuc)) deallocate (znuc)
        if (allocated(cent)) deallocate (cent)
    end subroutine deallocate_atom

end module atom

module ghostatom
    !> Arguments: newghostype, nghostcent

    integer :: newghostype
    integer :: nghostcent

    private
    public   :: newghostype, nghostcent
    save
end module ghostatom


