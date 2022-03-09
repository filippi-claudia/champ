module pseudo_mod
    !>Arguments : MPS_L, MPS_QUAD, MPS_GRID, MGAUSS
    integer, parameter :: MPS_L = 4
    integer, parameter :: MPS_QUAD = 86
    integer, parameter :: MPS_GRID = 1200
    ! for gauss ecp
    integer, parameter :: MGAUSS = 100

    private
    public :: MPS_L, MPS_QUAD, MPS_GRID, MGAUSS
    save
end module pseudo_mod

module pseudo
    !> Arguments: lpot, nloc, vps, vpso
    use precision_kinds, only: dp


    integer, dimension(:), allocatable :: lpot !(MCTYPE)
    integer :: nloc
    real(dp), dimension(:, :, :), allocatable :: vps !(MELEC,MCENT,MPS_L)
    real(dp), dimension(:, :, :, :), allocatable :: vpso !(MELEC,MCENT,MPS_L,MFORCE)

    private
    public :: lpot, nloc, vps, vpso
    public :: allocate_pseudo, deallocate_pseudo
    save
contains
    subroutine allocate_pseudo()
        use const, only: nelec
        use atom, only: nctype_tot
        use atom, only: ncent_tot
        use pseudo_mod, only: MPS_L
        use force_mod, only: MFORCE

        if (.not. allocated(lpot)) allocate (lpot(nctype_tot), source=0)
        if (.not. allocated(vps)) allocate (vps(nelec, ncent_tot, MPS_L))
        if (.not. allocated(vpso)) allocate (vpso(nelec, ncent_tot, MPS_L, MFORCE))
    end subroutine allocate_pseudo

    subroutine deallocate_pseudo()
        if (allocated(vpso)) deallocate (vpso)
        if (allocated(vps)) deallocate (vps)
        if (allocated(lpot)) deallocate (lpot)
    end subroutine deallocate_pseudo

end module pseudo

module pseudo_tm
    !> Arguments: arg, arg_ps, d2pot, nr_ps, r0, r0_ps, rmax, rmax_ps, vpseudo

    use precision_kinds, only: dp


    real(dp), dimension(:), allocatable :: arg !(MCTYPE)
    real(dp), dimension(:), allocatable :: arg_ps !(MCTYPE)
    real(dp), dimension(:, :, :), allocatable :: d2pot !(MPS_GRID,MCTYPE,MPS_L)
    integer, dimension(:), allocatable :: nr_ps !(MCTYPE)
    real(dp), dimension(:), allocatable :: r0 !(MCTYPE)
    real(dp), dimension(:), allocatable :: r0_ps !(MCTYPE)
    real(dp), dimension(:), allocatable :: rmax !(MCTYPE)
    real(dp), dimension(:), allocatable :: rmax_ps !(MCTYPE)
    real(dp), dimension(:, :, :), allocatable :: vpseudo !(MPS_GRID,MCTYPE,MPS_L)

    private
    public :: arg, arg_ps, d2pot, nr_ps, r0, r0_ps, rmax, rmax_ps, vpseudo
    public :: allocate_pseudo_tm, deallocate_pseudo_tm
    save
contains
    subroutine allocate_pseudo_tm()
        use atom, only: nctype_tot
        use pseudo_mod, only: MPS_L, MPS_GRID

        if (.not. allocated(arg)) allocate (arg(nctype_tot))
        if (.not. allocated(arg_ps)) allocate (arg_ps(nctype_tot))
        if (.not. allocated(d2pot)) allocate (d2pot(MPS_GRID, nctype_tot, MPS_L))
        if (.not. allocated(nr_ps)) allocate (nr_ps(nctype_tot), source=0)
        if (.not. allocated(r0)) allocate (r0(nctype_tot))
        if (.not. allocated(r0_ps)) allocate (r0_ps(nctype_tot))
        if (.not. allocated(rmax)) allocate (rmax(nctype_tot))
        if (.not. allocated(rmax_ps)) allocate (rmax_ps(nctype_tot))
        if (.not. allocated(vpseudo)) allocate (vpseudo(MPS_GRID, nctype_tot, MPS_L))
    end subroutine allocate_pseudo_tm

    subroutine deallocate_pseudo_tm()
        if (allocated(vpseudo)) deallocate (vpseudo)
        if (allocated(rmax_ps)) deallocate (rmax_ps)
        if (allocated(rmax)) deallocate (rmax)
        if (allocated(r0_ps)) deallocate (r0_ps)
        if (allocated(r0)) deallocate (r0)
        if (allocated(nr_ps)) deallocate (nr_ps)
        if (allocated(d2pot)) deallocate (d2pot)
        if (allocated(arg_ps)) deallocate (arg_ps)
        if (allocated(arg)) deallocate (arg)
    end subroutine deallocate_pseudo_tm

end module pseudo_tm

module m_pseudo
contains
subroutine allocate_m_pseudo()
    use pseudo, only: allocate_pseudo
    use pseudo_tm, only: allocate_pseudo_tm

    call allocate_pseudo()
    call allocate_pseudo_tm()
end subroutine allocate_m_pseudo

subroutine deallocate_m_pseudo()
    use pseudo, only: deallocate_pseudo
    use pseudo_tm, only: deallocate_pseudo_tm

    call deallocate_pseudo()
    call deallocate_pseudo_tm()
end subroutine deallocate_m_pseudo
end module
