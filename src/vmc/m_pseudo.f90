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
    use pseudo_mod, only: MPS_L
    use force_mod, only: MFORCE
    use precision_kinds, only: dp
    use vmc_mod, only: MELEC, MCENT, MCTYPE

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
        use forcepar, only: nforce
        use const, only: nelec
        use atom, only: nctype_tot
        use atom, only: ncent_tot
        use pseudo_mod, only: MPS_L
        use force_mod, only: MFORCE
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MCENT, MCTYPE
        if (.not. allocated(lpot)) allocate (lpot(nctype_tot))
        if (.not. allocated(vps)) allocate (vps(nelec, ncent_tot, MPS_L))
        if (.not. allocated(vpso)) allocate (vpso(nelec, ncent_tot, MPS_L, nforce))
    end subroutine allocate_pseudo

    subroutine deallocate_pseudo()
        if (allocated(vpso)) deallocate (vpso)
        if (allocated(vps)) deallocate (vps)
        if (allocated(lpot)) deallocate (lpot)
    end subroutine deallocate_pseudo

end module pseudo

module pseudo_champ
    !> Arguments: igrid_ps, rmax_coul, rmax_nloc
    use precision_kinds, only: dp
    use vmc_mod, only: MCTYPE

    integer, dimension(:), allocatable :: igrid_ps !(MCTYPE)
    real(dp), dimension(:), allocatable :: rmax_coul !(MCTYPE)
    real(dp), dimension(:), allocatable :: rmax_nloc !(MCTYPE)

    private
    public :: igrid_ps, rmax_coul, rmax_nloc
    public :: allocate_pseudo_champ, deallocate_pseudo_champ
    save
contains
    subroutine allocate_pseudo_champ()
        use atom, only: nctype_tot
        use precision_kinds, only: dp
        use vmc_mod, only: MCTYPE
        if (.not. allocated(igrid_ps)) allocate (igrid_ps(nctype_tot))
        if (.not. allocated(rmax_coul)) allocate (rmax_coul(nctype_tot))
        if (.not. allocated(rmax_nloc)) allocate (rmax_nloc(nctype_tot))
    end subroutine allocate_pseudo_champ

    subroutine deallocate_pseudo_champ()
        if (allocated(rmax_nloc)) deallocate (rmax_nloc)
        if (allocated(rmax_coul)) deallocate (rmax_coul)
        if (allocated(igrid_ps)) deallocate (igrid_ps)
    end subroutine deallocate_pseudo_champ

end module pseudo_champ

module pseudo_fahy
    !> Arguments: drad, dradl, nlrad, npotl, potl, ptnlc, rcmax
    use pseudo_mod, only: MPS_L, MPS_GRID
    use precision_kinds, only: dp
    use vmc_mod, only: MCTYPE

    real(dp), dimension(:), allocatable :: drad !(MCTYPE)
    real(dp), dimension(:), allocatable :: dradl !(MCTYPE)
    integer, dimension(:), allocatable :: nlrad !(MCTYPE)
    integer, dimension(:), allocatable :: npotl !(MCTYPE)
    real(dp), dimension(:, :), allocatable :: potl !(MPS_GRID,MCTYPE)
    real(dp), dimension(:, :, :), allocatable :: ptnlc !(MPS_GRID,MCTYPE,MPS_L)
    real(dp), dimension(:), allocatable :: rcmax !(MCTYPE)

    private
    public :: drad, dradl, nlrad, npotl, potl, ptnlc, rcmax
    public :: allocate_pseudo_fahy, deallocate_pseudo_fahy
    save
contains
    subroutine allocate_pseudo_fahy()
        use atom, only: nctype_tot
        use pseudo_mod, only: MPS_L, MPS_GRID
        use precision_kinds, only: dp
        use vmc_mod, only: MCTYPE
        if (.not. allocated(drad)) allocate (drad(nctype_tot))
        if (.not. allocated(dradl)) allocate (dradl(nctype_tot))
        if (.not. allocated(nlrad)) allocate (nlrad(nctype_tot))
        if (.not. allocated(npotl)) allocate (npotl(nctype_tot))
        if (.not. allocated(potl)) allocate (potl(MPS_GRID, nctype_tot))
        if (.not. allocated(ptnlc)) allocate (ptnlc(MPS_GRID, nctype_tot, MPS_L))
        if (.not. allocated(rcmax)) allocate (rcmax(nctype_tot))
    end subroutine allocate_pseudo_fahy

    subroutine deallocate_pseudo_fahy()
        if (allocated(rcmax)) deallocate (rcmax)
        if (allocated(ptnlc)) deallocate (ptnlc)
        if (allocated(potl)) deallocate (potl)
        if (allocated(npotl)) deallocate (npotl)
        if (allocated(nlrad)) deallocate (nlrad)
        if (allocated(dradl)) deallocate (dradl)
        if (allocated(drad)) deallocate (drad)
    end subroutine deallocate_pseudo_fahy

end module pseudo_fahy

module pseudo_tm
    !> Arguments: arg, arg_ps, d2pot, nr_ps, r0, r0_ps, rmax, rmax_ps, vpseudo
    use pseudo_mod, only: MPS_L, MPS_GRID
    use precision_kinds, only: dp
    use vmc_mod, only: MCTYPE

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
        use precision_kinds, only: dp
        use vmc_mod, only: MCTYPE
        if (.not. allocated(arg)) allocate (arg(nctype_tot))
        if (.not. allocated(arg_ps)) allocate (arg_ps(nctype_tot))
        if (.not. allocated(d2pot)) allocate (d2pot(MPS_GRID, nctype_tot, MPS_L))
        if (.not. allocated(nr_ps)) allocate (nr_ps(nctype_tot))
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

subroutine allocate_m_pseudo()
    use pseudo, only: allocate_pseudo
    use pseudo_champ, only: allocate_pseudo_champ
    use pseudo_fahy, only: allocate_pseudo_fahy
    use pseudo_tm, only: allocate_pseudo_tm

    call allocate_pseudo()
    call allocate_pseudo_champ()
    call allocate_pseudo_fahy()
    call allocate_pseudo_tm()
end subroutine allocate_m_pseudo


