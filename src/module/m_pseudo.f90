module pseudo_mod
    !>Arguments : MPS_L, MPS_QUAD, MPS_GRID, MGAUSS
    integer, parameter :: MPS_L = 5
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
      use multiple_geo, only: MFORCE
      use pseudo_mod, only: MPS_L
      use system,  only: ncent_tot,nctype_tot,nelec

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

module m_pseudo
contains
  subroutine allocate_m_pseudo()
    use pseudo,  only: allocate_pseudo
    call allocate_pseudo()
  end subroutine allocate_m_pseudo

  subroutine deallocate_m_pseudo()
    use pseudo,  only: deallocate_pseudo
    call deallocate_pseudo()
  end subroutine deallocate_m_pseudo
end module
