module config
    !> looks a lot like sampling stuff
    !> Arguments: delttn, enew, eold, nearestn, nearesto, pen, peo,
    !> psi2n, psi2o, psido, psijo, rminn, rminno, rmino, rminon, rvminn,
    !> rvminno, rvmino, rvminon, tjfn, tjfo, tjfoo, vnew, vold, xnew, xold
    use force_mod, only: MFORCE
    use precision_kinds, only: dp
    use vmc_mod, only: MELEC
    use mstates_mod, only: MSTATES

    real(dp), dimension(:), allocatable :: delttn !(MELEC)
    real(dp), dimension(:), allocatable :: enew !(MFORCE)
    real(dp), dimension(:, :), allocatable :: eold !(MSTATES, MFORCE)
    integer, dimension(:), allocatable :: nearestn !(MELEC)
    integer, dimension(:), allocatable :: nearesto !(MELEC)
    real(dp) :: pen
    real(dp), dimension(:), allocatable :: peo !(MSTATES)
    real(dp), dimension(:), allocatable :: psi2n !(MFORCE)
    real(dp), dimension(:, :), allocatable :: psi2o !(MSTATES, MFORCE)
    real(dp), dimension(:), allocatable :: psido !(MSTATES)
    real(dp) :: psijo
    real(dp), dimension(:), allocatable :: rminn !(MELEC)
    real(dp), dimension(:), allocatable :: rminno !(MELEC)
    real(dp), dimension(:), allocatable :: rmino !(MELEC)
    real(dp), dimension(:), allocatable :: rminon !(MELEC)
    real(dp), dimension(:, :), allocatable :: rvminn !(3, MELEC)
    real(dp), dimension(:, :), allocatable :: rvminno !(3, MELEC)
    real(dp), dimension(:, :), allocatable :: rvmino !(3, MELEC)
    real(dp), dimension(:, :), allocatable :: rvminon !(3, MELEC)
    real(dp) :: tjfn
    real(dp), dimension(:), allocatable :: tjfo !(MSTATES)
    real(dp) :: tjfoo
    real(dp), dimension(:, :), allocatable :: vnew !(3, MELEC)
    real(dp), dimension(:, :), allocatable :: vold !(3, MELEC)
    real(dp), dimension(:, :), allocatable :: xnew !(3, MELEC)
    real(dp), dimension(:, :), allocatable :: xold !(3, MELEC)

    private
    public   :: delttn, enew, eold, nearestn, nearesto, pen, peo, psi2n
    public   :: psi2o, psido, psijo, rminn, rminno, rmino, rminon, rvminn
    public   :: rvminno, rvmino, rvminon, tjfn, tjfo, tjfoo, vnew, vold, xnew, xold
    public :: allocate_config, deallocate_config
    save
contains
    subroutine allocate_config()
        use const, only: nelec
        use csfs, only: nstates
        use force_mod, only: MFORCE
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC
        use mstates_mod, only: MSTATES
        if (.not. allocated(delttn)) allocate (delttn(nelec))
        if (.not. allocated(enew)) allocate (enew(MFORCE))
        if (.not. allocated(eold)) allocate (eold(MSTATES, MFORCE))
        if (.not. allocated(nearestn)) allocate (nearestn(nelec))
        if (.not. allocated(nearesto)) allocate (nearesto(nelec))
        if (.not. allocated(peo)) allocate (peo(MSTATES))
        if (.not. allocated(psi2n)) allocate (psi2n(MFORCE))
        if (.not. allocated(psi2o)) allocate (psi2o(MSTATES, MFORCE))
        if (.not. allocated(psido)) allocate (psido(MSTATES))
        if (.not. allocated(rminn)) allocate (rminn(nelec))
        if (.not. allocated(rminno)) allocate (rminno(nelec))
        if (.not. allocated(rmino)) allocate (rmino(nelec))
        if (.not. allocated(rminon)) allocate (rminon(nelec))
        if (.not. allocated(rvminn)) allocate (rvminn(3, nelec))
        if (.not. allocated(rvminno)) allocate (rvminno(3, nelec))
        if (.not. allocated(rvmino)) allocate (rvmino(3, nelec))
        if (.not. allocated(rvminon)) allocate (rvminon(3, nelec))
        if (.not. allocated(tjfo)) allocate (tjfo(MSTATES))
        if (.not. allocated(vnew)) allocate (vnew(3, nelec))
        if (.not. allocated(vold)) allocate (vold(3, nelec))
        if (.not. allocated(xnew)) allocate (xnew(3, nelec))
        if (.not. allocated(xold)) allocate (xold(3, nelec))
    end subroutine allocate_config

    subroutine deallocate_config()
        if (allocated(xold)) deallocate (xold)
        if (allocated(xnew)) deallocate (xnew)
        if (allocated(vold)) deallocate (vold)
        if (allocated(vnew)) deallocate (vnew)
        if (allocated(tjfo)) deallocate (tjfo)
        if (allocated(rvminon)) deallocate (rvminon)
        if (allocated(rvmino)) deallocate (rvmino)
        if (allocated(rvminno)) deallocate (rvminno)
        if (allocated(rvminn)) deallocate (rvminn)
        if (allocated(rminon)) deallocate (rminon)
        if (allocated(rmino)) deallocate (rmino)
        if (allocated(rminno)) deallocate (rminno)
        if (allocated(rminn)) deallocate (rminn)
        if (allocated(psido)) deallocate (psido)
        if (allocated(psi2o)) deallocate (psi2o)
        if (allocated(psi2n)) deallocate (psi2n)
        if (allocated(peo)) deallocate (peo)
        if (allocated(nearesto)) deallocate (nearesto)
        if (allocated(nearestn)) deallocate (nearestn)
        if (allocated(eold)) deallocate (eold)
        if (allocated(enew)) deallocate (enew)
        if (allocated(delttn)) deallocate (delttn)
    end subroutine deallocate_config

end module config

module rnyucm
    !> I guess the random number generator
    !> used to move in the MC sampling
    !> Arguments: ll, lm

    integer  :: ll(4)
    integer  :: mm(4)
    data mm/502, 1521, 4071, 2107/
    data ll/0, 0, 0, 1/

    private
    public :: ll, mm
    save

end module rnyucm

module stats
    !> has to do with sampling
    !> Arguments: rejmax
    use precision_kinds, only: dp

    real(dp) :: rejmax

    private
    public :: rejmax
    save
end module stats

module step
    !> I guess has to do with the sampling
    !> Arguments: ekin, ekin2, rprob, suc, trunfb, try
    use precision_kinds, only: dp
    use vmc_mod, only: nrad

    real(dp), dimension(:), allocatable :: ekin !(nrad)
    real(dp), dimension(:), allocatable :: ekin2 !(nrad)
    real(dp), dimension(:), allocatable :: rprob !(nrad)
    real(dp), dimension(:), allocatable :: suc !(nrad)
    real(dp), dimension(:), allocatable :: trunfb !(nrad)
    real(dp), dimension(:), allocatable :: try !(nrad)

    private
    public :: ekin, ekin2, rprob, suc, trunfb, try
    public :: allocate_step, deallocate_step
    save
contains
    subroutine allocate_step()
        use precision_kinds, only: dp
        use vmc_mod, only: nrad
        if (.not. allocated(ekin)) allocate (ekin(nrad))
        if (.not. allocated(ekin2)) allocate (ekin2(nrad))
        if (.not. allocated(rprob)) allocate (rprob(nrad))
        if (.not. allocated(suc)) allocate (suc(nrad))
        if (.not. allocated(trunfb)) allocate (trunfb(nrad))
        if (.not. allocated(try)) allocate (try(nrad))
    end subroutine allocate_step

    subroutine deallocate_step()
        if (allocated(try)) deallocate (try)
        if (allocated(trunfb)) deallocate (trunfb)
        if (allocated(suc)) deallocate (suc)
        if (allocated(rprob)) deallocate (rprob)
        if (allocated(ekin2)) deallocate (ekin2)
        if (allocated(ekin)) deallocate (ekin)
    end subroutine deallocate_step

end module step

module tmpnode
    !> has to do with the sampling
    !> Arguments: distance_node_sum
    use precision_kinds, only: dp

    real(dp) :: distance_node_sum

    private
    public :: distance_node_sum
    save
end module tmpnode

module kinet
    !> kinetic energy ?
    !> only used in metropolis
    !> Arguments: dtdx2n, dtdx2o
    use precision_kinds, only: dp
    use vmc_mod, only: MELEC

    real(dp), dimension(:), allocatable :: dtdx2n !(MELEC)
    real(dp), dimension(:), allocatable :: dtdx2o !(MELEC)

    private
    public :: dtdx2n, dtdx2o
    public :: allocate_kinet, deallocate_kinet
    save
contains
    subroutine allocate_kinet()
        use const, only: nelec
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC
        if (.not. allocated(dtdx2n)) allocate (dtdx2n(nelec))
        if (.not. allocated(dtdx2o)) allocate (dtdx2o(nelec))
    end subroutine allocate_kinet

    subroutine deallocate_kinet()
        if (allocated(dtdx2o)) deallocate (dtdx2o)
        if (allocated(dtdx2n)) deallocate (dtdx2n)
    end subroutine deallocate_kinet

end module kinet

subroutine allocate_m_sampling()
    use config, only: allocate_config
    use step, only: allocate_step
    use kinet, only: allocate_kinet

    call allocate_config()
    call allocate_step()
    call allocate_kinet()
end subroutine allocate_m_sampling
