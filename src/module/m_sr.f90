module sr_mod
    !> Arguments:
    integer :: mparm
    integer :: mobs
    integer :: mconf

    private
    public :: mparm, mobs, mconf
    save
end module sr_mod

module sr_index
    !> Arguments: jelo, jelo2, jelohfj

    integer :: jelo
    integer :: jelo2
    integer :: jelohfj

    private
    public :: jelo, jelo2, jelohfj
    save
end module sr_index

module sr_mat_n
    !> Arguments: elocal, h_sr, jefj, jfj, jhfj, nconf_n, obs, s_diag, s_ii_inv, sr_ho, sr_o, wtg, obs_tot
    use sr_mod, only: mparm, mobs, mconf
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES

    real(dp), dimension(:, :), allocatable :: elocal !(mconf,MSTATES)
    real(dp), dimension(:), allocatable :: h_sr !(mparm)
    integer :: jefj
    integer :: jfj
    integer :: jhfj
    integer :: nconf_n
    real(dp), dimension(:, :), allocatable :: obs !(mobs,MSTATES)
    real(dp), dimension(:, :), allocatable :: s_diag !(mparm,MSTATES)
    real(dp), dimension(:), allocatable :: s_ii_inv !(mparm)
    real(dp), dimension(:, :), allocatable :: sr_ho !(mparm,mconf)
    real(dp), dimension(:, :), allocatable :: sr_o !(mparm,mconf)
    real(dp), dimension(:, :), allocatable :: wtg !(mconf,MSTATES)
    real(dp), dimension(:, :), allocatable :: obs_tot !(mobs,MSTATES)

    private
    public :: elocal, h_sr, jefj, jfj, jhfj, nconf_n, obs, s_diag, s_ii_inv, sr_ho, sr_o, wtg, obs_tot
    public :: allocate_sr_mat_n, deallocate_sr_mat_n
    save
contains
    subroutine allocate_sr_mat_n()
        use sr_mod, only: mparm, mobs, mconf
        use mstates_mod, only: MSTATES
        if (.not. allocated(elocal)) allocate (elocal(mconf, MSTATES), source=0.0_dp)
        if (.not. allocated(h_sr)) allocate (h_sr(mparm), source=0.0_dp)
        if (.not. allocated(obs)) allocate (obs(mobs, MSTATES), source=0.0_dp)
        if (.not. allocated(s_diag)) allocate (s_diag(mparm, MSTATES), source=0.0_dp)
        if (.not. allocated(s_ii_inv)) allocate (s_ii_inv(mparm), source=0.0_dp)
        if (.not. allocated(sr_ho)) allocate (sr_ho(mparm, mconf), source=0.0_dp)
        if (.not. allocated(sr_o)) allocate (sr_o(mparm, mconf), source=0.0_dp)
        if (.not. allocated(wtg)) allocate (wtg(mconf, MSTATES), source=0.0_dp)
        if (.not. allocated(obs_tot)) allocate (obs_tot(mobs, MSTATES), source=0.0_dp)
    end subroutine allocate_sr_mat_n

    subroutine deallocate_sr_mat_n()
        if (allocated(obs_tot)) deallocate (obs_tot)
        if (allocated(wtg)) deallocate (wtg)
        if (allocated(sr_o)) deallocate (sr_o)
        if (allocated(sr_ho)) deallocate (sr_ho)
        if (allocated(s_ii_inv)) deallocate (s_ii_inv)
        if (allocated(s_diag)) deallocate (s_diag)
        if (allocated(obs)) deallocate (obs)
        if (allocated(h_sr)) deallocate (h_sr)
        if (allocated(elocal)) deallocate (elocal)
    end subroutine deallocate_sr_mat_n

end module sr_mat_n

module m_sr
contains
subroutine allocate_m_sr()
    use sr_mat_n, only: allocate_sr_mat_n

    call allocate_sr_mat_n()
end subroutine allocate_m_sr


subroutine deallocate_m_sr()
    use sr_mat_n, only: deallocate_sr_mat_n

    call deallocate_sr_mat_n()
end subroutine deallocate_m_sr
end module 
