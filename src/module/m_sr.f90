module sr_mod
    !> Arguments:
    integer :: mparm
    integer :: mobs
    integer :: mconf
    integer :: i_sr_rescale, izvzb

    private
    public :: mparm, mobs, mconf
    public :: izvzb, i_sr_rescale
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
      use mstates_mod, only: MSTATES
      use precision_kinds, only: dp
      use sr_mod, only: mparm, mobs, mconf

    real(dp), dimension(:, :), allocatable :: elocal !(mconf,MSTATES)
    real(dp), dimension(:, :), allocatable :: h_sr !(mparm,MSTATES)
    real(dp), dimension(:, :), allocatable :: h_sr_penalty !(mparm,MSTATES)
    real(dp), dimension(:), allocatable :: isr_lambda !(MSTATES*(MSTATES-1)/2)
    real(dp), dimension(:, :), allocatable :: sr_lambda !(MSTATES,MSTATES)
    integer :: sr_state
    integer :: jefj
    integer :: jfj
    integer :: jhfj
    integer :: nconf_n
    real(dp), dimension(:, :), allocatable :: s_diag !(mparm,MSTATES)
    real(dp), dimension(:, :), allocatable :: s_ii_inv !(mparm,MSTATES)
    real(dp), dimension(:, :), allocatable :: sr_ho !(mparm,mconf)
    real(dp), dimension(:, :, :), allocatable :: sr_o !(mparm,mconf,MSTATES)
    real(dp), dimension(:, :), allocatable :: wtg !(mconf,MSTATES)
    real(dp), dimension(:, :), allocatable :: obs_tot !(mobs,MSTATES)

    private
    public :: elocal, h_sr, jefj, jfj, jhfj, nconf_n, s_diag, s_ii_inv, sr_ho, sr_o, wtg, obs_tot, isr_lambda, sr_lambda, sr_state, h_sr_penalty
    ! public :: obs
    public :: allocate_sr_mat_n, deallocate_sr_mat_n
    save
contains
    subroutine allocate_sr_mat_n()
      use mstates_mod, only: MSTATES
      use optwf_func, only: ifunc_omega
      use sr_mod, only: mparm, mobs, mconf, izvzb, i_sr_rescale
        if (.not. allocated(elocal)) allocate (elocal(mconf, MSTATES))
        if (.not. allocated(h_sr)) allocate (h_sr(mparm, MSTATES))
        if (.not. allocated(h_sr_penalty)) allocate (h_sr_penalty(mparm,MSTATES))
        if (.not. allocated(isr_lambda)) allocate (isr_lambda(MSTATES*(MSTATES-1)/2))
        if (.not. allocated(sr_lambda)) allocate (sr_lambda(MSTATES,MSTATES))
        if (.not. allocated(s_diag)) allocate (s_diag(mparm, MSTATES))
        if (.not. allocated(s_ii_inv)) allocate (s_ii_inv(mparm, MSTATES))
        if (.not. allocated(sr_ho)) allocate (sr_ho(mparm, mconf))
        if (.not. allocated(sr_o)) allocate (sr_o(mparm, mconf, MSTATES))
        if (.not. allocated(wtg)) allocate (wtg(mconf, MSTATES))
        if (.not. allocated(obs_tot)) allocate (obs_tot(mobs, MSTATES))
    end subroutine allocate_sr_mat_n

    subroutine deallocate_sr_mat_n()
        if (allocated(obs_tot)) deallocate (obs_tot)
        if (allocated(wtg)) deallocate (wtg)
        if (allocated(sr_o)) deallocate (sr_o)
        if (allocated(sr_ho)) deallocate (sr_ho)
        if (allocated(s_ii_inv)) deallocate (s_ii_inv)
        if (allocated(s_diag)) deallocate (s_diag)
        if (allocated(h_sr)) deallocate (h_sr)
        if (allocated(h_sr_penalty)) deallocate (h_sr_penalty)
        if (allocated(isr_lambda)) deallocate (isr_lambda)
        if (allocated(sr_lambda)) deallocate (sr_lambda)
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
