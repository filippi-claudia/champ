module optwf_func
    !> Arguments: ifunc_omega, omega, omega0, omega_hes, n_omegaf, n_omegat
      use precision_kinds, only: dp

    integer :: ifunc_omega
    real(dp) :: omega
    real(dp) :: omega0
    real(dp) :: omega_hes
    integer :: n_omegaf
    integer :: n_omegat

    private
    public :: ifunc_omega, omega, omega0, omega_hes, n_omegaf, n_omegat
    save
end module optwf_func

