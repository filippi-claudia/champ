module backflow_mod

contains

!> This function selects the backflow transformation to be used
!>
!> @details Backflow functions can be selected in the input file using the flag "backflow".
!> Flag 0: Blackflow disabled
!> Flag 1: Trivial backflow (no transformation)
!> @param[in]  x  Original electron coordinates (3, nelec)
!> @param[out] none
!>
!> @author Emiel Slootman
!> @date Nobember 2025
subroutine backflow(x)
    use precision_kinds, only: dp
    use error, only: fatal_error
    use m_backflow, only: ibackflow, quasi_x, dquasi_dx, d2quasi_dx2
    use system, only: nelec
    use qmckl_data, only: use_qmckl_orbitals
    implicit none
    real(dp), dimension(3, nelec), intent(in) :: x

    if (ibackflow == 0) return

    if (ibackflow == 1) then
        call trivial_backflow(x, quasi_x, dquasi_dx, d2quasi_dx2)
    else
        call fatal_error('Backflow type not recognized.')
    end if

    ! Only compute if not using QMCkl
    if (.not. use_qmckl_orbitals) then
        call bf_distances(quasi_x)
    end if

end

!> This subroutine computes distances between quasicoordinates and nuclei
!>
!> @param[in]  quasi_x  Quasicoordinates of electrons (3, nelec)
!>
!> @author Emiel Slootman
!> @date November 2025
subroutine bf_distances(quasi_x)
    use precision_kinds, only: dp
    use system, only: nelec, ncent_tot, cent
    use m_backflow, only: rvec_en_bf, r_en_bf

    implicit none
    real(dp), dimension(3, nelec), intent(in) :: quasi_x
    real(dp), dimension(3) :: rij
    integer :: iel, icent

    do iel = 1, nelec
        do icent = 1, ncent_tot
            rij = quasi_x(:, iel) - cent(:, icent)
            rvec_en_bf(:, iel, icent) = rij
            r_en_bf(iel, icent) = sqrt(dot_product(rij, rij))
        end do
    end do

end


subroutine trivial_backflow(x, quasi_x, dquasi_dx, d2quasi_dx2)
    use precision_kinds, only: dp
    use system, only: nelec
    implicit none
    real(dp), dimension(3, nelec), intent(in) :: x
    real(dp), dimension(3, nelec), intent(out) :: quasi_x
    real(dp), dimension(3, nelec, 3, nelec), intent(out) :: dquasi_dx
    real(dp), dimension(3, nelec, 3, 3, nelec), intent(out) :: d2quasi_dx2
    integer :: i

    ! Quasicoordinates equal to original coordinates
    quasi_x = x

    ! First derivatives are identity matrices
    dquasi_dx = 0.0d0
    do i = 1, nelec
        dquasi_dx(:, i, :, i) = 0.0d0
        dquasi_dx(1, i, 1, i) = 1.0d0
        dquasi_dx(2, i, 2, i) = 1.0d0
        dquasi_dx(3, i, 3, i) = 1.0d0
    end do

    ! Second derivatives are zero
    d2quasi_dx2 = 0.0d0
end


end module