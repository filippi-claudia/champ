module backflow_mod

contains

!> This function selects the backflow transformation to be used
!>
!> @details Backflow functions can be selected in the input file using the flag "backflow".
!> Flag 0: Blackflow disabled
!> Flag 1: Trivial backflow (no transformation)
!>
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
    else if (ibackflow == 2) then
        call gaussian_backflow(x, quasi_x, dquasi_dx, d2quasi_dx2)
    else if (ibackflow == 3) then
        call linear_backflow(x, quasi_x, dquasi_dx, d2quasi_dx2)
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
    integer :: iel, icent, k

    do iel = 1, nelec
        do icent = 1, ncent_tot
            do k=1,3
                rij(k) = quasi_x(k, iel) - cent(k, icent)
                rvec_en_bf(k, iel, icent) = rij(k)
            enddo
            r_en_bf(iel, icent) = sqrt(dot_product(rij(:), rij(:)))
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
    real(dp), dimension(3, nelec, nelec), intent(out) :: d2quasi_dx2
    integer :: i

    ! Quasicoordinates equal to original coordinates
    quasi_x = x

    ! First derivatives are identity matrices
    dquasi_dx = 0.0d0
    do i = 1, nelec
        dquasi_dx(1, i, 1, i) = 1.0d0
        dquasi_dx(2, i, 2, i) = 1.0d0
        dquasi_dx(3, i, 3, i) = 1.0d0
    end do

    ! Second derivatives are zero
    d2quasi_dx2 = 0.0d0
end

subroutine gaussian_backflow(x, quasi_x, dquasi_dx, d2quasi_dx2)
    use precision_kinds, only: dp
    use system, only: nelec
    implicit none
    real(dp), dimension(3, nelec), intent(in) :: x
    real(dp), dimension(3, nelec), intent(out) :: quasi_x
    real(dp), dimension(3, nelec, 3, nelec), intent(out) :: dquasi_dx
    real(dp), dimension(3, nelec, nelec), intent(out) :: d2quasi_dx2
    real(dp) :: sigma, mean, rij, G, dGdr, d2Gdr2, lambda
    real(dp) :: delta(3)
    integer :: i, j, k, a, b

    sigma=0.7d0
    mean=0.0d0
    lambda = 0.1d0


    quasi_x = 0.0_dp
    dquasi_dx = 0.0_dp
    d2quasi_dx2 = 0.0_dp

    ! Main loop
    do i = 1, nelec


        do a=1,3
            quasi_x(a,i) = x(a,i)
            dquasi_dx(a,i,a,i) = 1.0d0
        enddo

        do j = 1, nelec
            if (i == j) cycle

            ! delta_ij = r_i - r_j
            delta(:) = x(:, i) - x(:, j)

            rij = sqrt(delta(1)**2 + delta(2)**2 + delta(3)**2)

            ! Gaussian
            G = lambda * exp(-((rij-mean)/sigma)**2)

            ! dG/dr
            dGdr = G * ( - 2.0d0 * (rij-mean) / sigma**2 )

            ! d2G/dr2
            d2Gdr2 = G * ( (4.0d0*(rij-mean)**2)/sigma**4 - 2.0d0/sigma**2 - 4.0d0*(rij-mean)/sigma**2/rij)

            ! Update quasi_x


            !do k = 1, nelec
            !    do a = 1, 3
            !        do b = 1, 3

            !            dquasi_dx(a,i,b,k) = dquasi_dx(a,i,b,k) + &
            !                dGdr * delta(b)/rij * delta(a) * merge(1.0_dp, 0.0_dp, i==k) - &
            !                dGdr * delta(b)/rij * delta(a) * merge(1.0_dp, 0.0_dp, j==k)
            !            dquasi_dx(a,i,b,k) = dquasi_dx(a,i,b,k) + &
            !                G * merge(1.0_dp, 0.0_dp, a==b) * merge(1.0_dp, 0.0_dp, i==k) - &
            !                G * merge(1.0_dp, 0.0_dp, a==b) * merge(1.0_dp, 0.0_dp, j==k)                        

            !        end do
            !    end do
            !end do

            do a = 1, 3
                quasi_x(a, i) = quasi_x(a, i) + G * delta(a)
                do b = 1, 3
                    dquasi_dx(a,i,b,i) = dquasi_dx(a,i,b,i) + dGdr * delta(b)/rij * delta(a) 
                    dquasi_dx(a,i,b,j) = dquasi_dx(a,i,b,j) - dGdr * delta(b)/rij * delta(a) 
                enddo
                dquasi_dx(a,i,a,i) = dquasi_dx(a,i,a,i) + G
                dquasi_dx(a,i,a,j) = dquasi_dx(a,i,a,j) - G
            enddo


            do a = 1, 3
                d2quasi_dx2(a,i,i) = d2quasi_dx2(a,i,i) &
                    + delta(a) * d2Gdr2 &
                    + 2.0d0 * dGdr * delta(a)/rij
                d2quasi_dx2(a,i,j) = d2quasi_dx2(a,i,j) &
                    + delta(a) * d2Gdr2 &
                    + 2.0d0 * dGdr * delta(a)/rij
            end do

        end do
    end do


end subroutine gaussian_backflow



subroutine linear_backflow(x, quasi_x, dquasi_dx, d2quasi_dx2)
    use precision_kinds, only: dp
    use system, only: nelec
    implicit none
    real(dp), dimension(3, nelec), intent(in) :: x
    real(dp), dimension(3, nelec), intent(out) :: quasi_x
    real(dp), dimension(3, nelec, 3, nelec), intent(out) :: dquasi_dx
    real(dp), dimension(3, nelec, nelec), intent(out) :: d2quasi_dx2
    real(dp) :: sigma, mean, rij, G, dGdr, d2Gdr2, lambda
    real(dp) :: delta(3)
    integer :: i, j, k, a, b

    sigma=0.7d0
    mean=0.0d0
    lambda = 0.1d0


    quasi_x = 0.0_dp
    dquasi_dx = 0.0_dp
    d2quasi_dx2 = 0.0_dp

    ! Main loop
    do i = 1, nelec

        ! Start with q_i = r_i
        
        do a=1,3
            quasi_x(a, i) = x(a, i)

            dquasi_dx(a,i,a,i) = 1.0d0
        enddo

        do j = 1, nelec
            if (i == j) cycle

            do a = 1, 3
                delta(a) = x(a, i) - x(a, j)

                quasi_x(a, i) = quasi_x(a, i) + lambda * delta(a)
                dquasi_dx(a,i,a,i) = dquasi_dx(a,i,a,i) + lambda
                dquasi_dx(a,i,a,j) = dquasi_dx(a,i,a,j) - lambda
                    
            end do

        end do
    end do


end subroutine linear_backflow


end module