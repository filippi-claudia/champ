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
    use m_backflow, only: ibackflow, quasi_x, dquasi_dx, d2quasi_dx2, dquasi_dp
    use system, only: nelec
    use qmckl_data, only: use_qmckl_orbitals
    implicit none
    real(dp), dimension(3, nelec), intent(in) :: x

    if (ibackflow == 0) return

    if (ibackflow == 1) then
        call trivial_backflow(x, quasi_x, dquasi_dx, d2quasi_dx2)
    else if (ibackflow == 2) then
        call gaussian_backflow(x, quasi_x, dquasi_dx, d2quasi_dx2, dquasi_dp)
    else if (ibackflow == 3) then
        call linear_backflow(x, quasi_x, dquasi_dx, d2quasi_dx2, dquasi_dp)
    else if (ibackflow == 4) then
        call ee_backflow(x, quasi_x, dquasi_dx, d2quasi_dx2, dquasi_dp)
    else
        call fatal_error('Backflow type not recognized.')
    end if

    ! Only compute if not using QMCkl
    if (.not. use_qmckl_orbitals) then
        call bf_distances(quasi_x)
    end if

end

subroutine init_backflow()
    use m_backflow, only: ibackflow
    implicit none

    if (ibackflow == 0) return

    if (ibackflow == 2) then
        call init_gaussian_backflow()
    else if (ibackflow == 3) then
        continue
        !call init_linear_backflow()
    else if (ibackflow == 4) then
        call init_twobody_backflow(8,8)
    end if


end subroutine init_backflow

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

subroutine init_gaussian_backflow()
    use precision_kinds, only: dp
    use m_backflow, only: parm_bf, nparm_bf
    use m_backflow, only: allocate_m_backflow
    implicit none

    nparm_bf = 3
    call allocate_m_backflow

    parm_bf(1) = 0.0d0  ! lambda
    parm_bf(2) = 0.0d0  ! mean
    parm_bf(3) = 1.0d0  ! sigma

end subroutine init_gaussian_backflow

subroutine gaussian_backflow(x, quasi_x, dquasi_dx, d2quasi_dx2, dquasi_dp)
    use precision_kinds, only: dp
    use system, only: nelec
    use optwf_control, only: ioptci, ioptjas, ioptorb, ioptbf
    use m_backflow, only: parm_bf, nparm_bf
    implicit none
    real(dp), dimension(3, nelec), intent(in) :: x
    real(dp), dimension(3, nelec), intent(out) :: quasi_x
    real(dp), dimension(3, nelec, 3, nelec), intent(out) :: dquasi_dx
    real(dp), dimension(3, nelec, nelec), intent(out) :: d2quasi_dx2
    real(dp), dimension(3, nelec, nparm_bf), intent(out) :: dquasi_dp
    real(dp) :: sigma, mean, rij, G, dGdr, d2Gdr2, lambda
    real(dp) :: delta(3)
    integer :: i, j, k, a, b

    sigma=parm_bf(3)
    mean=parm_bf(2)
    lambda = parm_bf(1)


    quasi_x = 0.0_dp
    dquasi_dx = 0.0_dp
    d2quasi_dx2 = 0.0_dp
    dquasi_dp = 0.0_dp

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
            do a = 1, 3
                quasi_x(a, i) = quasi_x(a, i) + G * delta(a)
                do b = 1, 3
                    dquasi_dx(a,i,b,i) = dquasi_dx(a,i,b,i) + dGdr * delta(b)/rij * delta(a) 
                    dquasi_dx(a,i,b,j) = dquasi_dx(a,i,b,j) - dGdr * delta(b)/rij * delta(a) 
                enddo
                dquasi_dx(a,i,a,i) = dquasi_dx(a,i,a,i) + G
                dquasi_dx(a,i,a,j) = dquasi_dx(a,i,a,j) - G
                dquasi_dp(a,i,1) = dquasi_dp(a,i,1) + exp(-((rij-mean)/sigma)**2) * delta(a)
                dquasi_dp(a,i,2) = dquasi_dp(a,i,2) + G * 2.0d0 * (rij-mean) / sigma**2 * delta(a)
                dquasi_dp(a,i,3) = dquasi_dp(a,i,3) + G * 2.0d0 * (rij-mean)**2 / sigma**3 * delta(a)
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



subroutine linear_backflow(x, quasi_x, dquasi_dx, d2quasi_dx2, dquasi_dp)
    use precision_kinds, only: dp
    use system, only: nelec
    use m_backflow, only: parm_bf, nparm_bf
    implicit none
    real(dp), dimension(3, nelec), intent(in) :: x
    real(dp), dimension(3, nelec), intent(out) :: quasi_x
    real(dp), dimension(3, nelec, 3, nelec), intent(out) :: dquasi_dx
    real(dp), dimension(3, nelec, nelec), intent(out) :: d2quasi_dx2
    real(dp), dimension(3, nelec, 3), intent(out) :: dquasi_dp
    real(dp) :: sigma, mean, rij, G, dGdr, d2Gdr2, lambda
    real(dp) :: delta(3)
    integer :: i, j, k, a, b

    lambda = parm_bf(1)


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
                dquasi_dp(a,i,1) = dquasi_dp(a,i,1) + delta(a)
                    
            end do

        end do
    end do


end subroutine linear_backflow

subroutine init_twobody_backflow(orda, ordb)
    use precision_kinds, only: dp
    use m_backflow, only: parm_bf, nparm_bf, norda_bf, nordb_bf
    use system, only: nctype
    use m_backflow, only: allocate_m_backflow
    implicit none
    integer :: i, orda, ordb
 
    norda_bf = orda
    nordb_bf = ordb

    nparm_bf = 0
    if (norda_bf .gt. 0) then
        nparm_bf = nparm_bf + (2 + norda_bf) * nctype
    end if
    if (nordb_bf .gt. 0) then
        nparm_bf = nparm_bf + (2 + nordb_bf)
    end if

    call allocate_m_backflow
    
    parm_bf = 0.0d0
    parm_bf(1) = 3.0d0

    if (nordb_bf.gt.0) then
        do i = 1, nctype
            parm_bf(2 + nordb_bf + (2 + norda_bf)*(i-1) + 1) = 3.0d0
        end do
    end if
end subroutine init_twobody_backflow

subroutine ee_backflow(x, quasi_x, dquasi_dx, d2quasi_dx2, dquasi_dp)
    use precision_kinds, only: dp
    use system, only: nelec, iwctype, ncent, nctype, cent
    use optwf_control, only: ioptci, ioptjas, ioptorb, ioptbf
    use m_backflow, only: parm_bf, nparm_bf, norda_bf, nordb_bf
    implicit none
    real(dp), dimension(3, nelec), intent(in) :: x
    real(dp), dimension(3, nelec), intent(out) :: quasi_x
    real(dp), dimension(3, nelec, 3, nelec), intent(out) :: dquasi_dx
    real(dp), dimension(3, nelec, nelec), intent(out) :: d2quasi_dx2
    real(dp), dimension(3, nelec, nparm_bf), intent(out) :: dquasi_dp
    real(dp) :: rij, rr, C, cutoff
    real(dp) :: f, fp(3), fpp(3,3), eta, etap(3), etapp(3,3)
    real(dp) :: delta(3)
    integer :: i, j, k, a, b, offset_ee, offset_en, offset_een

    offset_ee = 0
    offset_en = offset_ee + 1+1+nordb_bf
    offset_een = offset_en + (1+1+norda_bf)*nctype

    quasi_x = 0.0_dp
    dquasi_dx = 0.0_dp
    d2quasi_dx2 = 0.0_dp
    dquasi_dp = 0.0_dp

    
    C = 3

    do i = 1, nelec
        do a=1,3
            quasi_x(a,i) = x(a,i)
            dquasi_dx(a,i,a,i) = 1.0d0
        enddo
    end do

    if (nordb_bf .eq. 0) goto 10
    cutoff = parm_bf(offset_ee+1)

    !parm_bf(2) = 0.0d0
    ! parm_bf(3) = 0.0d0

    do i = 1, nelec
        do j = 1, nelec
            if (i == j) cycle
            eta=0.0d0
            etap=0.0d0
            etapp=0.0d0

            ! delta_ij = r_i - r_j
            delta(:) = x(:, i) - x(:, j)

            rij = sqrt(delta(1)**2 + delta(2)**2 + delta(3)**2)
            if (rij > cutoff) cycle

            rr = 1.0d0

            f = ((cutoff - rij)/cutoff)**C


            do a = 1, 3
                fp(a) = -C/cutoff * ((cutoff - rij)/cutoff)**(C-1) * (delta(a)/rij)
                do b = 1, 3
                    fpp(a,b) = & 
                        (-C/cutoff) * ( - (C-1)/cutoff ) * ((cutoff - rij)/cutoff)**(C-2) * (delta(a)/rij) * (delta(b)/rij) + &
                        (-C/cutoff) * ((cutoff - rij)/cutoff)**(C-1) * ( (1/rij) - (delta(a)*delta(b)/(rij*rij*rij)) )
                end do
            end do

            do k = 0, nordb_bf
                eta = eta + parm_bf(offset_ee+k+2)*rr
                do a = 1, 3
                    dquasi_dp(a,i,offset_ee+k+2) = dquasi_dp(a,i,offset_ee+k+2) + rr * delta(a) * f
                    etap(a) = etap(a) + parm_bf(offset_ee+k+2) * (&
                        k * rr * delta(a)/(rij*rij) )
                    do b = 1, 3
                        etapp(a,b) = etapp(a,b) + parm_bf(offset_ee+k+2) * (&
                            k * (k-2) * rr * delta(a) * delta(b)/(rij*rij*rij*rij) )
                    enddo
                    etapp(a,a) = etapp(a,a) + parm_bf(offset_ee+k+2) * (&
                        k * rr/(rij*rij) )
                enddo
                rr = rr * rij
            end do

            do a = 1, 3
                quasi_x(a,i) = quasi_x(a,i) + eta * delta(a) * f
                ! dquasi_dp(a,i,1) = dquasi_dp(a,i,1) + eta * delta(a) * f * log((cutoff - rij)/cutoff)
                dquasi_dp(a,i,offset_ee+1) = dquasi_dp(a,i,offset_ee+1) + eta * delta(a) * C * ((cutoff - rij)/cutoff)**(C-1) * (rij/cutoff/cutoff)
                do b = 1, 3
                    dquasi_dx(a,i,b,i) = dquasi_dx(a,i,b,i) + (&
                        etap(b) * delta(a) * f + eta * delta(a) * fp(b) )
                    dquasi_dx(a,i,b,j) = dquasi_dx(a,i,b,j) - (&
                        etap(b) * delta(a) * f + eta * delta(a) * fp(b) )

                    d2quasi_dx2(a,i,i) = d2quasi_dx2(a,i,i) + (&
                        etapp(b,b) * delta(a) * f + &
                        2.0d0 * etap(b) * delta(a) * fp(b) + &
                        eta * delta(a) * fpp(b,b) )

                    d2quasi_dx2(a,i,j) = d2quasi_dx2(a,i,j) + (&
                        etapp(b,b) * delta(a) * f + &
                        2.0d0 * etap(b) * delta(a) * fp(b) + &
                        eta * delta(a) * fpp(b,b) )
                end do
                dquasi_dx(a,i,a,i) = dquasi_dx(a,i,a,i) + eta * f
                dquasi_dx(a,i,a,j) = dquasi_dx(a,i,a,j) - eta * f

                d2quasi_dx2(a,i,i) = d2quasi_dx2(a,i,i) + 2 * (eta * fp(a) + etap(a) * f)
                d2quasi_dx2(a,i,j) = d2quasi_dx2(a,i,j) + 2 * (eta * fp(a) + etap(a) * f)

            end do

        
        end do

    end do


10  continue

if (norda_bf .eq. 0) return

    ! do k = 1, nctype
    !    parm_bf(offset_en+(k-1)*(norda_bf+2)+2) = 0.0d0
    !    parm_bf(offset_en+(k-1)*(norda_bf+2)+3) = 0.0d0
    ! end do



    do j = 1, ncent
        cutoff = parm_bf(offset_en+(iwctype(j)-1)*(norda_bf+2)+1)

        do i = 1, nelec

             ! delta_ij = r_i - r_j
            delta(:) = x(:, i) - cent(:, j)

            rij = sqrt(delta(1)**2 + delta(2)**2 + delta(3)**2)
            if (rij > cutoff) cycle

            eta=0.0d0
            etap=0.0d0
            etapp=0.0d0

            rr = 1.0d0

            f = ((cutoff - rij)/cutoff)**C

            do a = 1, 3
                fp(a) = -C/cutoff * ((cutoff - rij)/cutoff)**(C-1) * (delta(a)/rij)
                do b = 1, 3
                    fpp(a,b) = & 
                        (-C/cutoff) * ( - (C-1)/cutoff ) * ((cutoff - rij)/cutoff)**(C-2) * (delta(a)/rij) * (delta(b)/rij) + &
                        (-C/cutoff) * ((cutoff - rij)/cutoff)**(C-1) * ( (1/rij) - (delta(a)*delta(b)/(rij*rij*rij)) )
                end do
            end do

            do k = 0, norda_bf
                eta = eta + parm_bf(offset_en + (iwctype(j)-1)*(norda_bf+2)+ k+2)*rr
                do a = 1, 3
                    dquasi_dp(a,i,offset_en + (iwctype(j)-1)*(norda_bf+2)+ k+2) = &
                         dquasi_dp(a,i,offset_en + (iwctype(j)-1)*(norda_bf+2)+ k+2) + rr * delta(a) * f
                    etap(a) = etap(a) + parm_bf(offset_en + (iwctype(j)-1)*(norda_bf+2)+ k+2) * (&
                        k * rr * delta(a)/(rij*rij) )
                    do b = 1, 3
                        etapp(a,b) = etapp(a,b) + parm_bf(offset_en + (iwctype(j)-1)*(norda_bf+2)+ k+2) * (&
                            k * (k-2) * rr * delta(a) * delta(b)/(rij*rij*rij*rij) )
                    enddo
                    etapp(a,a) = etapp(a,a) + parm_bf(offset_en + (iwctype(j)-1)*(norda_bf+2)+ k+2) * (&
                        k * rr/(rij*rij) )
                enddo
                rr = rr * rij
            end do

            do a = 1, 3
                quasi_x(a,i) = quasi_x(a,i) + eta * delta(a) * f
                dquasi_dp(a,i,offset_en + (iwctype(j)-1)*(norda_bf+2)+1) = dquasi_dp(a,i,offset_en + (iwctype(j)-1)*(norda_bf+2)+1) &
                    + eta * delta(a) * C * ((cutoff - rij)/cutoff)**(C-1) * (rij/cutoff/cutoff)
                do b = 1, 3
                    dquasi_dx(a,i,b,i) = dquasi_dx(a,i,b,i) + (&
                        etap(b) * delta(a) * f + eta * delta(a) * fp(b) )

                    d2quasi_dx2(a,i,i) = d2quasi_dx2(a,i,i) + (&
                        etapp(b,b) * delta(a) * f + &
                        2.0d0 * etap(b) * delta(a) * fp(b) + &
                        eta * delta(a) * fpp(b,b) )

                end do
                dquasi_dx(a,i,a,i) = dquasi_dx(a,i,a,i) + eta * f

                d2quasi_dx2(a,i,i) = d2quasi_dx2(a,i,i) + 2 * (eta * fp(a) + etap(a) * f)

            end do
            
        end do
    end do


end subroutine ee_backflow


end module