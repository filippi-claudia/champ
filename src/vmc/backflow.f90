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
        call rios_backflow(x, quasi_x, dquasi_dx, d2quasi_dx2, dquasi_dp)
    else
        call fatal_error('Backflow type not recognized.')
    end if

    ! Only compute if not using QMCkl
    if (.not. use_qmckl_orbitals) then
        call bf_distances(quasi_x)
    end if

end

subroutine single_backflow(iel, xold, xnew, quasi_x_new, dquasi_dx_new, d2quasi_dx2_new, indices)
    use precision_kinds, only: dp
    use m_backflow, only: ibackflow
    use system, only: nelec
    use qmckl_data, only: use_qmckl_orbitals
    use error, only: fatal_error
    implicit none
    integer, intent(in) :: iel
    real(dp), dimension(3, nelec), intent(in) :: xold
    real(dp), dimension(3), intent(in) :: xnew
    real(dp), dimension(3, nelec), intent(out) :: quasi_x_new
    real(dp), dimension(3, nelec, 3, nelec), intent(out) :: dquasi_dx_new
    real(dp), dimension(3, nelec, nelec), intent(out) :: d2quasi_dx2_new
    integer, dimension(nelec), intent(out) :: indices
    
    if (ibackflow == 0) return

    if (ibackflow == 4) then
        call single_rios_backflow(iel, xold, xnew, quasi_x_new, dquasi_dx_new, d2quasi_dx2_new, indices)
    else
        call fatal_error('Single backflow only implemented for Rios backflow.')
    end if

    ! Only compute if not using QMCkl
    if (.not. use_qmckl_orbitals) then
        call bf_distances(quasi_x_new)
    end if
end subroutine single_backflow

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
        call init_rios_backflow(5,5,5)
    end if
end subroutine init_backflow

subroutine init_backflow_arrays()
    use m_backflow, only: ibackflow
    implicit none

    if (ibackflow == 0) return

    if (ibackflow == 4) then
        call init_rios_backflow_arrays()
    end if
end subroutine init_backflow_arrays

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

subroutine init_rios_backflow(orda, ordb, ordc)
    use precision_kinds, only: dp
    use m_backflow, only: parm_bf, nparm_bf, norda_bf, nordb_bf, nordc_bf, cutoff_scale
    use system, only: nctype
    use m_backflow, only: allocate_m_backflow
    implicit none
    integer :: i, orda, ordb, ordc, tmpc, l, m, n
    intrinsic :: ceiling
 
    norda_bf = orda
    nordb_bf = ordb
    nordc_bf = ordc

    call init_backflow_arrays()
    
    parm_bf = 0.0d0
    parm_bf(1) = 3.0d0

    if (norda_bf.gt.0) then
        do i = 1, nctype
            parm_bf(1 + nordb_bf + (1 + norda_bf)*(i-1) + 1) = 3.0d0
        end do
    end if
    if (nordc_bf.gt.0) then
        tmpc = 0
        do l = 0, nordc_bf
            do m = 0, nordc_bf - l
                do n = 0, nordc_bf - l - m
                    tmpc = tmpc + 1
                end do
            end do
        end do

        do i = 1, 2 * nctype
            parm_bf(1 + nordb_bf + (1 + norda_bf)*nctype + (tmpc+1) * (i-1) + 1) = 3.0d0
        end do
    end if
    cutoff_scale = 3
end subroutine init_rios_backflow

subroutine init_rios_backflow_arrays()
    use m_backflow, only: allocate_m_backflow, ibackflow, norda_bf, nordb_bf, nordc_bf, nparm_bf, maxord
    use system, only: nctype
    implicit none
    integer :: i, tmpc, l, m, n
    intrinsic :: ceiling

    nparm_bf = 0
    if (norda_bf .gt. 0) then
        nparm_bf = nparm_bf + (1 + norda_bf) * nctype
    end if
    if (nordb_bf .gt. 0) then
        nparm_bf = nparm_bf + (1 + nordb_bf)
    end if
    if (nordc_bf .gt. 0) then
        tmpc = 0
        do l = 0, nordc_bf
            do m = 0, nordc_bf - l
                do n = 0, nordc_bf - l - m
                    tmpc = tmpc + 1
                end do
            end do
        end do

        nparm_bf = nparm_bf + (tmpc + 1) * nctype + (tmpc) * nctype
    end if

    maxord = max(norda_bf, nordb_bf, nordc_bf)

    call allocate_m_backflow

end subroutine init_rios_backflow_arrays


subroutine rios_distances(x)
    use precision_kinds, only: dp
    use system, only: nelec, ncent, cent, iwctype, nctype
    use m_backflow, only: parm_bf, nparm_bf, norda_bf, nordb_bf, nordc_bf, cutoff_scale
    use m_backflow, only: r_en, rvec_en, r_ee, rvec_ee, r_ee_gl, r_en_gl, maxord, cutoff_deriv
    implicit none
    real(dp), dimension(3, nelec), intent(in) :: x
    real(dp) :: r, inv_r, r_cutoff, inv_r_cutoff, cutoff, inv_cutoff
    integer :: i, j, nc, k, no, multb, multa, tmpc, l, m ,n, cc

    multb = 0
    if (nordb_bf .gt. 0) multb = 1
    multa = 0
    if (norda_bf .gt. 0) multa = 1

    tmpc = 0
    do l = 0, nordc_bf
        do m = 0, nordc_bf - l
            do n = 0, nordc_bf - l - m
                tmpc = tmpc + 1
            end do
        end do
    end do
    
    r_en = 0.0d0
    r_ee = 0.0d0
    rvec_en = 0.0d0
    rvec_ee = 0.0d0
    r_ee_gl = 0.0d0
    r_en_gl = 0.0d0

    do cc = 1, 2
        do nc = 1, ncent
            if (cc == 1) then
                r_cutoff = parm_bf((1+nordb_bf) + (iwctype(nc)-1)*(1+norda_bf) + 1)
            else
                r_cutoff = parm_bf((1+nordb_bf) + (1+norda_bf)*nctype + (tmpc+1)*(iwctype(nc)-1) + 1)
            end if

            inv_r_cutoff = 1.0d0 / r_cutoff  

            do i = 1, nelec
                do k = 1, 3
                    rvec_en(k, i, nc) = x(k, i) - cent(k, nc)
                end do
                r = sqrt( rvec_en(1,i,nc)**2 + rvec_en(2,i,nc)**2 + rvec_en(3,i,nc)**2 )

                inv_r = 1.0d0 / r 
                
                cutoff= ((r_cutoff - r) * inv_r_cutoff)**cutoff_scale
                inv_cutoff = 1/((r_cutoff - r) * inv_r_cutoff)

                r_en(i, nc, 0, cc) = cutoff
                if (r > r_cutoff) cycle

                cutoff_deriv(i, nc) = cutoff_scale * r/(r_cutoff * (r_cutoff - r))


                r_en_gl(i, 1, nc, 0, cc) = -cutoff_scale * inv_r_cutoff * inv_r * rvec_en(1, i, nc) * cutoff * inv_cutoff
                r_en_gl(i, 2, nc, 0, cc) = -cutoff_scale * inv_r_cutoff * inv_r * rvec_en(2, i, nc) * cutoff * inv_cutoff
                r_en_gl(i, 3, nc, 0, cc) = -cutoff_scale * inv_r_cutoff * inv_r * rvec_en(3, i, nc) * cutoff * inv_cutoff
                r_en_gl(i, 4, nc, 0, cc) = -cutoff_scale * inv_r_cutoff * 2.0d0 * inv_r * cutoff * inv_cutoff + &
                                            cutoff_scale * (cutoff_scale - 1) * inv_r_cutoff * inv_r_cutoff * inv_cutoff* inv_cutoff * cutoff 

                do no = 1, maxord
                    r_en(i, nc, no, cc) = r_en(i, nc, no-1, cc) * r 
                    r_en_gl(i, 1, nc, no, cc) = no * r_en(i, nc, no-1, cc) * inv_r * rvec_en(1, i, nc) - &
                                                cutoff_scale * inv_r_cutoff * inv_r * rvec_en(1, i, nc) * inv_cutoff * r_en(i, nc, no, cc)
                    r_en_gl(i, 2, nc, no, cc) = no * r_en(i, nc, no-1, cc) * inv_r * rvec_en(2, i, nc) - &
                                                cutoff_scale * inv_r_cutoff * inv_r * rvec_en(2, i, nc) * inv_cutoff * r_en(i, nc, no, cc)
                    r_en_gl(i, 3, nc, no, cc) = no * r_en(i, nc, no-1, cc) * inv_r * rvec_en(3, i, nc) - &
                                                cutoff_scale * inv_r_cutoff * inv_r * rvec_en(3, i, nc) * inv_cutoff * r_en(i, nc, no, cc)
                    !r_en_gl(i, 4, nc, no, cc) = r_en(i, nc, no-1, cc) * no * 2.0d0 * (&
                    !                            inv_r - cutoff_scale * inv_r_cutoff * inv_cutoff) + &
                    !                            r_en(i, nc, no, cc) * cutoff_scale * inv_r_cutoff * inv_cutoff * (&
                    !                            (cutoff_scale - 1) * inv_r_cutoff * inv_cutoff - &
                    !                            2.0d0 * inv_r)
                    r_en_gl(i, 4, nc, no, cc) = - cutoff_scale * inv_r_cutoff * (no+1) * 2.0d0 * r_en(i, nc, no-1, cc) * inv_cutoff +&
                                                no * (no+1) * r_en(i, nc, no-1, cc) * inv_r + &
                                                cutoff_scale * (cutoff_scale - 1) * inv_r_cutoff * inv_r_cutoff * inv_cutoff * inv_cutoff * r_en(i, nc, no, cc) 
                end do
                !do no = 2, maxord
                !    r_en_gl(i, 4, nc, no, cc) = r_en_gl(i, 4, nc, no, cc) + no * (no - 1) * r_en(i, nc, no-2, cc)
                !end do
            end do
        end do
    end do

    do j = 1, nelec
        do i = 1, nelec
            if (i == j) cycle
            do k = 1, 3
                rvec_ee(k, i, j) = x(k, i) - x(k, j)
            end do
            r = sqrt( rvec_ee(1,i,j)**2 + rvec_ee(2,i,j)**2 + rvec_ee(3,i,j)**2 )
            inv_r = 1.0d0 / r

            r_ee(i, j, 0) = 1.0d0
            r_ee_gl(i, 1, j, 0) = 0.0d0
            r_ee_gl(i, 2, j, 0) = 0.0d0
            r_ee_gl(i, 3, j, 0) = 0.0d0
            r_ee_gl(i, 4, j, 0) = 0.0d0

            do no = 1, maxord
                r_ee(i, j, no) = r_ee(i, j, no-1) * r
                r_ee_gl(i, 1, j, no) = no * r_ee(i, j, no-1) * inv_r * rvec_ee(1, i, j)
                r_ee_gl(i, 2, j, no) = no * r_ee(i, j, no-1) * inv_r * rvec_ee(2, i, j)
                r_ee_gl(i, 3, j, no) = no * r_ee(i, j, no-1) * inv_r * rvec_ee(3, i, j)
                r_ee_gl(i, 4, j, no) = no * r_ee(i, j, no-1) * 2.0d0 * inv_r
            end do
            do no = 2, maxord
                r_ee_gl(i, 4, j, no) = r_ee_gl(i, 4, j, no) + no * (no - 1) * r_ee(i, j, no-2)
            end do
        end do
    end do
end subroutine rios_distances

subroutine single_rios_distances(x, xnew, iel)
    use precision_kinds, only: dp
    use system, only: nelec, ncent, cent, iwctype, nctype
    use m_backflow, only: parm_bf, nparm_bf, norda_bf, nordb_bf, nordc_bf, cutoff_scale
    use m_backflow, only: r_en, rvec_en, r_ee, rvec_ee, r_ee_gl, r_en_gl, maxord
    implicit none
    real(dp), dimension(3, nelec), intent(in) :: x
    real(dp), dimension(3), intent(in) :: xnew
    integer, intent(in) :: iel
    real(dp) :: r, inv_r, r_cutoff, inv_r_cutoff, cutoff, inv_cutoff
    integer :: i, j, nc, k, no, multb, multa, tmpc, l, m ,n, cc

    multb = 0
    if (nordb_bf .gt. 0) multb = 1
    multa = 0
    if (norda_bf .gt. 0) multa = 1

    tmpc = 0
    do l = 0, nordc_bf
        do m = 0, nordc_bf - l
            do n = 0, nordc_bf - l - m
                tmpc = tmpc + 1
            end do
        end do
    end do
    

    do cc = 1, 2
        do nc = 1, ncent
            if (cc == 1) then
                r_cutoff = parm_bf((1+nordb_bf) + (iwctype(nc)-1)*(1+norda_bf) + 1)
            else
                r_cutoff = parm_bf((1+nordb_bf) + (1+norda_bf)*nctype + (tmpc+1)*(iwctype(nc)-1) + 1)
            end if

            inv_r_cutoff = 1.0d0 / r_cutoff  

    
            do k = 1, 3
                rvec_en(k, iel, nc) = xnew(k) - cent(k, nc)
            end do
            r = sqrt( rvec_en(1,iel,nc)**2 + rvec_en(2,iel,nc)**2 + rvec_en(3,iel,nc)**2 )

            inv_r = 1.0d0 / r 
            
            cutoff= ((r_cutoff - r) * inv_r_cutoff)**cutoff_scale
            inv_cutoff = 1/((r_cutoff - r) * inv_r_cutoff)

            r_en(iel, nc, 0, cc) = cutoff
            if (r > r_cutoff) cycle


            r_en_gl(iel, 1, nc, 0, cc) = -cutoff_scale * inv_r_cutoff * inv_r * rvec_en(1, iel, nc) * cutoff * inv_cutoff
            r_en_gl(iel, 2, nc, 0, cc) = -cutoff_scale * inv_r_cutoff * inv_r * rvec_en(2, iel, nc) * cutoff * inv_cutoff
            r_en_gl(iel, 3, nc, 0, cc) = -cutoff_scale * inv_r_cutoff * inv_r * rvec_en(3, iel, nc) * cutoff * inv_cutoff
            r_en_gl(iel, 4, nc, 0, cc) = -cutoff_scale * inv_r_cutoff * 2.0d0 * inv_r * cutoff * inv_cutoff + &
                                        cutoff_scale * (cutoff_scale - 1) * inv_r_cutoff * inv_r_cutoff * inv_cutoff* inv_cutoff * cutoff 

            do no = 1, maxord
                r_en(iel, nc, no, cc) = r_en(iel, nc, no-1, cc) * r 
                r_en_gl(iel, 1, nc, no, cc) = no * r_en(iel, nc, no-1, cc) * inv_r * rvec_en(1, iel, nc) - &
                                            cutoff_scale * inv_r_cutoff * inv_r * rvec_en(1, iel, nc) * inv_cutoff * r_en(iel, nc, no, cc)
                r_en_gl(iel, 2, nc, no, cc) = no * r_en(iel, nc, no-1, cc) * inv_r * rvec_en(2, iel, nc) - &
                                            cutoff_scale * inv_r_cutoff * inv_r * rvec_en(2, iel, nc) * inv_cutoff * r_en(iel, nc, no, cc)
                r_en_gl(iel, 3, nc, no, cc) = no * r_en(iel, nc, no-1, cc) * inv_r * rvec_en(3, iel, nc) - &
                                            cutoff_scale * inv_r_cutoff * inv_r * rvec_en(3, iel, nc) * inv_cutoff * r_en(iel, nc, no, cc)
                r_en_gl(iel, 4, nc, no, cc) = - cutoff_scale * inv_r_cutoff * (no+1) * 2.0d0 * r_en(iel, nc, no-1, cc) * inv_cutoff +&
                                            no * (no+1) * r_en(iel, nc, no-1, cc) * inv_r + &
                                            cutoff_scale * (cutoff_scale - 1) * inv_r_cutoff * inv_r_cutoff * inv_cutoff * inv_cutoff * r_en(iel, nc, no, cc) 
            end do
        end do
    end do

    do j = 1, nelec
        if (iel == j) cycle
        do k = 1, 3
            rvec_ee(k, iel, j) = xnew(k) - x(k, j)
            rvec_ee(k, j, iel) = -rvec_ee(k, iel, j)
        end do
        r = sqrt( rvec_ee(1,iel,j)**2 + rvec_ee(2,iel,j)**2 + rvec_ee(3,iel,j)**2 )
        inv_r = 1.0d0 / r

        r_ee(iel, j, 0) = 1.0d0
        r_ee_gl(iel, 1, j, 0) = 0.0d0
        r_ee_gl(iel, 2, j, 0) = 0.0d0
        r_ee_gl(iel, 3, j, 0) = 0.0d0
        r_ee_gl(iel, 4, j, 0) = 0.0d0

        r_ee(j, iel, 0) = 1.0d0
        r_ee_gl(j, 1, iel, 0) = 0.0d0
        r_ee_gl(j, 2, iel, 0) = 0.0d0
        r_ee_gl(j, 3, iel, 0) = 0.0d0
        r_ee_gl(j, 4, iel, 0) = 0.0d0

        do no = 1, maxord
            r_ee(iel, j, no) = r_ee(iel, j, no-1) * r
            r_ee(j, iel, no) = r_ee(iel, j, no)  
            r_ee_gl(iel, 1, j, no) = no * r_ee(iel, j, no-1) * inv_r * rvec_ee(1, iel, j)
            r_ee_gl(iel, 2, j, no) = no * r_ee(iel, j, no-1) * inv_r * rvec_ee(2, iel, j)
            r_ee_gl(iel, 3, j, no) = no * r_ee(iel, j, no-1) * inv_r * rvec_ee(3, iel, j)
            r_ee_gl(iel, 4, j, no) = no * r_ee(iel, j, no-1) * 2.0d0 * inv_r
            r_ee_gl(j, 1, iel, no) = -r_ee_gl(iel, 1, j, no)
            r_ee_gl(j, 2, iel, no) = -r_ee_gl(iel, 2, j, no)
            r_ee_gl(j, 3, iel, no) = -r_ee_gl(iel, 3, j, no)
        end do
        do no = 2, maxord
            r_ee_gl(iel, 4, j, no) = r_ee_gl(iel, 4, j, no) + no * (no - 1) * r_ee(iel, j, no-2)
            r_ee_gl(j, 4, iel, no) = r_ee_gl(iel, 4, j, no)
        end do
    end do
end subroutine single_rios_distances

subroutine rios_p()
    use precision_kinds, only: dp
    use system, only: nelec, ncent, nctype
    use m_backflow, only: nordb_bf, nordc_bf
    use m_backflow, only: r_en, r_ee, r_ee_gl, p, d_p
    implicit none
    integer :: l, m, nc, i, j, k

    do l = 0, nordc_bf
        do m = 0, nordc_bf - l
            do nc = 1, nctype
                do i = 1, nelec
                    do j = 1, nelec
                        p(i, nc, m, l) = r_en(j, nc, m, 2) * r_ee(j, i, l)
                        do k = 1, 4
                            d_p(i, k, nc, m, l) = r_en(j, nc, m, 2) * r_ee_gl(j, k, i, l)
                        end do
                    end do
                end do
            end do
        end do
    end do

end subroutine rios_p

subroutine rios_backflow(x, quasi_x, dquasi_dx, d2quasi_dx2, dquasi_dp)
    use precision_kinds, only: dp
    use system, only: nelec, iwctype, ncent, nctype, cent
    use optwf_control, only: ioptci, ioptjas, ioptorb, ioptbf
    use m_backflow, only: parm_bf, nparm_bf, norda_bf, nordb_bf, nordc_bf, cutoff_scale
    use m_backflow, only: r_en, rvec_en, r_ee, rvec_ee, r_ee_gl, r_en_gl, p, d_p, cutoff_deriv
    implicit none
    real(dp), dimension(3, nelec), intent(in) :: x
    real(dp), dimension(3, nelec), intent(out) :: quasi_x
    real(dp), dimension(3, nelec, 3, nelec), intent(out) :: dquasi_dx
    real(dp), dimension(3, nelec, nelec), intent(out) :: d2quasi_dx2
    real(dp), dimension(3, nelec, nparm_bf), intent(out) :: dquasi_dp
    real(dp) :: rij, rr, cutoff, b_one, a_one
    real(dp) :: f, fp(3), fpp, eta, etap(3), etapp
    real(dp) :: delta(3), delta_ril(3), delta_rjl(3)
    integer :: i, j, k, a, b, offset_ee, offset_en, offset_een, idx, C
    integer :: nc, l, m, n, tmpc, idx_phi, idx_theta, kk
    real(dp) :: ril, rjl, inv_ril, inv_rjl
    real(dp) :: inv_rij, inv_cutoff, tmp1, tmp2, cutoff1, cutoff2
    real(dp) :: fi, fpi(3), fppi(3), fj, fpj(3), fppj(3)
    real(dp) :: phi, theta, phipi(3), thetapi(3), phipj(3), thetapj(3), phippi, thetappi, phippj, thetappj

    offset_ee = 0
    if (nordb_bf .gt. 0) then
        offset_en = offset_ee + 1+nordb_bf
    else
        offset_en = offset_ee
    end if
    if (norda_bf .gt. 0) then
        offset_een = offset_en + (1+norda_bf)*nctype
    else
        offset_een = offset_en
    end if

    quasi_x = 0.0_dp
    dquasi_dx = 0.0_dp
    d2quasi_dx2 = 0.0_dp
    dquasi_dp = 0.0_dp

    
    C = cutoff_scale

    do i = 1, nelec
        do a=1,3
            quasi_x(a,i) = x(a,i)
            dquasi_dx(a,i,a,i) = 1.0d0
        enddo
    end do

    if (nordb_bf .eq. 0) goto 10
    cutoff = parm_bf(offset_ee+1)
    inv_cutoff = 1.0d0 / cutoff

    b_one = parm_bf(offset_ee+2) * cutoff/C


    do i = 1, nelec
        do j = 1, nelec
            if (i == j) cycle
            eta=0.0d0
            etap=0.0d0
            etapp=0.0d0

            delta(:) = x(:, i) - x(:, j)

            rij = sqrt(delta(1)**2 + delta(2)**2 + delta(3)**2)
            if (rij > cutoff) cycle

            f = ((cutoff - rij)*inv_cutoff)**C
            inv_rij = 1.0d0 / rij

            cutoff1 = f/((cutoff - rij)*inv_cutoff)
            cutoff2 = cutoff1/((cutoff - rij)*inv_cutoff)

            tmp1 = -C*inv_cutoff * cutoff1 * inv_rij
            do a = 1, 3
                fp(a) = tmp1 * delta(a)
            end do
            fpp = -2.0d0 * C * inv_cutoff * inv_rij * cutoff1 + &
                  + C * (C-1) * inv_cutoff * inv_cutoff * cutoff2 


            eta  = eta + b_one
            do a = 1, 3
                dquasi_dp(a,i,offset_ee+2) = dquasi_dp(a,i,offset_ee+2) + delta(a) * f * (1 + C * inv_cutoff * rij)
            end do

            rr = rij
            
            tmp2 = inv_rij*inv_rij
            do k = 1, nordb_bf
                eta = eta + parm_bf(offset_ee+k+1)*rr
                do a = 1, 3
                    dquasi_dp(a,i,offset_ee+k+1) = dquasi_dp(a,i,offset_ee+k+1) + rr * delta(a) * f
                    etap(a) = etap(a) + parm_bf(offset_ee+k+1) * delta(a) * tmp2 * k * rr
                enddo
                etapp = etapp + parm_bf(offset_ee+k+1) * k * rr * tmp2 * (k+1)
                rr = rr * rij
            end do

            tmp2 = eta * C * cutoff1 * (rij*inv_cutoff*inv_cutoff)
            do a = 1, 3
                quasi_x(a,i) = quasi_x(a,i) + eta * delta(a) * f
                ! dquasi_dp(a,i,1) = dquasi_dp(a,i,1) + eta * delta(a) * f * log((cutoff - rij)/cutoff)
                dquasi_dp(a,i,offset_ee+1) = dquasi_dp(a,i,offset_ee+1) + delta(a) * tmp2 + &
                    delta(a) * f * b_one * inv_cutoff
                do b = 1, 3
                    tmp1 = etap(b) * delta(a) * f + eta * delta(a) * fp(b)
                    dquasi_dx(a,i,b,i) = dquasi_dx(a,i,b,i) + tmp1
                    dquasi_dx(a,i,b,j) = dquasi_dx(a,i,b,j) - tmp1

                    tmp1 = 2.0d0 * etap(b) * delta(a) * fp(b)
                    d2quasi_dx2(a,i,i) = d2quasi_dx2(a,i,i) + tmp1
                    d2quasi_dx2(a,i,j) = d2quasi_dx2(a,i,j) + tmp1
                end do
                tmp1 = eta * f
                dquasi_dx(a,i,a,i) = dquasi_dx(a,i,a,i) + tmp1
                dquasi_dx(a,i,a,j) = dquasi_dx(a,i,a,j) - tmp1

                tmp1 = etapp * delta(a) * f + eta * delta(a) * fpp
                d2quasi_dx2(a,i,i) = d2quasi_dx2(a,i,i) + tmp1
                d2quasi_dx2(a,i,j) = d2quasi_dx2(a,i,j) + tmp1

                tmp1 = 2 * (eta * fp(a) + etap(a) * f)
                d2quasi_dx2(a,i,i) = d2quasi_dx2(a,i,i) + tmp1
                d2quasi_dx2(a,i,j) = d2quasi_dx2(a,i,j) + tmp1

            end do

        
        end do

    end do

10  continue

if (norda_bf .eq. 0) goto 20

    do j = 1, ncent
        idx = (iwctype(j)-1)*(norda_bf+1)
        cutoff = parm_bf(offset_en+idx+1)
        inv_cutoff = 1.0d0 / cutoff

        a_one = parm_bf(offset_en+idx+2) * cutoff/C

        do i = 1, nelec

            delta(:) = x(:, i) - cent(:, j)

            rij = sqrt(delta(1)**2 + delta(2)**2 + delta(3)**2)
            if (rij > cutoff) cycle
            inv_rij = 1.0d0 / rij

            eta=0.0d0
            etap=0.0d0
            etapp=0.0d0

            f = ((cutoff - rij)*inv_cutoff)**C

            cutoff1 = f/((cutoff - rij)*inv_cutoff)
            cutoff2 = cutoff1/((cutoff - rij)*inv_cutoff)

            tmp1 = -C*inv_cutoff * cutoff1 * inv_rij
            do a = 1, 3
                fp(a) = tmp1 * delta(a)
            end do
            fpp = -2.0d0 * C * inv_cutoff * inv_rij * cutoff1 + &
                  + C * (C-1) * inv_cutoff * inv_cutoff * cutoff2 

            eta  = eta + a_one
            do a = 1, 3
                dquasi_dp(a,i,offset_en + idx+2) = dquasi_dp(a,i,offset_en + idx+2) + delta(a) * f * (1 + C * inv_cutoff * rij)
            end do

            rr = rij

            tmp1 = inv_rij*inv_rij
            do k = 1, norda_bf
                eta = eta + parm_bf(offset_en + idx+ k+1)*rr
                do a = 1, 3
                    dquasi_dp(a,i,offset_en + idx+ k+1) = dquasi_dp(a,i,offset_en + idx+ k+1) + rr * delta(a) * f
                    etap(a) = etap(a) + parm_bf(offset_en + idx+ k+1) * k * rr * delta(a)*tmp1
                enddo
                etapp = etapp + parm_bf(offset_en + idx+ k+1) * k * rr * tmp1 * (k+1)
                rr = rr * rij
            end do


            tmp1 = eta * C * cutoff1 * (rij*inv_cutoff*inv_cutoff)
            do a = 1, 3 
                quasi_x(a,i) = quasi_x(a,i) + eta * delta(a) * f
                dquasi_dp(a,i,offset_en + idx+1) = dquasi_dp(a,i,offset_en + idx+1) + &
                    delta(a) * tmp1 + &
                    delta(a) * f * a_one*inv_cutoff
                do b = 1, 3
                    dquasi_dx(a,i,b,i) = dquasi_dx(a,i,b,i) + (&
                        etap(b) * delta(a) * f + eta * delta(a) * fp(b) )

                    d2quasi_dx2(a,i,i) = d2quasi_dx2(a,i,i) + 2.0d0 * etap(b) * delta(a) * fp(b) 

                end do
                dquasi_dx(a,i,a,i) = dquasi_dx(a,i,a,i) + eta * f

                d2quasi_dx2(a,i,i) = d2quasi_dx2(a,i,i) + etapp * delta(a) * f + eta * delta(a) * fpp

                d2quasi_dx2(a,i,i) = d2quasi_dx2(a,i,i) + 2 * (eta * fp(a) + etap(a) * f)

            end do
            
        end do
    end do

20  continue
    if (nordc_bf .eq. 0) return

    tmpc = 0
    do l = 0, nordc_bf
        do m = 0, nordc_bf - l
            do n = 0, nordc_bf - l - m
                tmpc = tmpc + 1
            end do
        end do
    end do

    call rios_distances(x)

     do i = 1, nelec
        do j = 1, nelec
            if (i == j) cycle
            do nc = 1, ncent
                idx_phi = (iwctype(nc)-1)*(tmpc + 1) + offset_een
                idx_theta = (tmpc + 1)*nctype + (iwctype(nc)-1)*(tmpc) + offset_een

                cutoff = parm_bf(idx_phi+1)
                inv_cutoff = 1.0d0 / cutoff

                if (r_en(i,nc,0,2) < 0 .or. r_en(j,nc,0,2) < 0) cycle

                phi = 0.00d0
                phipi = 0.0d0
                phipj = 0.0d0
                phippi = 0.0d0
                phippj = 0.0d0
                theta = 0.00d0
                thetapi = 0.0d0
                thetapj = 0.0d0
                thetappi = 0.0d0
                thetappj = 0.0d0


                k = idx_phi+2
                kk = idx_theta+1
                do l = 0, nordc_bf
                    do m = 0, nordc_bf - l
                        do n = 0, nordc_bf - l - m
                            phi = phi + parm_bf(k) * r_en(i,nc,l,2) * r_en(j,nc,m,2) * r_ee(i,j,n) 
                            theta = theta + parm_bf(kk) * r_en(i,nc,l,2) * r_en(j,nc,m,2) * r_ee(i,j,n) 
                            do a = 1, 3
                                tmp1 = r_en(i,nc,l,2) * r_en(j,nc,m,2) * r_ee(i,j,n) 
                                dquasi_dp(a,i,k) = dquasi_dp(a,i,k) + tmp1 * rvec_ee(a,i,j)
                                dquasi_dp(a,i,kk) = dquasi_dp(a,i,kk) + tmp1 * rvec_en(a,i,nc)

                                tmp1 = r_en_gl(i,a,nc,l,2) * r_en(j,nc,m,2) * r_ee(i,j,n)
                                phipi(a) = phipi(a) + parm_bf(k) * tmp1
                                thetapi(a) = thetapi(a) + parm_bf(kk) * tmp1

                                tmp1 = r_en(i,nc,l,2) * r_en_gl(j,a,nc,m,2) * r_ee(i,j,n)
                                phipj(a) = phipj(a) + parm_bf(k) * tmp1
                                thetapj(a) = thetapj(a) + parm_bf(kk) * tmp1

                                tmp1 = r_en(i,nc,l,2) * r_en(j,nc,m,2) * r_ee_gl(i,a,j,n) 
                                phipi(a) = phipi(a) + parm_bf(k) * tmp1
                                phipj(a) = phipj(a) - parm_bf(k) * tmp1
                                thetapi(a) = thetapi(a) + parm_bf(kk) * tmp1
                                thetapj(a) = thetapj(a) - parm_bf(kk) * tmp1

                                tmp1 = 2.0d0 * r_en_gl(i,a,nc,l,2) * r_en(j,nc,m,2) * r_ee_gl(i,a,j,n)
                                phippi = phippi + parm_bf(k) * tmp1
                                thetappi = thetappi + parm_bf(kk) * tmp1

                                tmp1 = 2.0d0 * r_en(i,nc,l,2) * r_en_gl(j,a,nc,m,2) * r_ee_gl(i,a,j,n)
                                phippj = phippj - parm_bf(k) * tmp1
                                thetappj = thetappj - parm_bf(kk) * tmp1
                            enddo

                            tmp1 = r_en_gl(i,4,nc,l,2) * r_en(j,nc,m,2) * r_ee(i,j,n)
                            phippi = phippi + parm_bf(k) * tmp1
                            thetappi = thetappi + parm_bf(kk) * tmp1

                            tmp1 = r_en(i,nc,l,2) * r_en_gl(j,4,nc,m,2) * r_ee(i,j,n) 
                            phippj = phippj + parm_bf(k) * tmp1
                            thetappj = thetappj + parm_bf(kk) * tmp1

                            tmp1 = r_en(i,nc,l,2) * r_en(j,nc,m,2) * r_ee_gl(i,4,j,n)
                            phippi = phippi + parm_bf(k) * tmp1
                            thetappi = thetappi + parm_bf(kk) * tmp1
                            phippj = phippj + parm_bf(k) * tmp1
                            thetappj = thetappj + parm_bf(kk) * tmp1

                            k = k + 1
                            kk = kk + 1
                        end do
                    end do
                end do

 
                do a = 1, 3 
                    quasi_x(a,i) = quasi_x(a,i) + phi * rvec_ee(a,i,j) + theta * rvec_en(a,i,nc)
                    dquasi_dp(a,i,idx_phi+1) = dquasi_dp(a,i,idx_phi+1) !+ phi * rvec_ee(a,i,j) * (cutoff_deriv(i,nc) + cutoff_deriv(j,nc)) &
                                                                        !+ theta * rvec_en(a,i,nc) * (cutoff_deriv(i,nc) + cutoff_deriv(j,nc))

                    do b = 1, 3
                        dquasi_dx(a,i,b,i) = dquasi_dx(a,i,b,i) + phipi(b) * rvec_ee(a,i,j) + thetapi(b) * rvec_en(a,i,nc)                        
                        dquasi_dx(a,i,b,j) = dquasi_dx(a,i,b,j) + phipj(b) * rvec_ee(a,i,j) + thetapj(b) * rvec_en(a,i,nc)
                    end do
                    dquasi_dx(a,i,a,i) = dquasi_dx(a,i,a,i) + phi + theta 
                    dquasi_dx(a,i,a,j) = dquasi_dx(a,i,a,j) - phi 

                    d2quasi_dx2(a,i,i) = d2quasi_dx2(a,i,i) + phippi * rvec_ee(a,i,j) + thetappi * rvec_en(a,i,nc)
                    d2quasi_dx2(a,i,j) = d2quasi_dx2(a,i,j) + phippj * rvec_ee(a,i,j) + thetappj * rvec_en(a,i,nc)

                    d2quasi_dx2(a,i,i) = d2quasi_dx2(a,i,i) + 2.0d0 * phipi(a) + 2.0d0 * thetapi(a)
                    d2quasi_dx2(a,i,j) = d2quasi_dx2(a,i,j) - 2.0d0 * phipj(a)
                end do
            end do
        end do
    end do



end subroutine rios_backflow


subroutine single_rios_backflow(iel, xold, xnew, quasi_x_new, dquasi_dx_new, d2quasi_dx2_new, indices)
    use precision_kinds, only: dp
    use system, only: nelec, iwctype, ncent, nctype, cent
    use optwf_control, only: ioptci, ioptjas, ioptorb, ioptbf
    use m_backflow, only: parm_bf, nparm_bf, norda_bf, nordb_bf, nordc_bf, cutoff_scale
    use m_backflow, only: quasi_x, dquasi_dx, d2quasi_dx2, r_ee, rvec_ee, r_en, rvec_en, r_ee_gl, r_en_gl
    implicit none
    integer, intent(in) :: iel
    real(dp), dimension(3, nelec), intent(in) :: xold
    real(dp), dimension(3), intent(in) :: xnew
    real(dp), dimension(3, nelec), intent(out) :: quasi_x_new
    real(dp), dimension(3, nelec, 3, nelec), intent(out) :: dquasi_dx_new
    real(dp), dimension(3, nelec, nelec), intent(out) :: d2quasi_dx2_new
    integer, dimension(nelec), intent(out) :: indices
    real(dp) :: rij, rr, cutoff, b_one, a_one
    real(dp) :: f, fp(3), fpp, eta, etap(3), etapp
    real(dp) :: delta(3)
    integer :: i, j, k, a, b, offset_ee, offset_en, offset_een, idx, C
    integer :: nc, l, m, n, tmpc, idx_phi, idx_theta, kk
    real(dp) :: inv_rij, inv_cutoff, tmp1, tmp2
    real(dp) :: cutoff1, cutoff2
    real(dp) :: fi, fpi(3), fppi(3), fj, fpj(3), fppj(3)
    real(dp) :: phi, theta, phipi(3), thetapi(3), phipj(3), thetapj(3), phippi, thetappi, phippj, thetappj
    real(dp) :: xtemp(3, nelec)


    offset_ee = 0
    if (nordb_bf .gt. 0) then
        offset_en = offset_ee + 1+nordb_bf
    else
        offset_en = offset_ee
    end if
    if (norda_bf .gt. 0) then
        offset_een = offset_en + (1+norda_bf)*nctype
    else
        offset_een = offset_en
    end if

    quasi_x_new = quasi_x
    dquasi_dx_new = dquasi_dx
    d2quasi_dx2_new = d2quasi_dx2

    indices = 0

    do k = 1, 3
        quasi_x_new(k, iel) = quasi_x_new(k, iel) - xold(k, iel) + xnew(k)
    end do
    indices(iel) = iel
    
    C = cutoff_scale

    if (nordb_bf .eq. 0) goto 10
    cutoff = parm_bf(offset_ee+1)
    inv_cutoff = 1.0d0 / cutoff

    b_one = parm_bf(offset_ee+2) * cutoff/C


    do j = 1, nelec
        if (j == iel) cycle

        delta(:) = xold(:, iel) - xold(:, j)
        rij = sqrt(delta(1)**2 + delta(2)**2 + delta(3)**2)

        if (rij < cutoff) then
            indices(j) = j

            eta=0.0d0
            etap=0.0d0
            etapp=0.0d0

            f = ((cutoff - rij)*inv_cutoff)**C
            inv_rij = 1.0d0 / rij

            cutoff1 = f/((cutoff - rij)*inv_cutoff)
            cutoff2 = cutoff1/((cutoff - rij)*inv_cutoff)

            tmp1 = -C*inv_cutoff * cutoff1 * inv_rij
            do a = 1, 3
                fp(a) = tmp1 * delta(a)
            end do
            fpp = -2.0d0 * C * inv_cutoff * inv_rij * cutoff1 + &
                  + C * (C-1) * inv_cutoff * inv_cutoff * cutoff2 


            eta = eta + b_one

            rr = rij
        
            tmp1 = inv_rij*inv_rij
            do k = 1, nordb_bf
                eta = eta + parm_bf(offset_ee+k+1)*rr
                do a = 1, 3
                    etap(a) = etap(a) + parm_bf(offset_ee+k+1) * delta(a) * tmp1 * k * rr
                enddo
                etapp = etapp + parm_bf(offset_ee+k+1) * k * rr * tmp1 * (k+1)
                rr = rr * rij
            end do


            do a = 1, 3
                quasi_x_new(a,iel) = quasi_x_new(a,iel) - eta * delta(a) * f 
                quasi_x_new(a,j) = quasi_x_new(a,j) + eta * delta(a) * f 

                do b = 1, 3
                    tmp1 = etap(b) * delta(a) * f + eta * delta(a) * fp(b)
                    dquasi_dx_new(a,iel,b,iel) = dquasi_dx_new(a,iel,b,iel) - tmp1
                    dquasi_dx_new(a,iel,b,j) = dquasi_dx_new(a,iel,b,j) + tmp1
                    dquasi_dx_new(a,j,b,iel) = dquasi_dx_new(a,j,b,iel) + tmp1
                    dquasi_dx_new(a,j,b,j) = dquasi_dx_new(a,j,b,j) - tmp1

                    tmp1 = 2.0d0 * etap(b) * delta(a) * fp(b)
                    d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) - tmp1
                    d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) - tmp1
                    d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) + tmp1
                    d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) + tmp1
                end do
                tmp1 = eta * f
                dquasi_dx_new(a,iel,a,iel) = dquasi_dx_new(a,iel,a,iel) - tmp1
                dquasi_dx_new(a,iel,a,j) = dquasi_dx_new(a,iel,a,j) + tmp1
                dquasi_dx_new(a,j,a,iel) = dquasi_dx_new(a,j,a,iel) + tmp1
                dquasi_dx_new(a,j,a,j) = dquasi_dx_new(a,j,a,j) - tmp1

                tmp1 = etapp * delta(a) * f + eta * delta(a) * fpp
                d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) - tmp1
                d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) - tmp1
                d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) + tmp1
                d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) + tmp1

                tmp1 = 2 * (eta * fp(a) + etap(a) * f)
                d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) - tmp1
                d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) - tmp1
                d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) + tmp1
                d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) + tmp1
            end do
        endif

        delta(:) = xnew(:) - xold(:, j)
        rij = sqrt(delta(1)**2 + delta(2)**2 + delta(3)**2)

        if (rij < cutoff) then
            indices(j) = j

            eta=0.0d0
            etap=0.0d0
            etapp=0.0d0

            f = ((cutoff - rij)*inv_cutoff)**C
            inv_rij = 1.0d0 / rij

            cutoff1 = f/((cutoff - rij)*inv_cutoff)
            cutoff2 = cutoff1/((cutoff - rij)*inv_cutoff)

            tmp1 = -C*inv_cutoff * cutoff1 * inv_rij
            do a = 1, 3
                fp(a) = tmp1 * delta(a)
            end do
            fpp = -2.0d0 * C * inv_cutoff * inv_rij * cutoff1 + &
                  + C * (C-1) * inv_cutoff * inv_cutoff * cutoff2 


            eta = eta + b_one

            rr = rij
        
            tmp1 = inv_rij*inv_rij
            do k = 1, nordb_bf
                eta = eta + parm_bf(offset_ee+k+1)*rr
                do a = 1, 3
                    etap(a) = etap(a) + parm_bf(offset_ee+k+1) * delta(a) * tmp1 * k * rr
                enddo
                etapp = etapp + parm_bf(offset_ee+k+1) * k * rr * tmp1 * (k+1)
                rr = rr * rij
            end do


            do a = 1, 3
                quasi_x_new(a,iel) = quasi_x_new(a,iel) + eta * delta(a) * f 
                quasi_x_new(a,j) = quasi_x_new(a,j) - eta * delta(a) * f 

                do b = 1, 3
                    tmp1 = etap(b) * delta(a) * f + eta * delta(a) * fp(b)
                    dquasi_dx_new(a,iel,b,iel) = dquasi_dx_new(a,iel,b,iel) + tmp1
                    dquasi_dx_new(a,iel,b,j) = dquasi_dx_new(a,iel,b,j) - tmp1
                    dquasi_dx_new(a,j,b,iel) = dquasi_dx_new(a,j,b,iel) - tmp1
                    dquasi_dx_new(a,j,b,j) = dquasi_dx_new(a,j,b,j) + tmp1

                    tmp1 = 2.0d0 * etap(b) * delta(a) * fp(b)
                    d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + tmp1
                    d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) + tmp1
                    d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) - tmp1
                    d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) - tmp1
                end do
                tmp1 = eta * f
                dquasi_dx_new(a,iel,a,iel) = dquasi_dx_new(a,iel,a,iel) + tmp1
                dquasi_dx_new(a,iel,a,j) = dquasi_dx_new(a,iel,a,j) - tmp1
                dquasi_dx_new(a,j,a,iel) = dquasi_dx_new(a,j,a,iel) - tmp1
                dquasi_dx_new(a,j,a,j) = dquasi_dx_new(a,j,a,j) + tmp1

                tmp1 = etapp * delta(a) * f + eta * delta(a) * fpp
                d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + tmp1
                d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) + tmp1
                d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) - tmp1
                d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) - tmp1

                tmp1 = 2 * (eta * fp(a) + etap(a) * f)
                d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + tmp1
                d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) + tmp1
                d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) - tmp1
                d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) - tmp1
            end do
        endif
    
    end do


10  continue

    if (norda_bf .eq. 0) goto 20

    do j = 1, ncent
        idx = (iwctype(j)-1)*(norda_bf+1)
        cutoff = parm_bf(offset_en+idx+1)
        inv_cutoff = 1.0d0 / cutoff

        a_one = parm_bf(offset_en+idx+2) * cutoff/C


        delta(:) = xold(:, iel) - cent(:, j)
        rij = sqrt(delta(1)**2 + delta(2)**2 + delta(3)**2)

        if (rij < cutoff) then
            indices(iel) = iel
            inv_rij = 1.0d0 / rij

            eta=0.0d0
            etap=0.0d0
            etapp=0.0d0

            f = ((cutoff - rij)*inv_cutoff)**C

            cutoff1 = f/((cutoff - rij)*inv_cutoff)
            cutoff2 = cutoff1/((cutoff - rij)*inv_cutoff)

            tmp1 = -C*inv_cutoff * cutoff1 * inv_rij
            do a = 1, 3
                fp(a) = tmp1 * delta(a)
            end do
            fpp = -2.0d0 * C * inv_cutoff * inv_rij * cutoff1 + &
                  + C * (C-1) * inv_cutoff * inv_cutoff * cutoff2 

            eta  = eta + a_one

            rr = rij

            tmp1 = inv_rij*inv_rij
            do k = 1, norda_bf
                eta = eta + parm_bf(offset_en + idx+ k+1)*rr
                do a = 1, 3
                    etap(a) = etap(a) + parm_bf(offset_en + idx+ k+1) * k * rr * delta(a)*tmp1
                enddo
                etapp = etapp + parm_bf(offset_en + idx+ k+1) * k * rr * tmp1 * (k+1)
                rr = rr * rij
            end do

            do a = 1, 3 
                quasi_x_new(a,iel) = quasi_x_new(a,iel) - eta * delta(a) * f
                do b = 1, 3
                    dquasi_dx_new(a,iel,b,iel) = dquasi_dx_new(a,iel,b,iel) - (&
                        etap(b) * delta(a) * f + eta * delta(a) * fp(b) )

                    d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) - (&
                         2.0d0 * etap(b) * delta(a) * fp(b)  )
                end do
                dquasi_dx_new(a,iel,a,iel) = dquasi_dx_new(a,iel,a,iel) - eta * f
                d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) - etapp * delta(a) * f - eta * delta(a) * fpp
                d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) - 2 * (eta * fp(a) + etap(a) * f)
            end do
        endif
        
        delta(:) = xnew(:) - cent(:, j)
        rij = sqrt(delta(1)**2 + delta(2)**2 + delta(3)**2)
        if (rij < cutoff) then
            indices(iel) = iel
            inv_rij = 1.0d0 / rij

            eta=0.0d0
            etap=0.0d0
            etapp=0.0d0

            f = ((cutoff - rij)*inv_cutoff)**C

            cutoff1 = f/((cutoff - rij)*inv_cutoff)
            cutoff2 = cutoff1/((cutoff - rij)*inv_cutoff)

            tmp1 = -C*inv_cutoff * cutoff1 * inv_rij
            do a = 1, 3
                fp(a) = tmp1 * delta(a)
            end do
            fpp = -2.0d0 * C * inv_cutoff * inv_rij * cutoff1 + &
                  + C * (C-1) * inv_cutoff * inv_cutoff * cutoff2 

            eta  = eta + a_one

            rr = rij

            tmp1 = inv_rij*inv_rij
            do k = 1, norda_bf
                eta = eta + parm_bf(offset_en + idx+ k+1)*rr
                do a = 1, 3
                    etap(a) = etap(a) + parm_bf(offset_en + idx+ k+1) * k * rr * delta(a)*tmp1
                enddo
                etapp = etapp + parm_bf(offset_en + idx+ k+1) * k * rr * tmp1 * (k+1)
                rr = rr * rij
            end do

            do a = 1, 3 
                quasi_x_new(a,iel) = quasi_x_new(a,iel) + eta * delta(a) * f
                do b = 1, 3
                    dquasi_dx_new(a,iel,b,iel) = dquasi_dx_new(a,iel,b,iel) + (&
                        etap(b) * delta(a) * f + eta * delta(a) * fp(b) )

                    d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + (&
                        2.0d0 * etap(b) * delta(a) * fp(b)  )
                end do
                dquasi_dx_new(a,iel,a,iel) = dquasi_dx_new(a,iel,a,iel) + eta * f
                d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + etapp * delta(a) * f + eta * delta(a) * fpp
                d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + 2 * (eta * fp(a) + etap(a) * f)
            end do
        endif

    end do


20  continue
    if (nordc_bf .eq. 0) return

    tmpc = 0
    do l = 0, nordc_bf
        do m = 0, nordc_bf - l
            do n = 0, nordc_bf - l - m
                tmpc = tmpc + 1
            end do
        end do
    end do

    !call single_rios_distances(xold, xold(:, iel), iel)
     call rios_distances(xold)
    do j = 1, nelec
        if (iel == j) cycle
        do nc = 1, ncent
            idx_phi = (iwctype(nc)-1)*(tmpc + 1) + offset_een
            idx_theta =  (tmpc + 1)*nctype + (iwctype(nc)-1)*(tmpc) + offset_een

            cutoff = parm_bf(idx_phi+1)
            inv_cutoff = 1.0d0 / cutoff

            if (r_en(iel,nc,0,2) > 0 .and. r_en(j,nc,0,2) > 0) then

                indices(j) = j

                phi = 0.00d0
                phipi = 0.0d0
                phipj = 0.0d0
                phippi = 0.0d0
                phippj = 0.0d0
                theta = 0.00d0
                thetapi = 0.0d0
                thetapj = 0.0d0
                thetappi = 0.0d0
                thetappj = 0.0d0

                k = idx_phi+2
                kk = idx_theta+1
                do l = 0, nordc_bf
                    do m = 0, nordc_bf - l
                        do n = 0, nordc_bf - l - m
                            phi = phi + parm_bf(k) * r_en(iel,nc,l,2) * r_en(j,nc,m,2) * r_ee(iel,j,n) 
                            theta = theta + parm_bf(kk) * r_en(iel,nc,l,2) * r_en(j,nc,m,2) * r_ee(iel,j,n) 
                            do a = 1, 3
                                tmp1 = r_en_gl(iel,a,nc,l,2) * r_en(j,nc,m,2) * r_ee(iel,j,n)
                                phipi(a) = phipi(a) + parm_bf(k) * tmp1
                                thetapi(a) = thetapi(a) + parm_bf(kk) * tmp1

                                tmp1 = r_en(iel,nc,l,2) * r_en_gl(j,a,nc,m,2) * r_ee(iel,j,n)
                                phipj(a) = phipj(a) + parm_bf(k) * tmp1
                                thetapj(a) = thetapj(a) + parm_bf(kk) * tmp1

                                tmp1 = r_en(iel,nc,l,2) * r_en(j,nc,m,2) * r_ee_gl(iel,a,j,n) 
                                phipi(a) = phipi(a) + parm_bf(k) * tmp1
                                phipj(a) = phipj(a) - parm_bf(k) * tmp1
                                thetapi(a) = thetapi(a) + parm_bf(kk) * tmp1
                                thetapj(a) = thetapj(a) - parm_bf(kk) * tmp1

                                tmp1 = 2.0d0 * r_en_gl(iel,a,nc,l,2) * r_en(j,nc,m,2) * r_ee_gl(iel,a,j,n)
                                phippi = phippi + parm_bf(k) * tmp1
                                thetappi = thetappi + parm_bf(kk) * tmp1

                                tmp1 = 2.0d0 * r_en(iel,nc,l,2) * r_en_gl(j,a,nc,m,2) * r_ee_gl(iel,a,j,n)
                                phippj = phippj - parm_bf(k) * tmp1
                                thetappj = thetappj - parm_bf(kk) * tmp1
                            enddo

                            tmp1 = r_en_gl(iel,4,nc,l,2) * r_en(j,nc,m,2) * r_ee(iel,j,n)
                            phippi = phippi + parm_bf(k) * tmp1
                            thetappi = thetappi + parm_bf(kk) * tmp1

                            tmp1 = r_en(iel,nc,l,2) * r_en_gl(j,4,nc,m,2) * r_ee(iel,j,n) 
                            phippj = phippj + parm_bf(k) * tmp1
                            thetappj = thetappj + parm_bf(kk) * tmp1

                            tmp1 = r_en(iel,nc,l,2) * r_en(j,nc,m,2) * r_ee_gl(iel,4,j,n)
                            phippi = phippi + parm_bf(k) * tmp1
                            thetappi = thetappi + parm_bf(kk) * tmp1
                            phippj = phippj + parm_bf(k) * tmp1
                            thetappj = thetappj + parm_bf(kk) * tmp1

                            k = k + 1
                            kk = kk + 1
                        end do
                    end do
                end do
                do a = 1, 3 
                    quasi_x_new(a,iel) = quasi_x_new(a,iel) - phi * rvec_ee(a,iel,j) - theta * rvec_en(a,iel,nc)

                    do b = 1, 3
                        dquasi_dx_new(a,iel,b,iel) = dquasi_dx_new(a,iel,b,iel) - phipi(b) * rvec_ee(a,iel,j) - thetapi(b) * rvec_en(a,iel,nc)                        
                        dquasi_dx_new(a,iel,b,j) = dquasi_dx_new(a,iel,b,j) - phipj(b) * rvec_ee(a,iel,j) - thetapj(b) * rvec_en(a,iel,nc)
                    end do
                    dquasi_dx_new(a,iel,a,iel) = dquasi_dx_new(a,iel,a,iel) - phi - theta 
                    dquasi_dx_new(a,iel,a,j) = dquasi_dx_new(a,iel,a,j) + phi 
  
                    d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) - phippi * rvec_ee(a,iel,j) - thetappi * rvec_en(a,iel,nc)
                    d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) - phippj * rvec_ee(a,iel,j) - thetappj * rvec_en(a,iel,nc)

                    d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) - 2.0d0 * phipi(a) - 2.0d0 * thetapi(a)
                    d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) + 2.0d0 * phipj(a)
                end do

                ! Now account for the contribution where j is the central electron
                phi = 0.0d0
                phipi = 0.0d0
                phipj = 0.0d0
                phippi = 0.0d0
                phippj = 0.0d0
                theta = 0.0d0
                thetapi = 0.0d0
                thetapj = 0.0d0
                thetappi = 0.0d0
                thetappj = 0.0d0

                k = idx_phi+2
                kk = idx_theta+1
                do l = 0, nordc_bf
                    do m = 0, nordc_bf - l
                        do n = 0, nordc_bf - l - m
                            phi = phi + parm_bf(k) * r_en(j,nc,l,2) * r_en(iel,nc,m,2) * r_ee(j,iel,n) 
                            theta = theta + parm_bf(kk) * r_en(j,nc,l,2) * r_en(iel,nc,m,2) * r_ee(j,iel,n) 
                            do a = 1, 3
                                tmp1 = r_en_gl(j,a,nc,l,2) * r_en(iel,nc,m,2) * r_ee(j,iel,n)
                                phipi(a) = phipi(a) + parm_bf(k) * tmp1
                                thetapi(a) = thetapi(a) + parm_bf(kk) * tmp1

                                tmp1 = r_en(j,nc,l,2) * r_en_gl(iel,a,nc,m,2) * r_ee(j,iel,n)
                                phipj(a) = phipj(a) + parm_bf(k) * tmp1
                                thetapj(a) = thetapj(a) + parm_bf(kk) * tmp1

                                tmp1 = r_en(j,nc,l,2) * r_en(iel,nc,m,2) * r_ee_gl(j,a,iel,n) 
                                phipi(a) = phipi(a) + parm_bf(k) * tmp1
                                phipj(a) = phipj(a) - parm_bf(k) * tmp1
                                thetapi(a) = thetapi(a) + parm_bf(kk) * tmp1
                                thetapj(a) = thetapj(a) - parm_bf(kk) * tmp1

                                tmp1 = 2.0d0 * r_en_gl(j,a,nc,l,2) * r_en(iel,nc,m,2) * r_ee_gl(j,a,iel,n)
                                phippi = phippi + parm_bf(k) * tmp1
                                thetappi = thetappi + parm_bf(kk) * tmp1

                                tmp1 = 2.0d0 * r_en(j,nc,l,2) * r_en_gl(iel,a,nc,m,2) * r_ee_gl(j,a,iel,n)
                                phippj = phippj - parm_bf(k) * tmp1
                                thetappj = thetappj - parm_bf(kk) * tmp1
                            enddo

                            tmp1 = r_en_gl(j,4,nc,l,2) * r_en(iel,nc,m,2) * r_ee(j,iel,n)
                            phippi = phippi + parm_bf(k) * tmp1
                            thetappi = thetappi + parm_bf(kk) * tmp1

                            tmp1 = r_en(j,nc,l,2) * r_en_gl(iel,4,nc,m,2) * r_ee(j,iel,n) 
                            phippj = phippj + parm_bf(k) * tmp1
                            thetappj = thetappj + parm_bf(kk) * tmp1

                            tmp1 = r_en(j,nc,l,2) * r_en(iel,nc,m,2) * r_ee_gl(j,4,iel,n)
                            phippi = phippi + parm_bf(k) * tmp1
                            thetappi = thetappi + parm_bf(kk) * tmp1
                            phippj = phippj + parm_bf(k) * tmp1
                            thetappj = thetappj + parm_bf(kk) * tmp1

                            k = k + 1
                            kk = kk + 1
                        end do
                    end do
                end do

                do a = 1, 3 
                    quasi_x_new(a,j) = quasi_x_new(a,j) - phi * rvec_ee(a,j,iel) - theta * rvec_en(a,j,nc)

                    do b = 1, 3
                        dquasi_dx_new(a,j,b,j) = dquasi_dx_new(a,j,b,j) - phipi(b) * rvec_ee(a,j,iel) - thetapi(b) * rvec_en(a,j,nc)                        
                        dquasi_dx_new(a,j,b,iel) = dquasi_dx_new(a,j,b,iel) - phipj(b) * rvec_ee(a,j,iel) - thetapj(b) * rvec_en(a,j,nc)
                    end do
                    dquasi_dx_new(a,j,a,j) = dquasi_dx_new(a,j,a,j) - phi - theta 
                    dquasi_dx_new(a,j,a,iel) = dquasi_dx_new(a,j,a,iel) + phi

                    d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) - phippi * rvec_ee(a,j,iel) - thetappi * rvec_en(a,j,nc)
                    d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) - phippj * rvec_ee(a,j,iel) - thetappj * rvec_en(a,j,nc)

                    d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) - 2.0d0 * phipi(a) - 2.0d0 * thetapi(a)
                    d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) + 2.0d0 * phipj(a)
                end do
            end if
        end do
    end do

    
    !call single_rios_distances(xold, xnew, iel)
    xtemp = xold
    xtemp(:, iel) = xnew
    call rios_distances(xtemp)


    do j = 1, nelec
        if (iel == j) cycle
        do nc = 1, ncent
            idx_phi = (iwctype(nc)-1)*(tmpc + 1) + offset_een
            idx_theta = (tmpc + 1)*nctype + (iwctype(nc)-1)*(tmpc) + offset_een

            cutoff = parm_bf(idx_phi+1)
            inv_cutoff = 1.0d0 / cutoff

            if (r_en(iel,nc,0,2) > 0 .and. r_en(j,nc,0,2) > 0) then

                indices(j) = j

                phi = 0.00d0
                phipi = 0.0d0
                phipj = 0.0d0
                phippi = 0.0d0
                phippj = 0.0d0
                theta = 0.00d0
                thetapi = 0.0d0
                thetapj = 0.0d0
                thetappi = 0.0d0
                thetappj = 0.0d0

                k = idx_phi+2
                kk = idx_theta+1
                do l = 0, nordc_bf
                    do m = 0, nordc_bf - l
                        do n = 0, nordc_bf - l - m
                            phi = phi + parm_bf(k) * r_en(iel,nc,l,2) * r_en(j,nc,m,2) * r_ee(iel,j,n) 
                            theta = theta + parm_bf(kk) * r_en(iel,nc,l,2) * r_en(j,nc,m,2) * r_ee(iel,j,n) 
                            do a = 1, 3
                                tmp1 = r_en_gl(iel,a,nc,l,2) * r_en(j,nc,m,2) * r_ee(iel,j,n)
                                phipi(a) = phipi(a) + parm_bf(k) * tmp1
                                thetapi(a) = thetapi(a) + parm_bf(kk) * tmp1

                                tmp1 = r_en(iel,nc,l,2) * r_en_gl(j,a,nc,m,2) * r_ee(iel,j,n)
                                phipj(a) = phipj(a) + parm_bf(k) * tmp1
                                thetapj(a) = thetapj(a) + parm_bf(kk) * tmp1

                                tmp1 = r_en(iel,nc,l,2) * r_en(j,nc,m,2) * r_ee_gl(iel,a,j,n) 
                                phipi(a) = phipi(a) + parm_bf(k) * tmp1
                                phipj(a) = phipj(a) - parm_bf(k) * tmp1
                                thetapi(a) = thetapi(a) + parm_bf(kk) * tmp1
                                thetapj(a) = thetapj(a) - parm_bf(kk) * tmp1

                                tmp1 = 2.0d0 * r_en_gl(iel,a,nc,l,2) * r_en(j,nc,m,2) * r_ee_gl(iel,a,j,n)
                                phippi = phippi + parm_bf(k) * tmp1
                                thetappi = thetappi + parm_bf(kk) * tmp1

                                tmp1 = 2.0d0 * r_en(iel,nc,l,2) * r_en_gl(j,a,nc,m,2) * r_ee_gl(iel,a,j,n)
                                phippj = phippj - parm_bf(k) * tmp1
                                thetappj = thetappj - parm_bf(kk) * tmp1
                            enddo

                            tmp1 = r_en_gl(iel,4,nc,l,2) * r_en(j,nc,m,2) * r_ee(iel,j,n)
                            phippi = phippi + parm_bf(k) * tmp1
                            thetappi = thetappi + parm_bf(kk) * tmp1

                            tmp1 = r_en(iel,nc,l,2) * r_en_gl(j,4,nc,m,2) * r_ee(iel,j,n) 
                            phippj = phippj + parm_bf(k) * tmp1
                            thetappj = thetappj + parm_bf(kk) * tmp1

                            tmp1 = r_en(iel,nc,l,2) * r_en(j,nc,m,2) * r_ee_gl(iel,4,j,n)
                            phippi = phippi + parm_bf(k) * tmp1
                            thetappi = thetappi + parm_bf(kk) * tmp1
                            phippj = phippj + parm_bf(k) * tmp1
                            thetappj = thetappj + parm_bf(kk) * tmp1

                            k = k + 1
                            kk = kk + 1
                        end do
                    end do
                end do
                do a = 1, 3 
                    quasi_x_new(a,iel) = quasi_x_new(a,iel) + phi * rvec_ee(a,iel,j) + theta * rvec_en(a,iel,nc)

                    do b = 1, 3
                        dquasi_dx_new(a,iel,b,iel) = dquasi_dx_new(a,iel,b,iel) + phipi(b) * rvec_ee(a,iel,j) + thetapi(b) * rvec_en(a,iel,nc)                        
                        dquasi_dx_new(a,iel,b,j) = dquasi_dx_new(a,iel,b,j) + phipj(b) * rvec_ee(a,iel,j) + thetapj(b) * rvec_en(a,iel,nc)
                    end do
                    dquasi_dx_new(a,iel,a,iel) = dquasi_dx_new(a,iel,a,iel) + phi + theta 
                    dquasi_dx_new(a,iel,a,j) = dquasi_dx_new(a,iel,a,j) - phi 

                    d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + phippi * rvec_ee(a,iel,j) + thetappi * rvec_en(a,iel,nc)
                    d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) + phippj * rvec_ee(a,iel,j) + thetappj * rvec_en(a,iel,nc)

                    d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + 2.0d0 * phipi(a) + 2.0d0 * thetapi(a)
                    d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) - 2.0d0 * phipj(a)
                end do

                ! Contribution with j as the central electron (new configuration)
                phi = 0.0d0
                phipi = 0.0d0
                phipj = 0.0d0
                phippi = 0.0d0
                phippj = 0.0d0
                theta = 0.0d0
                thetapi = 0.0d0
                thetapj = 0.0d0
                thetappi = 0.0d0
                thetappj = 0.0d0

                k = idx_phi+2
                kk = idx_theta+1
                do l = 0, nordc_bf
                    do m = 0, nordc_bf - l
                        do n = 0, nordc_bf - l - m
                            phi = phi + parm_bf(k) * r_en(j,nc,l,2) * r_en(iel,nc,m,2) * r_ee(j,iel,n) 
                            theta = theta + parm_bf(kk) * r_en(j,nc,l,2) * r_en(iel,nc,m,2) * r_ee(j,iel,n) 
                            do a = 1, 3
                                tmp1 = r_en_gl(j,a,nc,l,2) * r_en(iel,nc,m,2) * r_ee(j,iel,n)
                                phipi(a) = phipi(a) + parm_bf(k) * tmp1
                                thetapi(a) = thetapi(a) + parm_bf(kk) * tmp1

                                tmp1 = r_en(j,nc,l,2) * r_en_gl(iel,a,nc,m,2) * r_ee(j,iel,n)
                                phipj(a) = phipj(a) + parm_bf(k) * tmp1
                                thetapj(a) = thetapj(a) + parm_bf(kk) * tmp1

                                tmp1 = r_en(j,nc,l,2) * r_en(iel,nc,m,2) * r_ee_gl(j,a,iel,n) 
                                phipi(a) = phipi(a) + parm_bf(k) * tmp1
                                phipj(a) = phipj(a) - parm_bf(k) * tmp1
                                thetapi(a) = thetapi(a) + parm_bf(kk) * tmp1
                                thetapj(a) = thetapj(a) - parm_bf(kk) * tmp1

                                tmp1 = 2.0d0 * r_en_gl(j,a,nc,l,2) * r_en(iel,nc,m,2) * r_ee_gl(j,a,iel,n)
                                phippi = phippi + parm_bf(k) * tmp1
                                thetappi = thetappi + parm_bf(kk) * tmp1

                                tmp1 = 2.0d0 * r_en(j,nc,l,2) * r_en_gl(iel,a,nc,m,2) * r_ee_gl(j,a,iel,n)
                                phippj = phippj - parm_bf(k) * tmp1
                                thetappj = thetappj - parm_bf(kk) * tmp1
                            enddo

                            tmp1 = r_en_gl(j,4,nc,l,2) * r_en(iel,nc,m,2) * r_ee(j,iel,n)
                            phippi = phippi + parm_bf(k) * tmp1
                            thetappi = thetappi + parm_bf(kk) * tmp1

                            tmp1 = r_en(j,nc,l,2) * r_en_gl(iel,4,nc,m,2) * r_ee(j,iel,n) 
                            phippj = phippj + parm_bf(k) * tmp1
                            thetappj = thetappj + parm_bf(kk) * tmp1

                            tmp1 = r_en(j,nc,l,2) * r_en(iel,nc,m,2) * r_ee_gl(j,4,iel,n)
                            phippi = phippi + parm_bf(k) * tmp1
                            thetappi = thetappi + parm_bf(kk) * tmp1
                            phippj = phippj + parm_bf(k) * tmp1
                            thetappj = thetappj + parm_bf(kk) * tmp1

                            k = k + 1
                            kk = kk + 1
                        end do
                    end do
                end do

                do a = 1, 3 
                    quasi_x_new(a,j) = quasi_x_new(a,j) + phi * rvec_ee(a,j,iel) + theta * rvec_en(a,j,nc)

                    do b = 1, 3
                        dquasi_dx_new(a,j,b,j) = dquasi_dx_new(a,j,b,j) + phipi(b) * rvec_ee(a,j,iel) + thetapi(b) * rvec_en(a,j,nc)                        
                        dquasi_dx_new(a,j,b,iel) = dquasi_dx_new(a,j,b,iel) + phipj(b) * rvec_ee(a,j,iel) + thetapj(b) * rvec_en(a,j,nc)
                    end do
                    dquasi_dx_new(a,j,a,j) = dquasi_dx_new(a,j,a,j) + phi + theta 
                    dquasi_dx_new(a,j,a,iel) = dquasi_dx_new(a,j,a,iel) - phi

                    d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) + phippi * rvec_ee(a,j,iel) + thetappi * rvec_en(a,j,nc)
                    d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) + phippj * rvec_ee(a,j,iel) + thetappj * rvec_en(a,j,nc)

                    d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) + 2.0d0 * phipi(a) + 2.0d0 * thetapi(a)
                    d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) - 2.0d0 * phipj(a)
                end do
            end if
        end do
    end do
end subroutine single_rios_backflow


end module