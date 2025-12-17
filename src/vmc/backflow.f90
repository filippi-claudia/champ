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
        call init_rios_backflow(8,8,8)
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
    
    parm_bf = 0.001d0
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

    max_ord = max(norda_bf, nordb_bf, nordc_bf)

    call allocate_m_backflow

end subroutine init_rios_backflow_arrays


subroutine rios_distances(x)
    use precision_kinds, only: dp
    use system, only: nelec, ncent, cent, 
    use m_backflow, only: parm_bf, nparm_bf, norda_bf, nordb_bf, nordc_bf, cutoff_scale
    use m_backflow, only: r_en, rvec_en, r_ee, rvec_ee, r_ee_gl, r_en_gl, maxord
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

                r_en_gl(i, 1, nc, 0, cc) = -cutoff_scale * inv_r_cutoff * inv_r * rvec_en(1, i, nc) * cutoff * inv_r_cutoff
                r_en_gl(i, 2, nc, 0, cc) = -cutoff_scale * inv_r_cutoff * inv_r * rvec_en(2, i, nc) * cutoff * inv_r_cutoff
                r_en_gl(i, 3, nc, 0, cc) = -cutoff_scale * inv_r_cutoff * inv_r * rvec_en(3, i, nc) * cutoff * inv_r_cutoff
                r_en_gl(i, 4, nc, 0, cc) = -cutoff_scale * inv_r_cutoff * 2.0d0 * inv_r * cutoff * inv_r_cutoff + &
                                            cutoff_scale * (cutoff_scale - 1) * inv_r_cutoff * inv_r_cutoff * inv_cutoff* inv_cutoff * cutoff 


                do no = 1, maxord
                    r_en(i, nc, no, cc) = r_en(i, nc, no-1, cc) * r 
                    r_en_gl(i, 1, nc, no, cc) = no * r_en(i, nc, no-1, cc) * inv_r * rvec_en(1, i, nc) - &
                                                cutoff_scale * inv_r_cutoff * inv_r * rvec_en(1, i, nc) * inv_cutoff
                    r_en_gl(i, 2, nc, no, cc) = no * r_en(i, nc, no-1, cc) * inv_r * rvec_en(2, i, nc) - &
                                                cutoff_scale * inv_r_cutoff * inv_r * rvec_en(2, i, nc) * inv_cutoff
                    r_en_gl(i, 3, nc, no, cc) = no * r_en(i, nc, no-1, cc) * inv_r * rvec_en(3, i, nc) - &
                                                cutoff_scale * inv_r_cutoff * inv_r * rvec_en(3, i, nc) * inv_cutoff
                    r_en_gl(i, 4, nc, no, cc) = r_en(i, nc, no-1, cc) * no * 2.0d0 * (&
                                                inv_r - cutoff_scale * inv_r_cutoff * inv_cutoff) + &
                                                r_en(i, nc, no, cc) * cutoff_scale * inv_r_cutoff * inv_cutoff * (&
                                                (cutoff_scale - 1) * inv_r_cutoff * inv_cutoff - &
                                                2.0d0 * inv_r)
                end do
                do no = 2, maxord
                    r_en_gl(i, 4, nc, no, cc) = r_en_gl(i, 4, nc, no, cc) + no * (no - 1) * r_en(i, nc, no-2, cc)
                end do
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
            r_ee_gl(i, :, j, 0) = 0.0d0

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

subroutine rios_p()
    use precision_kinds, only: dp
    use system, only: nelec, ncent,, nctype
    use m_backflow, only: nordb_bf, nordc_bf
    use m_backflow, only: r_en, r_ee, r_ee_gl, p, dp
    implicit none
    integer :: l, m, nc, i, j, k

    do l = 0, nordc_bf
        do m = 0, nordc_bf - l
            do nc = 1, nctype
                do i = 1, nelec
                    do j = 1, nelec
                        p(i, nc, m, l) = r_en(j, nc, m, 2) * r_ee(j, i, l)
                        do k = 1, 4
                            dp(i, k, nc, m, l) = r_en(j, nc, m, 2) * r_ee_gl(j, k, i, l)
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
    use m_backflow, only: r_en, rvec_en, r_ee, rvec_ee, r_ee_gl, r_en_gl, p, dp
    implicit none
    real(dp), dimension(3, nelec), intent(in) :: x
    real(dp), dimension(3, nelec), intent(out) :: quasi_x
    real(dp), dimension(3, nelec, 3, nelec), intent(out) :: dquasi_dx
    real(dp), dimension(3, nelec, nelec), intent(out) :: d2quasi_dx2
    real(dp), dimension(3, nelec, nparm_bf), intent(out) :: dquasi_dp
    real(dp) :: rij, rr, cutoff, b_one, a_one
    real(dp) :: f, fp(3), fpp(3), eta, etap(3), etapp(3)
    real(dp) :: delta(3), delta_ril(3), delta_rjl(3)
    integer :: i, j, k, a, b, offset_ee, offset_en, offset_een, idx, C
    integer :: nc, l, m, n, tmpc, idx_phi, idx_theta, kk
    real(dp) :: ril, rjl, inv_ril, inv_rjl
    real(dp) :: inv_rij, inv_cutoff, tmp1, tmp2, cutoff1, cutoff2
    real(dp) :: fi, fpi(3), fppi(3), fj, fpj(3), fppj(3)
    real(dp) :: phi, theta, phipi(3), thetapi(3), phipj(3), thetapj(3), phipp(3), thetapp(3)

    offset_ee = 0
    offset_en = offset_ee + 1+nordb_bf
    offset_een = offset_en + (1+norda_bf)*nctype

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

    b_one = parm_bf(offset_ee+2) * C * inv_cutoff


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
            tmp2 = -C*inv_cutoff * ( - (C-1)*inv_cutoff ) * cutoff2 * inv_rij * inv_rij
            do a = 1, 3
                fp(a) = tmp1 * delta(a)
                fpp(a) = tmp2 * delta(a) * delta(a) + &
                         tmp1 * ( 1 - delta(a)*delta(a)*inv_rij*inv_rij )
            end do




            eta  = eta + parm_bf(offset_ee+2) + b_one*rij
            tmp1 = b_one * inv_rij
            tmp2 = b_one * inv_rij * inv_rij * inv_rij
            do a = 1, 3
                dquasi_dp(a,i,offset_ee+2) = dquasi_dp(a,i,offset_ee+2) + delta(a) * f * (1 + C * inv_cutoff * rij)
                etap(a) = etap(a) + tmp1 * delta(a)
                etapp(a) = etapp(a) + tmp1 - tmp2 *  delta(a) * delta(a)
            end do

            rr = rij*rij
            
            tmp2 = inv_rij*inv_rij
            do k = 2, nordb_bf
                eta = eta + parm_bf(offset_ee+k+1)*rr
                do a = 1, 3
                    dquasi_dp(a,i,offset_ee+k+1) = dquasi_dp(a,i,offset_ee+k+1) + rr * delta(a) * f
                    etap(a) = etap(a) + parm_bf(offset_ee+k+1) * delta(a) * tmp2 * k * rr
                    etapp(a) = etapp(a) + parm_bf(offset_ee+k+1) * tmp2 * rr * k * (1 + (k-2) * delta(a) * delta(a) * tmp2)
                enddo
                rr = rr * rij
            end do


            tmp2 = eta * C * cutoff1 * (rij*inv_cutoff*inv_cutoff)
            do a = 1, 3
                quasi_x(a,i) = quasi_x(a,i) + eta * delta(a) * f
                ! dquasi_dp(a,i,1) = dquasi_dp(a,i,1) + eta * delta(a) * f * log((cutoff - rij)/cutoff)
                dquasi_dp(a,i,offset_ee+1) = dquasi_dp(a,i,offset_ee+1) + delta(a) * tmp2 - &
                    delta(a) * f * b_one * inv_cutoff * rij
                do b = 1, 3
                    tmp1 = etap(b) * delta(a) * f + eta * delta(a) * fp(b)
                    dquasi_dx(a,i,b,i) = dquasi_dx(a,i,b,i) + tmp1
                    dquasi_dx(a,i,b,j) = dquasi_dx(a,i,b,j) - tmp1

                    tmp1 = etapp(b) * delta(a) * f + &
                           2.0d0 * etap(b) * delta(a) * fp(b) + &
                           eta * delta(a) * fpp(b)

                    d2quasi_dx2(a,i,i) = d2quasi_dx2(a,i,i) + tmp1

                    d2quasi_dx2(a,i,j) = d2quasi_dx2(a,i,j) + tmp1
                end do
                tmp1 = eta * f
                dquasi_dx(a,i,a,i) = dquasi_dx(a,i,a,i) + tmp1
                dquasi_dx(a,i,a,j) = dquasi_dx(a,i,a,j) - tmp1

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

        a_one = parm_bf(offset_en+idx+2) * C * inv_cutoff

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
            tmp2 = -C*inv_cutoff * ( - (C-1)*inv_cutoff ) * cutoff2 * inv_rij * inv_rij
            do a = 1, 3
                fp(a) = tmp1 * delta(a)
                fpp(a) = tmp2 * delta(a) * delta(a) + &
                         tmp1 * ( 1 - delta(a)*delta(a)*inv_rij*inv_rij )
            end do

            eta  = eta + parm_bf(offset_en+idx+2) + a_one*rij
            tmp1 = a_one * inv_rij
            tmp2 = a_one * inv_rij * inv_rij * inv_rij
            do a = 1, 3
                dquasi_dp(a,i,offset_en + idx+2) = dquasi_dp(a,i,offset_en + idx+2) + delta(a) * f * (1 + C * inv_cutoff * rij)
                etap(a) = etap(a) + tmp1 * delta(a)
                etapp(a) = etapp(a) + tmp1 - tmp2 *  delta(a) * delta(a)
            end do

            rr = rij*rij

            tmp1 = inv_rij*inv_rij
            do k = 2, norda_bf
                eta = eta + parm_bf(offset_en + idx+ k+1)*rr
                do a = 1, 3
                    dquasi_dp(a,i,offset_en + idx+ k+1) = dquasi_dp(a,i,offset_en + idx+ k+1) + rr * delta(a) * f
                    etap(a) = etap(a) + parm_bf(offset_en + idx+ k+1) * k * rr * delta(a)*tmp1
                    etapp(a) = etapp(a) + parm_bf(offset_en + idx+ k+1) * k * tmp1 * rr * (1 + (k-2) * delta(a) * delta(a)*tmp1)
                enddo
                rr = rr * rij
            end do


            tmp1 = eta * C * cutoff1 * (rij*inv_cutoff*inv_cutoff)
            do a = 1, 3 
                quasi_x(a,i) = quasi_x(a,i) + eta * delta(a) * f
                dquasi_dp(a,i,offset_en + idx+1) = dquasi_dp(a,i,offset_en + idx+1) + &
                    delta(a) * tmp1 - &
                    delta(a) * f * a_one*inv_cutoff * rij
                do b = 1, 3
                    dquasi_dx(a,i,b,i) = dquasi_dx(a,i,b,i) + (&
                        etap(b) * delta(a) * f + eta * delta(a) * fp(b) )

                    d2quasi_dx2(a,i,i) = d2quasi_dx2(a,i,i) + (&
                        etapp(b) * delta(a) * f + &
                        2.0d0 * etap(b) * delta(a) * fp(b) + &
                        eta * delta(a) * fpp(b) )

                end do
                dquasi_dx(a,i,a,i) = dquasi_dx(a,i,a,i) + eta * f

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
    call rios_p()

     do i = 1, nelec
        do j = 1, nelec
            if (i == j) cycle
            do nc = 1, nctype
                idx_phi = (iwctype(nc)-1)*(tmpc + 1) + offset_een
                idx_theta = idx_phi + (tmpc) * nctype

                cutoff = parm_bf(idx_phi+1)
                inv_cutoff = 1.0d0 / cutoff

                delta_ril(:) = x(:, i) - cent(:, nc)
                delta_rjl(:) = x(:, j) - cent(:, nc)
                delta(:) = x(:, i) - x(:, j)
                ril = sqrt(delta_ril(1)**2 + delta_ril(2)**2 + delta_ril(3)**2)
                rjl = sqrt(delta_rjl(1)**2 + delta_rjl(2)**2 + delta_rjl(3)**2)
                rij = sqrt(delta(1)**2 + delta(2)**2 + delta(3)**2)
                if (ril > cutoff .or. rjl > cutoff) cycle

                phi = 0.00d0
                phipi = 0.0d0
                phipj = 0.0d0
                phipp = 0.0d0
                theta = 0.00d0
                thetapi = 0.0d0
                thetapj = 0.0d0
                thetapp = 0.0d0

                inv_rij = 1.0d0 / rij
                inv_ril = 1.0d0 / ril
                inv_rjl = 1.0d0 / rjl

                fi = ((cutoff - ril)*inv_cutoff)**C
                fj = ((cutoff - rjl)*inv_cutoff)**C

                
                tmp1 = -C*inv_cutoff * ((cutoff - ril)*inv_cutoff)**(C-1) * inv_ril
                tmp2 = -C*inv_cutoff * ( - (C-1)*inv_cutoff ) * ((cutoff - ril)*inv_cutoff)**(C-2) * inv_ril * inv_ril
                do a = 1, 3
                    fpi(a) = tmp1 * delta_ril(a)
                    fppi(a) = tmp2 * delta_ril(a) * delta_ril(a) + &
                            tmp1 * ( 1 - delta_ril(a)*delta_ril(a)*inv_ril*inv_ril )
                end do
                tmp1 = -C*inv_cutoff * ((cutoff - rjl)*inv_cutoff)**(C-1) * inv_rjl
                tmp2 = -C*inv_cutoff * ( - (C-1)*inv_cutoff ) * ((cutoff - rjl)*inv_cutoff)**(C-2) * inv_rjl * inv_rjl
                do a = 1, 3
                    fpj(a) = tmp1 * delta_rjl(a)
                    fppj(a) = tmp2 * delta_rjl(a) * delta_rjl(a) + &
                            tmp1 * ( 1 - delta_rjl(a)*delta_rjl(a)*inv_rjl*inv_rjl )
                end do

                k = idx_phi+2
                kk = idx_theta+1
                do l = 0, nordc_bf
                    do m = 0, nordc_bf - l
                        do n = 0, nordc_bf - l - m
                            phi = phi + parm_bf(k) * ril**l * rjl**m * rij**n
                            theta = theta + parm_bf(kk) * ril**l * rjl**m * rij**n
                            do a = 1, 3
                                dquasi_dp(a,i,k) = dquasi_dp(a,i,k) + ril**l * rjl**m * rij**n * delta(a) * fi * fj
                                dquasi_dp(a,i,kk) = dquasi_dp(a,i,kk) + ril**l * rjl**m * rij**n * delta_ril(a) * fi * fj
                                if (l.gt.0) then
                                    phipi(a) = phipi(a) + parm_bf(k) * l * ril**(l-1) * rjl**m * rij**n * inv_ril * delta_ril(a)
                                    phipp(a) = phipp(a) + parm_bf(kk) * l * ril**(l-1) * rjl**m * rij**n * inv_ril * (1 + delta_ril(a)*delta_ril(a)*inv_ril*inv_ril)
                                    thetapi(a) = thetapi(a) + parm_bf(kk) * l * ril**(l-1) * rjl**m * rij**n * inv_ril * delta_ril(a)
                                    thetapp(a) = thetapp(a) + parm_bf(kk) * l * ril**(l-1) * rjl**m * rij**n * inv_ril * (1 + delta_ril(a)*delta_ril(a)*inv_ril*inv_ril)
                                end if
                                if (m.gt.0) then
                                    phipj(a) = phipj(a) + parm_bf(k) * m * ril**l * rjl**(m-1) * rij**n * inv_rjl * delta_rjl(a)
                                    thetapj(a) = thetapj(a) + parm_bf(kk) * m * ril**l * rjl**(m-1) * rij**n * inv_rjl * delta_rjl(a)
                                end if
                                if (l.gt.1) then
                                    phipp(a) = phipp(a) + parm_bf(k) * l * (l-1) * ril**(l-2) * rjl**m * rij**n * inv_ril * inv_ril * delta_ril(a) * delta_ril(a)
                                    thetapp(a) = thetapp(a) + parm_bf(kk) * l * (l-1) * ril**(l-2) * rjl**m * rij**n * inv_ril * inv_ril * delta_ril(a) * delta_ril(a)
                                end if
                                if (n.gt.0) then
                                    phipi(a) = phipi(a) + parm_bf(k) * n * ril**l * rjl**m * rij**(n-1) * inv_rij * delta(a)
                                    phipj(a) = phipj(a) - parm_bf(k) * n * ril**l * rjl**m * rij**(n-1) * inv_rij * delta(a)
                                    phipp(a) = phipp(a) + parm_bf(k) * n * ril**l * rjl**m * rij**(n-1) * inv_rij * (1 + delta(a)*delta(a)*inv_rij*inv_rij)
                                    thetapi(a) = thetapi(a) + parm_bf(kk) * n * ril**l * rjl**m * rij**(n-1) * inv_rij * delta(a)
                                    thetapj(a) = thetapj(a) - parm_bf(kk) * n * ril**l * rjl**m * rij**(n-1) * inv_rij * delta(a)
                                    thetapp(a) = thetapp(a) + parm_bf(kk) * n * ril**l * rjl**m * rij**(n-1) * inv_rij * (1 + delta(a)*delta(a)*inv_rij*inv_rij)
                                end if
                                if(n.gt.1) then
                                    phipp(a) = phipp(a) + parm_bf(k) * n * (n-1) * ril**l * rjl**m * rij**(n-2) * inv_rij * inv_rij * delta(a) * delta(a)
                                    thetapp(a) = thetapp(a) + parm_bf(kk) * n * (n-1) * ril**l * rjl**m * rij**(n-2) * inv_rij * inv_rij * delta(a) * delta(a)
                                end if
                                if (l.gt.0 .and. n.gt.0) then
                                    phipp(a) = phipp(a) + parm_bf(k) * 2 * l * n * ril**(l-1) * rjl**m * rij**(n-1) * inv_ril * inv_rij * delta_ril(a) * delta(a)
                                    thetapp(a) = thetapp(a) + parm_bf(kk) * 2 * l * n * ril**(l-1) * rjl**m * rij**(n-1) * inv_ril * inv_rij * delta_ril(a) * delta(a)
                                end if
                            enddo
                            k = k + 1
                            kk = kk + 1
                        end do
                    end do
                end do

 
                do a = 1, 3 
                    quasi_x(a,i) = quasi_x(a,i) + phi * delta(a) * fi * fj
                    quasi_x(a,i) = quasi_x(a,i) + theta * delta_ril(a) * fi * fj
                    dquasi_dp(a,i,idx_phi+1) = dquasi_dp(a,i,idx_phi+1) + &
                        delta(a) * phi * fj * C *  ((cutoff - ril)*inv_cutoff)**(C-1) * (ril*inv_cutoff*inv_cutoff) + &
                        delta(a) * phi * fi * C *  ((cutoff - rjl)*inv_cutoff)**(C-1) * (rjl*inv_cutoff*inv_cutoff) + &
                        delta_ril(a) * theta * fj * C *  ((cutoff - ril)*inv_cutoff)**(C-1) * (ril*inv_cutoff*inv_cutoff) + &
                        delta_ril(a) * theta * fi * C *  ((cutoff - rjl)*inv_cutoff)**(C-1) * (rjl*inv_cutoff*inv_cutoff) 

                    do b = 1, 3
                        dquasi_dx(a,i,b,i) = dquasi_dx(a,i,b,i) + &
                            phipi(b) * delta(a) * fi * fj + phi * delta(a) * fpi(b) * fj + &
                            thetapi(b) * delta_ril(a) * fi * fj + theta * delta_ril(a) * fpi(b) * fj

                        
                        dquasi_dx(a,i,b,j) = dquasi_dx(a,i,b,j) + &
                            phipj(b) * delta(a) * fi * fj + phi * delta(a) * fi * fpj(b) + &
                            thetapj(b) * delta_ril(a) * fi * fj + theta * delta_ril(a) * fi * fpj(b)


                        d2quasi_dx2(a,i,i) = d2quasi_dx2(a,i,i) + (&
                            phipp(b) * delta(a) * fi * fj + &
                            2.0d0 * phipi(b) * delta(a) * fpi(b) * fj + &
                            phi * delta(a) * fppi(b) * fj ) + (&
                            thetapp(b) * delta_ril(a) * fi * fj + &
                            2.0d0 * thetapi(b) * delta_ril(a) * fpi(b) * fj + &
                            theta * delta_ril(a) * fppi(b) * fj )

                        d2quasi_dx2(a,j,j) = d2quasi_dx2(a,j,j) + (&
                            phipp(b) * delta(a) * fi * fj + &
                            2.0d0 * phipj(b) * delta(a) * fi * fpj(b) + &
                            phi * delta(a) * fi * fppj(b) ) + (&
                            thetapp(b) * delta_ril(a) * fi * fj + &
                            2.0d0 * thetapj(b) * delta_ril(a) * fi * fpj(b) + &
                            theta * delta_ril(a) * fi * fppj(b) )
                    end do
                    dquasi_dx(a,i,a,i) = dquasi_dx(a,i,a,i) + phi * fi * fj + theta * fi * fj
                    dquasi_dx(a,i,a,j) = dquasi_dx(a,i,a,j) - phi * fi * fj

                    d2quasi_dx2(a,i,i) = d2quasi_dx2(a,i,i) + 2.0d0 * (phi * fpi(a) * fj + phipi(a) * fi * fj) + &
                                            2.0d0 * (theta * fpi(a) * fj + thetapi(a) * fi * fj)

                    d2quasi_dx2(a,j,j) = d2quasi_dx2(a,j,j) + 2.0d0 * (phi * fi * fpj(a) + phipj(a) * fi * fj)


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
    use m_backflow, only: quasi_x, dquasi_dx, d2quasi_dx2
    implicit none
    integer, intent(in) :: iel
    real(dp), dimension(3, nelec), intent(in) :: xold
    real(dp), dimension(3), intent(in) :: xnew
    real(dp), dimension(3, nelec), intent(out) :: quasi_x_new
    real(dp), dimension(3, nelec, 3, nelec), intent(out) :: dquasi_dx_new
    real(dp), dimension(3, nelec, nelec), intent(out) :: d2quasi_dx2_new
    integer, dimension(nelec), intent(out) :: indices
    real(dp) :: rij, rr, cutoff, b_one, a_one
    real(dp) :: f, fp(3), fpp(3), eta, etap(3), etapp(3)
    real(dp) :: delta(3)
    integer :: i, j, k, a, b, offset_ee, offset_en, offset_een, idx, C
    real(dp) :: inv_rij, inv_cutoff, tmp1, tmp2
    real(dp) :: cutoff1, cutoff2

    offset_ee = 0
    offset_en = offset_ee + 1+nordb_bf
    offset_een = offset_en + (1+norda_bf)*nctype

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

    b_one = parm_bf(offset_ee+2) * C * inv_cutoff


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
            tmp2 = -C*inv_cutoff * ( - (C-1)*inv_cutoff ) * cutoff2 * inv_rij * inv_rij
            do a = 1, 3
                fp(a) = tmp1 * delta(a)
                fpp(a) = tmp2 * delta(a) * delta(a) + &
                            tmp1 * ( 1 - delta(a)*delta(a)*inv_rij*inv_rij )
            end do


            eta = eta + parm_bf(offset_ee+2) + b_one*rij
            tmp1 = b_one * inv_rij
            tmp2 = b_one * inv_rij * inv_rij * inv_rij
            do a = 1, 3
                etap(a) = etap(a) + tmp1 * delta(a)
                etapp(a) = etapp(a) + tmp1 - tmp2 *  delta(a) * delta(a)
            end do


            rr = rij*rij
        
            tmp1 = inv_rij*inv_rij
            do k = 2, nordb_bf
                eta = eta + parm_bf(offset_ee+k+1)*rr
                do a = 1, 3
                    etap(a) = etap(a) + parm_bf(offset_ee+k+1) * delta(a) * tmp1 * k * rr
                    etapp(a) = etapp(a) + parm_bf(offset_ee+k+1) * tmp1 * rr * k * (1 + (k-2) * delta(a) * delta(a) * tmp1)
                enddo
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

                    tmp1 = etapp(b) * delta(a) * f + &
                            2.0d0 * etap(b) * delta(a) * fp(b) + &
                            eta * delta(a) * fpp(b)

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
            tmp2 = -C*inv_cutoff * ( - (C-1)*inv_cutoff ) * cutoff2 * inv_rij * inv_rij
            do a = 1, 3
                fp(a) = tmp1 * delta(a)
                fpp(a) = tmp2 * delta(a) * delta(a) + &
                            tmp1 * ( 1 - delta(a)*delta(a)*inv_rij*inv_rij )
            end do


            eta = eta + parm_bf(offset_ee+2) + b_one*rij
            tmp1 = b_one * inv_rij
            tmp2 = b_one * inv_rij * inv_rij * inv_rij
            do a = 1, 3
                etap(a) = etap(a) + tmp1 * delta(a)
                etapp(a) = etapp(a) + tmp1 - tmp2 *  delta(a) * delta(a)
            end do


            rr = rij*rij
        
            tmp1 = inv_rij*inv_rij
            do k = 2, nordb_bf
                eta = eta + parm_bf(offset_ee+k+1)*rr
                do a = 1, 3
                    etap(a) = etap(a) + parm_bf(offset_ee+k+1) * delta(a) * tmp1 * k * rr
                    etapp(a) = etapp(a) + parm_bf(offset_ee+k+1) * tmp1 * rr * k * (1 + (k-2) * delta(a) * delta(a) * tmp1)
                enddo
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

                    tmp1 = etapp(b) * delta(a) * f + &
                            2.0d0 * etap(b) * delta(a) * fp(b) + &
                            eta * delta(a) * fpp(b)

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

        a_one = parm_bf(offset_en+idx+2) * C * inv_cutoff


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
            tmp2 = -C*inv_cutoff * ( - (C-1)*inv_cutoff ) * cutoff2 * inv_rij * inv_rij
            do a = 1, 3
                fp(a) = tmp1 * delta(a)
                fpp(a) = tmp2 * delta(a) * delta(a) + &
                            tmp1 * ( 1 - delta(a)*delta(a)*inv_rij*inv_rij )
            end do

            eta  = eta + parm_bf(offset_en+idx+2) + a_one*rij

            tmp1 = a_one * inv_rij
            tmp2 = a_one * inv_rij * inv_rij * inv_rij
            do a = 1, 3
                etap(a) = etap(a) + tmp1 * delta(a)
                etapp(a) = etapp(a) + tmp1 - tmp2 *  delta(a) * delta(a)
            end do

            rr = rij*rij

            tmp1 = inv_rij*inv_rij
            do k = 2, norda_bf
                eta = eta + parm_bf(offset_en + idx+ k+1)*rr
                do a = 1, 3
                    etap(a) = etap(a) + parm_bf(offset_en + idx+ k+1) * k * rr * delta(a)*tmp1
                    etapp(a) = etapp(a) + parm_bf(offset_en + idx+ k+1) * k * tmp1 * rr * (1 + (k-2) * delta(a) * delta(a)*tmp1)
                enddo
                rr = rr * rij
            end do

            do a = 1, 3 
                quasi_x_new(a,iel) = quasi_x_new(a,iel) - eta * delta(a) * f
                do b = 1, 3
                    dquasi_dx_new(a,iel,b,iel) = dquasi_dx_new(a,iel,b,iel) - (&
                        etap(b) * delta(a) * f + eta * delta(a) * fp(b) )

                    d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) - (&
                        etapp(b) * delta(a) * f + &
                        2.0d0 * etap(b) * delta(a) * fp(b) + &
                        eta * delta(a) * fpp(b) )
                end do
                dquasi_dx_new(a,iel,a,iel) = dquasi_dx_new(a,iel,a,iel) - eta * f
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
            tmp2 = -C*inv_cutoff * ( - (C-1)*inv_cutoff ) * cutoff2 * inv_rij * inv_rij
            do a = 1, 3
                fp(a) = tmp1 * delta(a)
                fpp(a) = tmp2 * delta(a) * delta(a) + &
                            tmp1 * ( 1 - delta(a)*delta(a)*inv_rij*inv_rij )
            end do

            eta  = eta + parm_bf(offset_en+idx+2) + a_one*rij

            tmp1 = a_one * inv_rij
            tmp2 = a_one * inv_rij * inv_rij * inv_rij
            do a = 1, 3
                etap(a) = etap(a) + tmp1 * delta(a)
                etapp(a) = etapp(a) + tmp1 - tmp2 *  delta(a) * delta(a)
            end do

            rr = rij*rij

            tmp1 = inv_rij*inv_rij
            do k = 2, norda_bf
                eta = eta + parm_bf(offset_en + idx+ k+1)*rr
                do a = 1, 3
                    etap(a) = etap(a) + parm_bf(offset_en + idx+ k+1) * k * rr * delta(a)*tmp1
                    etapp(a) = etapp(a) + parm_bf(offset_en + idx+ k+1) * k * tmp1 * rr * (1 + (k-2) * delta(a) * delta(a)*tmp1)
                enddo
                rr = rr * rij
            end do

            do a = 1, 3 
                quasi_x_new(a,iel) = quasi_x_new(a,iel) + eta * delta(a) * f
                do b = 1, 3
                    dquasi_dx_new(a,iel,b,iel) = dquasi_dx_new(a,iel,b,iel) + (&
                        etap(b) * delta(a) * f + eta * delta(a) * fp(b) )

                    d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + (&
                        etapp(b) * delta(a) * f + &
                        2.0d0 * etap(b) * delta(a) * fp(b) + &
                        eta * delta(a) * fpp(b) )
                end do
                dquasi_dx_new(a,iel,a,iel) = dquasi_dx_new(a,iel,a,iel) + eta * f
                d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + 2 * (eta * fp(a) + etap(a) * f)
            end do
        endif

    end do


20  continue
    if (nordc_bf .eq. 0) return
end subroutine single_rios_backflow


end module