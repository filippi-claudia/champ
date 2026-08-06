module backflow_mod

contains

!> This function selects the backflow transformation to be used
!>
!> @details Backflow functions can be selected in the input file using the flag "backflow".
!> Flag 0: Blackflow disabled
!> Flag 1: Rios backflow
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

    if (ibackflow == 1) then
        call single_rios_backflow(iel, xold, xnew, quasi_x_new, dquasi_dx_new, d2quasi_dx2_new, indices)
    else
        call fatal_error('Single backflow only implemented for Rios backflow.')
    end if

    ! Only compute if not using QMCkl
    if (.not. use_qmckl_orbitals) then
        call bf_distances(quasi_x_new)
    end if
end subroutine single_backflow

subroutine init_backflow(iflag)
    use m_backflow, only: ibackflow
    implicit none
    integer, intent(in) :: iflag


    if (ibackflow == 0) return

    if (ibackflow == 1) then
        call init_rios_backflow(iflag,5,5,5)
    end if
end subroutine init_backflow

subroutine init_backflow_arrays()
    use m_backflow, only: ibackflow
    implicit none

    if (ibackflow == 0) return

    if (ibackflow == 1) then
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


subroutine init_rios_backflow(iflag, orda, ordb, ordc)
    use precision_kinds, only: dp
    use m_backflow, only: parm_bf, nparm_bf, norda_bf, nordb_bf, nspin_bf_ee, nordc_bf, cutoff_scale, ncparm_bf, cutoff_b_offset, cutoff_a_offset, cutoff_c_offset, ee_coeff_active, en_coeff_active
    use system, only: nctype
    use m_backflow, only: allocate_m_backflow
    implicit none
    integer, intent(in) :: iflag
    integer :: i, orda, ordb, ordc, l, m, n
    intrinsic :: ceiling

    if (iflag.eq.0) then
        norda_bf = orda
        nordb_bf = ordb
        nordc_bf = ordc
        if (nordb_bf.gt.0) then
            if (nspin_bf_ee.ne.1 .and. nspin_bf_ee.ne.2) nspin_bf_ee = 2
        else
            nspin_bf_ee = 1
        endif
    endif

    call init_backflow_arrays()

    if (iflag.eq.0) then
       
        parm_bf = 0.0d0
        if (ee_coeff_active.eq.1) parm_bf(cutoff_b_offset) = 3.0d0
        if (en_coeff_active.eq.1) then
            do i = 1, nctype
                parm_bf(cutoff_a_offset + i - 1) = 3.0d0
            end do
        end if
        if (nordc_bf.gt.0) then
            do i = 1, nctype
                parm_bf(cutoff_c_offset + i - 1) = 3.0d0
            end do
        end if
        cutoff_scale = 3

        call init_cusp()
        call fix_cusp()
    endif
end subroutine init_rios_backflow

subroutine init_rios_backflow_arrays()
    use m_backflow, only: allocate_m_backflow, ibackflow, norda_bf, nordb_bf, nspin_bf_ee, nordc_bf, nparm_bf, maxord, ncparm_bf, c_cuspconst, set_backflow_parameter_layout
    use system, only: nctype
    implicit none
    ! The c-section dimensions must be known before caching the complete layout.
    if (nordc_bf > 0) then
        ncparm_bf = (nordc_bf+1)*(nordc_bf+2)*(nordc_bf+3)/6
        c_cuspconst = 10*(nordc_bf + 1)
    else
        ncparm_bf = 0
        c_cuspconst = 0
    end if

    maxord = max(norda_bf, nordb_bf, nordc_bf)
    call set_backflow_parameter_layout()

    call allocate_m_backflow

end subroutine init_rios_backflow_arrays

subroutine init_cusp()
    use precision_kinds, only: dp
    use m_backflow, only: parm_bf, c_cuspconst, nspin_bf_ee, nordc_bf, ncparm_bf, cutoff_scale, B, dB_dcutoff, cusp_parameters, cusp_indices, inv_cusp_indices, inv_cusp_parameters, cusp_cutoff_deriv, basis_klmn, cutoff_c_offset, een_coeff_offset, een_component_block_size
    use system, only: nctype
    use control, only: ipr
    use contrl_file,    only: ounit
    implicit none
    integer :: k, l, m, n, alpha, idx, idx1, idx2, info, offset, idx_phi, idx_theta, i, j, linefound, eq_idx
    integer :: ee_spin, phi_channel, theta_channel
    integer :: pr, max_row, idx_pivot, cnt
    integer, dimension(c_cuspconst) :: ipiv
    real(dp) :: cutoff, pivot, max_val, factor, tmp, dtmp

    ! Initialize cusp dependency matrix B and its cutoff derivative to zero
    B = 0.0d0
    dB_dcutoff = 0.0d0
    cusp_cutoff_deriv = 0.0d0

    offset = een_coeff_offset

    basis_klmn = 0
    ! Store explicit absolute indices in the file order: phi(par), phi(anti), theta(par), theta(anti).
    do ee_spin = 1, nspin_bf_ee
        idx = 0
        do n = 1, nctype
            do k = 0, nordc_bf
                do l = 0, nordc_bf - k
                    do m = 0, nordc_bf - k - l
                        idx = idx + 1
                        basis_klmn(k,l,m,n,1,ee_spin) = offset + (ee_spin-1)*een_component_block_size + idx
                        basis_klmn(k,l,m,n,2,ee_spin) = offset + nspin_bf_ee*een_component_block_size + (ee_spin-1)*een_component_block_size + idx
                    end do
                end do
            end do
        end do
    end do

    cusp_parameters = 0.0d0
    cusp_indices = 0
    inv_cusp_parameters = 0.0d0
    inv_cusp_indices = 0

    eq_idx = 1

    do ee_spin = 1, nspin_bf_ee
        ! Define spin-specific cusp equations below: ee_spin=1 is parallel; ee_spin=2 is antiparallel.
        ! Branch on ee_spin around the B/dB_dcutoff updates here when their conditions differ.
        do n = 1, nctype
            phi_channel = n + (ee_spin - 1) * nctype
            theta_channel = nspin_bf_ee * nctype + n + (ee_spin - 1) * nctype
            idx_phi = offset + (ee_spin - 1) * een_component_block_size + ncparm_bf * (n - 1)
            idx_theta = offset + nspin_bf_ee * een_component_block_size + (ee_spin - 1) * een_component_block_size + ncparm_bf * (n - 1)
            cutoff = parm_bf(cutoff_c_offset + n - 1)

        do alpha=0,nordc_bf
            do k=0,nordc_bf
                do l=0,nordc_bf - k
                    do m=0,nordc_bf-k-l
                        ! if (k > l) then
                        ! else
                            idx1 = basis_klmn(k,l,m,n,1,ee_spin) - idx_phi
                            idx2 = basis_klmn(k,l,m,n,2,ee_spin) - idx_theta
                        ! endif

                        ! Parallel and antiparallel spin
                        if (k .eq. 0 .and. (l+m).eq.alpha) then
                            B(alpha+1, idx1, phi_channel) = B(alpha+1, idx1, phi_channel) - cutoff_scale/cutoff
                            ! d/d(cutoff) of -cutoff_scale/cutoff = cutoff_scale/cutoff^2
                            dB_dcutoff(alpha+1, idx1, phi_channel) = dB_dcutoff(alpha+1, idx1, phi_channel) + cutoff_scale/(cutoff*cutoff)
                        endif
                        ! Parallel and antiparallel spin
                        if (k .eq. 1 .and. (l+m).eq.alpha .and. (k+l+m).le.nordc_bf) then
                            B(alpha+1, idx1, phi_channel) = B(alpha+1, idx1, phi_channel) + 1.0d0
                            ! Constant term has zero derivative
                        endif
                        ! Parallel and antiparallel spin
                        if (l .eq. 0 .and. (k+m).eq.alpha) then
                            B(alpha+nordc_bf+1+1, idx1, phi_channel) = B(alpha+nordc_bf+1+1, idx1, phi_channel) - cutoff_scale/cutoff
                            B(alpha+1, idx2, theta_channel) = B(alpha+1, idx2, theta_channel) - cutoff_scale/cutoff
                            dB_dcutoff(alpha+nordc_bf+1+1, idx1, phi_channel) = dB_dcutoff(alpha+nordc_bf+1+1, idx1, phi_channel) + cutoff_scale/(cutoff*cutoff)
                            dB_dcutoff(alpha+1, idx2, theta_channel) = dB_dcutoff(alpha+1, idx2, theta_channel) + cutoff_scale/(cutoff*cutoff)
                        endif
                        ! Parallel and antiparallel spin
                        if (l .eq. 1 .and. (k+m).eq.alpha .and. (k+l+m).le.nordc_bf) then
                            B(alpha+nordc_bf+1+1, idx1, phi_channel) = B(alpha+nordc_bf+1+1, idx1, phi_channel) + 1.0d0
                            B(alpha+1, idx2, theta_channel) = B(alpha+1, idx2, theta_channel) + 1.0d0
                        endif
                        ! Only the phi contribution is parallel-spin-specific.
                        if (m .eq. 1 .and. (k+l).eq.alpha .and. (k+l+m).le.nordc_bf .and. ee_spin .eq. 1) then
                            B(alpha+(nordc_bf+1)*2+1, idx1, phi_channel) = B(alpha+(nordc_bf+1)*2+1, idx1, phi_channel) + 1.0d0
                        endif
                        ! The theta contribution applies to parallel and antiparallel spin.
                        if (m .eq. 1 .and. (k+l).eq.alpha .and. (k+l+m).le.nordc_bf) then
                            B(alpha+(nordc_bf+1)  +1, idx2, theta_channel) = B(alpha+(nordc_bf+1)  +1, idx2, theta_channel) + 1.0d0
                        endif

                        ! All these only for parallel spin
                        if (m .eq. 1 .and. (k+l).eq.alpha .and. (k+l+m).le.nordc_bf .and. ee_spin .eq. 1) then
                            B(alpha+(nordc_bf+1)*2+1, idx2, theta_channel) = B(alpha+(nordc_bf+1)*2+1, idx2, theta_channel) - cutoff_scale/cutoff
                            dB_dcutoff(alpha+(nordc_bf+1)*2+1, idx2, theta_channel) = dB_dcutoff(alpha+(nordc_bf+1)*2+1, idx2, theta_channel) + cutoff_scale/(cutoff*cutoff)
                            if (k .gt. 0) then
                                B(alpha+(nordc_bf+1)*2+1, idx2, theta_channel) = B(alpha+(nordc_bf+1)*2+1, idx2, theta_channel) - k/cutoff
                                dB_dcutoff(alpha+(nordc_bf+1)*2+1, idx2, theta_channel) = dB_dcutoff(alpha+(nordc_bf+1)*2+1, idx2, theta_channel) + k/(cutoff*cutoff)
                                B(alpha+(nordc_bf+1)*2  , idx2, theta_channel) = B(alpha+(nordc_bf+1)*2  , idx2, theta_channel) + k
                            end if

                            B(alpha+(nordc_bf+1)*3+1, idx2, theta_channel) = B(alpha+(nordc_bf+1)*3+1, idx2, theta_channel) - cutoff_scale/cutoff
                            dB_dcutoff(alpha+(nordc_bf+1)*3+1, idx2, theta_channel) = dB_dcutoff(alpha+(nordc_bf+1)*3+1, idx2, theta_channel) + cutoff_scale/(cutoff*cutoff)
                            if (l .gt. 0) then
                                B(alpha+(nordc_bf+1)*3+1, idx2, theta_channel) = B(alpha+(nordc_bf+1)*3+1, idx2, theta_channel) - l/cutoff
                                dB_dcutoff(alpha+(nordc_bf+1)*3+1, idx2, theta_channel) = dB_dcutoff(alpha+(nordc_bf+1)*3+1, idx2, theta_channel) + l/(cutoff*cutoff)
                                B(alpha+(nordc_bf+1)*3  , idx2, theta_channel) = B(alpha+(nordc_bf+1)*3  , idx2, theta_channel) + l
                            end if
                        end if
                    enddo
                enddo
            enddo
        enddo

        ! Solve Phi terms (B(:,:,phi_channel)) using Gaussian Elimination to RREF
        ! Apply same operations to dB_dcutoff to get derivative of RREF result
        pr = 0
        do j = 1, ncparm_bf ! Loop over variables (columns)
             if (pr >= c_cuspconst) exit
            
             ! -- Partial Pivoting --
             ! Find the row with the largest absolute value in the current column j
             ! starting from the current pivot row pr + 1
             max_val = 0.0d0
             max_row = -1
             do i = pr + 1, c_cuspconst
                 if (abs(B(i,j,phi_channel)) > max_val) then
                     max_val = abs(B(i,j,phi_channel))
                     max_row = i
                 end if
             end do
            
             ! If column is zero (or near zero), skip it -> this variable is a free variable
             if (max_val < 1.0d-12) cycle
            
             pr = pr + 1 ! We found a pivot for this column, move to next row
            
             ! Swap rows to bring the pivot to position (pr, j)
             ! Apply same swap to dB_dcutoff
             if (max_row /= pr) then
                 do k = 1, ncparm_bf
                     tmp = B(pr, k, phi_channel)
                     B(pr, k, phi_channel) = B(max_row, k, phi_channel)
                     B(max_row, k, phi_channel) = tmp
                     dtmp = dB_dcutoff(pr, k, phi_channel)
                     dB_dcutoff(pr, k, phi_channel) = dB_dcutoff(max_row, k, phi_channel)
                     dB_dcutoff(max_row, k, phi_channel) = dtmp
                 end do
             end if
            
             ! Normalize the pivot row so the leading coefficient is 1.0
             ! For derivative: d(B/pivot)/dL = dB/dL / pivot - B * dpivot/dL / pivot^2
             !                               = (dB/dL * pivot - B * dpivot/dL) / pivot^2
             pivot = B(pr, j, phi_channel)
             dtmp = dB_dcutoff(pr, j, phi_channel)  ! Save dpivot/dL before loop modifies it
             do k = j, ncparm_bf
                 ! First update derivative before updating B
                 dB_dcutoff(pr, k, phi_channel) = (dB_dcutoff(pr, k, phi_channel) * pivot - B(pr, k, phi_channel) * dtmp) / (pivot * pivot)
                 B(pr, k, phi_channel) = B(pr, k, phi_channel) / pivot
             end do
            
             ! Eliminate entries in this column for all other rows (above and below)
             ! This converts the matrix to Reduced Row Echelon Form (RREF)
             ! For derivative: d(B - factor*B_pr)/dL = dB/dL - dfactor/dL * B_pr - factor * dB_pr/dL
             do i = 1, c_cuspconst
                 if (i == pr) cycle
                 factor = B(i, j, phi_channel)
                 dtmp = dB_dcutoff(i, j, phi_channel)  ! Save dfactor/dL before loop modifies it
                 if (abs(factor) < 1.0d-12) then
                     ! Still apply derivative elimination even if factor is small
                     ! because dfactor/dL might not be small
                     do k = j, ncparm_bf
                         dB_dcutoff(i, k, phi_channel) = dB_dcutoff(i, k, phi_channel) - dtmp * B(pr, k, phi_channel) &
                                              - factor * dB_dcutoff(pr, k, phi_channel)
                     end do
                     cycle
                 end if
                 do k = j, ncparm_bf
                     dB_dcutoff(i, k, phi_channel) = dB_dcutoff(i, k, phi_channel) - dtmp * B(pr, k, phi_channel) &
                                          - factor * dB_dcutoff(pr, k, phi_channel)
                     B(i, k, phi_channel) = B(i, k, phi_channel) - factor * B(pr, k, phi_channel)
                 end do
             end do
        end do
       
        ! Extract dependencies from RREF matrix
        do i = 1, pr ! Loop over the pivot rows
             ! Find the pivot column for this row (first non-zero entry)
             idx_pivot = -1
             do j = 1, ncparm_bf
                 if (abs(B(i, j, phi_channel)) > 1.0d-12) then
                     idx_pivot = j
                     exit
                 end if
             end do
            
             if (idx_pivot == -1) cycle ! Should not happen if pr logic is correct
            
             ! The variable corresponding to idx_pivot is a DEPENDENT variable
             cusp_indices(eq_idx, 1) = idx_phi + idx_pivot
             cusp_parameters(eq_idx, 1) = 1.0d0
             ! Cutoff derivative of dependent parameter coefficient is from dB_dcutoff
             cusp_cutoff_deriv(eq_idx, 1) = 0.0d0  ! Coefficient of dependent var is always 1, so derivative is 0
            
             ! All other non-zero entries in this row correspond to FREE variables
             ! (or their linear combination) that the dependent variable depends on.
             ! Because of RREF, there is only one pivot per row, so this relationship is unique.
             cnt = 1
             do k = idx_pivot + 1, ncparm_bf
                 if (abs(B(i, k, phi_channel)) > 1.0d-12 .or. abs(dB_dcutoff(i, k, phi_channel)) > 1.0d-12) then
                     cnt = cnt + 1
                     cusp_indices(eq_idx, cnt) = idx_phi + k
                     ! Move terms to RHS: x_dep + c * x_free = 0  =>  x_dep = -c * x_free
                     cusp_parameters(eq_idx, cnt) = -B(i, k, phi_channel)
                     ! Derivative: d(-B)/dL = -dB/dL
                     cusp_cutoff_deriv(eq_idx, cnt) = -dB_dcutoff(i, k, phi_channel)
                    
                     ! Store inverse mapping for derivatives:
                     ! When optimizing x_free, we must also update the derivative wrt x_dep
                     inv_cusp_indices(idx_phi + k, 1) = idx_phi + k
                     do l = 2, ncparm_bf
                          if (inv_cusp_indices(idx_phi + k, l) .eq. 0) then
                              inv_cusp_indices(idx_phi + k, l) = cusp_indices(eq_idx, 1)
                              inv_cusp_parameters(idx_phi + k, l) = cusp_parameters(eq_idx, cnt)
                              exit
                          endif
                     end do
                 end if
             end do
             eq_idx = eq_idx + 1
        end do

        ! Solve Theta terms (B(:,:,theta_channel)) using Gaussian Elimination to RREF
        ! Apply same operations to dB_dcutoff to get derivative of RREF result
        pr = 0
        do j = 1, ncparm_bf ! Loop over variables (columns)
             if (pr >= c_cuspconst) exit
            
             ! -- Partial Pivoting --
             ! Find the row with the largest absolute value in the current column j
             ! starting from the current pivot row pr + 1
             max_val = 0.0d0
             max_row = -1
             do i = pr + 1, c_cuspconst
                 if (abs(B(i,j,theta_channel)) > max_val) then
                     max_val = abs(B(i,j,theta_channel))
                     max_row = i
                 end if
             end do
            
             ! If column is zero (or near zero), skip it -> this variable is a free variable
             if (max_val < 1.0d-12) cycle
            
             pr = pr + 1 ! We found a pivot for this column, move to next row
            
             ! Swap rows to bring the pivot to position (pr, j)
             ! Apply same swap to dB_dcutoff
             if (max_row /= pr) then
                 do k = 1, ncparm_bf
                     tmp = B(pr, k, theta_channel)
                     B(pr, k, theta_channel) = B(max_row, k, theta_channel)
                     B(max_row, k, theta_channel) = tmp
                     dtmp = dB_dcutoff(pr, k, theta_channel)
                     dB_dcutoff(pr, k, theta_channel) = dB_dcutoff(max_row, k, theta_channel)
                     dB_dcutoff(max_row, k, theta_channel) = dtmp
                 end do
             end if
            
             ! Normalize the pivot row so the leading coefficient is 1.0
             ! For derivative: d(B/pivot)/dL = dB/dL / pivot - B * dpivot/dL / pivot^2
             !                               = (dB/dL * pivot - B * dpivot/dL) / pivot^2
             pivot = B(pr, j, theta_channel)
             dtmp = dB_dcutoff(pr, j, theta_channel)  ! Save dpivot/dL before loop modifies it
             do k = j, ncparm_bf
                 ! First update derivative before updating B
                 dB_dcutoff(pr, k, theta_channel) = (dB_dcutoff(pr, k, theta_channel) * pivot - B(pr, k, theta_channel) * dtmp) / (pivot * pivot)
                 B(pr, k, theta_channel) = B(pr, k, theta_channel) / pivot
             end do
            
             ! Eliminate entries in this column for all other rows (above and below)
             ! This converts the matrix to Reduced Row Echelon Form (RREF)
             ! For derivative: d(B - factor*B_pr)/dL = dB/dL - dfactor/dL * B_pr - factor * dB_pr/dL
             do i = 1, c_cuspconst
                 if (i == pr) cycle
                 factor = B(i, j, theta_channel)
                 dtmp = dB_dcutoff(i, j, theta_channel)  ! Save dfactor/dL before loop modifies it
                 if (abs(factor) < 1.0d-12) then
                     ! Still apply derivative elimination even if factor is small
                     ! because dfactor/dL might not be small
                     do k = j, ncparm_bf
                         dB_dcutoff(i, k, theta_channel) = dB_dcutoff(i, k, theta_channel) - dtmp * B(pr, k, theta_channel) &
                                              - factor * dB_dcutoff(pr, k, theta_channel)
                     end do
                     cycle
                 end if
                 do k = j, ncparm_bf
                     dB_dcutoff(i, k, theta_channel) = dB_dcutoff(i, k, theta_channel) - dtmp * B(pr, k, theta_channel) &
                                          - factor * dB_dcutoff(pr, k, theta_channel)
                     B(i, k, theta_channel) = B(i, k, theta_channel) - factor * B(pr, k, theta_channel)
                 end do
             end do
        end do
       
        ! Extract dependencies from RREF matrix
        do i = 1, pr ! Loop over the pivot rows
             ! Find the pivot column for this row (first non-zero entry)
             idx_pivot = -1
             do j = 1, ncparm_bf
                 if (abs(B(i, j, theta_channel)) > 1.0d-12) then
                     idx_pivot = j
                     exit
                 end if
             end do
            
             if (idx_pivot == -1) cycle ! Should not happen if pr logic is correct
            
             ! The variable corresponding to idx_pivot is a DEPENDENT variable
             cusp_indices(eq_idx, 1) = idx_theta + idx_pivot
             cusp_parameters(eq_idx, 1) = 1.0d0
             cusp_cutoff_deriv(eq_idx, 1) = 0.0d0  ! Coefficient of dependent var is always 1, so derivative is 0
            
             ! All other non-zero entries in this row correspond to FREE variables
             ! (or their linear combination) that the dependent variable depends on.
             ! Because of RREF, there is only one pivot per row, so this relationship is unique.
             cnt = 1
             do k = idx_pivot + 1, ncparm_bf
                 if (abs(B(i, k, theta_channel)) > 1.0d-12 .or. abs(dB_dcutoff(i, k, theta_channel)) > 1.0d-12) then
                     cnt = cnt + 1
                     cusp_indices(eq_idx, cnt) = idx_theta + k
                     ! Move terms to RHS: x_dep + c * x_free = 0  =>  x_dep = -c * x_free
                     cusp_parameters(eq_idx, cnt) = -B(i, k, theta_channel)
                     ! Derivative: d(-B)/dL = -dB/dL
                     cusp_cutoff_deriv(eq_idx, cnt) = -dB_dcutoff(i, k, theta_channel)
                    
                     ! Store inverse mapping for derivatives:
                     ! When optimizing x_free, we must also update the derivative wrt x_dep
                     inv_cusp_indices(idx_theta + k, 1) = idx_theta + k
                     do l = 2, ncparm_bf
                          if (inv_cusp_indices(idx_theta + k, l) .eq. 0) then
                              inv_cusp_indices(idx_theta + k, l) = cusp_indices(eq_idx, 1)
                              inv_cusp_parameters(idx_theta + k, l) = cusp_parameters(eq_idx, cnt)
                              exit
                          endif
                     end do
                 end if
             end do
             eq_idx = eq_idx + 1
        end do

        end do
    end do

    if (ipr.gt.2) then
        write(ounit, *) "Cusp conditions for Rios backflow:"
        idx = 0
        do i=1, c_cuspconst*nctype*nspin_bf_ee
            if (cusp_indices(i,1) .eq. 0) cycle
            idx = idx + 1

            do j = 1, ncparm_bf
                if (cusp_indices(i,j) .eq. 0) exit
                write(ounit, '(I8, F20.12)') cusp_indices(i,j), cusp_parameters(i,j)
            end do
            write(ounit, *) "---------------------"
        end do

        write(ounit, *) "Total number of backflow cusp conditions:", idx
    endif

end subroutine init_cusp

subroutine fix_cusp()
    use precision_kinds, only: dp
    use system, only: nctype
    use m_backflow, only: parm_bf, ncparm_bf, nspin_bf_ee, cusp_parameters, cusp_indices, c_cuspconst
    implicit none

    integer :: k, kk

    ! Each row already contains absolute indices for one phi/theta/spin channel.
    do k = 1, c_cuspconst*nctype*nspin_bf_ee
        if (cusp_indices(k,1) .gt. 0) then
            parm_bf(cusp_indices(k,1)) = 0.0d0
            do kk = 2, ncparm_bf
                if (cusp_indices(k,kk) .eq. 0) exit
                parm_bf(cusp_indices(k,1)) = parm_bf(cusp_indices(k,1)) + &
                    cusp_parameters(k,kk) * parm_bf(cusp_indices(k,kk))
            end do
        end if
    end do

end subroutine fix_cusp


subroutine rios_distances(x)
    use precision_kinds, only: dp
    use system, only: nelec, ncent, cent, iwctype, nctype
    use m_backflow, only: parm_bf, nparm_bf, norda_bf, nordb_bf, nspin_bf_ee, nordc_bf, cutoff_scale, ncparm_bf, r_en, rvec_en, r_ee, rvec_ee, r_ee_gl, r_en_gl, maxord, cutoff_deriv, cutoff_a_offset, cutoff_c_offset, en_coeff_active, een_coeff_active
    use m_backflow, only: r_en, rvec_en, r_ee, rvec_ee, r_ee_gl, r_en_gl, maxord, cutoff_deriv
    implicit none
    real(dp), dimension(3, nelec), intent(in) :: x
    real(dp) :: r, inv_r, r_cutoff, inv_r_cutoff, cutoff, inv_cutoff
    integer :: i, j, nc, k, no, l, m, n, cc

   
    r_en = 0.0d0
    r_ee = 0.0d0
    rvec_en = 0.0d0
    rvec_ee = 0.0d0
    r_ee_gl = 0.0d0
    r_en_gl = 0.0d0
    cutoff_deriv = 0.0d0

    do cc = 1, 2
        do nc = 1, ncent
            if (cc == 1) then
                if (en_coeff_active .eq. 0) cycle
                r_cutoff = parm_bf(cutoff_a_offset + iwctype(nc) - 1)
            else
                if (een_coeff_active .eq. 0) cycle
                r_cutoff = parm_bf(cutoff_c_offset + iwctype(nc) - 1)
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
    use m_backflow, only: parm_bf, nparm_bf, norda_bf, nordb_bf, nspin_bf_ee, nordc_bf, cutoff_scale, ncparm_bf, r_en, rvec_en, r_ee, rvec_ee, r_ee_gl, r_en_gl, maxord, single_r_en, single_rvec_en, single_r_ee, single_rvec_ee, single_r_ee_gl, single_r_en_gl, cutoff_a_offset, cutoff_c_offset, en_coeff_active, een_coeff_active
    use m_backflow, only: r_en, rvec_en, r_ee, rvec_ee, r_ee_gl, r_en_gl, maxord
    use m_backflow, only: single_r_en, single_rvec_en, single_r_ee, single_rvec_ee, single_r_ee_gl, single_r_en_gl
    implicit none
    real(dp), dimension(3, nelec), intent(in) :: x
    real(dp), dimension(3), intent(in) :: xnew
    integer, intent(in) :: iel
    real(dp) :: r, inv_r, r_cutoff, inv_r_cutoff, cutoff, inv_cutoff
    integer :: i, j, nc, k, no, l, m, n, cc

   
    do cc = 1, 2
        do nc = 1, ncent
            if (cc == 1) then
                if (en_coeff_active .eq. 0) cycle
                r_cutoff = parm_bf(cutoff_a_offset + iwctype(nc) - 1)
            else
                if (een_coeff_active .eq. 0) cycle
                r_cutoff = parm_bf(cutoff_c_offset + iwctype(nc) - 1)
            end if

            inv_r_cutoff = 1.0d0 / r_cutoff 

   
            do k = 1, 3
                single_rvec_en(k, nc) = xnew(k) - cent(k, nc)
            end do
            r = sqrt( single_rvec_en(1,nc)**2 + single_rvec_en(2,nc)**2 + single_rvec_en(3,nc)**2 )

            inv_r = 1.0d0 / r
           
            cutoff= ((r_cutoff - r) * inv_r_cutoff)**cutoff_scale
            inv_cutoff = 1/((r_cutoff - r) * inv_r_cutoff)

            single_r_en(nc, 0, cc) = cutoff
            if (r > r_cutoff) cycle


            single_r_en_gl(1, nc, 0, cc) = -cutoff_scale * inv_r_cutoff * inv_r * single_rvec_en(1, nc) * cutoff * inv_cutoff
            single_r_en_gl(2, nc, 0, cc) = -cutoff_scale * inv_r_cutoff * inv_r * single_rvec_en(2, nc) * cutoff * inv_cutoff
            single_r_en_gl(3, nc, 0, cc) = -cutoff_scale * inv_r_cutoff * inv_r * single_rvec_en(3, nc) * cutoff * inv_cutoff
            ! single_r_en_gl(4, nc, 0, cc) = -cutoff_scale * inv_r_cutoff * 2.0d0 * inv_r * cutoff * inv_cutoff + &
            !                             cutoff_scale * (cutoff_scale - 1) * inv_r_cutoff * inv_r_cutoff * inv_cutoff* inv_cutoff * cutoff

            do no = 1, maxord
                single_r_en(nc, no, cc) = single_r_en(nc, no-1, cc) * r
                single_r_en_gl(1, nc, no, cc) = no * single_r_en(nc, no-1, cc) * inv_r * single_rvec_en(1, nc) - &
                                            cutoff_scale * inv_r_cutoff * inv_r * single_rvec_en(1, nc) * inv_cutoff * single_r_en(nc, no, cc)
                single_r_en_gl(2, nc, no, cc) = no * single_r_en(nc, no-1, cc) * inv_r * single_rvec_en(2, nc) - &
                                            cutoff_scale * inv_r_cutoff * inv_r * single_rvec_en(2, nc) * inv_cutoff * single_r_en(nc, no, cc)
                single_r_en_gl(3, nc, no, cc) = no * single_r_en(nc, no-1, cc) * inv_r * single_rvec_en(3, nc) - &
                                            cutoff_scale * inv_r_cutoff * inv_r * single_rvec_en(3, nc) * inv_cutoff * single_r_en(nc, no, cc)
                ! single_r_en_gl(4, nc, no, cc) = - cutoff_scale * inv_r_cutoff * (no+1) * 2.0d0 * single_r_en(nc, no-1, cc) * inv_cutoff +&
                !                             no * (no+1) * single_r_en(nc, no-1, cc) * inv_r + &
                !                             cutoff_scale * (cutoff_scale - 1) * inv_r_cutoff * inv_r_cutoff * inv_cutoff * inv_cutoff * single_r_en(nc, no, cc)
            end do
        end do
    end do

    do j = 1, nelec
        if (iel == j) cycle
        do k = 1, 3
            single_rvec_ee(k, j) = xnew(k) - x(k, j)
        end do
        r = sqrt( single_rvec_ee(1,j)**2 + single_rvec_ee(2,j)**2 + single_rvec_ee(3,j)**2 )
        inv_r = 1.0d0 / r

        single_r_ee(j, 0) = 1.0d0
        single_r_ee_gl(1, j, 0) = 0.0d0
        single_r_ee_gl(2, j, 0) = 0.0d0
        single_r_ee_gl(3, j, 0) = 0.0d0
        ! single_r_ee_gl(4, j, 0) = 0.0d0


        do no = 1, maxord
            single_r_ee(j, no) = single_r_ee(j, no-1) * r
            single_r_ee_gl(1, j, no) = no * single_r_ee(j, no-1) * inv_r * single_rvec_ee(1, j)
            single_r_ee_gl(2, j, no) = no * single_r_ee(j, no-1) * inv_r * single_rvec_ee(2, j)
            single_r_ee_gl(3, j, no) = no * single_r_ee(j, no-1) * inv_r * single_rvec_ee(3, j)
            ! single_r_ee_gl(4, j, no) = no * single_r_ee(j, no-1) * 2.0d0 * inv_r
        end do
        ! do no = 2, maxord
        !     single_r_ee_gl(4, j, no) = single_r_ee_gl(4, j, no) + no * (no - 1) * single_r_ee(j, no-2)
        ! end do
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
    use system, only: nelec, nup, iwctype, ncent, nctype, cent
    use optwf_control, only: ioptci, ioptjas, ioptorb, ioptbf
    use m_backflow, only: parm_bf, nparm_bf, norda_bf, nordb_bf, nspin_bf_ee, nordc_bf, cutoff_scale
    use m_backflow, only: r_en, rvec_en, r_ee, rvec_ee, r_ee_gl, r_en_gl, p, d_p, cutoff_deriv
    use m_backflow, only: inv_cusp_indices, inv_cusp_parameters, cusp_indices, c_cuspconst, ncparm_bf
    use m_backflow, only: cusp_cutoff_deriv, cusp_parameters, basis_klmn, cutoff_b_offset, cutoff_a_offset, cutoff_c_offset, ee_coeff_offset, en_coeff_offset, een_coeff_offset, ee_coeff_block_size, en_coeff_block_size, een_component_block_size
    implicit none
    real(dp), dimension(3, nelec), intent(in) :: x
    real(dp), dimension(3, nelec), intent(out) :: quasi_x
    real(dp), dimension(3, nelec, 3, nelec), intent(out) :: dquasi_dx
    real(dp), dimension(3, nelec, nelec), intent(out) :: d2quasi_dx2
    real(dp), dimension(3, nelec, nparm_bf), intent(out) :: dquasi_dp
    real(dp) :: rij, rr, cutoff, b_one, a_one
    real(dp) :: f, fp(3), fpp, eta, etap(3), etapp
    real(dp) :: delta(3), delta_ril(3), delta_rjl(3)
    integer :: i, j, k, a, b, idx, C
    integer :: ee_spin, offset_ee_ch, offset_een_ch
    integer :: nc, l, m, n, idx_phi, idx_theta, kk
    real(dp) :: ril, rjl, inv_ril, inv_rjl
    real(dp) :: inv_rij, inv_cutoff, tmp1, tmp2, cutoff1, cutoff2
    real(dp) :: fi, fpi(3), fppi(3), fj, fpj(3), fppj(3)
    real(dp) :: dp_dep_dcutoff
    real(dp) :: phi, theta, phipi(3), thetapi(3), phipj(3), thetapj(3), phippi, thetappi, phippj, thetappj

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

    do i = 1, nelec
        do j = 1, nelec
            if (i == j) cycle

            if (nspin_bf_ee.eq.1) then
                ee_spin = 1
            else
                if ((i.le.nup .and. j.le.nup) .or. (i.gt.nup .and. j.gt.nup)) then
                    ee_spin = 1
                else
                    ee_spin = 2
                endif
            endif
            offset_ee_ch = ee_coeff_offset + (ee_spin-1)*ee_coeff_block_size

            cutoff = parm_bf(cutoff_b_offset)
            inv_cutoff = 1.0d0 / cutoff
            if (ee_spin .eq. 1) then
                b_one = parm_bf(offset_ee_ch+2) * cutoff/C
                parm_bf(offset_ee_ch+1) = b_one
            else
                b_one = parm_bf(offset_ee_ch+1)
            endif



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
            if (ee_spin .eq. 1) then
                do a = 1, 3
                    dquasi_dp(a,i,offset_ee_ch+2) = dquasi_dp(a,i,offset_ee_ch+2) + delta(a) * f * cutoff / C
                end do
            else
                do a = 1, 3
                    dquasi_dp(a,i,offset_ee_ch+1) = dquasi_dp(a,i,offset_ee_ch+1) + f * delta(a)
                end do
            end if

            rr = rij
           
            tmp2 = inv_rij*inv_rij
            do k = 1, nordb_bf
                eta = eta + parm_bf(offset_ee_ch+k+1)*rr
                do a = 1, 3
                    dquasi_dp(a,i,offset_ee_ch+k+1) = dquasi_dp(a,i,offset_ee_ch+k+1) + rr * delta(a) * f
                    etap(a) = etap(a) + parm_bf(offset_ee_ch+k+1) * delta(a) * tmp2 * k * rr
                enddo
                etapp = etapp + parm_bf(offset_ee_ch+k+1) * k * rr * tmp2 * (k+1)
                rr = rr * rij
            end do

            tmp2 = eta * C * cutoff1 * (rij*inv_cutoff*inv_cutoff)
            do a = 1, 3
                quasi_x(a,i) = quasi_x(a,i) + eta * delta(a) * f
                dquasi_dp(a,i,cutoff_b_offset) = dquasi_dp(a,i,cutoff_b_offset) + delta(a) * tmp2
                if (ee_spin .eq. 1) then
                    dquasi_dp(a,i,cutoff_b_offset) = dquasi_dp(a,i,cutoff_b_offset) + delta(a) * f * b_one * inv_cutoff
                endif
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

10  if (norda_bf .eq. 0) goto 20

    do j = 1, ncent
        idx = (iwctype(j)-1)*(norda_bf+1)
        cutoff = parm_bf(cutoff_a_offset + iwctype(j) - 1)
        inv_cutoff = 1.0d0 / cutoff

        a_one = parm_bf(en_coeff_offset+idx+1)

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
                dquasi_dp(a,i,en_coeff_offset + idx+1) = dquasi_dp(a,i,en_coeff_offset + idx+1) + delta(a) * f
            end do

            rr = rij

            tmp1 = inv_rij*inv_rij
            do k = 1, norda_bf
                eta = eta + parm_bf(en_coeff_offset + idx+ k+1)*rr
                do a = 1, 3
                    dquasi_dp(a,i,en_coeff_offset + idx+ k+1) = dquasi_dp(a,i,en_coeff_offset + idx+ k+1) + rr * delta(a) * f
                    etap(a) = etap(a) + parm_bf(en_coeff_offset + idx+ k+1) * k * rr * delta(a)*tmp1
                enddo
                etapp = etapp + parm_bf(en_coeff_offset + idx+ k+1) * k * rr * tmp1 * (k+1)
                rr = rr * rij
            end do


            tmp1 = eta * C * cutoff1 * (rij*inv_cutoff*inv_cutoff)
            do a = 1, 3
                quasi_x(a,i) = quasi_x(a,i) + eta * delta(a) * f
                dquasi_dp(a,i,cutoff_a_offset + iwctype(j) - 1) = dquasi_dp(a,i,cutoff_a_offset + iwctype(j) - 1) + &
                    delta(a) * tmp1
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

20  if (nordc_bf .eq. 0) return

    call rios_distances(x)

     do i = 1, nelec
        do j = 1, nelec
            if (i == j) cycle
            if (nspin_bf_ee.eq.1) then
                ee_spin = 1
            else
                if ((i.le.nup .and. j.le.nup) .or. (i.gt.nup .and. j.gt.nup)) then
                    ee_spin = 1
                else
                    ee_spin = 2
                endif
            endif
            offset_een_ch = een_coeff_offset + (ee_spin-1)*een_component_block_size
            do nc = 1, ncent
                idx_phi = (iwctype(nc)-1)*ncparm_bf + offset_een_ch
                idx_theta = nspin_bf_ee*een_component_block_size + (iwctype(nc)-1)*ncparm_bf + offset_een_ch

                cutoff = parm_bf(cutoff_c_offset + iwctype(nc) - 1)
                inv_cutoff = 1.0d0 / cutoff

                ! TODO change all these cutoff checks, as they only work for odd C values.
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

                do l = 0, nordc_bf
                    do m = 0, nordc_bf - l
                        do n = 0, nordc_bf - l - m
                            ! if (l > m) then
                            ! else
                                k = basis_klmn(l,m,n,iwctype(nc),1,ee_spin)
                                kk = basis_klmn(l,m,n,iwctype(nc),2,ee_spin)
                            ! end if
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
                    idx = cutoff_c_offset + iwctype(nc) - 1
                    dquasi_dp(a,i,idx) = dquasi_dp(a,i,idx) + phi * rvec_ee(a,i,j) * (cutoff_deriv(i,nc) + cutoff_deriv(j,nc)) &
                                                                  + theta * rvec_en(a,i,nc) * (cutoff_deriv(i,nc) + cutoff_deriv(j,nc))

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


        if (ioptbf .gt. 0) then
            ! Here we first add the derivatives due to the cusp condition, and then zero the fixed variables.
            ! Ideally, we would not even have the fixed variables in the parameter list, but this is more complicated to implement.
            do a = 1, 3
                do k = 1, nparm_bf
                    if (inv_cusp_indices(k,1) .gt. 0) then
                        do kk = 2, ncparm_bf
                            if (inv_cusp_indices(k, kk) .eq. 0) exit
                            dquasi_dp(a,i,k) = dquasi_dp(a,i,k) + inv_cusp_parameters(k,kk) * dquasi_dp(a,i,inv_cusp_indices(k,kk))
                        end do
                    end if
                end do
               
                ! Add chain rule contribution for cutoff derivative:
                ! The cusp constraint is: p_dep = sum(cusp_parameters(k,kk) * p_free)
                ! When cutoff changes, the coefficients change:
                ! d(p_dep)/d(cutoff) = sum(cusp_cutoff_deriv(k,kk) * p_free)
                ! This contributes to dquasi/d(cutoff) via:
                ! dquasi/d(cutoff) += dquasi/d(p_dep) * d(p_dep)/d(cutoff)
                do k = 1, c_cuspconst*nctype*nspin_bf_ee
                    if (cusp_indices(k,1) .gt. 0) then
                        dp_dep_dcutoff = 0.0d0
                        do kk = 2, ncparm_bf
                            if (cusp_indices(k,kk) .eq. 0) exit
                            dp_dep_dcutoff = dp_dep_dcutoff + cusp_cutoff_deriv(k,kk) * parm_bf(cusp_indices(k,kk))
                        end do

                        ! All spin/component blocks share one cutoff per nuclear type.
                        idx = cutoff_c_offset + mod((cusp_indices(k,1) - een_coeff_offset - 1) / ncparm_bf, nctype)
                        dquasi_dp(a,i,idx) = dquasi_dp(a,i,idx) + dquasi_dp(a,i,cusp_indices(k,1)) * dp_dep_dcutoff
                    end if
                end do

                do k = 1, c_cuspconst*nctype*nspin_bf_ee
                    if (cusp_indices(k,1) .gt. 0) dquasi_dp(a,i,cusp_indices(k,1)) = 0.0d0
                end do
            end do
        endif
    end do

    ! rij = sqrt((x(1,1)-x(1,3))**2 + (x(2,1)-x(2,3))**2 + (x(3,1)-x(3,3))**2)
    ! do a=1,3
    !         print *, 'X:', (quasi_x(a,1)  - (quasi_x(a,3) ))/(rij)
    !         ! print *, 'd2X:', d2quasi_dx2(a,1,1), d2quasi_dx2(a,3,3)
    ! end do
    !     print *, 'done'




end subroutine rios_backflow



subroutine single_rios_backflow(iel, xold, xnew, quasi_x_new, dquasi_dx_new, d2quasi_dx2_new, indices)
    use precision_kinds, only: dp
    use system, only: nelec, nup, iwctype, ncent, nctype, cent
    use optwf_control, only: ioptci, ioptjas, ioptorb, ioptbf
    use m_backflow, only: parm_bf, nparm_bf, norda_bf, nordb_bf, nspin_bf_ee, nordc_bf, cutoff_scale, ncparm_bf, basis_klmn, cutoff_b_offset, cutoff_a_offset, cutoff_c_offset, ee_coeff_offset, en_coeff_offset, een_coeff_offset, ee_coeff_block_size, en_coeff_block_size, een_component_block_size
    use m_backflow, only: quasi_x, dquasi_dx, d2quasi_dx2, r_ee, rvec_ee, r_en, rvec_en, r_ee_gl, r_en_gl
    use m_backflow, only: single_r_ee, single_rvec_ee, single_r_en, single_rvec_en, single_r_ee_gl, single_r_en_gl
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
    integer :: i, j, k, a, b, idx, C
    integer :: ee_spin, offset_ee_ch, offset_een_ch
    integer :: nc, l, m, n, idx_phi, idx_theta, kk
    real(dp) :: inv_rij, inv_cutoff, tmp1, tmp2
    real(dp) :: cutoff1, cutoff2
    real(dp) :: fi, fpi(3), fppi(3), fj, fpj(3), fppj(3)
    real(dp) :: phi, theta, phipi(3), thetapi(3), phipj(3), thetapj(3), phippi, thetappi, phippj, thetappj
    real(dp) :: xtemp(3, nelec)

    quasi_x_new = quasi_x
    dquasi_dx_new = dquasi_dx
    !d2quasi_dx2_new = d2quasi_dx2  ! Not needed during single-electron loop

    indices = 0

    do k = 1, 3
        quasi_x_new(k, iel) = quasi_x_new(k, iel) - xold(k, iel) + xnew(k)
    end do
    indices(iel) = iel
   
    C = cutoff_scale

    if (nordb_bf .eq. 0) goto 30

    do j = 1, nelec
        if (j == iel) cycle

        if (nspin_bf_ee.eq.1) then
            ee_spin = 1
        else
            if ((iel.le.nup .and. j.le.nup) .or. (iel.gt.nup .and. j.gt.nup)) then
                ee_spin = 1
            else
                ee_spin = 2
            endif
        endif
        offset_ee_ch = ee_coeff_offset + (ee_spin-1)*ee_coeff_block_size

        cutoff = parm_bf(cutoff_b_offset)
        inv_cutoff = 1.0d0 / cutoff

        if (ee_spin .eq. 1) then
            b_one = parm_bf(offset_ee_ch+2) * cutoff/C
        else
            b_one = parm_bf(offset_ee_ch+1)
        endif

        delta(:) = xold(:, iel) - xold(:, j)
        rij = sqrt(delta(1)**2 + delta(2)**2 + delta(3)**2)

        if (rij < cutoff) then
            indices(j) = j

            eta=0.0d0
            etap=0.0d0
            ! etapp=0.0d0

            f = ((cutoff - rij)*inv_cutoff)**C
            inv_rij = 1.0d0 / rij

            cutoff1 = f/((cutoff - rij)*inv_cutoff)
            cutoff2 = cutoff1/((cutoff - rij)*inv_cutoff)

            tmp1 = -C*inv_cutoff * cutoff1 * inv_rij
            do a = 1, 3
                fp(a) = tmp1 * delta(a)
            end do
            ! fpp = -2.0d0 * C * inv_cutoff * inv_rij * cutoff1 + &
            !       + C * (C-1) * inv_cutoff * inv_cutoff * cutoff2


            eta = eta + b_one

            rr = rij
       
            tmp1 = inv_rij*inv_rij
            do k = 1, nordb_bf
                eta = eta + parm_bf(offset_ee_ch+k+1)*rr
                do a = 1, 3
                    etap(a) = etap(a) + parm_bf(offset_ee_ch+k+1) * delta(a) * tmp1 * k * rr
                enddo
                ! etapp = etapp + parm_bf(ee_coeff_offset+k+1) * k * rr * tmp1 * (k+1)
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

                    !tmp1 = 2.0d0 * etap(b) * delta(a) * fp(b)
                    !d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) - tmp1
                    !d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) - tmp1
                    !d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) + tmp1
                    !d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) + tmp1
                end do
                tmp1 = eta * f
                dquasi_dx_new(a,iel,a,iel) = dquasi_dx_new(a,iel,a,iel) - tmp1
                dquasi_dx_new(a,iel,a,j) = dquasi_dx_new(a,iel,a,j) + tmp1
                dquasi_dx_new(a,j,a,iel) = dquasi_dx_new(a,j,a,iel) + tmp1
                dquasi_dx_new(a,j,a,j) = dquasi_dx_new(a,j,a,j) - tmp1

                !tmp1 = etapp * delta(a) * f + eta * delta(a) * fpp
                !d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) - tmp1
                !d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) - tmp1
                !d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) + tmp1
                !d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) + tmp1

                !tmp1 = 2 * (eta * fp(a) + etap(a) * f)
                !d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) - tmp1
                !d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) - tmp1
                !d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) + tmp1
                !d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) + tmp1
            end do
        endif

        delta(:) = xnew(:) - xold(:, j)
        rij = sqrt(delta(1)**2 + delta(2)**2 + delta(3)**2)

        if (rij < cutoff) then
            indices(j) = j

            eta=0.0d0
            etap=0.0d0
            ! etapp=0.0d0

            f = ((cutoff - rij)*inv_cutoff)**C
            inv_rij = 1.0d0 / rij

            cutoff1 = f/((cutoff - rij)*inv_cutoff)
            cutoff2 = cutoff1/((cutoff - rij)*inv_cutoff)

            tmp1 = -C*inv_cutoff * cutoff1 * inv_rij
            do a = 1, 3
                fp(a) = tmp1 * delta(a)
            end do
            ! fpp = -2.0d0 * C * inv_cutoff * inv_rij * cutoff1 + &
            !       + C * (C-1) * inv_cutoff * inv_cutoff * cutoff2


            eta = eta + b_one

            rr = rij
       
            tmp1 = inv_rij*inv_rij
            do k = 1, nordb_bf
                eta = eta + parm_bf(offset_ee_ch+k+1)*rr
                do a = 1, 3
                    etap(a) = etap(a) + parm_bf(offset_ee_ch+k+1) * delta(a) * tmp1 * k * rr
                enddo
                ! etapp = etapp + parm_bf(ee_coeff_offset+k+1) * k * rr * tmp1 * (k+1)
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

                    !tmp1 = 2.0d0 * etap(b) * delta(a) * fp(b)
                    !d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + tmp1
                    !d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) + tmp1
                    !d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) - tmp1
                    !d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) - tmp1
                end do
                tmp1 = eta * f
                dquasi_dx_new(a,iel,a,iel) = dquasi_dx_new(a,iel,a,iel) + tmp1
                dquasi_dx_new(a,iel,a,j) = dquasi_dx_new(a,iel,a,j) - tmp1
                dquasi_dx_new(a,j,a,iel) = dquasi_dx_new(a,j,a,iel) - tmp1
                dquasi_dx_new(a,j,a,j) = dquasi_dx_new(a,j,a,j) + tmp1

                !tmp1 = etapp * delta(a) * f + eta * delta(a) * fpp
                !d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + tmp1
                !d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) + tmp1
                !d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) - tmp1
                !d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) - tmp1

                !tmp1 = 2 * (eta * fp(a) + etap(a) * f)
                !d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + tmp1
                !d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) + tmp1
                !d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) - tmp1
                !d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) - tmp1
            end do
        endif
   
    end do


30  if (norda_bf .eq. 0) goto 40

    do j = 1, ncent
        idx = (iwctype(j)-1)*(norda_bf+1)
        cutoff = parm_bf(cutoff_a_offset + iwctype(j) - 1)
        inv_cutoff = 1.0d0 / cutoff

        a_one = parm_bf(en_coeff_offset+idx+1)


        delta(:) = xold(:, iel) - cent(:, j)
        rij = sqrt(delta(1)**2 + delta(2)**2 + delta(3)**2)

        if (rij < cutoff) then
            indices(iel) = iel
            inv_rij = 1.0d0 / rij

            eta=0.0d0
            etap=0.0d0
            ! etapp=0.0d0

            f = ((cutoff - rij)*inv_cutoff)**C

            cutoff1 = f/((cutoff - rij)*inv_cutoff)
            cutoff2 = cutoff1/((cutoff - rij)*inv_cutoff)

            tmp1 = -C*inv_cutoff * cutoff1 * inv_rij
            do a = 1, 3
                fp(a) = tmp1 * delta(a)
            end do
            ! fpp = -2.0d0 * C * inv_cutoff * inv_rij * cutoff1 + &
            !       + C * (C-1) * inv_cutoff * inv_cutoff * cutoff2

            eta  = eta + a_one

            rr = rij

            tmp1 = inv_rij*inv_rij
            do k = 1, norda_bf
                eta = eta + parm_bf(en_coeff_offset + idx+ k+1)*rr
                do a = 1, 3
                    etap(a) = etap(a) + parm_bf(en_coeff_offset + idx+ k+1) * k * rr * delta(a)*tmp1
                enddo
                ! etapp = etapp + parm_bf(en_coeff_offset + idx+ k+2) * k * rr * tmp1 * (k+1)
                rr = rr * rij
            end do

            do a = 1, 3
                quasi_x_new(a,iel) = quasi_x_new(a,iel) - eta * delta(a) * f
                do b = 1, 3
                    dquasi_dx_new(a,iel,b,iel) = dquasi_dx_new(a,iel,b,iel) - (&
                        etap(b) * delta(a) * f + eta * delta(a) * fp(b) )

                    !d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) - (&
                    !     2.0d0 * etap(b) * delta(a) * fp(b)  )
                end do
                dquasi_dx_new(a,iel,a,iel) = dquasi_dx_new(a,iel,a,iel) - eta * f
                !d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) - etapp * delta(a) * f - eta * delta(a) * fpp
                !d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) - 2 * (eta * fp(a) + etap(a) * f)
            end do
        endif
       
        delta(:) = xnew(:) - cent(:, j)
        rij = sqrt(delta(1)**2 + delta(2)**2 + delta(3)**2)
        if (rij < cutoff) then
            indices(iel) = iel
            inv_rij = 1.0d0 / rij

            eta=0.0d0
            etap=0.0d0
            ! etapp=0.0d0

            f = ((cutoff - rij)*inv_cutoff)**C

            cutoff1 = f/((cutoff - rij)*inv_cutoff)
            cutoff2 = cutoff1/((cutoff - rij)*inv_cutoff)

            tmp1 = -C*inv_cutoff * cutoff1 * inv_rij
            do a = 1, 3
                fp(a) = tmp1 * delta(a)
            end do
            ! fpp = -2.0d0 * C * inv_cutoff * inv_rij * cutoff1 + &
            !       + C * (C-1) * inv_cutoff * inv_cutoff * cutoff2

            eta  = eta + a_one

            rr = rij

            tmp1 = inv_rij*inv_rij
            do k = 1, norda_bf
                eta = eta + parm_bf(en_coeff_offset + idx+ k+1)*rr
                do a = 1, 3
                    etap(a) = etap(a) + parm_bf(en_coeff_offset + idx+ k+1) * k * rr * delta(a)*tmp1
                enddo
                ! etapp = etapp + parm_bf(en_coeff_offset + idx+ k+2) * k * rr * tmp1 * (k+1)
                rr = rr * rij
            end do

            do a = 1, 3
                quasi_x_new(a,iel) = quasi_x_new(a,iel) + eta * delta(a) * f
                do b = 1, 3
                    dquasi_dx_new(a,iel,b,iel) = dquasi_dx_new(a,iel,b,iel) + (&
                        etap(b) * delta(a) * f + eta * delta(a) * fp(b) )

                    !d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + (&
                    !    2.0d0 * etap(b) * delta(a) * fp(b)  )
                end do
                dquasi_dx_new(a,iel,a,iel) = dquasi_dx_new(a,iel,a,iel) + eta * f
                !d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + etapp * delta(a) * f + eta * delta(a) * fpp
                !d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + 2 * (eta * fp(a) + etap(a) * f)
            end do
        endif

    end do


40  if (nordc_bf .eq. 0) return


    !call single_rios_distances(xold, xold(:, iel), iel)
    !call rios_distances(xold)
    do j = 1, nelec
        if (iel == j) cycle
        if (nspin_bf_ee.eq.1) then
            ee_spin = 1
        else
            if ((iel.le.nup .and. j.le.nup) .or. (iel.gt.nup .and. j.gt.nup)) then
                ee_spin = 1
            else
                ee_spin = 2
            endif
        endif
        offset_een_ch = een_coeff_offset + (ee_spin-1)*een_component_block_size
        do nc = 1, ncent
            idx_phi = (iwctype(nc)-1)*ncparm_bf + offset_een_ch
            idx_theta = nspin_bf_ee*een_component_block_size + (iwctype(nc)-1)*ncparm_bf + offset_een_ch

            cutoff = parm_bf(cutoff_c_offset + iwctype(nc) - 1)
            inv_cutoff = 1.0d0 / cutoff

            if (r_en(iel,nc,0,2) > 0 .and. r_en(j,nc,0,2) > 0) then

                indices(j) = j

                phi = 0.00d0
                phipi = 0.0d0
                phipj = 0.0d0
                ! phippi = 0.0d0
                ! phippj = 0.0d0
                theta = 0.00d0
                thetapi = 0.0d0
                thetapj = 0.0d0
                ! thetappi = 0.0d0
                ! thetappj = 0.0d0

                do l = 0, nordc_bf
                    do m = 0, nordc_bf - l
                        do n = 0, nordc_bf - l - m
                            ! if (l > m) then
                            ! else
                                k = basis_klmn(l,m,n,iwctype(nc),1,ee_spin)
                                kk = basis_klmn(l,m,n,iwctype(nc),2,ee_spin)
                            ! end if
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

                                ! tmp1 = 2.0d0 * r_en_gl(iel,a,nc,l,2) * r_en(j,nc,m,2) * r_ee_gl(iel,a,j,n)
                                ! phippi = phippi + parm_bf(k) * tmp1
                                ! thetappi = thetappi + parm_bf(kk) * tmp1

                                ! tmp1 = 2.0d0 * r_en(iel,nc,l,2) * r_en_gl(j,a,nc,m,2) * r_ee_gl(iel,a,j,n)
                                ! phippj = phippj - parm_bf(k) * tmp1
                                ! thetappj = thetappj - parm_bf(kk) * tmp1
                            enddo

                            ! tmp1 = r_en_gl(iel,4,nc,l,2) * r_en(j,nc,m,2) * r_ee(iel,j,n)
                            ! phippi = phippi + parm_bf(k) * tmp1
                            ! thetappi = thetappi + parm_bf(kk) * tmp1

                            ! tmp1 = r_en(iel,nc,l,2) * r_en_gl(j,4,nc,m,2) * r_ee(iel,j,n)
                            ! phippj = phippj + parm_bf(k) * tmp1
                            ! thetappj = thetappj + parm_bf(kk) * tmp1

                            ! tmp1 = r_en(iel,nc,l,2) * r_en(j,nc,m,2) * r_ee_gl(iel,4,j,n)
                            ! phippi = phippi + parm_bf(k) * tmp1
                            ! thetappi = thetappi + parm_bf(kk) * tmp1
                            ! phippj = phippj + parm_bf(k) * tmp1
                            ! thetappj = thetappj + parm_bf(kk) * tmp1
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
 
                    !d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) - phippi * rvec_ee(a,iel,j) - thetappi * rvec_en(a,iel,nc)
                    !d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) - phippj * rvec_ee(a,iel,j) - thetappj * rvec_en(a,iel,nc)

                    !d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) - 2.0d0 * phipi(a) - 2.0d0 * thetapi(a)
                    !d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) + 2.0d0 * phipj(a)
                end do

                ! Now account for the contribution where j is the central electron
                phi = 0.0d0
                phipi = 0.0d0
                phipj = 0.0d0
                ! phippi = 0.0d0
                ! phippj = 0.0d0
                theta = 0.0d0
                thetapi = 0.0d0
                thetapj = 0.0d0
                ! thetappi = 0.0d0
                ! thetappj = 0.0d0

                do l = 0, nordc_bf
                    do m = 0, nordc_bf - l
                        do n = 0, nordc_bf - l - m
                            ! if (l > m) then
                            ! else
                                k = basis_klmn(l,m,n,iwctype(nc),1,ee_spin)
                                kk = basis_klmn(l,m,n,iwctype(nc),2,ee_spin)
                            ! end if
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

                                ! tmp1 = 2.0d0 * r_en_gl(j,a,nc,l,2) * r_en(iel,nc,m,2) * r_ee_gl(j,a,iel,n)
                                ! phippi = phippi + parm_bf(k) * tmp1
                                ! thetappi = thetappi + parm_bf(kk) * tmp1

                                ! tmp1 = 2.0d0 * r_en(j,nc,l,2) * r_en_gl(iel,a,nc,m,2) * r_ee_gl(j,a,iel,n)
                                ! phippj = phippj - parm_bf(k) * tmp1
                                ! thetappj = thetappj - parm_bf(kk) * tmp1
                            enddo

                            ! tmp1 = r_en_gl(j,4,nc,l,2) * r_en(iel,nc,m,2) * r_ee(j,iel,n)
                            ! phippi = phippi + parm_bf(k) * tmp1
                            ! thetappi = thetappi + parm_bf(kk) * tmp1

                            ! tmp1 = r_en(j,nc,l,2) * r_en_gl(iel,4,nc,m,2) * r_ee(j,iel,n)
                            ! phippj = phippj + parm_bf(k) * tmp1
                            ! thetappj = thetappj + parm_bf(kk) * tmp1

                            ! tmp1 = r_en(j,nc,l,2) * r_en(iel,nc,m,2) * r_ee_gl(j,4,iel,n)
                            ! phippi = phippi + parm_bf(k) * tmp1
                            ! thetappi = thetappi + parm_bf(kk) * tmp1
                            ! phippj = phippj + parm_bf(k) * tmp1
                            ! thetappj = thetappj + parm_bf(kk) * tmp1
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

                    !d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) - phippi * rvec_ee(a,j,iel) - thetappi * rvec_en(a,j,nc)
                    !d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) - phippj * rvec_ee(a,j,iel) - thetappj * rvec_en(a,j,nc)

                    !d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) - 2.0d0 * phipi(a) - 2.0d0 * thetapi(a)
                    !d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) + 2.0d0 * phipj(a)
                end do
            end if
        end do
    end do

   
    call single_rios_distances(xold, xnew, iel)
    !xtemp = xold
    !xtemp(:, iel) = xnew
    !call rios_distances(xtemp)


    do j = 1, nelec
        if (iel == j) cycle
        if (nspin_bf_ee.eq.1) then
            ee_spin = 1
        else
            if ((iel.le.nup .and. j.le.nup) .or. (iel.gt.nup .and. j.gt.nup)) then
                ee_spin = 1
            else
                ee_spin = 2
            endif
        endif
        offset_een_ch = een_coeff_offset + (ee_spin-1)*een_component_block_size
        do nc = 1, ncent
            idx_phi = (iwctype(nc)-1)*ncparm_bf + offset_een_ch
            idx_theta = nspin_bf_ee*een_component_block_size + (iwctype(nc)-1)*ncparm_bf + offset_een_ch

            cutoff = parm_bf(cutoff_c_offset + iwctype(nc) - 1)
            inv_cutoff = 1.0d0 / cutoff

            if (single_r_en(nc,0,2) > 0 .and. r_en(j,nc,0,2) > 0) then

                indices(j) = j

                phi = 0.00d0
                phipi = 0.0d0
                phipj = 0.0d0
                ! phippi = 0.0d0
                ! phippj = 0.0d0
                theta = 0.00d0
                thetapi = 0.0d0
                thetapj = 0.0d0
                ! thetappi = 0.0d0
                ! thetappj = 0.0d0

                do l = 0, nordc_bf
                    do m = 0, nordc_bf - l
                        do n = 0, nordc_bf - l - m
                            ! if (l > m) then
                            ! else
                                k = basis_klmn(l,m,n,iwctype(nc),1,ee_spin)
                                kk = basis_klmn(l,m,n,iwctype(nc),2,ee_spin)
                            ! end if
                            phi = phi + parm_bf(k) * single_r_en(nc,l,2) * r_en(j,nc,m,2) * single_r_ee(j,n)
                            theta = theta + parm_bf(kk) * single_r_en(nc,l,2) * r_en(j,nc,m,2) * single_r_ee(j,n)
                            do a = 1, 3
                                tmp1 = single_r_en_gl(a,nc,l,2) * r_en(j,nc,m,2) * single_r_ee(j,n)
                                phipi(a) = phipi(a) + parm_bf(k) * tmp1
                                thetapi(a) = thetapi(a) + parm_bf(kk) * tmp1

                                tmp1 = single_r_en(nc,l,2) * r_en_gl(j,a,nc,m,2) * single_r_ee(j,n)
                                phipj(a) = phipj(a) + parm_bf(k) * tmp1
                                thetapj(a) = thetapj(a) + parm_bf(kk) * tmp1

                                tmp1 = single_r_en(nc,l,2) * r_en(j,nc,m,2) * single_r_ee_gl(a,j,n)
                                phipi(a) = phipi(a) + parm_bf(k) * tmp1
                                phipj(a) = phipj(a) - parm_bf(k) * tmp1
                                thetapi(a) = thetapi(a) + parm_bf(kk) * tmp1
                                thetapj(a) = thetapj(a) - parm_bf(kk) * tmp1

                                ! tmp1 = 2.0d0 * single_r_en_gl(a,nc,l,2) * r_en(j,nc,m,2) * single_r_ee_gl(a,j,n)
                                ! phippi = phippi + parm_bf(k) * tmp1
                                ! thetappi = thetappi + parm_bf(kk) * tmp1

                                ! tmp1 = 2.0d0 * single_r_en(nc,l,2) * r_en_gl(j,a,nc,m,2) * single_r_ee_gl(a,j,n)
                                ! phippj = phippj - parm_bf(k) * tmp1
                                ! thetappj = thetappj - parm_bf(kk) * tmp1
                            enddo

                            ! tmp1 = single_r_en_gl(4,nc,l,2) * r_en(j,nc,m,2) * single_r_ee(j,n)
                            ! phippi = phippi + parm_bf(k) * tmp1
                            ! thetappi = thetappi + parm_bf(kk) * tmp1

                            ! tmp1 = single_r_en(nc,l,2) * r_en_gl(j,4,nc,m,2) * single_r_ee(j,n)
                            ! phippj = phippj + parm_bf(k) * tmp1
                            ! thetappj = thetappj + parm_bf(kk) * tmp1

                            ! tmp1 = single_r_en(nc,l,2) * r_en(j,nc,m,2) * single_r_ee_gl(4,j,n)
                            ! phippi = phippi + parm_bf(k) * tmp1
                            ! thetappi = thetappi + parm_bf(kk) * tmp1
                            ! phippj = phippj + parm_bf(k) * tmp1
                            ! thetappj = thetappj + parm_bf(kk) * tmp1
                        end do
                    end do
                end do
                do a = 1, 3
                    quasi_x_new(a,iel) = quasi_x_new(a,iel) + phi * single_rvec_ee(a,j) + theta * single_rvec_en(a,nc)

                    do b = 1, 3
                        dquasi_dx_new(a,iel,b,iel) = dquasi_dx_new(a,iel,b,iel) + phipi(b) * single_rvec_ee(a,j) + thetapi(b) * single_rvec_en(a,nc)                       
                        dquasi_dx_new(a,iel,b,j) = dquasi_dx_new(a,iel,b,j) + phipj(b) * single_rvec_ee(a,j) + thetapj(b) * single_rvec_en(a,nc)
                    end do
                    dquasi_dx_new(a,iel,a,iel) = dquasi_dx_new(a,iel,a,iel) + phi + theta
                    dquasi_dx_new(a,iel,a,j) = dquasi_dx_new(a,iel,a,j) - phi

                    !d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + phippi * single_rvec_ee(a,j) + thetappi * single_rvec_en(a,nc)
                    !d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) + phippj * single_rvec_ee(a,j) + thetappj * single_rvec_en(a,nc)

                    !d2quasi_dx2_new(a,iel,iel) = d2quasi_dx2_new(a,iel,iel) + 2.0d0 * phipi(a) + 2.0d0 * thetapi(a)
                    !d2quasi_dx2_new(a,iel,j) = d2quasi_dx2_new(a,iel,j) - 2.0d0 * phipj(a)
                end do

                ! Contribution with j as the central electron (new configuration)
                phi = 0.0d0
                phipi = 0.0d0
                phipj = 0.0d0
                ! phippi = 0.0d0
                ! phippj = 0.0d0
                theta = 0.0d0
                thetapi = 0.0d0
                thetapj = 0.0d0
                ! thetappi = 0.0d0
                ! thetappj = 0.0d0

                do l = 0, nordc_bf
                    do m = 0, nordc_bf - l
                        do n = 0, nordc_bf - l - m
                            ! if (l > m) then
                            ! else
                                k = basis_klmn(l,m,n,iwctype(nc),1,ee_spin)
                                kk = basis_klmn(l,m,n,iwctype(nc),2,ee_spin)
                            ! end if
                            phi = phi + parm_bf(k) * r_en(j,nc,l,2) * single_r_en(nc,m,2) * single_r_ee(j,n)
                            theta = theta + parm_bf(kk) * r_en(j,nc,l,2) * single_r_en(nc,m,2) * single_r_ee(j,n)
                            do a = 1, 3
                                tmp1 = r_en_gl(j,a,nc,l,2) * single_r_en(nc,m,2) * single_r_ee(j,n)
                                phipi(a) = phipi(a) + parm_bf(k) * tmp1
                                thetapi(a) = thetapi(a) + parm_bf(kk) * tmp1

                                tmp1 = r_en(j,nc,l,2) * single_r_en_gl(a,nc,m,2) * single_r_ee(j,n)
                                phipj(a) = phipj(a) + parm_bf(k) * tmp1
                                thetapj(a) = thetapj(a) + parm_bf(kk) * tmp1

                                tmp1 = -r_en(j,nc,l,2) * single_r_en(nc,m,2) * single_r_ee_gl(a,j,n)
                                phipi(a) = phipi(a) + parm_bf(k) * tmp1
                                phipj(a) = phipj(a) - parm_bf(k) * tmp1
                                thetapi(a) = thetapi(a) + parm_bf(kk) * tmp1
                                thetapj(a) = thetapj(a) - parm_bf(kk) * tmp1

                                ! tmp1 = -2.0d0 * r_en_gl(j,a,nc,l,2) * single_r_en(nc,m,2) * single_r_ee_gl(a,j,n)
                                ! phippi = phippi + parm_bf(k) * tmp1
                                ! thetappi = thetappi + parm_bf(kk) * tmp1

                                ! tmp1 = -2.0d0 * r_en(j,nc,l,2) * single_r_en_gl(a,nc,m,2) * single_r_ee_gl(a,j,n)
                                ! phippj = phippj - parm_bf(k) * tmp1
                                ! thetappj = thetappj - parm_bf(kk) * tmp1
                            enddo

                            ! tmp1 = r_en_gl(j,4,nc,l,2) * single_r_en(nc,m,2) * single_r_ee(j,n)
                            ! phippi = phippi + parm_bf(k) * tmp1
                            ! thetappi = thetappi + parm_bf(kk) * tmp1

                            ! tmp1 = r_en(j,nc,l,2) * single_r_en_gl(4,nc,m,2) * single_r_ee(j,n)
                            ! phippj = phippj + parm_bf(k) * tmp1
                            ! thetappj = thetappj + parm_bf(kk) * tmp1

                            ! tmp1 = r_en(j,nc,l,2) * single_r_en(nc,m,2) * single_r_ee_gl(4,j,n)
                            ! phippi = phippi + parm_bf(k) * tmp1
                            ! thetappi = thetappi + parm_bf(kk) * tmp1
                            ! phippj = phippj + parm_bf(k) * tmp1
                            ! thetappj = thetappj + parm_bf(kk) * tmp1
                        end do
                    end do
                end do

                do a = 1, 3
                    quasi_x_new(a,j) = quasi_x_new(a,j) - phi * single_rvec_ee(a,j) + theta * rvec_en(a,j,nc)

                    do b = 1, 3
                        dquasi_dx_new(a,j,b,j) = dquasi_dx_new(a,j,b,j) - phipi(b) * single_rvec_ee(a,j) + thetapi(b) * rvec_en(a,j,nc)                       
                        dquasi_dx_new(a,j,b,iel) = dquasi_dx_new(a,j,b,iel) - phipj(b) * single_rvec_ee(a,j) + thetapj(b) * rvec_en(a,j,nc)
                    end do
                    dquasi_dx_new(a,j,a,j) = dquasi_dx_new(a,j,a,j) + phi + theta
                    dquasi_dx_new(a,j,a,iel) = dquasi_dx_new(a,j,a,iel) - phi

                    !d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) - phippi * single_rvec_ee(a,j) + thetappi * rvec_en(a,j,nc)
                    !d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) - phippj * single_rvec_ee(a,j) + thetappj * rvec_en(a,j,nc)

                    !d2quasi_dx2_new(a,j,j) = d2quasi_dx2_new(a,j,j) + 2.0d0 * phipi(a) + 2.0d0 * thetapi(a)
                    !d2quasi_dx2_new(a,j,iel) = d2quasi_dx2_new(a,j,iel) - 2.0d0 * phipj(a)
                end do
            end if
        end do
    end do
end subroutine single_rios_backflow


subroutine backflow_accept(iel)
    use precision_kinds, only: dp
    use system, only: nelec,ncent
    use m_backflow, only: quasi_x, dquasi_dx, d2quasi_dx2, r_ee, rvec_ee, r_en, rvec_en, r_ee_gl, r_en_gl, maxord
    use m_backflow, only: single_r_ee, single_rvec_ee, single_r_en, single_rvec_en, single_r_ee_gl, single_r_en_gl
    use m_backflow, only: quasi_x_new, dquasi_dx_new, d2quasi_dx2_new
    implicit none
    integer, intent(in) :: iel

    integer :: i, nc, cc, k, p

    quasi_x = quasi_x_new
    dquasi_dx = dquasi_dx_new
    !d2quasi_dx2 = d2quasi_dx2_new  ! Not needed during single-electron loop, maybe for forces later

    do cc = 1, 2
        do nc = 1, ncent
            do k = 1, 3
                rvec_en(k, iel, nc) = single_rvec_en(k, nc)
            enddo
            do p = 0, maxord
                r_en(iel, nc, p, 2) = single_r_en(nc, p, 2)
                r_en_gl(iel, 1, nc, p, 2) = single_r_en_gl(1, nc, p, 2)
                r_en_gl(iel, 2, nc, p, 2) = single_r_en_gl(2, nc, p, 2)
                r_en_gl(iel, 3, nc, p, 2) = single_r_en_gl(3, nc, p, 2)
                ! r_en_gl(iel, 4, nc, p, 2) = single_r_en_gl(4, nc, p, 2)
            end do
        end do
    end do
           


    do i = 1, nelec
        do k = 1, 3
            rvec_ee(k, iel, i) = single_rvec_ee(k, i)
            rvec_ee(k, i, iel) = -single_rvec_ee(k, i)
        enddo
        do p = 0, maxord
            r_ee(iel, i, p) = single_r_ee(i, p)
            r_ee(i, iel, p) = single_r_ee(i, p)
            r_ee_gl(iel, 1, i, p) = single_r_ee_gl(1, i, p)
            r_ee_gl(iel, 2, i, p) = single_r_ee_gl(2, i, p)
            r_ee_gl(iel, 3, i, p) = single_r_ee_gl(3, i, p)
            ! r_ee_gl(iel, 4, i, p) = single_r_ee_gl(4, i, p)
            r_ee_gl(i, 1, iel, p) = -single_r_ee_gl(1, i, p)
            r_ee_gl(i, 2, iel, p) = -single_r_ee_gl(2, i, p)
            r_ee_gl(i, 3, iel, p) = -single_r_ee_gl(3, i, p)
            ! r_ee_gl(i, 4, iel, p) = single_r_ee_gl(4, i, p)
        enddo
    end do

end subroutine backflow_accept

end module
