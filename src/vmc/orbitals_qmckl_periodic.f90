module orbitals_qmckl_periodic_mod

contains

subroutine orbitals_qmckl_periodic(x,rvec_en,r_en)

    use orbval, only: ddorb, dorb, nadorb, orb
    use precision_kinds, only: dp
    use slater, only: norb
    use system, only: ncent_tot, nelec
    use periodic, only : n_images, ell
    use find_pimage, only: find_image_pbc_no_norm

    use const
    use qmckl_data

    implicit none

    real(dp), allocatable :: mo_vgl_qmckl(:,:,:)
    integer :: rc, ictx
    integer*8 :: n8
    real(dp), dimension(3,nelec) :: ximage, ximage0
    real(dp), allocatable :: ao_vgl_qmckl(:,:,:)
    real(dp), allocatable :: ao_qmckl(:,:,:)
    integer*8 :: na8, i_image, ivgl, i_basis, i_elec
    character*(1024) :: err_message = ''

    integer :: i, iorb, k, m

    real(dp), dimension(3,*) :: x
    real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
    real(dp), dimension(nelec,ncent_tot) :: r_en

    if (use_qmckl_jastrow) then
        ictx = qmckl_no_ctx-1
    else
        ictx = qmckl_no_ctx
    end if

    
    ! get number MO's
    rc = qmckl_get_mo_basis_mo_num(qmckl_ctx(ictx), n8)
    if (rc /= QMCKL_SUCCESS) then
        print *, 'Error getting mo_num from QMCkl'
        stop
    end if

    k=1 ! until state specific orbitals can be used
    orb(1:nelec,1:norb+nadorb,k) = 0.d0
    dorb(1:norb+nadorb,1:nelec,:,k) = 0.d0
    ddorb(1:norb+nadorb,1:nelec,k) = 0.d0
    allocate(mo_vgl_qmckl(n8, 5, nelec))

    ! Apply PBC and translate coordinates to the central image
    ximage0 = x(:,1:nelec)
    do i=1,nelec
       call find_image_pbc_no_norm(ximage0(:,i))
    enddo



    do i_image=1,n_images          
        !get electron coordinates of the image  
           ximage(1,:) = ximage0(1,:) - ell(1,i_image)
           ximage(2,:) = ximage0(2,:) - ell(2,i_image)
           ximage(3,:) = ximage0(3,:) - ell(3,i_image)

        ! Send electron coordinates to QMCkl to compute the MOs at these positions
        rc = qmckl_set_point(qmckl_ctx(ictx), 'N', nelec*1_8, ximage, nelec*3_8)
        if (rc /= QMCKL_SUCCESS) then
            print *, 'Error setting electron coordinates in QMCkl'
            stop
        end if

        ! Compute the MOs
        rc = qmckl_get_mo_basis_mo_vgl_inplace(qmckl_ctx(ictx), mo_vgl_qmckl, n8*nelec*5_8)

        if (rc /= QMCKL_SUCCESS) then
            print *, 'Error getting MOs from QMCkl'
            stop
        end if

        ! pass computed qmckl orbitals back to champ
        do i=1,nelec
            do iorb=1,norb+nadorb
                orb  (  i,iorb,k) = orb  (  i,iorb,k) + mo_vgl_qmckl(iorb,1,i)
                dorb (iorb,i,1,k) = dorb (iorb,i,1,k) + mo_vgl_qmckl(iorb,2,i)
                dorb (iorb,i,2,k) = dorb (iorb,i,2,k) + mo_vgl_qmckl(iorb,3,i)
                dorb (iorb,i,3,k) = dorb (iorb,i,3,k) + mo_vgl_qmckl(iorb,4,i)
                ddorb(  iorb,i,k) = ddorb(  iorb,i,k) + mo_vgl_qmckl(iorb,5,i)
            end do
        end do
    end do
    deallocate(mo_vgl_qmckl)
    
return
end

subroutine orbitalse_qmckl_periodic(iel,x,rvec_en,r_en,iflag)

    use multislatern, only: ddorbn, dorbn, orbn
    use precision_kinds, only: dp
    use slater, only: norb
    use system, only: ncent_tot, nelec
    use periodic, only : n_images, ell
    use find_pimage, only: find_image_pbc_no_norm

    use qmckl_data


    implicit none

    integer :: iel, iflag, iorb, ictx, k

    real(dp), dimension(3,*) :: x
    real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
    real(dp), dimension(nelec,ncent_tot) :: r_en

    real(dp), allocatable :: mo_vgl_qmckl(:,:)
    integer :: rc, i_image
    integer*8 :: n8
    character*(1024) :: err_message = ''
    real(dp), dimension(3) :: ximage, ximage0

    if (iflag .eq. 0) then
        ictx = 1
    else
        ictx = 2
    end if

    rc = qmckl_get_mo_basis_mo_num(qmckl_ctx(ictx), n8)
    if (rc /= QMCKL_SUCCESS) then
        print *, 'Error getting mo_num from QMCkl'
        stop
    end if
    
    allocate(mo_vgl_qmckl(n8, 5))
    
    k=1  ! until state-specific orbitals can use QMCKL
    orbn(1:norb, k) = 0.d0
    dorbn(1:norb, 1:3, k) = 0.d0
    if(iflag.gt.0) ddorbn(1:norb, k) = 0.d0
    ! Apply PBC and translate coordinates to the central image
    ximage0 = x(:,iel)
    call find_image_pbc_no_norm(ximage0)

    do i_image=1, n_images
        !get electron coordinates of the image
        ximage(1) = ximage0(1) - ell(1,i_image)
        ximage(2) = ximage0(2) - ell(2,i_image)
        ximage(3) = ximage0(3) - ell(3,i_image)

        ! set one electron coordinates
        rc = qmckl_set_point(qmckl_ctx(ictx), 'N', 1_8, ximage, 3_8)
        if (rc /= QMCKL_SUCCESS) then
            print *, 'Error setting electron coords orbitalse'
            call qmckl_last_error(qmckl_ctx(ictx),err_message)
            print *, trim(err_message)
            call abort()
        end if

        ! Compute the MOs
        rc = qmckl_get_mo_basis_mo_vgl_inplace(qmckl_ctx(ictx), mo_vgl_qmckl, n8*5_8)

        do iorb=1,norb
            orbn(iorb,k)= orbn(iorb,k) + mo_vgl_qmckl(iorb,1)
            dorbn(iorb,1,k)= dorbn(iorb,1,k) + mo_vgl_qmckl(iorb,2)
            dorbn(iorb,2,k)= dorbn(iorb,2,k) + mo_vgl_qmckl(iorb,3)
            dorbn(iorb,3,k)= dorbn(iorb,3,k) + mo_vgl_qmckl(iorb,4)
            if(iflag.gt.0) ddorbn(iorb,k)= ddorbn(iorb,k) + mo_vgl_qmckl(iorb,5)
        enddo
    enddo
    deallocate(mo_vgl_qmckl)
return
end

subroutine orbitals_quad_qmckl_periodic(nxquad,xquad,rvec_en,r_en,orbn,dorbn,da_orbn,iwforb)

    use m_force_analytic, only: iforce_analy
    use multiple_geo, only: iwf
    use optwf_control, only: ioptorb
    use optwf_control, only: method
    use orbval,  only: nadorb
    use phifun,  only: dphin,n0_ibasis,n0_ic,n0_nbasis,phin
    use precision_kinds, only: dp
    use qua,     only: nquad
    use slater,  only: coef,norb
    use sr_mod,  only: i_sr_rescale
    use system,  only: iwctype,ncent,ncent_tot,nelec
    use vmc_mod, only: norb_tot, nwftypeorb
    use error,   only: fatal_error
    use contrl_file, only: ounit
    use periodic, only : n_images, ell, rlatt, rlatt_inv
    use find_pimage, only: find_image_pbc_no_norm
    
    use qmckl_data

    implicit none

    integer :: ic, iq, ictx, i, j
    integer :: iorb, k, nxquad, iwforb
    integer :: nadorb_sav

    real(dp), dimension(3,*) :: xquad
    real(dp), dimension(nquad*nelec*2, ncent_tot) :: r_en
    real(dp), dimension(3,nquad*nelec*2, ncent_tot) :: rvec_en
    real(dp), dimension(norb_tot, *) :: orbn
    real(dp), dimension(norb_tot, nquad*nelec*2, 3) :: dorbn
    real(dp), dimension(norb,3,nxquad,ncent_tot) :: da_orbn

    real(dp), allocatable :: mo_qmckl(:,:)
    integer :: rc, i_image
    integer*8 :: n8
    real(dp), dimension(3, nxquad) :: ximage, ximage0
    real(dp), dimension(3, nxquad) :: xquad_work
    real(dp), dimension(3) :: s_frac

    nadorb_sav=nadorb

    ictx = 2

    if(ioptorb.eq.0.or.(method(1:3).ne.'lin'.and.i_sr_rescale.eq.0)) then
       ictx=1
       nadorb = 0
    end if
    
    rc = qmckl_get_mo_basis_mo_num(qmckl_ctx(ictx), n8)
    if (rc /= QMCKL_SUCCESS) then
        print *, 'orbitals quad Error getting mo_num from QMCkl'
        print *, "n8", n8
        stop
    end if

    allocate(mo_qmckl(n8, nxquad))
    orbn(1:norb+nadorb,1:nxquad) = 0.d0
    dorbn(1:norb,1:nxquad,1:3) = 0.d0

    !transform quadrature points to the central image 
    ximage0 = xquad(:,1:nxquad)
    do i=1,nxquad
       call find_image_pbc_no_norm(ximage0(:,i))
    enddo

    do i_image=1,n_images
        !get image quadrature points coordinates
        ximage(1,:) = ximage0(1,:) - ell(1,i_image)
        ximage(2,:) = ximage0(2,:) - ell(2,i_image)
        ximage(3,:) = ximage0(3,:) - ell(3,i_image)

        ! Send electron coordinates to QMCkl to compute the MOs at these positions
        rc = qmckl_set_point(qmckl_ctx(ictx), 'N', nxquad*1_8, ximage, nxquad*3_8)
        if (rc /= QMCKL_SUCCESS) then
            print *, 'orbitals quad Error setting electron coordinates in QMCkl'
            print *, "nxquad", nxquad
            stop
        end if

        if(iforce_analy.gt.0) then

            if (ictx .eq. 2) then
                rc = qmckl_set_point(qmckl_ctx(1), 'N', nxquad*1_8, ximage, nxquad*3_8)
                if (rc /= QMCKL_SUCCESS) then
                    stop
                end if
            end if

            rc = qmckl_get_forces_mo_value_inplace(qmckl_ctx(1), da_orbn, nxquad*norb*3_8*ncent_tot)
            if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl MO forces.')

            do ic=1,ncent
                do iq=1,nxquad
                    do k =1,3
                        do iorb=1,norb
                            dorbn(iorb,iq,k)=dorbn(iorb,iq,k)-da_orbn(iorb,k,iq,ic)
                        enddo
                    enddo
                enddo
            enddo
        endif


        ! Compute the MOs
        rc = qmckl_get_mo_basis_mo_value_inplace(qmckl_ctx(ictx), mo_qmckl, nxquad*n8)

        if (rc /= QMCKL_SUCCESS) then
            print *, 'Error orbitals quad getting MOs from QMCkl'
            stop
        end if

        orbn(1:norb+nadorb,1:nxquad) = orbn(1:norb+nadorb,1:nxquad) + mo_qmckl(1:norb+nadorb,1:nxquad)
    end do

    deallocate(mo_qmckl)

    nadorb = nadorb_sav
return
end

subroutine init_context_qmckl(update_coef)

    use qmckl_data
    use slater,  only: norb, coef
    use system,  only: ncent, znuc, iwctype
    use orbval,  only: nadorb
    use vmc_mod, only: norb_tot
    use precision_kinds, only: dp
    use coefs,   only: nbasis
    use pseudo,  only: nloc
    use contrl_file, only: ounit
    use error,   only: fatal_error

    implicit none

    integer(qmckl_exit_code)   :: rc
    integer*8                  :: n8
    integer*8                  :: ncheck, ictx, ictx_lim
    integer*8                  :: norb_qmckl(qmckl_no_ctx_max)
    integer, allocatable       :: keep(:)
    character*(1024)           :: err_message = ''
    integer                    :: k
  
    logical                    :: do_nucl_fitcusp
    real(dp), allocatable      :: nucl_fitcusp_radius(:)
    real(dp), parameter        :: a_cusp = 1.74891d0
    real(dp), parameter        :: b_cusp = 0.126057d0
    character(len=100)         :: int_format     = '(A, T40, ":: ", T50, I0)'
    character(len=100)         :: array_format   = '(A, "(",I0,")", T40, ":: ", T42, F25.16)'
    logical                    :: update_coef

    if (use_qmckl_jastrow) then
        ictx_lim = qmckl_no_ctx-1
    else
        ictx_lim = qmckl_no_ctx
    end if

    do ictx = 1, ictx_lim
        rc = qmckl_set_numprec_precision(qmckl_ctx(ictx), 53) ! 24
        if (rc .ne. QMCKL_SUCCESS) call fatal_error('INPUT: QMCkl error: Unable to set precision')
    end do

    write(ounit, *) " QMCkl precision set to 53 bits"


    do ictx=1,ictx_lim
        rc = qmckl_get_mo_basis_mo_num(qmckl_ctx(ictx), n8)
        if (rc /= QMCKL_SUCCESS) call fatal_error('INPUT: QMCkl getting mo_num from trexio file')

        write(ounit,int_format) "QMCkl number mo found", n8
        if(n8.ne.norb_tot) call fatal_error('INPUT: QMCkl getting wrong number of orbitals')
        
       if (update_coef) then
          rc = qmckl_set_mo_basis_coefficient(qmckl_ctx(ictx), coef(:,:, 1), nbasis*norb_tot*1_8)
          if (rc /= QMCKL_SUCCESS) call fatal_error('INPUT: QMCkl setting orbital coefficients')
       end if
    enddo


    n8 = norb_tot*1_8

    allocate(keep(n8))

    ! maximum mo's number to occupied in determinants and for optimization
    norb_qmckl(1)=norb
    norb_qmckl(2)=norb+nadorb
    norb_qmckl(2)=min(norb_qmckl(2),norb_tot) ! perhaps done before

    do ictx=1,ictx_lim

      if (n8 > norb_qmckl(ictx)) then

       ! selecting range of orbitals to compute qith QMCkl
       keep(1:norb_qmckl(ictx)) = 1
       keep((norb_qmckl(ictx)+1):n8) = 0

       rc = qmckl_mo_basis_select_mo(qmckl_ctx(ictx), keep, n8)
       if (rc /= QMCKL_SUCCESS) write(ounit,*) 'Error 01 selecting MOs in verify orbitals'

       ! getting new number of orbitals to be computed
       rc = qmckl_get_mo_basis_mo_num(qmckl_ctx(ictx), ncheck)
       if (rc /= QMCKL_SUCCESS) call fatal_error('INPUT: QMCkl mo_num from verify orbitals')
       write(ounit,int_format) "QMCkl number of orbitals after mo's selec", ncheck
       write(ounit,int_format) "QMCkl norb_qmckl after mo's selec", norb_qmckl(ictx)

       if (ncheck /= norb_qmckl(ictx)) call fatal_error('INPUT: Problem in MO selection in QMCkl verify orb')
      endif

    enddo
    ! deallocate keep
    deallocate(keep)

    if(nloc.eq.0) then
      allocate(nucl_fitcusp_radius(ncent))
      do_nucl_fitcusp = .true.

      if(.not. do_nucl_fitcusp) then
         nucl_fitcusp_radius = 0.d0
       else
        do k=1,ncent
          nucl_fitcusp_radius(k) = 1.d0/(a_cusp*znuc(iwctype(k))+b_cusp)
          write(ounit, array_format) "Radius fit cusps for atom", k, nucl_fitcusp_radius(k)
        enddo

        ! Avoid dummy atoms
        do k=1,ncent
         if (znuc(iwctype(k)) < 5.d-1) then
           nucl_fitcusp_radius(k) = 0.d0
         endif
        enddo

        do ictx=1,ictx_lim
          write(ounit, *) "Context for QMCKl set mo basis r cusp  ", qmckl_ctx(ictx)
          rc = qmckl_set_mo_basis_r_cusp(qmckl_ctx(ictx),dble(nucl_fitcusp_radius(:)), int(ncent,8))
          write(ounit, *) "Status QMCKl set mo basis r cusp  ", rc

          if (rc /= QMCKL_SUCCESS) call fatal_error('PARSER: QMCkl error: Unable to set cusp parameters')
        enddo
      endif
    endif

return
end


end module
