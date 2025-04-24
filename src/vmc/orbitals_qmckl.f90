module orbitals_qmckl_mod

contains

subroutine orbitals_qmckl(x,rvec_en,r_en)

    use orbval, only: ddorb, dorb, nadorb, orb
    use precision_kinds, only: dp
    use slater, only: norb
    use system, only: ncent_tot, nelec

    use const
    use qmckl_data

    implicit none

    real(dp), allocatable :: mo_vgl_qmckl(:,:,:)
    integer :: rc
    integer*8 :: n8
    real(dp), dimension(3,nelec) :: xqmckl
    real(dp), dimension(3,nelec) :: xelec
    real(dp), allocatable :: ao_vgl_qmckl(:,:,:)
    real(dp), allocatable :: ao_qmckl(:,:,:)
    integer*8 :: na8, i_image, ivgl, i_basis, i_elec
    real(dp), dimension(3) :: r_image
    real(dp) :: rnorm
    character*(1024) :: err_message = ''

    integer :: i, iorb, k

    real(dp), dimension(3,*) :: x
    real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
    real(dp), dimension(nelec,ncent_tot) :: r_en

    
    ! get number MO's
    rc = qmckl_get_mo_basis_mo_num(qmckl_ctx(qmckl_no_ctx-1), n8)
    if (rc /= QMCKL_SUCCESS) then
        print *, 'Error getting mo_num from QMCkl'
        stop
    end if

    allocate(mo_vgl_qmckl(n8, 5, nelec))
    ! Send electron coordinates to QMCkl to compute the MOs at these positions
    rc = qmckl_set_point(qmckl_ctx(qmckl_no_ctx-1), 'N', nelec*1_8, x, nelec*3_8)
    if (rc /= QMCKL_SUCCESS) then
        print *, 'Error setting electron coordinates in QMCkl'
        stop
    end if

    ! Compute the MOs
    rc = qmckl_get_mo_basis_mo_vgl_inplace(qmckl_ctx(qmckl_no_ctx-1), mo_vgl_qmckl, n8*nelec*5_8)

    if (rc /= QMCKL_SUCCESS) then
        print *, 'Error getting MOs from QMCkl'
        stop
    end if

    ! pass computed qmckl orbitals back to champ
    k=1 ! until state specific orbitals can be used
    do i=1,nelec
        do iorb=1,norb+nadorb
            orb  (  i,iorb,k) = mo_vgl_qmckl(iorb,1,i)
            dorb (iorb,i,1,k) = mo_vgl_qmckl(iorb,2,i)
            dorb (iorb,i,2,k) = mo_vgl_qmckl(iorb,3,i)
            dorb (iorb,i,3,k) = mo_vgl_qmckl(iorb,4,i)
            ddorb(  iorb,i,k) = mo_vgl_qmckl(iorb,5,i)
        end do
    end do

    deallocate(mo_vgl_qmckl)
return
end

subroutine orbitalse_qmckl(iel,x,rvec_en,r_en,iflag)

    use multislatern, only: ddorbn, dorbn, orbn
    use precision_kinds, only: dp
    use slater, only: norb
    use system, only: ncent_tot, nelec

    use qmckl_data


    implicit none

    integer :: h, iel, iflag, iorb, ictx, k

    real(dp), dimension(3,*) :: x
    real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
    real(dp), dimension(nelec,ncent_tot) :: r_en

    real(dp), allocatable :: mo_vgl_qmckl(:,:)
    integer :: rc
    integer*8 :: n8, na8, i_image, ivgl, i_basis
    character*(1024) :: err_message = ''
    real(dp), allocatable :: ao_vgl_qmckl(:,:,:)
    real(dp), allocatable :: ao_qmckl(:,:)
    real(dp), dimension(3) :: xiel
    real(dp) :: rnorm

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

    ! set one electron coordinates
    rc = qmckl_set_point(qmckl_ctx(ictx), 'N', 1_8, x(1:3,iel), 3_8)
    if (rc /= QMCKL_SUCCESS) then
        print *, 'Error setting electron coords orbitalse'
        call qmckl_last_error(qmckl_ctx(ictx),err_message)
        print *, trim(err_message)
        call abort()
    end if

    allocate(mo_vgl_qmckl(n8, 5))

    ! Compute the MOs
    rc = qmckl_get_mo_basis_mo_vgl_inplace(qmckl_ctx(ictx), mo_vgl_qmckl, n8*5_8)

    k=1  ! until state-specific orbitals can use QMCKL

    if(iflag.gt.0) then
        do iorb=1,norb
!     orbn(iorb)=mo_vgl_qmckl(iorb,1,1)
!     dorbn(iorb,1)=mo_vgl_qmckl(iorb,2,1)
!     dorbn(iorb,2)=mo_vgl_qmckl(iorb,3,1)
!     dorbn(iorb,3)=mo_vgl_qmckl(iorb,4,1)
!     ddorbn(iorb)=mo_vgl_qmckl(iorb,5,1)

            orbn(iorb,k)=mo_vgl_qmckl(iorb,1)
            dorbn(iorb,1,k)=mo_vgl_qmckl(iorb,2)
            dorbn(iorb,2,k)=mo_vgl_qmckl(iorb,3)
            dorbn(iorb,3,k)=mo_vgl_qmckl(iorb,4)
            ddorbn(iorb,k)=mo_vgl_qmckl(iorb,5)
        enddo
    else
        do iorb=1,norb
!     orbn(iorb)=mo_vgl_qmckl(iorb,1,1)
!     dorbn(iorb,1)=mo_vgl_qmckl(iorb,2,1)
!     dorbn(iorb,2)=mo_vgl_qmckl(iorb,3,1)
!     dorbn(iorb,3)=mo_vgl_qmckl(iorb,4,1)

            orbn(iorb,k)=mo_vgl_qmckl(iorb,1)
            dorbn(iorb,1,k)=mo_vgl_qmckl(iorb,2)
            dorbn(iorb,2,k)=mo_vgl_qmckl(iorb,3)
            dorbn(iorb,3,k)=mo_vgl_qmckl(iorb,4)

        enddo
    endif

    deallocate(mo_vgl_qmckl)
return
end


subroutine orbitals_quad_qmckl(nxquad,xquad,rvec_en,r_en,orbn,dorbn,da_orbn,iwforb)

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

    use qmckl_data

    implicit none

    integer :: ic, iel, ier, ii, iq
    integer :: iorb, k, m, m0, nxquad, iwforb
    integer :: nadorb_sav

    real(dp), dimension(3,*) :: xquad
    real(dp), dimension(nquad*nelec*2, ncent_tot) :: r_en
    real(dp), dimension(3,nquad*nelec*2, ncent_tot) :: rvec_en
    real(dp), dimension(norb_tot, *) :: orbn
    real(dp), dimension(norb_tot, nquad*nelec*2, 3) :: dorbn
    real(dp), dimension(3,ncent_tot, norb_tot, *) :: da_orbn
    real(dp), dimension(:,:,:,:),allocatable :: da_orbn_temp
    real(dp), dimension(3) :: dtmp
    real(dp) :: ddtmp

    real(dp), allocatable :: mo_qmckl(:,:)
    integer :: rc
    real(dp), allocatable :: ao_qmckl(:,:,:)
    real(dp), allocatable :: ao_vgl_qmckl(:,:,:)
    integer*8 :: n8, na8, i_image, ivgl, i_basis
    character*(1024) :: err_message = ''
    real(dp), allocatable :: xqmckl(:,:)
    real(dp), allocatable :: xqmckl_i(:,:)
    real(dp) :: rnorm
    real(dp), dimension(3) :: r_image

    nadorb_sav=nadorb

    if(ioptorb.eq.0.or.(method(1:3).ne.'lin'.and.i_sr_rescale.eq.0)) nadorb=0


    ! Send electron coordinates to QMCkl to compute the MOs at these positions
    rc = qmckl_set_point(qmckl_ctx(1), 'N', nxquad*1_8, xquad, nxquad*3_8)
    if (rc /= QMCKL_SUCCESS) then
        print *, 'orbitals quad Error setting electron coordinates in QMCkl'
        print *, "nxquad", nxquad
        stop
    end if

    rc = qmckl_get_mo_basis_mo_num(qmckl_ctx(1), n8)
    if (rc /= QMCKL_SUCCESS) then
        print *, 'orbitals quad Error getting mo_num from QMCkl'
        print *, "n8", n8
        stop
    end if

    allocate(mo_qmckl(n8, nxquad))

    ! Compute the MOs
    rc = qmckl_get_mo_basis_mo_value_inplace(qmckl_ctx(1), mo_qmckl, nxquad*n8)

    if (rc /= QMCKL_SUCCESS) then
        print *, 'Error orbitals quad getting MOs from QMCkl'
        stop
    end if

    orbn(1:norb+nadorb,1:nxquad) = mo_qmckl(1:norb+nadorb,1:nxquad)

    deallocate(mo_qmckl)

    ! To fix - QMCkl does not give da_orbitals
    if(iforce_analy.gt.0) then
        allocate(da_orbn_temp(norb,3,nxquad,ncent))  

        rc = qmckl_get_forces_mo_value_inplace(qmckl_ctx(1), da_orbn_temp, nxquad*norb*3_8*ncent)
        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl MO forces.')
        do iq=1,nxquad

            do iorb=1,norb
                do ic=1,ncent
                    do k=1,3
                        da_orbn(k,ic,iorb,iq)=da_orbn_temp(iorb,k,iq,ic)
                    enddo
                enddo

                do k=1,3
                    dorbn(iorb,iq,k)=0.d0
                enddo
                do ic=1,ncent
                    do k=1,3
                        dorbn(iorb,iq,k)=dorbn(iorb,iq,k)-da_orbn(k,ic,iorb,iq)
                    enddo
                enddo
            enddo
        enddo
        ! enddo nxquad
        deallocate(da_orbn_temp)

    endif
    ! endif iforce

    nadorb = nadorb_sav
return
end

end module
