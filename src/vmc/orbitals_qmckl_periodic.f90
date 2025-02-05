module orbitals_qmckl_periodic_mod

contains

subroutine orbitals_qmckl_periodic(x,rvec_en,r_en)

    use coefs, only: nbasis
    use find_pimage, only: find_image_pbc
    use multiple_geo, only: iwf
    use orbval, only: ddorb, dorb, nadorb, orb
    use precision_kinds, only: dp
    use phifun, only: phin, dphin, d2phin, n0_ibasis, n0_nbasis
    use periodic, only : n_images, ell
    use slater, only: norb, coef
    use system, only: ncent_tot, nelec
    use vmc_mod, only: nwftypeorb

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

    integer :: i, iorb, k, m
    integer :: m0, j

    real(dp), dimension(3,*) :: x
    real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
    real(dp), dimension(nelec,ncent_tot) :: r_en
    real(dp), dimension(:), allocatable :: auxorb !(norb+nadorb)
    real(dp), dimension(:, :), allocatable :: auxdorb !(norb+nadorb)
    real(dp), dimension(:), allocatable :: auxddorb !(norb+nadorb)
    if (.not. allocated(auxorb)) allocate (auxorb(norb+nadorb))
    if (.not. allocated(auxdorb)) allocate (auxdorb(norb+nadorb,3))
    if (.not. allocated(auxddorb)) allocate (auxddorb(norb+nadorb))

    ! get number of atomic orbitals
    rc = qmckl_get_ao_basis_ao_num(qmckl_ctx(qmckl_no_ctx-1), na8)
    if (rc /= QMCKL_SUCCESS) then
        print *, 'Error getting mo_num from QMCkl'
        stop
    end if

    if (nbasis.ne.na8) then
        print *, 'Error getting ao_num from QMCkl'
        stop
    end if

    ! first image zero

    ! allocate ao_vlg array
    allocate(ao_qmckl(nbasis, 5, nelec))
    ao_qmckl=0.d0

    ! set initial coordinates for electrons zero image
    xelec=x(1:3,1:nelec)

    !setting pbc conditions
    do i_elec=1, nelec
    call find_image_pbc(xelec(1:3,i_elec),rnorm)
    enddo

    ! Send electron coordinates to QMCkl to compute the MOs at these positions
    rc = qmckl_set_point(qmckl_ctx(qmckl_no_ctx-1), 'N', nelec*1_8, xelec, nelec*3_8)
    if (rc /= QMCKL_SUCCESS) then
        print *, 'Error setting electron coordinates in QMCkl'
        stop
    end if

    ! computing aos zero image
    rc = qmckl_get_ao_basis_ao_vgl_inplace(qmckl_ctx(qmckl_no_ctx-1), ao_qmckl, nbasis*5_8*nelec)
    if (rc /= QMCKL_SUCCESS) then
        print *, 'Error getting AOs from QMCkl zero image'
    endif

    ! computing images distance for iel
    if(n_images.gt.0) then
    ! allocate ao_vgl array for all images
        allocate(ao_vgl_qmckl(nbasis, 5, nelec))

        do i_image=1, n_images
            ! set electron images distances
            xqmckl=0.d0
            r_image=ell(1:3,i_image)
            do i_elec=1, nelec
                xqmckl(1:3,i_elec)=xelec(1:3,i_elec)-r_image(1:3)
            enddo

            ! send electron images coordinates
            rc = qmckl_set_point(qmckl_ctx(qmckl_no_ctx-1), 'N', 1_8*nelec, xqmckl, 3_8*nelec)
            if (rc /= QMCKL_SUCCESS) then
                print *, 'Error setting electron coords orbitalse'
                call qmckl_last_error(qmckl_ctx(qmckl_no_ctx-1),err_message)
                print *, trim(err_message)
                call abort()
            end if

            ao_vgl_qmckl=0.d0
            ! computing aos for the given image
            rc = qmckl_get_ao_basis_ao_vgl_inplace(qmckl_ctx(qmckl_no_ctx-1), ao_vgl_qmckl, nbasis*5_8*nelec)
            if (rc /= QMCKL_SUCCESS) then
                print *, 'Error getting AOs from QMCkl zero image'
            endif

            ! adding contribution of the given image
            do i_elec=1, nelec
                do ivgl=1, 5
                    do i_basis=1, nbasis
                        ao_qmckl(i_basis,ivgl, i_elec)=ao_qmckl(i_basis,ivgl,i_elec)+ao_vgl_qmckl(i_basis,ivgl,i_elec)
                    enddo
                enddo
            enddo

        enddo
        ! enddo loop images

    endif
    ! endif images

    ! pass QMCkl ao's to champ
    do i_elec=1, nelec
        do i_basis=1, nbasis
            phin(i_basis,i_elec)=ao_qmckl(i_basis,1,i_elec)
            dphin(i_basis,i_elec,1)=ao_qmckl(i_basis,2,i_elec)
            dphin(i_basis,i_elec,2)=ao_qmckl(i_basis,3,i_elec)
            dphin(i_basis,i_elec,3)=ao_qmckl(i_basis,4,i_elec)
            d2phin(i_basis,i_elec)=ao_qmckl(i_basis,5,i_elec)
        enddo
    enddo

    ! deallocate QMCkl champ arrays
    if(allocated(ao_qmckl)) deallocate(ao_qmckl)
    if(allocated(ao_vgl_qmckl)) deallocate(ao_vgl_qmckl)

! Vectorization dependent code selection
#ifdef VECTORIZATION

    ! Following loop changed for better vectorization AVX512/AVX2
    if(nwftypeorb.gt.1) then
        do k=1,nwftypeorb
            do i=1,nelec
                auxorb=0.d0
                auxdorb=0.d0
                auxddorb=0.d0
                do iorb=1,norb+nadorb
                    orb(i,iorb,k)=0.d0
                    dorb(iorb,i,1,k)=0.d0
                    dorb(iorb,i,2,k)=0.d0
                    dorb(iorb,i,3,k)=0.d0
                    ddorb(iorb,i,k)=0.d0
                    do m=1,nbasis
!     orb  (  i,iorb,k)=orb  (  i,iorb,k)+coef(m,iorb,k)*phin  ( m,i)
!     dorb (iorb,i,1,k)=dorb (iorb,i,1,k)+coef(m,iorb,k)*dphin (m,i,1)
!     dorb (iorb,i,2,k)=dorb (iorb,i,2,k)+coef(m,iorb,k)*dphin (m,i,2)
!     dorb (iorb,i,3,k)=dorb (iorb,i,3,k)+coef(m,iorb,k)*dphin (m,i,3)
!     ddorb(  iorb,i,k)=ddorb(iorb,i,k)+coef(m,iorb,k)*d2phin( m,i)
                        auxorb  (iorb)=auxorb  (iorb)+coef(m,iorb,k)*phin  ( m,i)
                        auxdorb (iorb,1)=auxdorb (iorb,1)+coef(m,iorb,k)*dphin (m,i,1)
                        auxdorb (iorb,2)=auxdorb (iorb,2)+coef(m,iorb,k)*dphin (m,i,2)
                        auxdorb (iorb,3)=auxdorb (iorb,3)+coef(m,iorb,k)*dphin (m,i,3)
                        auxddorb(  iorb)=auxddorb(iorb)+coef(m,iorb,k)*d2phin( m,i)
                    enddo
                enddo
                orb(i,1:(norb+nadorb),k)=auxorb(1:(norb+nadorb))
                dorb(1:(norb+nadorb),i,1:3,k)=auxdorb(1:(norb+nadorb),1:3)
                ddorb(1:(norb+nadorb),i,k)=auxddorb(1:(norb+nadorb))
            enddo
        enddo

    else

        do i=1,nelec
            auxorb=0.d0
            auxdorb=0.d0
            auxddorb=0.d0
            do iorb=1,norb+nadorb
                orb(i,iorb,1)=0.d0
                dorb(iorb,i,1,1)=0.d0
                dorb(iorb,i,2,1)=0.d0
                dorb(iorb,i,3,1)=0.d0
                ddorb(iorb,i,1)=0.d0
                do m=1,nbasis
!     orb  (  i,iorb,1)=orb  (  i,iorb,1)+coef(m,iorb,iwf)*phin  ( m,i)
!     dorb (iorb,i,1,1)=dorb (iorb,i,1,1)+coef(m,iorb,iwf)*dphin (m,i,1)
!     dorb (iorb,i,2,1)=dorb (iorb,i,2,1)+coef(m,iorb,iwf)*dphin (m,i,2)
!     dorb (iorb,i,3,1)=dorb (iorb,i,3,1)+coef(m,iorb,iwf)*dphin (m,i,3)
!     ddorb(  iorb,i,1)=ddorb(iorb,i,1)+coef(m,iorb,iwf)*d2phin( m,i)
                    auxorb  (iorb)=auxorb  (iorb)+coef(m,iorb,iwf)*phin  ( m,i)
                    auxdorb (iorb,1)=auxdorb (iorb,1)+coef(m,iorb,iwf)*dphin (m,i,1)
                    auxdorb (iorb,2)=auxdorb (iorb,2)+coef(m,iorb,iwf)*dphin (m,i,2)
                    auxdorb (iorb,3)=auxdorb (iorb,3)+coef(m,iorb,iwf)*dphin (m,i,3)
                    auxddorb(  iorb)=auxddorb(iorb)+coef(m,iorb,iwf)*d2phin( m,i)
                enddo
            enddo
            orb(i,1:(norb+nadorb),1)=auxorb(1:(norb+nadorb))
            dorb(1:(norb+nadorb),i,1:3,1)=auxdorb(1:(norb+nadorb),1:3)
            ddorb(1:(norb+nadorb),i,1)=auxddorb(1:(norb+nadorb))
        enddo

    endif
    ! nwftype endif

#else
! keep the old localization code if no vectorization instructions available

    if(nwftypeorb.gt.1) then

        do k=1,nwftypeorb
            do i=1,nelec
                auxorb=0.d0
                auxdorb=0.d0
                auxddorb=0.d0
                do iorb=1,norb+nadorb
                    orb(i,iorb,k)=0.d0
                    dorb(iorb,i,1,k)=0.d0
                    dorb(iorb,i,2,k)=0.d0
                    dorb(iorb,i,3,k)=0.d0
                    ddorb(iorb,i,k)=0.d0
                    do m0=1,n0_nbasis(i)
                        m=n0_ibasis(m0,i)
!     orb  (  i,iorb,k)=orb  (  i,iorb,k)+coef(m,iorb,k)*phin  ( m,i)
!     dorb (iorb,i,1,k)=dorb (iorb,i,1,k)+coef(m,iorb,k)*dphin (m,i,1)
!     dorb (iorb,i,2,k)=dorb (iorb,i,2,k)+coef(m,iorb,k)*dphin (m,i,2)
!     dorb (iorb,i,3,k)=dorb (iorb,i,3,k)+coef(m,iorb,k)*dphin (m,i,3)
!     ddorb(iorb,i,k)=ddorb(iorb,i,k)+coef(m,iorb,k)*d2phin( m,i)
                        auxorb  (iorb)=auxorb  (iorb)+coef(m,iorb,k)*phin  ( m,i)
                        auxdorb (iorb,1)=auxdorb (iorb,1)+coef(m,iorb,k)*dphin (m,i,1)
                        auxdorb (iorb,2)=auxdorb (iorb,2)+coef(m,iorb,k)*dphin (m,i,2)
                        auxdorb (iorb,3)=auxdorb (iorb,3)+coef(m,iorb,k)*dphin (m,i,3)
                        auxddorb(  iorb)=auxddorb(iorb)+coef(m,iorb,k)*d2phin( m,i)
                    enddo
                enddo
                orb(i,1:(norb+nadorb),k)=auxorb(1:(norb+nadorb))
                dorb(1:(norb+nadorb),i,1:3,k)=auxdorb(1:(norb+nadorb),1:3)
                ddorb(1:(norb+nadorb),i,k)=auxddorb(1:(norb+nadorb))
            enddo
        enddo

    else

        do i=1,nelec
            auxorb=0.d0
            auxdorb=0.d0
            auxddorb=0.d0
            do iorb=1,norb+nadorb
                orb(i,iorb,1)=0.d0
                dorb(iorb,i,1,1)=0.d0
                dorb(iorb,i,2,1)=0.d0
                dorb(iorb,i,3,1)=0.d0
                ddorb(iorb,i,1)=0.d0
                do m0=1,n0_nbasis(i)
                    m=n0_ibasis(m0,i)
!     orb  (  i,iorb,1)=orb  (  i,iorb,1)+coef(m,iorb,iwf)*phin  ( m,i)
!     dorb (iorb,i,1,1)=dorb (iorb,i,1,1)+coef(m,iorb,iwf)*dphin (m,i,1)
!     dorb (iorb,i,2,1)=dorb (iorb,i,2,1)+coef(m,iorb,iwf)*dphin (m,i,2)
!     dorb (iorb,i,3,1)=dorb (iorb,i,3,1)+coef(m,iorb,iwf)*dphin (m,i,3)
!     ddorb(iorb,i,1)=ddorb(iorb,i,1)+coef(m,iorb,iwf)*d2phin( m,i)
                    auxorb  (iorb)=auxorb  (iorb)+coef(m,iorb,iwf)*phin  ( m,i)
                    auxdorb (iorb,1)=auxdorb (iorb,1)+coef(m,iorb,iwf)*dphin (m,i,1)
                    auxdorb (iorb,2)=auxdorb (iorb,2)+coef(m,iorb,iwf)*dphin (m,i,2)
                    auxdorb (iorb,3)=auxdorb (iorb,3)+coef(m,iorb,iwf)*dphin (m,i,3)
                    auxddorb(  iorb)=auxddorb(iorb)+coef(m,iorb,iwf)*d2phin( m,i)
                enddo
            enddo
            orb(i,1:(norb+nadorb),1)=auxorb(1:(norb+nadorb))
            dorb(1:(norb+nadorb),i,1:3,1)=auxdorb(1:(norb+nadorb),1:3)
            ddorb(1:(norb+nadorb),i,1)=auxddorb(1:(norb+nadorb))
        enddo

    endif
    ! nwftype endif

#endif
! vectorization endif
return
end

subroutine orbitalse_qmckl_periodic(iel,x,rvec_en,r_en,iflag)

    use coefs, only: nbasis
    use find_pimage, only: find_image_pbc
    use multiple_geo, only: iwf
    use multislatern, only: ddorbn, dorbn, orbn
    use periodic, only : n_images, ell
    use precision_kinds, only: dp
    use phifun,  only: d2phin, dphin, phin
    use slater, only: norb, coef
    use system, only: ncent_tot, nelec
    use vmc_mod, only: nwftypeorb

    use qmckl_data

    implicit none

    integer :: iel, iflag, iorb, m
    integer :: k, j, ictx

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
    real(dp), dimension(3,n_images) :: xqmckl
    real(dp) :: rnorm

    if (iflag .eq. 0) then
        ictx = 2
    else
        ictx = 3
    end if

    ! get number of atomic orbitals
    rc = qmckl_get_ao_basis_ao_num(qmckl_ctx(ictx), na8)
    if (rc /= QMCKL_SUCCESS) then
       print *, 'Error getting mo_num from QMCkl'
       stop
    end if


    if (nbasis.ne.na8) then
       print *, 'Error getting ao_num from QMCkl'
       stop
    end if

    ! allocate ao_vlg array
    allocate(ao_qmckl(nbasis, 5))
    ao_qmckl=0.d0

    xiel=x(1:3,iel)
    ! applying pbc (wrapping inside the box)
    call find_image_pbc(xiel,rnorm)

    ! set one electron coordinates
    rc = qmckl_set_point(qmckl_ctx(ictx), 'N', 1_8, xiel, 3_8)
    if (rc /= QMCKL_SUCCESS) then
       print *, 'Error setting electron coords orbitalse'
       call qmckl_last_error(qmckl_ctx(ictx),err_message)
       print *, trim(err_message)
       call abort()
    end if

    ! computing aos zero image
    rc = qmckl_get_ao_basis_ao_vgl_inplace(qmckl_ctx(ictx), ao_qmckl, nbasis*5_8)
    if (rc /= QMCKL_SUCCESS) then
        print *, 'Error getting AOs from QMCkl zero image'
    endif

    ! intialize before computation just in case garbage appears
    xqmckl=0.d0
    ! computing images distance for iel
    if(n_images.gt.0) then

        ! allocate ao_vgl array for all images
        allocate(ao_vgl_qmckl(nbasis, 5, n_images))
        ao_vgl_qmckl=0.d0

        do i_image=1, n_images
            xqmckl(1:3,i_image)=xiel(1:3)-ell(1:3,i_image)
        enddo

        ! send coordinates
        rc = qmckl_set_point(qmckl_ctx(ictx), 'N', 1_8*n_images, xqmckl, 3_8*n_images)
        if (rc /= QMCKL_SUCCESS) then
            print *, 'Error setting electron coords orbitalse'
            call qmckl_last_error(qmckl_ctx(ictx),err_message)
            print *, trim(err_message)
            call abort()
        end if

        rc = qmckl_get_ao_basis_ao_vgl_inplace(qmckl_ctx(ictx), ao_vgl_qmckl, nbasis*n_images*5_8)
        if (rc /= QMCKL_SUCCESS) then
            print *, 'Error getting AOs from QMCkl'
        endif


        ! summing it all over the image0
        do i_image=1, n_images
            do ivgl=1, 5
                do i_basis=1, nbasis
                    ao_qmckl(i_basis,ivgl)=ao_qmckl(i_basis,ivgl)+ao_vgl_qmckl(i_basis,ivgl,i_image)
                enddo
            enddo
        enddo

    endif

    ! copying QMCkl ao's back to champ
    do i_basis=1, nbasis
        phin(i_basis,iel)=ao_qmckl(i_basis,1)
        dphin(i_basis,iel,1)=ao_qmckl(i_basis,2)
        dphin(i_basis,iel,2)=ao_qmckl(i_basis,3)
        dphin(i_basis,iel,3)=ao_qmckl(i_basis,4)
        d2phin(i_basis,iel)=ao_qmckl(i_basis,5)
    enddo

    if(allocated(ao_qmckl)) deallocate(ao_qmckl)
    if(allocated(ao_vgl_qmckl)) deallocate(ao_vgl_qmckl)

    if(iflag.gt.0) then

        if(nwftypeorb.gt.1) then

            do k=1,nwftypeorb
                do iorb=1,norb
                    orbn(iorb,k)=0.d0
                    dorbn(iorb,1,k)=0.d0
                    dorbn(iorb,2,k)=0.d0
                    dorbn(iorb,3,k)=0.d0
                    ddorbn(iorb,k)=0.d0
                    do m=1,nbasis
                        orbn(iorb,k)=orbn(iorb,k)+coef(m,iorb,k)*phin(m,iel)
                        dorbn(iorb,1,k)=dorbn(iorb,1,k)+coef(m,iorb,k)*dphin(m,iel,1)
                        dorbn(iorb,2,k)=dorbn(iorb,2,k)+coef(m,iorb,k)*dphin(m,iel,2)
                        dorbn(iorb,3,k)=dorbn(iorb,3,k)+coef(m,iorb,k)*dphin(m,iel,3)
                        ddorbn(iorb,k)=ddorbn(iorb,k)+coef(m,iorb,k)*d2phin(m,iel)
                    enddo
                enddo
            enddo

        else

            do iorb=1,norb
                orbn(iorb,1)=0.d0
                dorbn(iorb,1,1)=0.d0
                dorbn(iorb,2,1)=0.d0
                dorbn(iorb,3,1)=0.d0
                ddorbn(iorb,1)=0.d0
                do m=1,nbasis
                    orbn(iorb,1)=orbn(iorb,1)+coef(m,iorb,iwf)*phin(m,iel)
                    dorbn(iorb,1,1)=dorbn(iorb,1,1)+coef(m,iorb,iwf)*dphin(m,iel,1)
                    dorbn(iorb,2,1)=dorbn(iorb,2,1)+coef(m,iorb,iwf)*dphin(m,iel,2)
                    dorbn(iorb,3,1)=dorbn(iorb,3,1)+coef(m,iorb,iwf)*dphin(m,iel,3)
                    ddorbn(iorb,1)=ddorbn(iorb,1)+coef(m,iorb,iwf)*d2phin(m,iel)
                enddo
            enddo


        endif
        ! endif nwftype
    else
    ! else iflag

        if(nwftypeorb.gt.1) then

            do k=1,nwftypeorb
                do iorb=1,norb
                    orbn(iorb,k)=0.d0
                    dorbn(iorb,1,k)=0.d0
                    dorbn(iorb,2,k)=0.d0
                    dorbn(iorb,3,k)=0.d0
                    do m=1,nbasis
                    orbn(iorb,k)=orbn(iorb,k)+coef(m,iorb,k)*phin(m,iel)
                    dorbn(iorb,1,k)=dorbn(iorb,1,k)+coef(m,iorb,k)*dphin(m,iel,1)
                    dorbn(iorb,2,k)=dorbn(iorb,2,k)+coef(m,iorb,k)*dphin(m,iel,2)
                    dorbn(iorb,3,k)=dorbn(iorb,3,k)+coef(m,iorb,k)*dphin(m,iel,3)
                    enddo
                enddo
            enddo

        else

            do iorb=1,norb
                orbn(iorb,1)=0.d0
                dorbn(iorb,1,1)=0.d0
                dorbn(iorb,2,1)=0.d0
                dorbn(iorb,3,1)=0.d0
                do m=1,nbasis
                    orbn(iorb,1)=orbn(iorb,1)+coef(m,iorb,iwf)*phin(m,iel)
                    dorbn(iorb,1,1)=dorbn(iorb,1,1)+coef(m,iorb,iwf)*dphin(m,iel,1)
                    dorbn(iorb,2,1)=dorbn(iorb,2,1)+coef(m,iorb,iwf)*dphin(m,iel,2)
                    dorbn(iorb,3,1)=dorbn(iorb,3,1)+coef(m,iorb,iwf)*dphin(m,iel,3)
                enddo
            enddo

        endif
        ! endif nwftype

    endif
    ! endif iflag

return
end

subroutine orbitals_quad_qmckl_periodic(nxquad,xquad,rvec_en,r_en,orbn,dorbn,da_orbn,iwforb)

    use coefs, only: nbasis
    use find_pimage, only: find_image_pbc
    use m_force_analytic, only: iforce_analy
    use multiple_geo, only: iwf
    use numbas2, only: ibas0,ibas1
    use optwf_control, only: ioptorb
    use optwf_control, only: method
    use orbval,  only: nadorb
    use phifun,  only: dphin,n0_ibasis,n0_ic,n0_nbasis,phin
    use periodic, only : n_images, ell
    use precision_kinds, only: dp
    use qua,     only: nquad
    use slater,  only: coef,norb
    use sr_mod,  only: i_sr_rescale
    use system,  only: ncent,ncent_tot,nelec
    use vmc_mod, only: norb_tot, nwftypeorb

    use qmckl_data

    implicit none

    integer :: ic, iq
    integer :: iorb, k, m, m0, nxquad, iwforb
    integer :: nadorb_sav

    real(dp), dimension(3,*) :: xquad
    real(dp), dimension(nquad*nelec*2, ncent_tot) :: r_en
    real(dp), dimension(3,nquad*nelec*2, ncent_tot) :: rvec_en
    real(dp), dimension(norb_tot, *) :: orbn
    real(dp), dimension(norb_tot, nquad*nelec*2, 3) :: dorbn
    real(dp), dimension(3,ncent_tot, norb_tot, *) :: da_orbn
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

    if(nwftypeorb.gt.1) iwf=1

    ! get number of atomic orbitals
    rc = qmckl_get_ao_basis_ao_num(qmckl_ctx(1), na8)
    if (rc /= QMCKL_SUCCESS) then
        print *, 'Error getting mo_num from QMCkl'
        stop
    end if
    
    
    if (nbasis.ne.na8) then
        print *, 'Error getting ao_num from QMCkl'
        stop
    end if
    
    ! Image zero calculation
    allocate(ao_qmckl(nbasis, 5, nxquad))
    ao_qmckl=0.d0

    allocate(xqmckl(3,nxquad))
    xqmckl=xquad(1:3,1:nxquad)
    ! apply pbc (wraping inside the box)
    do iq=1,nxquad
        call find_image_pbc(xqmckl(1:3,iq),rnorm)
    enddo
    
    ! Send electron coordinates to QMCkl to compute the MOs at these positions
    rc = qmckl_set_point(qmckl_ctx(1), 'N', nxquad*1_8, xqmckl, nxquad*3_8)
    if (rc /= QMCKL_SUCCESS) then
        print *, 'Error setting electron coordinates QMCkl orbitals_quad'
        stop
    end if
    
    ! computing ao's zero image
    rc = qmckl_get_ao_basis_ao_vgl_inplace(qmckl_ctx(1), ao_qmckl, nxquad*5_8*nbasis)
    if (rc /= QMCKL_SUCCESS) then
        print *, 'Error getting AOs from QMCkl zero image'
    endif
    
    ! computing images distance for nxquad points
    if(n_images.gt.0) then
    
        ! allocate ao_vgl array for all images
        allocate(ao_vgl_qmckl(nbasis, 5, nxquad))
        allocate(xqmckl_i(3,nxquad))

        do i_image=1, n_images
            ! initilialize xqmckl_i
            ao_vgl_qmckl=0.d0
            xqmckl_i=0.d0

            r_image=ell(1:3,i_image)
            do iq=1, nxquad
            xqmckl_i(1:3,iq)=xqmckl(1:3,iq)-r_image(1:3)
            enddo
    
            ! send coordinates of xquad image
            rc = qmckl_set_point(qmckl_ctx(1), 'N', 1_8*nxquad, xqmckl_i, 3_8*nxquad)
            if (rc /= QMCKL_SUCCESS) then
                print *, 'Error electron coords orbitals quad'
                call qmckl_last_error(qmckl_ctx(1),err_message)
                print *, trim(err_message)
                call abort()
            end if
    
            ! computing aos for the given image
            rc = qmckl_get_ao_basis_ao_vgl_inplace(qmckl_ctx(1),ao_vgl_qmckl, nxquad*5_8*nbasis)
            if (rc /= QMCKL_SUCCESS) then
                print *, 'Error getting AOs from QMCkl zero image'
            endif
    
    
            ! add contribution of the given image
            do iq=1, nxquad
                do ivgl=1, 5
                    do i_basis=1, nbasis
                        ao_qmckl(i_basis,ivgl, iq)=ao_qmckl(i_basis,ivgl,iq)+ao_vgl_qmckl(i_basis,ivgl,iq)
                    enddo
                enddo
            enddo

        enddo
        !enddo images    
    
    endif
    ! endif periodic images
    
    ! passing the qmckl ao's back to champ
    do iq=1, nxquad
        do i_basis=1, nbasis
            phin(i_basis,iq)=ao_qmckl(i_basis,1,iq)
            dphin(i_basis,iq,1)=ao_qmckl(i_basis,2,iq)
            dphin(i_basis,iq,2)=ao_qmckl(i_basis,3,iq)
            dphin(i_basis,iq,3)=ao_qmckl(i_basis,4,iq)
        enddo
    enddo

    if(allocated(ao_qmckl)) deallocate(ao_qmckl)
    if(allocated(ao_vgl_qmckl)) deallocate(ao_vgl_qmckl)
    if(allocated(xqmckl)) deallocate(xqmckl)
    if(allocated(xqmckl_i)) deallocate(xqmckl_i)

    if(nwftypeorb.gt.1) iwf=iwforb

    do iq=1,nxquad
    
! Vectorization dependent code selection
#ifdef VECTORIZATION
    ! The following loop changed for better vectorization AVX512/AVX2
        do iorb=1,norb+nadorb
            orbn(iorb,iq)=0.d0
            do m=1,nbasis
                orbn(iorb,iq)=orbn(iorb,iq)+coef(m,iorb,iwf)*phin(m,iq)
            enddo
        enddo
#else
        do iorb=1,norb+nadorb
            orbn(iorb,iq)=0.d0
            do m0=1,n0_nbasis(iq)
                m=n0_ibasis(m0,iq)
                orbn(iorb,iq)=orbn(iorb,iq)+coef(m,iorb,iwf)*phin(m,iq)
            enddo
        enddo
#endif
    
        if(iforce_analy.gt.0) then
            do iorb=1,norb
                do ic=1,ncent
                    do k=1,3
                        da_orbn(k,ic,iorb,iq)=0.d0
                    enddo
                enddo
#ifdef VECTORIZATION
                do ic=1,ncent
                    do k=1,3
                        do m=ibas0(ic),ibas1(ic)
                            da_orbn(k,ic,iorb,iq)=da_orbn(k,ic,iorb,iq)-coef(m,iorb,iwf)*dphin(m,iq,k)
                        enddo
                    enddo
                enddo
#else
                do m0=1,n0_nbasis(iq)
                    m=n0_ibasis(m0,iq)
                    ic=n0_ic(m0,iq)
                    do k=1,3
                        da_orbn(k,ic,iorb,iq)=da_orbn(k,ic,iorb,iq)-coef(m,iorb,iwf)*dphin(m,iq,k)
                    enddo
                enddo
#endif
                do k=1,3
                    dorbn(iorb,iq,k)=0.d0
                enddo
                do ic=1,ncent
                    do k=1,3
                        dorbn(iorb,iq,k)=dorbn(iorb,iq,k)-da_orbn(k,ic,iorb,iq)
                    enddo
                enddo
            enddo
    
        endif

    enddo
    ! enddo nxquad

    nadorb = nadorb_sav

return
end

end module