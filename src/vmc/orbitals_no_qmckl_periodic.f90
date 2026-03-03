module orbitals_no_qmckl_periodic_mod
    interface !LAPACK interface
    SUBROUTINE dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
! *  -- Reference BLAS level3 routine --
! *  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
! *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
    END SUBROUTINE
    SUBROUTINE dcopy(N,DX,INCX,DY,INCY)
      INTEGER INCX,INCY,N
      DOUBLE PRECISION DX(*),DY(*)
    END SUBROUTINE
  end interface
contains

subroutine orbitals_no_qmckl_periodic(x,rvec_en,r_en)

    use basis_fns_mod, only: basis_fns 
    use coefs, only: nbasis
    use multiple_geo, only: iwf !Current wavefunction index?
    use m_force_analytic, only: iforce_analy !input whether to calculate the hessian matrix
    use orbval, only: ddorb, dorb, nadorb, orb !Allocates arrays for the molecular orbitals and it's derivative
    use phifun, only: phin, dphin, d2phin, n0_ibasis, n0_nbasis !array used for the atomic orbitals
    use precision_kinds, only: dp 
    use slater, only: norb, coef  !allocates arrays for slater matrices
    use system, only: ncent_tot, nelec 
    use vmc_mod, only: nwftypeorb ! Number of orbital wavefunction types
    use periodic, only : n_images, ell
    implicit none
        
    integer :: i, ider, iorb, k, m
    integer :: m0, j
    integer :: i_image

    real(dp), dimension(3,*) :: x
    real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
    real(dp), dimension(nelec,ncent_tot) :: r_en
    real(dp), dimension(3, nelec, ncent_tot) :: rvec_en_image
    real(dp), dimension(nelec, ncent_tot) :: r_en_image
    real(dp), dimension(:), allocatable :: auxorb !(norb+nadorb)
    real(dp), dimension(:, :), allocatable :: auxdorb !(norb+nadorb)
    real(dp), dimension(:), allocatable :: auxddorb !(norb+nadorb)
    if (.not. allocated(auxorb)) allocate (auxorb(norb+nadorb))
    if (.not. allocated(auxdorb)) allocate (auxdorb(norb+nadorb,3))
    if (.not. allocated(auxddorb)) allocate (auxddorb(norb+nadorb))
    ! get basis functions for all electrons
    ider=2
    if(iforce_analy.eq.1) ider=3
    !Calculate rvec_en_image before the loop to avoid redundant calculations. This is used in the basis_fns subroutine to calculate the basis functions for each image.

    ! do i_image=1,n_images
    !     rvec_en_image(i_image, 1, :, :) = rvec_en(1, :, :) - ell(1, i_image)
    !     rvec_en_image(i_image, 2, :, :) = rvec_en(2, :, :) - ell(2, i_image)
    !     rvec_en_image(i_image, 3, :, :) = rvec_en(3, :, :) - ell(3, i_image)
    !     r_en_image(i_image, :, :) = sqrt(sum(rvec_en_image(i_image, :, :, :)**2, dim=2))
    ! enddo

    
    orb(1:nelec,1:(norb+nadorb),:)=0.d0
    dorb(1:(norb+nadorb),1:nelec,:,:)=0.d0
    ddorb(1:(norb+nadorb),1:nelec,:)=0.d0
#ifdef VECTORIZATION
    if(nwftypeorb.gt.1) then
        do i_image=1,n_images
            rvec_en_image(1, :, :) = rvec_en(1, :, :) - ell(1,i_image)
            rvec_en_image(2, :, :) = rvec_en(2, :, :) - ell(2,i_image)
            rvec_en_image(3, :, :) = rvec_en(3, :, :) - ell(3,i_image)
            r_en_image(:,:) = sqrt(rvec_en_image(1,:,:)**2 + rvec_en_image(2,:,:)**2 + rvec_en_image(3,:,:)**2)
            call basis_fns(1,nelec,nelec,rvec_en_image,r_en_image,ider)
            do k=1,nwftypeorb
                do i=1,nelec
                    auxorb=0.d0
                    auxdorb=0.d0
                    auxddorb=0.d0
                    do iorb=1,norb+nadorb
                        do m=1,nbasis
                            auxorb(iorb)    =auxorb(iorb)+coef(m,iorb,k)*phin(m,i)
                            auxdorb(iorb,1) =auxdorb(iorb,1)+coef(m,iorb,k)*dphin(m,i,1)
                            auxdorb(iorb,2) =auxdorb(iorb,2)+coef(m,iorb,k)*dphin(m,i,2)
                            auxdorb(iorb,3) =auxdorb(iorb,3)+coef(m,iorb,k)*dphin(m,i,3)
                            auxddorb(iorb)  =auxddorb(iorb)+coef(m,iorb,k)*d2phin(m,i)
                        enddo ! loop over basis functions(atomic orbitals)
                    enddo ! loop over molecular orbitals     
                    orb(i,1:(norb+nadorb),k) = orb(i,1:(norb+nadorb),k) + auxorb(1:(norb+nadorb))
                    dorb(1:(norb+nadorb),i,1:3,k) = dorb(1:(norb+nadorb),i,1:3,k) + auxdorb(1:(norb+nadorb),1:3)
                    ddorb(1:(norb+nadorb),i,k) = ddorb(1:(norb+nadorb),i,k) + auxddorb(1:(norb+nadorb))
                enddo! loop over electrons
            enddo ! loop over wavefunction types
        enddo !images loop
    else
        do i_image=1,n_images
            rvec_en_image(1, :, :) = rvec_en(1, :, :) - ell(1,i_image)
            rvec_en_image(2, :, :) = rvec_en(2, :, :) - ell(2,i_image)
            rvec_en_image(3, :, :) = rvec_en(3, :, :) - ell(3,i_image)
            r_en_image(:,:) = sqrt(rvec_en_image(1,:,:)**2 + rvec_en_image(2,:,:)**2 + rvec_en_image(3,:,:)**2)
            call basis_fns(1,nelec,nelec,rvec_en_image,r_en_image,ider)
            do i=1,nelec
                auxorb=0.d0
                auxdorb=0.d0
                auxddorb=0.d0
                do iorb=1,norb+nadorb
                    do m=1,nbasis
                        auxorb(iorb)    =auxorb(iorb)+coef(m,iorb,iwf)*phin(m,i)
                        auxdorb(iorb,1) =auxdorb(iorb,1)+coef(m,iorb,iwf)*dphin(m,i,1)
                        auxdorb(iorb,2) =auxdorb(iorb,2)+coef(m,iorb,iwf)*dphin(m,i,2)
                        auxdorb(iorb,3) =auxdorb(iorb,3)+coef(m,iorb,iwf)*dphin(m,i,3)
                        auxddorb(iorb)  =auxddorb(iorb)+coef(m,iorb,iwf)*d2phin(m,i)
                    enddo ! loop over basis functions(atomic orbitals)
                enddo ! loop over molecular orbitals     
                orb(i,1:(norb+nadorb),1) = orb(i,1:(norb+nadorb),1) + auxorb(1:(norb+nadorb))
                dorb(1:(norb+nadorb),i,1:3,1) = dorb(1:(norb+nadorb),i,1:3,1) + auxdorb(1:(norb+nadorb),1:3)
                ddorb(1:(norb+nadorb),i,1) = ddorb(1:(norb+nadorb),i,1) + auxddorb(1:(norb+nadorb))
            enddo! loop over electrons
        enddo !images loop

    endif!nwftype endif
#else
    if(nwftypeorb.gt.1) then
        do image=1,n_images
            rvec_en_image(1, :, :) = rvec_en(1, :, :) - ell(1,i_image)
            rvec_en_image(2, :, :) = rvec_en(2, :, :) - ell(2,i_image)
            rvec_en_image(3, :, :) = rvec_en(3, :, :) - ell(3,i_image)
            r_en_image(:,:) = sqrt(rvec_en_image(1,:,:)**2 + rvec_en_image(2,:,:)**2 + rvec_en_image(3,:,:)**2)
            call basis_fns(1,nelec,nelec,rvec_en_image,r_en_image,ider)
            do k=1=nwftypeorb
                do i=1,nelec
                    auxorb=0.d0
                    auxdorb=0.d0
                    auxddorb=0.d0
                    do iorb=1,norb+nadorb
                        do m0=1,n0_nbasis(i)
                            m=n0_ibasis(m0,i)
                            auxorb(iorb)    =auxorb(iorb)+coef(m,iorb,k)*phin(m,i)
                            auxdorb(iorb,1) =auxdorb(iorb,1)+coef(m,iorb,k)*dphin(m,i,1)
                            auxdorb(iorb,2) =auxdorb(iorb,2)+coef(m,iorb,k)*dphin(m,i,2)
                            auxdorb(iorb,3) =auxdorb(iorb,3)+coef(m,iorb,k)*dphin(m,i,3)
                            auxddorb(iorb)  =auxddorb(iorb)+coef(m,iorb,k)*d2phin(m,i)
                        enddo ! loop over basis functions(atomic orbitals)
                    enddo ! loop over molecular orbitals
                    orb(i,1:(norb+nadorb),k) = orb(i,1:(norb+nadorb),k) + auxorb(1:(norb+nadorb))
                    dorb(1:(norb+nadorb),i,1:3,k) = dorb(1:(norb+nadorb),i,1:3,k) + auxdorb(1:(norb+nadorb),1:3)
                    ddorb(1:(norb+nadorb),i,k) = ddorb(1:(norb+nadorb),i,k) + auxddorb(1:(norb+nadorb))
                enddo ! loop over electrons
            enddo ! loop over wavefunction types
        enddo !images loop
    else
        do i_image=1,n_images
            rvec_en_image(1, :, :) = rvec_en(1, :, :) - ell(1,i_image)
            rvec_en_image(2, :, :) = rvec_en(2, :, :) - ell(2,i_image)
            rvec_en_image(3, :, :) = rvec_en(3, :, :) - ell(3,i_image)
            r_en_image(:,:) = sqrt(rvec_en_image(1,:,:)**2 + rvec_en_image(2,:,:)**2 + rvec_en_image(3,:,:)**2)
            call basis_fns(1,nelec,nelec,rvec_en_image,r_en_image,ider)
            do i=1,nelec
                auxorb=0.d0
                auxdorb=0.d0
                auxddorb=0.d0
                do iorb=1,norb+nadorb
                    do m0=1,n0_basis(i)
                        m=n0_ibasis(m0,i)
                        auxorb(iorb)    =auxorb(iorb)+coef(m,iorb,iwf)*phin(m,i)
                        auxdorb(iorb,1) =auxdorb(iorb,1)+coef(m,iorb,iwf)*dphin(m,i,1)
                        auxdorb(iorb,2) =auxdorb(iorb,2)+coef(m,iorb,iwf)*dphin(m,i,2)
                        auxdorb(iorb,3) =auxdorb(iorb,3)+coef(m,iorb,iwf)*dphin(m,i,3)
                        auxddorb(iorb)  =auxddorb(iorb)+coef(m,iorb,iwf)*d2phin(m,i)
                    enddo ! loop over basis functions(atomic orbitals)
                enddo ! loop over molecular orbitals     
                orb(i,1:(norb+nadorb),1) = orb(i,1:(norb+nadorb),1) + auxorb(1:(norb+nadorb))
                dorb(1:(norb+nadorb),i,1:3,1) = dorb(1:(norb+nadorb),i,1:3,1) + auxdorb(1:(norb+nadorb),1:3)
                ddorb(1:(norb+nadorb),i,1) = ddorb(1:(norb+nadorb),i,1) + auxddorb(1:(norb+nadorb))
            enddo! loop over electrons
        enddo !images loop

    endif !nwftype endif
#endif 
!vectorization endif
return
end


subroutine orbitalse_no_qmckl_periodic(iel,x,rvec_en,r_en,iflag)
  
    use basis_fns_mod, only: basis_fns
    use coefs, only: nbasis
    use multiple_geo, only: iwf
    use multislatern, only: ddorbn, dorbn, orbn
    use phifun, only: d2phin, dphin, n0_ibasis, n0_nbasis
    use phifun, only: phin
    use precision_kinds, only: dp
    use slater, only: norb, coef
    use system, only: ncent_tot, nelec
    use vmc_mod, only: nwftypeorb
    use periodic, only : n_images, ell

    implicit none

    integer :: iel, ider, iflag, iorb, m
    integer :: m0, k, j
    integer :: i_image

    real(dp), dimension(3,*) :: x
    real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
    real(dp), dimension(nelec,ncent_tot) :: r_en
    real(dp), dimension(3, nelec, ncent_tot) :: rvec_en_image
    real(dp), dimension(nelec, ncent_tot) :: r_en_image

    ider=1
    if(iflag.gt.0) ider=2


    orbn(1:norb,:)=0.d0
    dorbn(1:norb,:,:)=0.d0
    ddorbn(1:norb,:)=0.d0
! Vectorization dependent code. useful for AVX512 and AVX2
#ifdef VECTORIZATION
    if(iflag.gt.0) then
        if(nwftypeorb.gt.1) then
            do i_image=1,n_images
                rvec_en_image(1, iel, :) = rvec_en(1, iel, :) - ell(1,i_image)
                rvec_en_image(2, iel, :) = rvec_en(2, iel, :) - ell(2,i_image)
                rvec_en_image(3, iel, :) = rvec_en(3, iel, :) - ell(3,i_image)
                r_en_image(iel,:) = sqrt(rvec_en_image(1,iel,:)**2 + rvec_en_image(2,iel,:)**2 + rvec_en_image(3,iel,:)**2)
                call basis_fns(iel,iel,nelec,rvec_en_image,r_en_image,ider)
                do k=1,nwftypeorb
                    do iorb=1,norb
                        do m=1,nbasis
                            orbn(iorb,k)=orbn(iorb,k)+coef(m,iorb,k)*phin(m,iel)
                            dorbn(iorb,1,k)=dorbn(iorb,1,k)+coef(m,iorb,k)*dphin(m,iel,1)
                            dorbn(iorb,2,k)=dorbn(iorb,2,k)+coef(m,iorb,k)*dphin(m,iel,2)
                            dorbn(iorb,3,k)=dorbn(iorb,3,k)+coef(m,iorb,k)*dphin(m,iel,3)
                            ddorbn(iorb,k)=ddorbn(iorb,k)+coef(m,iorb,k)*d2phin(m,iel)
                        enddo ! loop over basis functions(atomic orbitals)
                    enddo ! loop over molecular orbitals
                enddo ! loop over wavefunction types
            enddo ! loop over images
        else
            do i_image=1,n_images
                rvec_en_image(1, iel, :) = rvec_en(1, iel, :) - ell(1,i_image)
                rvec_en_image(2, iel, :) = rvec_en(2, iel, :) - ell(2,i_image)
                rvec_en_image(3, iel, :) = rvec_en(3, iel, :) - ell(3,i_image)
                r_en_image(iel,:) = sqrt(rvec_en_image(1,iel,:)**2 + rvec_en_image(2,iel,:)**2 + rvec_en_image(3,iel,:)**2)
                call basis_fns(iel,iel,nelec,rvec_en_image,r_en_image,ider)
                do iorb=1,norb
                    do m=1,nbasis
                        orbn(iorb,1)=orbn(iorb,1)+coef(m,iorb,iwf)*phin(m,iel)
                        dorbn(iorb,1,1)=dorbn(iorb,1,1)+coef(m,iorb,iwf)*dphin(m,iel,1)
                        dorbn(iorb,2,1)=dorbn(iorb,2,1)+coef(m,iorb,iwf)*dphin(m,iel,2)
                        dorbn(iorb,3,1)=dorbn(iorb,3,1)+coef(m,iorb,iwf)*dphin(m,iel,3)
                        ddorbn(iorb,1)=ddorbn(iorb,1)+coef(m,iorb,iwf)*d2phin(m,iel)
                    enddo ! loop over basis functions(atomic orbitals)
                enddo ! loop over molecular orbitals
            enddo ! loop over images
        endif !endif nwftype

    else! else iflag
        if(nwftypeorb.gt.1) then
            do i_image=1,n_images
                rvec_en_image(1, iel, :) = rvec_en(1, iel, :) - ell(1,i_image)
                rvec_en_image(2, iel, :) = rvec_en(2, iel, :) - ell(2,i_image)
                rvec_en_image(3, iel, :) = rvec_en(3, iel, :) - ell(3,i_image)
                r_en_image(iel,:) = sqrt(rvec_en_image(1,iel,:)**2 + rvec_en_image(2,iel,:)**2 + rvec_en_image(3,iel,:)**2)
                call basis_fns(iel,iel,nelec,rvec_en_image,r_en_image,ider)
                do k=1,nwftypeorb
                    do iorb=1,norb
                        do m=1,nbasis
                            orbn(iorb,k)=orbn(iorb,k)+coef(m,iorb,k)*phin(m,iel)
                            dorbn(iorb,1,k)=dorbn(iorb,1,k)+coef(m,iorb,k)*dphin(m,iel,1)
                            dorbn(iorb,2,k)=dorbn(iorb,2,k)+coef(m,iorb,k)*dphin(m,iel,2)
                            dorbn(iorb,3,k)=dorbn(iorb,3,k)+coef(m,iorb,k)*dphin(m,iel,3)
                        enddo ! loop over basis functions(atomic orbitals)
                    enddo ! loop over molecular orbitals
                enddo !loop over wavefunction types
            enddo ! loop over images
        else
            do i_image=1,n_images
                rvec_en_image(1, iel, :) = rvec_en(1, iel, :) - ell(1,i_image)
                rvec_en_image(2, iel, :) = rvec_en(2, iel, :) - ell(2,i_image)
                rvec_en_image(3, iel, :) = rvec_en(3, iel, :) - ell(3,i_image)
                r_en_image(iel,:) = sqrt(rvec_en_image(1,iel,:)**2 + rvec_en_image(2,iel,:)**2 + rvec_en_image(3,iel,:)**2)
                call basis_fns(iel,iel,nelec,rvec_en_image,r_en_image,ider)
                do iorb=1,norb
                    do m=1,nbasis
                        orbn(iorb,1)=orbn(iorb,1)+coef(m,iorb,iwf)*phin(m,iel)
                        dorbn(iorb,1,1)=dorbn(iorb,1,1)+coef(m,iorb,iwf)*dphin(m,iel,1)
                        dorbn(iorb,2,1)=dorbn(iorb,2,1)+coef(m,iorb,iwf)*dphin(m,iel,2)
                        dorbn(iorb,3,1)=dorbn(iorb,3,1)+coef(m,iorb,iwf)*dphin(m,iel,3)
                    enddo ! loop over basis functions(atomic orbitals)
                enddo ! loop over molecular orbitals
            enddo ! loop over images
    
        endif ! nwftype 


    endif !iflag

#else
! Keep the localization for the non-vectorized code
       
    if(iflag.gt.0) then
        if(nwftypeorb.gt.1) then
            do i_image = 1,n_images
                rvec_en_image(1, iel, :) = rvec_en(1, iel, :) - ell(1,i_image)
                rvec_en_image(2, iel, :) = rvec_en(2, iel, :) - ell(2,i_image)
                rvec_en_image(3, iel, :) = rvec_en(3, iel, :) - ell(3,i_image)
                r_en_image(iel,:) = sqrt(rvec_en_image(1,iel,:)**2 + rvec_en_image(2,iel,:)**2 + rvec_en_image(3,iel,:)**2)
                call basis_fns(iel,iel,nelec,rvec_en_image,r_en_image,ider)
                do k=1,nwftypeorb
                    do iorb=1,norb
                        do m0=1,n0_nbasis(iel)
                            m=n0_ibasis(m0,iel)
                            orbn(iorb,k)=orbn(iorb,k)+coef(m,iorb,k)*phin(m,iel)
                            dorbn(iorb,1,k)=dorbn(iorb,1,k)+coef(m,iorb,k)*dphin(m,iel,1)
                            dorbn(iorb,2,k)=dorbn(iorb,2,k)+coef(m,iorb,k)*dphin(m,iel,2)
                            dorbn(iorb,3,k)=dorbn(iorb,3,k)+coef(m,iorb,k)*dphin(m,iel,3)
                            ddorbn(iorb,k)=ddorbn(iorb,k)+coef(m,iorb,k)*d2phin(m,iel)
                        enddo
                    enddo
                enddo
            enddo
        else
            do image=1,n_images
                rvec_en_image(1, iel, :) = rvec_en(1, iel, :) - ell(1,i_image)
                rvec_en_image(2, iel, :) = rvec_en(2, iel, :) - ell(2,i_image)
                rvec_en_image(3, iel, :) = rvec_en(3, iel, :) - ell(3,i_image)
                r_en_image(iel,:) = sqrt(rvec_en_image(1,iel,:)**2 + rvec_en_image(2,iel,:)**2 + rvec_en_image(3,iel,:)**2)
                call basis_fns(iel,iel,nelec,rvec_en_image,r_en_image,ider)
                do iorb=1,norb
                    do m0=1,n0_nbasis(iel)
                        m=n0_ibasis(m0,iel)
                        orbn(iorb,1)=orbn(iorb,1)+coef(m,iorb,iwf)*phin(m,iel)
                        dorbn(iorb,1,1)=dorbn(iorb,1,1)+coef(m,iorb,iwf)*dphin(m,iel,1)
                        dorbn(iorb,2,1)=dorbn(iorb,2,1)+coef(m,iorb,iwf)*dphin(m,iel,2)
                        dorbn(iorb,3,1)=dorbn(iorb,3,1)+coef(m,iorb,iwf)*dphin(m,iel,3)
                        ddorbn(iorb,1)=ddorbn(iorb,1)+coef(m,iorb,iwf)*d2phin(m,iel)
                    enddo
                enddo
            enddo
        endif
        ! endif nwftype

    else
    ! else iflag
        if(nwftypeorb.gt.1) then
            do i_image=1,n_images
                rvec_en_image(1, iel, :) = rvec_en(1, iel, :) - ell(1,i_image)
                rvec_en_image(2, iel, :) = rvec_en(2, iel, :) - ell(2,i_image)
                rvec_en_image(3, iel, :) = rvec_en(3, iel, :) - ell(3,i_image)
                r_en_image(iel,:) = sqrt(rvec_en_image(1,iel,:)**2 + rvec_en_image(2,iel,:)**2 + rvec_en_image(3,iel,:)**2)
                call basis_fns(iel,iel,nelec,rvec_en_image,r_en_image,ider)
                do k=1,nwftypeorb
                    do iorb=1,norb
                        do m0=1,n0_nbasis(iel)
                            m=n0_ibasis(m0,iel)
                            orbn(iorb,k)=orbn(iorb,k)+coef(m,iorb,k)*phin(m,iel)
                            dorbn(iorb,1,k)=dorbn(iorb,1,k)+coef(m,iorb,k)*dphin(m,iel,1)
                            dorbn(iorb,2,k)=dorbn(iorb,2,k)+coef(m,iorb,k)*dphin(m,iel,2)
                            dorbn(iorb,3,k)=dorbn(iorb,3,k)+coef(m,iorb,k)*dphin(m,iel,3)
                        enddo
                    enddo
                enddo
            enddo
        else
            do i_image=1,n_images
                rvec_en_image(1, iel, :) = rvec_en(1, iel, :) - ell(1,i_image)
                rvec_en_image(2, iel, :) = rvec_en(2, iel, :) - ell(2,i_image)
                rvec_en_image(3, iel, :) = rvec_en(3, iel, :) - ell(3,i_image)
                r_en_image(iel,:) = sqrt(rvec_en_image(1,iel,:)**2 + rvec_en_image(2,iel,:)**2 + rvec_en_image(3,iel,:)**2)
                call basis_fns(iel,iel,nelec,rvec_en_image,r_en_image,ider)
                do iorb=1,norb
                    do m0=1,n0_nbasis(iel)
                        m=n0_ibasis(m0,iel)
                        orbn(iorb,1)=orbn(iorb,1)+coef(m,iorb,iwf)*phin(m,iel)
                        dorbn(iorb,1,1)=dorbn(iorb,1,1)+coef(m,iorb,iwf)*dphin(m,iel,1)
                        dorbn(iorb,2,1)=dorbn(iorb,2,1)+coef(m,iorb,iwf)*dphin(m,iel,2)
                        dorbn(iorb,3,1)=dorbn(iorb,3,1)+coef(m,iorb,iwf)*dphin(m,iel,3)
                    enddo
                enddo
            enddo
        endif
        ! endif nwftype

    endif
    ! endif iflag
       
#endif
!  endif vectorization
return  
end

subroutine orbitals_quad_no_qmckl_periodic(nxquad,xquad,rvec_en,r_en,orbn,dorbn,da_orbn,iwforb)

    use basis_fns_mod, only: basis_fns
    use coefs,   only: nbasis
    use multiple_geo, only: iwf
    use numbas2, only: ibas0,ibas1
    use m_force_analytic, only: iforce_analy
    use optwf_control, only: ioptorb
    use optwf_control, only: method
    use orbval,  only: nadorb
    use phifun,  only: dphin,n0_ibasis,n0_ic,n0_nbasis,phin
    use precision_kinds, only: dp
    use qua,     only: nquad
    use slater,  only: coef,norb
    use sr_mod,  only: i_sr_rescale
    use system,  only: ncent,ncent_tot,nelec
    use vmc_mod, only: norb_tot, nwftypeorb
    use contrl_file, only: ounit
    use periodic, only : n_images, ell

    implicit none

    integer :: ic, ider, iq
    integer :: iorb, k, m, m0, nxquad, iwforb
    integer :: nadorb_sav, i_image

    real(dp), dimension(3,*) :: xquad
    real(dp), dimension(nquad*nelec*2, ncent_tot) :: r_en
    real(dp), dimension(3,nquad*nelec*2, ncent_tot) :: rvec_en
    real(dp), dimension(nquad*nelec*2, ncent_tot) :: r_en_image
    real(dp), dimension(3,nquad*nelec*2, ncent_tot) :: rvec_en_image
    real(dp), dimension(norb_tot, *) :: orbn
    real(dp), dimension(norb_tot, nquad*nelec*2, 3) :: dorbn
    real(dp), dimension(norb,3,nxquad,ncent_tot) :: da_orbn
    real(dp), dimension(3) :: dtmp
    real(dp) :: ddtmp

    nadorb_sav=nadorb

    if(ioptorb.eq.0.or.(method(1:3).ne.'lin'.and.i_sr_rescale.eq.0)) nadorb=0

    ! get basis functions for electron iel
    ider=0
    if(iforce_analy.gt.0) then
         ider=1
        !
        do iorb=1,norb
            do ic=1,ncent
                do k=1,3
                    da_orbn(iorb,k,iq,ic)=0.d0
                enddo
            enddo
            do k=1,3
                dorbn(iorb,iq,k)=0.d0
            enddo
        enddo
    end if
    
    orbn(1:(norb+nadorb),1:nxquad)=0.d0

    do i_image=1,n_images
        rvec_en_image(1, :, :) = rvec_en(1, :, :) - ell(1, i_image)
        rvec_en_image(2, :, :) = rvec_en(2, :, :) - ell(2, i_image)
        rvec_en_image(3, :, :) = rvec_en(3, :, :) - ell(3, i_image)
        r_en_image(:,:) = sqrt(rvec_en_image(1,:,:)**2 + rvec_en_image(2,:,:)**2 + rvec_en_image(3,:,:)**2)

        if(nwftypeorb.gt.1) iwf=1
        call basis_fns(1,nxquad,nquad*nelec*2,rvec_en_image,r_en_image,ider)
        if(nwftypeorb.gt.1) iwf=iwforb

        do iq=1,nxquad

    ! Vectorization dependent code selection
    #ifdef VECTORIZATION
    ! The following loop changed for better vectorization AVX512/AVX2
            do iorb=1,norb+nadorb
                do m=1,nbasis
                    orbn(iorb,iq)=orbn(iorb,iq)+coef(m,iorb,iwf)*phin(m,iq)
                enddo
            enddo
    #else
            do iorb=1,norb+nadorb
                do m0=1,n0_nbasis(iq)
                    m=n0_ibasis(m0,iq)
                    orbn(iorb,iq)=orbn(iorb,iq)+coef(m,iorb,iwf)*phin(m,iq)
                enddo
            enddo
    #endif

            if(iforce_analy.gt.0) then
                do iorb=1,norb
    #ifdef VECTORIZATION
                    do ic=1,ncent
                        do k=1,3
                            do m=ibas0(ic),ibas1(ic)
                                da_orbn(iorb,k,iq,ic)=da_orbn(iorb,k,iq,ic)-coef(m,iorb,iwf)*dphin(m,iq,k)
                            enddo
                        enddo
                    enddo
    #else
                    do m0=1,n0_nbasis(iq)
                        m=n0_ibasis(m0,iq)
                        ic=n0_ic(m0,iq)
                        do k=1,3
                            da_orbn(iorb,k,iq,ic)=da_orbn(iorb,k,iq,ic)-coef(m,iorb,iwf)*dphin(m,iq,k)
                        enddo
                    enddo
    #endif
                    do ic=1,ncent
                        do k=1,3
                            dorbn(iorb,iq,k)=dorbn(iorb,iq,k)-da_orbn(iorb,k,iq,ic)
                        enddo
                    enddo
                enddo
            endif
            ! endiff iforce
        enddo
        ! enddo nxquad

        nadorb = nadorb_sav
    enddo
return
end

end module