module orbitals_no_qmckl_mod
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

subroutine orbitals_no_qmckl(x,rvec_en,r_en)

    use basis_fns_mod, only: basis_fns
    use coefs, only: nbasis
    use multiple_geo, only: iwf
    use m_force_analytic, only: iforce_analy
    use orbval, only: ddorb, dorb, nadorb, orb
    use phifun, only: phin, dphin, d2phin, n0_ibasis, n0_nbasis
    use precision_kinds, only: dp
    use slater, only: norb, coef
    use system, only: ncent_tot, nelec
    use vmc_mod, only: nwftypeorb
        
    implicit none
        
    integer :: i, ider, iorb, k, m
    integer :: m0, j

    real(dp), dimension(3,*) :: x
    real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
    real(dp), dimension(nelec,ncent_tot) :: r_en
!     real(dp), dimension(nelec,nbasis) :: bhin
!     real(dp), dimension(3*nelec,nbasis) :: dbhin
!     real(dp), dimension(nelec,nbasis) :: d2bhin
    real(dp), dimension(:), allocatable :: auxorb !(norb+nadorb)
    real(dp), dimension(:, :), allocatable :: auxdorb !(norb+nadorb)
    real(dp), dimension(:), allocatable :: auxddorb !(norb+nadorb)
    if (.not. allocated(auxorb)) allocate (auxorb(norb+nadorb))
    if (.not. allocated(auxdorb)) allocate (auxdorb(norb+nadorb,3))
    if (.not. allocated(auxddorb)) allocate (auxddorb(norb+nadorb))

    ! get basis functions for all electrons
    ider=2
    if(iforce_analy.eq.1) ider=3

    call basis_fns(1,nelec,nelec,rvec_en,r_en,ider)

    !     in alternativa al loop 26
    !     do jbasis=1,nbasis
    !     i=0
    !     do ielec=1,nelec
    !     bhin(ielec,jbasis)=phin(jbasis,ielec)
    !     do l=1,3
    !     i=i+1
    !     dbhin(i,jbasis)=dphin(jbasis,ielec,l)
    !     enddo
    !     d2bhin(ielec,jbasis)=d2phin(jbasis,ielec)
    !     enddo
    !     enddo
    !     call dgemm('n','n',  nelec,norb,nbasis,1.d0,bhin,   nelec,  coef(1,1,iwf),nbasis,0.d0,orb,   nelec)
    !     call dgemm('n','n',3*nelec,norb,nbasis,1.d0,dbhin,3*nelec,  coef(1,1,iwf),nbasis,0.d0,dorb,3*nelec)
    !     call dgemm('n','n',  nelec,norb,nbasis,1.d0,d2bhin, nelec,  coef(1,1,iwf),nbasis,0.d0,ddorb, nelec)
    
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
    !  nwftype endif
    
#endif
! vectorization endif
return
end


subroutine orbitalse_no_qmckl(iel,x,rvec_en,r_en,iflag)

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

    implicit none

    integer :: iel, ider, iflag, iorb, m
    integer :: m0, k, j

    real(dp), dimension(3,*) :: x
    real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
    real(dp), dimension(nelec,ncent_tot) :: r_en


    ider=1
    if(iflag.gt.0) ider=2

    call basis_fns(iel,iel,nelec,rvec_en,r_en,ider)

! Vectorization dependent code. useful for AVX512 and AVX2
#ifdef VECTORIZATION
       
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
       !endif nwftype

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


    endif
    ! endif nwtype
       
#else
! Keep the localization for the non-vectorized code
       
    if(iflag.gt.0) then
        if(nwftypeorb.gt.1) then
            do k=1,nwftypeorb
                do iorb=1,norb
                    orbn(iorb,k)=0.d0
                    dorbn(iorb,1,k)=0.d0
                    dorbn(iorb,2,k)=0.d0
                    dorbn(iorb,3,k)=0.d0
                    ddorbn(iorb,k)=0.d0
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
        else
            do iorb=1,norb
                orbn(iorb,1)=0.d0
                dorbn(iorb,1,1)=0.d0
                dorbn(iorb,2,1)=0.d0
                dorbn(iorb,3,1)=0.d0
                ddorbn(iorb,1)=0.d0
                do m0=1,n0_nbasis(iel)
                    m=n0_ibasis(m0,iel)
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
                    do m0=1,n0_nbasis(iel)
                        m=n0_ibasis(m0,iel)
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
                do m0=1,n0_nbasis(iel)
                    m=n0_ibasis(m0,iel)
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
       
#endif
!  endif vectorization
return
end

subroutine orbitals_quad_no_qmckl(nxquad,xquad,rvec_en,r_en,orbn,dorbn,da_orbn,iwforb)

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

    implicit none

    integer :: ic, ider, iq
    integer :: iorb, k, m, m0, nxquad, iwforb
    integer :: nadorb_sav

    real(dp), dimension(3,*) :: xquad
    real(dp), dimension(nquad*nelec*2, ncent_tot) :: r_en
    real(dp), dimension(3,nquad*nelec*2, ncent_tot) :: rvec_en
    real(dp), dimension(norb_tot, *) :: orbn
    real(dp), dimension(norb_tot, nquad*nelec*2, 3) :: dorbn
    real(dp), dimension(norb,3,nxquad,ncent_tot) :: da_orbn
    real(dp), dimension(3) :: dtmp
    real(dp) :: ddtmp

    nadorb_sav=nadorb

    if(ioptorb.eq.0.or.(method(1:3).ne.'lin'.and.i_sr_rescale.eq.0)) nadorb=0

    ! get basis functions for electron iel
    ider=0
    if(iforce_analy.gt.0) ider=1

    if(nwftypeorb.gt.1) iwf=1
    call basis_fns(1,nxquad,nquad*nelec*2,rvec_en,r_en,ider)
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
                        da_orbn(iorb,k,iq,ic)=0.d0
                    enddo
                enddo
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
                do k=1,3
                    dorbn(iorb,iq,k)=0.d0
                enddo
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

return
end

end module
