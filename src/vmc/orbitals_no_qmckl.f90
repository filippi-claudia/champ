module orbitals_no_qmckl_mod
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
      integer :: kcoef, m0

      real(dp), dimension(3,*) :: x
      real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
      real(dp), dimension(nelec,ncent_tot) :: r_en
      real(dp), dimension(norb+nadorb) :: auxorb
      real(dp), dimension(norb+nadorb,3) :: auxdorb
      real(dp), dimension(norb+nadorb) :: auxddorb

      ! get basis functions for all electrons
      ider=2
      if(iforce_analy.eq.1) ider=3

      call basis_fns(1,nelec,nelec,rvec_en,r_en,ider)

      ! Vectorization dependent code selection
#ifdef VECTORIZATION
      ! Following loop changed for better vectorization AVX512/AVX2
      do k=1,nwftypeorb
        kcoef=k
        if(nwftypeorb.eq.1) kcoef=iwf
        do i=1,nelec
          auxorb=0.d0
          auxdorb=0.d0
          auxddorb=0.d0
          do iorb=1,norb+nadorb
            do m=1,nbasis
              auxorb(iorb)=auxorb(iorb)+coef(m,iorb,kcoef)*phin(m,i)
              auxdorb(iorb,1)=auxdorb(iorb,1)+coef(m,iorb,kcoef)*dphin(m,i,1)
              auxdorb(iorb,2)=auxdorb(iorb,2)+coef(m,iorb,kcoef)*dphin(m,i,2)
              auxdorb(iorb,3)=auxdorb(iorb,3)+coef(m,iorb,kcoef)*dphin(m,i,3)
              auxddorb(iorb)=auxddorb(iorb)+coef(m,iorb,kcoef)*d2phin(m,i)
            enddo
          enddo
          orb(i,1:(norb+nadorb),k)=auxorb(1:(norb+nadorb))
          dorb(1:(norb+nadorb),i,1:3,k)=auxdorb(1:(norb+nadorb),1:3)
          ddorb(1:(norb+nadorb),i,k)=auxddorb(1:(norb+nadorb))
        enddo
      enddo

#else
      ! Use localized basis lists when vectorization is not enabled.

      do k=1,nwftypeorb
        kcoef=k
        if(nwftypeorb.eq.1) kcoef=iwf
        do i=1,nelec
          auxorb=0.d0
          auxdorb=0.d0
          auxddorb=0.d0
          do iorb=1,norb+nadorb
            do m0=1,n0_nbasis(i)
              m=n0_ibasis(m0,i)
              auxorb(iorb)=auxorb(iorb)+coef(m,iorb,kcoef)*phin(m,i)
              auxdorb(iorb,1)=auxdorb(iorb,1)+coef(m,iorb,kcoef)*dphin(m,i,1)
              auxdorb(iorb,2)=auxdorb(iorb,2)+coef(m,iorb,kcoef)*dphin(m,i,2)
              auxdorb(iorb,3)=auxdorb(iorb,3)+coef(m,iorb,kcoef)*dphin(m,i,3)
              auxddorb(iorb)=auxddorb(iorb)+coef(m,iorb,kcoef)*d2phin(m,i)
            enddo
          enddo
          orb(i,1:(norb+nadorb),k)=auxorb(1:(norb+nadorb))
          dorb(1:(norb+nadorb),i,1:3,k)=auxdorb(1:(norb+nadorb),1:3)
          ddorb(1:(norb+nadorb),i,k)=auxddorb(1:(norb+nadorb))
        enddo
      enddo

#endif
      ! endif vectorization
      end subroutine orbitals_no_qmckl

      subroutine orbitalse_no_qmckl(iel,x,rvec_en,r_en,iflag)

      use basis_fns_mod, only: basis_fns
      use coefs, only: nbasis
      use multiple_geo, only: iwf
      use multislatern, only: ddorbn, dorbn, orbn
      use phifun, only: d2phin, dphin, n0_ibasis, n0_nbasis, phin
      use precision_kinds, only: dp
      use slater, only: norb, coef
      use system, only: ncent_tot, nelec
      use vmc_mod, only: nwftypeorb

      implicit none

      integer :: iel, ider, iflag, iorb, m
      integer :: k, kcoef, m0

      real(dp), dimension(3,*) :: x
      real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
      real(dp), dimension(nelec,ncent_tot) :: r_en

      ider=1
      if(iflag.gt.0) ider=2

      call basis_fns(iel,iel,nelec,rvec_en,r_en,ider)

      ! Vectorization dependent code, useful for AVX512 and AVX2.
#ifdef VECTORIZATION

      do k=1,nwftypeorb
        kcoef=k
        if(nwftypeorb.eq.1) kcoef=iwf
        orbn(1:norb,k)=0.d0
        dorbn(1:norb,1,k)=0.d0
        dorbn(1:norb,2,k)=0.d0
        dorbn(1:norb,3,k)=0.d0
        do iorb=1,norb
          do m=1,nbasis
            orbn(iorb,k)=orbn(iorb,k)+coef(m,iorb,kcoef)*phin(m,iel)
            dorbn(iorb,1,k)=dorbn(iorb,1,k)+coef(m,iorb,kcoef)*dphin(m,iel,1)
            dorbn(iorb,2,k)=dorbn(iorb,2,k)+coef(m,iorb,kcoef)*dphin(m,iel,2)
            dorbn(iorb,3,k)=dorbn(iorb,3,k)+coef(m,iorb,kcoef)*dphin(m,iel,3)
          enddo
        enddo
        if(iflag.gt.0) then
          ddorbn(1:norb,k)=0.d0
          do iorb=1,norb
            do m=1,nbasis
              ddorbn(iorb,k)=ddorbn(iorb,k)+coef(m,iorb,kcoef)*d2phin(m,iel)
            enddo
          enddo
        endif
      enddo

#else
      ! Use localized basis lists when vectorization is not enabled.

      do k=1,nwftypeorb
        kcoef=k
        if(nwftypeorb.eq.1) kcoef=iwf
        orbn(1:norb,k)=0.d0
        dorbn(1:norb,1,k)=0.d0
        dorbn(1:norb,2,k)=0.d0
        dorbn(1:norb,3,k)=0.d0
        do iorb=1,norb
          do m0=1,n0_nbasis(iel)
            m=n0_ibasis(m0,iel)
            orbn(iorb,k)=orbn(iorb,k)+coef(m,iorb,kcoef)*phin(m,iel)
            dorbn(iorb,1,k)=dorbn(iorb,1,k)+coef(m,iorb,kcoef)*dphin(m,iel,1)
            dorbn(iorb,2,k)=dorbn(iorb,2,k)+coef(m,iorb,kcoef)*dphin(m,iel,2)
            dorbn(iorb,3,k)=dorbn(iorb,3,k)+coef(m,iorb,kcoef)*dphin(m,iel,3)
          enddo
        enddo
        if(iflag.gt.0) then
          ddorbn(1:norb,k)=0.d0
          do iorb=1,norb
            do m0=1,n0_nbasis(iel)
              m=n0_ibasis(m0,iel)
              ddorbn(iorb,k)=ddorbn(iorb,k)+coef(m,iorb,kcoef)*d2phin(m,iel)
            enddo
          enddo
        endif
      enddo

#endif
      ! endif vectorization
      end subroutine orbitalse_no_qmckl

      subroutine orbitals_quad_no_qmckl(nxquad,xquad,rvec_en,r_en,orbn,dorbn,da_orbn,iwforb)

      use basis_fns_mod, only: basis_fns
      use coefs,   only: nbasis
      use multiple_geo, only: iwf
      use numbas2, only: ibas0,ibas1
      use m_force_analytic, only: iforce_analy
      use optwf_control, only: ioptorb, method
      use orbval,  only: nadorb
      use phifun,  only: dphin,n0_ibasis,n0_ic,n0_nbasis,phin
      use precision_kinds, only: dp
      use qua,     only: nquad
      use slater,  only: coef,norb
      use sr_mod,  only: i_sr_rescale
      use system,  only: ncent,ncent_tot,nelec
      use vmc_mod, only: norb_tot, nwftypeorb

      implicit none

      integer :: ic, ider, iq
      integer :: iorb, k, m, m0, nxquad, iwforb
      integer :: norb_eval

      real(dp), dimension(3,*) :: xquad
      real(dp), dimension(nquad*nelec*2, ncent_tot) :: r_en
      real(dp), dimension(3,nquad*nelec*2, ncent_tot) :: rvec_en
      real(dp), dimension(norb_tot, *) :: orbn
      real(dp), dimension(norb_tot, nquad*nelec*2, 3) :: dorbn
      real(dp), dimension(norb,3,nxquad,ncent_tot) :: da_orbn

      norb_eval=norb+nadorb
      if(ioptorb.eq.0.or.(method(1:3).ne.'lin'.and.i_sr_rescale.eq.0)) norb_eval=norb

      ! get basis functions for quadrature points
      ider=0
      if(iforce_analy.gt.0) ider=1

      if(nwftypeorb.gt.1) iwf=1
      call basis_fns(1,nxquad,nquad*nelec*2,rvec_en,r_en,ider)
      if(nwftypeorb.gt.1) iwf=iwforb

      do iq=1,nxquad

      ! Vectorization dependent code selection
#ifdef VECTORIZATION
      ! The following loop changed for better vectorization AVX512/AVX2
        orbn(1:norb_eval,iq)=0.d0
        do iorb=1,norb_eval
          do m=1,nbasis
            orbn(iorb,iq)=orbn(iorb,iq)+coef(m,iorb,iwf)*phin(m,iq)
          enddo
        enddo
#else
        orbn(1:norb_eval,iq)=0.d0
        do iorb=1,norb_eval
          do m0=1,n0_nbasis(iq)
            m=n0_ibasis(m0,iq)
            orbn(iorb,iq)=orbn(iorb,iq)+coef(m,iorb,iwf)*phin(m,iq)
          enddo
        enddo
#endif

        if(iforce_analy.gt.0) then
          do iorb=1,norb
            do ic=1,ncent
              da_orbn(iorb,1:3,iq,ic)=0.d0
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
            dorbn(iorb,iq,1:3)=0.d0
            do ic=1,ncent
              dorbn(iorb,iq,1:3)=dorbn(iorb,iq,1:3)-da_orbn(iorb,1:3,iq,ic)
            enddo
          enddo
        endif
        ! endif iforce
      enddo
      ! enddo nxquad

      end subroutine orbitals_quad_no_qmckl

end module orbitals_no_qmckl_mod
