module orbitals_mod
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
      subroutine orbitals(x,rvec_en,r_en)
! Written by Cyrus Umrigar starting from Kevin Schmidt's routine
! Modified by A. Scemama

      use control, only: ipr
      use contrl_per, only: iperiodic
      use contrl_file, only: ounit
      use m_force_analytic, only: iforce_analy
      use orbitals_no_qmckl_mod, only: orbitals_no_qmckl
      use orbval, only: ddorb, dorb, nadorb, orb
      use precision_kinds, only: dp
      use slater, only: norb, coef
      use system, only: ncent_tot, nelec
      use vmc_mod, only: nwftypeorb

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
      use orbitals_qmckl_periodic_mod, only: orbitals_qmckl_periodic
      use orbitals_qmckl_mod, only: orbitals_qmckl
#endif

      implicit none

      integer :: i, j, iorb, k

      real(dp), dimension(3,*) :: x
      real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
      real(dp), dimension(nelec,ncent_tot) :: r_en

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
      if (iperiodic.eq.0) then
         call orbitals_qmckl(x,rvec_en,r_en)
      else
         call orbitals_qmckl_periodic(x,rvec_en,r_en)
      endif
#else
      call orbitals_no_qmckl(x,rvec_en,r_en)
#endif

      if(iforce_analy.eq.1) call da_orbitals

      if(ipr.ge.0) then
         do j=1,nwftypeorb
            do iorb=1,norb+nadorb
               write(ounit,'(''orb set,iorb,orb='',2i4,1000f15.11)') j,iorb,(orb(i,iorb,j),i=1,nelec)
            enddo
            do iorb=1,norb+nadorb
               write(ounit,'(''orb set,iorb,d2orb='',2i4,1000f15.11)') j,iorb,(ddorb(iorb,i,j),i=1,nelec)
            enddo
            do k=1,3
               do iorb=1,norb+nadorb
                  write(ounit,'(''orb set,dir,iorb,dorb='',3i4,1000f12.8)') j,k,iorb,(dorb(iorb,i,k,j),i=1,nelec)
               enddo
            enddo
         enddo
      endif

      return
      end
!------------------------------------------------------------------------------------

      subroutine da_orbitals

      use system, only: ncent, nelec
      use da_orbval, only: da_d2orb, da_dorb, da_orb
      use numbas2, only: ibas0, ibas1
      use phifun, only: d2phin_all, d3phin, dphin
      use multiple_geo, only: iwf
      use coefs, only: nbasis
      use slater, only: norb, coef
      use precision_kinds, only: dp
      use slater,  only: coef,norb
      use system,  only: ncent,nelec
      use contrl_file, only: ounit
      use error,   only: fatal_error


#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND)
      use qmckl_data
#endif

      implicit none

      integer :: ibasis, i, ic, ielec, j, k
      integer :: l, m, n

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND)

      integer(qmckl_exit_code) :: rc

      real(dp), dimension(:,:,:,:),allocatable :: da_orb_two
      real(dp), dimension(:,:,:,:,:),allocatable :: da_dorb_two
      real(dp), dimension(:,:,:,:),allocatable :: da_d2orb_two

      double precision, parameter :: alpha = 1.0d0, beta = 0.0d0
      double precision :: identity(3,3)

      allocate (da_dorb_two(norb,3,  nelec, 3, ncent))
      allocate (da_orb_two(norb,3,nelec,ncent))
      allocate (da_d2orb_two(norb,nelec,3,ncent))

      rc = qmckl_get_forces_mo_value_inplace(qmckl_ctx(1), da_orb_two, nelec*norb*3_8*ncent)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCKL forces of MO values.')

      rc = qmckl_get_forces_mo_g_inplace(qmckl_ctx(1), da_dorb_two, 3*nelec*norb*3_8*ncent)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCKL forces of MO gradients/')

      rc = qmckl_get_forces_mo_l(qmckl_ctx(1), da_d2orb_two, nelec*norb*3_8*ncent)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCKL forces of MO laplacian.')

      ! do j = 1, norb
      !       do i = 1, nelec
      !             do k = 1, 3
      !                   do ic = 1, ncent
      !                         da_orb(k, i, j, ic) = da_orb_two(j, i, k, ic)
      !                         da_d2orb(k, i, j, ic) = da_d2orb_two(j, i, k, ic)
      !                         do l = 1, 3
      !                               da_dorb(k, l, i, j, ic) = da_dorb_two(j, l, i, k, ic)
      !                         enddo
      !                   enddo
      !             enddo
      !       enddo
      ! enddo

      do ic = 1, ncent
            do k = 1,3
                  do i = 1, nelec
                        do j = 1, norb
                              da_orb(k, i, j, ic) = da_orb_two(j, k, i, ic)
                              da_d2orb(k, i, j, ic) = da_d2orb_two(j, i, k, ic)
                        enddo
                        do l = 1, 3
                              do j = 1, norb
                              da_dorb(k, l, i, j, ic) = da_dorb_two(j, l, i, k, ic)
                              enddo
                       enddo
                  enddo
            enddo
      enddo

  
      ! Create a 3x3 identity matrix
      ! identity = reshape([1.0d0, 0.0d0, 0.0d0, &
      !                      0.0d0, 1.0d0, 0.0d0, &
      !                      0.0d0, 0.0d0, 1.0d0], [3,3])
  
      ! do ic = 1, ncent
      !     ! Transpose da_orb using dgemm
      !     call dgemm('T', 'N', 3, nelec * norb, 3, alpha, da_orb_two(1, 1, 1, ic), nelec, &
      !                identity, 3, beta, da_orb(1, 1, 1, ic), 3)
  
      !     ! Transpose da_d2orb using dgemm
      !     call dgemm('T', 'N', 3, nelec * norb, 3, alpha, da_d2orb_two(1, 1, 1, ic), 3, &
      !                identity, 3, beta, da_d2orb(1, 1, 1, ic), 3)
  
      !     ! Transpose da_dorb using dgemm
      !     call dgemm('T', 'N', 3 * 3, nelec * norb, 3, alpha, da_dorb_two(1, 1, 1, 1, ic), 3, &
      !                identity, 3, beta, da_dorb(1, 1, 1, 1, ic), 3 * 3)
      ! end do

      deallocate(da_orb_two,da_dorb_two,da_d2orb_two)
#else
      
      real(dp), dimension(3*nelec,nbasis) :: tphin
      real(dp), dimension(3*3*nelec,nbasis) :: t2phin_all
      real(dp), dimension(3*nelec,nbasis) :: t3phin

      do ibasis=1,nbasis
       i=0
       j=0
       do ielec=1,nelec
        do l=1,3
         i=i+1
         tphin(i,ibasis)=dphin(ibasis,ielec,l)
         t3phin(i,ibasis)=d3phin(l,ibasis,ielec)
         do k=1,3
          j=j+1
          t2phin_all(j,ibasis)=d2phin_all(k,l,ibasis,ielec)
         enddo
        enddo
       enddo
      enddo
      n=3*nelec
      m=3*nelec
      do ic=1,ncent
        k=ibas1(ic)-ibas0(ic)+1
        j=ibas0(ic)
        call dgemm('n','n',  n,norb,k,-1.d0,tphin(1,j)     ,  m,coef(j,1,iwf),nbasis,0.d0,da_orb(1,1,1,ic)   ,  m)
        call dgemm('n','n',  n,norb,k,-1.d0,t3phin(1,j)    ,  m,coef(j,1,iwf),nbasis,0.d0,da_d2orb(1,1,1,ic) ,  m)
        call dgemm('n','n',3*n,norb,k,-1.d0,t2phin_all(1,j),3*m,coef(j,1,iwf),nbasis,0.d0,da_dorb(1,1,1,1,ic),3*m)
      enddo

#endif

      ! do ic = 1, ncent
      !       do i = 1, nelec
      !             write(ounit, *), 'da_orb', (da_orb(k,i,1,ic), k=1, 3)
      !             write(ounit, *), 'da_d2orb', (da_d2orb(k,i,1,ic), k=1, 3)
      !             write(ounit, *), 'da_dorb', (da_dorb(k,1,i,1,ic), k=1, 3)

      !       enddo
      ! enddo
      return
      end
!------------------------------------------------------------------------------------
      subroutine orbitalse(iel,x,rvec_en,r_en,iflag)

      use contrl_per, only: iperiodic
      use orbitals_no_qmckl_mod, only: orbitalse_no_qmckl
      use precision_kinds, only: dp
      use system, only: ncent_tot, nelec

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
      use orbitals_qmckl_periodic_mod, only: orbitalse_qmckl_periodic
      use orbitals_qmckl_mod, only: orbitalse_qmckl
#endif

      implicit none

      integer :: iel, iflag

      real(dp), dimension(3,*) :: x
      real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
      real(dp), dimension(nelec,ncent_tot) :: r_en

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
      if (iperiodic.eq.0) then
         call orbitalse_qmckl(iel,x,rvec_en,r_en,iflag)
      else
         call orbitalse_qmckl_periodic(iel,x,rvec_en,r_en,iflag)
      endif
#else
      call orbitalse_no_qmckl(iel,x,rvec_en,r_en,iflag)
#endif

      return
      end
!------------------------------------------------------------------------------------
end module
