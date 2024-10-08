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
      use multiple_geo, only: iwf
      use phifun, only: phin, dphin, d2phin, n0_ibasis, n0_nbasis
      use coefs, only: nbasis
      use slater, only: norb, coef
      use contrl_per, only: iperiodic
      use m_force_analytic, only: iforce_analy
      use grid3dflag, only: i3dlagorb, i3dsplorb
      use system, only: ncent_tot, nelec
      use orbval, only: ddorb, dorb, nadorb, orb
      use precision_kinds, only: dp
      use contrl_file,    only: ounit
      use grid3d_orbitals, only: spline_mo
      use grid3d_orbitals, only: lagrange_mos, lagrange_mos_grad, lagrange_mos_2
      use basis_fns_mod, only: basis_fns
      use coefs,   only: nbasis
      use contrl_file, only: ounit
      use contrl_per, only: iperiodic
      use control, only: ipr
      use grid3d_orbitals, only: lagrange_mos,lagrange_mos_2
      use grid3d_orbitals, only: lagrange_mos_grad,spline_mo
      use grid3dflag, only: i3dlagorb,i3dsplorb
      use m_force_analytic, only: iforce_analy
      use multiple_geo, only: iwf
      use orbval,  only: ddorb,dorb,nadorb,orb
      use phifun,  only: d2phin,dphin,n0_ibasis,n0_nbasis,phin
      use vmc_mod, only: nwftypeorb
      use periodic, only : n_images, ell
      use find_pimage, only: find_image_pbc
      use precision_kinds, only: dp

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
      use const
      use qmckl_data
#endif

      implicit none

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
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
#endif


      integer :: i, ier, ider, iorb, k, m
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

      ier=1

!     spline interpolation
      if(i3dsplorb.eq.2) then
         do k=1,nwftypeorb
            do i=1,nelec
               ier = 0.d0
               do iorb=1,norb+nadorb
                  ddorb(iorb,i,k)=1.d0 ! compute the laplacian
                  dorb(iorb,i,1,k)=1.d0 ! compute the gradients
                  dorb(iorb,i,2,k)=1.d0 ! compute the gradients
                  dorb(iorb,i,3,k)=1.d0 ! compute the gradients
                  call spline_mo(x(1,i),iorb,orb(i,iorb,k),dorb(iorb,i,:,k),ddorb(iorb,i,k),ier)
               enddo
               if(ier.eq.1) then

                  call basis_fns(i,i,nelec,rvec_en,r_en,2)

                  do iorb=1,norb+nadorb
                     orb(i,iorb,k)=0.d0
                     dorb(iorb,i,1,k)=0.d0
                     dorb(iorb,i,2,k)=0.d0
                     dorb(iorb,i,3,k)=0.d0
                     ddorb(iorb,i,k)=0.d0
                     do m0=1,n0_nbasis(i)
                        m=n0_ibasis(m0,i)
                        orb(i,iorb,k)=orb(i,iorb,k)+coef(m,iorb,iwf)*phin(m,i)
                        dorb(iorb,i,1,k)=dorb(iorb,i,1,k)+coef(m,iorb,iwf)*dphin(m,i,1)
                        dorb(iorb,i,2,k)=dorb(iorb,i,2,k)+coef(m,iorb,iwf)*dphin(m,i,2)
                        dorb(iorb,i,3,k)=dorb(iorb,i,3,k)+coef(m,iorb,iwf)*dphin(m,i,3)
                        ddorb(iorb,i,k)=ddorb(iorb,i,k)+coef(m,iorb,iwf)*d2phin(m,i)
                     enddo
                  enddo
               endif
            enddo
         enddo
!
!     Lagrange interpolation, did not inclue multiorb here yet
      elseif(i3dlagorb.eq.2) then
         do k=1,nwftypeorb
            do i=1,nelec
               ier=0
               call lagrange_mos(1,x(1,i),orb(1,1,k),i,ier)
               call lagrange_mos_grad(2,x(1,i),dorb(1,1,1,k),i,ier)
               call lagrange_mos_grad(3,x(1,i),dorb(1,1,1,k),i,ier)
               call lagrange_mos_grad(4,x(1,i),dorb(1,1,1,k),i,ier)
               call lagrange_mos_2(5,x(1,i),ddorb(1,1,k),i,ier)

               if(ier.eq.1) then

                  if(nwftypeorb.gt.1) iwf=1
                  call basis_fns(i,i,nelec,rvec_en,r_en,2)
                  if(nwftypeorb.gt.1) iwf=k

                  do iorb=1,norb+nadorb
                     orb(i,iorb,k)=0.d0
                     dorb(iorb,i,1,k)=0.d0
                     dorb(iorb,i,2,k)=0.d0
                     dorb(iorb,i,3,k)=0.d0
                     ddorb(iorb,i,k)=0.d0
                     do m0=1,n0_nbasis(i)
                        m=n0_ibasis(m0,i)
                        orb(i,iorb,k)=orb(i,iorb,k)+coef(m,iorb,iwf)*phin(m,i)
                        dorb(iorb,i,1,k)=dorb(iorb,i,1,k)+coef(m,iorb,iwf)*dphin(m,i,1)
                        dorb(iorb,i,2,k)=dorb(iorb,i,2,k)+coef(m,iorb,iwf)*dphin(m,i,2)
                        dorb(iorb,i,3,k)=dorb(iorb,i,3,k)+coef(m,iorb,iwf)*dphin(m,i,3)
                        ddorb(iorb,i,k)=ddorb(iorb,i,k)+coef(m,iorb,iwf)*d2phin(m,i)
                     enddo
                  enddo
               endif
            enddo
         enddo

!     no 3d interpolation
      else

!     get basis functions for all electrons
         ider=2
         if(iforce_analy.eq.1) ider=3

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 

         if(iperiodic.eq.0) then

!     get number MO's
         rc = qmckl_get_mo_basis_mo_num(qmckl_ctx, n8)
         if (rc /= QMCKL_SUCCESS) then
            print *, 'Error getting mo_num from QMCkl'
            stop
         end if


         allocate(mo_vgl_qmckl(n8, 5, nelec))


!     Send electron coordinates to QMCkl to compute the MOs at these positions
         rc = qmckl_set_point(qmckl_ctx, 'N', nelec*1_8, x, nelec*3_8)

         if (rc /= QMCKL_SUCCESS) then
            print *, 'Error setting electron coordinates in QMCkl'
            stop
         end if


!     Compute the MOs
         rc = qmckl_get_mo_basis_mo_vgl_inplace( &
              qmckl_ctx, &
              mo_vgl_qmckl, &
              n8*nelec*5_8)

         if (rc /= QMCKL_SUCCESS) then
            print *, 'Error getting MOs from QMCkl'
            stop
         end if


!     print*, "inside qmckl"

!     pass computed qmckl orbitals back to champ

         k=1                    ! until state specific orbitals can be used
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


      else
! else iperiodic

! call to be replaced with QMCkl
!     call basis_fns(1,nelec,nelec,rvec_en,r_en,ider)

!     ! get number of atomic orbitals
         rc = qmckl_get_ao_basis_ao_num(qmckl_ctx, na8)
         if (rc /= QMCKL_SUCCESS) then
            print *, 'Error getting mo_num from QMCkl'
            stop
         end if

!     print*,"na8",na8
!     print*,"nbasis",nbasis


         if (nbasis.ne.na8) then
            print *, 'Error getting ao_num from QMCkl'
            stop
         end if


!     first image zero

!     allocate ao_vlg array
         allocate(ao_qmckl(nbasis, 5, nelec))
         ao_qmckl=0.d0

!     ! set initial coordinates for electrons zero image
         xelec=x(1:3,1:nelec)

!     !setting pbc conditions
         do i_elec=1, nelec
            call find_image_pbc(xelec(1:3,i_elec),rnorm)
         enddo

!     Send electron coordinates to QMCkl to compute the MOs at these positions
         rc = qmckl_set_point(qmckl_ctx, 'N', nelec*1_8, xelec, nelec*3_8)
         if (rc /= QMCKL_SUCCESS) then
            print *, 'Error setting electron coordinates in QMCkl'
            stop
         end if


!     computing aos zero image
         rc = qmckl_get_ao_basis_ao_vgl_inplace(qmckl_ctx, &
              ao_qmckl, nbasis*5_8*nelec)
         if (rc /= QMCKL_SUCCESS) then
            print *, 'Error getting AOs from QMCkl zero image'
         endif

!     computing images distance for iel
         if(n_images.gt.0) then
!     allocate ao_vgl array for all images
            allocate(ao_vgl_qmckl(nbasis, 5, nelec))

            do i_image=1, n_images

! set electron images distances
               xqmckl=0.d0
               r_image=ell(1:3,i_image)
               do i_elec=1, nelec
                  xqmckl(1:3,i_elec)=xelec(1:3,i_elec)-r_image(1:3)
               enddo

!     send electron images coordinates
               rc = qmckl_set_point(qmckl_ctx, 'N', 1_8*nelec, xqmckl, 3_8*nelec)
               if (rc /= QMCKL_SUCCESS) then
                  print *, 'Error setting electron coords orbitalse'
                  call qmckl_last_error(qmckl_ctx,err_message)
                  print *, trim(err_message)
                  call abort()
               end if


               ao_vgl_qmckl=0.d0
!     computing aos for the given image
               rc = qmckl_get_ao_basis_ao_vgl_inplace(qmckl_ctx, &
                    ao_vgl_qmckl, nbasis*5_8*nelec)
               if (rc /= QMCKL_SUCCESS) then
                  print *, 'Error getting AOs from QMCkl zero image'
               endif


!     adding contribution of the given image
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
         !endif images

!! pass QMCkl ao's to champ
         do i_elec=1, nelec
            do i_basis=1, nbasis
               phin(i_basis,i_elec)=ao_qmckl(i_basis,1,i_elec)
               dphin(i_basis,i_elec,1)=ao_qmckl(i_basis,2,i_elec)
               dphin(i_basis,i_elec,2)=ao_qmckl(i_basis,3,i_elec)
               dphin(i_basis,i_elec,3)=ao_qmckl(i_basis,4,i_elec)
               d2phin(i_basis,i_elec)=ao_qmckl(i_basis,5,i_elec)
            enddo
         enddo


!     ! deallocate QMCkl champ arrays
         if(allocated(ao_qmckl)) deallocate(ao_qmckl)
         if(allocated(ao_vgl_qmckl)) deallocate(ao_vgl_qmckl)




!     Vectorization dependent code selection
#ifdef VECTORIZATION
!     Following loop changed for better vectorization AVX512/AVX2
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
!     nwftype endif


#else
!     keep the old localization code if no vectorization instructions available

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
!     nwftype endif


#endif
!     vectorization endif




      endif
!nedif iperiodic


#else
!     else qmckl

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

!     Vectorization dependent code selection
#ifdef VECTORIZATION
!     Following loop changed for better vectorization AVX512/AVX2
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
!     nwftype endif


#else
!     keep the old localization code if no vectorization instructions available

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
!     nwftype endif


#endif
!     vectorization endif


#endif
!     endif qmckl usage or not

      endif


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

      implicit none

      integer :: ibasis, i, ic, ielec, j, k
      integer :: l, m, n

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

      return
      end
!------------------------------------------------------------------------------------
      subroutine orbitalse(iel,x,rvec_en,r_en,iflag)

      use phifun, only: d2phin, dphin, n0_ibasis, n0_nbasis
      use phifun, only: phin
      use multiple_geo, only: iwf
      use coefs, only: nbasis
      use slater, only: norb, coef
      use contrl_per, only: iperiodic
      use system, only: ncent_tot, nelec
      use grid3dflag, only: i3dlagorb, i3dsplorb
      use multislatern, only: ddorbn, dorbn, orbn
      use precision_kinds, only: dp
      use grid3d_orbitals, only: spline_mo, lagrange_mose
      use grid3d_orbitals, only: lagrange_mos_grade
      use basis_fns_mod, only: basis_fns
      use coefs,   only: nbasis
      use contrl_per, only: iperiodic
      use grid3d_orbitals, only: lagrange_mos_grade,lagrange_mose
      use grid3d_orbitals, only: spline_mo
      use grid3dflag, only: i3dlagorb,i3dsplorb
      use multiple_geo, only: iwf
      use multislatern, only: ddorbn,dorbn,orbn
      use phifun,  only: d2phin,dphin,n0_ibasis,n0_nbasis,phin
      use precision_kinds, only: dp
      use vmc_mod, only: nwftypeorb
      use control, only: ipr
      use contrl_file, only: ounit
      use periodic, only : n_images, ell
      use find_pimage, only: find_image_pbc

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
      use qmckl_data
#endif


      implicit none

      integer :: iel, ier, ider, iflag, iorb, m
      integer :: m0, k, j

      real(dp), dimension(3,*) :: x
      real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
      real(dp), dimension(nelec,ncent_tot) :: r_en


#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
      real(dp), allocatable :: mo_vgl_qmckl(:,:)
      integer :: rc
      integer*8 :: n8, na8, i_image, ivgl, i_basis
      character*(1024) :: err_message = ''
      real(dp), allocatable :: ao_vgl_qmckl(:,:,:)
      real(dp), allocatable :: ao_qmckl(:,:)
      real(dp), dimension(3) :: xiel
      real(dp), dimension(3,n_images) :: xqmckl
      real(dp) :: rnorm
#endif




!     get the value and gradients from the 3d-interpolated orbitals
      ier=0
!     spline interplolation
      if(i3dsplorb.ge.1) then
         do k=1,nwftypeorb
            do iorb=1,norb
               ddorbn(iorb,k)=0.0d0 ! Don't compute the laplacian
               dorbn(iorb,1,k)=1.0d0 ! compute the gradients
               dorbn(iorb,2,k)=1.0d0 ! compute the gradients
               dorbn(iorb,3,k)=1.0d0 ! compute the gradients
               call spline_mo(x(1,iel),iorb,orbn(iorb,k),dorbn(iorb,:,k),ddorbn(iorb,k),ier)
            enddo
         enddo

!     Lagrange interpolation
      elseif(i3dlagorb.ge.1) then
         do k=1,nwftypeorb
            call lagrange_mose(1,x(1,iel),orbn(:,k),ier)
            call lagrange_mos_grade(2,x(1,iel),dorbn(:,:,k),ier)
            call lagrange_mos_grade(3,x(1,iel),dorbn(:,:,k),ier)
            call lagrange_mos_grade(4,x(1,iel),dorbn(:,:,k),ier)
         enddo
      else
         ier=1
      endif

      if(ier.eq.1) then
! get basis functions for electron iel

         ider=1
         if(iflag.gt.0) ider=2

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
         if(iperiodic.eq.0) then

            rc = qmckl_get_mo_basis_mo_num(qmckl_ctx, n8)
            if (rc /= QMCKL_SUCCESS) then
               print *, 'Error getting mo_num from QMCkl'
               stop
            end if



!     Send electron coordinates to QMCkl to compute the MOs at these positions
!     rc = qmckl_set_point(qmckl_ctx, 'N', 1_8, x(:,iel), 3_8)
!! set one electron coordinates
            rc = qmckl_set_point(qmckl_ctx, 'N', 1_8, x(1:3,iel), 3_8)
            if (rc /= QMCKL_SUCCESS) then
               print *, 'Error setting electron coords orbitalse'
               call qmckl_last_error(qmckl_ctx,err_message)
               print *, trim(err_message)
               call abort()
            end if

!     !allocate mo_vlg array
!     allocate(mo_vgl_qmckl(n8, 5, 1))
            allocate(mo_vgl_qmckl(n8, 5))

!     Compute the MOs
            rc = qmckl_get_mo_basis_mo_vgl_inplace( &
                 qmckl_ctx, &
                 mo_vgl_qmckl, &
                 n8*5_8)

            k=1                 ! until state-specific orbitals can use QMCKL

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


      else
!     iperiodic else


!     original statement to replace
!     call basis_fns(iel,iel,nelec,rvec_en,r_en,ider)


!     ! setting test to verify qmckl-ao calculation

!! get number of atomic orbitals
         rc = qmckl_get_ao_basis_ao_num(qmckl_ctx, na8)
            if (rc /= QMCKL_SUCCESS) then
               print *, 'Error getting mo_num from QMCkl'
               stop
            end if


            if (nbasis.ne.na8) then
               print *, 'Error getting ao_num from QMCkl'
               stop
            end if


!     !allocate ao_vlg array
            allocate(ao_qmckl(nbasis, 5))
            ao_qmckl=0.d0

!     first for image zero
!            xiel=0.d0
            xiel=x(1:3,iel)
! applying pbc (wrapping inside the box)
            call find_image_pbc(xiel,rnorm)

!     Send electron coordinates to QMCkl to compute the MOs at these positions
!     rc = qmckl_set_point(qmckl_ctx, 'N', 1_8, x(:,iel), 3_8)
!! set one electron coordinates
            rc = qmckl_set_point(qmckl_ctx, 'N', 1_8, xiel, 3_8)
            if (rc /= QMCKL_SUCCESS) then
               print *, 'Error setting electron coords orbitalse'
               call qmckl_last_error(qmckl_ctx,err_message)
               print *, trim(err_message)
               call abort()
            end if


!computing aos zero image
            rc = qmckl_get_ao_basis_ao_vgl_inplace(qmckl_ctx, &
                 ao_qmckl, nbasis*5_8)
            if (rc /= QMCKL_SUCCESS) then
               print *, 'Error getting AOs from QMCkl zero image'
            endif


! intialize before computation just in case garbage appears
            xqmckl=0.d0
!     ! computing images distance for iel
            if(n_images.gt.0) then

!!allocate ao_vgl array for all images
               allocate(ao_vgl_qmckl(nbasis, 5, n_images))
               ao_vgl_qmckl=0.d0

               do i_image=1, n_images
                  xqmckl(1:3,i_image)=xiel(1:3)-ell(1:3,i_image)
               enddo



            ! send coordinates
               rc = qmckl_set_point(qmckl_ctx, 'N', 1_8*n_images, xqmckl, 3_8*n_images)
               if (rc /= QMCKL_SUCCESS) then
                  print *, 'Error setting electron coords orbitalse'
                  call qmckl_last_error(qmckl_ctx,err_message)
                  print *, trim(err_message)
                  call abort()
               end if


               rc = qmckl_get_ao_basis_ao_vgl_inplace(qmckl_ctx, &
                    ao_vgl_qmckl, nbasis*n_images*5_8)
               if (rc /= QMCKL_SUCCESS) then
                  print *, 'Error getting AOs from QMCkl'
               endif


!     ! summing it all over the image0
               do i_image=1, n_images
                  do ivgl=1, 5
                     do i_basis=1, nbasis
                        ao_qmckl(i_basis,ivgl)=ao_qmckl(i_basis,ivgl)+ao_vgl_qmckl(i_basis,ivgl,i_image)
                     enddo
                  enddo
               enddo



            endif

!! copying QMCkl ao's back to champ
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
!     endif nwftype


         else
!else iflag

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
!     endif nwftype

         endif
!     endif iflag


      endif
!     iperiodic endif


#else
!!else qmckl

         call basis_fns(iel,iel,nelec,rvec_en,r_en,ider)

!     Vectorization dependent code. useful for AVX512 and AVX2
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
!     else iflag

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
!endif nwtype

#else
!     Keep the localization for the non-vectorized code


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
!endif nwftype




         else
!else iflag


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
!endif nwftype


         endif
!endif iflag

#endif
!     endif vectorization

#endif
!     endif qmckl

         endif
! endif for ier


      return
      end
!------------------------------------------------------------------------------------
end module
