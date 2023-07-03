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
c Written by Cyrus Umrigar starting from Kevin Schmidt's routine
c Modified by A. Scemama

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
      use precision_kinds, only: dp
      use pw_orbitals, only: orbitals_pw
      use vmc_mod, only: nwftypeorb

#ifdef QMCKL_FOUND
      use const
      use qmckl_data                                               
#endif 
      
      implicit none

      integer :: i, ier, ider, iorb, k, m
      integer :: m0, j

      real(dp), dimension(3,*) :: x
      real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
      real(dp), dimension(nelec,ncent_tot) :: r_en
c     real(dp), dimension(nelec,nbasis) :: bhin
c     real(dp), dimension(3*nelec,nbasis) :: dbhin
c     real(dp), dimension(nelec,nbasis) :: d2bhin


#ifdef QMCKL_FOUND
      real(dp), allocatable :: mo_vgl_qmckl(:,:,:)            
      integer :: rc                                           
      integer*8 :: n8                                         
#endif 

      
      ier=1
      if(iperiodic.eq.0) then

c spline interpolation
        if(i3dsplorb.eq.2) then
          do k=1,nwftypeorb
            do i=1,nelec
            ier = 0.d0
              do iorb=1,norb+nadorb
                ddorb(iorb,i,k)=1.d0    ! compute the laplacian
                dorb(iorb,i,1,k)=1.d0   ! compute the gradients
                dorb(iorb,i,2,k)=1.d0   ! compute the gradients
                dorb(iorb,i,3,k)=1.d0   ! compute the gradients
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
! Lagrange interpolation, did not inclue multiorb here yet
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

c no 3d interpolation
      else
         
c     get basis functions for all electrons
         ider=2
         if(iforce_analy.eq.1) ider=3

#ifdef QMCKL_FOUND


! get number MO's
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
         rc = qmckl_get_mo_basis_mo_vgl_inplace(
     &        qmckl_ctx,
     &        mo_vgl_qmckl,
     &        n8*nelec*5_8)
         
         if (rc /= QMCKL_SUCCESS) then
            print *, 'Error getting MOs from QMCkl'
            stop
         end if
         

!            print*, "inside qmckl"
         
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
         
         
#else

         call basis_fns(1,nelec,nelec,rvec_en,r_en,ider)



c in alternativa al loop 26
c     do jbasis=1,nbasis
c     i=0
c     do ielec=1,nelec
c     bhin(ielec,jbasis)=phin(jbasis,ielec)
c     do l=1,3
c     i=i+1
c     dbhin(i,jbasis)=dphin(jbasis,ielec,l)
c     enddo
c     d2bhin(ielec,jbasis)=d2phin(jbasis,ielec)
c     enddo
c     enddo
c     call dgemm('n','n',  nelec,norb,nbasis,1.d0,bhin,   nelec,  coef(1,1,iwf),nbasis,0.d0,orb,   nelec)
c     call dgemm('n','n',3*nelec,norb,nbasis,1.d0,dbhin,3*nelec,  coef(1,1,iwf),nbasis,0.d0,dorb,3*nelec)
c     call dgemm('n','n',  nelec,norb,nbasis,1.d0,d2bhin, nelec,  coef(1,1,iwf),nbasis,0.d0,ddorb, nelec)
         
!     Vectorization dependent code selection
#ifdef VECTORIZATION
!     Following loop changed for better vectorization AVX512/AVX2
         if(nwftypeorb.gt.1) then
            
            do k=1,nwftypeorb
               do i=1,nelec
                  do iorb=1,norb+nadorb
                     orb(i,iorb,k)=0.d0
                     dorb(iorb,i,1,k)=0.d0
                     dorb(iorb,i,2,k)=0.d0
                     dorb(iorb,i,3,k)=0.d0
                     ddorb(iorb,i,k)=0.d0
                     do m=1,nbasis
                        orb  (  i,iorb,k)=orb  (  i,iorb,k)+coef(m,iorb,k)*phin  ( m,i)
                        dorb (iorb,i,1,k)=dorb (iorb,i,1,k)+coef(m,iorb,k)*dphin (m,i,1)
                        dorb (iorb,i,2,k)=dorb (iorb,i,2,k)+coef(m,iorb,k)*dphin (m,i,2)
                        dorb (iorb,i,3,k)=dorb (iorb,i,3,k)+coef(m,iorb,k)*dphin (m,i,3)
                        ddorb(  iorb,i,k)=ddorb(iorb,i,k)+coef(m,iorb,k)*d2phin( m,i)
                     enddo
                  enddo
               enddo
            enddo
                        
         else
            
            do i=1,nelec
               do iorb=1,norb+nadorb
                  orb(i,iorb,1)=0.d0
                  dorb(iorb,i,1,1)=0.d0
                  dorb(iorb,i,2,1)=0.d0
                  dorb(iorb,i,3,1)=0.d0
                  ddorb(iorb,i,1)=0.d0
                  do m=1,nbasis
                     orb  (  i,iorb,1)=orb  (  i,iorb,1)+coef(m,iorb,iwf)*phin  ( m,i)
                     dorb (iorb,i,1,1)=dorb (iorb,i,1,1)+coef(m,iorb,iwf)*dphin (m,i,1)
                     dorb (iorb,i,2,1)=dorb (iorb,i,2,1)+coef(m,iorb,iwf)*dphin (m,i,2)
                     dorb (iorb,i,3,1)=dorb (iorb,i,3,1)+coef(m,iorb,iwf)*dphin (m,i,3)
                     ddorb(  iorb,i,1)=ddorb(iorb,i,1)+coef(m,iorb,iwf)*d2phin( m,i)
                  enddo
               enddo
            enddo
            
            
         endif
            
         
         

         
#else
!     keep the old localization code if no vectorization instructions available

         if(nwftypeorb.gt.1) then

            do k=1,nwftypeorb
               do i=1,nelec
                  do iorb=1,norb+nadorb
                     orb(i,iorb,k)=0.d0
                     dorb(iorb,i,1,k)=0.d0
                     dorb(iorb,i,2,k)=0.d0
                     dorb(iorb,i,3,k)=0.d0
                     ddorb(iorb,i,k)=0.d0
                     do m0=1,n0_nbasis(i)
                        m=n0_ibasis(m0,i)
                        orb  (  i,iorb,k)=orb  (  i,iorb,k)+coef(m,iorb,k)*phin  ( m,i)
                        dorb (iorb,i,1,k)=dorb (iorb,i,1,k)+coef(m,iorb,k)*dphin (m,i,1)
                        dorb (iorb,i,2,k)=dorb (iorb,i,2,k)+coef(m,iorb,k)*dphin (m,i,2)
                        dorb (iorb,i,3,k)=dorb (iorb,i,3,k)+coef(m,iorb,k)*dphin (m,i,3)
                        ddorb(iorb,i,k)=ddorb(iorb,i,k)+coef(m,iorb,k)*d2phin( m,i)
                     enddo
                  enddo
               enddo
            enddo
                        
         else
            
            do i=1,nelec
               do iorb=1,norb+nadorb
                  orb(i,iorb,1)=0.d0
                  dorb(iorb,i,1,1)=0.d0
                  dorb(iorb,i,2,1)=0.d0
                  dorb(iorb,i,3,1)=0.d0
                  ddorb(iorb,i,1)=0.d0
                  do m0=1,n0_nbasis(i)
                     m=n0_ibasis(m0,i)
                     orb  (  i,iorb,1)=orb  (  i,iorb,1)+coef(m,iorb,iwf)*phin  ( m,i)
                     dorb (iorb,i,1,1)=dorb (iorb,i,1,1)+coef(m,iorb,iwf)*dphin (m,i,1)
                     dorb (iorb,i,2,1)=dorb (iorb,i,2,1)+coef(m,iorb,iwf)*dphin (m,i,2)
                     dorb (iorb,i,3,1)=dorb (iorb,i,3,1)+coef(m,iorb,iwf)*dphin (m,i,3)
                     ddorb(iorb,i,1)=ddorb(iorb,i,1)+coef(m,iorb,iwf)*d2phin( m,i)
                  enddo
               enddo
            enddo
            
         endif

         
          
#endif



#endif
!endif qmckl usage or not

       endif


       if(iforce_analy.eq.1) call da_orbitals

       else
        call orbitals_pw(x,orb,dorb,ddorb)
      endif

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
c------------------------------------------------------------------------------------

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
c------------------------------------------------------------------------------------
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
      use pw_orbitals, only: orbitals_pw_grade
      use vmc_mod, only: nwftypeorb
      use control, only: ipr
      use contrl_file, only: ounit

#ifdef QMCKL_FOUND
      use qmckl_data
#endif

      
      implicit none

      integer :: iel, ier, ider, iflag, iorb, m
      integer :: m0, k, j

      real(dp), dimension(3,*) :: x
      real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
      real(dp), dimension(nelec,ncent_tot) :: r_en


#ifdef QMCKL_FOUND
c     real(dp), allocatable :: mo_vgl_qmckl(:,:,:)
      real(dp), allocatable :: mo_vgl_qmckl(:,:)
      integer :: rc
      integer*8 :: n8
      character*(1024) :: err_message = ''
#endif


      
      if(iperiodic.eq.0) then

c get the value and gradients from the 3d-interpolated orbitals
         ier=0
c     spline interplolation
         if(i3dsplorb.ge.1) then
            do k=1,nwftypeorb
              do iorb=1,norb
                 ddorbn(iorb,k)=0.0d0   ! Don't compute the laplacian
                 dorbn(iorb,1,k)=1.0d0  ! compute the gradients
                 dorbn(iorb,2,k)=1.0d0  ! compute the gradients
                 dorbn(iorb,3,k)=1.0d0  ! compute the gradients
                 call spline_mo(x(1,iel),iorb,orbn(iorb,k),dorbn(iorb,:,k),ddorbn(iorb,k),ier)
              enddo
            enddo

c     Lagrange interpolation
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
c get basis functions for electron iel

            ider=1
            if(iflag.gt.0) ider=2
            
#ifdef QMCKL_FOUND

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

!!allocate mo_vlg array
c     allocate(mo_vgl_qmckl(n8, 5, 1))
            allocate(mo_vgl_qmckl(n8, 5))
              
!     Compute the MOs
            rc = qmckl_get_mo_basis_mo_vgl_inplace(
     &           qmckl_ctx,
     &           mo_vgl_qmckl,
     &           n8*5_8)

            k=1 ! until state-specific orbitals can use QMCKL

            if(iflag.gt.0) then
               do iorb=1,norb
c     orbn(iorb)=mo_vgl_qmckl(iorb,1,1)
c     dorbn(iorb,1)=mo_vgl_qmckl(iorb,2,1)
c     dorbn(iorb,2)=mo_vgl_qmckl(iorb,3,1)
c     dorbn(iorb,3)=mo_vgl_qmckl(iorb,4,1)
c     ddorbn(iorb)=mo_vgl_qmckl(iorb,5,1)

                  orbn(iorb,k)=mo_vgl_qmckl(iorb,1)
                  dorbn(iorb,1,k)=mo_vgl_qmckl(iorb,2)
                  dorbn(iorb,2,k)=mo_vgl_qmckl(iorb,3)
                  dorbn(iorb,3,k)=mo_vgl_qmckl(iorb,4)
                  ddorbn(iorb,k)=mo_vgl_qmckl(iorb,5)
               enddo
            else
               do iorb=1,norb
c     orbn(iorb)=mo_vgl_qmckl(iorb,1,1)
c     dorbn(iorb,1)=mo_vgl_qmckl(iorb,2,1)
c     dorbn(iorb,2)=mo_vgl_qmckl(iorb,3,1)
c     dorbn(iorb,3)=mo_vgl_qmckl(iorb,4,1)

                  orbn(iorb,k)=mo_vgl_qmckl(iorb,1)
                  dorbn(iorb,1,k)=mo_vgl_qmckl(iorb,2)
                  dorbn(iorb,2,k)=mo_vgl_qmckl(iorb,3)
                  dorbn(iorb,3,k)=mo_vgl_qmckl(iorb,4)

               enddo
            endif
           
            deallocate(mo_vgl_qmckl)
            
                      

#else
            

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

               
            else

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
                  enddo

               endif

               
              


            else


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

               
               
            endif


#endif


#endif
         endif
c endif for ier
      else
         do k=1,nwftypeorb
           call orbitals_pw_grade(iel,x(1,iel),orbn(:,k),dorbn(:,:,k),ddorbn(:,k))
         enddo
      endif

      return
      end
c------------------------------------------------------------------------------------
      end module
