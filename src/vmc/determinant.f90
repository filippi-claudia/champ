module determinant_mod
      interface !LAPACK interface
        DOUBLE PRECISION FUNCTION ddot(N,DX,INCX,DY,INCY)
          INTEGER incx,incy,n
          DOUBLE PRECISION dx(*),dy(*)
        END FUNCTION
      end interface
contains
      subroutine determinant(ipass,x,rvec_en,r_en)
! Written by Cyrus Umrigar starting from Kevin Schmidt's routine
! Modified by A. Scemama

      use contrl_file, only: ounit
      use control, only: ipr,mode
      use dorb_m,  only: iworbd
      use matinv_mod, only: matinv
      use multidet, only: kchange,kref_fixed
      use multislater, only: allocate_multislater,detiab
      use optorb_f_mod, only: optorb_define
      use optwf_control, only: ioptorb, ioptbf
      use optwf_handle_wf, only: dcopy
      use orbitals_mod, only: orbitals
      use orbval,  only: ddorb,dorb,nadorb,orb
      use precision_kinds, only: dp
      use set_input_data, only: multideterminants_define
      use slater,  only: d2dx2,ddx,kref,ndet,norb,slmi
      use system,  only: ncent_tot,ndn,nelec,nup
      use vmc_mod, only: nmat_dim, norb_tot, nwftypeorb
      use m_backflow, only: quasi_x, dquasi_dx, d2quasi_dx2, ibackflow, rvec_en_bf, r_en_bf
      use m_backflow, only: dslm, d2slm, d2orb, deriv_parm_bf, nparm_bf, dquasi_dp
      use backflow_mod, only: backflow

      implicit none

      integer, parameter :: max_ref_checks = 10
      integer :: i, iab, icheck, ii, ik, p
      integer :: index, ipass, ish, j, k, l, kk, jj, m, n
      integer :: jk, jorb, nel, spin_to_change
      real(dp) :: tmp
      real(dp), dimension(3, nelec) :: x
      real(dp), dimension(3, nelec, ncent_tot) :: rvec_en
      real(dp), dimension(nelec, ncent_tot) :: r_en
      real(dp), dimension(3, nmat_dim, 2, nwftypeorb) :: fp
      real(dp), dimension(nmat_dim, 2, nwftypeorb) :: fpp
      real(dp), allocatable :: tmp1(:,:,:), tmp2_diag(:,:,:)

! Compute orbital values and derivatives.
      if (ibackflow .gt. 0) then
          call backflow(x)
          call orbitals(quasi_x, rvec_en_bf, r_en_bf)
          ddx = 0.0d0
          d2dx2 = 0.0d0
          if (ioptbf.gt.0) then
            deriv_parm_bf = 0.0d0
          endif
      else
          call orbitals(x,rvec_en,r_en)
      endif

      kchange=0
      icheck=0
      10 continue
      call allocate_multislater()



  if (ibackflow.gt.0) then


      do iab=1,2

        if(iab.eq.1) then
          ish=0
          nel=nup
        else
          ish=nup
          nel=ndn
        endif

        allocate(tmp1(nel, nel,3))
        allocate(tmp2_diag(nel,3,3))
        

        do k=1,nwftypeorb
          detiab(kref,iab,k)=1.d0
          jk=-nel
          do j=1,nel
            jorb=iworbd(j+ish,kref)
            jk=jk+nel
! Build the reference Slater matrix and derivative buffers.
            call dcopy(nel,orb(1+ish,jorb,k),1,slmi(1+jk,iab,k),1)
            call dcopy(nel,dorb(jorb,1+ish,1,k),norb_tot,dslm(1,1+jk,iab,k),3)
            call dcopy(nel,dorb(jorb,1+ish,2,k),norb_tot,dslm(2,1+jk,iab,k),3)
            call dcopy(nel,dorb(jorb,1+ish,3,k),norb_tot,dslm(3,1+jk,iab,k),3)
            do kk=1,3
              do jj=1,3
                call dcopy(nel,d2orb(kk,jj,jorb,1+ish,k),3*3*norb_tot,d2slm (kk,jj,1+jk,iab,k),3*3)
              enddo 
            enddo
          enddo

! Invert the reference Slater matrix in place and save its determinant.
! After matinv, slmi holds the inverse matrix for this spin/orbital set.
          if(nel.gt.0) call matinv(slmi(1,iab,k),nel,detiab(kref,iab,k))


          tmp1=0.0d0
          tmp2_diag=0.0d0
          do kk = 1,3
            do l = 1, nel
              do j = 1, nel
                do m = 1, nel
                  ! B7 in Kwon
                  tmp1(j,m,kk) = tmp1(j,m,kk) + slmi((j-1)*nel+l,iab,k) * dslm(kk,(l-1)*nel+m,iab,k) 
                enddo
              enddo
            enddo
          enddo
          do jj = 1,3
            do kk = 1,3
              do l = 1, nel
                do j = 1, nel
                  ! Only the diagonal tmp2(j,j,kk,jj) contributes downstream
                  tmp2_diag(j,kk,jj) = tmp2_diag(j,kk,jj) + slmi((j-1)*nel+l,iab,k) * d2slm(kk,jj,(l-1)*nel+j,iab,k)
                enddo
              enddo
            enddo
          enddo

          if (ioptbf.gt.0) then
            do p =1,nparm_bf
              do j=1,nel
                do kk=1,3
                  deriv_parm_bf(p)=deriv_parm_bf(p) + tmp1(j,j,kk) * dquasi_dp(kk,j+ish,p)
                enddo
              enddo
            enddo
          endif

          do i=1,nelec
            do j=1,nel
              do kk=1,3
                do jj=1,3
                  ddx(kk,i,k)=ddx(kk,i,k)+tmp1(j,j,jj)*dquasi_dx(jj,j+ish,kk,i)
                enddo
              enddo
            enddo
          enddo


          do i = 1, nelec
            do ii = 1,3
              do kk = 1,3
                do jj = 1,3
                  do m = 1, nel
                    do j = 1, nel
                      d2dx2(i,k) = d2dx2(i,k) - tmp1(j,m,kk) * tmp1(m,j,jj) * dquasi_dx(kk,m+ish,ii,i) * dquasi_dx(jj,j+ish,ii,i)
                    enddo
                  enddo
                end do
              end do
            end do
          end do
          do i=1,nelec
            do j=1,nel
              do kk=1,3
                do ii=1,3
                  do jj=1,3
                    d2dx2(i,k)=d2dx2(i,k) + tmp2_diag(j,ii,jj) * dquasi_dx(ii,j+ish,kk,i) * dquasi_dx(jj,j+ish,kk,i)
                  enddo
                enddo
              enddo
            enddo
          enddo
          do i=1,nelec
            do j=1,nel
              do kk=1,3
                d2dx2(i,k) = d2dx2(i,k) + tmp1(j,j,kk) * d2quasi_dx2(kk,j+ish,i)
              enddo
            enddo
          enddo
        enddo
        deallocate(tmp1,tmp2_diag)
      enddo



      do k=1,nwftypeorb
        do i=1,nelec
          d2dx2(i,k)=d2dx2(i,k)+&
          ddx(1,i,k)*ddx(1,i,k)+&
          ddx(2,i,k)*ddx(2,i,k)+&
          ddx(3,i,k)*ddx(3,i,k)
        enddo
      enddo



  else
      do iab=1,2

        if(iab.eq.1) then
          ish=0
          nel=nup
        else
          ish=nup
          nel=ndn
        endif

        call allocate_multislater() ! properly accessing array elements

        do k=1,nwftypeorb
          detiab(kref,iab,k)=1.d0

          jk=-nel
          do j=1,nel
            jorb=iworbd(j+ish,kref)

            jk=jk+nel
! Build the reference Slater matrix and derivative buffers.
            call dcopy(nel,orb(1+ish,jorb,k),1,slmi(1+jk,iab,k),1)
            call dcopy(nel,dorb(jorb,1+ish,1,k),norb_tot,fp(1,j,iab,k),nel*3)
            call dcopy(nel,dorb(jorb,1+ish,2,k),norb_tot,fp(2,j,iab,k),nel*3)
            call dcopy(nel,dorb(jorb,1+ish,3,k),norb_tot,fp(3,j,iab,k),nel*3)
            call dcopy(nel,ddorb (jorb,1+ish,k),norb_tot,fpp(j,iab,k),nel)
          enddo

! Invert the reference Slater matrix in place and save its determinant.
! After matinv, slmi holds the inverse matrix for this spin/orbital set.
          if(nel.gt.0) call matinv(slmi(1,iab,k),nel,detiab(kref,iab,k))

! Contract the inverse matrix with the orbital derivative buffers to get
! grad(det)/det and laplacian(det)/det for each electron.
          ik=-nel
          do i=1,nel
            ik=ik+nel
            ddx(1,i+ish,k)=ddot(nel,slmi(1+ik,iab,k),1,fp(1,1+ik,iab,k),3)
            ddx(2,i+ish,k)=ddot(nel,slmi(1+ik,iab,k),1,fp(2,1+ik,iab,k),3)
            ddx(3,i+ish,k)=ddot(nel,slmi(1+ik,iab,k),1,fp(3,1+ik,iab,k),3)
            d2dx2(i+ish,k)=ddot(nel,slmi(1+ik,iab,k),1,fpp(1+ik,iab,k),1)
          enddo
        enddo
      enddo

  endif

          


      if(ipr.ge.4) then
        do k=1,nwftypeorb
          write(ounit,'(A,i4)') "Dets for orbital set ", k
          write(ounit,'(''detu,detd'',9d12.5)') detiab(kref,1,k),detiab(kref,2,k)
        enddo
      endif

! DMC would need per-walker determinant kref before allowing kref changes.
      if(index(mode,'dmc').eq.0 .and. kref_fixed.eq.0) then ! allow if kref is allowed to vary
         icheck=icheck+1
         if(ndet.gt.1.and.kref.lt.ndet.and.icheck.le.max_ref_checks) then
            call check_detref(ipass,icheck,spin_to_change)
            if(spin_to_change.gt.0) goto 10

! Reshuffle determinants only if the new kref was accepted.
            if(spin_to_change.eq.0 .and. kchange.gt.0) then
               call multideterminants_define(kchange)
               if (ioptorb.ne.0) then
                  norb=norb+nadorb
!                 write(ounit, *) norb
                  call optorb_define(kchange)
               endif
            endif

         endif

! Reshuffle determinants if the maximum number of reference checks was exhausted.
         if (kchange.eq.max_ref_checks) then
            call multideterminants_define(kchange)
            if (ioptorb.ne.0) then
               norb=norb+nadorb
!              write(ounit, *) norb
               call optorb_define(kchange)
            endif
            write(ounit, *) "kref changed but it is not optimal"
         endif

      endif

      end subroutine determinant
!-----------------------------------------------------------------------
      subroutine check_detref(ipass,icheck,spin_to_change)

      use control, only: ipr
      use error, only: fatal_error
      use estpsi, only: detref
      use contrl_file, only: ounit
      use multidet, only: k_aux,kchange,kref_old,ndetiab
      use multideterminant_mod, only: idiff
      use multislater, only: allocate_multislater,detiab
      use precision_kinds, only: dp
      use slater,  only: kref
      implicit none

      real(dp), parameter :: detref_log10_threshold = 6.d0
      integer :: iab, icheck, ipass, spin_to_change
      real(dp) :: dcheck, dlogdet

      spin_to_change=0

      if(ipass.le.2) return

      call allocate_multislater() !access elements after allocating

      do iab=1,2
        dlogdet=dlog10(dabs(detiab(kref,iab,1)))
        dcheck=detref(iab,1)/ipass-dlogdet
        if(iab.eq.1.and.dcheck.gt.detref_log10_threshold) spin_to_change=1
        if(iab.eq.2.and.dcheck.gt.detref_log10_threshold) spin_to_change=2
        if(ipr.ge.2) write(ounit,*) 'check',dlogdet,detref(iab,1)/ipass
      enddo

      if(ipr.ge.2) write(ounit,*) 'check detref',spin_to_change

! to change kref if required
      if(spin_to_change.gt.0) then

         if (kref .gt. 1 .and. icheck .eq. 1) then
            kref = 1
         endif

         if (idiff(kref_old, kref, spin_to_change) .eq. 0) then
            if (icheck .gt. ndetiab(spin_to_change)) &
                call fatal_error('MULTIDET_DEFINE: kref > ndet')
            kref=k_aux(icheck,spin_to_change)
         endif

         kref_old = kref
         kchange = kchange + 1
      endif

      end subroutine check_detref
!-----------------------------------------------------------------------
      subroutine compute_bmatrices_kin

      use Bloc, only: b_da, b_dj, bkin
      use constants, only: hb
      use da_jastrow, only: da_vj
      use da_orbval, only: da_d2orb,da_dorb
      use derivjas, only: g
      use m_force_analytic, only: iforce_analy
      use optwf_control, only: ioptjas
      use optwf_handle_wf, only: dcopy
      use optwf_parms, only: nparmj
      use orbval,  only: ddorb,dorb,nadorb
      use precision_kinds, only: dp
      use slater,  only: norb
      use sr_more, only: daxpy
      use system, only: ncent, nelec
      use velocity_jastrow, only: vj
      use vmc_mod, only: bjxtoj, bjxtoo, nbjx
      implicit none

      integer :: i, ic, iorb, iparm, l, k, xo, xj


! resize ddorb and dorb if necessary
! call resize_matrix(ddorb, norb+nadorb, 2)
! call resize_matrix(b, norb+nadorb, 1)
! call resize_tensor(dorb, norb+nadorb, 3)

      do k=1,nbjx
        xo=bjxtoo(k)
        xj=bjxtoj(k)
! compute kinetic contribution of B+Btilde to compute Eloc
        do i=1,nelec
          do iorb=1,norb+nadorb
            bkin(iorb,i,k)=-hb*(ddorb(iorb,i,xo)+2*(vj(1,i,xj)*dorb(iorb,i,1,xo) &
            +vj(2,i,xj)*dorb(iorb,i,2,xo)+vj(3,i,xj)*dorb(iorb,i,3,xo)))
          enddo
        enddo
      enddo
! compute derivative of kinetic contribution of B+Btilde wrt jastrow parameters
      if(ioptjas.gt.0) then
        do k=1,nbjx
          xo=bjxtoo(k)
          xj=bjxtoj(k)
          do iparm=1,nparmj
            do i=1,nelec
              do iorb=1,norb
                b_dj(iorb,i,iparm,k)=-2*hb*(g(1,i,iparm,xj)*dorb(iorb,i,1,xo) &
                                           +g(2,i,iparm,xj)*dorb(iorb,i,2,xo) &
                                           +g(3,i,iparm,xj)*dorb(iorb,i,3,xo))
              enddo
            enddo
          enddo
        enddo
      endif

! compute derivative of kinetic contribution of B+Btilde wrt nuclear coordinates
      if(iforce_analy.eq.1) then
        do ic=1,ncent
          do iorb=1,norb
            call dcopy(3*nelec,da_d2orb(1,1,iorb,ic),1,b_da(1,1,iorb,ic),1)
            !write(ounit, *), 'da_d2orb', da_d2orb(1,1,iorb,ic)
          enddo
        enddo
        do ic=1,ncent
          do i=1,nelec
            do l=1,3
              call daxpy(norb,2*vj(1,i,1),da_dorb(l,1,i,1,ic),9*nelec,b_da(l,i,1,ic),3*nelec)
              call daxpy(norb,2*vj(2,i,1),da_dorb(l,2,i,1,ic),9*nelec,b_da(l,i,1,ic),3*nelec)
              call daxpy(norb,2*vj(3,i,1),da_dorb(l,3,i,1,ic),9*nelec,b_da(l,i,1,ic),3*nelec)
              call daxpy(norb,2*da_vj(l,1,i,ic),dorb(1:norb,i,1,1),1,b_da(l,i,1,ic),3*nelec)
              call daxpy(norb,2*da_vj(l,2,i,ic),dorb(1:norb,i,2,1),1,b_da(l,i,1,ic),3*nelec)
              call daxpy(norb,2*da_vj(l,3,i,ic),dorb(1:norb,i,3,1),1,b_da(l,i,1,ic),3*nelec)
              !write(ounit, *), 'vj', vj(1,i,1), vj(2,i,1)
              !write(ounit, *), 'da_dorb', da_dorb(l,1,i,1,ic), da_dorb(l,2,i,1,ic)
              !write(ounit, *), 'da_vj', da_vj(l,1,i,ic), da_vj(l,2,i,ic)
              !write(ounit, *), 'dorb', dorb(1,i,1,1), dorb(1,i,2,1)
              do iorb=1,norb
                b_da(l,i,iorb,ic)=-hb*b_da(l,i,iorb,ic)
              enddo
            enddo
          enddo
        enddo
      endif
!     do 10 ic=1,ncent
!       do 10 iorb=1,norb
!         do 10 i=1,nelec
!           do 10 l=1,3
! 10          db(l,i,iorb,ic)=da_d2orb(l,i,iorb,ic)+two*(
!    &           vj(1,i)*da_dorb(l,1,i,iorb,ic)
!    &          +vj(2,i)*da_dorb(l,2,i,iorb,ic)
!    &          +vj(3,i)*da_dorb(l,3,i,iorb,ic)
!    &          +da_vj(l,1,i,ic)*dorb(iorb,i,1)
!    &          +da_vj(l,2,i,ic)*dorb(iorb,i,2)
!    &          +da_vj(l,3,i,ic)*dorb(iorb,i,3))

      end subroutine compute_bmatrices_kin
end module determinant_mod
