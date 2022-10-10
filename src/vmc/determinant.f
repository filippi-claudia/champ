      module determinant_mod
      interface !LAPACK interface
        DOUBLE PRECISION FUNCTION ddot(N,DX,INCX,DY,INCY)
          INTEGER incx,incy,n
          DOUBLE PRECISION dx(*),dy(*)
        END FUNCTION
      end interface
      contains
      subroutine determinant(ipass,x,rvec_en,r_en)
c Written by Cyrus Umrigar starting from Kevin Schmidt's routine
c Modified by A. Scemama

      use control, only: ipr
      use dets, only: ndet
      use system, only: ndn, nup
      use multidet, only: kref, kchange, kref_fixed
      use dorb_m, only: iworbd
      use control, only: mode

      use orbval, only: ddorb, dorb, orb
      use slater, only: d2dx2, ddx, fp, fpp, slmi
      use system, only: nelec

      use multislater, only: detiab, allocate_multislater
      use system, only: ncent_tot
      use precision_kinds, only: dp
      use contrl_file, only: ounit
      use optwf_handle_wf, only: dcopy
      use matinv_mod, only: matinv
      use orbitals_mod, only: orbitals
      use set_input_data, only: multideterminants_define
      use optorb_f_mod, only: optorb_define
      use optwf_contrl, only: ioptorb
      use coefs, only: norb
      use orbval, only: nadorb
      use vmc_mod, only: norb_tot

      implicit none

      integer :: i, iab, icheck, ii, ik
      integer :: index, ipass, ish, j
      integer :: jk, jorb, nel, newref
      real(dp), dimension(3, *) :: x
      real(dp), dimension(3, nelec, ncent_tot) :: rvec_en
      real(dp), dimension(nelec, ncent_tot) :: r_en
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: half = 0.5d0


c compute orbitals
      call orbitals(x,rvec_en,r_en)

      kchange=0
      icheck=0
  10  continue

      do iab=1,2

         if(iab.eq.1) then
            ish=0
            nel=nup
         else
            ish=nup
            nel=ndn
         endif

         call allocate_multislater() ! properly accessing array elements
         detiab(kref,iab)=1.d0

         jk=-nel
         do j=1,nel
            jorb=iworbd(j+ish,kref)

            jk=jk+nel

            call dcopy(nel,orb(1+ish,jorb),1,slmi(1+jk,iab),1)
            call dcopy(nel,dorb(jorb,1+ish:nel+ish,1),1,fp(1,j,iab),nel*3)
            call dcopy(nel,dorb(jorb,1+ish:nel+ish,2),1,fp(2,j,iab),nel*3)
            call dcopy(nel,dorb(jorb,1+ish:nel+ish,3),1,fp(3,j,iab),nel*3)
            call dcopy(nel,ddorb (jorb,1+ish:nel+ish),1,fpp (j,iab),nel)
         enddo

c     calculate the inverse transpose matrix and itsdeterminant
         if(nel.gt.0) call matinv(slmi(1,iab),nel,detiab(kref,iab))

c     loop through up spin electrons
c     take inner product of transpose inverse with derivative
c     vectors to get (1/detup)*d(detup)/dx and (1/detup)*d2(detup)/dx**2
         ik=-nel
         do i=1,nel
            ik=ik+nel
            ddx(1,i+ish)=ddot(nel,slmi(1+ik,iab),1,fp(1,1+ik,iab),3)
            ddx(2,i+ish)=ddot(nel,slmi(1+ik,iab),1,fp(2,1+ik,iab),3)
            ddx(3,i+ish)=ddot(nel,slmi(1+ik,iab),1,fp(3,1+ik,iab),3)
            d2dx2(i+ish)=ddot(nel,slmi(1+ik,iab),1,fpp( 1+ik,iab),1)
         enddo

         if(ipr.ge.4) then
            ik=-nel
            do i=1,nel
               ik=ik+nel
               write(ounit,*) 'slmi',iab,'M',(slmi(ii+ik,iab),ii=1,nel)
            enddo
         endif
      enddo

      if(ipr.ge.4) write(ounit,'(''detu,detd'',9d12.5)') detiab(kref,1),detiab(kref,2)

c     for dmc must be implemented: for each iw, must save not only kref,kref_old but also cdet etc.
      if(index(mode,'dmc').eq.0 .and. kref_fixed.eq.0) then ! allow if kref is allowed to vary
         icheck=icheck+1
         if(ndet.gt.1.and.kref.lt.ndet.and.icheck.le.10) then
            call check_detref(ipass,icheck,newref)
            if(newref.gt.0) goto 10

c reshuffling determinants just if the new kref was accepted
            if(newref.eq.0 .and. kchange.gt.0) then
               call multideterminants_define(kchange,icheck)
               if (ioptorb.ne.0) then
                  norb=norb+nadorb
                  write(ounit, *) norb
                  call optorb_define
               endif
            endif

         endif

c reshuffling determinants if the maximum number of iterations looking for kref was exhausted
         if (kchange.eq.10) then
            call multideterminants_define(kchange,icheck)
            if (ioptorb.ne.0) then
               norb=norb+nadorb
               write(ounit, *) norb
               call optorb_define
            endif
            write(ounit, *) "kref changed but it is not optimal"
         endif



      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine check_detref(ipass,icheck,iflag)

      use control, only: ipr
      use estpsi, only: detref
      use multidet, only: kref, kref_old, kchange
      use multislater, only: detiab, allocate_multislater
      use precision_kinds, only: dp
      use contrl_file, only: ounit
      use dets, only: ndet
      use multideterminant_mod, only: idiff
      use error, only: fatal_error
      implicit none

      integer :: iab, icheck, iflag, ipass
      real(dp) :: dcheck, dlogdet


      iflag=0

      if(ipass.le.2) return

      call allocate_multislater() !access elements after allocating
      do iab=1,2
        dlogdet=dlog10(dabs(detiab(kref,iab)))
c     dcheck=dabs(dlogdet-detref(iab)/ipass)
c     if(iab.eq.1.and.dcheck.gt.6) iflag=1
c     if(iab.eq.2.and.dcheck.gt.6) iflag=2
        dcheck=detref(iab)/ipass-dlogdet
        if(iab.eq.1.and.dcheck.gt.6) iflag=1
        if(iab.eq.2.and.dcheck.gt.6) iflag=2
        if(ipr.ge.2) write(ounit,*) 'check',dlogdet,detref(iab)/ipass
      enddo

      if(ipr.ge.2) write(ounit,*) 'check detref',iflag
      if(iflag.gt.0) then


c     block of code decoupled from multideterminants_define
c to change kref if the change is accepted or required
         if (kref .gt. 1 .and. icheck .eq. 1) then
            kref = 1
         endif



         if (idiff(kref_old, kref, iflag) .eq. 0) then
            kref = kref + 1
            if (kref .gt. ndet) then
               call fatal_error('MULTIDET_DEFINE: kref > ndet')
            endif
         endif

         write (ounit, *) 'kref change', iflag, kref_old, kref

         kref_old = kref

         kchange = kchange + 1


      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_bmatrices_kin

      use system, only: ncent
      use system, only: nelec
      use const, only: hb
      use da_jastrow4val, only: da_vj
      use da_orbval, only: da_d2orb, da_dorb
      use derivjas, only: g
      use optwf_contrl, only: ioptjas
      use optwf_parms, only: nparmj
      use Bloc, only: b_da
      use Bloc, only: b_dj
      use coefs, only: norb
      use Bloc, only: b
      use force_analy, only: iforce_analy
      use velocity_jastrow, only: vj
      use orbval, only: ddorb, dorb, nadorb
      use precision_kinds, only: dp
      use optwf_handle_wf, only: dcopy
      use sr_more, only: daxpy
      implicit none

      integer :: i, ic, iorb, iparm, l

      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: half = 0.5d0


      ! resize ddor and dorb if necessary
      ! call resize_matrix(ddorb, norb+nadorb, 2)
      ! call resize_matrix(b, norb+nadorb, 1)
      ! call resize_tensor(dorb, norb+nadorb, 3)

c compute kinetic contribution of B+Btilde to compute Eloc
      do i=1,nelec
        do iorb=1,norb+nadorb
          b(iorb,i)=-hb*(ddorb(iorb,i)+2*(vj(1,i)*dorb(iorb,i,1)+vj(2,i)*dorb(iorb,i,2)+vj(3,i)*dorb(iorb,i,3)))
        enddo
      enddo

c compute derivative of kinetic contribution of B+Btilde wrt jastrow parameters
      if(ioptjas.gt.0) then
        do iparm=1,nparmj
          do i=1,nelec
            do iorb=1,norb
              b_dj(iorb,i,iparm)=-2*hb*(g(1,i,iparm)*dorb(iorb,i,1)+g(2,i,iparm)*dorb(iorb,i,2)+g(3,i,iparm)*dorb(iorb,i,3))
            enddo
          enddo
        enddo
      endif

c compute derivative of kinetic contribution of B+Btilde wrt nuclear coordinates
      if(iforce_analy.eq.1) then
        do ic=1,ncent
          do iorb=1,norb
            call dcopy(3*nelec,da_d2orb(1,1,iorb,ic),1,b_da(1,1,iorb,ic),1)
          enddo
        enddo
        do ic=1,ncent
          do i=1,nelec
            do l=1,3
              call daxpy(norb,2*vj(1,i),da_dorb(l,1,i,1,ic),9*nelec,b_da(l,i,1,ic),3*nelec)
              call daxpy(norb,2*vj(2,i),da_dorb(l,2,i,1,ic),9*nelec,b_da(l,i,1,ic),3*nelec)
              call daxpy(norb,2*vj(3,i),da_dorb(l,3,i,1,ic),9*nelec,b_da(l,i,1,ic),3*nelec)
              call daxpy(norb,2*da_vj(l,1,i,ic),dorb(1:norb,i,1),1,b_da(l,i,1,ic),3*nelec)
              call daxpy(norb,2*da_vj(l,2,i,ic),dorb(1:norb,i,2),1,b_da(l,i,1,ic),3*nelec)
              call daxpy(norb,2*da_vj(l,3,i,ic),dorb(1:norb,i,3),1,b_da(l,i,1,ic),3*nelec)
              do iorb=1,norb
                b_da(l,i,iorb,ic)=-hb*b_da(l,i,iorb,ic)
              enddo
            enddo
          enddo
        enddo
      endif
c     do 10 ic=1,ncent
c       do 10 iorb=1,norb
c         do 10 i=1,nelec
c           do 10 l=1,3
c 10          db(l,i,iorb,ic)=da_d2orb(l,i,iorb,ic)+two*(
c    &           vj(1,i)*da_dorb(l,1,i,iorb,ic)
c    &          +vj(2,i)*da_dorb(l,2,i,iorb,ic)
c    &          +vj(3,i)*da_dorb(l,3,i,iorb,ic)
c    &          +da_vj(l,1,i,ic)*dorb(iorb,i,1)
c    &          +da_vj(l,2,i,ic)*dorb(iorb,i,2)
c    &          +da_vj(l,3,i,ic)*dorb(iorb,i,3))

      return
      end
      end module
