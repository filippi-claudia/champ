      module nonloc_mod
      use error,   only: fatal_error
      contains

      subroutine nonloc(x,rshift,rvec_en,r_en,vpsp_det,dvpsp_dj,t_vpsp,i_vpsp)
c Written by Claudia Filippi, modified by Cyrus Umrigar and A. Scemama
      use Bloc,    only: b,bkin,b_dj
      use b_tmove, only: b_t,iskip
      use contrl_file, only: ounit, errunit
      use control, only: ipr,mode
      use deriv_nonloc, only: deriv_nonlocj_quad
      use jastrow_update, only: fso
      use m_force_analytic, only: alfgeo,iforce_analy,iuse_zmat
      use multislater, only: detiab
      use optwf_control, only: ioptjas
      use optwf_parms, only: nparmj
      use orbval,  only: nadorb
      use precision_kinds, only: dp
      use pseudo,  only: lpot,vps
      use pseudo_mod, only: MPS_QUAD
      use qua,     only: nquad,wq,xq,yq,zq
      use scale_dist_mod, only: scale_dist,scale_dist1
      use slater,  only: norb,slmi
      use system,  only: cent,iwctype,ncent,ncent_tot,nelec,nup
      use vmc_mod, only: norb_tot, nwftypeorb, nwftypejas
      use contrl_per, only: iperiodic
      use csfs,    only: nstates
      use vmc_mod, only: stoj, stoo, nbjx, bjxtoo, bjxtoj

      implicit none

      ! variables in subroutine call
      integer :: i_vpsp
      real(dp), dimension(3,*) :: x
      real(dp), dimension(nelec,ncent_tot) :: r_en
      real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
      real(dp), dimension(3,nelec,ncent_tot) :: rshift
      real(dp), dimension(2,nbjx) :: vpsp_det
      real(dp), dimension(nparmj,nbjx) :: dvpsp_dj
      real(dp), dimension(ncent_tot,MPS_QUAD,*) :: t_vpsp

      ! local variables
      integer :: i, i1, i2, iab, istate, auxy
      integer :: ic, ict, iel, index
      integer :: iorb, iparm, iq, iqq
      integer :: jc, k, l, nxquad, ndim, iwforb, iwfjas, ibjx, xo, xj
      real(dp) :: ri, term1, term2
      real(dp), parameter :: one = 1.d0
      
      ! local, allocatable arrays
      integer, allocatable :: iequad(:)
      integer, allocatable :: icquad(:)
      integer, allocatable :: iqquad(:)

      real(dp), allocatable :: costh(:)
      real(dp), allocatable :: term_radial(:)
      real(dp), allocatable :: term_radial_jas(:,:)
      real(dp), allocatable :: xquad(:,:)
      real(dp), allocatable :: det_ratio(:,:)
      real(dp), allocatable :: psij_ratio(:,:)
      real(dp), allocatable :: dpsij_ratio(:,:,:)
      real(dp), allocatable :: da_psij_ratio(:,:,:)
      real(dp), allocatable :: r_en_quad(:,:)
      real(dp), allocatable :: rvec_en_quad(:,:,:)
      real(dp), allocatable :: orbn(:,:,:)
      real(dp), allocatable :: dorbn(:,:,:,:)
      real(dp), allocatable :: da_orbn(:,:,:,:)
      real(dp), allocatable :: vjn(:,:)

      ndim = nquad*nelec*2

      ! allocating local arrays
        ! integer arrays
      if(allocated(iequad)) deallocate(iequad)
      allocate(iequad(ndim))
      if(allocated(icquad)) deallocate(icquad)
      allocate(icquad(ndim))
      if(allocated(iqquad)) deallocate(iqquad)
      allocate(iqquad(ndim))
      if(allocated(costh)) deallocate(costh)
      allocate(costh(ndim))
        ! dp real arrays
      if(allocated(term_radial)) deallocate(term_radial)
      allocate(term_radial(ndim))
      if(allocated(term_radial_jas)) deallocate(term_radial_jas)
      allocate(term_radial_jas(ndim,nwftypejas))
      if(allocated(xquad)) deallocate(xquad)
      allocate(xquad(3,ndim))
      if(allocated(det_ratio)) deallocate(det_ratio)
      allocate(det_ratio(ndim,nwftypeorb))
      if(allocated(psij_ratio)) deallocate(psij_ratio)
      allocate(psij_ratio(ndim,nwftypejas))
      if(allocated(dpsij_ratio)) deallocate(dpsij_ratio)
      allocate(dpsij_ratio(nparmj,ndim,nwftypejas))
      if(allocated(da_psij_ratio)) deallocate(da_psij_ratio)
      allocate(da_psij_ratio(3,ncent_tot,ndim))
      if(allocated(r_en_quad)) deallocate(r_en_quad)
      allocate(r_en_quad(ndim,ncent_tot))
      if(allocated(rvec_en_quad)) deallocate(rvec_en_quad)
      allocate(rvec_en_quad(3,ndim,ncent_tot))
      if(allocated(orbn)) deallocate(orbn)
      allocate(orbn(norb_tot,ndim,nwftypeorb))
      if(allocated(dorbn)) deallocate(dorbn)
      allocate(dorbn(norb_tot,ndim,3,nwftypeorb))
      if(allocated(da_orbn)) deallocate(da_orbn)
      allocate(da_orbn(3,ncent_tot,norb_tot,ndim))
      if(allocated(vjn)) deallocate(vjn)
      allocate(vjn(3,ndim))

      do k=1,nbjx
        vpsp_det(1,k)=0.d0
        vpsp_det(2,k)=0.d0
        do iparm=1,nparmj
          dvpsp_dj(iparm,k)=0.d0
        enddo
      enddo

      if(i_vpsp.gt.0)then
        i1=i_vpsp
        i2=i_vpsp
       else
        i1=1
        i2=nelec
      endif

      nxquad=0
      do i=i1,i2

        do k=1,nbjx
          do iorb=1,norb+nadorb
            b(iorb,i,k)=bkin(iorb,i,k)
          enddo
        enddo

        do ic=1,ncent
          ict=iwctype(ic)

c vps was calculated by calling getvps_tm from nonloc_pot
          iskip(i,ic)=1
          do l=1,lpot(ict)-1
            if(dabs(vps(i,ic,l)).gt.1.d-4) iskip(i,ic)=0
          enddo

c skip if non-local components are zero
          if(iskip(i,ic).eq.0) then

             ri=one/r_en(i,ic)

c loop quadrature points
             do iq=1,nquad
                
                nxquad=nxquad+1
                iequad(nxquad)=i
                icquad(nxquad)=ic
                iqquad(nxquad)=iq
                
                costh(nxquad)=rvec_en(1,i,ic)*xq(iq)+rvec_en(2,i,ic)*yq(iq)+rvec_en(3,i,ic)*zq(iq)
                costh(nxquad)=costh(nxquad)*ri
                
                if(iperiodic.eq.0) then
                   xquad(1,nxquad)=r_en(i,ic)*xq(iq)+cent(1,ic)
                   xquad(2,nxquad)=r_en(i,ic)*yq(iq)+cent(2,ic)
                   xquad(3,nxquad)=r_en(i,ic)*zq(iq)+cent(3,ic)
                else
                   xquad(1,nxquad)=r_en(i,ic)*xq(iq)+cent(1,ic)+rshift(1,i,ic)
                   xquad(2,nxquad)=r_en(i,ic)*yq(iq)+cent(2,ic)+rshift(2,i,ic)
                   xquad(3,nxquad)=r_en(i,ic)*zq(iq)+cent(3,ic)+rshift(3,i,ic)
                endif
                
                r_en_quad(nxquad,ic)=r_en(i,ic)
                call distance_quad(nxquad,ic,xquad(1,nxquad),r_en_quad,rvec_en_quad,rshift)
                
             enddo
          elseif(i_vpsp.ne.0)then
             if(index(mode,'dmc').ne.0) then
                do iq=1,nquad
                   t_vpsp(ic,iq,i)=0.d0
                   do iorb=1,norb
                      b_t(iorb,iq,ic,i)=0.d0
                   enddo
                enddo
             else
                do iq=1,nquad
                   t_vpsp(ic,iq,i)=0.d0
                enddo
             endif
c     endif iskip
          endif
          
       enddo
      enddo

      if(nxquad.eq.0) return

      do iwforb=1,nwftypeorb
        call orbitals_quad(nxquad,xquad,rvec_en_quad,r_en_quad,orbn(1,1,iwforb),
     &                   dorbn(1,1,1,iwforb),da_orbn,iwforb)
        call nonlocd_quad(nxquad,iequad,orbn(1,1,iwforb),det_ratio(1,iwforb),iwforb)
      enddo
      if(ioptjas.eq.0) then
        do iwfjas=1,nwftypejas 
            call nonlocj_quad(nxquad,xquad,iequad,x,rshift,r_en,
     &         rvec_en_quad,r_en_quad,psij_ratio(1,iwfjas),vjn,da_psij_ratio,
     &         fso(1,1,iwfjas),iwfjas)
        enddo
      else
        do iwfjas=1,nwftypejas
            call deriv_nonlocj_quad(nxquad,xquad,iequad,x,rshift,r_en,
     &         rvec_en_quad,r_en_quad,psij_ratio(1,iwfjas),dpsij_ratio(1,1,iwfjas),vjn,
     &         da_psij_ratio,iwfjas)
        enddo
      endif

      do iq=1,nxquad

        iel=iequad(iq)
        ic=icquad(iq)
        iqq=iqquad(iq)

        ict=iwctype(ic)

        iab=1
        if(iel.gt.nup) iab=2
        
        term_radial(iq)=0.d0
        do l=1,lpot(ict)-1
          term_radial(iq)=term_radial(iq)+yl0(l,costh(iq))*vps(iel,ic,l)
        enddo
        do iwfjas=1,nwftypejas
          term_radial_jas(iq,iwfjas)=term_radial(iq)*wq(iqq)*exp(psij_ratio(iq,iwfjas))
        enddo
c vpsp_det  = vnl(D_kref J)/(D_kref J)
        do ibjx=1,nbjx
          xj=bjxtoj(ibjx)
          xo=bjxtoo(ibjx)
          vpsp_det(iab,ibjx)=vpsp_det(iab,ibjx)+term_radial_jas(iq,xj)*det_ratio(iq,xo)
c pseudopotential contribution to B_eloc matrix
          do iorb=1,norb+nadorb
            b(iorb,iel,ibjx)=b(iorb,iel,ibjx)+term_radial_jas(iq,xj)*orbn(iorb,iq,xo)
          enddo
        enddo

c dvpsp_dj  = vnl(D_kref dJ)/(D_kref J)
        if(ioptjas.gt.0) then
          do ibjx=1,nbjx
            xj=bjxtoj(ibjx)
            xo=bjxtoo(ibjx)
            term2=term_radial_jas(iq,xj)*det_ratio(iq,xo)
            do iparm=1,nparmj
              dvpsp_dj(iparm,ibjx)=dvpsp_dj(iparm,ibjx)+term2*dpsij_ratio(iparm,iq,xj)
            enddo
            do iparm=1,nparmj
              do iorb=1,norb
                b_dj(iorb,iel,iparm,ibjx)=b_dj(iorb,iel,iparm,ibjx)
     &          +orbn(iorb,iq,xo)*term_radial_jas(iq,xj)*dpsij_ratio(iparm,iq,xj)
              enddo
            enddo
          enddo
        endif

c transition probabilities for Casula's moves in DMC
        do istate=1,nstates
c         swap order of this loop and if statement        
          if(index(mode,'dmc').ne.0) then
            t_vpsp(ic,iqq,iel)=det_ratio(iq,1)*term_radial_jas(iq,1)
            do iorb=1,norb
              b_t(iorb,iqq,ic,iel)=orbn(iorb,iq,istate)*term_radial_jas(iq,1)
            enddo
          endif
        enddo


      enddo !loop over nquad

      if(iforce_analy.gt.0) call compute_da_bnl(nxquad,iequad,icquad,iqquad,r_en,rvec_en,costh,term_radial_jas(1,1)
     &,orbn(1,1,1),dorbn(1,1,1,1),da_orbn,psij_ratio(1,1),vjn,da_psij_ratio)

      if(ipr.ge.4) then 
        write(ounit,'(''vpsp_det,det,r_en(1)='',100d12.4)')
     & (vpsp_det(iab,1),detiab(1,iab,1),iab=1,2),r_en(1,1)
        if(nxquad.eq.0) write(errunit,*) "warning nxquad zero", nxquad
      endif

      return
      end
c-----------------------------------------------------------------------
      function yl0(l,costh)

      use precision_kinds, only: dp
      implicit none

      integer :: l
      real(dp) :: costh, yl0

      yl0 = 0.0

      if(l.eq.1) then
        yl0=1.d0
       elseif(l.eq.2) then
        yl0=3.d0*costh
       elseif(l.eq.3) then
        yl0=2.5d0*(3*costh*costh-1)
       elseif(l.eq.4) then
        yl0=3.5d0*(5*costh*costh*costh-3*costh)
       else
        call fatal_error('YL0: implemented to l=4 only')
      endif

      return
      end

c-----------------------------------------------------------------------

      function dyl0(l,costh)
      use precision_kinds, only: dp
      implicit none

      integer :: l
      real(dp) :: costh, dyl0

      dyl0 = 0.0

      if(l.eq.1) then
        dyl0=0.d0
       elseif(l.eq.2) then
        dyl0=3.d0
       elseif(l.eq.3) then
        dyl0=15.d0*costh
       elseif(l.eq.4) then
        dyl0=10.5d0*(5*costh*costh-1)
       else
        call fatal_error('YL0: implemented to l=4 only')
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine distance_quad(iq,ic,x,r_en_quad,rvec_en_quad,rshift)

      use contrl_per, only: iperiodic
      use m_force_analytic, only: iforce_analy
      use precision_kinds, only: dp
      use pw_find_image, only: find_image4
      use qua,     only: xq,yq,zq
      use scale_dist_mod, only: scale_dist,scale_dist1
      use system,  only: cent,ncent,ncent_tot,nelec
      use qua,     only: nquad

      implicit none

      integer :: iq, ic, jc, k

      real(dp), dimension(3) :: x
      real(dp), dimension(3,nelec,ncent_tot) :: rshift
      real(dp), dimension(3,nquad*nelec*2,ncent_tot) :: rvec_en_quad
      real(dp), dimension(nquad*nelec*2,ncent_tot) :: r_en_quad
      real(dp), parameter :: one = 1.d0

      do jc=1,ncent
        do k=1,3
          rvec_en_quad(k,iq,jc)=x(k)-cent(k,jc)
        enddo

        if(jc.ne.ic) then
          if(iperiodic.eq.0) then
            r_en_quad(iq,jc)=0
            do k=1,3
              r_en_quad(iq,jc)=r_en_quad(iq,jc)+rvec_en_quad(k,iq,jc)**2
            enddo
            r_en_quad(iq,jc)=dsqrt(r_en_quad(iq,jc))
           else
            call find_image4(rshift(1,iq,jc),rvec_en_quad(1,iq,jc),r_en_quad(iq,jc))
          endif

        endif

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine orbitals_quad(nxquad,xquad,rvec_en,r_en,orbn,dorbn,da_orbn,iwforb)
c Written by Claudia Filippi, modified by Cyrus Umrigar and A. Scemama

      use m_force_analytic, only: iforce_analy
      use basis_fns_mod, only: basis_fns
      use coefs,   only: nbasis
      use contrl_per, only: iperiodic
      use grid3d_orbitals, only: lagrange_mose,spline_mo
      use grid3dflag, only: i3dlagorb,i3dsplorb
      use multiple_geo, only: iwf
      use optwf_control, only: ioptorb
      use optwf_control, only: method
      use orbval,  only: ddorb,nadorb
      use phifun,  only: dphin,n0_ibasis,n0_ic,n0_nbasis,phin
      use pw_orbitals_e, only: orbitals_pwe
      use slater,  only: coef,norb
      use sr_mod,  only: i_sr_rescale
      use system,  only: iwctype,ncent,ncent_tot,nelec
      use vmc_mod, only: norb_tot, nwftypeorb
      use qua,     only: nquad
      use precision_kinds, only: dp

      use qua,     only: nquad
      use precision_kinds, only: dp

#ifdef QMCKL_FOUND
      use qmckl_data
#endif

      
      implicit none

      integer :: ic, iel, ider, ier, ii, iq
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


#ifdef QMCKL_FOUND
      real(dp), allocatable :: mo_qmckl(:,:)
      integer :: rc
      integer*8 :: n8
#endif  
      
      
      nadorb_sav=nadorb

      if(ioptorb.eq.0.or.(method(1:3).ne.'lin'.and.i_sr_rescale.eq.0)) nadorb=0


      if(iperiodic.eq.0) then

c     get the value from the 3d-interpolated orbitals
         ier=0
         if(i3dsplorb.ge.1) then
            do iq=1,nxquad
               do iorb=1,norb+nadorb
                  ddtmp=0       ! Don't compute the laplacian
                  dtmp(1)=0     ! Don't compute the gradients
                  dtmp(2)=0     ! Don't compute the gradients
                  dtmp(3)=0     ! Don't compute the gradients
                  call spline_mo(xquad(1,iq),iorb,orbn(iorb,iq),dtmp,ddtmp,ier)
               enddo
            enddo
         elseif(i3dlagorb.ge.1) then
            do iq=1,nxquad
               call lagrange_mose(1,xquad(1,iq),orbn(iorb,iq),ier)
            enddo
         else
            ier=1
         endif

         if(ier.eq.1) then
c     get basis functions for electron iel
            ider=0
            if(iforce_analy.gt.0) ider=1

#ifdef QMCKL_FOUND

!     Send electron coordinates to QMCkl to compute the MOs at these positions
            rc = qmckl_set_point(qmckl_ctx, 'N', nxquad*1_8, xquad, nxquad*3_8)
            if (rc /= QMCKL_SUCCESS) then
                print *, 'orbitals quad Error setting electron coordinates in QMCkl'
                print *, "nxquad", nxquad
                stop
             end if
             
             rc = qmckl_get_mo_basis_mo_num(qmckl_ctx, n8)
             if (rc /= QMCKL_SUCCESS) then
                print *, 'orbitals quad Error getting mo_num from QMCkl'
                print *, "n8", n8 
                stop
             end if


             allocate(mo_qmckl(n8, nxquad))

!     Compute the MOs
             rc = qmckl_get_mo_basis_mo_value_inplace(
     &            qmckl_ctx,
     &            mo_qmckl,
     &            nxquad*n8)
             
             if (rc /= QMCKL_SUCCESS) then
                print *, 'Error orbitals quad getting MOs from QMCkl'
                stop
             end if
             
             orbn(1:norb+nadorb,1:nxquad) = mo_qmckl(1:norb+nadorb,1:nxquad)
                
             deallocate(mo_qmckl)
             
             if(iforce_analy.gt.0) then
                do iq=1,nxquad

                   do iorb=1,norb
                      do ic=1,ncent
                        do k=1,3
                           da_orbn(k,ic,iorb,iq)=0.d0
                        enddo
                     enddo
                     do m0=1,n0_nbasis(iq)
                        m=n0_ibasis(m0,iq)
                        ic=n0_ic(m0,iq)
                        ii=iwctype(ic)
                        do k=1,3
                           da_orbn(k,ic,iorb,iq)=da_orbn(k,ic,iorb,iq)-coef(m,iorb,iwf)*dphin(m,iq,k)
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
             endif

#else

          if(nwftypeorb.gt.1) iwf=1
          call basis_fns(1,nxquad,nquad*nelec*2,rvec_en,r_en,ider)
          if(nwftypeorb.gt.1) iwf=iwforb

            do iq=1,nxquad

!     Vectorization dependent code selection
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
              do m0=1,n0_nbasis(iq)
                m=n0_ibasis(m0,iq)
                ic=n0_ic(m0,iq)
                ii=iwctype(ic)
                do k=1,3
                  da_orbn(k,ic,iorb,iq)=da_orbn(k,ic,iorb,iq)-coef(m,iorb,iwf)*dphin(m,iq,k)
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
          endif
c         write(ounit,*)'orb_quad iel,ren',iel,rvec_en(1,iel,1),rvec_en(1,iel,2)
c         write(ounit,*)'orb_quad da_orb', da_orbn(1,1,1),dphin(1,iel,1)

        enddo

#endif

        endif

            enddo

#endif
            
         endif

      else
         
         call orbitals_pwe(iel,xquad,orbn)
         
      endif

      nadorb = nadorb_sav

      return
      end

c-----------------------------------------------------------------------
      subroutine nonlocd_quad(nxquad,iequad,orb,ratio,iwforb)
c Written by Claudia Filippi, modified by Cyrus Umrigar and A. Scemama

      use dorb_m,  only: iworbd
      use precision_kinds, only: dp
      use slater,  only: kref
      use system,  only: ndn,nup
      use vmc_mod, only: nmat_dim
      use slater, only: slmi
      use contrl_file,    only: ounit
      use vmc_mod, only: norb_tot

      implicit none

      integer :: nxquad, iq, iel, ikel, j, iwforb
      integer, dimension(*) :: iequad
      real(dp), dimension(*) :: ratio
      real(dp), dimension(norb_tot,*) :: orb

      do iq=1,nxquad

      iel=iequad(iq)

      if(iel.le.nup) then

        ikel=nup*(iel-1)

        ratio(iq)=0.d0
        do j=1,nup
          ratio(iq)=ratio(iq)+slmi(j+ikel,1,iwforb)*orb(iworbd(j,kref),iq)
        enddo

       else

        ikel=ndn*(iel-nup-1)

        ratio(iq)=0.d0
        do j=1,ndn
          ratio(iq)=ratio(iq)+slmi(j+ikel,2,iwforb)*orb(iworbd(j+nup,kref),iq)
        enddo

      endif

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine nonlocj_quad(nxquad,xquad,iequad,x,rshift,r_en,rvec_en_quad,r_en_quad,ratio_jn,vjn,da_psij_ratio,fso,iwfjas)

c Written by Claudia Filippi, modified by Cyrus Umrigar

      use bparm,   only: nocuspb,nspin2b
      use contrl_per, only: iperiodic
      use da_jastrow4val, only: da_j
      use jastrow, only: isc,sspinn
      use m_force_analytic, only: iforce_analy
      use nonlpsi, only: dpsianl,dpsibnl,psianl,psibnl,psinl
      use precision_kinds, only: dp
      use pw_find_image, only: find_image3
      use scale_dist_mod, only: scale_dist,scale_dist1
      use system,  only: iwctype,ncent,ncent_tot,nelec,nup
      use optwf_control, only: ioptjas
      use qua,     only: nquad

      implicit none

      integer :: i, ic, iel, ipar, isb, iwfjas
      integer :: iq, it, j, jj, k, nxquad
      integer, dimension(*) :: iequad

      real(dp) :: dd1u, dum, dumk, fsumn
      real(dp) :: rij, u
      real(dp), dimension(nelec,*) :: fso
      real(dp), dimension(3,*) :: x
      real(dp), dimension(3,*) :: xquad
      real(dp), dimension(3,nelec,ncent_tot) :: rshift
      real(dp), dimension(nelec,ncent_tot) :: r_en
      real(dp), dimension(3,nquad*nelec*2,*) :: rvec_en_quad
      real(dp), dimension(nquad*nelec*2,ncent_tot) :: r_en_quad
      real(dp), dimension(nelec,ncent_tot) :: rr_en
      real(dp), dimension(nelec,ncent_tot) :: rr_en2
      real(dp), dimension(ncent_tot) :: rr_en_quad
      real(dp), dimension(ncent_tot) :: rr_en2_quad
      real(dp), dimension(nelec,nelec) :: fsn
      real(dp), dimension(3) :: dx
      real(dp), dimension(nelec,ncent_tot) :: dd1
      real(dp), dimension(ncent_tot) :: dd1_quad
      real(dp), dimension(3,*) :: vjn
      real(dp), dimension(*) :: ratio_jn
      real(dp), dimension(3,ncent_tot,*) :: da_psij_ratio
      real(dp), parameter :: half = .5d0

      if(iforce_analy.eq.0) then
        do ic=1,ncent
          do i=1,nelec
            call scale_dist(r_en(i,ic),rr_en(i,ic),1)
            call scale_dist(r_en(i,ic),rr_en2(i,ic),2)
          enddo
        enddo

       else
        do ic=1,ncent
          do i=1,nelec
            call scale_dist1(r_en(i,ic),rr_en(i,ic),dd1(i,ic),1)
cJF added to see what happens --> gives same as iforce_analy = 0
c           call scale_dist(r_en(i,ic),rr_en2(i,ic),2)
            if(ioptjas.gt.0) call scale_dist(r_en(i,ic),rr_en2(i,ic),2)
          enddo
        enddo
      endif

      do iq=1,nxquad

      iel=iequad(iq)

      if(iforce_analy.eq.0) then
        do ic=1,ncent
          call scale_dist(r_en_quad(iq,ic),rr_en_quad(ic),1)
          call scale_dist(r_en_quad(iq,ic),rr_en2_quad(ic),2)
        enddo
       else
        do ic=1,ncent
          call scale_dist1(r_en_quad(iq,ic),rr_en_quad(ic),dd1_quad(ic),1)
cJF added to see what happens --> gives same as iforce_analy = 0
c         call scale_dist(r_en_quad(iq,ic),rr_en2_quad(ic),2)
          if(ioptjas.gt.0) call scale_dist(r_en_quad(iq,ic),rr_en2_quad(ic),2)
        enddo
      endif

      fsumn=0
      do k=1,3
         vjn(k,iq)=0.d0
      enddo

      if (nelec.lt.2) goto 47

      do jj=1,nelec

        if(jj.eq.iel) goto 45
        if(jj.lt.iel) then
          i=iel
          j=jj
         else
          i=jj
          j=iel
        endif

        sspinn=1
        ipar=0
        isb=1
        if(i.le.nup .or. j.gt.nup) then
          if(nspin2b.eq.2) then
            isb=2
           elseif(nocuspb.eq.0) then
            sspinn=half
          endif
          ipar=1
        endif

        do k=1,3
          dx(k)=x(k,jj)-xquad(k,iq)
        enddo

        if(iperiodic.eq.0) then
          rij=0
          do k=1,3
            rij=rij+dx(k)**2
          enddo
          rij=dsqrt(rij)
         else
          call find_image3(dx,rij)
        endif

c e-e terms
        if(iforce_analy.eq.0) then
          call scale_dist(rij,u,1)
         else
          call scale_dist1(rij,u,dd1u,1)
          dum=dpsibnl(u,isb,ipar,iwfjas)*dd1u/rij
          do k=1,3
            dumk=-dum*dx(k)
            vjn(k,iq)=vjn(k,iq)+dumk
          enddo
        endif

        fsn(i,j)=psibnl(u,isb,ipar,iwfjas)

c e-e-n terms
c The scaling is switched in psinl, so do not do it here.
      if(isc.ge.12) call scale_dist(rij,u,3)

        do ic=1,ncent
          it=iwctype(ic)
          fsn(i,j)=fsn(i,j) +
     &    psinl(u,rshift(1,i,ic),rshift(1,j,ic),rr_en2_quad(ic),rr_en2(jj,ic),it,iwfjas)
        enddo

        fsumn=fsumn+fsn(i,j)-fso(i,j)
   45 continue
      enddo

c e-n terms
   47 fsn(iel,iel)=0

      do ic=1,ncent
        it=iwctype(ic)
        fsn(iel,iel)=fsn(iel,iel)+psianl(rr_en_quad(ic),it,iwfjas)
      enddo

      fsumn=fsumn+fsn(iel,iel)-fso(iel,iel)
      ratio_jn(iq)=fsumn

      if(iforce_analy.gt.0) then

       do ic=1,ncent
        it=iwctype(ic)
        dum=dpsianl(rr_en_quad(ic),it,iwfjas)*dd1_quad(ic)/r_en_quad(iq,ic)
        do k=1,3
          dumk=dum*rvec_en_quad(k,iq,ic)
          vjn(k,iq)=vjn(k,iq)+dumk
          da_psij_ratio(k,ic,iq)=-dumk-da_j(k,iel,ic)
        enddo
       enddo

      endif

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_da_bnl(nxquad,iequad,icquad,iqquad,r_en,rvec_en,costh,term_radial
     &                                ,orbn,dorbn,da_orbn,psij_ratio,vjn,da_psij_ratio)

      use Bloc,    only: b_da
      use contrl_file,    only: ounit
      use da_pseudo, only: da_vps
      use m_force_analytic, only: iforce_analy
      use precision_kinds, only: dp
      use pseudo,  only: lpot,vps
      use qua,     only: nquad,wq,xq,yq,zq
      use slater,  only: norb
      use system,  only: iwctype,ncent,ncent_tot,nelec
      use vmc_mod, only: norb_tot

      implicit none

      integer :: i, ic, ict, iorb, nxquad
      integer :: iq, iqq, jc, k, l
      integer, dimension(nquad*nelec*2) :: iequad
      integer, dimension(nquad*nelec*2) :: icquad
      integer, dimension(nquad*nelec*2) :: iqquad

      real(dp) :: da_term_radial, db_tmp1, db_tmp2, db_tmp3
      real(dp) :: dum, r_eni, r_eni2, sav_db
      real(dp), dimension(nquad*nelec*2) :: costh
      real(dp), dimension(nelec,ncent_tot) :: r_en
      real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
      real(dp), dimension(norb_tot,nquad*nelec*2) :: orbn
      real(dp), dimension(norb_tot,nquad*nelec*2,3) :: dorbn
      real(dp), dimension(3,ncent_tot,norb_tot,nquad*nelec*2) :: da_orbn
      real(dp), dimension(nquad*nelec*2) :: psij_ratio
      real(dp), dimension(nquad*nelec*2) :: term_radial
      real(dp), dimension(3,*) :: vjn
      real(dp), dimension(3,ncent_tot,*) :: da_psij_ratio
      real(dp), dimension(3) :: term_radial_da_vps
      real(dp), parameter :: one = 1.d0

      if(iforce_analy.eq.0) return

      do iq=1,nxquad

      i=iequad(iq)
      ic=icquad(iq)
      iqq=iqquad(iq)

      ict=iwctype(ic)
c     sav_db=b_da(1,i,1,ic)

      da_term_radial=0.d0
      do k=1,3
        term_radial_da_vps(k)=0.d0
      enddo
      do l=1,lpot(ict)-1
        da_term_radial=da_term_radial+dyl0(l,costh(iq))*vps(i,ic,l)
        do k=1,3
          term_radial_da_vps(k)=term_radial_da_vps(k)+yl0(l,costh(iq))*da_vps(k,i,ic,l)
        enddo
      enddo
      da_term_radial=da_term_radial*wq(iqq)*exp(psij_ratio(iq))
      do k=1,3
        term_radial_da_vps(k)=term_radial_da_vps(k)*wq(iqq)*exp(psij_ratio(iq))
      enddo

      r_eni=1.d0/r_en(i,ic)
      r_eni2=r_eni*r_eni
      do iorb=1,norb
        b_da(1,i,iorb,ic)=b_da(1,i,iorb,ic)+term_radial_da_vps(1)*orbn(iorb,iq)
     &                   +da_term_radial*(-xq(iqq)*r_eni+costh(iq)*rvec_en(1,i,ic)*r_eni2)*orbn(iorb,iq)
        b_da(2,i,iorb,ic)=b_da(2,i,iorb,ic)+term_radial_da_vps(2)*orbn(iorb,iq)
     &                   +da_term_radial*(-yq(iqq)*r_eni+costh(iq)*rvec_en(2,i,ic)*r_eni2)*orbn(iorb,iq)
        b_da(3,i,iorb,ic)=b_da(3,i,iorb,ic)+term_radial_da_vps(3)*orbn(iorb,iq)
     &                   +da_term_radial*(-zq(iqq)*r_eni+costh(iq)*rvec_en(3,i,ic)*r_eni2)*orbn(iorb,iq)

         db_tmp1=term_radial(iq)*(dorbn(iorb,iq,1)+orbn(iorb,iq)*vjn(1,iq))
         db_tmp2=term_radial(iq)*(dorbn(iorb,iq,2)+orbn(iorb,iq)*vjn(2,iq))
         db_tmp3=term_radial(iq)*(dorbn(iorb,iq,3)+orbn(iorb,iq)*vjn(3,iq))

         dum=xq(iqq)*db_tmp1+yq(iqq)*db_tmp2+zq(iqq)*db_tmp3

         b_da(1,i,iorb,ic)=b_da(1,i,iorb,ic)-dum*rvec_en(1,i,ic)*r_eni+db_tmp1
         b_da(2,i,iorb,ic)=b_da(2,i,iorb,ic)-dum*rvec_en(2,i,ic)*r_eni+db_tmp2
         b_da(3,i,iorb,ic)=b_da(3,i,iorb,ic)-dum*rvec_en(3,i,ic)*r_eni+db_tmp3

         do jc=1,ncent
c          if(jc.ne.ic) then
             do k=1,3
               b_da(k,i,iorb,jc)=b_da(k,i,iorb,jc)+term_radial(iq)*(da_orbn(k,jc,iorb,iq)+orbn(iorb,iq)*da_psij_ratio(k,jc,iq))
             enddo
c          endif
         enddo
      enddo

      enddo

c     write(ounit,*) 'AFT',iel,ic,iq,b_da(1,iel,1,ic)-sav_db
      return
      end

      end module
