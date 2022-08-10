      module nonloc_mod
      use error,   only: fatal_error
      contains
      subroutine nonloc(x,rshift,rvec_en,r_en,vpsp_det,dvpsp_dj,t_vpsp,i_vpsp)
c Written by Claudia Filippi, modified by Cyrus Umrigar and A. Scemama
      use Bloc,    only: b,b_dj
      use b_tmove, only: b_t,iskip
      use contrl_file, only: ounit
      use control, only: ipr,mode
      use deriv_nonloc, only: deriv_nonlocj
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
      use system,  only: iwctype,ncent,ncent_tot,nelec,nup
      use vmc_mod, only: norb_tot

      implicit none

      integer :: i, i1, i2, i_vpsp, iab
      integer :: ic, ict, iel, index
      integer :: iorb, iparm, iq, is
      integer :: jc, k, l
      real(dp) :: costh
      real(dp) :: det_ratio, gives
      real(dp) :: psij_ratio, ri, see
      real(dp) :: term, term_radial
      real(dp), dimension(3,*) :: x
      real(dp), dimension(3,nelec,ncent_tot) :: rshift
      real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
      real(dp), dimension(nelec,ncent_tot) :: r_en
      real(dp), dimension(nelec,ncent_tot) :: rr_en
      real(dp), dimension(nelec,ncent_tot) :: rr_en2
      real(dp), dimension(ncent_tot) :: rr_en_sav
      real(dp), dimension(ncent_tot) :: rr_en2_sav
      real(dp), dimension(3) :: xsav
      real(dp), dimension(3,ncent_tot) :: rshift_sav
      real(dp), dimension(3,ncent_tot) :: rvec_en_sav
      real(dp), dimension(ncent_tot) :: r_en_sav
      real(dp), dimension(*) :: vpsp_det
      real(dp), dimension(*) :: dvpsp_dj
      real(dp), dimension(ncent_tot,MPS_QUAD,*) :: t_vpsp
      real(dp), dimension(nparmj) :: dpsij_ratio
      real(dp), dimension(norb_tot) :: orbn
      real(dp), dimension(norb_tot,3) :: dorbn
      real(dp), dimension(3,ncent_tot,norb_tot) :: da_orbn
      real(dp), dimension(3) :: term_radial_da_vps
      real(dp), dimension(3) :: vjn
      real(dp), dimension(3,ncent_tot) :: da_ratio_jn
      real(dp), dimension(nelec,ncent_tot) :: dd1
      real(dp), dimension(ncent_tot) :: dd1_sav
      real(dp), parameter :: one = 1.d0

      ! call resize_matrix(b, norb+nadorb, 1)

      do ic=1,ncent
cJF this is the culprit
        if(iforce_analy.eq.0) then
          do i=1,nelec
            call scale_dist(r_en(i,ic),rr_en(i,ic),1)
            call scale_dist(r_en(i,ic),rr_en2(i,ic),2)
          enddo
         else
          do i=1,nelec
            call scale_dist1(r_en(i,ic),rr_en(i,ic),dd1(i,ic),1)
cJF added to see what happens --> gives same as iforce_analy = 0
c           call scale_dist(r_en(i,ic),rr_en2(i,ic),2)
            if(ioptjas.gt.0) call scale_dist(r_en(i,ic),rr_en2(i,ic),2)
          enddo
        endif
      enddo

      vpsp_det(1)=0
      vpsp_det(2)=0
        do iparm=1,nparmj
          dvpsp_dj(iparm)=0
        enddo

      if(i_vpsp.gt.0)then
        i1=i_vpsp
        i2=i_vpsp
       else
        i1=1
        i2=nelec
      endif



      do i=i1,i2

        iab=1
        if(i.gt.nup) iab=2

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

            do k=1,3
              xsav(k)=x(k,i)
            enddo
            do jc=1,ncent
              r_en_sav(jc)=r_en(i,jc)
              rr_en_sav(jc)=rr_en(i,jc)
              rr_en2_sav(jc)=rr_en2(i,jc)
              dd1_sav(jc)=dd1(i,jc)
              do k=1,3
                rshift_sav(k,jc)=rshift(k,i,jc)
                rvec_en_sav(k,jc)=rvec_en(k,i,jc)
              enddo
            enddo

c loop quadrature points
            do iq=1,nquad

              costh=rvec_en_sav(1,ic)*xq(iq)
     &             +rvec_en_sav(2,ic)*yq(iq)
     &             +rvec_en_sav(3,ic)*zq(iq)
              costh=costh*ri

              call dist_quad(i,ic,iq,x(1,i),r_en,rvec_en,rshift,rr_en,rr_en2,dd1)

              iel=i
              call orbitals_quad(iel,x,rvec_en,r_en,orbn,dorbn,da_orbn,iforce_analy)

              call nonlocd(iel,orbn,detiab(1,1),detiab(1,2),slmi(1,1),slmi(1,2),det_ratio)
              if(ioptjas.gt.0) then
                call deriv_nonlocj(iel,x,rshift,rvec_en,r_en,rr_en,rr_en2,dd1,psij_ratio,dpsij_ratio,vjn,da_ratio_jn)
               else
                call nonlocj(iel,x,rshift,rvec_en,r_en,rr_en,rr_en2,dd1,fso,psij_ratio,vjn,da_ratio_jn)
              endif

c             write(ounit,*) 'PSI',psij_ratio
c             write(ounit,*) 'HELLO',(vjn(k),k=1,3)
c             do ll=1,ncent
c               write(ounit,*) 'HELLO',(da_ratio_jn(k,ll),k=1,3)
c             enddo

              term_radial=0.d0
              do l=1,lpot(ict)-1
                term_radial=term_radial+yl0(l,costh)*vps(i,ic,l)
              enddo
              term_radial=term_radial*wq(iq)*exp(psij_ratio)

c vpsp_det  = vnl(D_kref J)/(D_kref J)
              vpsp_det(iab)=vpsp_det(iab)+det_ratio*term_radial

c pseudopotential contribution to B_eloc matrix
               do iorb=1,norb+nadorb
                 b(iorb,i)=b(iorb,i)+orbn(iorb)*term_radial
               enddo

c dvpsp_dj  = vnl(D_kref dJ)/(D_kref J)
              if(ioptjas.gt.0) then
                term=term_radial*det_ratio
                do iparm=1,nparmj
                  dvpsp_dj(iparm)=dvpsp_dj(iparm)+term*dpsij_ratio(iparm)

                  do iorb=1,norb
                    b_dj(iorb,i,iparm)=b_dj(iorb,i,iparm)+orbn(iorb)*term_radial*dpsij_ratio(iparm)
                  enddo
                enddo
              endif

c transition probabilities for Casula's moves in DMC
              if(index(mode,'dmc').ne.0) then
                t_vpsp(ic,iq,i)=det_ratio*term_radial
                do iorb=1,norb
                  b_t(iorb,iq,ic,i)=orbn(iorb)*term_radial
                enddo
              endif

              if(iforce_analy.gt.0)
     &        call compute_da_bnl(iel,ic,ict,iq,r_en_sav,rvec_en_sav,costh,
     &                                     term_radial,orbn,dorbn,da_orbn,psij_ratio,vjn,da_ratio_jn)

c end loop quadrature points
            enddo

            do k=1,3
              x(k,i)=xsav(k)
            enddo
            do jc=1,ncent
              r_en(i,jc)=r_en_sav(jc)
              rr_en(i,jc)=rr_en_sav(jc)
              rr_en2(i,jc)=rr_en2_sav(jc)
              dd1(i,jc)=dd1_sav(jc)
              do k=1,3
                rshift(k,i,jc)=rshift_sav(k,jc)
                rvec_en(k,i,jc)=rvec_en_sav(k,jc)
              enddo
            enddo

c elseif iskip
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
c endif dmc
            endif
c endif iskip
          endif
c end loop nelec, ncent
        enddo
      enddo

      if(ipr.ge.4) write(ounit,'(''vpsp_det,det,r_en(1)='',100d12.4)')
     &,(vpsp_det(iab),detiab(1,iab),iab=1,2),r_en(1,1)

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
      subroutine dist_quad(i,ic,iq,x,r_en,rvec_en,rshift,rr_en,rr_en2,dd1)

      use contrl_per, only: iperiodic
      use m_force_analytic, only: iforce_analy
      use precision_kinds, only: dp
      use pw_find_image, only: find_image4
      use qua,     only: xq,yq,zq
      use scale_dist_mod, only: scale_dist,scale_dist1
      use system,  only: cent,ncent,ncent_tot,nelec

      implicit none

      integer :: i, ic, iq, jc, k

      real(dp), dimension(3) :: x
      real(dp), dimension(3,nelec,ncent_tot) :: rshift
      real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
      real(dp), dimension(nelec,ncent_tot) :: r_en
      real(dp), dimension(nelec,ncent_tot) :: rr_en
      real(dp), dimension(nelec,ncent_tot) :: rr_en2
      real(dp), dimension(nelec,ncent_tot) :: dd1
      real(dp), parameter :: one = 1.d0




      if(iperiodic.eq.0) then
        x(1)=r_en(i,ic)*xq(iq)+cent(1,ic)
        x(2)=r_en(i,ic)*yq(iq)+cent(2,ic)
        x(3)=r_en(i,ic)*zq(iq)+cent(3,ic)
       else
        x(1)=r_en(i,ic)*xq(iq)+cent(1,ic)+rshift(1,i,ic)
        x(2)=r_en(i,ic)*yq(iq)+cent(2,ic)+rshift(2,i,ic)
        x(3)=r_en(i,ic)*zq(iq)+cent(3,ic)+rshift(3,i,ic)
      endif

      do jc=1,ncent
        do k=1,3
          rvec_en(k,i,jc)=x(k)-cent(k,jc)
        enddo

        if(jc.ne.ic) then
          if(iperiodic.eq.0) then
            r_en(i,jc)=0
            do k=1,3
              r_en(i,jc)=r_en(i,jc)+rvec_en(k,i,jc)**2
            enddo
            r_en(i,jc)=dsqrt(r_en(i,jc))
           else
            call find_image4(rshift(1,i,jc),rvec_en(1,i,jc),r_en(i,jc))
          endif

          if(iforce_analy.eq.0) then
            call scale_dist(r_en(i,jc),rr_en(i,jc),1)
            call scale_dist(r_en(i,jc),rr_en2(i,jc),2)
           else
            call scale_dist1(r_en(i,jc),rr_en(i,jc),dd1(i,jc),1)
            call scale_dist(r_en(i,jc),rr_en2(i,jc),2)
          endif
        endif

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine orbitals_quad(iel,x,rvec_en,r_en,orbn,dorbn,da_orbn,iforce_analy)
c Written by Claudia Filippi, modified by Cyrus Umrigar and A. Scemama

      use basis_fns_mod, only: basis_fns
      use coefs,   only: nbasis
      use contrl_per, only: iperiodic
      use grid3d_orbitals, only: lagrange_mose,spline_mo
      use grid3dflag, only: i3dlagorb,i3dsplorb
      use multiple_geo, only: iwf
      use optwf_control, only: ioptorb,method
      use orbval,  only: ddorb,nadorb
      use phifun,  only: dphin,n0_ibasis,n0_ic,n0_nbasis,phin
      use precision_kinds, only: dp
      use pw_orbitals_e, only: orbitals_pwe
      use slater,  only: coef,norb
      use system,  only: iwctype,ncent,ncent_tot,nelec
      use trexio_basis_fns_mod, only: trexio_basis_fns
      use trexio_read_data, only: trexio_has_group_orbitals
      use vmc_mod, only: norb_tot
#if defined(TREXIO_FOUND)
#endif

      implicit none

      integer :: ic, iel, ider, ier, iforce_analy, ii
      integer :: iorb, k, m, m0
      integer :: nadorb_sav

      real(dp), dimension(3) :: x
      real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
      real(dp), dimension(nelec,ncent_tot) :: r_en
      real(dp), dimension(*) :: orbn
      real(dp), dimension(norb_tot, 3) :: dorbn
      real(dp), dimension(3,ncent_tot,*) :: da_orbn
      real(dp), dimension(3) :: dtmp

      nadorb_sav=nadorb

      if(ioptorb.eq.0.or.method(1:3).ne.'lin') nadorb=0

      ! call resize_tensor(coef, norb+nadorb, 2)

      if(iperiodic.eq.0) then

c get the value from the 3d-interpolated orbitals
        ier=0
        if(i3dsplorb.ge.1) then
          do iorb=1,norb+nadorb
            ddorb=0     ! Don't compute the laplacian
            dtmp(1)=0   ! Don't compute the gradients
            dtmp(2)=0   ! Don't compute the gradients
            dtmp(3)=0   ! Don't compute the gradients
            call spline_mo(x,iorb,orbn(iorb),dtmp,ddorb(iorb,iel),ier)
          enddo
         elseif(i3dlagorb.ge.1) then
          call lagrange_mose(1,x,orbn,ier)
         else
          ier=1
        endif

        if(ier.eq.1) then
c get basis functions for electron iel
          ider=0
          if(iforce_analy.gt.0) ider=1
#if defined(TREXIO_FOUND)
          if (trexio_has_group_orbitals) then
            call trexio_basis_fns(iel,iel,rvec_en,r_en,ider)
          else
            call basis_fns(iel,iel,rvec_en,r_en,ider)
          endif
#else
          call basis_fns(iel,iel,rvec_en,r_en,ider)
#endif          


! Vectorization dependent code selection
#ifdef VECTORIZATION
          ! The following loop changed for better vectorization AVX512/AVX2
          do iorb=1,norb+nadorb
             orbn(iorb)=0.d0
             do m=1,nbasis
                orbn(iorb)=orbn(iorb)+coef(m,iorb,iwf)*phin(m,iel)
             enddo
          enddo
#else
          do iorb=1,norb+nadorb
             orbn(iorb)=0.d0
             do m0=1,n0_nbasis(iel)
                m=n0_ibasis(m0,iel)
                orbn(iorb)=orbn(iorb)+coef(m,iorb,iwf)*phin(m,iel)
             enddo
          enddo
#endif

          if(iforce_analy.gt.0) then
            do iorb=1,norb
              do ic=1,ncent
                do k=1,3
                  da_orbn(k,ic,iorb)=0.d0
                enddo
              enddo
              do m0=1,n0_nbasis(iel)
                m=n0_ibasis(m0,iel)
                ic=n0_ic(m0,iel)
                ii=iwctype(ic)
                do k=1,3
                  da_orbn(k,ic,iorb)=da_orbn(k,ic,iorb)-coef(m,iorb,iwf)*dphin(m,iel,k)
                enddo
              enddo
              do k=1,3
                dorbn(iorb,k)=0.d0
              enddo
              do ic=1,ncent
                do k=1,3
                   dorbn(iorb,k)=dorbn(iorb,k)-da_orbn(k,ic,iorb)
                enddo
              enddo
            enddo
          endif
c         write(ounit,*)'orb_quad iel,ren',iel,rvec_en(1,iel,1),rvec_en(1,iel,2)
c         write(ounit,*)'orb_quad da_orb', da_orbn(1,1,1),dphin(1,iel,1)
        endif

       else

        call orbitals_pwe(iel,x,orbn)

      endif

      nadorb = nadorb_sav

      return
      end

c-----------------------------------------------------------------------
      subroutine nonlocd(iel,orb,detu,detd,slmui,slmdi,ratio)
c Written by Claudia Filippi, modified by Cyrus Umrigar and A. Scemama

      use dorb_m,  only: iworbd
      use precision_kinds, only: dp
      use slater,  only: kref
      use system,  only: ndn,nup
      use vmc_mod, only: nmat_dim

      implicit none

      integer :: iel, ikel, j
      real(dp) :: ratio
      real(dp), dimension(3) :: x
      real(dp), dimension(*) :: detu
      real(dp), dimension(*) :: detd
      real(dp), dimension(nmat_dim) :: slmui
      real(dp), dimension(nmat_dim) :: slmdi
      real(dp), dimension(*) :: orb



      if(iel.le.nup) then

        ikel=nup*(iel-1)

        ratio=0.d0
        do j=1,nup
          ratio=ratio+slmui(j+ikel)*orb(iworbd(j,kref))
        enddo

       else

        ikel=ndn*(iel-nup-1)

        ratio=0.d0
        do j=1,ndn
          ratio=ratio+slmdi(j+ikel)*orb(iworbd(j+nup,kref))
        enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine nonlocj(iel,x,rshift,rvec_en,r_en,rr_en,rr_en2,dd1,fso,ratio_jn,vjn,da_ratio_jn)
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


      implicit none

      integer :: i, ic, iel, ipar, isb
      integer :: it, j, jj, k
      real(dp) :: dd1u, dum, dumk, fsumn
      real(dp) :: ratio_jn, rij, u
      real(dp), dimension(nelec,*) :: fso
      real(dp), dimension(3,*) :: x
      real(dp), dimension(3,nelec,ncent_tot) :: rshift
      real(dp), dimension(3,nelec,*) :: rvec_en
      real(dp), dimension(nelec,ncent_tot) :: r_en
      real(dp), dimension(nelec,ncent_tot) :: rr_en
      real(dp), dimension(nelec,ncent_tot) :: rr_en2
      real(dp), dimension(nelec,nelec) :: fsn
      real(dp), dimension(3) :: dx
      real(dp), dimension(nelec,ncent_tot) :: dd1
      real(dp), dimension(3) :: vjn
      real(dp), dimension(3,ncent_tot) :: da_ratio_jn
      real(dp), parameter :: half = .5d0






      fsumn=0

       do k=1,3
         vjn(k)=0.d0
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
          dx(k)=x(k,jj)-x(k,iel)
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
          dum=dpsibnl(u,isb,ipar)*dd1u/rij
          do k=1,3
            dumk=-dum*dx(k)
            vjn(k)=vjn(k)+dumk
          enddo
        endif

        fsn(i,j)=psibnl(u,isb,ipar)

c e-e-n terms
c The scaling is switched in psinl, so do not do it here.
      if(isc.ge.12) call scale_dist(rij,u,3)

        do ic=1,ncent
          it=iwctype(ic)
          fsn(i,j)=fsn(i,j) +
     &    psinl(u,rshift(1,i,ic),rshift(1,j,ic),rr_en2(i,ic),rr_en2(j,ic),it)
        enddo
        fsumn=fsumn+fsn(i,j)-fso(i,j)
   45 continue
      enddo

c e-n terms
   47 fsn(iel,iel)=0

      do ic=1,ncent
        it=iwctype(ic)
        fsn(iel,iel)=fsn(iel,iel)+psianl(rr_en(iel,ic),it)
      enddo

      fsumn=fsumn+fsn(iel,iel)-fso(iel,iel)
      ratio_jn=fsumn

      if(iforce_analy.eq.0) return

       do ic=1,ncent
        it=iwctype(ic)
        dum=dpsianl(rr_en(iel,ic),it)*dd1(iel,ic)/r_en(iel,ic)
        do k=1,3
          dumk=dum*rvec_en(k,iel,ic)
          vjn(k)=vjn(k)+dumk
          da_ratio_jn(k,ic)=-dumk-da_j(k,iel,ic)
        enddo
       enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_da_bnl(i,ic,ict,iq,r_en_sav,rvec_en_sav,costh,
     &                                   term_radial,orbn,dorbn,da_orbn,psij_ratio,vjn,da_ratio_jn)

      use Bloc,    only: b_da
      use da_pseudo, only: da_vps
      use m_force_analytic, only: iforce_analy
      use precision_kinds, only: dp
      use pseudo,  only: lpot,vps
      use qua,     only: wq,xq,yq,zq
      use slater,  only: norb
      use system,  only: ncent,ncent_tot
      use vmc_mod, only: norb_tot

      implicit none

      integer :: i, ic, ict, iel, iorb
      integer :: iq, jc, k, l
      real(dp) :: costh, da_term_radial, db_tmp1, db_tmp2, db_tmp3
      real(dp) :: dum, psij_ratio, r_en_savi
      real(dp) :: r_en_savi2, sav_db, term_radial
      real(dp), dimension(3,ncent_tot) :: rvec_en_sav
      real(dp), dimension(ncent_tot) :: r_en_sav
      real(dp), dimension(norb_tot) :: orbn
      real(dp), dimension(norb_tot, 3) :: dorbn
      real(dp), dimension(3) :: vjn
      real(dp), dimension(3,ncent_tot,norb_tot) :: da_orbn
      real(dp), dimension(3,ncent_tot) :: da_ratio_jn
      real(dp), dimension(3) :: term_radial_da_vps
      real(dp), parameter :: one = 1.d0




      if(iforce_analy.eq.0) return

      iel=i
      sav_db=b_da(1,iel,1,ic)

      da_term_radial=0.d0
      do k=1,3
        term_radial_da_vps(k)=0.d0
      enddo

      do l=1,lpot(ict)-1
        da_term_radial=da_term_radial+dyl0(l,costh)*vps(i,ic,l)
        do k=1,3
          term_radial_da_vps(k)=term_radial_da_vps(k)+yl0(l,costh)*da_vps(k,i,ic,l)*wq(iq)*exp(psij_ratio)
        enddo
      enddo
      da_term_radial=da_term_radial*wq(iq)*exp(psij_ratio)

      r_en_savi=1.d0/r_en_sav(ic)
      r_en_savi2=r_en_savi*r_en_savi
      do iorb=1,norb
        b_da(1,iel,iorb,ic)=b_da(1,iel,iorb,ic)+term_radial_da_vps(1)*orbn(iorb)
     &                   +da_term_radial*(-xq(iq)*r_en_savi+costh*rvec_en_sav(1,ic)*r_en_savi2)*orbn(iorb)
        b_da(2,iel,iorb,ic)=b_da(2,iel,iorb,ic)+term_radial_da_vps(2)*orbn(iorb)
     &                   +da_term_radial*(-yq(iq)*r_en_savi+costh*rvec_en_sav(2,ic)*r_en_savi2)*orbn(iorb)
        b_da(3,iel,iorb,ic)=b_da(3,iel,iorb,ic)+term_radial_da_vps(3)*orbn(iorb)
     &                   +da_term_radial*(-zq(iq)*r_en_savi+costh*rvec_en_sav(3,ic)*r_en_savi2)*orbn(iorb)


         db_tmp1=term_radial*(dorbn(iorb,1)+orbn(iorb)*vjn(1))
         db_tmp2=term_radial*(dorbn(iorb,2)+orbn(iorb)*vjn(2))
         db_tmp3=term_radial*(dorbn(iorb,3)+orbn(iorb)*vjn(3))

         dum=xq(iq)*db_tmp1+yq(iq)*db_tmp2+zq(iq)*db_tmp3

         b_da(1,iel,iorb,ic)=b_da(1,iel,iorb,ic)-dum*rvec_en_sav(1,ic)*r_en_savi+db_tmp1
         b_da(2,iel,iorb,ic)=b_da(2,iel,iorb,ic)-dum*rvec_en_sav(2,ic)*r_en_savi+db_tmp2
         b_da(3,iel,iorb,ic)=b_da(3,iel,iorb,ic)-dum*rvec_en_sav(3,ic)*r_en_savi+db_tmp3

         do jc=1,ncent
c          if(jc.ne.ic) then
             do k=1,3
               b_da(k,iel,iorb,jc)=b_da(k,iel,iorb,jc)+term_radial*(da_orbn(k,jc,iorb)+orbn(iorb)*da_ratio_jn(k,jc))
             enddo
c          endif
         enddo
      enddo

c     write(ounit,*) 'AFT',iel,ic,iq,b_da(1,iel,1,ic)-sav_db
      return
      end
      end module
