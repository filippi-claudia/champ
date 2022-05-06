      subroutine nonloc(x,rshift,rvec_en,r_en,vpsp_det,dvpsp_dj,t_vpsp,i_vpsp)
c     Written by Claudia Filippi, modified by Cyrus Umrigar, A. Scemama and RLPB
      use pseudo_mod, only: MPS_QUAD
      use optjas, only: MPARMJ
      use vmc_mod, only: MELEC, MORB, MDET, MCENT
      use vmc_mod, only: MMAT_DIM
      use mstates_mod, only: MSTATES
      use atom, only: iwctype, ncent
      use const, only: nelec, ipr
      use elec, only: nup
      use jaso, only: fso
      use optwf_contrl, only: ioptjas
      use optwf_parms, only: nparmj
      use Bloc, only: b_dj
      use coefs, only: norb
      use csfs, only: nstates
      use contr3, only: mode
      use Bloc, only: b
      use force_analy, only: iforce_analy, iuse_zmat, alfgeo
      use pseudo, only: lpot, vps
      use b_tmove, only: b_t, iskip
      use Bloc, only: b
      use force_analy, only: iforce_analy, iuse_zmat, alfgeo
      use pseudo, only: lpot, vps
      use b_tmove, only: b_t, iskip
      use qua, only: nquad, wq, xq, yq, zq
      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb
      use slater, only: d2dx2, ddx, fp, fpp, slmi
      use multislater, only: detiab

      implicit real*8(a-h,o-z)

      parameter (one=1.d0)
      dimension x(3,*),rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
      dimension rr_en(MELEC,MCENT),rr_en2(MELEC,MCENT),rr_en_sav(MCENT),rr_en2_sav(MCENT)
     &     ,xsav(3),rshift_sav(3,MCENT),rvec_en_sav(3,MCENT),r_en_sav(MCENT)
      dimension vpsp_det(2,MSTATES),dvpsp_dj(MPARMJ,MSTATES),t_vpsp(MCENT,MPS_QUAD,*)
      dimension dpsij_ratio(MPARMJ)
      dimension orbn(MORB,MSTATES),dorbn(3,MORB,MSTATES),
     &     da_orbn(3,MCENT,MORB,MSTATES),term_radial_da_vps(3)
      dimension vjn(3),da_ratio_jn(3,MCENT),dd1(MELEC,MCENT),dd1_sav(MCENT)

      do ic=1,ncent
         if(iforce_analy.eq.0) then
            do i=1,nelec
               call scale_dist(r_en(i,ic),rr_en(i,ic),1)
               call scale_dist(r_en(i,ic),rr_en2(i,ic),2)
            enddo
         else
            do i=1,nelec
               call scale_dist1(r_en(i,ic),rr_en(i,ic),dd1(i,ic),1)
               if(ioptjas.gt.0) call scale_dist(r_en(i,ic),rr_en2(i,ic),2)
            enddo
         endif
      enddo

      vpsp_det(:,:)=0.0d0
      dvpsp_dj(:,:)=0.0d0

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
c     vps was calculated by calling getvps_tm from nonloc_pot
            iskip(i,ic)=1
            do l=1,lpot(ict)-1
               if(dabs(vps(i,ic,l)).gt.1.d-4) iskip(i,ic)=0
            enddo
c     skip if non-local components are zero
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
                  do  k=1,3
                     rshift_sav(k,jc)=rshift(k,i,jc)
                     rvec_en_sav(k,jc)=rvec_en(k,i,jc)
                  enddo
               enddo
c     loop quadrature points
               do iq=1,nquad
                  costh=rvec_en_sav(1,ic)*xq(iq)
     &                 +rvec_en_sav(2,ic)*yq(iq)
     &                 +rvec_en_sav(3,ic)*zq(iq)
                  costh=costh*ri

                  call dist_quad(i,ic,iq,x(1,i),r_en,rvec_en,rshift,rr_en,rr_en2,dd1)

                  iel=i
                  do istate=1,nstates
                     call orbitals_quad(iel,x,rvec_en,r_en,
     &                    orbn(1,istate),dorbn(1,1,istate),
     &                    da_orbn(1,1,1,istate),iforce_analy,istate)
                     call nonlocd(iel,orbn(1,istate),detiab(1,1,istate),detiab(1,2,istate),
     &                    slmi(1,1,istate),slmi(1,2,istate),det_ratio)

                     if(ioptjas.gt.0) then
                        call deriv_nonlocj(iel,x,rshift,rvec_en,r_en,rr_en,
     &                       rr_en2,dd1,psij_ratio,dpsij_ratio,vjn,da_ratio_jn,istate)
                     else
                        call nonlocj(iel,x,rshift,rvec_en,r_en,rr_en,
     &                       rr_en2,dd1,fso(:,:,istate),psij_ratio,vjn,da_ratio_jn,istate)
                     endif

                     term_radial=0.0d0
                     do l=1,lpot(ict)-1
                        term_radial=term_radial+yl0(l,costh)*vps(i,ic,l)
                     enddo
                     term_radial=term_radial*wq(iq)*exp(psij_ratio)
                     
c     vpsp_det  = vnl(D_kref J)/(D_kref J)
                     vpsp_det(iab,istate)=vpsp_det(iab,istate)+det_ratio*term_radial

c     pseudopotential contribution to B_eloc matrix
                     do iorb=1,norb+nadorb
                        b(iorb,i,istate)=b(iorb,i,istate)+orbn(iorb,istate)*term_radial
                     enddo

c     dvpsp_dj  = vnl(D_kref dJ)/(D_kref J)
                     if(ioptjas.gt.0) then
                        term=term_radial*det_ratio
                        do iparm=1,nparmj
                           dvpsp_dj(iparm,istate)=dvpsp_dj(iparm,istate)+term*dpsij_ratio(iparm)
                           do iorb=1,norb
                              b_dj(iorb,i,iparm,istate)=b_dj(iorb,i,iparm,istate)
     &                     +orbn(iorb,istate)*term_radial*dpsij_ratio(iparm)
                           enddo
                        enddo
                     endif

c     transition probabilities for Casula's moves in DMC
                     if(index(mode,'dmc').ne.0) then
                        t_vpsp(ic,iq,i)=det_ratio*term_radial 
                        do iorb=1,norb
                           b_t(iorb,iq,ic,i,istate)=orbn(iorb,istate)*term_radial
                        enddo
                     endif

                     if(iforce_analy.gt.0) then
                        call compute_da_bnl(iel,ic,ict,iq,r_en_sav,rvec_en_sav,costh,
     &                       term_radial,orbn(1,istate),dorbn(1,1,istate),
     &                       da_orbn(1,1,1,istate),psij_ratio,vjn,da_ratio_jn,istate)
                     endif
                  enddo
                  
c     end loop quadrature points
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
c     elseif iskip
            elseif(i_vpsp.ne.0)then
               do istate=1,nstates
                  do iq=1,nquad
                     t_vpsp(ic,iq,i)=0.0d0
                     do iorb=1,norb
                        b_t(iorb,iq,ic,i,istate)=0.0d0
                     enddo
                  enddo
               enddo
c     endif iskip
            endif
c     end loop nelec, ncent
         enddo
      enddo

      if(ipr.ge.4) then
         do istate=1,nstates
            write(6,*) "STATE", istate
            write(6,'(''vpsp_det,det,r_en(1)='',100d12.4)')
     &           ,(vpsp_det(iab,istate),detiab(1,iab,istate),iab=1,2),r_en(1,1)
         enddo
      endif

      end subroutine

c-----------------------------------------------------------------------

      function yl0(l,costh)

      implicit real*8(a-h,o-z)

      if(l.eq.1) then
         yl0=1.0d0
      elseif(l.eq.2) then
         yl0=3.0d0*costh
      elseif(l.eq.3) then
         yl0=2.5d0*(3*costh*costh-1)
      else
         call fatal_error('YL0: implemented to l=3 only')
      endif

      end function

c-----------------------------------------------------------------------

      function dyl0(l,costh)

      implicit real*8(a-h,o-z)

      if(l.eq.1) then
         dyl0=0.0d0
      elseif(l.eq.2) then
         dyl0=3.0d0
      elseif(l.eq.3) then
         dyl0=15.0d0*costh
      else
         call fatal_error('YL0: implemented to l=3 only')
      endif

      end function

c-----------------------------------------------------------------------

      subroutine dist_quad(i,ic,iq,x,r_en,rvec_en,rshift,rr_en,rr_en2,dd1)
      use vmc_mod, only: MELEC, MCENT
      use atom, only: cent, ncent
      use contrl_per, only: iperiodic
      use force_analy, only: iforce_analy
      use qua, only: xq, yq, zq

      implicit real*8(a-h,o-z)

      parameter (one=1.d0)
      dimension x(3),rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
      dimension rr_en(MELEC,MCENT),rr_en2(MELEC,MCENT),dd1(MELEC,MCENT)

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

      end subroutine

c-----------------------------------------------------------------------

      subroutine orbitals_quad(iel,x,rvec_en,r_en,orbn,dorbn,da_orbn,iforce_analy,istate)
c     Written by Claudia Filippi, modified by Cyrus Umrigar and A. Scemama
      use vmc_mod, only: MELEC, MORB, MCENT
      use atom, only: iwctype, ncent
      use phifun, only: dphin, n0_ibasis, n0_ic, n0_nbasis
      use phifun, only: phin
      use wfsec, only: iwf
      use coefs, only: coef, nbasis, norb
      use contrl_per, only: iperiodic
      use grid3dflag, only: i3dlagorb, i3dsplorb
      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb

      implicit real*8(a-h,o-z)

      dimension x(3),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
      dimension orbn(*),dorbn(3,*),da_orbn(3,MCENT,*),dtmp(3)

      if(iperiodic.eq.0) then
c     get the value from the 3d-interpolated orbitals
         ier=0
         if(i3dsplorb.ge.1) then
            do iorb=1,norb+nadorb
               ddorb=0          ! Don't compute the laplacian
               dtmp(1)=0        ! Don't compute the gradients
               dtmp(2)=0        ! Don't compute the gradients
               dtmp(3)=0        ! Don't compute the gradients
               call spline_mo(x,iorb,orbn(iorb),dtmp,ddorb,ier)
            enddo
         elseif(i3dlagorb.ge.1) then
            call lagrange_mose(1,x,orbn,ier)
         else
            ier=1
         endif 
         if(ier.eq.1) then
c     get basis functions for electron iel
            call basis_fnse_v(iel,rvec_en,r_en)

            do iorb=1,norb+nadorb
               orbn(iorb)=0.d0
               do  m0=1,n0_nbasis(iel)
                  m=n0_ibasis(m0,iel)
                  orbn(iorb)=orbn(iorb)+coef(m,iorb,istate,iwf)*phin(m,iel)
               enddo
            enddo
            if(iforce_analy.gt.0) then
               do iorb=1,norb
                  do  ic=1,ncent
                     do  k=1,3
                        da_orbn(k,ic,iorb)=0.0d0
                     enddo
                  enddo
                  do  m0=1,n0_nbasis(iel)
                     m=n0_ibasis(m0,iel)
                     ic=n0_ic(m0,iel)
                     ii=iwctype(ic)
                     do  k=1,3
                        da_orbn(k,ic,iorb)=da_orbn(k,ic,iorb)
     &                       -coef(m,iorb,istate,iwf)*dphin(k,m,iel)
                     enddo
                  enddo
                  do  k=1,3
                     dorbn(k,iorb)=0.0d0
                  enddo
                  do ic=1,ncent
                     do k=1,3
                        dorbn(k,iorb)=dorbn(k,iorb)-da_orbn(k,ic,iorb)
                     enddo
                  enddo
               enddo
            endif
         endif
      else
         call orbitals_pwe(iel,x,orbn)
      endif

      end subroutine

c-----------------------------------------------------------------------

      subroutine nonlocd(iel,orb,detu,detd,slmui,slmdi,ratio)
c     Written by Claudia Filippi, modified by Cyrus Umrigar and A. Scemama
      use vmc_mod, only: MMAT_DIM
      use elec, only: ndn, nup
      use multidet, only: kref
      use dorb_m, only: iworbd

      implicit real*8(a-h,o-z)

      dimension x(3)
      dimension detu(*),detd(*),slmui(MMAT_DIM),slmdi(MMAT_DIM)
      dimension orb(*),dorb(3)

      if(iel.le.nup) then
         ikel=nup*(iel-1)
         ratio=0.0d0
         do j=1,nup
            ratio=ratio+slmui(j+ikel)*orb(iworbd(j,kref))
         enddo
      else
         ikel=ndn*(iel-nup-1)
         ratio=0.0d0
         do j=1,ndn
            ratio=ratio+slmdi(j+ikel)*orb(iworbd(j+nup,kref))
         enddo
      endif

      end subroutine

c-----------------------------------------------------------------------

      subroutine nonlocj(iel,x,rshift,rvec_en,r_en,rr_en,
     &rr_en2,dd1,fso,ratio_jn,vjn,da_ratio_jn,istate)
c     Written by Claudia Filippi, modified by Cyrus Umrigar
      use vmc_mod, only: MELEC, MCENT
      use atom, only: iwctype, ncent
      use jaspar, only: sspinn, is
      use const, only: nelec
      use da_jastrow4val, only: da_j
      use elec, only: nup
      use bparm, only: nocuspb, nspin2b
      use contr2, only: isc
      use contrl_per, only: iperiodic
      use force_analy, only: iforce_analy

      implicit real*8(a-h,o-z)

      parameter (half=.5d0)
      dimension fso(MELEC,*)
      dimension x(3,*),rshift(3,MELEC,MCENT),rvec_en(3,MELEC,*)
      dimension r_en(MELEC,MCENT),rr_en(MELEC,MCENT),rr_en2(MELEC,MCENT)
     &     ,fsn(MELEC,MELEC),dx(3),dd1(MELEC,MCENT),vjn(3),da_ratio_jn(3,MCENT)

      fsumn=0.0d0
      do k=1,3
         vjn(k)=0.d0
      enddo

      if (nelec.lt.2) goto 47
      do jj=1,nelec
         if(jj.eq.iel) cycle
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

c     e-e terms
         if(iforce_analy.eq.0) then
            call scale_dist(rij,u,1)
         else
            call scale_dist1(rij,u,dd1u,1)
            dum=dpsibnl(u,isb,ipar,istate)*dd1u/rij
            do k=1,3
               dumk=-dum*dx(k)
               vjn(k)=vjn(k)+dumk
            enddo
         endif
         fsn(i,j)=psibnl(u,isb,ipar,istate)

c     e-e-n terms
c     The scaling is switched in psinl, so do not do it here.
         if(isc.ge.12) call scale_dist(rij,u,3)

         do ic=1,ncent
            it=iwctype(ic)
            fsn(i,j)=fsn(i,j) +
     &               psinl(u,rshift(1,i,ic),rshift(1,j,ic),rr_en2(i,ic),rr_en2(j,ic),it,istate)
         enddo
         fsumn=fsumn+fsn(i,j)-fso(i,j)
      enddo

c     e-n terms
 47   fsn(iel,iel)=0.0d0

      do ic=1,ncent
         it=iwctype(ic)
         fsn(iel,iel)=fsn(iel,iel)+psianl(rr_en(iel,ic),it,istate)
      enddo

      fsumn=fsumn+fsn(iel,iel)-fso(iel,iel)
      ratio_jn=fsumn

      if(iforce_analy.eq.0) return

      do  ic=1,ncent
         it=iwctype(ic)
         dum=dpsianl(rr_en(iel,ic),it,istate)*dd1(iel,ic)/r_en(iel,ic)
         do  k=1,3
            dumk=dum*rvec_en(k,iel,ic)
            vjn(k)=vjn(k)+dumk
            da_ratio_jn(k,ic)=-dumk-da_j(k,iel,ic,istate)
         enddo
      enddo

      end subroutine

c-----------------------------------------------------------------------

      subroutine compute_da_bnl(i,ic,ict,iq,r_en_sav,rvec_en_sav,costh,
     &     term_radial,orbn,dorbn,da_orbn,psij_ratio,vjn,da_ratio_jn,istate)

      use vmc_mod, only: MORB, MCENT
      use atom, only: ncent
      use Bloc, only: b_da
      use coefs, only: norb
      use force_analy, only: iforce_analy
      use pseudo, only: lpot, vps
      use da_pseudo, only: da_vps
      use qua, only: wq, xq, yq, zq

      implicit real*8(a-h,o-z)

      parameter (one=1.d0)
      dimension rvec_en_sav(3,MCENT),r_en_sav(MCENT)
      dimension orbn(MORB),dorbn(3,MORB),vjn(3)
      dimension da_orbn(3,MCENT,MORB),da_ratio_jn(3,MCENT)
      dimension term_radial_da_vps(3)

      if(iforce_analy.eq.0) return

      iel=i
      sav_db=b_da(1,iel,1,ic,istate)

      da_term_radial=0.0d0
      do k=1,3
         term_radial_da_vps(k)=0.d0
      enddo

      do l=1,lpot(ict)-1
         da_term_radial=da_term_radial+dyl0(l,costh)*vps(i,ic,l)
         do k=1,3
            term_radial_da_vps(k)=term_radial_da_vps(k)
     &           +yl0(l,costh)*da_vps(k,i,ic,l)*wq(iq)*exp(psij_ratio)
         enddo
      enddo
      da_term_radial=da_term_radial*wq(iq)*exp(psij_ratio)

      r_en_savi=1.0d0/r_en_sav(ic)
      r_en_savi2=r_en_savi*r_en_savi
      do iorb=1,norb
         b_da(1,iel,iorb,ic,istate)=b_da(1,iel,iorb,ic,istate)+term_radial_da_vps(1)*orbn(iorb)
     &        +da_term_radial*(-xq(iq)*r_en_savi
     &        +costh*rvec_en_sav(1,ic)*r_en_savi2)*orbn(iorb)
         b_da(2,iel,iorb,ic,istate)=b_da(2,iel,iorb,ic,istate)+term_radial_da_vps(2)*orbn(iorb)
     &        +da_term_radial*(-yq(iq)*r_en_savi
     &        +costh*rvec_en_sav(2,ic)*r_en_savi2)*orbn(iorb)
         b_da(3,iel,iorb,ic,istate)=b_da(3,iel,iorb,ic,istate)+term_radial_da_vps(3)*orbn(iorb)
     &        +da_term_radial*(-zq(iq)*r_en_savi
     &        +costh*rvec_en_sav(3,ic)*r_en_savi2)*orbn(iorb)
         
         db_tmp1=term_radial*(dorbn(1,iorb)+orbn(iorb)*vjn(1))
         db_tmp2=term_radial*(dorbn(2,iorb)+orbn(iorb)*vjn(2))
         db_tmp3=term_radial*(dorbn(3,iorb)+orbn(iorb)*vjn(3))

         dum=xq(iq)*db_tmp1+yq(iq)*db_tmp2+zq(iq)*db_tmp3

         b_da(1,iel,iorb,ic,istate)=b_da(1,iel,iorb,ic,istate)-dum*rvec_en_sav(1,ic)*r_en_savi+db_tmp1
         b_da(2,iel,iorb,ic,istate)=b_da(2,iel,iorb,ic,istate)-dum*rvec_en_sav(2,ic)*r_en_savi+db_tmp2
         b_da(3,iel,iorb,ic,istate)=b_da(3,iel,iorb,ic,istate)-dum*rvec_en_sav(3,ic)*r_en_savi+db_tmp3

         do jc=1,ncent
            do k=1,3
               b_da(k,iel,iorb,jc,istate)=b_da(k,iel,iorb,jc,istate)
     &              +term_radial*(da_orbn(k,jc,iorb)+orbn(iorb)*da_ratio_jn(k,jc))
            enddo
         enddo
      enddo
      
      end subroutine
