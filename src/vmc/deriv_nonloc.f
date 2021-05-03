      subroutine deriv_nonlocj(istate,iel,x,rshift,rvec_en,r_en,rr_en,
     &     rr_en2,dd1,value,gn,vjn,da_ratio_jn)
c     Written by Claudia Filippi, modified by Cyrus Umrigar
      use vmc_mod, only: MELEC, MCENT
      use atom, only: iwctype, nctype, ncent
      use jaspar, only: nspin2, sspinn, is
      use const, only: nelec
      use da_jastrow4val, only: da_j
      use derivjas, only: go, gvalue
      use elec, only: nup
      use jaso, only: fso
      use jaspointer, only: npoint, npointa
      use optwf_nparmj, only: nparma, nparmb, nparmc
      use optwf_parms, only: nparmj
      use bparm, only: nocuspb, nspin2b
      use contr2, only: ijas
      use contr2, only: isc
      use contrl_per, only: iperiodic
      use force_analy, only: iforce_analy

      implicit real*8(a-h,o-z)

      parameter (half=.5d0)
      dimension x(3,*),rshift(3,MELEC,MCENT),rvec_en(3,MELEC,*)
      dimension r_en(MELEC,MCENT),rr_en(MELEC,MCENT),rr_en2(MELEC,MCENT)
     &     ,gn(*),fsn(MELEC,MELEC),dx(3),dd1(MELEC,MCENT),vjn(3),da_ratio_jn(3,MCENT)

      fsumn=0.0d0

      do k=1,3
         vjn(k)=0.0d0
      enddo

      do iparm=1,nparmj
         gn(iparm)=0.0d0
      enddo

      if (nelec.lt.2) goto 47

      ipara=nparma(1)
      if(ijas.ge.4.and.ijas.le.6) then
         do it=2,nctype
            ipara=ipara+nparma(it)
         enddo
      endif

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
         if(i.le.nup .or. j.gt.nup) then
            if(nspin2.ge.2) then
               is=2
               isb=is
               if(nspin2.eq.3 .and. j.gt.nup) then
                  is=3
                  isb=is
               endif
            else
               is=1
               isb=is
               if(nspin2b.eq.2) then
                  isb=2
               elseif(nocuspb.eq.0) then
                  sspinn=half
               endif
            endif
            ipar=1
         else
            is=1
            isb=1
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
            dum=dpsibnl(u,isb,ipar)*dd1u/rij
            do k=1,3
               dumk=-dum*dx(k)
               vjn(k)=vjn(k)+dumk
            enddo
         endif

         iparm0=ipara
         if(isb.eq.2) iparm0=iparm0+nparmb(1)
         fsn(i,j)=deriv_psibnl(u,gn(iparm0+1),isb,ipar)

         do jparm=1,nparmb(isb)
            iparm=iparm0+jparm
            gn(iparm)=gn(iparm)-go(i,j,iparm,istate)
         enddo

c     e-e-n terms
c     The scaling is switched in deriv_psinl, so do not do it here.
         if(isc.ge.12) call scale_dist(rij,u,3)

         do ic=1,ncent
            it=iwctype(ic)
            iparm0=npoint(it)
            fsn(i,j)=fsn(i,j) +
     &           deriv_psinl(u,rshift(1,i,ic),rshift(1,j,ic),
     &           rr_en2(i,ic),rr_en2(j,ic),gn(iparm0+1),it)
         enddo

         do it=1,nctype
            iparm0=npoint(it)
            do jparm=1,nparmc(it)
               iparm=iparm0+jparm
               gn(iparm)=gn(iparm)-go(i,j,iparm,istate)
            enddo
         enddo

         fsumn=fsumn+fsn(i,j)-fso(i,j,istate)
      enddo

c     e-n terms
 47   fsn(iel,iel)=0.0d0

      if(ijas.ge.4.and.ijas.le.6) then
         do ic=1,ncent
            it=iwctype(ic)
            iparm0=npointa(it)
            fsn(iel,iel)=fsn(iel,iel)+
     &           deriv_psianl(rr_en(iel,ic),gn(iparm0+1),it)
         enddo
         do it=1,nctype
            iparm0=npointa(it)
            do jparm=1,nparma(it)
               iparm=iparm0+jparm
               gn(iparm)=gn(iparm)-go(iel,iel,iparm,istate)
            enddo
         enddo
      endif

      fsumn=fsumn+fsn(iel,iel)-fso(iel,iel,istate)
      value=fsumn

      if(iforce_analy.eq.0) return

      do ic=1,ncent
         it=iwctype(ic)
         dum=dpsianl(rr_en(iel,ic),it)*dd1(iel,ic)/r_en(iel,ic)
         do k=1,3
            dumk=dum*rvec_en(k,iel,ic)
            vjn(k)=vjn(k)+dumk
            da_ratio_jn(k,ic)=-dumk-da_j(k,iel,ic,istate)
         enddo
      enddo

      end subroutine
