module deriv_nonloc
contains
      subroutine deriv_nonlocj_quad1(nxquad,xquad,ielquad,x,r_en,rvec_en_quad,r_en_quad, &
                                    psij_ratio,dpsij_ratio,vjn,da_psij_ratio,iwfjas)

! Written by Claudia Filippi, modified by Cyrus Umrigar
      use bparm,   only: nocuspb,nspin2b
      use contrl_per, only: iperiodic
      use da_jastrow4val, only: da_j
      use deriv_nonlpsi, only: deriv_psianl,deriv_psibnl,deriv_psinl
      use derivjas, only: go
      use jaspointer, only: npoint,npointa
      use jastrow, only: ijas,is,isc,nspin2,sspinn
      use jastrow_update, only: fso
      use m_force_analytic, only: iforce_analy
      use nonlpsi, only: dpsianl,dpsibnl
      use optwf_control, only: ioptjas
      use optwf_nparmj, only: nparma,nparmb,nparmc
      use optwf_parms, only: nparmj
      use precision_kinds, only: dp
      use find_pimage, only: find_image_pbc
      use qua, only: nquad
      use system,  only: iwctype,ncent,ncent_tot,nctype,nelec,nup
      use contrl_file, only: ounit

      implicit none

      integer :: i, ic, iel, ipar, ipara
      integer :: iparm, iparm0, isb, it, nxquad
      integer :: j, jj, iq, jparm, k, iwfjas
      integer, dimension(nquad*nelec*2) :: ielquad


      real(dp) :: dum, dumk, fsumn, rij
      real(dp), dimension(3,*) :: x
      real(dp), dimension(3,*) :: xquad
      real(dp), dimension(nelec,ncent_tot) :: r_en
      real(dp), dimension(3,nquad*nelec*2,*) :: rvec_en_quad
      real(dp), dimension(nquad*nelec*2,ncent_tot) :: r_en_quad
      real(dp), dimension(*) :: psij_ratio
      real(dp), dimension(nparmj,*) :: dpsij_ratio
      real(dp), dimension(3,ncent_tot,*) :: da_psij_ratio
      real(dp), dimension(nelec,nelec) :: fsn
      real(dp), dimension(3) :: dx
      real(dp), dimension(3,*) :: vjn
      real(dp), parameter :: half = .5d0


      do iq=1,nxquad

      iel=ielquad(iq)

      fsumn=0

      do k=1,3
        vjn(k,iq)=0.d0
      enddo

! TMP
!     do 5 iparm=1,nparmj
!   5   dpsij_ratio(iparm)=gvalue(iparm)
      do iparm=1,nparmj
        dpsij_ratio(iparm,iq)=0
      enddo

      if (nelec.lt.2) goto 47

      ipara=nparma(1)
!      write(ounit,*) 'nparmj,nparma(1),ijas', nparmj, nparma(1), ijas
        do it=2,nctype
          ipara=ipara+nparma(it)
!          write(ounit,*) 'it,nparma(it)', it, nparma(it)
        enddo


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
            sspinn=half
            if(nspin2b.eq.2) then
              isb=2
             elseif(nocuspb.gt.0) then
              sspinn=1
            endif
          endif
          ipar=1
         else
          is=1
          isb=1
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
          call find_image_pbc(dx,rij)
        endif

! e-e terms
        if(iforce_analy.gt.0) then
           dum=dpsibnl(rij,isb,ipar,iwfjas)/rij
           do k=1,3
              dumk=-dum*dx(k)
              vjn(k,iq)=vjn(k,iq)+dumk
           enddo
        endif

        iparm0=ipara
        if(isb.eq.2) iparm0=iparm0+nparmb(1)
        fsn(i,j)=deriv_psibnl(rij,dpsij_ratio(iparm0+1,iq),isb,ipar,iwfjas)

        do jparm=1,nparmb(isb)
          iparm=iparm0+jparm
!          write(ounit,*) 'jparm,iparm0,iparm,nparmb(isb),isb', jparm,iparm0,iparm,nparmb(isb),isb
          dpsij_ratio(iparm,iq)=dpsij_ratio(iparm,iq)-go(i,j,iparm,iwfjas)
        enddo

! e-e-n terms
        do ic=1,ncent
          it=iwctype(ic)
          if(nparmc(it).gt.0) then
            iparm0=npoint(it)
            fsn(i,j)=fsn(i,j) + &
            deriv_psinl(rij,r_en_quad(iq,ic),r_en(jj,ic),dpsij_ratio(iparm0+1,iq),it,iwfjas)
          endif
        enddo

        do it=1,nctype
          iparm0=npoint(it)
          do jparm=1,nparmc(it)
            iparm=iparm0+jparm
            dpsij_ratio(iparm,iq)=dpsij_ratio(iparm,iq)-go(i,j,iparm,iwfjas)
          enddo
        enddo

        fsumn=fsumn+fsn(i,j)-fso(i,j,iwfjas)
   45 continue
      enddo

! e-n terms
   47 fsn(iel,iel)=0

      do ic=1,ncent
         it=iwctype(ic)
         iparm0=npointa(it)
         fsn(iel,iel)=fsn(iel,iel)+ &
              deriv_psianl(r_en_quad(iq,ic),dpsij_ratio(iparm0+1,iq),it,iwfjas)
!          write(ounit,*) 'ic,it,iwctype(ic),iparm0,iparm0+1', ic,it,iwctype(ic),iparm0,iparm0+1
      enddo
      do it=1,nctype
         iparm0=npointa(it)
         do jparm=1,nparma(it)
            iparm=iparm0+jparm
            dpsij_ratio(iparm,iq)=dpsij_ratio(iparm,iq)-go(iel,iel,iparm,iwfjas)
!     write(ounit,*) 'it,npointa(it),jparm,iparm0,iparm', it,npointa(it),jparm,iparm0,iparm
         enddo
      enddo

      fsumn=fsumn+fsn(iel,iel)-fso(iel,iel,iwfjas)
      psij_ratio(iq)=fsumn

      if(iforce_analy.gt.0) then

       do ic=1,ncent
        it=iwctype(ic)
        dum=dpsianl(r_en_quad(iq,ic),it,iwfjas)/r_en_quad(iq,ic)
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

      subroutine deriv_nonlocj_quad4(nxquad,xquad,ielquad,x,r_en,rvec_en_quad,r_en_quad, &
          psij_ratio,dpsij_ratio,vjn,da_psij_ratio,iwfjas)

! Written by Claudia Filippi, modified by Cyrus Umrigar
      use bparm,   only: nocuspb,nspin2b
      use contrl_per, only: iperiodic
      use da_jastrow4val, only: da_j
      use deriv_nonlpsi, only: deriv_psianl,deriv_psibnl,deriv_psinl
      use derivjas, only: go
      use jaspointer, only: npoint,npointa
      use jastrow, only: ijas,is,isc,nspin2,sspinn
      use jastrow_update, only: fso
      use m_force_analytic, only: iforce_analy
      use nonlpsi, only: dpsianl,dpsibnl
      use optwf_control, only: ioptjas
      use optwf_nparmj, only: nparma,nparmb,nparmc
      use optwf_parms, only: nparmj
      use precision_kinds, only: dp
      use find_pimage, only: find_image_pbc
      use qua, only: nquad
      use scale_dist_mod, only: scale_dist,scale_dist1
      use system,  only: iwctype,ncent,ncent_tot,nctype,nelec,nup
      use contrl_file, only: ounit

      implicit none

      integer :: i, ic, iel, ipar, ipara
      integer :: iparm, iparm0, isb, it, nxquad
      integer :: j, jj, iq, jparm, k, iwfjas
      integer, dimension(nquad*nelec*2) :: ielquad


      real(dp) :: dd1u, dum
      real(dp) :: dumk, fsumn, rij, u
      real(dp), dimension(3,*) :: x
      real(dp), dimension(3,*) :: xquad
      real(dp), dimension(nelec,ncent_tot) :: r_en
      real(dp), dimension(3,nquad*nelec*2,*) :: rvec_en_quad
      real(dp), dimension(nquad*nelec*2,ncent_tot) :: r_en_quad
      real(dp), dimension(nelec,ncent_tot) :: rr_en
      real(dp), dimension(nelec,ncent_tot) :: rr_en2
      real(dp), dimension(ncent_tot) :: rr_en_quad
      real(dp), dimension(ncent_tot) :: rr_en2_quad
      real(dp), dimension(*) :: psij_ratio
      real(dp), dimension(nparmj,*) :: dpsij_ratio
      real(dp), dimension(3,ncent_tot,*) :: da_psij_ratio
      real(dp), dimension(nelec,nelec) :: fsn
      real(dp), dimension(3) :: dx
      real(dp), dimension(nelec,ncent_tot) :: dd1
      real(dp), dimension(ncent_tot) :: dd1_quad
      real(dp), dimension(3,*) :: vjn
      real(dp), parameter :: half = .5d0

      do ic=1,ncent
!JF this is the culprit
        if(iforce_analy.eq.0) then
          do i=1,nelec
            call scale_dist(r_en(i,ic),rr_en(i,ic),1)
            call scale_dist(r_en(i,ic),rr_en2(i,ic),2)
          enddo
         else
          do i=1,nelec
            call scale_dist1(r_en(i,ic),rr_en(i,ic),dd1(i,ic),1)
!JF added to see what happens --> gives same as iforce_analy = 0
!           call scale_dist(r_en(i,ic),rr_en2(i,ic),2)
            if(ioptjas.gt.0) call scale_dist(r_en(i,ic),rr_en2(i,ic),2)
          enddo
        endif
      enddo

      do iq=1,nxquad

      iel=ielquad(iq)

      if(iforce_analy.eq.0) then
        do ic=1,ncent
          call scale_dist(r_en_quad(iq,ic),rr_en_quad(ic),1)
          call scale_dist(r_en_quad(iq,ic),rr_en2_quad(ic),2)
        enddo
       else
        do ic=1,ncent
          call scale_dist1(r_en_quad(iq,ic),rr_en_quad(ic),dd1_quad(ic),1)
!JF added to see what happens --> gives same as iforce_analy = 0
!         call scale_dist(r_en_quad(iq,ic),rr_en2_quad(ic),2)
          if(ioptjas.gt.0) call scale_dist(r_en_quad(iq,ic),rr_en2_quad(ic),2)
        enddo
      endif

      fsumn=0

      do k=1,3
        vjn(k,iq)=0.d0
      enddo

! TMP
!     do 5 iparm=1,nparmj
!   5   dpsij_ratio(iparm)=gvalue(iparm)
      do iparm=1,nparmj
        dpsij_ratio(iparm,iq)=0
      enddo

      if (nelec.lt.2) goto 47

      ipara=nparma(1)
!      write(ounit,*) 'nparmj,nparma(1),ijas', nparmj, nparma(1), ijas
      if(ijas.ge.4.and.ijas.le.6) then
        do it=2,nctype
          ipara=ipara+nparma(it)
!          write(ounit,*) 'it,nparma(it)', it, nparma(it)
        enddo
      endif

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
          dx(k)=x(k,jj)-xquad(k,iq)
        enddo
        if(iperiodic.eq.0) then
          rij=0
          do k=1,3
            rij=rij+dx(k)**2
          enddo
          rij=dsqrt(rij)
         else
          call find_image_pbc(dx,rij)
        endif

! e-e terms
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

        iparm0=ipara
        if(isb.eq.2) iparm0=iparm0+nparmb(1)
        fsn(i,j)=deriv_psibnl(u,dpsij_ratio(iparm0+1,iq),isb,ipar,iwfjas)

        do jparm=1,nparmb(isb)
          iparm=iparm0+jparm
!          write(ounit,*) 'jparm,iparm0,iparm,nparmb(isb),isb', jparm,iparm0,iparm,nparmb(isb),isb
          dpsij_ratio(iparm,iq)=dpsij_ratio(iparm,iq)-go(i,j,iparm,iwfjas)
        enddo

! e-e-n terms
! The scaling is switched in deriv_psinl, so do not do it here.
      if(isc.ge.12) call scale_dist(rij,u,3)

        do ic=1,ncent
          it=iwctype(ic)
          if(nparmc(it).gt.0) then
            iparm0=npoint(it)
            fsn(i,j)=fsn(i,j) + deriv_psinl(u,rr_en2_quad(ic),rr_en2(jj,ic),dpsij_ratio(iparm0+1,iq),it,iwfjas)
          endif
        enddo

        do it=1,nctype
          iparm0=npoint(it)
          do jparm=1,nparmc(it)
            iparm=iparm0+jparm
            dpsij_ratio(iparm,iq)=dpsij_ratio(iparm,iq)-go(i,j,iparm,iwfjas)
          enddo
        enddo

        fsumn=fsumn+fsn(i,j)-fso(i,j,iwfjas)
      45 continue
      enddo

! e-n terms
      47 fsn(iel,iel)=0

      if(ijas.ge.4.and.ijas.le.6) then
        do ic=1,ncent
          it=iwctype(ic)
          iparm0=npointa(it)
          fsn(iel,iel)=fsn(iel,iel)+ &
          deriv_psianl(rr_en_quad(ic),dpsij_ratio(iparm0+1,iq),it,iwfjas)
!          write(ounit,*) 'ic,it,iwctype(ic),iparm0,iparm0+1', ic,it,iwctype(ic),iparm0,iparm0+1
        enddo
        do it=1,nctype
          iparm0=npointa(it)
          do jparm=1,nparma(it)
            iparm=iparm0+jparm
            dpsij_ratio(iparm,iq)=dpsij_ratio(iparm,iq)-go(iel,iel,iparm,iwfjas)
!            write(ounit,*) 'it,npointa(it),jparm,iparm0,iparm', it,npointa(it),jparm,iparm0,iparm
          enddo
        enddo
      endif

      fsumn=fsumn+fsn(iel,iel)-fso(iel,iel,iwfjas)
      psij_ratio(iq)=fsumn

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
end module
