module deriv_nonloc
contains
      subroutine deriv_nonlocj_quad1(nxquad,xquad,ielquad,x,r_en,rvec_en_quad,r_en_quad, &
                                    psij_ratio,dpsij_ratio,vjn,da_psij_ratio,iwfjas)

! Written by Claudia Filippi, modified by Cyrus Umrigar
      use bparm,   only: nocuspb,nspin2b
      use contrl_per, only: iperiodic
      use da_jastrow, only: da_j
      use deriv_nonlpsi, only: deriv_psianl,deriv_psibnl,deriv_psinl
      use derivjas, only: go
      use jaspointer, only: npoint,npointa
      use jastrow, only: is,nspin2,sspinn
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
!      write(ounit,*) 'nparmj,nparma(1),', nparmj, nparma(1)
      do it=2,nctype
        ipara=ipara+nparma(it)
!        write(ounit,*) 'it,nparma(it)', it, nparma(it)
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
          da_psij_ratio(k,ic,iq)=-dumk-da_j(k,iel,iel,ic)
        enddo
       enddo

      endif

      enddo

      return
      end

      subroutine deriv_nonlocj_quad4(nxquad,xquad,ielquad,x,rvec_en,r_en,rvec_en_quad,r_en_quad, &
          psij_ratio,dpsij_ratio,vjn,da_psij_ratio,iwfjas)

! Written by Claudia Filippi, modified by Cyrus Umrigar
      use bparm,   only: nocuspb,nspin2b
      use contrl_file, only: ounit
      use contrl_per, only: iperiodic
      use da_jastrow, only: da_j
      use deriv_nonlpsi, only: deriv_psianl,deriv_psibnl,deriv_psinl
      use derivjas, only: go
      use find_pimage, only: find_image_pbc
      use jaspointer, only: npoint,npointa
      use jastrow, only: is,nspin2,sspinn,nordc
      use jastrow_update, only: fso, fjo
      use m_force_analytic, only: iforce_analy
      use nonlpsi, only: dpsianl,dpsibnl,dpsinl
      use optwf_control, only: ioptjas
      use optwf_nparmj, only: nparma,nparmb,nparmc
      use optwf_parms, only: nparmj
      use precision_kinds, only: dp
      use qua, only: nquad
      use scale_dist_mod, only: scale_dist,scale_dist1
      use system,  only: iwctype,ncent,ncent_tot,nctype,nelec,nup
#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND)
      use qmckl
      use qmckl_data
      use error, only: fatal_error
      use deriv_jastrow_qmckl_mod, only: deriv_jastrowe_qmckl
#endif
      implicit none

      integer :: i, ic, iel, ipar, ipara
      integer :: iparm, iparm0, isb, it, nxquad
      integer :: j, jj, iq, jparm, k, iwfjas
      integer, dimension(nquad*nelec*2) :: ielquad


      real(dp) :: dd1i, dd1ij, dd1j, dd1u, dum, dumk
      real(dp) :: fi, fj, fsumn, fu, rij, u
      real(dp), dimension(3,*) :: x
      real(dp), dimension(3,*) :: xquad
      real(dp), dimension(3,nelec,ncent_tot) :: rvec_en
      real(dp), dimension(nelec,ncent_tot) :: r_en
      real(dp), dimension(3,nquad*nelec*2,*) :: rvec_en_quad
      real(dp), dimension(nquad*nelec*2,ncent_tot) :: r_en_quad
      real(dp), dimension(nelec,ncent_tot) :: rr_en
      real(dp), dimension(ncent_tot) :: rr_en_quad
      real(dp), dimension(*) :: psij_ratio
      real(dp), dimension(nparmj,*) :: dpsij_ratio
      real(dp), dimension(3,ncent_tot,*) :: da_psij_ratio
      real(dp), dimension(nelec,nelec) :: fsn
      real(dp), dimension(3) :: dx
      real(dp), dimension(nelec,ncent_tot) :: dd1
      real(dp), dimension(ncent_tot) :: dd1_quad
      real(dp), dimension(3,*) :: vjn
      real(dp), dimension(3, nelec) :: fjn
      real(dp), dimension(3, ncent) :: da_single_en, da_single_een
      real(dp) :: d2n
      real(dp), parameter :: half = .5d0
      integer(qmckl_exit_code) :: rc
      !for testing
      !real(dp), dimension(nparmj,nxquad) :: dpsij_ratio_new
      do ic=1,ncent
        if(iforce_analy.eq.0) then
          do i=1,nelec
            call scale_dist(r_en(i,ic),rr_en(i,ic))
          enddo
         else
          do i=1,nelec
            call scale_dist1(r_en(i,ic),rr_en(i,ic),dd1(i,ic))
          enddo
        endif
      enddo
      
      do iq=1,nxquad

      iel=ielquad(iq)

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND)
      
      if (iforce_analy .gt. 0) then
        call deriv_jastrowe_qmckl(iel,xquad(:,iq),fjn(:,:),d2n,fsumn,dpsij_ratio(:,iq),1)   

        rc = qmckl_get_forces_jastrow_single_en(qmckl_ctx(qmckl_no_ctx), da_single_en, 3_8*ncent)
        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl Jastrow single en force.')
        rc = qmckl_get_forces_jastrow_single_een(qmckl_ctx(qmckl_no_ctx), da_single_een, 3_8*ncent)
        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl Jastrow single een force.')
        do ic =1, ncent
          do k = 1, 3
            da_psij_ratio(k,ic,iq)=da_single_en(k,ic)+da_single_een(k,ic)
            ! write(ounit, *), 'da_psij', da_psij_ratio(k,ic,iq), da_single_en(k,ic), da_single_een(k,ic)
          enddo
        enddo
        do k=1,3
          vjn(k,iq)=fjn(k,iel)+fjo(k,iel,1)
          !write(ounit, *), 'vjn', vjn(k,iq)
        enddo

      else
        call deriv_jastrowe_qmckl(iel,xquad(:,iq),fjn(:,:),d2n,fsumn,dpsij_ratio(:,iq),0)
      endif
      
      ! print*, "begin deriv_jastrowe_qmckl"
      ! print*, "iel", iel
      ! print*, "xquad", xquad(:,iq)
      ! print*, "fjn", fjn(:,:)
      ! print*, "d2n", d2n
      ! print*, "fsumn", fsumn
      ! print*, "dpsij_ratio", dpsij_ratio(:,iq)
      ! print*, "end deriv_jastrowe_qmckl"

      psij_ratio(iq) = fsumn
#else
      if(iforce_analy.eq.0) then
        do ic=1,ncent
          call scale_dist(r_en_quad(iq,ic),rr_en_quad(ic))
        enddo
       else
        da_psij_ratio(:,:,iq)=0.d0
        do ic=1,ncent
          call scale_dist1(r_en_quad(iq,ic),rr_en_quad(ic),dd1_quad(ic))
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
!      write(ounit,*) 'nparmj,nparma(1)', nparmj, nparma(1)
      do it=2,nctype
        ipara=ipara+nparma(it)
!        write(ounit,*) 'it,nparma(it)', it, nparma(it)
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
          call scale_dist(rij,u)
         else
          call scale_dist1(rij,u,dd1u)
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

        if(nordc.gt.1) then
          do ic=1,ncent
            it=iwctype(ic)
            if(nparmc(it).gt.0) then
              iparm0=npoint(it)
              fsn(i,j)=fsn(i,j) + deriv_psinl(u,rr_en_quad(ic),rr_en(jj,ic),dpsij_ratio(iparm0+1,iq),it,iwfjas)
            endif
          enddo

          do it=1,nctype
            iparm0=npoint(it)
            do jparm=1,nparmc(it)
              iparm=iparm0+jparm
              dpsij_ratio(iparm,iq)=dpsij_ratio(iparm,iq)-go(i,j,iparm,iwfjas)
            enddo
          enddo

          if(iforce_analy.gt.0) then
            do ic=1,ncent
              it=iwctype(ic)
              dd1ij=dd1u
              dd1i=dd1_quad(ic)
              dd1j=dd1(jj,ic)
              dum=dpsinl(u,rr_en_quad(ic),rr_en(jj,ic),fu,fi,fj,dd1ij,dd1i,dd1j,it,iwfjas)
              fu=fu*dd1ij/rij
              fi=fi*dd1i/r_en_quad(iq,ic)
              fj=fj*dd1j/r_en(jj,ic)

              do k=1,3
                dumk=fi*rvec_en_quad(k,iq,ic)
                vjn(k,iq)=vjn(k,iq)+dumk-fu*dx(k)
                dumk=-dumk-fj*rvec_en(k,jj,ic)
                da_psij_ratio(k,ic,iq)=da_psij_ratio(k,ic,iq)+dumk-da_j(k,i,j,ic)
              enddo
            enddo
          endif
        endif

        fsumn=fsumn+fsn(i,j)-fso(i,j,iwfjas)

      45 continue
      enddo

! e-n terms
      47 fsn(iel,iel)=0

      do ic=1,ncent
        it=iwctype(ic)
        iparm0=npointa(it)
        fsn(iel,iel)=fsn(iel,iel)+ &
        deriv_psianl(rr_en_quad(ic),dpsij_ratio(iparm0+1,iq),it,iwfjas)
!        write(ounit,*) 'ic,it,iwctype(ic),iparm0,iparm0+1', ic,it,iwctype(ic),iparm0,iparm0+1
      enddo
      do it=1,nctype
        iparm0=npointa(it)
        do jparm=1,nparma(it)
          iparm=iparm0+jparm
          dpsij_ratio(iparm,iq)=dpsij_ratio(iparm,iq)-go(iel,iel,iparm,iwfjas)
!          write(ounit,*) 'it,npointa(it),jparm,iparm0,iparm', it,npointa(it),jparm,iparm0,iparm
        enddo
      enddo

      fsumn=fsumn+fsn(iel,iel)-fso(iel,iel,iwfjas)
      psij_ratio(iq)=fsumn

      if(iforce_analy.gt.0) then
        do ic=1,ncent
          it=iwctype(ic)
          dum=dpsianl(rr_en_quad(ic),it,iwfjas)*dd1_quad(ic)/r_en_quad(iq,ic)
          do k=1,3
            dumk=dum*rvec_en_quad(k,iq,ic)
            vjn(k,iq)=vjn(k,iq)+dumk
            da_psij_ratio(k,ic,iq)=da_psij_ratio(k,ic,iq)-dumk-da_j(k,iel,iel,ic)
          enddo
        enddo
      endif
#endif
      ! print*, "begin deriv_jastrowe old"
      ! print*, "iel", iel
      ! print*, "xquad", xquad(:,iq)
      ! print*, "fjn", fjn(:,:)
      ! print*, "d2n", d2n
      ! print*, "fsumn", fsumn
      ! print*, "dpsij_ratio", dpsij_ratio(:,iq)
      ! print*, "end deriv_jastrowe old"

      enddo

      return
      end
end module
