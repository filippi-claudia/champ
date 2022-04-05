      module basis_fns_mod
      contains
      subroutine basis_fns(ie1,ie2,rvec_en,r_en,ider)
c calculate the values of the basis functions and their derivatives
c ider = 0 -> value
c ider = 1 -> value, gradient
c ider = 2 -> value, gradient, laplacian
c ider = 3 -> value, gradient, laplacian, forces

      use numbas_mod, only: MRWF
      use atom, only: iwctype, ncent, ncent_tot
      use ghostatom, only: nghostcent
      use const, only: nelec
      use numbas, only: iwrwf, nrbas
      use numbas1, only: iwlbas, nbastyp
      use phifun, only: phin, dphin, d2phin, d2phin_all, d3phin, n0_nbasis
      use wfsec, only: iwf
      use force_analy, only: iforce_analy
      use splfit_mod, only: splfit
      use slm_mod, only: slm

      use precision_kinds, only: dp
      implicit none

      integer :: it, ic, ider, irb
      integer :: iwlbas0, j, k
      integer :: ie1, ie2, k0, l, l0, ll
      real(dp) :: y, ddy_lap
      real(dp) :: r, r2, ri, ri2
      real(dp), dimension(3, nelec, ncent_tot) :: rvec_en
      real(dp), dimension(nelec, ncent_tot) :: r_en
      real(dp), dimension(3) :: dy
      real(dp), dimension(3, 3) :: ddy
      real(dp), dimension(3) :: dlapy
      real(dp), dimension(4, MRWF) :: wfv
      real(dp), dimension(3) :: xc
      real(dp), parameter :: one = 1.d0

      do j=ie1,ie2
        n0_nbasis(j)=0
      enddo

      l=0
c loop through centers
      do ic=1,ncent+nghostcent

        it=iwctype(ic)

        ll=0
        k0=0
        l0=l

c numerical atomic orbitals
        do k=ie1,ie2

c get distance to center
          xc(1)=rvec_en(1,k,ic)
          xc(2)=rvec_en(2,k,ic)
          xc(3)=rvec_en(3,k,ic)

          r=r_en(k,ic)
          r2=r*r
          ri=one/r
          ri2=ri*ri

          do irb=1,nrbas(it)
            call splfit(r,irb,it,iwf,wfv(1,irb),ider)
          enddo

c compute sml and combine to generate molecular orbitals
          l=l0
          ll=0
          iwlbas0=0
          do j=1,nbastyp(it)
            l=l+1
            ll=ll+1
            irb=iwrwf(ll,it)
            if(iwlbas(ll,it).ne.iwlbas0.or.k.ne.k0) then
              k0=k
              iwlbas0=iwlbas(ll,it)
              call slm(iwlbas0,xc,r2,y,dy,ddy,ddy_lap,dlapy,ider)
            endif

            call phi_combine(iwlbas0,xc,ri,ri2,wfv(1,irb),y,dy,ddy,ddy_lap,dlapy,
     &           phin(l,k),dphin(l,k,:),d2phin(l,k),d2phin_all(1,1,l,k),d3phin(1,l,k),ider)
                       
            
            call n0_inc(l,k,ic)
          enddo

        enddo

c loop over all atoms
      enddo

      return
      end
c-------------------------------------------------------------------
      subroutine phi_combine(l,xc,ri,ri2,wfv,y,dy,ddy,ddy_lap,dlapy,phi,dphi,d2phi,d2phi_all,d3phi,ider)
      use precision_kinds, only: dp
      implicit none

      integer :: iforce_analy, ider, ii, jj, l
      real(dp) :: d2phi, ddy_lap, dum, dum1, phi
      real(dp) :: prod, ri, ri2, ri3
      real(dp) :: y
      real(dp), dimension(3) :: xc
      real(dp), dimension(3) :: xcri
      real(dp), dimension(4) :: wfv
      real(dp), dimension(4) :: wfvn
      real(dp), dimension(3) :: dy
      real(dp), dimension(3, 3) :: ddy
      real(dp), dimension(3) :: dphi
      real(dp), dimension(3, 3) :: d2phi_all
      real(dp), dimension(3) :: d3phi
      real(dp), dimension(3) :: dlapy
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: three = 3.d0
      real(dp), parameter :: four = 4.d0
      real(dp), parameter :: five = 5.d0
      real(dp), parameter :: six = 6.d0

      if(ider.eq.0) then

        if(l.eq.1) then
          wfvn(1)=wfv(1)
        elseif(l.lt.5) then
          wfvn(1)=wfv(1)*ri
        elseif(l.lt.10) then
          wfvn(1)=wfv(1)*ri2
        elseif(l.lt.20) then
          ri3=ri*ri2
          wfvn(1)=wfv(1)*ri3
        else
         stop 'to fix for >f functions'
        endif

        phi=y*wfvn(1)

      elseif(ider.eq.1) then

        xcri(1)=xc(1)*ri
        xcri(2)=xc(2)*ri
        xcri(3)=xc(3)*ri
        if(l.eq.1) then
          wfvn(1)=wfv(1)
          wfvn(2)=wfv(2)
        elseif(l.lt.5) then
          wfvn(2)=-wfv(1)*ri2+wfv(2)*ri
          wfvn(1)=wfv(1)*ri
        elseif(l.lt.10) then
          wfvn(2)=(-two*wfv(1)*ri+wfv(2))*ri2
          wfvn(1)=wfv(1)*ri2
        elseif(l.lt.20) then
          ri3=ri*ri2
          wfvn(2)=(-three*wfv(1)*ri+wfv(2))*ri3
          wfvn(1)=wfv(1)*ri3
        else
         stop 'to fix for >f functions'
        endif

        phi=y*wfvn(1)
        do jj=1,3
          dphi(jj)=y*xcri(jj)*wfvn(2)+dy(jj)*wfvn(1)
        enddo

      elseif(ider.ge.2) then

        xcri(1)=xc(1)*ri
        xcri(2)=xc(2)*ri
        xcri(3)=xc(3)*ri
        if(l.eq.1) then
          wfvn(1)=wfv(1)
          wfvn(2)=wfv(2)
          wfvn(3)=wfv(3)
        elseif(l.lt.5) then
          wfvn(3)=ri*(wfv(3)+two*ri*(wfv(1)*ri-wfv(2)))
          wfvn(2)=-wfv(1)*ri2+wfv(2)*ri
          wfvn(1)=wfv(1)*ri
        elseif(l.lt.10) then
          wfvn(3)=ri2*(wfv(3)+two*ri*(three*wfv(1)*ri-two*wfv(2)))
          wfvn(2)=(-two*wfv(1)*ri+wfv(2))*ri2
          wfvn(1)=wfv(1)*ri2
        elseif(l.lt.20) then
          ri3=ri*ri2
          wfvn(3)=ri3*(wfv(3)+six*ri*(two*wfv(1)*ri-wfv(2)))
          wfvn(2)=(-three*wfv(1)*ri+wfv(2))*ri3
          wfvn(1)=wfv(1)*ri3
        else
         stop 'to fix for >f functions'
        endif

        phi=y*wfvn(1)
        d2phi=y*wfvn(3)+y*two*ri*wfvn(2)+ddy_lap*wfvn(1)
        dum=0.d0
        do jj=1,3
          dphi(jj)=y*xcri(jj)*wfvn(2)+dy(jj)*wfvn(1)
          dum=dum+dy(jj)*xcri(jj)
        enddo
        d2phi=d2phi+two*dum*wfvn(2)

        if(ider.eq.3) then

          if(l.eq.1) then
            wfvn(4)=wfv(4)
          elseif(l.lt.5) then
            wfvn(4)=-ri*wfvn(3)+ri*(wfv(4)-two*ri*(two*wfv(1)*ri2-two*wfv(2)*ri+wfv(3)))
          elseif(l.lt.10) then
            wfvn(4)=-two*ri*wfvn(3)+ri2*(wfv(4)-two*ri*(six*wfv(1)*ri2-five*wfv(2)*ri+two*wfv(3)))
          elseif(l.lt.20) then
            wfvn(4)=-three*ri*wfvn(3)+ri3*(wfv(4)-six*ri*(four*wfv(1)*ri2-three*wfv(2)*ri+wfv(3)))
          else
           stop 'to fix for >f functions'
          endif

          do jj=1,3
            dum1=0
            do ii=1,3
              dum1=dum1+ddy(jj,ii)*xcri(ii)
            enddo
            d3phi(jj)=wfvn(4)*y*xcri(jj)
     &               +wfvn(3)*(dy(jj)+two*xcri(jj)*(y*ri+dum))
     &               +wfvn(2)*(xcri(jj)*(ddy_lap-two*ri*(dum+y*ri))+two*(dum1+two*dy(jj)*ri))
     &               +wfvn(1)*dlapy(jj)
          enddo

          do jj=1,3
            do ii=jj,3
              prod=xcri(jj)*xcri(ii)
              d2phi_all(ii,jj)=ddy(ii,jj)*wfvn(1)+wfvn(2)*(dy(ii)*xcri(jj)+dy(jj)*xcri(ii)-y*ri*prod)+wfvn(3)*y*prod
              d2phi_all(jj,ii)=d2phi_all(ii,jj)
            enddo
            d2phi_all(jj,jj)=d2phi_all(jj,jj)+y*ri*wfvn(2)
          enddo
        endif

      endif

      return
      end
c-------------------------------------------------------------------
      subroutine n0_inc(l,k,ic)

      use phifun, only: phin, dphin, n0_ibasis, n0_ic, n0_nbasis
      implicit none

      integer :: ic, k, l

      if(abs(phin(l,k))+abs(dphin(l,k,1))+abs(dphin(l,k,2))+abs(dphin(l,k,3)).gt.1.d-20)then
       n0_nbasis(k)=n0_nbasis(k)+1
       n0_ibasis(n0_nbasis(k),k)=l
       n0_ic(n0_nbasis(k),k)=ic
      endif

      return
      end
      end module
