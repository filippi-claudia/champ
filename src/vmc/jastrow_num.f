      module jastrow_num_mod
      contains
      subroutine jastrow_num(x,v,d2,value)
c Written by Cyrus Umrigar

c **Warning** This routine needs to be upgraded to calculate distances
c correctly for periodic systems if we add in capability to use
c numerical Laplacian for periodic systems.

      use bparm,   only: nocuspb,nspin2b
      use distance_mod, only: r_ee,r_en
      use jaspar6, only: c1_jas6
      use jastrow, only: ijas,is,nspin2,scalek,sspin,sspinn
      use multiple_geo, only: iwf
      use precision_kinds, only: dp
      use psi_mod, only: psi,psia,psib
      use system,  only: cent,iwctype,ncent,ncent_tot,nelec,nup
      use vmc_mod, only: nmat_dim2

      implicit none

      integer :: i, ic, ij, im1, ipar
      integer :: isb, it, j, jk
      integer :: k, m, ncentt
      real(dp) :: c3, d2, ftmp, psiac
      real(dp) :: psiami, psiami2, psiapi, psiapi2
      real(dp) :: psibc, psibmi, psibmi2, psibpi
      real(dp) :: psibpi2, psic, psimi, psimi2
      real(dp) :: psimj, psimj2, psipi, psipi2
      real(dp) :: psipj, psipj2, term, usum
      real(dp) :: value, wtj
      real(dp), dimension(3, *) :: x
      real(dp), dimension(3, *) :: v
      real(dp), dimension(3, nelec, ncent_tot) :: rp
      real(dp), dimension(3, nelec, ncent_tot) :: rm
      real(dp), dimension(3, nelec, ncent_tot) :: rp2
      real(dp), dimension(3, nelec, ncent_tot) :: rm2
      real(dp), dimension(3, nmat_dim2) :: rrp
      real(dp), dimension(3, nmat_dim2) :: rrm
      real(dp), dimension(3, nmat_dim2) :: rrp2
      real(dp), dimension(3, nmat_dim2) :: rrm2
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: half = .5d0
      real(dp), parameter :: eps = .5d-4
      real(dp), parameter :: eps2 = 2.d0*eps
      real(dp), parameter :: eps4 = 4.d0*eps
      real(dp), parameter :: epssq = eps**2
      real(dp), parameter :: eps2sq = eps2**2
      real(dp), parameter :: d1b12 = 8.333333333333333d-2
      real(dp), parameter :: d2b3 = 0.666666666666667d0
      real(dp), parameter :: d4b3 = 1.333333333333333d0

c subroutine to calculate jastrow factor,its derivatives
c and the potential

      do i=1,nelec
      v(1,i)=zero
      v(2,i)=zero
      v(3,i)=zero
      enddo
      d2=zero
      usum=zero

      ncentt=ncent

c Calculate e-N and e-e inter-particle distances
      do ic=1,ncent
        ij=0
        do i=1,nelec
          r_en(i,ic)=zero
          do m=1,3
            r_en(i,ic)=r_en(i,ic)+(x(m,i)-cent(m,ic))**2
          enddo
          do m=1,3
            rp(m,i,ic)=r_en(i,ic)+(x(m,i)-cent(m,ic))*eps2+epssq
            rm(m,i,ic)=r_en(i,ic)-(x(m,i)-cent(m,ic))*eps2+epssq
            rp2(m,i,ic)=r_en(i,ic)+(x(m,i)-cent(m,ic))*eps4+eps2sq
            rm2(m,i,ic)=r_en(i,ic)-(x(m,i)-cent(m,ic))*eps4+eps2sq
            rp(m,i,ic)=dsqrt(rp(m,i,ic))
            rm(m,i,ic)=dsqrt(rm(m,i,ic))
            rp2(m,i,ic)=dsqrt(rp2(m,i,ic))
            rm2(m,i,ic)=dsqrt(rm2(m,i,ic))
          enddo
          r_en(i,ic)=dsqrt(r_en(i,ic))

          do j=1,i-1
            ij=ij+1
            r_ee(ij)=(x(1,i)-x(1,j))**2 + (x(2,i)-x(2,j))**2
     &      + (x(3,i)-x(3,j))**2
            do m=1,3
              rrp(m,ij)=r_ee(ij)+(x(m,i)-x(m,j))*eps2+epssq
              rrm(m,ij)=r_ee(ij)-(x(m,i)-x(m,j))*eps2+epssq
              rrp2(m,ij)=r_ee(ij)+(x(m,i)-x(m,j))*eps4+eps2sq
              rrm2(m,ij)=r_ee(ij)-(x(m,i)-x(m,j))*eps4+eps2sq
              rrp(m,ij)=dsqrt(rrp(m,ij))
              rrm(m,ij)=dsqrt(rrm(m,ij))
              rrp2(m,ij)=dsqrt(rrp2(m,ij))
              rrm2(m,ij)=dsqrt(rrm2(m,ij))
            enddo
            r_ee(ij)=dsqrt(r_ee(ij))
          enddo
        enddo
      enddo

c if nelec is < 2 then no pairs so exit
      if (nelec.lt.2) goto 60

      ij=0
      do i=2,nelec
c3      jk=0
        im1=i-1
        do j=1,im1
          ij=ij+1
          if(i.le.nup .or. j.gt.nup) then
c           parallel spins
            if(nspin2.ge.2) then
              is=2
              sspinn=one
              if(nspin2.eq.3 .and. j.gt.nup) is=3
             else
              is=1
              sspinn=half
            endif
           else
c           anti-parallel spins
            sspinn=one
            is=1
          endif

          if(ijas.ge.3.and.ijas.le.6) then
            sspinn=one
            ipar=0
            if(i.le.nup .or. j.gt.nup) then
              isb=is
              if(nspin2b.eq.2) then
                isb=2
               elseif(nocuspb.eq.0) then
                sspinn=half
              endif
              ipar=1
             else
              isb=1
            endif

            psibc=psib(r_ee(ij),isb,ipar)
            usum=usum+psibc

            do m=1,3
              psibpi=psib(rrp(m,ij),isb,ipar)
              psibmi=psib(rrm(m,ij),isb,ipar)
              psibpi2=psib(rrp2(m,ij),isb,ipar)
              psibmi2=psib(rrm2(m,ij),isb,ipar)

              ftmp=(-d1b12*(psibpi2-psibmi2)+d2b3*(psibpi-psibmi))/eps
              v(m,i)=v(m,i)+ftmp
              v(m,j)=v(m,j)-ftmp

              d2=d2+two*(-d1b12*((psibmi2-psibc)+(psibpi2-psibc))
     &                   +d4b3 *((psibmi -psibc)+(psibpi -psibc)))
            enddo
          endif

          do ic=1,ncentt
          it=iwctype(ic)

          wtj=one
          if(ijas.eq.2) wtj=one/ncentt

          sspin=sspinn*wtj
          psic=psi(r_ee(ij),r_en(i,ic),r_en(j,ic),it)
          usum=usum+psic

          do m=1,3

            psipi=psi(rrp(m,ij),rp(m,i,ic),r_en(j,ic),it)
            psimi=psi(rrm(m,ij),rm(m,i,ic),r_en(j,ic),it)
            psipj=psi(rrm(m,ij),r_en(i,ic),rp(m,j,ic),it)
            psimj=psi(rrp(m,ij),r_en(i,ic),rm(m,j,ic),it)
            psipi2=psi(rrp2(m,ij),rp2(m,i,ic),r_en(j,ic),it)
            psimi2=psi(rrm2(m,ij),rm2(m,i,ic),r_en(j,ic),it)
            psipj2=psi(rrm2(m,ij),r_en(i,ic),rp2(m,j,ic),it)
            psimj2=psi(rrp2(m,ij),r_en(i,ic),rm2(m,j,ic),it)

            v(m,i)=v(m,i)+(-d1b12*(psipi2-psimi2)
     &      + d2b3*(psipi-psimi))/eps
            v(m,j)=v(m,j)+(-d1b12*(psipj2-psimj2)
     &      + d2b3*(psipj-psimj))/eps

            d2=d2-d1b12*((psimi2-psic)+(psipi2-psic))
     &           +d4b3 *((psimi -psic)+(psipi -psic))
     &           -d1b12*((psimj2-psic)+(psipj2-psic))
     &           +d4b3 *((psimj -psic)+(psipj -psic))
          enddo

          enddo
        enddo
      enddo

      if(ijas.ge.3.and.ijas.le.6) then
        do i=1,nelec
          is=1

          do ic=1,ncentt
            it=iwctype(ic)

            psiac=psia(r_en(i,ic),it)
            usum=usum+psiac

            do m=1,3
              psiapi=psia(rp(m,i,ic),it)
              psiami=psia(rm(m,i,ic),it)
              psiapi2=psia(rp2(m,i,ic),it)
              psiami2=psia(rm2(m,i,ic),it)

              v(m,i)=v(m,i)+(-d1b12*(psiapi2-psiami2)
     &        + d2b3*(psiapi-psiami))/eps

              d2=d2-d1b12*((psiami2-psiac)+(psiapi2-psiac))
     &             +d4b3 *((psiami -psiac)+(psiapi -psiac))
            enddo
          enddo
        enddo
      endif

      d2=d2/eps**2

   60 continue

      if(ijas.eq.6) then
        term=1/(c1_jas6*scalek(iwf))
        usum=term*usum
        d2=term*d2
        do i=1,nelec
          do k=1,3
            v(k,i)=term*v(k,i)
          enddo
        enddo
      endif

      value=usum

      return
      end
      end module
