      module psi_mod
      use scale_dist_mod, only: scale_dist,switch_scale
      contains
      function psi(rij,ri,rj,it)

      use jastrow, only: nordc
      use jaspar6, only: cutjas
      use jastrow, only: c,ijas,nordj
      use multiple_geo, only: iwf
      use precision_kinds, only: dp
      implicit none

      integer :: it, jp, k, l, l_hi
      integer :: ll, m, n
      real(dp) :: psi, ri, rij, rj
      real(dp) :: rri, rrj, s, t
      real(dp) :: u
      real(dp), dimension(0:nordj) :: uu
      real(dp), dimension(0:nordj) :: ss
      real(dp), dimension(0:nordj) :: tt
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0


















      psi=0

      u=rij
      s=ri+rj
      t=dabs(ri-rj)

      call scale_dist(rij,u,3)
      call scale_dist(ri,rri,2)
      call scale_dist(rj,rrj,2)

      if(nordc.le.1) return

      if(ri.gt.cutjas .or. rj.gt.cutjas) return

      if(ijas.eq.4.or.ijas.eq.5) then
        call switch_scale(u)
        call switch_scale(rri)
        call switch_scale(rrj)
      endif

      uu(0)=one
      ss(0)=2
      tt(0)=one
      do jp=1,nordc
        uu(jp)=u**jp
        ss(jp)=rri**jp+rrj**jp
        tt(jp)=(rri*rrj)**jp
      enddo

      ll=0
      do n=2,nordc
        do k=n-1,0,-1
          if(k.eq.0) then
            l_hi=n-k-2
           else
            l_hi=n-k
          endif
          do l=l_hi,0,-1
            m=(n-k-l)/2
            if(2*m.eq.n-k-l) then
              ll=ll+1
              psi=psi+c(ll,it,iwf)*uu(k)*ss(l)*tt(m)
            endif
          enddo
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      function psia(ri,it)


      use jastrow, only: norda
      use jaspar6, only: cutjas
      use jastrow, only: a4,asymp_jasa,ijas
      use multiple_geo, only: iwf
      use precision_kinds, only: dp
      implicit none

      integer :: i, it
      real(dp) :: ri, rri
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp) :: psia










      psia=zero
      if(ijas.lt.4.or.ijas.gt.6) return

      if(ri.gt.cutjas) return

      call scale_dist(ri,rri,1)

      if(ijas.eq.4.or.ijas.eq.5) then
        psia=a4(1,it,iwf)*rri/(one+a4(2,it,iwf)*rri)-asymp_jasa(it)
        do i=2,norda
          psia=psia+a4(i+1,it,iwf)*rri**i
        enddo
       elseif(ijas.eq.6) then
        psia=a4(1,it,iwf)*rri/(one+a4(2,it,iwf)*(1-rri))
        do i=2,norda
          psia=psia+a4(i+1,it,iwf)*rri**i
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      function psib(rij,isb,ipar)


      use jastrow, only: nordb
      use jaspar6, only: cutjas
      use jastrow, only: asymp_jasb,b,ijas,sspinn
      use multiple_geo, only: iwf
      use precision_kinds, only: dp
      implicit none

      integer :: i, ipar, isb
      real(dp) :: rij, u
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp) :: psib










      psib=zero
      if(ijas.lt.4.or.ijas.gt.6) return

      if(rij.gt.cutjas) return

      u=rij
      call scale_dist(rij,u,1)

      if(ijas.eq.4) then
        psib=sspinn*b(1,isb,iwf)*u/(one+b(2,isb,iwf)*u)-asymp_jasb(ipar+1)
        do i=2,nordb
          psib=psib+b(i+1,isb,iwf)*u**i
        enddo
       elseif(ijas.eq.5) then
        psib=b(1,isb,iwf)*u/(one+b(2,isb,iwf)*u)-asymp_jasb(ipar+1)
        do i=2,nordb
          psib=psib+b(i+1,isb,iwf)*u**i
        enddo
        psib=sspinn*psib
       elseif(ijas.eq.6) then
        psib=b(1,isb,iwf)*u/(one+b(2,isb,iwf)*(1-u))
        do i=2,nordb
          psib=psib+b(i+1,isb,iwf)*u**i
        enddo
        psib=sspinn*psib
      endif

      return
      end
      end module
