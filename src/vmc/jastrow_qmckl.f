      module jastrow_qmckl_mod
      contains

      subroutine jastrow_init_qmckl(iwft)

      use jastrow, only: norda,nordb,nordc,a4,b,c,scalek
      use jastrow4_mod, only: nterms4
      use precision_kinds, only: dp
      use scale_dist_mod, only: scale_dist2,switch_scale2
      use system,  only: iwctype,ncent,nelec,nup,nctype

      use qmckl_data

      implicit none
      integer, intent(in) :: iwft
      integer(qmckl_exit_code) :: rc, rc2
      double precision :: scalek_en(nctype)
      integer*8 :: dimc, norda_l, nordb_l

      norda_l = max(norda, 1)
      nordb_l = max(nordb, 1)

      rc = qmckl_set_jastrow_aord_num (qmckl_ctx(iwft), norda_l*1_8)
      rc2 = qmckl_check(qmckl_ctx(iwft), rc)
      if (rc /= QMCKL_SUCCESS) stop __LINE__

      rc = qmckl_set_jastrow_bord_num (qmckl_ctx(iwft), (nordb_l+1)*1_8)
      rc2 = qmckl_check(qmckl_ctx(iwft), rc)
      if (rc /= QMCKL_SUCCESS) stop __LINE__

      rc = qmckl_set_jastrow_cord_num (qmckl_ctx(iwft), nordc*1_8)
      rc2 = qmckl_check(qmckl_ctx(iwft), rc)
      if (rc /= QMCKL_SUCCESS) stop __LINE__

      rc = qmckl_set_jastrow_type_nucl_num (qmckl_ctx(iwft), nctype*1_8)
      rc2 = qmckl_check(qmckl_ctx(iwft), rc)
      if (rc /= QMCKL_SUCCESS) stop __LINE__

      rc = qmckl_set_jastrow_type_nucl_vector (qmckl_ctx(iwft),
     &      int(iwctype(:),8), size(iwctype,1)*1_8)
      rc2 = qmckl_check(qmckl_ctx(iwft), rc)
      if (rc /= QMCKL_SUCCESS) stop __LINE__

      rc = qmckl_set_jastrow_a_vector (qmckl_ctx(iwft),
     &      a4(1:norda_l+1,1:nctype,iwft), (norda_l+1)*nctype)
      rc2 = qmckl_check(qmckl_ctx(iwft), rc)
      if (rc /= QMCKL_SUCCESS) stop __LINE__

      rc = qmckl_set_jastrow_b_vector (qmckl_ctx(iwft),
     &      b(1:(nordb_l+1),1:nctype,iwft), (nordb_l+1)*nctype*1_8)
      rc2 = qmckl_check(qmckl_ctx(iwft), rc)
      if (rc /= QMCKL_SUCCESS) stop __LINE__

      dimc = nterms4(nordc)
      rc = qmckl_set_jastrow_c_vector (qmckl_ctx(iwft),
     &      c(1:dimc,1:nctype,iwft), dimc*nctype*1_8)
      rc2 = qmckl_check(qmckl_ctx(iwft), rc)
      if (rc /= QMCKL_SUCCESS) stop __LINE__

      rc = qmckl_set_jastrow_rescale_factor_ee (qmckl_ctx(iwft), scalek(iwft))
      rc2 = qmckl_check(qmckl_ctx(iwft), rc)
      if (rc /= QMCKL_SUCCESS) stop __LINE__

      scalek_en(:) = scalek(iwft)
      rc = qmckl_set_jastrow_rescale_factor_en (qmckl_ctx(iwft),
     &      scalek_en, nctype*1_8)
      rc2 = qmckl_check(qmckl_ctx(iwft), rc)
      if (rc /= QMCKL_SUCCESS) stop __LINE__

      end subroutine


!---------------------------------------------------------------------------------

      subroutine jastrow_qmckl_common(x,v,d2,value)

      use qmckl_data
      use multiple_geo, only: iwf
      use system,  only: iwctype,ncent,nelec,nup,nctype

      implicit none
      double precision, intent(in)  :: x(3,*)
      double precision, intent(out) :: v(3,*)
      double precision, intent(out) :: d2
      double precision, intent(out) :: value

      integer(qmckl_exit_code) :: rc, rc2
      double precision :: jen(1), jee(1), jeen(1)
!      double precision :: ven(, vee, veen

      rc = qmckl_set_electron_coord(qmckl_ctx(iwf), 'N', 1_8, x, nelec*3_8)
      rc2 = qmckl_check(qmckl_ctx(iwf), rc)
      if (rc /= QMCKL_SUCCESS) stop __LINE__

      rc = qmckl_get_jastrow_factor_ee(qmckl_ctx(iwf), jee, 1_8)
      rc2 = qmckl_check(qmckl_ctx(iwf), rc)
      if (rc /= QMCKL_SUCCESS) stop __LINE__

      rc = qmckl_get_jastrow_factor_en(qmckl_ctx(iwf), jen, 1_8)
      rc2 = qmckl_check(qmckl_ctx(iwf), rc)
      if (rc /= QMCKL_SUCCESS) stop __LINE__

!     rc = qmckl_get_jastrow_factor_een(qmckl_ctx(iwf), jeen, 1_8)
!     rc2 = qmckl_check(qmckl_ctx(iwf), rc)
!     if (rc /= QMCKL_SUCCESS) stop __LINE__

      value = jee(1) + jen(1) !+ jeen(1)
      print *, '-------------- Fee:', jee(1), jen(1), value

      end subroutine

!---------------------------------------------------------------------------------

      subroutine jastrow_qmckl(x,v,d2,div_vj,value)

      use qmckl_data
      use multiple_geo, only: iwf
      use jastrow4_mod, only: jastrow4

      use jastrow, only: asymp_jasa,asymp_jasb

      use system,  only: iwctype,ncent,nelec,nup,nctype
      implicit none
      double precision, intent(in)  :: x(3,*)
      double precision, intent(out) :: v(3,*)
      double precision, intent(out) :: d2
      double precision, intent(out) :: div_vj(*)
      double precision, intent(out) :: value

      integer(qmckl_exit_code) :: rc, rc2

      integer :: it,i 
      double precision :: xx(10000)

      ! TODO CHECK
!      print *, 'scalek = ', scalek(1)
!      write(*,'(''asympa='',10f10.6)') (asymp_jasa(it),it=1,nctype)
!      write(*,'(''asympb='',10f10.6)') (asymp_jasb(it),it=1,2)
!
!      rc = qmckl_get_jastrow_asymp_jasa(qmckl_ctx(1), xx, nctype)
!      rc2 = qmckl_check(qmckl_ctx(1), rc)
!      if (rc /= QMCKL_SUCCESS) stop __LINE__
!      write(*,'(''asympa='',10f10.6)') (xx(it),it=1,nctype)
!
!      rc = qmckl_get_jastrow_asymp_jasb(qmckl_ctx(1), xx, 2)
!      rc2 = qmckl_check(qmckl_ctx(1), rc)
!      if (rc /= QMCKL_SUCCESS) stop __LINE__
!      write(*,'(''asympb='',10f10.6)') (xx(it),it=1,2)
!

      call jastrow4(x,v,d2,div_vj,value)
      print *, 'before', value

      call jastrow_qmckl_common(x,v,d2,value)
      print *, 'after', value
      stop

      end subroutine

!---------------------------------------------------------------------------------

      subroutine jastrowe_qmckl(iel,x,v,d2,value,iflag)

      use qmckl_data
      use multiple_geo, only: iwf
      use jastrow4e_mod, only: jastrow4e

      implicit none
      integer, intent(in)           :: iel, iflag
      double precision, intent(in)  :: x(3,*)
      double precision, intent(out) :: v(3,*)
      double precision, intent(out) :: d2
      double precision, intent(out) :: value

      integer(qmckl_exit_code) :: rc, rc2
!      double precision :: ven(, vee, veen

      if (iflag /= 0) then
        print *, 'jastrowe_qmckl only implemented for iflag=0'
      endif

      ! TODO CHECK
      call jastrow4e(iel,x,v,d2,value,iflag)
      print *, 'before', value

      call jastrow_qmckl_common(x,v,d2,value)
      print *, 'after', value

      end subroutine

!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

!      integer :: i, ic, ij, im1, iord
!      integer :: ipar, isb, it, j
!      integer :: k, l, l_hi, ll
!      integer :: m, n
!      real(dp) :: bot, bot2, boti, botii, botu
!      real(dp) :: botuu, d2, dd1, dd10
!      real(dp) :: dd2, dd7, dd8, dd9
!      real(dp) :: fc, fee, feeu, feeuu
!      real(dp) :: fen, feni, feni_save, fenii
!      real(dp) :: fenii_save, fi, fii, fj
!      real(dp) :: fjj, fsum, fu, fui
!      real(dp) :: fuj, fuu, ri, rij
!      real(dp) :: rj, s, t, term
!      real(dp) :: top, topi, topii, topu
!      real(dp) :: topuu, u2mst, u2pst, value
!      real(dp), dimension(3, *) :: x
!      real(dp), dimension(3, *) :: v
!      real(dp), dimension(*) :: div_vj
!      real(dp), dimension(-2:nordj) :: uu
!      real(dp), dimension(-2:nordj) :: ss
!      real(dp), dimension(-2:nordj) :: tt
!      real(dp), dimension(-2:nordj) :: rri
!      real(dp), dimension(-2:nordj) :: rrj
!      real(dp), parameter :: half = .5d0
!      real(dp), parameter :: eps = 1.d-12
!
!      fsum=0
!      do i=-2,-1
!        uu(i)=0
!        ss(i)=0
!        tt(i)=0
!        rri(i)=0
!        rrj(i)=0
!      enddo
!      uu(0)=1
!      ss(0)=2
!      tt(0)=1
!      rri(0)=1
!      rrj(0)=1
!
!      if (nelec.lt.2) goto 65
!
!c e-e and e-e-n terms
!      ij=0
!      do i=2,nelec
!      im1=i-1
!      do j=1,im1
!      ij=ij+1
!
!      fso(i,j)=0
!      d2ijo(i,j)=0
!      do k=1,3
!        fijo(k,i,j)=0
!        fijo(k,j,i)=0
!      enddo
!
!      sspinn=1
!      isb=1
!      ipar=0
!      if(i.le.nup .or. j.gt.nup) then
!        if(nspin2b.eq.2) then
!          isb=2
!         elseif(nocuspb.eq.0) then
!          sspinn=half
!        endif
!        ipar=1
!      endif
!
!      rij=r_ee(ij)
!
!      call scale_dist2(rij,uu(1),dd1,dd2,1)
!c     write(ounit,'(''rij,u in ee'',2f9.5)') rij,uu(1)
!
!c Check rij after scaling because uu(1) used in e-e-n terms too
!      if(rij.gt.cutjas) goto 22
!
!      top=sspinn*b(1,isb,iwf)*uu(1)
!      topu=sspinn*b(1,isb,iwf)
!      topuu=0
!
!      if(ijas.eq.4.or.ijas.eq.5) then
!        bot=1+b(2,isb,iwf)*uu(1)
!        botu=b(2,isb,iwf)
!       elseif(ijas.eq.6) then
!        bot=1+b(2,isb,iwf)*(1-uu(1))
!        botu=-b(2,isb,iwf)
!      endif
!      botuu=0
!      bot2=bot*bot
!
!      fee=top/bot-asymp_jasb(ipar+1)
!      feeu=topu/bot-botu*top/bot2
!      feeuu=topuu-(botuu*top+2*botu*topu)/bot+2*botu**2*top/bot2
!      feeuu=feeuu/bot
!
!      do iord=2,nordb
!        uu(iord)=uu(1)*uu(iord-1)
!        if(ijas.eq.4) then
!          fee=fee+b(iord+1,isb,iwf)*uu(iord)
!          feeu=feeu+b(iord+1,isb,iwf)*iord*uu(iord-1)
!          feeuu=feeuu+b(iord+1,isb,iwf)*iord*(iord-1)*uu(iord-2)
!         elseif(ijas.eq.5.or.ijas.eq.6) then
!          fee=fee+sspinn*b(iord+1,isb,iwf)*uu(iord)
!          feeu=feeu+sspinn*b(iord+1,isb,iwf)*iord*uu(iord-1)
!          feeuu=feeuu+sspinn*b(iord+1,isb,iwf)*iord*(iord-1)*uu(iord-2)
!        endif
!      enddo
!
!      feeuu=feeuu*dd1*dd1+feeu*dd2
!      feeu=feeu*dd1/rij
!
!      fso(i,j)=fso(i,j)+fee
!      do k=1,3
!        fijo(k,i,j)= fijo(k,i,j) + feeu*rvec_ee(k,ij)
!        fijo(k,j,i)= fijo(k,j,i) - feeu*rvec_ee(k,ij)
!      enddo
!      d2ijo(i,j)=d2ijo(i,j)+2*(feeuu+2*feeu)
!
!c There are no C terms to order 1.
!   22 if(nordc.le.1) goto 55
!
!      if(isc.ge.12) call scale_dist2(rij,uu(1),dd1,dd2,3)
!      if(ijas.eq.4.or.ijas.eq.5) then
!        call switch_scale2(uu(1),dd1,dd2)
!        do iord=2,nordc
!          uu(iord)=uu(1)*uu(iord-1)
!        enddo
!      endif
!c     write(ounit,'(''rij,u in een'',2f12.9)') rij,uu(1)
!
!      do ic=1,ncent
!        it=iwctype(ic)
!
!        ri=r_en(i,ic)
!        rj=r_en(j,ic)
!
!        if(ri.gt.cutjas .or. rj.gt.cutjas) goto 50
!        do k=1,3
!          if(abs(rshift(k,i,ic)-rshift(k,j,ic)).gt.eps) goto 50
!        enddo
!
!        call scale_dist2(ri,rri(1),dd7,dd9,2)
!        call scale_dist2(rj,rrj(1),dd8,dd10,2)
!
!        if(ijas.eq.4.or.ijas.eq.5) then
!          call switch_scale2(rri(1),dd7,dd9)
!          call switch_scale2(rrj(1),dd8,dd10)
!        endif
!c     write(ounit,'(''ri,rri in een'',2f12.9)') ri,rri(1)
!
!        s=ri+rj
!        t=ri-rj
!c       u2mt2=rij*rij-t*t
!        u2pst=rij*rij+s*t
!        u2mst=rij*rij-s*t
!c       s2mu2=s*s-rij*rij
!c       s2mt2=s*s-t*t
!
!        do iord=1,nordc
!          rri(iord)=rri(1)*rri(iord-1)
!          rrj(iord)=rrj(1)*rrj(iord-1)
!          ss(iord)=rri(iord)+rrj(iord)
!          tt(iord)=rri(iord)*rrj(iord)
!        enddo
!
!        fc=0
!        fu=0
!        fuu=0
!        fi=0
!        fii=0
!        fj=0
!        fjj=0
!        fui=0
!        fuj=0
!        ll=0
!        do n=2,nordc
!          do k=n-1,0,-1
!            if(k.eq.0) then
!              l_hi=n-k-2
!             else
!              l_hi=n-k
!            endif
!            do l=l_hi,0,-1
!              m=(n-k-l)/2
!              if(2*m.eq.n-k-l) then
!                ll=ll+1
!                fc=fc+c(ll,it,iwf)*uu(k)*ss(l)*tt(m)
!                fu=fu+c(ll,it,iwf)*k*uu(k-1)*ss(l)*tt(m)
!                fuu=fuu+c(ll,it,iwf)*k*(k-1)*uu(k-2)*ss(l)*tt(m)
!                fi=fi+c(ll,it,iwf)*uu(k)
!     &          *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
!                fii=fii+c(ll,it,iwf)*uu(k)
!     &          *((l+m)*(l+m-1)*rri(l+m-2)*rrj(m)
!     &          +m*(m-1)*rri(m-2)*rrj(l+m))
!                fj=fj+c(ll,it,iwf)*uu(k)
!     &          *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))
!                fjj=fjj+c(ll,it,iwf)*uu(k)
!     &          *((l+m)*(l+m-1)*rrj(l+m-2)*rri(m)
!     &          +m*(m-1)*rrj(m-2)*rri(l+m))
!                fui=fui+c(ll,it,iwf)*k*uu(k-1)
!     &          *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
!                fuj=fuj+c(ll,it,iwf)*k*uu(k-1)
!     &          *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))
!              endif
!c     write(ounit,'(''rij,ri,rj'',9f10.5)') rij,ri,rj,uu(1),rri(1),rrj(1)
!            enddo
!          enddo
!        enddo
!
!        fuu=fuu*dd1*dd1+fu*dd2
!        fu=fu*dd1/rij
!
!        fui=fui*dd1*dd7
!        fuj=fuj*dd1*dd8
!
!        fii=fii*dd7*dd7+fi*dd9
!        fjj=fjj*dd8*dd8+fj*dd10
!        fi=fi*dd7/ri
!        fj=fj*dd8/rj
!
!        fso(i,j)=fso(i,j) + fc
!
!        fijo(1,i,j)=fijo(1,i,j) + fi*rvec_en(1,i,ic)+fu*rvec_ee(1,ij)
!        fijo(2,i,j)=fijo(2,i,j) + fi*rvec_en(2,i,ic)+fu*rvec_ee(2,ij)
!        fijo(3,i,j)=fijo(3,i,j) + fi*rvec_en(3,i,ic)+fu*rvec_ee(3,ij)
!        fijo(1,j,i)=fijo(1,j,i) + fj*rvec_en(1,j,ic)-fu*rvec_ee(1,ij)
!        fijo(2,j,i)=fijo(2,j,i) + fj*rvec_en(2,j,ic)-fu*rvec_ee(2,ij)
!        fijo(3,j,i)=fijo(3,j,i) + fj*rvec_en(3,j,ic)-fu*rvec_ee(3,ij)
!c       write(ounit,'(''i,j,fijo2='',2i5,9d12.4)') i,j,(fijo(k,i,j),k=1,3)
!
!        d2ijo(i,j)=d2ijo(i,j) + 2*(fuu + 2*fu) + fui*u2pst/(ri*rij)
!     &  + fuj*u2mst/(rj*rij) + fii + 2*fi + fjj + 2*fj
!
!  50  continue
!      enddo
!
!  55  fsum=fsum+fso(i,j)
!      v(1,i)=v(1,i)+fijo(1,i,j)
!      v(2,i)=v(2,i)+fijo(2,i,j)
!      v(3,i)=v(3,i)+fijo(3,i,j)
!      v(1,j)=v(1,j)+fijo(1,j,i)
!      v(2,j)=v(2,j)+fijo(2,j,i)
!      v(3,j)=v(3,j)+fijo(3,j,i)
!      div_vj(i)=div_vj(i)+d2ijo(i,j)/2
!      div_vj(j)=div_vj(j)+d2ijo(i,j)/2
!      d2=d2+d2ijo(i,j)
!      enddo
!      enddo
!
!c e-n terms
!  65  do i=1,nelec
!
!        fso(i,i)=0
!        fijo(1,i,i)=0
!        fijo(2,i,i)=0
!        fijo(3,i,i)=0
!        d2ijo(i,i)=0
!
!        do ic=1,ncent
!          it=iwctype(ic)
!
!          ri=r_en(i,ic)
!          if(ri.gt.cutjas) goto 80
!
!          call scale_dist2(ri,rri(1),dd7,dd9,1)
!c     write(ounit,'(''ri,rri in en'',2f9.5)') ri,rri(1)
!
!          top=a4(1,it,iwf)*rri(1)
!          topi=a4(1,it,iwf)
!          topii=0
!
!          if(ijas.eq.4.or.ijas.eq.5) then
!            bot=1+a4(2,it,iwf)*rri(1)
!            boti=a4(2,it,iwf)
!           elseif(ijas.eq.6) then
!            bot=1+a4(2,it,iwf)*(1-rri(1))
!            boti=-a4(2,it,iwf)
!          endif
!          botii=0
!          bot2=bot*bot
!
!          fen=top/bot-asymp_jasa(it)
!          feni=topi/bot-boti*top/bot2
!          fenii=topii-(botii*top+2*boti*topi)/bot+2*boti**2*top/bot2
!          fenii=fenii/bot
!
!        do iord=2,norda
!          rri(iord)=rri(1)**iord
!          fen=fen+a4(iord+1,it,iwf)*rri(iord)
!          feni=feni+a4(iord+1,it,iwf)*iord*rri(iord-1)
!          fenii=fenii+a4(iord+1,it,iwf)*iord*(iord-1)*rri(iord-2)
!        enddo
!
!          feni_save=feni
!          fenii_save=fenii
!          fenii=fenii*dd7*dd7+feni*dd9
!          feni=feni*dd7/ri
!
!          fso(i,i)=fso(i,i)+fen
!
!          fijo(1,i,i)=fijo(1,i,i) + feni*rvec_en(1,i,ic)
!          fijo(2,i,i)=fijo(2,i,i) + feni*rvec_en(2,i,ic)
!          fijo(3,i,i)=fijo(3,i,i) + feni*rvec_en(3,i,ic)
!c         write(ounit,'(''fijo='',9d12.4)') (fijo(k,i,i),k=1,3),feni,rvec_en(1,i,ic)
!
!          d2ijo(i,i) = d2ijo(i,i) + fenii + 2*feni
!
!          if(iforce_analy.eq.1) call da_jastrow4(iwf,i,ic,it,rvec_en(1,i,ic),ri,rri,feni_save,fenii_save,dd7,dd9)
!   80   continue
!        enddo
!
!        fsum=fsum+fso(i,i)
!        v(1,i)=v(1,i)+fijo(1,i,i)
!        v(2,i)=v(2,i)+fijo(2,i,i)
!        v(3,i)=v(3,i)+fijo(3,i,i)
!c       write(ounit,'(''v='',9d12.4)') (v(k,i),k=1,3)
!        div_vj(i)=div_vj(i)+d2ijo(i,i)
!        d2=d2+d2ijo(i,i)
!      enddo
!
!      if(ijas.eq.6) then
!        term=1/(c1_jas6*scalek(iwf))
!        fsum=term*fsum
!        d2=term*d2
!        do i=1,nelec
!          div_vj(i)=term*div_vj(i)
!          do k=1,3
!            v(k,i)=term*v(k,i)
!          enddo
!          do j=1,nelec
!            d2ijo(i,j)=term*d2ijo(i,j)
!            do k=1,3
!              fijo(k,i,j)=term*fijo(k,i,j)
!            enddo
!          enddo
!        enddo
!      endif
!
!      fsumo=fsum
!      d2o=d2
!      do i=1,nelec
!        fjo(1,i)=v(1,i)
!        fjo(2,i)=v(2,i)
!        fjo(3,i)=v(3,i)
!      enddo
!
!      value=fsum
!
!      return
!      end

      end module
