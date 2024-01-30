module properties_mod
contains
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     additional properties
!     Friedemann Schautz
!
!
!     properties so far:
!     1   2   3   4      5      6
!     <x> <y> <z> <x**2> <y**2> <z**2>
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prop_compute(coord)
      use precision_kinds, only: dp
      use prp000,  only: iprop,nprop
      use prp001,  only: vprop
      use system,  only: nelec
      use periodic, only: ngnorm_sim,ngvec_sim
      implicit none

      integer :: i, jprop, m

      real(dp), dimension(3,*) :: coord

!     electron coordinates

      if(iprop.eq.0) return
      do i=1,nprop
       vprop(i)=0.d0
      enddo

      do i=1,nelec
       do m=1,3
        vprop(m)  = vprop(m)+coord(m,i)
        vprop(3+m)= vprop(3+m) + coord(m,i)**2
       enddo
      enddo

      jprop=6

      call sofk(jprop)

      jprop=6+(ngvec_sim-1)

      call rhok(jprop)

      end

!-----------------------------------------------------------------------
      subroutine prop_init(iflg)
      use prp000,  only: iprop,nprop
      use prp003,  only: vprop_cm2,vprop_cum,vprop_sum

      implicit none

      integer :: i, iflg

      if(iprop.eq.0) return

      do i=1,nprop
       vprop_sum(i)=0.d0
      enddo

! $ iflg = 0: init *cum, *cm2 as well
      if(iflg.gt.0) return
      do i=1,nprop
       vprop_cum(i)=0.d0
       vprop_cm2(i)=0.d0
      enddo
      end

!-----------------------------------------------------------------------
      subroutine prop_cum(w)
      use precision_kinds, only: dp
      use prp000,  only: iprop,nprop
      use prp003,  only: vprop_cm2,vprop_cum,vprop_sum

      implicit none

      integer :: i
      real(dp) :: vprop_now, w

      if(iprop.eq.0) return
      do i=1,nprop
       vprop_now = vprop_sum(i)/w
       vprop_cm2(i)=vprop_cm2(i)+ vprop_sum(i)*vprop_now
       vprop_cum(i)=vprop_cum(i)+ vprop_sum(i)
      enddo
      end

!-----------------------------------------------------------------------
      subroutine prop_avrg(wcum,iblk,pav,perr)
      use precision_kinds, only: dp
      use prp000,  only: iprop,nprop
      use prp003,  only: vprop_cm2,vprop_cum

      implicit none

      integer :: i, iblk
      real(dp) :: wcum, x, x2
      real(dp), dimension(nprop) :: pav
      real(dp), dimension(nprop) :: perr

      if(iprop.eq.0) return
      do i=1,nprop
       perr(i)=err(vprop_cum(i),vprop_cm2(i))
       pav(i)=vprop_cum(i)/wcum
      enddo
contains
        elemental pure function err(x,x2)
        implicit none
        real(dp) err
        real(dp), intent(in) :: x, x2
        err = dsqrt(abs(x2/wcum-(x/wcum)**2)/iblk)
        end function
      end
!-----------------------------------------------------------------------
      subroutine prop_dump(iu)
      use prp000,  only: iprop,nprop
      use prp003,  only: vprop_cm2,vprop_cum

      implicit none

      integer :: i, iu



      if(iprop.eq.0) return
      write(iu) nprop
      write(iu) (vprop_cum(i),vprop_cm2(i),i=1,nprop)
      end
      subroutine prop_rstrt(iu)
      use prp000,  only: iprop,nprop
      use prp003,  only: vprop_cm2,vprop_cum

      implicit none

      integer :: i, iu



      if(iprop.eq.0) return
      read(iu) nprop
      read(iu) (vprop_cum(i),vprop_cm2(i),i=1,nprop)
      end
!-----------------------------------------------------------------------
      subroutine prop_fin(passes,iblk,efin,eerr)
      use contrl_file, only: ounit
      use precision_kinds, only: dp
      use prp000,  only: iprop,ipropprt
      implicit none

      integer :: iblk, ipropprt_sav
      real(dp) :: eerr, efin, passes


      if(iprop.eq.0) return
      write(ounit,'(''--- additional properties ---'')')
      ipropprt_sav=ipropprt
      ipropprt=-1
      call prop_prt(passes,iblk,ounit)
      ipropprt=ipropprt_sav

      end
!-----------------------------------------------------------------------
      subroutine prop_prt(w,iblk,iu)
      use m_icount, only: icount_prop
      use periodic, only: ngnorm_sim, gvec_sim, ngvec_sim
      use precision_kinds, only: dp
      use prp000,  only: iprop,ipropprt,nprop
      use prp003,  only: cc_nuc
      use system,  only: nelec

      implicit none

      integer :: i, iblk, iu
      real(dp) :: dble, dip, diperr, dipx
      real(dp) :: dipy, dipz, norm_aux, w
      real(dp), dimension(nprop) :: pav
      real(dp), dimension(nprop) :: perr

! compute averages and print then out

! ipropprt 0 no printout
!          1 each iteration full printout
!         >1 after ipropprt iterations reduced printout
!         -1 force printout

      if(iprop.eq.0.or.ipropprt.eq.0) return

      if(ipropprt.gt.0.and.icount_prop.ne.ipropprt) then
        icount_prop=icount_prop+1
        return
      endif

      icount_prop=1

      call prop_avrg(w,iblk,pav(1),perr(1))


      write(iu,10)
      write(iu,20) 'X  ',pav(1),perr(1),pav(1)/dble(nelec),perr(1) &
           /dble(nelec)
      write(iu,20) 'Y  ',pav(2),perr(2),pav(2)/dble(nelec),perr(2) &
           /dble(nelec)
      write(iu,20) 'Z  ',pav(3),perr(3),pav(3)/dble(nelec),perr(3) &
           /dble(nelec)
      write(iu,20) 'XX ',pav(4),perr(4),pav(4)/dble(nelec),perr(4) &
           /dble(nelec)
      write(iu,20) 'YY ',pav(5),perr(5),pav(5)/dble(nelec),perr(5) &
           /dble(nelec)
      write(iu,20) 'ZZ ',pav(6),perr(6),pav(6)/dble(nelec),perr(6) &
           /dble(nelec)

!....dipole
      write(iu,50) 'center of nuclear charge: ',cc_nuc
      write(iu,30)
      dipx=cc_nuc(1)*nelec*2.5417 - pav(1) *2.5417
      dipy=cc_nuc(2)*nelec*2.5417 - pav(2) *2.5417
      dipz=cc_nuc(3)*nelec*2.5417 - pav(3) *2.5417
      dip=dsqrt(dipx**2+dipy**2+dipz**2)
      diperr=dabs (perr(1)*2.5417 * dipx / dip) + &
             dabs (perr(2)*2.5417 * dipy / dip) + &
             dabs (perr(3)*2.5417 * dipz / dip)
      write(iu,40) 'Dip X ',dipx,perr(1)*2.5417
      write(iu,40) 'Dip Y ',dipy,perr(2)*2.5417
      write(iu,40) 'Dip Z ',dipz,perr(3)*2.5417
      write(iu,40) 'Dip   ',dip,diperr

      do i=6+1,6+ngvec_sim-1
        call gnormf(3,gvec_sim(1,i-5), norm_aux)
        write(iu,'(''s(k)     '',t17,f12.7,f12.7,'' +-'' &
        ,f12.7,f12.7)') norm_aux,pav(i),perr(i), &
        pav(i)-(pav(i+ngvec_sim-1)**2)-(pav(i+2*(ngvec_sim-1))**2)
      enddo
      do i=6+ngvec_sim,6+2*(ngvec_sim-1)
        call gnormf(3,gvec_sim(1,i-5-ngvec_sim+1), norm_aux)
        write(iu,'(''cos(kr)  '',t17,f12.7,f12.7,'' +-'' &
        ,f12.7)') norm_aux,pav(i),perr(i)
      enddo
      do i=6+2*(ngvec_sim-1)+1,nprop
        call gnormf(3,gvec_sim(1,i-5-2*(ngvec_sim-1)), norm_aux)
        write(iu,'(''sin(kr)  '',t17,f12.7,f12.7,'' +-'' &
        ,f12.7)') norm_aux,pav(i),perr(i)
      enddo

      10 format('-------- property operator averages  ----------')
      20 format(a3,'   ',f16.8,' +- ',f16.8,' ( ',f16.8,' +- ',f16.8,' )')
      30 format('-------- dipole operator averages  ----------')
      40 format(a6,'   ',f16.8,' +- ',f16.8)
      50 format(a25,' ',3f12.8)

      end

!*********************************************************************
        subroutine prop_cc_nuc(znuc,cent,iwctype,mctype,mcent, &
        ncent,cc_nuc)
!*********************************************************************


      implicit none

        integer  mctype,mcent,ncent
        integer  iwctype(mcent)
        real(8) znuc(mctype),cent(3,mcent)
        real(8) cc_nuc(3), tmp
        integer i,id

        cc_nuc(1)=0.d0
        cc_nuc(2)=0.d0
        cc_nuc(3)=0.d0
        tmp=0.d0
        do i=1,ncent
          id=znuc(iwctype(i))
          if(znuc(iwctype(i)).gt.2) id=znuc(iwctype(i))+2
          cc_nuc(1)=cc_nuc(1)+znuc(iwctype(i))*cent(1,i)
          cc_nuc(2)=cc_nuc(2)+znuc(iwctype(i))*cent(2,i)
          cc_nuc(3)=cc_nuc(3)+znuc(iwctype(i))*cent(3,i)
          tmp=tmp+znuc(iwctype(i))
        enddo
        cc_nuc(1)=cc_nuc(1)/tmp
        cc_nuc(2)=cc_nuc(2)/tmp
        cc_nuc(3)=cc_nuc(3)/tmp
!        write (*,*) 'Center of nuclear charge:', cc_nuc(:),tmp


        return
        end

!-----------------------------------------------------------------------

      subroutine sofk(ip)
! Written by Edgar Landinez and Saverio Moroni

      use prp001,   only: vprop
      use ewald,    only: cos_e_sum_sim, sin_e_sum_sim
      use periodic, only: igmult_sim, ngnorm_sim, ngvec_sim

      use system, only: nelec
      use precision_kinds, only: dp
      implicit none

      integer :: im, ip, ivec, k,ik
      real(dp) :: skcum
      real(dp), dimension(ngvec_sim) :: cos1_sum
      real(dp), dimension(ngvec_sim) :: cos2_sum
      real(dp), dimension(ngvec_sim) :: sin1_sum
      real(dp), dimension(ngvec_sim) :: sin2_sum

      cos1_sum=cos_e_sum_sim
      cos2_sum=cos_e_sum_sim
      sin1_sum=sin_e_sum_sim
      sin2_sum=sin_e_sum_sim

!     print*,"sofk, ngnorm", ngnorm_sim

      !vprop(ip+1)=(cos1_sum(1)*cos2_sum(1)+sin1_sum(1)*sin2_sum(1))/nelec

!      ivec=1
!      do k=2,ngnorm_sim
!         skcum=0.d0
!         do im=1,igmult_sim(k)
!            ivec=ivec+1
!            skcum=skcum+(cos1_sum(ivec)*cos2_sum(ivec)+sin1_sum(ivec)*sin2_sum(ivec))
!         enddo
!         vprop(ip+k-1)=skcum/(nelec*igmult_sim(k))
!      enddo

      ivec=1

      ik=0
      do k=2,ngnorm_sim
         do im=1,igmult_sim(k)
            ik=ik+1
            ivec=ivec+1
            vprop(ip+ik)=(cos1_sum(ivec)*cos2_sum(ivec)+sin1_sum(ivec)*sin2_sum(ivec))
         enddo
      enddo


      
      return
      end

! ------------------------------------------------------------------------------------------------
      subroutine rhok(ip)
! Written by Edgar Landinez and Saverio Moroni

      use prp001,   only: vprop
      use ewald,    only: cos_e_sum_sim, sin_e_sum_sim
      use periodic, only: igmult_sim, ngnorm_sim, ngvec_sim

      use system, only: nelec
      use precision_kinds, only: dp
      implicit none

      integer :: im, ip, ivec, k, ik 
      real(dp) :: skcum

!     print*,"sofk, ngnorm", ngnorm_sim

      !vprop(ip+1)=(cos1_sum(1)*cos2_sum(1)+sin1_sum(1)*sin2_sum(1))/nelec

      ivec=1
      ik=2
      do k=2,ngvec_sim
         vprop(ip+ik-1)=cos_e_sum_sim(k)
         ik=ik+1
      enddo

      do k=2,ngvec_sim
         vprop(ip+ik-1)=sin_e_sum_sim(k)
         ik=ik+1
      enddo
      return
      end
! ------------------------------------------------------------------------------------------------
      subroutine gnormf(n,gvec,gnorm)
      use precision_kinds, only: dp
      implicit none
      integer :: i, n
      real(dp),  dimension(n)::gvec
      real(dp) :: gnorm
      gnorm=0.d0
      do i=1,n
         gnorm=gnorm+gvec(i)*gvec(i)
      enddo
      gnorm=dsqrt(gnorm)

      return
      end

end module
