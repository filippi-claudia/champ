module prop_dmc
contains
!----------------------------------------------
! Additional properties (<X> etc )
! DMC related routines
! F.Schautz
!----------------------------------------------

      subroutine prop_prt_dmc(iblk,ifinal,wgcum,wgcm2)

      use contrl_file, only: ounit
      use control_dmc, only: dmc_nconf,dmc_nstep
      use multiple_geo, only: MFORCE
      use precision_kinds, only: dp
      use properties, only: MAXPROP
      use prp000,  only: iprop,ipropprt,nprop
      use prp003,  only: cc_nuc,vprop_cm2,vprop_cum
      use system,  only: nelec
!      use contrl, only: nconf, nstep
      implicit none

      integer :: i, iblk, icount, ifinal, iperr
      real(dp) :: dip, diperr, dipx, dipy
      real(dp) :: dipz, evalg_eff
      real(dp) :: rtevalg_eff1
      real(dp), dimension(MFORCE) :: wgcum
      real(dp), dimension(MFORCE) :: wgcm2
      real(dp), dimension(MAXPROP) :: perr
      real(dp), dimension(MAXPROP) :: pav

      character (len=3) pnames(MAXPROP)
      data pnames /'X  ','Y  ','Z  ','XX ','YY ','ZZ '/
      data icount /1/
      save icount

      if(iprop.eq.0.or.(ipropprt.eq.0.and.ifinal.ne.1)) return

! ipropprt 0 no printout
!          1 each iteration full printout
!         >1 after ipropprt iterations reduced printout

      if(icount.ne.ipropprt.and.ifinal.eq.0) then
        icount=icount+1
        return
      endif

      icount=1

      do i=1,nprop
        pav(i)=vprop_cum(i)/wgcum(1)
        if(iblk.eq.1) then
          perr(i)=0
         else
          perr(i)=errg(vprop_cum(i),vprop_cm2(i),1)
          iperr=nint(100000*perr(i))
        endif

        if(ifinal.eq.0) then
          write(ounit,'(a3,25x,f10.5,''  ('',i5,'')'')') pnames(i),pav(i),iperr
         else
          evalg_eff=dmc_nconf*dmc_nstep*rn_eff(wgcum(1),wgcm2(1))
          rtevalg_eff1=dsqrt(evalg_eff-1)
          write(ounit,'(''property '',a3,t17,f12.7,'' +-'',f11.7,f9.5)') pnames(i),pav(i),perr(i),perr(i)*rtevalg_eff1
        endif
      enddo
!....dipole
      write(ounit,50) 'center of nuclear charge: ',cc_nuc
      write(ounit,30)
      dipx=cc_nuc(1)*nelec*2.5417 - pav(1) *2.5417
      dipy=cc_nuc(2)*nelec*2.5417 - pav(2) *2.5417
      dipz=cc_nuc(3)*nelec*2.5417 - pav(3) *2.5417
      dip=dsqrt(dipx**2+dipy**2+dipz**2)
      diperr=dabs (perr(1)*2.5417 * dipx / dip) + &
             dabs (perr(2)*2.5417 * dipy / dip) + &
             dabs (perr(3)*2.5417 * dipz / dip)
      write(ounit,40) 'Dip X ',dipx,perr(1)*2.5417
      write(ounit,40) 'Dip Y ',dipy,perr(2)*2.5417
      write(ounit,40) 'Dip Z ',dipz,perr(3)*2.5417
      write(ounit,40) 'Dip   ',dip,diperr


      30   format('-------- dipole operator averages  ----------')
      40   format(a6,'   ',f16.8,' +- ',f16.8)
      50   format(a25,' ',3f12.8)

      return
contains
        elemental pure function rn_eff(w,w2)
          implicit none
          real(dp), intent(in) :: w, w2
          real(dp)             :: rn_eff
          rn_eff=w**2/w2
        end function
        elemental pure function error(x,x2,w,w2)
          implicit none
          real(dp), intent(in) :: x, x2,w,w2
          real(dp)             :: error
          error=dsqrt(max((x2/w-(x/w)**2)/(rn_eff(w,w2)-1),0.d0))
        end function
        elemental pure function errg(x,x2,i)
          implicit none
          real(dp), intent(in) :: x, x2
          integer, intent(in)  :: i
          real(dp)             :: errg
          errg=error(x,x2,wgcum(i),wgcm2(i))
        end function
      end
!----------------------------------------------------------------------
      subroutine prop_save_dmc(iw)

      use prp000,  only: iprop,nprop
      use prp001,  only: vprop
      use prp002,  only: vprop_old

      implicit none

      integer :: i, iw


      if(iprop.eq.0) return
      do i=1,nprop
       vprop_old(i,iw)=vprop(i)
      enddo
      end
!----------------------------------------------------------------------
      subroutine prop_sum_dmc(p,q,iw)

      use precision_kinds, only: dp
      use prp000,  only: iprop,nprop
      use prp001,  only: vprop
      use prp002,  only: vprop_old
      use prp003,  only: vprop_sum

      implicit none

      integer :: i, iw
      real(dp) :: p, q

      if(iprop.eq.0) return

      do i=1,nprop
       vprop_sum(i)=vprop_sum(i)+p*vprop(i)+q*vprop_old(i,iw)
      enddo
      end
!----------------------------------------------------------------------
      subroutine prop_splitj(iw,iw2)

      use prp000,  only: nprop
      use prp002,  only: vprop_old

      implicit none

      integer :: i, iw, iw2


      do i=1,nprop
       vprop_old(i,iw2)=vprop_old(i,iw)
      enddo
      end
end module
