      subroutine efficiency_sample(ipass,determ_s,determ_psig)
      use mstates_ctrl, only: iefficiency, nstates_psig
      use mstates2, only: effcm2, effcum

      implicit real*8(a-h,o-z)

      dimension determ_s(*)

      if(iefficiency.eq.0) return

      determ_psigi=1.0d0/determ_psig

      do j=1,nstates_psig
         ratio=determ_s(j)*determ_psigi
         wi=ratio*ratio
         effcum(j)=effcum(j)+wi
         effcm2(j)=effcm2(j)+wi*wi
      enddo

      end subroutine

c----------------------------------------------------------------------

      subroutine efficiency_init
      use mstates_ctrl, only: nstates_psig
      use mstates2, only: effcm2, effcum

      implicit real*8(a-h,o-z)

      do j=1,nstates_psig
         effcum(j)=0.0d0
         effcm2(j)=0.0d0
      enddo

      end subroutine

c----------------------------------------------------------------------

      subroutine efficiency_prt(passes)
      use mstates_ctrl, only: iefficiency, nstates_psig
      use mstates2, only: effcm2, effcum

      implicit real*8(a-h,o-z)

      if(iefficiency.eq.0) return

      write(6,*)
      write(6,'(''efficiency for multiple states'')')
      do j=1,nstates_psig
         efficiency=effcum(j)*effcum(j)/effcm2(j)/passes
         write(6,'(''efficiency state '',i4,f8.3)') j,efficiency
      enddo

      end subroutine

c----------------------------------------------------------------------

      subroutine efficiency_dump(iu)
      use mstates_ctrl, only: iefficiency, nstates_psig
      use mstates2, only: effcm2, effcum

      implicit real*8(a-h,o-z)

      if(iefficiency.eq.0) return

      write(iu) nstates_psig,(effcum(i),effcm2(i),i=1,nstates_psig)

      end subroutine

c-----------------------------------------------------------------------

      subroutine efficiency_rstrt(iu)
      use mstates_ctrl, only: iefficiency, nstates_psig
      use mstates2, only: effcm2, effcum
      ! nstates below is undefined and 
      ! it' also the case in the master branch
      ! one replacement is:
      use csfs, only: nstates
      ! no idea if that would be correct
      implicit real*8(a-h,o-z)

      if(iefficiency.eq.0) return

      read(iu) nstates_psig,(effcum(i),effcm2(i),i=1,nstates_psig)
      if(nstates_psig.ne.nstates) then
         call fatal('EFFICIENCY: different nstates')
      endif

      end subroutine
