      module multiple_states
      use error, only: fatal_error
      contains
c----------------------------------------------------------------------
      subroutine efficiency_sample(ipass,determ_s,determ_psig)

      use mstates_ctrl, only: iefficiency, nstates_psig
      use mstates2, only: effcm2, effcum
      use precision_kinds, only: dp
      implicit none

      integer :: ipass, j
      real(dp) :: determ_psig, determ_psigi, ratio, wi
      real(dp), dimension(*) :: determ_s




      if(iefficiency.eq.0) return

      determ_psigi=1.d0/determ_psig
c     write(ounit,*) ((determ_s(j)*determ_psigi)**2,j=1,nstates_psig)

      do j=1,nstates_psig
        ratio=determ_s(j)*determ_psigi
        wi=ratio*ratio
        effcum(j)=effcum(j)+wi
        effcm2(j)=effcm2(j)+wi*wi
      enddo

c     write(88,*) (effcum(j),effcm2(j),j=1,nstates_psig),(determ_s(j),j=1,nstates_psig),determ_psig
c     write(88,*) (effcum(j)*effcum(j)/effcm2(j)/ipass,j=1,nstates_psig)

      end
c----------------------------------------------------------------------
      subroutine efficiency_init
      use mstates_ctrl, only: nstates_psig
      use mstates2, only: effcm2, effcum
      implicit none

      integer :: j



      do j=1,nstates_psig
        effcum(j)=0
        effcm2(j)=0
      enddo

      end
c----------------------------------------------------------------------
      subroutine efficiency_prt(passes)
      use mstates_ctrl, only: iefficiency, nstates_psig
      use mstates2, only: effcm2, effcum
      use contrl_file,    only: ounit, errunit
      use precision_kinds, only: dp
      implicit none

      integer :: j
      real(dp) :: efficiency, passes




      if(iefficiency.eq.0) return

      write(ounit,*)
      write(ounit,'(''efficiency for multiple states'')')
      do j=1,nstates_psig
        efficiency=effcum(j)*effcum(j)/effcm2(j)/passes
c       write(ounit,*) effcum(j)*effcum(j)/passes,effcm2(j)
        write(ounit,'(''efficiency state '',i4,f8.3)') j,efficiency
      enddo

      end
c----------------------------------------------------------------------
      subroutine efficiency_dump(iu)
      use mstates_ctrl, only: iefficiency, nstates_psig
      use mstates2, only: effcm2, effcum
      implicit none

      integer :: i, iu





      if(iefficiency.eq.0) return

      write(iu) nstates_psig,(effcum(i),effcm2(i),i=1,nstates_psig)

      end
c-----------------------------------------------------------------------
      subroutine efficiency_rstrt(iu)
      use mstates_ctrl, only: iefficiency, nstates_psig
      use mstates2, only: effcm2, effcum

      ! nstates below is undefined and
      ! it' also the case in the master branch
      ! one replacement is
      ! use csfs, only: nstates
      ! no idea if that would be correct
      use precision_kinds, only: dp

      use csfs, only: nstates
      implicit none

      integer :: i, iu
      real(dp) :: different


      if(iefficiency.eq.0) return

      read(iu) nstates_psig,(effcum(i),effcm2(i),i=1,nstates_psig)
      if(nstates_psig.ne.nstates) call fatal_error('EFFICIENCY: different nstates')

      end
c-----------------------------------------------------------------------
      end module
