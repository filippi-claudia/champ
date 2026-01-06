module multiple_states
      use error,   only: fatal_error
contains
!----------------------------------------------------------------------
      subroutine efficiency_sample(ipass,determ_s,psij,psi2g)

      use mstates2, only: effcm2,effcum
      use mstates_ctrl, only: iefficiency,nstates_psig
      use precision_kinds, only: dp
      use vmc_mod, only: stoj
      implicit none

      integer :: ipass, j, istoj1
      real(dp) :: psi2g, psi2gi, ratio, wi
      real(dp), dimension(*) :: determ_s
      real(dp), dimension(*) :: psij

      if(iefficiency.eq.0) return

      psi2gi=1.d0/psi2g

      istoj1=stoj(1)
      wi=psi2gi
      effcum(1)=effcum(1)+wi
      effcm2(1)=effcm2(1)+wi*wi

      do j=2,nstates_psig
        ratio=determ_s(j)/determ_s(1)*exp(psij(stoj(j))-psij(istoj1))
        wi=ratio*ratio*psi2gi
        effcum(j)=effcum(j)+wi
        effcm2(j)=effcm2(j)+wi*wi
      enddo

!     write(88,*) (effcum(j),effcm2(j),j=1,nstates_psig),(determ_s(j),j=1,nstates_psig),determ_psig
!     write(88,*) (effcum(j)*effcum(j)/effcm2(j)/ipass,j=1,nstates_psig)

      end
!----------------------------------------------------------------------
      subroutine efficiency_init
      use mstates2, only: effcm2,effcum
      use mstates_ctrl, only: nstates_psig
      implicit none

      integer :: j



      do j=1,nstates_psig
        effcum(j)=0.0d0
        effcm2(j)=0.0d0
      enddo

      end
!----------------------------------------------------------------------
      subroutine efficiency_prt(passes)
      use mstates_ctrl, only: iefficiency, nstates_psig
      use mstates2, only: effcm2, effcum
      use contrl_file,    only: ounit, errunit
      use csfs,           only: anormo
      use precision_kinds, only: dp

      implicit none

      integer :: j
      real(dp) :: efficiency, passes




      if(iefficiency.eq.0) return

      write(ounit,*)
      write(ounit,'(''efficiency for multiple states'')')
      do j=1,nstates_psig
        efficiency=effcum(j)*effcum(j)/effcm2(j)/passes
!       write(ounit,*) effcum(j)*effcum(j)/passes,effcm2(j)
        write(ounit,'(''efficiency state '',i4,f8.3)') j,efficiency
!        write(ounit,'(''anorm correction '',i4, E15.8)') j,anormo(j)
      enddo

      end
!----------------------------------------------------------------------
      subroutine efficiency_dump(iu)
      use mstates2, only: effcm2,effcum
      use mstates_ctrl, only: iefficiency,nstates_psig
      implicit none

      integer :: i, iu





      if(iefficiency.eq.0) return

      write(iu) nstates_psig,(effcum(i),effcm2(i),i=1,nstates_psig)

      end
!-----------------------------------------------------------------------
      subroutine efficiency_rstrt(iu)
      use csfs,    only: nstates
      use mstates2, only: effcm2,effcum
      use mstates_ctrl, only: iefficiency,nstates_psig
      use precision_kinds, only: dp

! nstates below is undefined and
! it' also the case in the master branch
! one replacement is
! use csfs, only: nstates
! no idea if that would be correct

      implicit none

      integer :: i, iu
      real(dp) :: different


      if(iefficiency.eq.0) return

      read(iu) nstates_psig,(effcum(i),effcm2(i),i=1,nstates_psig)
      if(nstates_psig.ne.nstates) call fatal_error('EFFICIENCY: different nstates')

      end
!-----------------------------------------------------------------------
end module
