c----------------------------------------------------------------------
      subroutine efficiency_sample(ipass,determ_s,determ_psig)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'optci.h'
      include 'mstates.h'

      dimension determ_s(*)

      if(iefficiency.eq.0) return

      determ_psigi=1.d0/determ_psig
c     write(6,*) ((determ_s(j)*determ_psigi)**2,j=1,nstates_psig)

      do 100 j=1,nstates_psig
        ratio=determ_s(j)*determ_psigi
        wi=ratio*ratio
        effcum(j)=effcum(j)+wi
  100   effcm2(j)=effcm2(j)+wi*wi

c     write(88,*) (effcum(j),effcm2(j),j=1,nstates_psig),(determ_s(j),j=1,nstates_psig),determ_psig
c     write(88,*) (effcum(j)*effcum(j)/effcm2(j)/ipass,j=1,nstates_psig)

      end
c----------------------------------------------------------------------
      subroutine efficiency_init
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'optci.h'
      include 'mstates.h'
      do 100 j=1,nstates_psig
        effcum(j)=0
  100   effcm2(j)=0

      end
c----------------------------------------------------------------------
      subroutine efficiency_prt(passes)
      use contrl, only: idump, irstar, isite, n_conf, nblk, nblkeq, nconf_new, nstep
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'optci.h'
      include 'mstates.h'

      if(iefficiency.eq.0) return

      write(6,*)
      write(6,'(''efficiency for multiple states'')')
      do 200 j=1,nstates_psig
        efficiency=effcum(j)*effcum(j)/effcm2(j)/passes
c       write(6,*) effcum(j)*effcum(j)/passes,effcm2(j)
  200   write(6,'(''efficiency state '',i4,f8.3)') j,efficiency

      end
c----------------------------------------------------------------------
      subroutine efficiency_dump(iu)
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'optci.h'
      include 'mstates.h'

      if(iefficiency.eq.0) return

      write(iu) nstates_psig,(effcum(i),effcm2(i),i=1,nstates_psig)

      end
c-----------------------------------------------------------------------
      subroutine efficiency_rstrt(iu)
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'optci.h'
      include 'mstates.h'

      if(iefficiency.eq.0) return

      read(iu) nstates_psig,(effcum(i),effcm2(i),i=1,nstates_psig)
      if(nstates_psig.ne.nstates) call fatal('EFFICIENCY: different nstates')

      end
c-----------------------------------------------------------------------
