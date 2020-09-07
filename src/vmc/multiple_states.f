c----------------------------------------------------------------------
      subroutine efficiency_sample(ipass,determ_s,determ_psig)

      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use mstates_ctrl, only: iefficiency, nstates_psig
      use mstates2, only: effcm2, effcum
      implicit real*8(a-h,o-z)


      include 'optci.h'

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
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use mstates_ctrl, only: nstates_psig
      use mstates2, only: effcm2, effcum
      implicit real*8(a-h,o-z)


      include 'optci.h'
      do 100 j=1,nstates_psig
        effcum(j)=0
  100   effcm2(j)=0

      end
c----------------------------------------------------------------------
      subroutine efficiency_prt(passes)
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use mstates_ctrl, only: iefficiency, nstates_psig
      use mstates2, only: effcm2, effcum
      implicit real*8(a-h,o-z)



      include 'optci.h'

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
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use mstates_ctrl, only: iefficiency, nstates_psig
      use mstates2, only: effcm2, effcum
      implicit real*8(a-h,o-z)



      include 'optci.h'

      if(iefficiency.eq.0) return

      write(iu) nstates_psig,(effcum(i),effcm2(i),i=1,nstates_psig)

      end
c-----------------------------------------------------------------------
      subroutine efficiency_rstrt(iu)
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use mstates_ctrl, only: iefficiency, nstates_psig
      use mstates2, only: effcm2, effcum
      implicit real*8(a-h,o-z)



      include 'optci.h'

      if(iefficiency.eq.0) return

      read(iu) nstates_psig,(effcum(i),effcm2(i),i=1,nstates_psig)
      if(nstates_psig.ne.nstates) call fatal('EFFICIENCY: different nstates')

      end
c-----------------------------------------------------------------------
