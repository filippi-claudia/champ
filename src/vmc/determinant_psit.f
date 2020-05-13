      subroutine determinant_psit(determ,istate)

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use dets, only: cdet, ndet
      use elec, only: ndn, nup
      use wfsec, only: iwf, iwftype, nwftype
      use contrl_per, only: iperiodic, ibasis

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      parameter (one=1.d0,half=0.5d0)

      common /multislater/ detu(MDET),detd(MDET)

      determ=0

      do 110 k=1,ndet
  110   determ=determ+detu(k)*detd(k)*cdet(k,istate,iwf)


      return
      end
