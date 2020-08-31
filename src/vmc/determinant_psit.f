      subroutine determinant_psit(determ,istate)

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use dets, only: cdet, ndet
      use wfsec, only: iwf

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      common /multislater/ detiab(MDET,2)

      determ=0
      do 110 k=1,ndet
  110   determ=determ+detiab(k,1)*detiab(k,2)*cdet(k,istate,iwf)


      return
      end
