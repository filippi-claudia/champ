      subroutine determinant_psit(determ,istate)

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use dets, only: cdet, ndet
      use elec, only: ndn, nup
      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      parameter (one=1.d0,half=0.5d0)

      common /contrl_per/ iperiodic,ibasis

      common /multislater/ detu(MDET),detd(MDET)

      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      determ=0

      do 110 k=1,ndet
  110   determ=determ+detu(k)*detd(k)*cdet(k,istate,iwf)


      return
      end
