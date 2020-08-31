c-----------------------------------------------------------------------
      subroutine determinante_psit(iel,determ,istate)

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use dets, only: cdet, ndet
      use elec, only: ndn, nup
      use wfsec, only: iwf, iwftype, nwftype
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      common /multislater/ detiab(MDET,2)
      common /multislatern/ detn(MDET)
     &,orbn(MORB),dorbn(3,MORB),ddorbn(MORB)


      determ=0

      do 110 k=1,ndet
        if(iel.le.nup) then
          det=detn(k)*detiab(k,2)
         else
          det=detiab(k,1)*detn(k)
        endif

        determ=determ+det*cdet(k,istate,iwf)
  110 continue

      return
      end
