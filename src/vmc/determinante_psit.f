c-----------------------------------------------------------------------
      subroutine determinante_psit(iel,determ,istate)

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /elec/ nup,ndn

      common /multislater/ detu(MDET),detd(MDET)
      common /multislatern/ detn(MDET)
     &,orb(MORB),dorb(3,MORB),ddorb(MORB)

      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      determ=0

      do 110 k=1,ndet
        if(iel.le.nup) then
          det=detn(k)*detd(k)
         else
          det=detu(k)*detn(k)
        endif

        determ=determ+det*cdet(k,istate,iwf)
  110 continue

      return
      end
