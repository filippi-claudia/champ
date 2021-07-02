c-----------------------------------------------------------------------
      subroutine determinante_psit(iel,determ,istate)

      use dets, only: cdet, ndet
      use elec, only: nup
      use wfsec, only: iwf
      use multislatern, only: detn

      use multislater, only: detiab
      use precision_kinds, only: dp
      implicit none

      integer :: iel, istate, k
      real(dp) :: det, determ






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
