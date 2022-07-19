      module determinante_psit_mod
      contains
c-----------------------------------------------------------------------
      subroutine determinante_psit(iel,determ,istate)

      use dets, only: cdet, ndet
      use wfsec, only: iwf
      use multislatern, only: detn

      use multislater, only: detiab
      use precision_kinds, only: dp
      use system, only: nup
      implicit none

      integer :: iel, istate, k
      real(dp) :: det, determ






      determ=0

      if(iel.le.nup) then
       do k=1,ndet
          determ=determ+detn(k)*detiab(k,2)*cdet(k,istate,iwf)
       enddo
      else
         do k=1,ndet
            determ=determ+detn(k)*detiab(k,1)*cdet(k,istate,iwf)
         enddo
      endif
      

      return
      end
      end module
