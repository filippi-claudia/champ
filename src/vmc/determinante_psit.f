      module determinante_psit_mod
      contains
c-----------------------------------------------------------------------
      subroutine determinante_psit(iel,determ,istate)

      use multiple_geo, only: iwf
      use multislater, only: detiab
      use multislatern, only: detn
      use precision_kinds, only: dp
      use slater, only: ndet, cdet
      use system, only: nup
      implicit none

      integer :: iel, istate, k
      real(dp) :: det, determ

      !STU add istate to orbital mapping here, detn, detiab
      determ=0
      if(iel.le.nup) then
       do k=1,ndet
          determ=determ+detn(k,istate)*detiab(k,2,istate)*cdet(k,istate,iwf)
       enddo
      else
         do k=1,ndet
            determ=determ+detn(k,istate)*detiab(k,1,istate)*cdet(k,istate,iwf)
         enddo
      endif
      

      return
      end
      end module
