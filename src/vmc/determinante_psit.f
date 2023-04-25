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
      use vmc_mod, only: stoo, nwftypeorb
      implicit none

      integer :: iel, istate, k, iwf_save, o
      real(dp) :: det, determ

      o=stoo(istate)
      iwf_save=iwf
      
      if(nwftypeorb.gt.1) iwf=1
      determ=0.d0
      if(iel.le.nup) then
       do k=1,ndet
          determ=determ+detn(k,o)*detiab(k,2,o)*cdet(k,istate,iwf)
       enddo
      else
         do k=1,ndet
            determ=determ+detn(k,o)*detiab(k,1,o)*cdet(k,istate,iwf)
         enddo
      endif
      iwf=iwf_save

      return
      end
      end module
