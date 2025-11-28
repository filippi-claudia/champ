module determinante_psit_mod
contains
!-----------------------------------------------------------------------
      subroutine determinante_psit(iel,determ,istate)

      use multiple_geo, only: iwf
      use multislater, only: detiab
      use multislatern, only: detn
      use precision_kinds, only: dp
      use slater, only: ndet, cdet
      use system, only: nup
      use vmc_mod, only: stoo, nwftypeorb
      use m_backflow, only: ibackflow, detn_bf
      implicit none

      integer :: iel, istate, k, iwf_save, o
      real(dp) :: det, determ

      o=stoo(istate)
      iwf_save=iwf

      if(nwftypeorb.gt.1) iwf=1
      determ=0.d0
      
      if (ibackflow.gt.0) then
         determ = detn_bf(1,iwf) * detn_bf(2,iwf)*cdet(1,istate,iwf)
      else
         if(iel.le.nup) then
         do k=1,ndet
            determ=determ+detn(k,o)*detiab(k,2,o)*cdet(k,istate,iwf)
         enddo
         else
            do k=1,ndet
               determ=determ+detn(k,o)*detiab(k,1,o)*cdet(k,istate,iwf)
            enddo
         endif
      endif
      iwf=iwf_save

      return
      end
end module
