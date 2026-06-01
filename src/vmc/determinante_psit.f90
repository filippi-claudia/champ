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
      implicit none

      integer :: iel, istate, iab_other, k, kcdet, o
      real(dp) :: determ

      o=stoo(istate)

      kcdet=iwf
      if(nwftypeorb.gt.1) kcdet=1

      if(iel.le.nup) then
        iab_other=2
      else
        iab_other=1
      endif

      determ=0.d0
      do k=1,ndet
         determ=determ+detn(k,o)*detiab(k,iab_other,o)*cdet(k,istate,kcdet)
      enddo

      return
      end
end module
