module determinant_psit_mod
contains
      subroutine determinant_psit(determ,istate)

      use slater, only: ndet, cdet
      use multiple_geo, only: iwf
      use multislater, only: detiab
      use precision_kinds, only: dp
      use vmc_mod, only: stoo, nwftypeorb
      implicit none

      integer :: istate, k, iwf_save, o
      real(dp) :: determ

      determ=0.0d0
      iwf_save=iwf
      o=stoo(istate)
      if(nwftypeorb.gt.1) iwf=1
      do k=1,ndet
        determ=determ+detiab(k,1,o)*detiab(k,2,o)*cdet(k,istate,iwf)
      enddo
      iwf=iwf_save


      return
      end
end module
