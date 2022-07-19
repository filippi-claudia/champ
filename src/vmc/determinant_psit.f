      module determinant_psit_mod
      contains
      subroutine determinant_psit(determ,istate)

      use multiple_geo, only: iwf
      use multislater, only: detiab
      use precision_kinds, only: dp
      use slater, only: ndet
      use slater, only: cdet
      implicit none

      integer :: istate, k
      real(dp) :: determ


      determ=0.0d0
      do k=1,ndet
        determ=determ+detiab(k,1)*detiab(k,2)*cdet(k,istate,iwf)
      enddo


      return
      end
      end module
