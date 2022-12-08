      module determinant_psit_mod
      contains
      subroutine determinant_psit(determ,istate)

      use slater, only: ndet, cdet
      use multiple_geo, only: iwf

      use multislater, only: detiab
      use precision_kinds, only: dp
      implicit none

      integer :: istate, k
      real(dp) :: determ

      !STU use state orb mapping here
      determ=0.0d0
      do k=1,ndet
        determ=determ+detiab(k,1,istate)*detiab(k,2,istate)*cdet(k,istate,iwf)
      enddo


      return
      end
      end module
