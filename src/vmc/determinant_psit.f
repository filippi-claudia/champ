      subroutine determinant_psit(determ,istate)

      use dets, only: cdet, ndet
      use wfsec, only: iwf

      use multislater, only: detiab
      use precision_kinds, only: dp
      implicit none

      integer :: istate, k
      real(dp) :: determ




      determ=0
      do 110 k=1,ndet
  110   determ=determ+detiab(k,1)*detiab(k,2)*cdet(k,istate,iwf)


      return
      end
