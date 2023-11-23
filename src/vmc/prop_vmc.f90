module prop_vmc
contains
!-----------------------------------------------------------------------
      subroutine prop_save
      use prp000,  only: iprop,nprop
      use prp001,  only: vprop
      use prp002,  only: vprop_old2

      implicit none

      integer :: i



      if(iprop.eq.0) return
      do i=1,nprop
       vprop_old2(i)=vprop(i)
      enddo
      end
!-----------------------------------------------------------------------
      subroutine prop_sum(p,q)
      use precision_kinds, only: dp
      use prp000,  only: iprop,nprop
      use prp001,  only: vprop
      use prp002,  only: vprop_old2
      use prp003,  only: vprop_sum


      implicit none

      integer :: i
      real(dp) :: p, q



      if(iprop.eq.0) return
      do i=1,nprop
       vprop_sum(i)=vprop_sum(i)+p*vprop(i)+q*vprop_old2(i)
      enddo
      end
end module
