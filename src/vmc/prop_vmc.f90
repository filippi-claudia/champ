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
      subroutine prop_sum(wtg)
      use precision_kinds, only: dp
      use mstates_mod, only: MSTATES
      use csfs,    only: nstates
      use prp000,  only: iprop,npropps
      use prp001,  only: vprop
      use prp003,  only: vprop_sum
      use mpiconf, only: wid

      implicit none

      integer :: i,i0,istate
      real(dp), dimension(MSTATES) :: wtg

      if(iprop.eq.0) return

      do istate=1,nstates
       i0=(istate-1)*npropps
       do i=i0+1,i0+npropps
        vprop_sum(i)=vprop_sum(i)+wtg(istate)*vprop(i)
       enddo
      enddo
      end
!-----------------------------------------------------------------------
      subroutine prop_sum_moveall(p,q)
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
