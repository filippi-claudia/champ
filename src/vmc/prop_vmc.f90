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
      subroutine prop_sum(wtg,wtg_sqrt)
      use contrl_per, only: iperiodic
      use csfs,    only: nstates
      use mpiconf, only: wid
      use precision_kinds, only: dp
      use prp000,  only: iprop,npropps
      use prp001,  only: vprop
      use prp003,  only: vprop_sum

      implicit none

      integer :: i,i0,istate,j0,jstate,k0
      real(dp), dimension(*) :: wtg, wtg_sqrt

      if(iprop.eq.0) return

      k0=3
      if(iperiodic.eq.1) k0=k0*2

      do istate=1,nstates
        i0=(istate-1)*npropps
        do jstate=1,nstates
          j0=i0+k0*(jstate-1)
          do i=j0+1,j0+k0
            vprop_sum(i)=vprop_sum(i)+wtg_sqrt(istate)*wtg_sqrt(jstate)*vprop(i)
          enddo
        enddo
      enddo

      do istate=1,nstates
        i0=(istate-1)*npropps+6*nstates
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
