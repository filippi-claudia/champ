c-----------------------------------------------------------------------
      subroutine prop_save
      use prp000, only: iprop, ipropprt, nprop
      use prp001, only: vprop
      implicit real*8(a-h,o-z)


      include 'properties.h'
      include 'prop_vmc.h'

      if(iprop.eq.0) return
      do i=1,nprop
       vprop_old(i)=vprop(i)
      enddo
      end
c-----------------------------------------------------------------------
      subroutine prop_sum(p,q)
      use prp000, only: iprop, ipropprt, nprop
      use prp001, only: vprop
      implicit real*8(a-h,o-z)


      include 'properties.h'
      include 'prop_vmc.h'

      if(iprop.eq.0) return
      do i=1,nprop
       vprop_sum(i)=vprop_sum(i)+p*vprop(i)+q*vprop_old(i)
      enddo
      end
