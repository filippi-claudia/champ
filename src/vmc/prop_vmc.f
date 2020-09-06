c-----------------------------------------------------------------------
      subroutine prop_save
      use prp000, only: iprop, ipropprt, nprop
      use prp001, only: vprop
      use prp002, only: vprop_old2

      implicit real*8(a-h,o-z)

      include 'properties.h'

      if(iprop.eq.0) return
      do i=1,nprop
       vprop_old2(i)=vprop(i)
      enddo
      end
c-----------------------------------------------------------------------
      subroutine prop_sum(p,q)
      use prp000, only: iprop, ipropprt, nprop
      use prp001, only: vprop
      use prp002, only: vprop_old2

      use prp003, only: cc_nuc, vprop_cm2, vprop_cum, vprop_sum

      implicit real*8(a-h,o-z)


      include 'properties.h'

      if(iprop.eq.0) return
      do i=1,nprop
       vprop_sum(i)=vprop_sum(i)+p*vprop(i)+q*vprop_old2(i)
      enddo
      end
