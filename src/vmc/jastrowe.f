      subroutine jastrowe(iel,x,v,d2,value,iflag,istate)
c Written by Claudia Filippi by modifying jastrow
      use const, only: nelec
      implicit real*8(a-h,o-z)

      dimension x(3,*),v(3,*)

      do i=1,nelec
         v(1,i)=0.0d0
         v(2,i)=0.0d0
         v(3,i)=0.0d0
      enddo

      call jastrow4e(iel,x,v,d2,value,iflag,istate)

      end subroutine
