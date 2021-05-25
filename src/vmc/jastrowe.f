      subroutine jastrowe(iel,x,v,d2,value,iflag)
c Written by Claudia Filippi by modifying jastrow

      use const, only: nelec
      use precision_kinds, only: dp

      implicit none

      integer :: i, iel, iflag
      real(dp) :: d2, value
      real(dp), dimension(3, *) :: x
      real(dp), dimension(3, *) :: v
      real(dp), parameter :: zero = 0.d0


      do 10 i=1,nelec
        v(1,i)=zero
        v(2,i)=zero
   10   v(3,i)=zero

      call jastrow4e(iel,x,v,d2,value,iflag)

      return
      end
