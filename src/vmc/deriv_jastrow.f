      module deriv_jastrow_mod
      contains
      subroutine deriv_jastrow(x,v,d2,div_vj,value)
c Written by Claudia Filippi

      use system, only: nelec

      use precision_kinds, only: dp

      use deriv_jastrow4_mod, only: deriv_jastrow4
      implicit none

      integer :: i
      real(dp) :: d2, value
      real(dp), dimension(3, *) :: x
      real(dp), dimension(3, *) :: v
      real(dp), dimension(nelec) :: div_vj
      real(dp), parameter :: zero = 0.d0


      do i=1,nelec
        v(1,i)=zero
        v(2,i)=zero
        v(3,i)=zero
      enddo
      d2=zero

      call deriv_jastrow4(x,v,d2,value)

      return
      end
      end module
