      module random_mod
      private
        interface ! xoroshiro c bindings
          function xoroshiro_next() result(res)
     &          bind(c,name="next_uniform_double")
            use iso_c_binding, only: c_double
            real(c_double) res
          end function
          subroutine xoroshiro_set_state(s)
     &               bind(c,name="xoroshiro_set_state")
            use iso_c_binding, only: c_int32_t
            integer(c_int32_t), intent(in) :: s(4)
          end subroutine
          subroutine xoroshiro_get_state(s)
     &               bind(c,name="xoroshiro_get_state")
            use iso_c_binding, only: c_int32_t
            integer(c_int32_t), intent(out) :: s(4)
          end subroutine
        end interface

        interface ! transparent interface for CHAMP
          function next_rn_if() 
            use precision_kinds, only: dp
            real(dp) :: next_rn_if
          end function
        end interface
        procedure(next_rn_if), pointer :: random_dp => NULL()
      public :: random_dp
      public :: setrn, savern
      contains
      subroutine setrn(iseed)

      use contrl_file, only: ounit
      use precision_kinds, only: dp
      use random,  only: ll,switch_rng
      implicit none

      integer       , intent(in) :: iseed(4)

      select case(switch_rng)
        case(0)
          random_dp => rannyu_reference_wrap
          call setrn_rannyu_reference(iseed)
        case default
          random_dp => xoroshiro_wrap
          call xoroshiro_set_state(iseed)
      end select
      endsubroutine

c--------------------------XOROSHIRO----------

      function xoroshiro_wrap()
      use precision_kinds, only: dp
      implicit none
        real(dp)                          :: xoroshiro_wrap
        xoroshiro_wrap = xoroshiro_next()
      end function


c--------------------------RANNYU-----------

      subroutine setrn_rannyu_reference(iseed)
c NYU linear congruential random number generator.
c Uses 48 bits, rather than the usual 32 bits.

      use contrl_file, only: ounit
      use precision_kinds, only: dp
      use random,  only: ll
      implicit none

      integer :: i
      integer, intent(in) :: iseed(4)

      do i=1,4
         ll(i)=iseed(i)
      enddo
      ll(4)=2*(ll(4)/2)+1
      return
      end

      subroutine rannyu_reference(val, idum)
      use random, only: ll, mm
      use precision_kinds, only: dp
      implicit none

      real(dp), intent(inout) :: val
      integer, intent(in) :: idum
      integer :: i1, i2, i3, i4
      integer, parameter :: itwo12 = 4096

      real(dp), parameter :: two12i = 2.44140625d-4


c     On bat it is more efficient to not precompute 2**12 and 2**-12
c     Therefore use original verion of code instead of this one.
      i1=ll(1)*mm(4)+ll(2)*mm(3)+ll(3)*mm(2)+ll(4)*mm(1)
      i2=ll(2)*mm(4)+ll(3)*mm(3)+ll(4)*mm(2)
      i3=ll(3)*mm(4)+ll(4)*mm(3)
      i4=ll(4)*mm(4)
      ll(4)=mod(i4,itwo12)
      i3=i3+i4/itwo12
      ll(3)=mod(i3,itwo12)
      i2=i2+i3/itwo12
      ll(2)=mod(i2,itwo12)
      ll(1)=mod(i1+i2/itwo12,itwo12)
      val=two12i *(dfloat(ll(1))+
     &       two12i*(dfloat(ll(2))+
     &       two12i*(dfloat(ll(3))+
     &       two12i*(dfloat(ll(4))))))
      return
      end

      function rannyu_reference_wrap()
      use precision_kinds, only: dp
      implicit none
        real(dp)                          :: rannyu_reference_wrap
        integer                           :: idum = 1
        call rannyu_reference(rannyu_reference_wrap, idum)
      end function

c-----------------------------------------------------------------------

      subroutine savern(iseed)
      use random, only: ll, switch_rng
      implicit none
      integer iseed(4)

      select case(switch_rng)
        case(0) ! Rannyu
          iseed=ll
        case default ! xoroshiro
          call xoroshiro_get_state(iseed)
      end select
      end

      end module
