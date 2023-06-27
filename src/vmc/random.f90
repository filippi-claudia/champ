module random_mod
      private
        interface ! xoroshiro c bindings
          function xoroshiro_next() result(res) &
                bind(c,name="next_uniform_double")
            use iso_c_binding, only: c_double
            real(c_double) res
          end function
          subroutine xoroshiro_set_state(s) &
                     bind(c,name="xoroshiro_set_state")
            use iso_c_binding, only: c_int32_t
            integer(c_int32_t), intent(in) :: s(8)
          end subroutine
          subroutine xoroshiro_get_state(s) &
                     bind(c,name="xoroshiro_get_state")
            use iso_c_binding, only: c_int32_t
            integer(c_int32_t), intent(out) :: s(8)
          end subroutine
          subroutine jumprn() &
                     bind(c,name="jump")
          end subroutine
        end interface

      public :: random_dp
      public :: setrn, savern, jumprn
contains
      subroutine setrn(iseed)

      use contrl_file, only: ounit
      use precision_kinds, only: dp
      use random,  only: ll,switch_rng
      implicit none

      integer       , intent(in) :: iseed(8) ! 256-bits (8x32bits integer)

        call xoroshiro_set_state(iseed)

      endsubroutine

      function random_dp()
      use precision_kinds, only: dp
      implicit none
        real(dp)                          :: random_dp
        random_dp = xoroshiro_next()
      end function

      subroutine savern(iseed)
      use random, only: ll, switch_rng
      implicit none
      integer iseed(8)

        call xoroshiro_get_state(iseed)
      end

end module
