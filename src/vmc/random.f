      module random_mod
      private
        interface ! xoshiro c bindings
          function xoshiro_next() result(res)
     &          bind(c,name="next_uniform_double")
            use iso_c_binding, only: c_double
            real(c_double) res
          end function
          subroutine xoshiro_seed(s) bind(c,name="xoshiro_seed")
            use iso_c_binding, only: c_int64_t
            integer(c_int64_t), intent(in) :: s(4)
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
      subroutine setrn(legaseed)

      use iso_c_binding, only: c_int64_t
      use contrl_file, only: ounit
      use precision_kinds, only: dp
      use rnyucm,  only: ll,switch_rng
      implicit none

      integer(c_int64_t)         :: iseed(4)
      integer       , intent(in) :: legaseed(4)

      iseed(1) = splitmix64(legaseed(1))
      iseed(2) = splitmix64(legaseed(2))
      iseed(3) = splitmix64(legaseed(3))
      iseed(4) = splitmix64(legaseed(4))

      select case(switch_rng)
        case(0)
          random_dp => rannyu_reference_wrap
          call setrn_rannyu_reference(legaseed)
        case(-1)
          random_dp => std_reference_wrap
          call setrn_rannyu_std(legaseed)
        case default
          random_dp => xoshiro_wrap
          call xoshiro_seed(iseed)
      end select
      endsubroutine

c--------------------------XOSHIRO----------

      function xoshiro_wrap()
      use precision_kinds, only: dp
      implicit none
        real(dp)                          :: xoshiro_wrap
        xoshiro_wrap = xoshiro_next()
      end function


c--------------------------RANNYU-----------

      subroutine setrn_rannyu_reference(iseed)
c NYU linear congruential random number generator.
c Uses 48 bits, rather than the usual 32 bits.

      use contrl_file, only: ounit
      use precision_kinds, only: dp
      use rnyucm,  only: ll
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
      use rnyucm, only: ll, mm
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

c------------------------------COMPILER SPECIFIC----------------
      subroutine setrn_rannyu_std(iseed)
      use contrl_file, only: ounit
      use precision_kinds, only: dp
      implicit none

      integer, intent(in) :: iseed(4)
      integer size
      integer, allocatable :: new_seed(:)
      integer i

      ! determine the seed size
      call random_seed(size = size)

      ! allocate the new seed
      allocate(new_seed(size))

      ! input the seed
      do i=1,min(size,4)
            new_seed(i) = iseed(i)
      end do

      ! check the new seed
      write(ounit, *) 'seed size', size
      do i=1,size
            write(ounit, *) 'seed     ', i, new_seed(i)
      end do

      ! set the new seed
      call random_seed(put=new_seed)

      ! deallocate
      deallocate(new_seed)
      return
      end

      function std_reference_wrap()
      use precision_kinds, only: dp
      implicit none
        real(dp)                          :: std_reference_wrap
        call random_number(std_reference_wrap)
      end function

c-----------------------------------------------------------------------

      subroutine savern(iseed)
      use rnyucm, only: ll
      implicit none

      integer :: i


      integer iseed(4)
      do i=1,4
         iseed(i)=ll(i)
      enddo
      return
      end

c------------------------------------SPLITMIX HELPER
c helps with generating a quality seed (64bit) with only 32bit of input
c needed because of legacy reasons the seed is 4x32=128bit, but
c xoshiro needs 256bits of seed (preferably not half being zero) which happens with a simple cast
c Victor Azizi
      function splitmix64(seed) result(res)
        use iso_fortran_env, only: int64
        integer       , intent(in) :: seed
        integer(int64)             :: res

        integer(INT64),parameter::cnst_step1= -7046029254386353131_int64
        integer(INT64),parameter::cnst_step2= -4658895280553007687_int64 
        integer(INT64),parameter::cnst_step3= -7723592293110705685_int64
        integer(INT64) :: state

        state   = seed + cnst_step1
        state   = ieor_shiftr( i= state , shift= 30 ) * cnst_step2
        state   = ieor_shiftr( i= state , shift= 27 ) * cnst_step3
        res     = ieor_shiftr( i= state , shift= 31 )
      contains
      pure elemental function ieor_shiftr ( i , shift )

          ! argument(s) for this <function>
          integer (INT64) , intent(in) :: i
          integer         , intent(in) :: shift

          ! return value of this <function>
          integer(INT64) :: ieor_shiftr

          ieor_shiftr = ieor( i= i, j= shiftr( i= i , shift= shift ) )

      end function ieor_shiftr
      end function
      end module
