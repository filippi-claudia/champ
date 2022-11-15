module test_random_mod
  use fortutf
  use precision_kinds, only: dp
  implicit none
  contains
subroutine test_random
  use random_mod, only: setrn, random_dp
  use rnyucm, only: switch_rng

  integer :: n                     

  call setrn([1,2,3,4])

  call tag_test("test if random is between 0 and 1, xoshiro")
  do n = 1,100
    if (random_dp() > 1. .or. random_dp() < 0.) then
      call FAIL()
    end if
  end do
  call tag_test("test if random is between 0 and 1, legacy")
  switch_rng = 0
  call setrn([1,2,3,4])
  do n = 1,100
    if (random_dp() > 1. .or. random_dp() < 0.) then
      call FAIL()
    end if
  end do

  !call assert_equal(reshape(a,[9]),a_expected)

end subroutine test_random

end module
