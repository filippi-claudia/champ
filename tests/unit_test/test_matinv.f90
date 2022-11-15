module test_matinv_mod
  use fortutf
  use precision_kinds, only: dp
  implicit none
  contains
subroutine test_matinv
  use system, only: nelec
  use matinv_mod, only: matinv

  real(dp), dimension(3,3) :: a       = reshape([ 1, 0, 0,  &
                                                  0, 1, 0,  &
                                                  0, 0, 1], [3,3])
  real(dp), dimension(  9) :: a_expected      = [ 1, 0, 0   &
                                                , 0, 1, 0   &
                                                , 0, 0, 1]
  real(dp)                 :: a_expected_det = 1.

  real(dp), dimension(3,3) :: b       = reshape([ 1, 0, 0   &
                                                , 5, 1, 0   &
                                                , 0, 0, 1], [3,3])
  real(dp), dimension(  9) :: b_expected      = [ 1, 0, 0   &
                                                ,-5, 1, 0   &
                                                , 0, 0, 1]
  real(dp)                 :: b_expected_det  = 1.

  real(dp)                 :: det             = -1.


  nelec = 6 !TODO matinv could be independent of nelec

  call tag_test("test matinv_a")
  call matinv(a,3,det)
  call assert_equal(reshape(a,[9]),a_expected)
  call assert_equal(det, a_expected_det)

  call tag_test("test matinv_b")
  call matinv(b,3,det)
  call assert_equal(reshape(b,[9]),b_expected)
  call assert_equal(det, b_expected_det)
end subroutine test_matinv

end module
