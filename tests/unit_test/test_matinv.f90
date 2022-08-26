module test_matinv_mod
  use fortutf
  use precision_kinds, only: dp
  implicit none
  contains
subroutine test_matinv_id

    use system, only: nelec
    use matinv_mod, only: matinv
    real(dp), dimension(3,3) :: a               = reshape([1,0,0,0,1,0,0,0,1],[3,3])
    real(dp), dimension(3,3) :: expected_output = reshape([1,0,0,0,1,0,0,0,1],[3,3])
    real(dp)                 :: determinant     = -1.
    real(dp)                 :: expected_determinant = 1.

    call tag_test("test matinv")

    nelec = 3 ! TODO is this really necessary?

    call matinv(a,3,determinant)

    call assert_equal(reshape(a,[9]),reshape(expected_output,[9]))
    call assert_equal(determinant,expected_determinant)

end subroutine test_matinv_id

subroutine test_matinv_2

    use system, only: nelec
    use matinv_mod, only: matinv
    real(dp), dimension(3,3) :: a               = reshape([1,0,0,5,1,0,0,0,1],[3,3])
    real(dp), dimension(3,3) :: expected_output = reshape([1,0,0,-5,1,0,0,0,1],[3,3])
    real(dp)                 :: determinant     = -1
    real(dp)                 :: expected_determinant = 1.

    call tag_test("test matinv")

    nelec = 3 ! TODO is this really necessary?

    call matinv(a,3,determinant)

    call assert_equal(reshape(a,[9]),reshape(expected_output,[9]))
    call assert_equal(determinant,expected_determinant)

end subroutine test_matinv_2

end module
