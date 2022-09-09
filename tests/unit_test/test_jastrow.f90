module test_jastrow_mod
  use fortutf
  use precision_kinds, only: dp
  implicit none
  contains
subroutine test_jastrow
  use main_mod, only: initialize_main, finalize_main
  use jastrow_mod, only: jastrow
  use contrl_file, only: file_input

  real(dp), dimension(3,3) :: a       = reshape([ 1, 0, 0,  &
                                                  0, 1, 0,  &
                                                  0, 0, 1], [3,3])
  real(dp), dimension(  9) :: a_expected      = [ 1, 0, 0   &
                                                , 0, 1, 0   &
                                                , 0, 0, 1]
  real(dp)                 :: a_expected_det = 1.

  real(dp), dimension(3,3) :: v
  real(dp)                 :: value, d2        
  integer                  :: ifr            = 1.
  integer                  :: stat
  
  file_input = './jastrow/jastrow.inp'

  call initialize_main()

  call tag_test("test jastrow_a")

  call jastrow(a,v,d2,value, ifr)

  call finalize_main()

end subroutine test_jastrow

end module
