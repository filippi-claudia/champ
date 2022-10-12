module test_pe_mod
  use fortutf
  use precision_kinds, only: dp
  implicit none
  contains
subroutine test_pe
  use main_mod, only: initialize_main, finalize_main
  use contrl_file, only: file_input

  implicit none
  integer                  :: ifr            = 1.
  integer                  :: stat


  call initialize_main()

  call tag_test("test pe")

  call finalize_main()

end subroutine test_pe

end module test_pe_mod
