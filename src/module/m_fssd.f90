module fssd
  use precision_kinds, only: dp
  use atom, only: ncent_tot


  implicit none

  integer :: ifssd, norm_fssd, div_fssd, iset
  real(dp) :: alfa_dfssd, errbar_fssd
  real(dp), dimension(:, :), allocatable :: d_fssd

  private
  public   :: d_fssd, ifssd, errbar_fssd, iset
  public   :: norm_fssd, div_fssd, alfa_dfssd
  public   :: allocate_fssd, deallocate_fssd

  save
contains
  subroutine allocate_fssd()
    if (.not. allocated(d_fssd)) allocate (d_fssd(3, ncent_tot))
  end subroutine allocate_fssd

  subroutine deallocate_fssd()
    if (allocated(d_fssd)) deallocate (d_fssd)
  end subroutine deallocate_fssd

end module fssd
