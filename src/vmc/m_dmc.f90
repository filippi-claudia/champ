module dmc_mod
    !> Arguments: MWALK, MFPROD, MFPRD1, MPATH
    !> MWALK: Maximum number of walkers
    integer, parameter :: MWALK = 600
    integer, parameter :: MFPROD=3201
    integer, parameter :: MFPRD1=MFPROD-1
    integer, parameter :: MPATH=999

    private
    public :: MWALK, MFPROD, MFPRD1, MPATH
    save
end module  dmc_mod

module age
  !> Arguments: iage, ioldest, ioldestmx
  use precision_kinds, only: dp
  use dmc_mod, only: MWALK 

  integer, dimension(:), allocatable:: iage
  integer  :: ioldest
  integer  :: ioldestmx

  private
  public :: iage, ioldest, ioldestmx
  public :: allocate_iage, deallocate_iage
  save
contains

  subroutine allocate_iage()
      if (.not. allocated(iage)) allocate (iage(MWALK))
  end subroutine allocate_iage

  subroutine deallocate_iage()
      if (allocated(iage)) deallocate (iage)
  end subroutine deallocate_iage
end module age

module iterat
  !> Arguments: iblk, ipass

   integer  :: iblk
   integer  :: ipass

   private
   public :: iblk, ipass
   save
end module iterat
