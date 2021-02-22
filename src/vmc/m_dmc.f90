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

module branch
  !> Arguments: eest, eigv, eold, ff, fprod, nwalk, pwt, wdsumo, wgdsumo, wt, wtgen,
  !> wthist
   use precision_kinds, only: dp
   use force_mod, only: MFORCE, MFORCE_WT_PRD
   use dmc_mod, only: MWALK, MFPRD1

   real(dp) :: eest
   real(dp) :: eigv
   real(dp), dimension(:,:), allocatable :: eold !(MWALK,MFORCE)
   real(dp), dimension(:), allocatable :: ff !(MFPRD1)
   real(dp) :: fprod
   integer  :: nwalk
   real(dp), dimension(:,:), allocatable:: pwt !(MWALK,MFORCE)
   real(dp) :: wdsumo
   real(dp) :: wgdsumo
   real(dp), dimension(:), allocatable :: wt !(MWALK)
   real(dp), dimension(:), allocatable :: wtgen !(MFPRD1)
   real(dp), dimension(:,:,:), allocatable :: wthist!(MWALK,MFORCE_WT_PRD,MFORCE)

   private
   public :: eest, eigv, eold, ff, fprod, nwalk, pwt, wdsumo, wgdsumo, wt
   public :: wtgen, wthist
   public :: allocate_branch, deallocate_branch
   save

contains
   subroutine allocate_branch()
      if (.not. allocated(eold)) allocate(eold(MWALK,MFORCE))
      if (.not. allocated(ff)) allocate(ff(MFPRD1))
      if (.not. allocated(pwt)) allocate(pwt(MWALK,MFORCE))
      if (.not. allocated(wt)) allocate(wt(MWALK))
      if (.not. allocated(wtgen)) allocate(wtgen(MFPRD1))
      if (.not. allocated(wthist)) allocate(wthist(MWALK,MFORCE_WT_PRD,MFORCE))
   end subroutine allocate_branch

   subroutine deallocate_branch()
      if (allocated(eold)) deallocate(eold)
      if (allocated(ff)) deallocate(ff)
      if (allocated(pwt)) deallocate(pwt)
      if (allocated(wt)) deallocate(wt)
      if (allocated(wtgen)) deallocate(wtgen)
      if (allocated(wthist)) deallocate(wthist)
   end subroutine deallocate_branch
end module branch
