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

module c_averages 
  !> Arguments: mprop, prop, wprop, cum_av, cum_av2, cum_w
   use precision_kinds, only: dp

   integer, parameter :: mprop = 100
   real(dp), dimension(:), allocatable :: prop !(mprop)
   real(dp), dimension(:), allocatable :: wprop !(mprop)
   real(dp), dimension(:), allocatable :: cum_av !(mprop)
   real(dp), dimension(:), allocatable :: cum_av2 !(mprop)
   real(dp), dimension(:), allocatable :: cum_w !(mprop)

   private
   public :: mprop, prop, wprop, cum_av, cum_av2, cum_w
   public :: allocate_c_averages, deallocate_c_averages
   save

contains
   subroutine allocate_c_averages()
      if (.not. allocated(prop)) allocate(prop(mprop))
      if (.not. allocated(wprop)) allocate(wprop(mprop))
      if (.not. allocated(cum_av)) allocate(cum_av(mprop))
      if (.not. allocated(cum_av2)) allocate(cum_av2(mprop))
      if (.not. allocated(cum_w)) allocate(cum_w(mprop))
   end subroutine allocate_c_averages

   subroutine deallocate_c_averages()
      if (allocated(prop)) deallocate(prop)
      if (allocated(wprop)) deallocate(wprop)
      if (allocated(cum_av)) deallocate(cum_av)
      if (allocated(cum_av2)) deallocate(cum_av2)
      if (allocated(cum_w)) deallocate(cum_w)
   end subroutine deallocate_c_averages
end module c_averages 

module c_averages_index 
   !> Arguments: jeloc, jderiv 
   use force_mod, only: MFORCE

   integer :: jeloc
   integer, dimension(:,:), allocatable :: jderiv !(3,MFORCE) 

   private
   public :: jeloc, jderiv  
   public :: allocate_c_averages_index, deallocate_c_averages_index 
   save

contains
   subroutine allocate_c_averages_index()
      if (.not. allocated(jderiv)) allocate(jderiv(3, MFORCE))
   end subroutine allocate_c_averages_index

   subroutine deallocate_c_averages_index()
      if (allocated(jderiv)) deallocate(jderiv)
   end subroutine deallocate_c_averages_index
end module c_averages_index 

module jacobsave
   !> Arguments: ajacob, ajacold
   use dmc_mod, only: MWALK
   use force_mod, only: MFORCE
   use precision_kinds, only: dp

    real(dp) :: ajacob
    real(dp), dimension(:,:), allocatable :: ajacold

    private
    public :: ajacob, ajacold
    public :: allocate_jacobsave, deallocate_jacobsave
    save

contains
   subroutine allocate_jacobsave()
      if (.not. allocated(ajacold)) allocate(ajacold(MWALK, MFORCE))
   end subroutine allocate_jacobsave

   subroutine deallocate_jacobsave()
      if (allocated(ajacold)) deallocate(ajacold)
   end subroutine deallocate_jacobsave

end module jacobsave
