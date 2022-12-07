module dmc_mod
    !> Arguments: mwalk, MFPROD, MFPRD1
    !> mwalk: Maximum number of walkers

    implicit none

    integer            :: mwalk
    integer, parameter :: MFPROD=200
    integer, parameter :: MFPRD1=MFPROD-1

    private
    public :: mwalk, set_mwalk, MFPROD, MFPRD1
    save
    contains
    subroutine set_mwalk
      use control_dmc, only: dmc_nconf
      mwalk = dmc_nconf + int(0.5*dmc_nconf)  ! maximum number of walkers are 50% more than nwalk
    end subroutine set_mwalk
end module  dmc_mod

module age
  !> Arguments: iage, ioldest, ioldestmx
      use dmc_mod, only: mwalk

  implicit none

  integer, dimension(:), allocatable:: iage
  integer  :: ioldest
  integer  :: ioldestmx

  private
  public :: iage, ioldest, ioldestmx
  public :: allocate_iage, deallocate_iage
  save
contains

  subroutine allocate_iage()
      if (.not. allocated(iage)) allocate (iage(mwalk), source=0)
  end subroutine allocate_iage

  subroutine deallocate_iage()
      if (allocated(iage)) deallocate (iage)
  end subroutine deallocate_iage
end module age

module branch
  !> Arguments: eest, eigv, eold, ff, fprod, nwalk, pwt, wdsumo, wgdsumo, wt, wtgen,
  !> wthist
      use dmc_mod, only: MFPRD1,mwalk
      use multiple_geo, only: MFORCE,MFORCE_WT_PRD
      use precision_kinds, only: dp

  implicit none

   real(dp) :: eest
   real(dp) :: esigma
   real(dp) :: eigv
   real(dp), dimension(:,:), allocatable :: eold !(mwalk,MFORCE)
   real(dp), dimension(:), allocatable :: ff !(MFPRD1)
   real(dp) :: fprod
   integer  :: nwalk
   real(dp), dimension(:,:), allocatable:: pwt !(mwalk,MFORCE)
   real(dp) :: wdsumo
   real(dp) :: wgdsumo
   real(dp), dimension(:), allocatable :: wt !(mwalk)
   real(dp), dimension(:), allocatable :: wtgen !(MFPRD1)
   real(dp), dimension(:,:,:), allocatable :: wthist!(mwalk,MFORCE_WT_PRD,MFORCE)

   private
   public :: eest, esigma, eigv, eold, ff, fprod, nwalk, pwt, wdsumo, wgdsumo, wt
   public :: wtgen, wthist
   public :: allocate_branch, deallocate_branch
   save

contains
   subroutine allocate_branch()
      if (.not. allocated(eold)) allocate(eold(mwalk,MFORCE))
      if (.not. allocated(ff)) allocate(ff(0:MFPRD1))
      if (.not. allocated(pwt)) allocate(pwt(mwalk,MFORCE))
      if (.not. allocated(wt)) allocate(wt(mwalk))
      if (.not. allocated(wtgen)) allocate(wtgen(0:MFPRD1))
      if (.not. allocated(wthist)) allocate(wthist(mwalk,0:MFORCE_WT_PRD,MFORCE))
   end subroutine allocate_branch

!splitj.f:      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(mwalk,MFORCE),
!splitj.f:     &pwt(mwalk,MFORCE),wthist(mwalk,0:MFORCE_WT_PRD,MFORCE),
!splitj.f-     &wt(mwalk),eigv,eest,wdsumo,wgdsumo,fprod,nwalk
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

   implicit none

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
      use multiple_geo, only: MFORCE

   implicit none

   integer :: jeloc
   integer, dimension(:,:), allocatable :: jderiv !(3,MFORCE)

   private
   public :: jeloc, jderiv
   public :: allocate_c_averages_index, deallocate_c_averages_index
   save

contains
   subroutine allocate_c_averages_index()
      if (.not. allocated(jderiv)) allocate(jderiv(3, MFORCE), source=0)
   end subroutine allocate_c_averages_index

   subroutine deallocate_c_averages_index()
      if (allocated(jderiv)) deallocate(jderiv)
   end subroutine deallocate_c_averages_index
end module c_averages_index

module jacobsave
   !> Arguments: ajacob, ajacold
      use dmc_mod, only: mwalk
      use multiple_geo, only: MFORCE
      use precision_kinds, only: dp

   implicit none

    real(dp) :: ajacob
    real(dp), dimension(:,:), allocatable :: ajacold

    private
    public :: ajacob, ajacold
    public :: allocate_jacobsave, deallocate_jacobsave
    save

contains
   subroutine allocate_jacobsave()
      if (.not. allocated(ajacold)) allocate(ajacold(mwalk, MFORCE))
   end subroutine allocate_jacobsave

   subroutine deallocate_jacobsave()
      if (allocated(ajacold)) deallocate(ajacold)
   end subroutine deallocate_jacobsave

end module jacobsave

module velratio
   !> Arguments: fratio, xdrifted
      use dmc_mod, only: mwalk
      use multiple_geo, only: MFORCE
      use precision_kinds, only: dp

   implicit none

   real(dp), dimension(:,:), allocatable :: fratio !(mwalk,MFORCE)
   real(dp), dimension(:,:,:,:), allocatable :: xdrifted !(3,MELEC,mwalk,MFORCE)

   private
   public :: fratio, xdrifted
   public :: allocate_velratio, deallocate_velratio
   save

contains
   subroutine allocate_velratio()
      use dmc_mod, only: mwalk
      use multiple_geo, only: MFORCE
      use system,  only: nelec
      if (.not. allocated(fratio)) allocate(fratio(mwalk, MFORCE))
      if (.not. allocated(xdrifted)) allocate(xdrifted(3, nelec, mwalk, MFORCE))
   end subroutine allocate_velratio

   subroutine deallocate_velratio()
      if (allocated(fratio)) deallocate(fratio)
      if (allocated(xdrifted)) deallocate(xdrifted)
   end subroutine deallocate_velratio
 end module velratio
