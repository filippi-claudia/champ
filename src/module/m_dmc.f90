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
   use dmc_mod, only: mwalk, MFPRD1
   use multiple_geo, only: MFORCE, MFORCE_WT_PRD
   use precision_kinds, only: dp
   use system, only: nelec, ncent

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

module fragments 
   use dmc_mod, only: mwalk
   use multiple_geo, only: MFORCE
   use precision_kinds, only: dp
   use system, only: nelec, ncent

   implicit none
   ! electron fragments
   real(dp), dimension(:), allocatable :: eloc_i !(nelec)
   real(dp), dimension(:,:,:), allocatable :: eloco_i !(nelec,mwalk,MFORCE)
   real(dp), dimension(:), allocatable :: eest_i !(nelec)
   real(dp), dimension(:), allocatable :: esum_i !(nelec)
   real(dp), dimension(:), allocatable :: ecum_i !(nelec)
   real(dp), dimension(:), allocatable :: v2_i !(nelec)
   real(dp), dimension(:), allocatable :: vav2_i !(nelec)
   real(dp), dimension(:,:,:), allocatable :: fratio_i !(nelec,mwalk,MFORCE)
   real(dp), dimension(:), allocatable :: fration_i !(nelec)
   ! fragments
   integer :: nfrag
   integer, dimension(:), allocatable :: ifragcent !(ncent)
   integer, dimension(:), allocatable :: ifragelec !(nelec)
   real(dp), dimension(:), allocatable :: elocfrag !(nfrag)
   real(dp), dimension(:,:,:), allocatable :: elocofrag !(nfrag,mwalk,MFORCE)
   real(dp), dimension(:), allocatable :: eestfrag !(nfrag)
   real(dp), dimension(:), allocatable :: esumfrag !(nfrag)
   real(dp), dimension(:), allocatable :: ecumfrag !(nfrag)
   real(dp), dimension(:), allocatable :: v2frag !(nfrag)
   real(dp), dimension(:), allocatable :: vav2frag !(nfrag)
   real(dp), dimension(:,:,:), allocatable :: fratiofrag !(nfrag,mwalk,MFORCE)
   real(dp), dimension(:), allocatable :: frationfrag !(nfrag)
   real(dp), dimension(:), allocatable :: potnnfrag !(nfrag)
   real(dp), dimension(:), allocatable :: ibranching_cfrag !(nfrag)
   real(dp), dimension(:), allocatable :: sqrt_nelecfrag !(nfrag)
   real(dp), dimension(:), allocatable :: tauefffrag !(nfrag)
   real(dp), dimension(:), allocatable :: dfus2acfrag !(nfrag)
   real(dp), dimension(:), allocatable :: dfus2unfrag !(nfrag)
   real(dp), dimension(:), allocatable :: etrialfrag !(nfrag)
   ! Accumulators for the fragments
   real(dp), dimension(:), allocatable :: egsum1frag !(nfrag)
   real(dp), dimension(:), allocatable :: egsumfrag !(nfrag)
   real(dp), dimension(:), allocatable :: egcum1frag !(nfrag)
   real(dp), dimension(:), allocatable :: egcumfrag !(nfrag)
   real(dp), dimension(:), allocatable :: egcm21frag !(nfrag)
   real(dp), dimension(:), allocatable :: egcm2frag !(nfrag)

   private
   ! electron fragments
   public :: eest_i, eloc_i, esum_i, ecum_i, eloco_i
   public :: v2_i, vav2_i, fratio_i, fration_i
   ! full fragments
   public :: nfrag, ifragcent, ifragelec, potnnfrag, sqrt_nelecfrag
   public :: eestfrag, elocfrag, esumfrag, ecumfrag, elocofrag
   public :: v2frag, vav2frag, fratiofrag, frationfrag, ibranching_cfrag
   public :: dfus2acfrag, dfus2unfrag, tauefffrag, etrialfrag
   public :: egsum1frag, egsumfrag, egcum1frag, egcumfrag, egcm21frag, egcm2frag

   public :: allocate_fragments, deallocate_fragments
   save

contains
   subroutine allocate_fragments
      ! electron fragments
      if (.not. allocated(eloc_i)) allocate(eloc_i(nelec))
      if (.not. allocated(eloco_i)) allocate(eloco_i(nelec,mwalk,MFORCE))
      if (.not. allocated(eest_i)) allocate(eest_i(nelec))
      if (.not. allocated(esum_i)) allocate(esum_i(nelec))
      if (.not. allocated(ecum_i)) allocate(ecum_i(nelec))
      if (.not. allocated(v2_i)) allocate(v2_i(nelec))
      if (.not. allocated(vav2_i)) allocate(vav2_i(nelec))
      if (.not. allocated(fration_i)) allocate(fration_i(nelec))
      if (.not. allocated(fratio_i)) allocate(fratio_i(nelec,mwalk,MFORCE))
      ! full fragments
      if (.not. allocated(ifragelec)) allocate(ifragelec(nelec))
      if (.not. allocated(ifragcent)) allocate(ifragcent(ncent))
      if (.not. allocated(sqrt_nelecfrag)) allocate(sqrt_nelecfrag(nfrag))
      if (.not. allocated(elocfrag)) allocate(elocfrag(nfrag))
      if (.not. allocated(elocofrag)) allocate(elocofrag(nfrag,mwalk,MFORCE))
      if (.not. allocated(eestfrag)) allocate(eestfrag(nfrag))
      if (.not. allocated(esumfrag)) allocate(esumfrag(nfrag))
      if (.not. allocated(ecumfrag))allocate(ecumfrag(nfrag))
      if (.not. allocated(v2frag)) allocate(v2frag(nfrag))
      if (.not. allocated(vav2frag)) allocate(vav2frag(nfrag))
      if (.not. allocated(fratiofrag)) allocate(fratiofrag(nfrag,mwalk,MFORCE))
      if (.not. allocated(frationfrag)) allocate(frationfrag(nfrag))
      if (.not. allocated(potnnfrag)) allocate(potnnfrag(nfrag))
      if (.not. allocated(ibranching_cfrag)) allocate(ibranching_cfrag(nfrag))
      if (.not. allocated(tauefffrag)) allocate(tauefffrag(nfrag))
      if (.not. allocated(dfus2acfrag)) allocate(dfus2acfrag(nfrag))
      if (.not. allocated(dfus2unfrag)) allocate(dfus2unfrag(nfrag))
      if (.not. allocated(etrialfrag)) allocate(etrialfrag(nfrag))
      if (.not. allocated(egsum1frag)) allocate(egsum1frag(nfrag))
      if (.not. allocated(egsumfrag)) allocate(egsumfrag(nfrag))
      if (.not. allocated(egcum1frag)) allocate(egcum1frag(nfrag))
      if (.not. allocated(egcumfrag)) allocate(egcumfrag(nfrag))
      if (.not. allocated(egcm21frag)) allocate(egcm21frag(nfrag))
      if (.not. allocated(egcm2frag)) allocate(egcm2frag(nfrag))
   end subroutine allocate_fragments

   subroutine deallocate_fragments
      ! electron fragments
      if (allocated(eloc_i)) deallocate(eloc_i)
      if (allocated(eest_i)) deallocate(eest_i)
      if (allocated(esum_i)) deallocate(esum_i)
      if (allocated(ecum_i)) deallocate(ecum_i)
      if (allocated(v2_i)) deallocate(v2_i)
      if (allocated(vav2_i)) deallocate(vav2_i)
      if (allocated(fration_i)) deallocate(fration_i)
      if (allocated(fratio_i)) deallocate(fratio_i)
      ! full fragments
      if (allocated(ifragelec)) deallocate(ifragelec)
      if (allocated(ifragcent)) deallocate(ifragcent)
      if (allocated(sqrt_nelecfrag)) deallocate(sqrt_nelecfrag)
      ! for some odd reason deallocating elocfrag gives a segmentation fault.
      if (allocated(elocfrag)) deallocate(elocfrag)
      if (allocated(elocofrag)) deallocate(elocofrag)
      if (allocated(eestfrag)) deallocate(eestfrag)
      if (allocated(esumfrag)) deallocate(esumfrag)
      if (allocated(ecumfrag)) deallocate(ecumfrag)
      if (allocated(v2frag)) deallocate(v2frag)
      if (allocated(vav2frag)) deallocate(vav2frag)
      if (allocated(fratiofrag)) deallocate(fratiofrag)
      if (allocated(frationfrag)) deallocate(frationfrag)
      if (allocated(potnnfrag)) deallocate(potnnfrag)
      if (allocated(ibranching_cfrag)) deallocate(ibranching_cfrag)
      if (allocated(tauefffrag)) deallocate(tauefffrag)
      if (allocated(dfus2acfrag)) deallocate(dfus2acfrag)
      if (allocated(dfus2unfrag)) deallocate(dfus2unfrag)
      if (allocated(etrialfrag)) deallocate(etrialfrag)
      if (allocated(egsum1frag)) deallocate(egsum1frag)
      if (allocated(egsumfrag)) deallocate(egsumfrag)
      if (allocated(egcum1frag)) deallocate(egcum1frag)
      if (allocated(egcumfrag)) deallocate(egcumfrag)
      if (allocated(egcm21frag)) deallocate(egcm21frag)
      if (allocated(egcm2frag)) deallocate(egcm2frag)
   end subroutine deallocate_fragments
end module fragments

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
      use system, only: nelec
      if (.not. allocated(fratio)) allocate(fratio(mwalk, MFORCE))
      if (.not. allocated(xdrifted)) allocate(xdrifted(3, nelec, mwalk, MFORCE))
   end subroutine allocate_velratio

   subroutine deallocate_velratio()
      if (allocated(fratio)) deallocate(fratio)
      if (allocated(xdrifted)) deallocate(xdrifted)
   end subroutine deallocate_velratio
 end module velratio
