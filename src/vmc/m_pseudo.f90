 module pseudo_mod
  !>Arguments : MPS_L, MPS_QUAD, MPS_GRID, MGAUSS
  integer, parameter :: MPS_L=4
  integer, parameter :: MPS_QUAD=86
  integer, parameter :: MPS_GRID=1200
  ! for gauss ecp
  integer, parameter :: MGAUSS=100 

  private 
  public :: MPS_L, MPS_QUAD, MPS_GRID, MGAUSS
  save
 end module pseudo_mod
 
module pseudo
   !> Arguments: lpot, nloc, vps, vpso
   use pseudo_mod, only: MPS_L
   use force_mod, only: MFORCE
   use precision_kinds, only: dp
   use vmc, only: MELEC, MCENT, MCTYPE

   integer  :: lpot(MCTYPE)
   integer  :: nloc
   real(dp) :: vps(MELEC,MCENT,MPS_L)
   real(dp) :: vpso(MELEC,MCENT,MPS_L,MFORCE)

   private 
   public :: lpot, nloc, vps, vpso 
   save
end module pseudo

 module pseudo_champ
   !> Arguments: igrid_ps, rmax_coul, rmax_nloc
   use precision_kinds, only: dp
   use vmc, only: MCTYPE

   integer  :: igrid_ps(MCTYPE)
   real(dp) :: rmax_coul(MCTYPE)
   real(dp) :: rmax_nloc(MCTYPE)

   private
   public :: igrid_ps, rmax_coul, rmax_nloc
   save
 end module pseudo_champ

 module pseudo_fahy
   !> Arguments: drad, dradl, nlrad, npotl, potl, ptnlc, rcmax
   use pseudo_mod, only: MPS_L, MPS_GRID
   use precision_kinds, only: dp
   use vmc, only: MCTYPE

   real(dp) :: drad(MCTYPE)
   real(dp) :: dradl(MCTYPE)
   integer  :: nlrad(MCTYPE)
   integer  :: npotl(MCTYPE)
   real(dp) :: potl(MPS_GRID,MCTYPE)
   real(dp) :: ptnlc(MPS_GRID,MCTYPE,MPS_L)
   real(dp) :: rcmax(MCTYPE)

   private
   public :: drad, dradl, nlrad, npotl, potl, ptnlc, rcmax
   save
 end module pseudo_fahy

 module pseudo_tm
   !> Arguments: arg, arg_ps, d2pot, nr_ps, r0, r0_ps, rmax, rmax_ps, vpseudo
   use pseudo_mod, only: MPS_L, MPS_GRID
   use precision_kinds, only: dp
   use vmc, only: MCTYPE

    real(dp) :: arg(MCTYPE)
    real(dp) :: arg_ps(MCTYPE)
    real(dp) :: d2pot(MPS_GRID,MCTYPE,MPS_L)
    integer  :: nr_ps(MCTYPE)
    real(dp) :: r0(MCTYPE)
    real(dp) :: r0_ps(MCTYPE)
    real(dp) :: rmax(MCTYPE)
    real(dp) :: rmax_ps(MCTYPE)
    real(dp) :: vpseudo(MPS_GRID,MCTYPE,MPS_L)

    private
    public :: arg, arg_ps, d2pot, nr_ps, r0, r0_ps, rmax, rmax_ps, vpseudo
    save
 end module pseudo_tm