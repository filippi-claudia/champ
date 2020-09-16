 module ewald_mod
  ! NP         is the number of factors of (r/cutr-1) we multiply polynomial by
  ! IVOL_RATIO is the ratio of the simulation to primitive cell volumes
  ! IBIG_RATIO is the ratio of the number of vectors or norms before the optimal Ewald separation
  !            to the number after the separation
  ! NSYM       is the ratio of the number of vectors to the number of norms
  !            and depends on the symmetry of the lattice.
  integer, parameter :: NCOEFX=10, NPX=4, IVOL_RATIO=10, IBIG_RATIO=15, NSYM=8
  integer, parameter :: NGNORMX=1000, NGVECX=NGNORMX*NSYM, NG1DX=60
  integer, parameter :: NGNORM_SIMX=NGNORMX*IVOL_RATIO, NGVEC_SIMX=NGVECX*IVOL_RATIO
  integer, parameter :: NGNORM_BIGX=IBIG_RATIO*NGNORMX, NGVEC_BIGX=IBIG_RATIO*NGVECX
  integer, parameter :: NGNORM_SIM_BIGX=IBIG_RATIO*NGNORM_SIMX, NGVEC_SIM_BIGX=IBIG_RATIO*NGVEC_SIMX

  private 
  public :: NCOEFX, NPX, IVOL_RATIO, IBIG_RATIO, NSYM
  public :: NGNORMX, NGVECX, NG1DX
  public :: NGNORM_SIMX, NGVEC_SIMX
  public :: NGNORM_BIGX, NGVEC_BIGX
  public :: NGNORM_SIM_BIGX, NGVEC_SIM_BIGX
  save 
 end module ewald_mod
    
 module ewald
   !> Arguments: b_coul, b_coul_sim, y_coul, y_coul_sim
   use precision_kinds, only: dp
   use ewald_mod, only: NCOEFX, NGNORMX, NGNORM_SIMX

   real(dp) :: b_coul(NCOEFX)
   real(dp) :: b_coul_sim(NCOEFX)
   real(dp) :: y_coul(NGNORMX)
   real(dp) :: y_coul_sim(NGNORM_SIMX)

   private
   public   ::  b_coul, b_coul_sim, y_coul, y_coul_sim
   save
 end module ewald

 module ewald_basis
   !> Arguments: vps_basis_fourier
   use precision_kinds, only: dp
   use ewald_mod, only: NGNORM_BIGX

   real(dp) :: vps_basis_fourier(NGNORM_BIGX)

   private
   public   :: vps_basis_fourier
   save
 end module ewald_basis


 module periodic
   !> Arguments: cutg, cutg_big, cutg_sim, cutg_sim_big, cutr, cutr_sim, glatt, glatt_inv, glatt_sim, gnorm, gnorm_sim, gvec, gvec_sim, igmult, igmult_sim, igvec, igvec_sim, ireal_imag, isrange, k_inv, kvec, nband, ncoef_per, ng1d, ng1d_sim, ngnorm, ngnorm_big, ngnorm_orb, ngnorm_sim, ngnorm_sim_big, ngvec, ngvec_big, ngvec_orb, ngvec_sim, ngvec_sim_big, nkvec, np, npoly, rknorm, rkvec, rkvec_shift, rlatt, rlatt_inv, rlatt_sim, rlatt_sim_inv, vcell, vcell_sim, znuc2_sum, znuc_sum
   use ewald_mod, only: IVOL_RATIO
   use ewald_mod, only: NGNORM_BIGX, NGVEC_BIGX
   use ewald_mod, only: NGNORM_SIM_BIGX, NGVEC_SIM_BIGX
   use precision_kinds, only: dp
   use vmc, only: MORB

   real(dp) :: cutg
   real(dp) :: cutg_big
   real(dp) :: cutg_sim
   real(dp) :: cutg_sim_big
   real(dp) :: cutr
   real(dp) :: cutr_sim
   real(dp) :: glatt(3,3)
   real(dp) :: glatt_inv(3,3)
   real(dp) :: glatt_sim(3,3)
   real(dp) :: gnorm(NGNORM_BIGX)
   real(dp) :: gnorm_sim(NGNORM_SIM_BIGX)
   real(dp) :: gvec(3,NGVEC_BIGX)
   real(dp) :: gvec_sim(3,NGVEC_SIM_BIGX)
   integer  :: igmult(NGNORM_BIGX)
   integer  :: igmult_sim(NGNORM_SIM_BIGX)
   integer  :: igvec(3,NGVEC_BIGX)
   integer  :: igvec_sim(3,NGVEC_SIM_BIGX)
   integer  :: ireal_imag(MORB)
   integer  :: isrange
   integer  :: k_inv(IVOL_RATIO)
   integer  :: kvec(3,IVOL_RATIO)
   integer  :: nband(IVOL_RATIO)
   integer  :: ncoef_per
   integer  :: ng1d(3)
   integer  :: ng1d_sim(3)
   integer  :: ngnorm
   integer  :: ngnorm_big
   integer  :: ngnorm_orb
   integer  :: ngnorm_sim
   integer  :: ngnorm_sim_big
   integer  :: ngvec
   integer  :: ngvec_big
   integer  :: ngvec_orb
   integer  :: ngvec_sim
   integer  :: ngvec_sim_big
   integer  :: nkvec
   integer  :: np
   integer  :: npoly
   real(dp) :: rknorm(IVOL_RATIO)
   real(dp) :: rkvec(3,IVOL_RATIO)
   real(dp) :: rkvec_shift(3)
   real(dp) :: rlatt(3,3)
   real(dp) :: rlatt_inv(3,3)
   real(dp) :: rlatt_sim(3,3)
   real(dp) :: rlatt_sim_inv(3,3)
   real(dp) :: vcell
   real(dp) :: vcell_sim
   real(dp) :: znuc2_sum
   real(dp) :: znuc_sum

   private
   public :: cutg, cutg_big, cutg_sim, cutg_sim_big, cutr, cutr_sim, glatt, glatt_inv
   public :: glatt_sim, gnorm, gnorm_sim, gvec, gvec_sim, igmult, igmult_sim, igvec
   public :: igvec_sim, ireal_imag, isrange, k_inv, kvec, nband, ncoef_per, ng1d, ng1d_sim
   public :: ngnorm, ngnorm_big, ngnorm_orb, ngnorm_sim, ngnorm_sim_big, ngvec, ngvec_big
   public :: ngvec_orb, ngvec_sim, ngvec_sim_big, nkvec, np, npoly, rknorm, rkvec, rkvec_shift
   public :: rlatt, rlatt_inv, rlatt_sim, rlatt_sim_inv, vcell, vcell_sim, znuc2_sum, znuc_sum
   save
 end module periodic

  module pworbital
   !> Arguments: c_im, c_ip, c_rm, c_rp, icmplx, isortg, isortk, ngorb
   use ewald_mod, only: IVOL_RATIO
   use ewald_mod, only: NGVECX
   use precision_kinds, only: dp
   use vmc, only: MORB

   real(dp) :: c_im(NGVECX,MORB)
   real(dp) :: c_ip(NGVECX,MORB)
   real(dp) :: c_rm(NGVECX,MORB)
   real(dp) :: c_rp(NGVECX,MORB)
   integer  :: icmplx
   integer  :: isortg(NGVECX,MORB)
   integer  :: isortk(IVOL_RATIO)
   integer  :: ngorb(MORB)

   private
   public :: c_im, c_ip, c_rm, c_rp, icmplx, isortg, isortk, ngorb
   save
 end module pworbital