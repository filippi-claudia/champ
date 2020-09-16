 module optorb_mod
  ! flags and dimensions for orbital optimization
  ! maximal number of terms, max dim of reduced matrices
  integer, parameter :: MXORBOP=8000
  integer, parameter :: MXREDUCED=1
  integer, parameter :: MXMATDIM=MXREDUCED*(MXREDUCED+1)
  integer, parameter :: MXMATDIM2=MXMATDIM/2
  integer, parameter :: MXREP=10
  private
  public :: MXORBOP, MXREDUCED, MXMATDIM, MXMATDIM2, MXREP
  save
 end module optorb_mod

module orb_mat_001
   !> Arguments: orb_o, orb_oe, orb_ho
   use optorb_mod, only: MXORBOP
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: orb_ho(MXORBOP,MSTATES)
   real(dp) :: orb_o(MXORBOP,MSTATES)
   real(dp) :: orb_oe(MXORBOP,MSTATES)

   private
   public :: orb_o, orb_oe, orb_ho
   save
 end module orb_mat_001

 module orb_mat_002
   !> Arguments: orb_ho_old, orb_o_old, orb_oe_old
   use optorb_mod, only: MXORBOP
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: orb_ho_old(MXORBOP,MSTATES)
   real(dp) :: orb_o_old(MXORBOP,MSTATES)
   real(dp) :: orb_oe_old(MXORBOP,MSTATES)

   private
   public :: orb_ho_old, orb_o_old, orb_oe_old
   save
 end module orb_mat_002

 module orb_mat_003
   !> Arguments: orb_o_sum, orb_o_cum
   use optorb_mod, only: MXORBOP
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: orb_o_cum(MXORBOP,MSTATES)
   real(dp) :: orb_o_sum(MXORBOP,MSTATES)

   private
   public :: orb_o_sum, orb_o_cum
   save
 end module orb_mat_003

 module orb_mat_004
   !> Arguments: orb_oe_sum, orb_oe_cum
   use optorb_mod, only: MXORBOP
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: orb_oe_cum(MXORBOP,MSTATES)
   real(dp) :: orb_oe_sum(MXORBOP,MSTATES)

   private
   public :: orb_oe_sum, orb_oe_cum
   save
 end module orb_mat_004

 module orb_mat_005
   !> Arguments: orb_ho_cum
   use optorb_mod, only: MXORBOP
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: orb_ho_cum(MXORBOP,MSTATES)

   private
   public :: orb_ho_cum
   save
 end module orb_mat_005

 module orb_mat_006
   !> Arguments: orb_oo_cum
   use optorb_mod, only: MXMATDIM2
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: orb_oo_cum(MXMATDIM2,MSTATES)

   private
   public :: orb_oo_cum
   save
 end module orb_mat_006

 module orb_mat_007
   !> Arguments: orb_oho_cum
   use optorb_mod, only: MXMATDIM
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: orb_oho_cum(MXMATDIM,MSTATES)

   private
   public :: orb_oho_cum
   save
 end module orb_mat_007

 module orb_mat_022
   use optorb_mod, only: MXORBOP
   !> Arguments: ideriv

   integer  :: ideriv(2,MXORBOP)

   private
   public :: ideriv
   save
 end module orb_mat_022

 module orb_mat_024
   !> Arguments: orb_oe_bsum, orb_f_bcum, orb_e_bsum, orb_w_bsum, orb_o_bsum, orb_f_bcm2
   use optorb_mod, only: MXORBOP
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: orb_e_bsum(MSTATES)
   real(dp) :: orb_f_bcm2(MXORBOP,MSTATES)
   real(dp) :: orb_f_bcum(MXORBOP,MSTATES)
   real(dp) :: orb_o_bsum(MXORBOP,MSTATES)
   real(dp) :: orb_oe_bsum(MXORBOP,MSTATES)
   real(dp) :: orb_w_bsum(MSTATES)

   private
   public :: orb_oe_bsum, orb_f_bcum, orb_e_bsum, orb_w_bsum, orb_o_bsum, orb_f_bcm2
   save
 end module orb_mat_024

 module orb_mat_030
   !> Arguments: orb_wcum, orb_ecum
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: orb_ecum(MSTATES)
   real(dp) :: orb_wcum(MSTATES)

   private
   public :: orb_wcum, orb_ecum
   save
 end module orb_mat_030

 module orb_mat_033
   use optorb_mod, only: MXORBOP
   !> Arguments: irepcol_ref, ideriv_ref, ideriv_iab

   integer  :: ideriv_iab(MXORBOP)
   integer  :: ideriv_ref(MXORBOP,2)
   integer  :: irepcol_ref(MXORBOP,2)

   private
   public :: irepcol_ref, ideriv_ref, ideriv_iab
   save
 end module orb_mat_033

 module optorb
   !> Arguments: dmat_diag, irrep, orb_energy
   use precision_kinds, only: dp
   use vmc, only: MORB

   real(dp) :: dmat_diag(MORB)
   integer  :: irrep(MORB)
   real(dp) :: orb_energy(MORB)

   private 
   public :: dmat_diag, irrep, orb_energy 
   save
 end module optorb

 module optorb_cblock   ! from optorb.h 
   ! norbterm: number of terms (possibly after a transformation)
   ! norbprim: number of primitive terms (determinant ratios)
   integer :: norbterm
   integer :: norbprim

   ! PLT: From old common block /orb004/
   integer :: nefp_blocks
   integer :: nb_current
   integer :: norb_f_bcum

   ! reduced correlation matrix pointers
   ! threshold in terms of std dev. , limit for keeping operators
   ! if iuse_trafo: linearly transformed operators sampled instead of primitive
   !     replacement operators 
   ! PLT: From old common blocks /orb006/ and /orb008/.
   integer :: isample_cmat
   integer :: nreduced 
   integer :: iuse_trafiuse_trafoo 

   ! Dumping block averages for error analysis. 
   ! PLT: From old common blocks /orb009/ and /orb010/.
   integer :: idump_blockav 
   integer :: iorbsample 
   integer :: ns_current

   ! Printing flags:
   integer :: iorbprt
   integer :: iorbprt_sav

   private 
   public :: norbterm, norbprim
   public :: nefp_blocks, nb_current, norb_f_bcum
   public :: isample_cmat, nreduced, iuse_trafiuse_trafoo 
   public :: idump_blockav, iorbsample, ns_current 
   public :: iorbprt, iorbprt_sav
   save
 end module optorb_cblock 

 module optorb_mix
   !> Arguments: iwmix_virt, norbopt, norbvirt
   use vmc, only: MORB

   integer  :: iwmix_virt(MORB,MORB)
   integer  :: norbopt
   integer  :: norbvirt

   private
   public :: iwmix_virt, norbopt, norbvirt
   save
 end module optorb_mix