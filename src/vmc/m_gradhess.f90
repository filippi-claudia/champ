module gradhess_all
   !> Arguments: MPARMALL, grad, h, s
   use optorb_mod, only: MXREDUCED
   use optjas, only: MPARMJ
   use optci, only: MXCIREDUCED
   use precision_kinds, only: dp

   integer, parameter :: MPARMALL=MPARMJ+MXCIREDUCED+MXREDUCED
   real(dp) :: grad(MPARMALL)
   real(dp) :: h(MPARMALL,MPARMALL)
   real(dp) :: s(MPARMALL,MPARMALL)

   private
   public :: MPARMALL, grad, h, s
   save
 end module gradhess_all

 module gradhessj
   !> Arguments: d2j, d2j_e, de, de_de, de_e, dj, dj_de, dj_dj, dj_dj_e, dj_e, dj_e2, e2
   use optjas, only: MPARMJ
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: d2j(MPARMJ,MPARMJ,MSTATES)
   real(dp) :: d2j_e(MPARMJ,MPARMJ,MSTATES)
   real(dp) :: de(MPARMJ,MSTATES)
   real(dp) :: de_de(MPARMJ,MPARMJ,MSTATES)
   real(dp) :: de_e(MPARMJ,MSTATES)
   real(dp) :: dj(MPARMJ,MSTATES)
   real(dp) :: dj_de(MPARMJ,MPARMJ,MSTATES)
   real(dp) :: dj_dj(MPARMJ,MPARMJ,MSTATES)
   real(dp) :: dj_dj_e(MPARMJ,MPARMJ,MSTATES)
   real(dp) :: dj_e(MPARMJ,MSTATES)
   real(dp) :: dj_e2(MPARMJ,MSTATES)
   real(dp) :: e2(MPARMJ,MSTATES)

   private 
   public :: d2j, d2j_e, de, de_de, de_e, dj, dj_de, dj_dj, dj_dj_e, dj_e, dj_e2, e2 
   save
 end module gradhessj

 module gradhessjo
   !> Arguments: d1d2a_old, d1d2b_old, d2d2a_old, d2d2b_old, denergy_old, gvalue_old
   use optjas, only: MPARMJ
   use precision_kinds, only: dp
   use vmc_mod, only: MCTYPE
   use mstates_mod, only: MSTATES

   real(dp) :: d1d2a_old(MCTYPE)
   real(dp) :: d1d2b_old(2)
   real(dp) :: d2d2a_old(MCTYPE)
   real(dp) :: d2d2b_old(2)
   real(dp) :: denergy_old(MPARMJ,MSTATES)
   real(dp) :: gvalue_old(MPARMJ)

   private
   public   ::  d1d2a_old, d1d2b_old, d2d2a_old, d2d2b_old, denergy_old, gvalue_old
   save
 end module gradhessjo

 module gradhess_ci
   !> Arguments: grad_ci, h_ci, s_ci
   use optci, only: MXCITERM, MXCIREDUCED
   use precision_kinds, only: dp

   real(dp) :: grad_ci(MXCITERM)
   real(dp) :: h_ci(MXCITERM,MXCIREDUCED)
   real(dp) :: s_ci(MXCITERM,MXCIREDUCED)

   private
   public   ::  grad_ci, h_ci, s_ci
   save
 end module gradhess_ci

 module gradhess_jas
   !> Arguments: grad_jas, h_jas, s_jas
   use optjas, only: MPARMJ
   use precision_kinds, only: dp
 
   real(dp) :: grad_jas(MPARMJ)
   real(dp) :: h_jas(MPARMJ,MPARMJ)
   real(dp) :: s_jas(MPARMJ,MPARMJ)
 
   private
   public   ::  grad_jas, h_jas, s_jas
   save
 end module gradhess_jas
 
 module gradhess_mix_jas_ci
   !> Arguments: h_mix_jas_ci, s_mix_jas_ci
   use optjas, only: MPARMJ
   use optci, only: MXCITERM
   use precision_kinds, only: dp

   real(dp) :: h_mix_jas_ci(2*MPARMJ,MXCITERM)
   real(dp) :: s_mix_jas_ci(MPARMJ,MXCITERM)

   private
   public   ::  h_mix_jas_ci, s_mix_jas_ci
   save
 end module gradhess_mix_jas_ci

 module gradhess_mix_jas_orb
   !> Arguments: h_mix_jas_orb, s_mix_jas_orb
   use optorb_mod, only: MXREDUCED
   use optjas, only: MPARMJ
   use precision_kinds, only: dp

   real(dp) :: h_mix_jas_orb(2*MPARMJ,MXREDUCED)
   real(dp) :: s_mix_jas_orb(MPARMJ,MXREDUCED)

   private
   public   ::  h_mix_jas_orb, s_mix_jas_orb
   save
 end module gradhess_mix_jas_orb

 module gradhess_mix_orb_ci
   !> Arguments: h_mix_ci_orb, s_mix_ci_orb
   use optorb_mod, only: MXREDUCED
   use optci, only: MXCITERM
   use precision_kinds, only: dp

   real(dp) :: h_mix_ci_orb(2*MXCITERM,MXREDUCED)
   real(dp) :: s_mix_ci_orb(MXCITERM,MXREDUCED)

   private
   public   ::  h_mix_ci_orb, s_mix_ci_orb
   save
 end module gradhess_mix_orb_ci

 module gradjerr
   !> Arguments: dj_bsum, dj_e_bsum, dj_e_save, dj_save, e_bsum, grad_jas_bcm2, grad_jas_bcum
   use optjas, only: MPARMJ
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: dj_bsum(MPARMJ,MSTATES)
   real(dp) :: dj_e_bsum(MPARMJ,MSTATES)
   real(dp) :: dj_e_save(MPARMJ,MSTATES)
   real(dp) :: dj_save(MPARMJ,MSTATES)
   real(dp) :: e_bsum(MSTATES)
   real(dp) :: grad_jas_bcm2(MPARMJ,MSTATES)
   real(dp) :: grad_jas_bcum(MPARMJ,MSTATES)

   private
   public   ::  dj_bsum, dj_e_bsum, dj_e_save, dj_save, e_bsum, grad_jas_bcm2, grad_jas_bcum
   save
 end module gradjerr