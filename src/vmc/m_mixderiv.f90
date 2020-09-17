
 module mix_jas_ci
   !> Arguments: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci
   use optjas, only: MPARMJ
   use precision_kinds, only: dp
   use vmc, only: MDET

   real(dp) :: de_o_ci(MPARMJ,MDET)
   real(dp) :: dj_de_ci(MPARMJ,MDET)
   real(dp) :: dj_o_ci(MPARMJ,MDET)
   real(dp) :: dj_oe_ci(MPARMJ,MDET)

   private 
   public :: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci 
   save
 end module mix_jas_ci

 module mix_jas_orb
   !> Arguments: de_o, dj_ho, dj_o, dj_oe
   use optorb_mod, only: MXREDUCED
   use optjas, only: MPARMJ
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES
 
   real(dp) :: de_o(MPARMJ,MXREDUCED,MSTATES)
   real(dp) :: dj_ho(MPARMJ,MXREDUCED,MSTATES)
   real(dp) :: dj_o(MPARMJ,MXREDUCED,MSTATES)
   real(dp) :: dj_oe(MPARMJ,MXREDUCED,MSTATES)
 
   private 
   public :: de_o, dj_ho, dj_o, dj_oe 
   save
 end module mix_jas_orb

 module mix_orb_ci
   !> Arguments: ci_de_o, ci_o_ho, ci_o_o, ci_o_oe
   use optorb_mod, only: MXREDUCED
   use optci, only: MXCITERM
   use precision_kinds, only: dp

   real(dp) :: ci_de_o(MXCITERM,MXREDUCED)
   real(dp) :: ci_o_ho(MXCITERM,MXREDUCED)
   real(dp) :: ci_o_o(MXCITERM,MXREDUCED)
   real(dp) :: ci_o_oe(MXCITERM,MXREDUCED)

   private 
   public :: ci_de_o, ci_o_ho, ci_o_o, ci_o_oe 
   save
 end module mix_orb_ci