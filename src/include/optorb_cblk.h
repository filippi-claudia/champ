!     Common blocks (matrices) for efpci 
! 
! ratios  determinant/twf (efpci operators)
! old ratios

! sum, cum, square of operators
      common /orb_mat_003/ orb_o_sum(MXORBOP,MSTATES),orb_o_cum(MXORBOP,MSTATES)
      common /orb_mat_003/ orb_o_cm2(MXORBOP,MSTATES)

! sum, cum, square of operator * energy
      common /orb_mat_004/ orb_oe_sum(MXORBOP,MSTATES),orb_oe_cum(MXORBOP,MSTATES)
      common /orb_mat_004/ orb_oe_cm2(MXORBOP,MSTATES)

! sum, cum of operator products (no error on products anymore)
      common /orb_mat_005/ orb_ho_cum(MXORBOP,MSTATES)
      common /orb_mat_006/ orb_oo_cum(MXMATDIM2,MSTATES)
      common /orb_mat_007/ orb_oho_cum(MXMATDIM,MSTATES)

! block average of 'forces' : number of blocks in av, sum and sum**2
      common /orb_mat_024/ orb_e_bsum(MSTATES),orb_w_bsum(MSTATES)
      common /orb_mat_024/ orb_o_bsum(MXORBOP,MSTATES),orb_oe_bsum(MXORBOP,MSTATES)
      common /orb_mat_024/ orb_f_bcum(MXORBOP,MSTATES),orb_f_bcm2(MXORBOP,MSTATES)
      common /orb_mat_021/ i_have_cmat(MXREDUCED)
      common /orb_mat_022/ ideriv(2,MXORBOP)
      common /orb_mat_030/ orb_wcum(MSTATES),orb_ecum(MSTATES)
      common /orb_mat_033/ ideriv_iab(MXORBOP),ideriv_ref(MXORBOP,2),irepcol_ref(MXORBOP,2)
