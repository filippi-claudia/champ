!     Common blocks (matrices) for ci 
! 
!     ratios  determinant/twf (ci operators)
!     old ratios
      common /ci002_blk/ ci_o_old(MXCITERM),ci_oe_old(MXCITERM,MXCIREDUCED)
!
      common /ci003_blk/ ci_e(MXCITERM),ci_e_old(MXCITERM)
!
      common /ci004_blk/ ci_de(MXCITERM),ci_de_old(MXCITERM)
!     sum, cum of operator 
      common /ci005_blk/ ci_o_sum(MXCITERM),ci_o_cum(MXCITERM)
!
      common /ci006_blk/ ci_de_sum(MXCITERM),ci_de_cum(MXCITERM)
!     sum, cum, square of operator * energy
      common /ci008_blk/ ci_oe_sum(MXCITERM,MXCIREDUCED),ci_oe_cum(MXCITERM,MXCIREDUCED)
      common /ci008_blk/ ci_oe_cm2(MXCITERM,MXCIREDUCED)
!     sum, cum of operator products (no error on products anymore)
      common /ci009_blk/ ci_oo_sum(MXCIMATDIM),ci_oo_cum(MXCIMATDIM)
      common /ci009_blk/ ci_oo_cm2(MXCIMATDIM)
!     sum, cum of operator products (no error on products anymore)
      common /ci010_blk/ ci_ooe_sum(MXCIMATDIM),ci_ooe_cum(MXCIMATDIM)
