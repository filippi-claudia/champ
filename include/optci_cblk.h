c     Common blocks (matrices) for ci 
c 
c     ratios  determinant/twf (ci operators)
      common /ci001_blk/ ci_o(MXCITERM),ci_oe(MXCITERM,MXCIREDUCED)
c     old ratios
      common /ci002_blk/ ci_o_old(MXCITERM),ci_oe_old(MXCITERM,MXCIREDUCED)
c
      common /ci003_blk/ ci_e(MXCITERM),ci_e_old(MXCITERM)
c
      common /ci004_blk/ ci_de(MXCITERM),ci_de_old(MXCITERM)
c     sum, cum of operator 
      common /ci005_blk/ ci_o_sum(MXCITERM),ci_o_cum(MXCITERM)
c
      common /ci006_blk/ ci_de_sum(MXCITERM),ci_de_cum(MXCITERM)
c     sum, cum, square of operator * energy
      common /ci008_blk/ ci_oe_sum(MXCITERM,MXCIREDUCED),ci_oe_cum(MXCITERM,MXCIREDUCED),
     &               ci_oe_cm2(MXCITERM,MXCIREDUCED)
c     sum, cum of operator products (no error on products anymore)
      common /ci009_blk/ ci_oo_sum(MXCIMATDIM),ci_oo_cum(MXCIMATDIM),
     &               ci_oo_cm2(MXCIMATDIM)
c     sum, cum of operator products (no error on products anymore)
      common /ci010_blk/ ci_ooe_sum(MXCIMATDIM),ci_ooe_cum(MXCIMATDIM)
