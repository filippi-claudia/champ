c     flags and dimensions for generalized CI expectation values
c
c     maximal number of terms, max dim of reduced matrices
      parameter (MXCITERM=1100,MXCIREDUCED=1,MXCIMATDIM=MXCITERM*(MXCIREDUCED+1)/2)
c     ici  : flag whether CI should be sampled
c     nciterm: number of terms (possibly after a transformation)
c     nciterm2:number of terms for Nightingale's method (EL(i) O(k))
c     nciprim: number of primitive terms (determinant ratios)
c     iciprt : print flag 
      common /ci000/ nciterm, nciprim, iciprt
