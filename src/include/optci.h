!     flags and dimensions for generalized CI expectation values
!
!     maximal number of terms, max dim of reduced matrices
      parameter (MXCITERM=1000,MXCIREDUCED=1,MXCIMATDIM=MXCITERM*(MXCIREDUCED+1)/2)
!     ici  : flag whether CI should be sampled
!     nciterm: number of terms (possibly after a transformation)
!     nciterm2:number of terms for Nightingale's method (EL(i) O(k))
!     nciprim: number of primitive terms (determinant ratios)
!     iciprt : print flag 
      common /ci000/ nciterm, nciprim, iciprt
