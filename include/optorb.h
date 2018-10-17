c flags and dimensions for orbital optimization
c
c maximal number of terms, max dim of reduced matrices
      parameter (MXORBOP=3000,MXREDUCED=1,
     &     MXMATDIM=MXREDUCED*(MXREDUCED+1),MXMATDIM2=MXMATDIM/2,
     &     MXREP=10)
c norbterm: number of terms (possibly after a transformation)
c norbprim: number of primitive terms (determinant ratios)
c iorbprt : print flag 
      common /orb000/ norbterm, norbprim
      common /orb004/ nefp_blocks,nb_current,norb_f_bcum
c reduced correlation matrix pointers
c threshold in terms of std dev. , limit for keeping operators
c if iuse_trafo: linearly transformed operators sampled instead of primitive
c     replacement operators 
      common /orb006/ isample_cmat,nreduced
      common /orb008/ iuse_trafo
c dumping block averages for error analysis
      common /orb009/ idump_blockav
      common /orb010/ iorbsample,ns_current
