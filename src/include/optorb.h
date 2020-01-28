! flags and dimensions for orbital optimization
!
! maximal number of terms, max dim of reduced matrices
      parameter( MXORBOP=8000,MXREDUCED=8000)
      parameter( MXMATDIM=MXREDUCED*(MXREDUCED+1),MXMATDIM2=MXMATDIM/2)
      parameter( MXREP=10)
! norbterm: number of terms (possibly after a transformation)
! norbprim: number of primitive terms (determinant ratios)
! iorbprt : print flag 
      common /orb000/ norbterm, norbprim
      common /orb004/ nefp_blocks,nb_current,norb_f_bcum
! reduced correlation matrix pointers
! threshold in terms of std dev. , limit for keeping operators
! if iuse_trafo: linearly transformed operators sampled instead of primitive
!     replacement operators 
      common /orb006/ isample_cmat,nreduced
      common /orb008/ iuse_trafo
! dumping block averages for error analysis
      common /orb009/ idump_blockav
      common /orb010/ iorbsample,ns_current
