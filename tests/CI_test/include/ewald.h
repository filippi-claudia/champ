! NP         is the number of factors of (r/cutr-1) we multiply polynomial by
! IVOL_RATIO is the ratio of the simulation to primitive cell volumes
! IBIG_RATIO is the ratio of the number of vectors or norms before the optimal Ewald separation
!            to the number after the separation
! NSYM       is the ratio of the number of vectors to the number of norms
!            and depends on the symmetry of the lattice.
      parameter (NCOEFX=10, NPX=4, IVOL_RATIO=10, IBIG_RATIO=15, NSYM=8)
      parameter (NGNORMX=1000, NGVECX=NGNORMX*NSYM, NG1DX=60)
      parameter (NGNORM_SIMX=NGNORMX*IVOL_RATIO, NGVEC_SIMX=NGVECX*IVOL_RATIO)
      parameter (NGNORM_BIGX=IBIG_RATIO*NGNORMX, NGVEC_BIGX=IBIG_RATIO*NGVECX)
      parameter (NGNORM_SIM_BIGX=IBIG_RATIO*NGNORM_SIMX, NGVEC_SIM_BIGX=IBIG_RATIO*NGVEC_SIMX)

!    &,NGNORMX=22, NGVECX=266, NG1DX=15
!    &,NGNORM_SIMX=54, NGVEC_SIMX=1052
!    &,NGNORM_BIGX=86, NGVEC_BIGX=2140
!    &,NGNORM_SIM_BIGX=213, NGVEC_SIM_BIGX=8440)

!    &,NGNORMX=100000, NGVECX=100000, NG1DX=15
!    &,NGNORM_SIMX=100000, NGVEC_SIMX=100000
!    &,NGNORM_BIGX=100000, NGVEC_BIGX=100000
!    &,NGNORM_SIM_BIGX=100000, NGVEC_SIM_BIGX=100000)

!cNP         is the number of factors of (r/cutr-1) we multiply polynomial by
!cIVOL_RATIO is the ratio of the simulation to primitive cell volumes
!cNSYM       is the ratio of the number of vectors to the number of norms
!c           and depends on the symmetry of the lattice.
!c    parameter (NCOEFX=10, NP=5, IVOL_RATIO=4, NSYM=48
!     parameter (NCOEFX=10, IVOL_RATIO=100, NSYM=48
!    &,NGNORMX=300, NGVECX=NGNORMX*NSYM, NG1DX=50
!    &,NGNORM_SIMX=NGNORMX*IVOL_RATIO, NGVEC_SIMX=NGVECX*IVOL_RATIO
!    &,NGNORM_BIGX=8*NGNORMX, NGVEC_BIGX=8*NGVECX
!    &,NGNORM_SIM_BIGX=8*NGNORM_SIMX, NGVEC_SIM_BIGX=8*NGVEC_SIMX)

!c   &,NGNORMX=22, NGVECX=266, NG1DX=15
!c   &,NGNORM_SIMX=54, NGVEC_SIMX=1052
!c   &,NGNORM_BIGX=86, NGVEC_BIGX=2140
!c   &,NGNORM_SIM_BIGX=213, NGVEC_SIM_BIGX=8440)

!c   &,NGNORMX=100000, NGVECX=100000, NG1DX=15
!c   &,NGNORM_SIMX=100000, NGVEC_SIMX=100000
!c   &,NGNORM_BIGX=100000, NGVEC_BIGX=100000
!c   &,NGNORM_SIM_BIGX=100000, NGVEC_SIM_BIGX=100000)
