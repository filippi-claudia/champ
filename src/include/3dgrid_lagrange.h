      integer LAGMAX, LAGSTART, LAGEND, MORB_OCC
!     Number of Lagrange interpolation points/axis
      parameter (LAGMAX=4)
      parameter (LAGSTART=-LAGMAX/2)
      parameter (LAGEND  =LAGSTART+LAGMAX-1)
      parameter (MORB_OCC=MELEC/2)

!     Spline fits of the orbitals
!     and boundary conditions (for the creation of the fit)

      real*4 orb_num_lag(5,MXNSTEP,MXNSTEP,MXNSTEP,MORB_OCC)
      common /orbital_num_lag/ denom(LAGSTART:LAGEND,3), step_inv(3,3)

