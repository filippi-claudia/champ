      integer LAGMAX, LAGSTART, LAGEND, MORB_OCC
c     Number of Lagrange interpolation points/axis
      parameter (LAGMAX=4)
      parameter (LAGSTART=-LAGMAX/2)
      parameter (LAGEND  =LAGSTART+LAGMAX-1)
      parameter (MORB_OCC=MELEC/2)

c     Spline fits of the orbitals
c     and boundary conditions (for the creation of the fit)

      real*4 orb_num_lag(5,MXNSTEP,MXNSTEP,MXNSTEP,MORB_OCC)
      common /orbital_num_lag/ orb_num_lag
     &,denom(LAGSTART:LAGEND,3), step_inv(3,3)

