c     flags and dimensions for the 3d grid objects

      integer MXNSTEP,MXNSTEP2,MXNSTEP3
      integer IUNDEFINED
      real*8  UNDEFINED, SHIFT

c     parameter (MXNSTEP=50)
      parameter (MXNSTEP=1)
      parameter (UNDEFINED = -1234567890.d0)
      parameter (IUNDEFINED = -1234567890)
      parameter (SHIFT = 2.d0)

      parameter (MXNSTEP2=MXNSTEP*MXNSTEP)
      parameter (MXNSTEP3=MXNSTEP2*MXNSTEP)

      real*4 grid3d(MXNSTEP,MXNSTEP,MXNSTEP), cart_from_int(MXNSTEP,3)
      common /grid3d_param/ origin(3), endpt(3), step3d(3), nstep3d(3)
      common /grid3d_array/ cart_from_int
      common /grid3d_data/  grid3d



