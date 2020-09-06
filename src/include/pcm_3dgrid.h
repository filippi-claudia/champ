!     flags and dimensions for the 3d grid objects

      integer MGRID_PCM,MGRID_PCM2,MGRID_PCM3
      integer IUNDEFINED
      real*8  UNDEFINED, SHIFT

      parameter (MGRID_PCM=1)
      parameter (UNDEFINED = -1234567890.d0)
      parameter (IUNDEFINED = -1234567890)

      parameter (MGRID_PCM2=MGRID_PCM*MGRID_PCM)
      parameter (MGRID_PCM3=MGRID_PCM2*MGRID_PCM)


