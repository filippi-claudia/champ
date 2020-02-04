c     flags and dimensions for the 3d grid objects

      integer MGRID_PCM,MGRID_PCM2,MGRID_PCM3
      integer IUNDEFINED
      real*8  UNDEFINED, SHIFT

      parameter (MGRID_PCM=1)
      parameter (UNDEFINED = -1234567890.d0)
      parameter (IUNDEFINED = -1234567890)

      parameter (MGRID_PCM2=MGRID_PCM*MGRID_PCM)
      parameter (MGRID_PCM3=MGRID_PCM2*MGRID_PCM)

      real*8 pcm_cart_from_int(MGRID_PCM,3)
      real*8 pcm_num_spl(8,MGRID_PCM,MGRID_PCM,MGRID_PCM)
      common /pcm_grid3d_contrl/ ipcm_3dgrid
      common /pcm_grid3d_param/ pcm_origin(3), pcm_endpt(3), pcm_step3d(3), ipcm_nstep3d(3)
      common /pcm_num_spl/ pcm_num_spl
      common /pcm_grid3d_array/ pcm_cart_from_int

