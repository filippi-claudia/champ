      DOUBLE PRECISION FUNCTION RANF(ISEED)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c     generates pseudo-random numbers in the interval (0,1)           ::
c     replacing the 1st 2 lines by the 3rd saves less than 10% of the ::
c     time on the vax (UNIX BSD4.2).  However replacing the function  ::
c     call by inline code results in more than a factor of 2 saving.  ::
c     The function call takes 10**-4 secs on a VAX 750 and 10**-5 secs::
c     on a GOULD 9000.                                                ::
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c     ISEED=ISEED*314159269+453806245
c     IF(ISEED.LT.0) ISEED=ISEED+2147483647+1
      ISEED=iand(ISEED*314159269+453806245,2147483647)
      RANF=DFLOAT(ISEED)*0.465661287524579692D-9
      RETURN
      END

