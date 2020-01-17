      SUBROUTINE SHELL(D,N)
      IMPLICIT REAL*8(A-H,O-Z)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C:::  SHELL-METZGER SORT IN ASCENDING ORDER.          ...CYRUS 1979  :::
C:::  MODIFIED SLIGHTLY FOR READIBILITY.          ...CYRUS 7 DEC 83  :::
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      DIMENSION D(N)
      LOGNB2=INT(DLOG(DFLOAT(N))/DLOG(2.D0)+1.D-14)
      M=N
      DO 20 NN=1,LOGNB2
        M=M/2
        K=N-M
        DO 20 J=1,K
          DO 10 I=J,1,-M
            L=I+M
            IF (D(L).GT.D(I)) GOTO 20
            T=D(I)
            D(I)=D(L)
            D(L)=T
   10       CONTINUE
   20     CONTINUE
      RETURN
      END
