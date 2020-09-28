      subroutine jassav(iel,iflag)
c Written by Claudia Filippi

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use const, only: nelec
      use jaso, only: d2ijo, d2o, fijo, fjo, fso, fsumo

      use velocity_jastrow, only: vj, vjn
      use jasn, only: d2ijn, d2n, fijn, fjn, fsn, fsumn

      implicit real*8(a-h,o-z)









      fsumo=fsumn
      do 10 i=1,nelec
        fjo(1,i)=fjn(1,i)
        fjo(2,i)=fjn(2,i)
  10    fjo(3,i)=fjn(3,i)

      do 20 j=1,iel
  20    fso(iel,j)=fsn(iel,j)

      do 30 j=iel+1,nelec
  30    fso(j,iel)=fsn(j,iel)

      do 40 j=1,nelec
        fijo(1,iel,j)=fijn(1,iel,j)
        fijo(2,iel,j)=fijn(2,iel,j)
        fijo(3,iel,j)=fijn(3,iel,j)
        fijo(1,j,iel)=fijn(1,j,iel)
        fijo(2,j,iel)=fijn(2,j,iel)
  40    fijo(3,j,iel)=fijn(3,j,iel)

      if(iflag.gt.0) then
        d2o=d2n
        do 50 j=1,iel
  50      d2ijo(iel,j)=d2ijn(iel,j)

        do 60 j=iel+1,nelec
  60      d2ijo(j,iel)=d2ijn(j,iel)
      endif

      do 100 i=1,nelec
        vj(1,i)=vjn(1,i)
        vj(2,i)=vjn(2,i)
 100    vj(3,i)=vjn(3,i)
      return
      end
