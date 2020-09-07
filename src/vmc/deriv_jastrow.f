      subroutine deriv_jastrow(x,v,d2,div_vj,value)
c Written by Claudia Filippi

      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use const, only: nelec

      implicit real*8(a-h,o-z)

      parameter (zero=0.d0)
      include 'force.h'
      include 'pseudo.h'
      dimension x(3,*),v(3,*),div_vj(MELEC)

      do 10 i=1,nelec
        v(1,i)=zero
        v(2,i)=zero
   10   v(3,i)=zero
      d2=zero

      call deriv_jastrow4(x,v,d2,value)

      return
      end
