      subroutine jastrowe(iel,x,v,d2,value,iflag)
c Written by Claudia Filippi by modifying jastrow

      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use const, only: nelec
      implicit real*8(a-h,o-z)




      parameter (zero=0.d0)


      include 'pseudo.h'

      include 'force.h'


      dimension x(3,*),v(3,*)

      do 10 i=1,nelec
        v(1,i)=zero
        v(2,i)=zero
   10   v(3,i)=zero

      call jastrow4e(iel,x,v,d2,value,iflag)

      return
      end
