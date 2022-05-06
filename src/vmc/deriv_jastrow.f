      subroutine deriv_jastrow(x,v,d2,div_vj,value,istate)
c Written by Claudia Filippi

      use vmc_mod, only: MELEC
      use const, only: nelec

      implicit real*8(a-h,o-z)

      parameter (zero=0.d0)
      dimension x(3,*),v(3,*),div_vj(MELEC)

      do 10 i=1,nelec
        v(1,i)=zero
        v(2,i)=zero
   10   v(3,i)=zero
      d2=zero

      call deriv_jastrow4(x,v,d2,value,istate)

      return
      end
