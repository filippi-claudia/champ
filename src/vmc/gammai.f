      module gammai_mod
      contains
      function gammai(a,x,xae,iflag)
c Written by Cyrus Umrigar
c Calculate the incomplete gamma function

c To speed up the evaluation several changes have been made.
c 1) The series is used for x < a+4 rather that x < a+1
c 2) Some of the abs() have been removed. OK for a>0
c 3) xae is passed in rather than being re-evaluated here.
c    xae = x**a * exp(-x)
c 4) Warning: To speed up the do loops, the test has been removed
c and the number of times the loop is performed is fixed by formulae
c that are appropriate for a=3/2 and 5/2, but not for general a.
c The formulae are:
c a) n=21+3*(x-a) for series, b) max(20+a-x,7) for cont frac
c--------------------------
c   x     5/2        3/2
c      need form  need form
c--------------------------
c  .01   6   14     6   17
c  .5   12   15    13   18
c 1.5   17   18    18   21
c 3.5   24   24    25   27
c 4.5   27   27    28   30
c 5.5   30   30    31   33
c 6.5   32   33
c-------------------------
c 5.5              16   16
c 7.5   13   15    13   14
c 9.5   11   13    12   12
c12.5    9   10    10    9
c25.5    7    7     7    7

c We have solved the problem of round-off as follows:
c The problem occurs for very large x. If both branches of this routine
c were changed to return the complementary fn gammac, then the problem
c would occur for very small x.  The problem is that for
c x large gammai is calculated by calculating its complement gammac and
c subtracting it from gamma.  However for x very large, gammai is almost
c zero and for x > 43 all significant digits are lost when gammac is
c subtracted from gamma. Since the expression for the area in metrop
c involves the difference of gammai's, the gamma's cancel giving area=0.
c The solution is to return gammai if x < a+4 and gammac otherwise.
c This requires changing all the metrop routines that call gammai.
c Returns   gammai, iflag=1,    x <  a+4
c          -gammac, iflag=-1,   x >= a+4 even though it is called gammai

      use contrl_file, only: ounit
      use precision_kinds, only: dp
      implicit none

      integer :: iflag, iter, itmax, large, n
      real(dp) :: a, a0, a1, an, ana
      real(dp) :: anf, ap, b0, b1
      real(dp) :: cc, ci, d3b2, d5b2
      real(dp) :: del, eps, fac, g
      real(dp) :: g3b2, g5b2, gamm, gammai, gammcf
      real(dp) :: gold, sum
      real(dp) :: x, xae

c     parameter (itmax=100,eps=1.d-14, d3b2=1.5d0,d5b2=2.5d0,
c    & g3b2=.886226925452758d0,g5b2=1.329340388179137d0)

      if(x.lt.a+4)then
c       write(ounit,*)'series'
        iflag=1
        ap=a
        sum=1/a
        del=sum
        iter=int(21+3*(x-a))
ci      do 10 n=1,itmax
        do n=1,iter
          ap=ap+1
          del=del*x/ap
          sum=sum+del
ci        if(del.lt.sum*eps)go to 20
        enddo
ci      call fatal_error ('a too large, itmax too small')
ci 20   gammai=sum*dexp(-x+a*dlog(x))
        gammai=sum*xae
       else
c       write(ounit,*) 'contin frac'
        iflag=-1
c       gold=0
        a0=1
        a1=x
        b0=0
        b1=1
        fac=1
        iter=max(20+nint(a-x),7)
ci      do 30 n=1,itmax
        do n=1,iter
          an=dble(n)
          ana=an-a
          a0=(a1+a0*ana)*fac
          b0=(b1+b0*ana)*fac
          anf=an*fac
          a1=x*a0+anf*a1
          b1=x*b0+anf*b1
            fac=1/a1
            g=b1*fac
c           if(abs((g-gold)/g).lt.eps)go to 40
ci          if(abs((g-gold)).lt.eps*g)go to 40
c           gold=g
        enddo
ci      call fatal_error ('a too large, itmax too small')
ci 40   gammcf=g*dexp(-x+a*dlog(x))
        gammcf=g*xae
        gammai=-gammcf
cc      if(a.eq.d3b2) then
cc        gammai=g3b2-gammcf
cc       elseif(a.eq.d5b2) then
cc        gammai=g5b2-gammcf
cc       else
cc        gammai=gamm(a)-gammcf
cc      endif
      endif
c     write(ounit,*) xae,dexp(-x+a*dlog(x))
      return
      end
c-----------------------------------------------------------------------

      function gamm(a)
c Written by Cyrus Umrigar
c Evaluate the Gamma function.
c For a=3/2 the error of this formula is 5.d-11
      use precision_kinds, only: dp
      implicit none

      integer :: j
      real(dp) :: a, fpf, gamm, half, one
      real(dp) :: ser, stp, tmp, x
      real(dp), dimension(6) :: cof
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     * -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data half,one,fpf/.5d0,1.d0,5.5d0/
      x=a-one
      tmp=x+fpf
c     tmp=(x+half)*dlog(tmp)-tmp
      tmp=tmp**(x+half)*dexp(-tmp)
      ser=one
      do j=1,6
        x=x+one
        ser=ser+cof(j)/x
      enddo
c     gammln=tmp+dlog(stp*ser)
c     gamm=dexp(gammln)
      gamm=stp*ser*tmp
      return
      end
      end module
