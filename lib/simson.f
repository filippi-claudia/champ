      subroutine simson(y,s,h,n)
      implicit real*8(a-h,o-z)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c:: Simpson integration routine wts 1/3, 4/3, 2/3, 4/3 ... 4/3, 1/3   ::
c:: Uses Newton's  4 pt integ. for the end pts        if n is even    ::
c:: Wts for this are 3/8, 9/8, 9/8, 3/8                               ::
c:: Hence if n is even correction for wts is:                         ::
c:: c3=(1/3+3/8-2/3), c2=(9/8-4/3), c1=(9/8-2/3), c0=(3/8-4/3)        ::
c:: ( =1/8=.1125,       =-5/8=-.625   =11/8=1.375   =-23/8=-2.875) /3 ::
c:: Error is (y'''')*(h**4)/90                        if n is odd     ::
c::   "   "          "        +(y'''')*(h**5)*(3/80)  if n is even    ::
c:: If the function is not analytic, the error seems to go as h**1.5, ::
c:: atleast for Log it does.                                          ::
c:: For non-analytic function one should probably use Teter's grater. ::
c::                                        ...Cyrus  Aug 82           ::
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      dimension y(*)
c...  c3=1/8, c2=-5/8, c1=11/8, c0=-23/8 
      data c3,c2,c1,c0/0.125d0,-0.625d0,1.375d0,-2.875d0/
      data zero,two,thrd/0.d0,2.d0,0.333333333333333d0/
      h3=h*thrd
      n1=n-1
      odd=zero
      eve=zero
      do 10 i=1,n1,2
        odd=odd+y(i)
   10   eve=eve+y(i+1)
      s=two*(odd+two*eve) - y(1)
      if(mod(n,2).eq.1) then
          s=(s+y(n)) * h3
        else
          s=(s+c3*y(n-3)+c2*y(n-2)+c1*y(n-1)+c0*y(n)) * h3
      endif
      return
      end
