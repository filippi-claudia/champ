      function bode(y,h,n)
      implicit real*8(a-h,o-z)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c:: Bode integration routine wts (14, 64, 24, 64, 14)/45 repeated     ::
c:: Error is (y'''''')*(h**6)/90                                      ::
c:: Assumes that number of pts is 4*n+1                               ::
c::                                           Cyrus Feb 92            ::
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      parameter(oneb45=1.d0/45.d0)
      parameter(c0=28.d0,c1=64.d0,c2=24.d0,half=.5d0)
      dimension y(*)

      if(mod(n-1,4).ne.0) stop 'n must be 4*n+1 in bode'

      n1=n-1
      s=half*c0*(y(n)-y(1))
      do 10 i=1,n1,4
   10   s=s+c0*y(i)+c1*(y(i+1)+y(i+3))+c2*y(i+2)
      bode=s*h*oneb45
      return
      end
