      function gauss()
c Written by Cyrus Umrigar
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Generate a gaussian random number using Box-Mueller method.
c Rewrote inline function as a function so that there would be no ambiquity as to
c the order in which the 2 rannyu's are evaluated.
c Could generate 2 numbers for almost the same price, but for
c backward compatibility generate just 1.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr

      gauss=dcos(2*pi*rannyu(0))
      gauss=gauss*sqrt(-2*dlog(rannyu(0)))

      return
      end
