      module gauss_mod
      contains
      function gauss()
c Written by Cyrus Umrigar
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Generate a gaussian random number using Box-Mueller method.
c Rewrote inline function as a function so that there would be no ambiquity as to
c the order in which the 2 rannyu's are evaluated.
c Could generate 2 numbers for almost the same price, but for
c backward compatibility generate just 1.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      use const, only: pi
      use precision_kinds, only: dp
      implicit none

      interface
         function rannyu(idum)
          use precision_kinds, only: dp
         implicit none
         integer,intent(in) :: idum
         real(dp) :: rannyu
         end function rannyu
      end interface

      real(dp) :: gauss

      gauss=dcos(2*pi*rannyu(0))
      gauss=gauss*sqrt(-2*dlog(rannyu(0)))

      return
      end
      end module
