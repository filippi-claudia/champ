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
      use constants, only: pi
      use precision_kinds, only: dp
      use random_mod, only: random_dp
      implicit none

      real(dp) :: gauss

      gauss=dcos(2*pi*random_dp())
      gauss=gauss*sqrt(-2*dlog(random_dp()))

      return
      end
      end module
