module gauss_mod
contains
      function gauss()
!> @author Cyrus Umrigar
!>
!> Generate a gaussian random number using Box-Mueller method.
!> Rewrote inline function as a function so that there would be no ambiquity as to
!> the order in which the 2 rannyu's are evaluated.
!> Could generate 2 numbers for almost the same price, but for
!> backward compatibility generate just 1.
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
