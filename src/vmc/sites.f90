module sites_mod
      use error,   only: fatal_error
contains
      subroutine sites(x,nelec,nsite)
! Written by Cyrus Umrigar
      use contrl_file, only: ounit
      use precision_kinds, only: dp
      use random_mod, only: random_dp
      use system,  only: cent,iwctype,ncent,znuc
      implicit none

      integer :: i, ic, ispin, j, ju
      integer :: l, nelec
      integer, dimension(*) :: nsite
      real(dp) :: site, sitsca
      real(dp), dimension(3,*) :: x
      real(dp), parameter :: half = 0.5d0


! routine to put electrons down around centers for an initial
! configuration if nothing else is available


! loop over spins and centers. If odd number of electrons on all
! atoms then the up-spins have an additional electron.

      l=0
      do ispin=1,2
        do i=1,ncent
          ju=(nsite(i)+2-ispin)/2
          do j=1,ju
            l=l+1
            if (l.gt.nelec) return
            if(j.eq.1) then
              sitsca=1/znuc(iwctype(i))
             elseif(j.le.5) then
              sitsca=2/(znuc(iwctype(i))-2)
             else
              sitsca=3/(znuc(iwctype(i))-10)
            endif
            do ic=1,3

! sample position from exponentials around center

             site=-dlog(random_dp())
             site=sign(site,(random_dp()-half))
             x(ic,l)=sitsca*site+cent(ic,i)
            enddo
          enddo
        enddo
      enddo
      write(ounit,'(''number of electrons placed ='',i5)') l
      if (l.lt.nelec) call fatal_error('SITES: bad input')
      return
      end
end module
