      subroutine sites(x,nelec,nsite)
c Written by Cyrus Umrigar
      use atom, only: znuc, cent, iwctype, ncent
      use contrl_file,    only: ounit
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

      integer :: i, ic, ispin, j, ju
      integer :: l, nelec
      integer, dimension(*) :: nsite
      real(dp) :: site, sitsca
      real(dp), dimension(3,*) :: x
      real(dp), parameter :: half = 0.5d0


c routine to put electrons down around centers for an initial
c configuration if nothing else is available


c loop over spins and centers. If odd number of electrons on all
c atoms then the up-spins have an additional electron.

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

c sample position from exponentials around center

             site=-dlog(rannyu(0))
             site=sign(site,(rannyu(0)-half))
             x(ic,l)=sitsca*site+cent(ic,i)
            enddo
          enddo
        enddo
      enddo
      write(ounit,'(''number of electrons placed ='',i5)') l
      if (l.lt.nelec) call fatal_error('SITES: bad input')
      return
      end
