      module verify_orbitals_mod
      contains
      subroutine verify_orbitals

      use dorb_m,  only: iworbd
      use error,   only: fatal_error
      use optwf_control, only: ioptorb
      use orbval,  only: nadorb,ndetorb,orb
      use slater,  only: ndet
      use system,  only: nelec
      use slater, only: norb

      implicit none

      integer :: i, j


c orbital indices in determinants of trial wave function
      ndetorb=0

      do i=1,ndet
       do j=1,nelec
        if(iworbd(j,i).gt.norb) then
         write (10,*) i,j,iworbd(j,i),norb
         call fatal_error('VERIFY: orbital index out of range')
        endif
        if(iworbd(j,i).gt.ndetorb) then
         ndetorb=iworbd(j,i)
        endif
       enddo
      enddo
 10   format('Det ',i4,' column ',i4,' orb index ',i4,' norb ',i4)

      if(ioptorb.eq.0) then
        norb=ndetorb
        nadorb=0
      endif

      return
      end
      end module
