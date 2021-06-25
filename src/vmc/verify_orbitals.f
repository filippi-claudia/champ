      subroutine verify_orbitals

      use vmc_mod, only: MELEC, MORB
      use const, only: nelec
      use dets, only: ndet
      use optwf_contrl, only: ioptorb
      use coefs, only: norb
      use dorb_m, only: iworbd
      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb

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

