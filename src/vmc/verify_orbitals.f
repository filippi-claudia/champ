      subroutine verify_orbitals
      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use dets, only: cdet, ndet
      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'force.h'
      include 'basis.h'
      include 'embed.h'
      include 'optorb.h'
      include 'optci.h'
      include 'mstates.h'
      include 'inputflags.h'

      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump
     &     ,irstar
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb
      common /dorb/ iworbd(MELEC,MDET)

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

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

      call p2gtid('optwf:ioptorb',ioptorb,0,1)

      if(ioptorb.eq.0) then
        norb=ndetorb
        nadorb=0
      endif

      return
      end
