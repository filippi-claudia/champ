      subroutine verify_orbitals
      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use dets, only: cdet, ndet
      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      use coefs, only: coef, nbasis, norb
      use contrl, only: idump, irstar, isite, n_conf, nblk, nblkeq, nconf_new, nstep
      use dorb_m, only: iworbd
      use basis, only: zex, betaq, n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz,
     & n4s, n4p, n4fxxx, n4fyyy, n4fzzz, n4fxxy, n4fxxz, n4fyyx, n4fyyz,
     & n4fzzx, n4fzzy, n4fxyz, nsa, npa, ndzra, ndz2a, ndxya, ndxza, ndyza

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'embed.h'
      include 'optorb.h'
      include 'optci.h'
      include 'mstates.h'
      include 'inputflags.h'

      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb


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
