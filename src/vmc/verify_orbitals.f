      subroutine verify_orbitals
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use embed, only: MXTABI, MXTAB, MXHTAB, MXXTAB
      use embed, only: MXHGRID, MXXGRID
      use const, only: nelec
      use dets, only: ndet
      use optwf_contrl, only: ioptorb
      use coefs, only: norb
      use dorb_m, only: iworbd

      implicit real*8(a-h,o-z)

      include 'force.h'
      include 'optorb.h'
      include 'optci.h'

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

