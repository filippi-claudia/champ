      module detsav_mod
      contains
      subroutine detsav(iel,iflag)
c Written by Claudia Filippi

      use csfs, only: nstates

      use slater, only: ndet
      use system, only: ndn, nup, nelec
      use multidet, only: ivirt, numrep_det, ndetiab, ndetsingle
      use slater, only: kref
      use slatn, only: slmin
      use ycompact, only: ymat
      use ycompactn, only: ymatn
      use slater, only: norb
      use dorb_m, only: iworbd
      use multimat, only: aa, wfmat
      use multimatn, only: aan, wfmatn
      use multislatern, only: detn, dorbn, orbn
      use orbval, only: dorb, orb
      use slater, only: fp, slmi
      use multislater, only: detiab
      use vmc_mod, only: MEXCIT
      use precision_kinds, only: dp

      implicit none

      integer :: i, iab, iel, iflag, ikel
      integer :: iorb, ish, istate, j
      integer :: k, kk, ndim, nel, ndim2, kn
      integer, dimension(ndet) :: ku
      integer, dimension(ndet) :: auxdim

      
      if(iel.le.nup) then
         iab=1
         nel=nup
         ish=0
      else
         iab=2
         nel=ndn
         ish=nup
      endif
      
      ikel=nel*(iel-ish-1)
      do istate=1,nstates
        do j=1,nel*nel
          slmi(j,iab,istate)=slmin(j,istate) !STU also here
        enddo
        do j=ivirt(iab),norb
          do i=1,nel
            ymat(j,i,iab,istate)=ymatn(j,i,istate)
            aa(i,j,iab,istate)=aan(i,j,istate) !STU check state orb mapping here
          enddo
        enddo
      

!     This loop should run just over unique or unequivalent determinants
! single excitations
        do k=1,ndetsingle(iab)
          wfmat(k,1,iab,istate)=wfmatn(k,1,istate) !STU also here
        enddo
! multiple excitations
        do k=ndetsingle(iab)+1,ndetiab(iab)
          ndim=numrep_det(k,iab)
          ndim2=ndim*ndim 
          wfmat(k,1:ndim2,iab,istate)=wfmatn(k,1:ndim2,istate) !STU also here
        enddo
      
      
      
      
        do j=1,nel !STU also here if all below just orbital dependent, close istate loop above and just loop over nwftypeorb
          fp(1,j+ikel,iab,istate)=dorbn(iworbd(j+ish,kref),1,istate)
          fp(2,j+ikel,iab,istate)=dorbn(iworbd(j+ish,kref),2,istate)
          fp(3,j+ikel,iab,istate)=dorbn(iworbd(j+ish,kref),3,istate)
        enddo
        do k=1,ndet !STU also here
          detiab(k,iab,istate)=detn(k,istate)
        enddo

      
        do iorb=1,norb !STU also here
          orb(iel,iorb,istate)=orbn(iorb,istate)
          dorb(iorb,iel,1,istate)=dorbn(iorb,1,istate)
          dorb(iorb,iel,1,istate)=dorbn(iorb,1,istate)
          dorb(iorb,iel,1,istate)=dorbn(iorb,1,istate)
        enddo
      enddo
      
      
      return
      end
      end module
