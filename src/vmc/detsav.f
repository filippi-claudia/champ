      module detsav_mod
      contains
      subroutine detsav(iel,iflag)
c Written by Claudia Filippi

      use csfs,    only: nstates
      use dorb_m,  only: iworbd
      use multidet, only: ivirt,ndetiab,ndetsingle,numrep_det, ndetdouble
      use multimat, only: aa,wfmat
      use multimatn, only: aan,wfmatn
      use multislater, only: detiab
      use multislatern, only: detn,dorbn,orbn
      use orbval,  only: dorb,orb
      use precision_kinds, only: dp
      use slater,  only: fp,kref,ndet,norb,slmi
      use slatn,   only: slmin
      use system,  only: ndn,nelec,nup
      use vmc_mod, only: MEXCIT
      use ycompact, only: ymat
      use ycompactn, only: ymatn

      implicit none

      integer :: i, iab, iel, iflag, ikel
      integer :: iorb, ish, istate, j
      integer :: k, kk, ndim, nel, ndim2, kn,kcum
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
      do j=1,nel*nel
         slmi(j,iab)=slmin(j)
      enddo
      do j=ivirt(iab),norb
         do i=1,nel
            do istate=1,nstates
               ymat(j,i,iab,istate)=ymatn(j,i,istate)
            enddo
            aa(i,j,iab)=aan(i+nelec*(j-1))
         enddo
      enddo
      

!     This loop should run just over unique or unequivalent determinants
!     single excitations
      if(ndetsingle(iab).ge.1)then
         do k=1,ndetsingle(iab)
            wfmat(k,1,iab)=wfmatn(k,1)
         enddo
      endif

!     double excitations       
      kcum=ndetsingle(iab)+ndetdouble(iab)
      if(ndetdouble(iab).ge.1)then
         do k=ndetsingle(iab)+1,kcum
            wfmat(k,1:4,iab)=wfmatn(k,1:4)
         enddo
      endif
      
!     multiple excitations
      if(kcum.lt.ndetiab(iab))then
         do k=ndetdouble(iab)+1,ndetiab(iab)
            ndim=numrep_det(k,iab)
            ndim2=ndim*ndim
            wfmat(k,1:ndim2,iab)=wfmatn(k,1:ndim2)
         enddo
      endif
      
      
      
      do j=1,nel
         fp(1,j+ikel,iab)=dorbn(iworbd(j+ish,kref),1)
         fp(2,j+ikel,iab)=dorbn(iworbd(j+ish,kref),2)
         fp(3,j+ikel,iab)=dorbn(iworbd(j+ish,kref),3)
      enddo
      do k=1,ndet
         detiab(k,iab)=detn(k)
      enddo

      
      do iorb=1,norb
         orb(iel,iorb)=orbn(iorb)
         dorb(iorb,iel,:)=dorbn(iorb,:)
      enddo
      
      
      return
      end
      end module
