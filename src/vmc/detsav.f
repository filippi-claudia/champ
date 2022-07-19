      module detsav_mod
      contains
      subroutine detsav(iel,iflag)
c Written by Claudia Filippi

      use csfs, only: nstates

      use dets, only: ndet
      use multidet, only: ivirt, kref, numrep_det, ndetiab, ndetsingle

      use slatn, only: slmin
      use ycompact, only: ymat
      use ycompactn, only: ymatn
      use coefs, only: norb
      use dorb_m, only: iworbd
      use multimat, only: aa, wfmat
      use multimatn, only: aan, wfmatn
      use multislatern, only: detn, dorbn, orbn

      use orbval, only: dorb, orb
      use slater, only: fp, slmi

      use multislater, only: detiab

      use vmc_mod, only: MEXCIT


      use precision_kinds, only: dp
      use system, only: nelec
      use system, only: nup
      use system, only: ndn
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
! single excitations
      do k=1,ndetsingle(iab)
          wfmat(k,1,iab)=wfmatn(k,1)
       enddo
! multiple excitations
      do k=ndetsingle(iab)+1,ndetiab(iab)
          ndim=numrep_det(k,iab)
          ndim2=ndim*ndim
          wfmat(k,1:ndim2,iab)=wfmatn(k,1:ndim2)
       enddo
      
      
      
      
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
