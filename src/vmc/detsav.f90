module detsav_mod
contains
      subroutine detsav(iel,iflag)
! Written by Claudia Filippi

      use csfs, only: nstates
      use dorb_m, only: iworbd
      use multidet, only: ivirt, ndetiab, numrep_det, ndetsingle, ndetdouble
      use multimat, only: aa, wfmat
      use multimatn, only: aan, wfmatn
      use multislater, only: detiab
      use multislatern, only: detn, dorbn, orbn
      use orbval, only: dorb, orb
      use precision_kinds, only: dp
      use slater, only: fp, slmi
      use slater, only: kref, ndet, norb
      use slatn, only: slmin
      use system, only: ndn, nup, nelec
      use ycompact, only: ymat
      use ycompactn, only: ymatn
      use vmc_mod, only: stoo

      implicit none

      integer :: i, iab, iel, iflag, ikel
      integer :: iorb, ish, istate, j
      integer :: k, kk, ndim, nel, ndim2, kcum

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
          slmi(j,iab,stoo(istate))=slmin(j,stoo(istate))
        enddo
        do j=ivirt(iab),norb
          do i=1,nel
            ymat(j,i,iab,istate)=ymatn(j,i,istate)
            aa(i,j,iab,stoo(istate))=aan(i,j,stoo(istate))
          enddo
        enddo

! loop over unique or unequivalent determinants
! single excitations
        if(ndetsingle(iab).ge.1)then
          do k=1,ndetsingle(iab)
            wfmat(k,1,iab,stoo(istate))=wfmatn(k,1,stoo(istate))
          enddo
        endif

! double excitations  
        kcum=ndetsingle(iab)+ndetdouble(iab)
        if(ndetdouble(iab).ge.1)then
           do k=ndetsingle(iab)+1,kcum
              wfmat(k,1:4,iab,stoo(istate))=wfmatn(k,1:4,stoo(istate))
           enddo
        endif

! multiple excitations
        if(kcum.lt.ndetiab(iab))then
           do k=ndetdouble(iab)+1,ndetiab(iab)
              ndim=numrep_det(k,iab)
              ndim2=ndim*ndim
              wfmat(k,1:ndim2,iab,stoo(istate))=wfmatn(k,1:ndim2,stoo(istate))
           enddo
        endif

        do j=1,nel
          fp(1,j+ikel,iab,stoo(istate))=dorbn(iworbd(j+ish,kref),1,stoo(istate))
          fp(2,j+ikel,iab,stoo(istate))=dorbn(iworbd(j+ish,kref),2,stoo(istate))
          fp(3,j+ikel,iab,stoo(istate))=dorbn(iworbd(j+ish,kref),3,stoo(istate))
        enddo
        do k=1,ndet
          detiab(k,iab,stoo(istate))=detn(k,stoo(istate))
        enddo

        do iorb=1,norb
          orb(iel,iorb,stoo(istate))=orbn(iorb,stoo(istate))
          dorb(iorb,iel,1,stoo(istate))=dorbn(iorb,1,stoo(istate))
          dorb(iorb,iel,2,stoo(istate))=dorbn(iorb,2,stoo(istate))
          dorb(iorb,iel,3,stoo(istate))=dorbn(iorb,3,stoo(istate))
        enddo
      enddo

      return
      end
end module
