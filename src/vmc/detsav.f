      subroutine detsav(iel,iflag)
c Written by Claudia Filippi

      use csfs, only: nstates

      use dets, only: ndet
      use elec, only: ndn, nup
      use multidet, only: ivirt, kref, numrep_det

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
      implicit none

      integer :: i, iab, iel, iflag, ikel
      integer :: iorb, ish, istate, j
      integer :: k, kk, ndim, nel








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
        aa(i,j,iab)=aan(i,j)
        enddo
      enddo

      do k=1,kref-1
     
         ndim=numrep_det(k,iab)
          do i=1,ndim*ndim
             wfmat(i,k,iab)=wfmatn(i,k)
          enddo

      enddo


      do k=kref+1,ndet
         
         ndim=numrep_det(k,iab)
         do i=1,ndim*ndim
            wfmat(i,k,iab)=wfmatn(i,k)
         enddo
         
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
           do kk=1,3
             dorb(kk,iel,iorb)=dorbn(iorb,kk)
           enddo
         enddo

      return
      end
