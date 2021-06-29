      subroutine detsav(iel,iflag)
c Written by Claudia Filippi

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
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
      use multislatern, only: ddorbn, detn, dorbn, orbn

      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb
      use slater, only: d2dx2, ddx, fp, fpp, slmi

      use multislater, only: detiab
      use precision_kinds, only: dp
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
      do 15 j=1,nel*nel
   15   slmi(j,iab)=slmin(j)
      do 30 j=ivirt(iab),norb
        do 30 i=1,nel
          do 20 istate=1,nstates
   20       ymat(j,i,iab,istate)=ymatn(j,i,istate)
   30   aa(i,j,iab)=aan(i,j)

        do 50 k=1,ndet
          if(k.eq.kref) go to 50
          ndim=numrep_det(k,iab)
          do 40 i=1,ndim*ndim
   40         wfmat(i,k,iab)=wfmatn(i,k)
   50   continue

        do 60 j=1,nel
          fp(1,j+ikel,iab)=dorbn(1,iworbd(j+ish,kref))
          fp(2,j+ikel,iab)=dorbn(2,iworbd(j+ish,kref))
   60     fp(3,j+ikel,iab)=dorbn(3,iworbd(j+ish,kref))
        do 70 k=1,ndet
   70     detiab(k,iab)=detn(k)

         do 80 iorb=1,norb
           orb(iel,iorb)=orbn(iorb)
           do 80 kk=1,3
   80        dorb(kk,iel,iorb)=dorbn(kk,iorb)

      return
      end
