      subroutine detsav(iel,iflag)
c Written by Claudia Filippi

      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use dets, only: cdet, ndet
      use elec, only: ndn, nup
      use multidet, only: iactv, irepcol_det, ireporb_det, ivirt, iwundet, kref, numrep_det

      use slatn, only: slmin
      use ycompact, only: dymat, ymat
      use ycompactn, only: ymatn
      use coefs, only: coef, nbasis, norb
      use dorb_m, only: iworbd
      use multimat, only: aa, wfmat
      use multimatn, only: aan, wfmatn
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      common /slater/ slmui(MMAT_DIM),slmdi(MMAT_DIM)
     &,fpu(3,MMAT_DIM),fpd(3,MMAT_DIM)
     &,fppu(MMAT_DIM),fppd(MMAT_DIM)
     &,ddx(3,MELEC),d2dx2(MELEC)
      common /multislater/ detu(MDET),detd(MDET)
      common /multislatern/ detn(MDET)
     &,orb(MORB),dorb(3,MORB),ddorb(MORB)
      common /orbval/ orbo(MELEC,MORB),dorbo(3,MELEC,MORB),ddorbo(MELEC,MORB),ndetorb,nadorb

      if(iel.le.nup) then
        ikel=nup*(iel-1)
        do 35 l=1,nup*nup
   35     slmui(l)=slmin(l)
       do 46 j=ivirt(1),norb
        do 46 i=1,nup
          do 40 istate=1,nstates
   40       ymat(j,i,1,istate)=ymatn(j,i,istate)
   46     aa(i,j,1)=aan(i,j)
        do 48 k=1,ndet
          if(k.eq.kref) go to 48
          ndim=numrep_det(k,1)
          do 47 i=1,ndim*ndim
   47         wfmat(i,k,1)=wfmatn(i,k)
   48   continue

        do 50 j=1,nup
          fpu(1,j+ikel)=dorb(1,iworbd(j,kref))
          fpu(2,j+ikel)=dorb(2,iworbd(j,kref))
   50     fpu(3,j+ikel)=dorb(3,iworbd(j,kref))
        do 60 k=1,ndet
   60     detu(k)=detn(k)

         do 63 iorb=1,norb
           orbo(iel,iorb)=orb(iorb)
           do 63 kk=1,3
   63        dorbo(kk,iel,iorb)=dorb(kk,iorb)
       else
        ikel=ndn*(iel-nup-1)
        do 65 j=1,ndn*ndn
   65     slmdi(j)=slmin(j)
       do 76 j=ivirt(2),norb
        do 76 i=1,ndn
          do 70 istate=1,nstates   
   70       ymat(j,i,2,istate)=ymatn(j,i,istate)
   76     aa(i,j,2)=aan(i,j)
        do 78 k=1,ndet
          if(k.eq.kref) go to 78
          ndim=numrep_det(k,2)
          do 77 i=1,ndim*ndim
   77         wfmat(i,k,2)=wfmatn(i,k)
   78   continue

        do 80 j=1,ndn
          fpd(1,j+ikel)=dorb(1,iworbd(j+nup,kref))
          fpd(2,j+ikel)=dorb(2,iworbd(j+nup,kref))
   80     fpd(3,j+ikel)=dorb(3,iworbd(j+nup,kref))
        do 90 k=1,ndet
   90     detd(k)=detn(k)

         do 93 iorb=1,norb
           orbo(iel,iorb)=orb(iorb)
           do 93 kk=1,3
   93        dorbo(kk,iel,iorb)=dorb(kk,iorb)
      endif

      return
      end
