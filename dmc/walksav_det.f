      subroutine walksav_det(iw)
c Written by Claudia Filippi

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'
      include 'mstates.h'

      parameter (MEXCIT=10)

      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /forcest/ fgcum(MFORCE),fgcm2(MFORCE)
      common /force_dmc/ itausec,nwprod

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /elec/ nup,ndn
      common /dorb/ iworbd(MELEC,MDET)

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb

      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      common /multidet/ kref,numrep_det(MDET,2),irepcol_det(MELEC,MDET,2),ireporb_det(MELEC,MDET,2)
     & ,iwundet(MDET,2),iactv(2),ivirt(2)

      common /slater/ slmui(MMAT_DIM),slmdi(MMAT_DIM)
     &,fpu(3,MMAT_DIM),fpd(3,MMAT_DIM)
     &,fppu(MMAT_DIM),fppd(MMAT_DIM)
     &,ddx(3,MELEC),d2dx2(MELEC)
      common /multislater/ detu(MDET),detd(MDET)

      common /multimat/ aa(MELEC,MORB,2),wfmat(MEXCIT**2,MDET,2)

      common /ycompact/ ymat(MORB,MELEC,2,MSTATES),dymat(MORB,MELEC,2,MSTATES)

      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb

      dimension krefw(MWALK),slmuiw(MMAT_DIM,MWALK)
     &,slmdiw(MMAT_DIM,MWALK)
     &,fpuw(3,MMAT_DIM,MWALK),fpdw(3,MMAT_DIM,MWALK)
     &,fppuw(MMAT_DIM,MWALK),fppdw(MMAT_DIM,MWALK)
     &,ddxw(3,MELEC,MWALK),d2dx2w(MELEC,MWALK)
     &,detuw(MDET,MWALK),detdw(MDET,MWALK)

      dimension aaw(MELEC,MORB,MWALK,2),wfmatw(MEXCIT**2,MDET,MWALK,2),ymatw(MORB,MELEC,MWALK,2,MSTATES)

      dimension orbw(MELEC,MORB,MWALK),dorbw(3,MELEC,MORB,MWALK)

      save krefw,slmuiw,slmdiw,fpuw,fpdw,fppuw,fppdw,detuw,detdw,ddxw,d2dx2w

      save aaw,wfmatw,ymatw,orbw,dorbw

       do 20 k=1,ndet
         detuw(k,iw)=detu(k)
   20    detdw(k,iw)=detd(k)

       krefw(iw)=kref
       do 40 j=1,nup*nup
         slmuiw(j,iw)=slmui(j)
         fpuw(1,j,iw)=fpu(1,j)
         fpuw(2,j,iw)=fpu(2,j)
   40    fpuw(3,j,iw)=fpu(3,j)
       do 50 j=1,ndn*ndn
         slmdiw(j,iw)=slmdi(j)
         fpdw(1,j,iw)=fpd(1,j)
         fpdw(2,j,iw)=fpd(2,j)
   50    fpdw(3,j,iw)=fpd(3,j)
       do 55 i=1,nelec
         ddxw(1,i,iw)=ddx(1,i)
         ddxw(2,i,iw)=ddx(2,i)
   55    ddxw(3,i,iw)=ddx(3,i)

       do 62 iab=1,2
         nel=nup
         if(iab.eq.2) nel=ndn
         do 58 j=ivirt(iab),norb
          do 58 i=1,nel
            do 57 istate=1,nstates
   57         ymatw(j,i,iw,iab,istate)=ymat(j,i,iab,istate)
   58       aaw(i,j,iw,iab)=aa(i,j,iab)
          do 62 k=1,ndet
            if(k.ne.kref) then
              ndim=numrep_det(k,iab)
              do 60 i=1,ndim*ndim
   60           wfmatw(i,k,iw,iab)=wfmat(i,k,iab)
            endif
   62  continue

       do 63 i=1,nelec
         do 63 iorb=1,norb
           orbw(i,iorb,iw)=orb(i,iorb)
           do 63 kk=1,3
   63        dorbw(kk,i,iorb,iw)=dorb(kk,i,iorb)

      return

      entry walkstrdet(iw)

      do 70 k=1,ndet
        detu(k)=detuw(k,iw)
   70   detd(k)=detdw(k,iw)

      kref=krefw(iw)
      do 80 j=1,nup*nup
        slmui(j)=slmuiw(j,iw)
        fpu(1,j)=fpuw(1,j,iw)
        fpu(2,j)=fpuw(2,j,iw)
   80   fpu(3,j)=fpuw(3,j,iw)
      do 90 j=1,ndn*ndn
        slmdi(j)=slmdiw(j,iw)
        fpd(1,j)=fpdw(1,j,iw)
        fpd(2,j)=fpdw(2,j,iw)
   90   fpd(3,j)=fpdw(3,j,iw)
      do 95 i=1,nelec
        ddx(1,i)=ddxw(1,i,iw)
        ddx(2,i)=ddxw(2,i,iw)
   95   ddx(3,i)=ddxw(3,i,iw)

       do 102 iab=1,2
         nel=nup
         if(iab.eq.2) nel=ndn
         do 98 j=ivirt(iab),norb
          do 98 i=1,nel
            do 97 istate=1,nstates
   97         ymat(j,i,iab,istate)=ymatw(j,i,iw,iab,istate)
   98       aa(i,j,iab)=aaw(i,j,iw,iab)
          do 102 k=1,ndet
            if(k.ne.kref) then
              ndim=numrep_det(k,iab)
              do 100 i=1,ndim*ndim
  100           wfmat(i,k,iab)=wfmatw(i,k,iw,iab)
            endif
  102  continue

       do 103 i=1,nelec
         do 103 iorb=1,norb
           orb(i,iorb)=orbw(i,iorb,iw)
           do 103 kk=1,3
  103        dorb(kk,i,iorb)=dorbw(kk,i,iorb,iw)

      return

      entry splitjdet(iw,iw2)

      do 110 k=1,ndet
        detuw(k,iw2)=detuw(k,iw)
  110   detdw(k,iw2)=detdw(k,iw)

      krefw(iw2)=krefw(iw)
      do 120 j=1,nup*nup
        slmuiw(j,iw2)=slmuiw(j,iw)
        fpuw(1,j,iw2)=fpuw(1,j,iw)
        fpuw(2,j,iw2)=fpuw(2,j,iw)
  120   fpuw(3,j,iw2)=fpuw(3,j,iw)
      do 130 j=1,ndn*ndn
        slmdiw(j,iw2)=slmdiw(j,iw)
        fpdw(1,j,iw2)=fpdw(1,j,iw)
        fpdw(2,j,iw2)=fpdw(2,j,iw)
  130   fpdw(3,j,iw2)=fpdw(3,j,iw)
      do 135 i=1,nelec
        ddxw(1,i,iw2)=ddxw(1,i,iw)
        ddxw(2,i,iw2)=ddxw(2,i,iw)
  135   ddxw(3,i,iw2)=ddxw(3,i,iw)

       do 142 iab=1,2
         nel=nup
         if(iab.eq.2) nel=ndn
         do 138 j=ivirt(iab),norb
          do 138 i=1,nel
            do 137 istate=1,nstates
  137         ymatw(j,i,iw2,iab,istate)=ymatw(j,i,iw,iab,istate)
  138       aaw(i,j,iw2,iab)=aaw(i,j,iw,iab)
          do 142 k=1,ndet
            if(k.ne.kref) then
              ndim=numrep_det(k,iab)
              do 140 i=1,ndim*ndim
  140           wfmatw(i,k,iw2,iab)=wfmatw(i,k,iw,iab)
            endif
  142  continue

       do 143 i=1,nelec
         do 143 iorb=1,norb
           orbw(i,iorb,iw2)=orbw(i,iorb,iw)
           do 143 kk=1,3
  143        dorbw(kk,i,iorb,iw2)=dorbw(kk,i,iorb,iw)

      return
      end
