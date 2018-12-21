      subroutine multideterminant_tmove(psid,iel_move)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'pseudo.h'
      include 'mstates.h'

      parameter (one=1.d0,half=0.5d0)
      parameter (MEXCIT=10)

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /elec/ nup,ndn
      common /dorb/ iworbd(MELEC,MDET)

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,lpot(MCTYPE),nloc
      common /qua/ xq0(MPS_QUAD),yq0(MPS_QUAD),zq0(MPS_QUAD)
     &,xq(MPS_QUAD),yq(MPS_QUAD),zq(MPS_QUAD),wq(MPS_QUAD),nquad

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb

      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      common /multidet/ kref,numrep_det(MDET,2),irepcol_det(MELEC,MDET,2),ireporb_det(MELEC,MDET,2)
     & ,iwundet(MDET,2),iactv(2),ivirt(2)

      common /slater/ slmi(MMAT_DIM,2)
     &,fpu(3,MMAT_DIM),fpd(3,MMAT_DIM)
     &,fppu(MMAT_DIM),fppd(MMAT_DIM)
     &,ddx(3,MELEC),d2dx2(MELEC)
      common /multislater/ detu(MDET),detd(MDET)

      common /multimat/ aa(MELEC,MORB,2),wfmat(MEXCIT**2,MDET,2)

      common /ycompact/ ymat(MORB,MELEC,2,MSTATES),dymat(MORB,MELEC,2,MSTATES)

      common /casula/ t_vpsp(MCENT,MPS_QUAD,MELEC),icasula,i_vpsp
      common /b_tmove/ b_t(MORB,MPS_QUAD,MCENT,MELEC),iskip(MELEC,MCENT)

      dimension gmat(MELEC,MORB)

      if(icasula.gt.0)then
        i1=iel_move
        i2=iel_move
       else
        i1=1
        i2=nelec
      endif

      do iel=i1,i2

      do ic=1,ncent
        
      if(iskip(iel,ic).eq.0) then

      if(iel.le.nup) then
        iab=1
        nel=nup
        ish=0
       else
        iab=2
        nel=ndn
        ish=nup
      endif

      detratio=detu(kref)*detd(kref)/psid

      jel=iel-ish

      do iq=1,nquad

        do jrep=ivirt(iab),norb
          dum=0
          do j=1,nel
            dum=dum+b_t(iworbd(j+ish,kref),iq,ic,iel)*aa(j,jrep,iab)
          enddo
          dum=b_t(jrep,iq,ic,iel)-dum

          do irep=iactv(iab),nel
            gmat(irep,jrep)=dum*slmi(irep+(jel-1)*nel,iab)
          enddo
        enddo

c     t_vpsp(ic,iq,iel)=t_vpsp_ref

      dum=0
      do jrep=ivirt(iab),norb
        do irep=iactv(iab),nel
          dum=dum+ymat(jrep,irep,iab,1)*gmat(irep,jrep)
        enddo
      enddo
      t_vpsp(ic,iq,iel)=t_vpsp(ic,iq,iel)+dum*detratio

      enddo

      endif

      enddo
      enddo

      return
      end
