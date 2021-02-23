      subroutine multideterminant_tmove(psid,iel_move)

      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X,
     &NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20,
     &radmax, delri, NEQSX, MTERMS, MCENT3, NCOEF, MEXCIT
      use const, only: delta, deltai, etrial, fbias, hb, imetro, ipr, nelec, pi
      use atom, only: cent, iwctype, ncent, nctype, pecent, znuc
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use mstates_mod, only: MSTATES, MDETCSFX
      use pseudo_mod, only: MPS_L, MPS_QUAD, MPS_GRID, MGAUSS
      use qua, only: nquad, wq, xq, xq0, yq, yq0, zq, zq0
      use b_tmove, only: b_t, iskip
      use casula, only: i_vpsp, icasula, t_vpsp
      use pseudo, only: lpot, nloc, vps, vpso
      use slater, only: d2dx2, ddx, fpd, fppd, fppu, fpu, slmi, slmui, slmdi
      use dets, only: cdet, ndet
      use elec, only: ndn, nup

      implicit real*8(a-h,o-z)



      parameter (one=1.d0,half=0.5d0)

      common /dorb/ iworbd(MELEC,MDET)
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /wfsec/ iwftype(MFORCE),iwf,nwftype
      common /multidet/ kref,numrep_det(MDET,2),irepcol_det(MELEC,MDET,2),ireporb_det(MELEC,MDET,2)
     & ,iwundet(MDET,2),iactv(2),ivirt(2)
      common /multislater/ detu(MDET),detd(MDET)
      common /multimat/ aa(MELEC,MORB,2),wfmat(MEXCIT**2,MDET,2)
      common /ycompact/ ymat(MORB,MELEC,2,MSTATES),dymat(MORB,MELEC,2,MSTATES)


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
