      subroutine optorb_deriv(psid,denergy,zmat,dzmat,emz,aaz,orbprim,eorbprim)

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use dets, only: cdet, ndet
      use elec, only: ndn, nup
      use multidet, only: iactv, irepcol_det, ireporb_det, ivirt, iwundet, kref, numrep_det

      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      implicit real*8(a-h,o-z)






      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'
      include 'optorb_cblk.h'

      parameter (MEXCIT=10)

      common /slater/ slmui(MMAT_DIM),slmdi(MMAT_DIM)
     &,fpu(3,MMAT_DIM),fpd(3,MMAT_DIM)
     &,fppu(MMAT_DIM),fppd(MMAT_DIM)
     &,ddx(3,MELEC),d2dx2(MELEC)
      common /multislater/ detu(MDET),detd(MDET)

      common /dorb/ iworbd(MELEC,MDET)

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb


      common /multimat/ aa(MELEC,MORB,2),wfmat(MEXCIT**2,MDET,2)


      common /Bloc/ b(MORB,MELEC),xmatu(MELEC**2),xmatd(MELEC**2)
     & ,tildem(MELEC,MORB,2)


      dimension zmat(MORB,MELEC,2),dzmat(MORB,MELEC,2),emz(MELEC,MELEC,2),aaz(MELEC,MELEC,2)
      dimension orbprim(*),eorbprim(*)

      if(ioptorb.eq.0) return

c     ns_current=ns_current+1
c     if(ns_current.ne.iorbsample) return
c ns_current reset in optorb_sum

      detratio=detu(kref)*detd(kref)/psid
      do 200 iterm=1,norbterm

        io=ideriv(1,iterm)
        jo=ideriv(2,iterm)

        dorb_psi_ref=0
        dorb_energy_ref=0.d0

        dorb_psi=0.d0
        dorb_energy=0.d0
        do iab=1,2

          if(iab.eq.1) then
            ish=0
            nel=nup
           else
            ish=nup
            nel=ndn
          endif

          if(io.ge.ivirt(iab)) then
            do i=1,nel
              dorb_psi=dorb_psi+zmat(io,i,iab)*orb(i+ish,jo)
              dorb_energy=dorb_energy+dzmat(io,i,iab)*orb(i+ish,jo)+zmat(io,i,iab)*b(jo,i+ish)
            enddo
          endif
          if(ideriv_ref(iterm,iab).gt.0) then
            irep=irepcol_ref(iterm,iab)

            dorb_psi_ref=dorb_psi_ref+aa(irep,jo,iab)
            dorb_energy_ref=dorb_energy_ref+tildem(irep,jo,iab)

            do i=1,nel
              dorb_psi=dorb_psi-aaz(irep,i,iab)*orb(i+ish,jo)
              dorb_energy=dorb_energy-emz(irep,i,iab)*orb(i+ish,jo)-aaz(irep,i,iab)*b(jo,i+ish)
            enddo
          endif

        enddo

        orbprim(iterm)=dorb_psi*detratio
        eorbprim(iterm)=dorb_energy*detratio+dorb_energy_ref-denergy*orbprim(iterm)

        orbprim(iterm)=orbprim(iterm)+dorb_psi_ref

 200  continue
          
      return
      end
c-----------------------------------------------------------------------
      subroutine optorb_compute(psid,eloc,deloc)

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'
      include 'optorb_cblk.h'




      common /zcompact/ zmat(MORB,MELEC,2,MSTATES),dzmat(MORB,MELEC,2,MSTATES)
     & ,emz(MELEC,MELEC,2,MSTATES),aaz(MELEC,MELEC,2,MSTATES)

      dimension psid(*),eloc(*),deloc(*)

      if(ioptorb.eq.0) return

      do 20 istate=1,nstates

        call optorb_deriv(psid(istate),deloc(istate)
     &   ,zmat(1,1,1,istate),dzmat(1,1,1,istate),emz(1,1,1,istate),aaz(1,1,1,istate)
     &   ,orb_o(1,istate),orb_ho(1,istate))
        
        do 20 i=1,norbterm
            orb_oe(i,istate)=orb_o(i,istate)*eloc(istate)
  20        orb_ho(i,istate)=orb_ho(i,istate)+eloc(istate)*orb_o(i,istate)

c     do iterm=1,norbterm
c        write(6,*) 'HELLO 1',iterm,orb_o(iterm,1),orb_ho(iterm,1),orb_oe(iterm,1)
c        write(6,*) 'HELLO 2',iterm,orb_o(iterm,2),orb_ho(iterm,2),orb_oe(iterm,2)
c     enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine optorb_sum(wtg_new,wtg_old,enew,eold,iflag)
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'
      include 'optorb_cblk.h'



      dimension wtg_new(*),wtg_old(*),enew(*),eold(*)

      if(ioptorb.eq.0) return

c     if(ns_current.ne.iorbsample) return
c ns_current reset
c     ns_current=0

      idiag_only=0
      call p2gtid('optwf:approx',iapprox,0,1)
      if(iapprox.gt.0) idiag_only=1

      do 200 istate=1,nstates

      p=wtg_new(istate)

      do 10 i=1,norbterm
       orb_o_sum(i,istate)=orb_o_sum(i,istate)+p*orb_o(i,istate)
       orb_oe_sum(i,istate) =orb_oe_sum(i,istate)+p*orb_oe(i,istate)
  10   orb_ho_cum(i,istate) =orb_ho_cum(i,istate)+p*orb_ho(i,istate)

      orb_wcum(istate)=orb_wcum(istate)+p
      orb_ecum(istate)=orb_ecum(istate)+p*enew(istate)

      if(isample_cmat.eq.0) go to 200

      if(idiag_only.eq.0) then
        idx=0
        do 20 i=1,nreduced
         ie=i
         do 20 j=1,i
          idx=idx+1
          je=j
  20      orb_oo_cum(idx,istate)=orb_oo_cum(idx,istate)+p*orb_o(ie,istate)*orb_o(je,istate)

        idx=0
        do 21 i=1,nreduced
         ie=i
         do 21 j=1,nreduced
          idx=idx+1
          je=j
  21      orb_oho_cum(idx,istate)=orb_oho_cum(idx,istate)+p*orb_o(je,istate)*orb_ho(ie,istate)
       else
        do 25 i=1,nreduced
          ie=i
          orb_oo_cum(i,istate)=orb_oo_cum(i,istate)+p*orb_o(ie,istate)*orb_o(ie,istate)
  25      orb_oho_cum(i,istate)=orb_oho_cum(i,istate)+p*orb_o(ie,istate)*orb_ho(ie,istate)
      endif

  200 continue

      if(iflag.eq.0) return

      do 300 istate=1,nstates

      q=wtg_old(istate)

      do 30 i=1,norbterm
       orb_o_sum(i,istate)=orb_o_sum(i,istate)+q*orb_o_old(i,istate)
  30   orb_oe_sum(i,istate) =orb_oe_sum(i,istate)+q*orb_oe_old(i,istate)

      orb_wcum(istate)=orb_wcum(istate)+q
      orb_ecum(istate)=orb_ecum(istate)+q*eold(istate)

      if(isample_cmat.eq.0) go to 300

      if(idiag_only.eq.0) then
        idx=0
        do 40 i=1,nreduced
         ie=i
         do 40 j=1,i
          idx=idx+1
          je=j
  40      orb_oo_cum(idx,istate)=orb_oo_cum(idx,istate)+q*orb_o_old(ie,istate)*orb_o_old(je,istate)

        idx=0
        do 41 i=1,nreduced
         ie=i
         do 41 j=1,nreduced
          idx=idx+1
          je=j
  41      orb_oho_cum(idx,istate)=orb_oho_cum(idx,istate)+q*orb_o_old(je,istate)*orb_ho_old(ie,istate)
       else
        do 45 i=1,nreduced
          ie=i
          orb_oo_cum(i,istate)=orb_oo_cum(i,istate)+q*orb_o_old(ie,istate)*orb_o_old(ie,istate)
  45      orb_oho_cum(i,istate)=orb_oho_cum(i,istate)+q*orb_o_old(ie,istate)*orb_ho_old(ie,istate)
      endif

  300 continue
      end
c-----------------------------------------------------------------------
      subroutine optorb_cum(wsum,esum)
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'
      include 'optorb_cblk.h'



      dimension wsum(*),esum(*)

      if(ioptorb.eq.0) return

      nb_current=nb_current+1

      do 200 istate=1,nstates

      orb_e_bsum(istate)=orb_e_bsum(istate)+esum(istate)
      orb_w_bsum(istate)=orb_w_bsum(istate)+wsum(istate)
      do 10 i=1,norbterm
       orb_o_bsum(i,istate)=orb_o_bsum(i,istate)+orb_o_sum(i,istate)
   10  orb_oe_bsum(i,istate)=orb_oe_bsum(i,istate)+orb_oe_sum(i,istate)

      if(nb_current.eq.nefp_blocks)then
       eb=orb_e_bsum(istate)/orb_w_bsum(istate)

       do 40 i=1,norbterm
         fnow=orb_oe_bsum(i,istate)/orb_w_bsum(istate)-orb_o_bsum(i,istate)/orb_w_bsum(istate)*eb
         orb_f_bcum(i,istate)=orb_f_bcum(i,istate)+fnow
   40    orb_f_bcm2(i,istate)=orb_f_bcm2(i,istate)+fnow**2

       orb_e_bsum(istate)=0.d0
       orb_w_bsum(istate)=0.d0
       do 50 i=1,norbterm
        orb_o_bsum(i,istate)=0.d0
   50   orb_oe_bsum(i,istate)=0.d0
      endif

      do 60 i=1,norbterm
       orb_o_cum(i,istate)=orb_o_cum(i,istate) + orb_o_sum(i,istate)
   60  orb_oe_cum(i,istate)=orb_oe_cum(i,istate) + orb_oe_sum(i,istate)

  200 continue

      if(nb_current.eq.nefp_blocks) then
        nb_current=0
        norb_f_bcum=norb_f_bcum+1
      endif

      if(idump_blockav.ne.0)then
       write(idump_blockav) esum(1)/wsum(1),(orb_o_sum(i,1)/wsum(1),orb_oe_sum(i,1)/wsum(1),i=1,norbterm)
      endif

      end

c-----------------------------------------------------------------------
      subroutine optorb_init(iflg)
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'
      include 'optorb_cblk.h'



      if(ioptorb.eq.0) return

      idiag_only=0
      call p2gtid('optwf:approx',iapprox,0,1)
      if(iapprox.gt.0) idiag_only=1

      do 100 istate=1,nstates

      do 10 i=1,norbterm
       orb_o_sum(i,istate)=0.d0
       orb_oe_sum(i,istate) =0.d0
       orb_o_bsum(i,istate)=0.d0
  10   orb_oe_bsum(i,istate)=0.d0
      orb_e_bsum(istate)=0.d0
      orb_w_bsum(istate)=0.d0

  100 continue
C$ iflg = 0: init *cum, *cm2 as well
      if(iflg.gt.0) return

      ns_current=0
      nb_current=0
      norb_f_bcum=0

      do 200 istate=1,nstates

      do 20 i=1,norbterm
       orb_o_cum(i,istate)=0.d0
       orb_oe_cum(i,istate) =0.d0
       orb_ho_cum(i,istate) =0.d0
       orb_f_bcum(i,istate)=0.d0
  20   orb_f_bcm2(i,istate)=0.d0
      orb_wcum(istate)=0.d0
      orb_ecum(istate)=0.d0

      if(isample_cmat.ne.0) then
       if(idiag_only.eq.0) then
         idx=0
         do 30 i=1,nreduced
          do 30 j=1,i
           idx=idx+1
  30       orb_oo_cum(idx,istate)=0.d0

         idx=0
         do 40 i=1,nreduced
          do 40 j=1,nreduced
           idx=idx+1
  40       orb_oho_cum(idx,istate)=0.d0
       else
         do 50 i=1,nreduced
           orb_oo_cum(i,istate)=0.d0
  50       orb_oho_cum(i,istate)=0.d0
       endif
      endif

  200 continue

      end
c-----------------------------------------------------------------------
      subroutine optorb_save
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'
      include 'optorb_cblk.h'



      if(ioptorb.eq.0) return

      do 200 istate=1,nstates

      do 10 i=1,norbterm
       orb_o_old(i,istate)=orb_o(i,istate)
       orb_oe_old(i,istate)=orb_oe(i,istate)
  10   orb_ho_old(i,istate)=orb_ho(i,istate)

  200 continue

      end
c-----------------------------------------------------------------------
      subroutine optorb_restore
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'
      include 'optorb_cblk.h'



      if(ioptorb.eq.0) return

      do 200 istate=1,nstates

      do 10 i=1,norbterm
       orb_o(i,istate)=orb_o_old(i,istate)
       orb_oe(i,istate)=orb_oe_old(i,istate)
  10   orb_ho(i,istate)=orb_ho_old(i,istate)

  200 continue

      end
c-----------------------------------------------------------------------
      subroutine optorb_avrg(wcum,eave,oav,eoav,fo,foerr,istate)
      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'mstates.h'
      include 'optorb.h'
      include 'optorb_cblk.h'
      dimension oav(*),eoav(*),fo(*),foerr(*)


      errn(x,x2,n)=dsqrt(dabs(x2/dble(n)-(x/dble(n))**2)/dble(n))

      if(ioptorb.eq.0) return

      do 30 i=1,norbterm
        oav(i)=orb_o_cum(i,istate)/wcum
        eoav(i)=orb_oe_cum(i,istate)/wcum
        fo(i)=eoav(i)-eave*oav(i)
   30   foerr(i)=errn(orb_f_bcum(i,istate),orb_f_bcm2(i,istate),norb_f_bcum)

      write(6,'(''ORB-PT: forces collected'',i4)') norb_f_bcum

      end
c-----------------------------------------------------------------------
      subroutine optorb_dump(iu)
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'
      include 'optorb_cblk.h'



      if(ioptorb.eq.0) return

      matdim=nreduced*(nreduced+1)/2
      call p2gtid('optwf:approx',iapprox,0,1)
      if(iapprox.gt.0) matdim=nreduced

      write(iu) norbprim,norbterm,nreduced
      write(iu) nefp_blocks,norb_f_bcum
      do 200 istate=1,nstates
      write(iu) (orb_o_cum(i,istate),i=1,norbterm)
      write(iu) (orb_oe_cum(i,istate),i=1,norbterm)
      write(iu) (orb_ho_cum(i,istate),i=1,norbterm)
      write(iu) (orb_f_bcum(i,istate),orb_f_bcm2(i,istate),i=1,norbterm)
      write(iu) (orb_oo_cum(i,istate),i=1,matdim)
      write(iu) (orb_oho_cum(i,istate),i=1,nreduced*nreduced)
      write(iu) orb_wcum(istate),orb_ecum(istate)
  200 continue

      end
c-----------------------------------------------------------------------
      subroutine optorb_rstrt(iu)
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'
      include 'optorb_cblk.h'



      if(ioptorb.eq.0) return
      read(iu) morbprim,morbterm,mreduced
      if(morbprim.ne.norbprim) then
       write (6,*) 'wrong number of primitive orb terms!'
       write (6,*) 'old ',morbprim,' new ',norbprim
       call fatal_error('OPTORB_RSTRT: Restart, inconsistent ORB information')
      endif
      if(morbterm.ne.norbterm) then
       write (6,*) 'wrong number of orb terms!'
       write (6,*) 'old ',morbterm,' new ',norbterm
       call fatal_error('OPTORB_RSTRT: Restart, inconsistent ORB information')
      endif

c nreduced has to be set since it will only be known for non-continuation runs
      nreduced=mreduced
      matdim=nreduced*(nreduced+1)/2 
      call p2gtid('optwf:approx',iapprox,0,1)
      if(iapprox.gt.0) matdim=nreduced

      read(iu) nefp_blocks,norb_f_bcum

      do 200 istate=1,nstates
      read(iu) (orb_o_cum(i,istate),i=1,norbterm)
      read(iu) (orb_oe_cum(i,istate),i=1,norbterm)
      read(iu) (orb_ho_cum(i,istate),i=1,norbterm)
      read(iu) (orb_f_bcum(i,istate),orb_f_bcm2(i,istate),i=1,norbterm)
      read(iu) (orb_oo_cum(i,istate),i=1,matdim)
      read(iu) (orb_oho_cum(i,istate),i=1,nreduced*nreduced)
      read(iu) orb_wcum,orb_ecum
  200 continue
      end
c-----------------------------------------------------------------------
      subroutine optorb_fin(wcum,ecum)
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      use optwf_parms, only: nparmd, nparme, nparmg, nparmj, nparml, nparms
      use sa_weights, only: iweight, nweight, weights
      implicit real*8(a-h,o-z)




      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'
      include 'optorb_cblk.h'
      include 'optci.h'
      include 'optjas.h'

      parameter(MPARMALL=MPARMJ+MXCIREDUCED+MXREDUCED)



      common /gradhess_all/ grad(MPARMALL),h(MPARMALL,MPARMALL),s(MPARMALL,MPARMALL)



      dimension oav(MXORBOP),eoav(MXORBOP),fo(MXORBOP),foerr(MXORBOP)
      dimension wcum(*),ecum(*)

      if(ioptorb.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      call p2gtid('optwf:iuse_orbeigv',iuse_orbeigv,0,1)
      call p2gtid('optwf:approx',iapprox,0,1)

      nparmd=max(nciterm-1,0)
      ish=nparmj+nparmd
      if(method.eq.'linear') ish=ish+1

      s(1,1)=0
      h(1,1)=0
      do 1 j=1,nreduced
        grad(j+ish)=0
        s(j+ish,1)=0
        h(j+ish,1)=0
        s(1,j+ish)=0
        h(1,j+ish)=0
        do 1 i=1,nreduced
          s(i+ish,j+ish)=0
   1      h(i+ish,j+ish)=0

      do 200 istate=1,nstates

      wts=weights(istate)

      passes=wcum(istate)
      passesi=1/passes
      eave=ecum(istate)*passesi

c     if(iorbsample.ne.1) then
c       passes=orb_wcum(istate)
c       passesi=1/passes
c       eave=orb_ecum(istate)*passesi
c     endif

      call optorb_avrg(passes,eave,oav(1),eoav(1),fo(1),foerr(1),istate)

c Hessian method
      if(method.eq.'hessian') then

        if(iuse_orbeigv.eq.0) then
c Formulas for exact orbital hessian not implemented
          call fatal_error('OPTORB_FIN: formulas for exact hessian not implemented')
        endif

c Linear method
       elseif(method.eq.'linear') then

        s(1,1)=1
        h(1,1)=h(1,1)+wts*eave
c Exact Hamiltonian 
        if(iuse_orbeigv.eq.0) then

c Hamiltonian on semi-orthogonal basis
        idx=0
        do 30 i=1,nreduced
          s(i+ish,1)=0
          s(1,i+ish)=0
          h(i+ish,1)=h(i+ish,1)+wts*(eoav(i)-eave*oav(i))
          h(1,i+ish)=h(1,i+ish)+wts*(orb_ho_cum(i,istate)*passesi-eave*oav(i))
c         write(6,*) 'H',wts,eoav(i)-eave*oav(i),orb_ho_cum(i,istate)*passesi-eave*oav(i)
          i0=1
          if(iapprox.gt.0) i0=i
          do 30 j=i0,i
            idx=idx+1
            orb_oo=orb_oo_cum(idx,istate)*passesi-oav(i)*oav(j)
            s(i+ish,j+ish)=s(i+ish,j+ish)+wts*orb_oo
   30       s(j+ish,i+ish)=s(i+ish,j+ish)

        i0=1
        i1=nreduced
        idx=0
        do 40 i=1,nreduced
          if(iapprox.gt.0) then
            i0=i
            i1=i
          endif
          do 40 j=i0,i1
            idx=idx+1
            orb_oho=(orb_oho_cum(idx,istate)-oav(j)*orb_ho_cum(i,istate))*passesi
     &             -oav(i)*eoav(j)+eave*oav(i)*oav(j)
   40       h(j+ish,i+ish)=h(j+ish,i+ish)+wts*orb_oho

       endif

c Perturbative method
       elseif(method.eq.'perturbative') then
            
        if(iuse_orbeigv.eq.0) then
c Formulas for exact orbital perturbative not implemented
          call fatal_error('OPTORB_FIN: formulas for exact perturbative not implemented')
         else
          do 60 i=1,nreduced
   60       grad(i)=grad(i)+wts*fo(i)
          idx=0
          do 70 i=1,nreduced
            do 70 j=1,i
              idx=idx+1
              s(i,j)=s(i,j)+wts*(orb_oo_cum(idx,istate)*passesi-oav(i)*oav(j))
   70         s(j,i)=s(i,j)
        endif
      endif

  200 continue

c Approximations on matrix elements
      if(method.eq.'linear') then
        if(iapprox.gt.0) then
          do 230 i=1,nreduced
            do 230 j=1,i-1
              s(i+ish,j+ish)=0
              s(j+ish,i+ish)=0
              h(i+ish,j+ish)=0
  230         h(j+ish,i+ish)=0
          if(iapprox.eq.2) then
            do 240 i=1,nreduced
  240         h(1,i+ish)=h(i+ish,1)
          endif
         elseif(iapprox.lt.0) then
          if(iapprox.eq.-1) then
            do 250 i=1,nreduced
  250         h(1,i+ish)=h(i+ish,1)
           elseif(iapprox.eq.-2) then
            do 260 i=1,nreduced
              h(1,i+ish)=h(i+ish,1)
              do 260 j=1,i-1
                h(i+ish,j+ish)=0.5*(h(i+ish,j+ish)+h(j+ish,i+ish))
  260           h(j+ish,i+ish)=h(i+ish,j+ish)
           elseif(iapprox.eq.-3) then
            do 270 i=1,nreduced
              h(1,i+ish)=0.5*(h(i+ish,1)+h(1,i+ish))
              h(i+ish,1)=h(1,i+ish)
              do 270 j=1,i-1
                h(i+ish,j+ish)=0.5*(h(i+ish,j+ish)+h(j+ish,i+ish))
  270           h(j+ish,i+ish)=h(i+ish,j+ish)
          endif
        endif
       elseif(method.eq.'perturbative') then
c Approximation: diagonal perturbative approach
        if(iapprox.gt.0) then
          do 280 i=1,nreduced
            do 280 j=1,i-1
              s(j,i)=0
  280         s(i,j)=0
        endif
      endif

      if(idump_blockav.ne.0) close(idump_blockav)

      end
c-----------------------------------------------------------------------
      subroutine detratio_col(nel,orb,icol,sinvt,ratio,isltnew)
      implicit real*8(a-h,o-z)
c values of new orbital
      dimension orb(nel)
c inverse transposed slater matrix (first index electron, 2nd orbital)
      dimension sinvt(nel,nel)
c compute ratio of new and old determinant, if isltnew is
c not zero, update inverse slater matrix as well
c the new determinant differs from the old by replacing column icol
c with the orbital values in orb

      ratio=0.d0
      do ie=1,nel
       ratio=ratio+sinvt(icol,ie)*orb(ie)
      enddo
      if(isltnew.gt.0) then
c matrix except replaced column
       do jcol=1,nel
        if(jcol.ne.icol) then
         sum=0.d0
         do je=1,nel
          sum=sum+orb(je)*sinvt(jcol,je)
         enddo
         sum=sum/ratio
         do je=1,nel
          sinvt(jcol,je)=sinvt(jcol,je)-sum*sinvt(icol,je)
         enddo
        endif
       enddo
c replaced column
       do ie=1,nel
        sinvt(icol,ie)=sinvt(icol,ie)/ratio
       enddo
      endif

      end
c-----------------------------------------------------------------------
      subroutine optorb_define
      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use dets, only: cdet, ndet
      use elec, only: ndn, nup
      use multidet, only: iactv, irepcol_det, ireporb_det, ivirt, iwundet, kref, numrep_det

      use optorb_mix, only: iwmix_virt, norbopt, norbvirt
      implicit real*8(a-h,o-z)






      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'
      include 'optorb_cblk.h'
      include 'inputflags.h'

      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb

      common /dorb/ iworbd(MELEC,MDET)
      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb

      common /optorb/ orb_energy(MORB),dmat_diag(MORB),irrep(MORB)

      data icount_orbdef /1/

      dimension iodet(2,MDET),iopos(2,MDET),iflag(2,MORB)

      dimension ne(2),m(2)

      save icount_orbdef

      iprt=3

      call p2gti('electrons:nelec',nelec,1)
      call p2gti('electrons:nup',nup,1)
      ndn=nelec-nup

      ne(1)=nup
      ne(2)=nelec
c orbital indices in determinants of trial wave function
      ndetorb=0

      do i=1,ndet
       do j=1,nelec
        if(iworbd(j,i).gt.norb) then
         write(6,1) i,j,iworbd(j,i),norb
         call fatal_error('VERIFY: orbital index out of range')
        endif
        if(iworbd(j,i).gt.ndetorb) then
         ndetorb=iworbd(j,i)
        endif
       enddo
      enddo
  1   format('Det ',i4,' column ',i4,' orb index ',i4,' norb ',i4)

c Number of external orbitals for orbital optimization
      call p2gtid('optwf:ncore',ncore,0,1)
      next_max=norb-ndetorb
      call p2gtid('optwf:nextorb',nadorb,next_max,1)
      if(nadorb.gt.next_max) nadorb=next_max
      if(iprt.gt.0) then
       write(6,'(''Determinantal orbitals in orbital optimization: '',i4)') ndetorb
       write(6,'(''External orbitals in orbital optimization: '',i4)') nadorb
       write(6,'(''Total orbitals in orbital optimization: '',i4)') nadorb+ndetorb-ncore
      endif
      norb=ndetorb

c Omit doubly occupied in all input determinants
      do 5 i=1,ndetorb
        iflag(1,i)=0
        do 3 k=1,ndet
          iocc=0
          do 2 j=1,nelec
   2        if(iworbd(j,k).eq.i) iocc=iocc+1
          if(iocc.ne.2) then
            iflag(1,i)=1
            goto 5
          endif
   3    continue
   5  continue
c Omit empty orbitals
      do 6 i=1,ndetorb
       iflag(2,i)=0
       do 6 k=1,ndet
        do 6 j=1,nelec
   6      if(iworbd(j,k).eq.i) iflag(2,i)=1
      do 8 i=ndetorb+1,ndetorb+nadorb
       iflag(1,i)=1
   8   iflag(2,i)=0

      if(norbopt.eq.0.or.norbvirt.eq.0) then
        do 9 io=1,ndetorb
         do 9 jo=ncore+1,ndetorb+nadorb
   9      iwmix_virt(io,jo)=jo
      elseif(norbopt.ne.ndetorb.or.norbvirt.lt.nadorb) then
       write(6,'(''OPTORB_DEFINE: norbopt,ndetorb'',2i6)') norbopt,ndetorb
       write(6,'(''OPTORB_DEFINE: noptvirt,nadorb'',2i6)') norbvirt,nadorb
       call fatal_error('OPTORB_DEFINE: Mixvirt block, inconsistent')
      endif

c Orbital variation io -> io+a*jo
c io: occupied orbitals in twf
c jo: all orbitals
c omitted if not same symmetry, or io empty, or both doubly occupied
      noporb=0
      iterm=0

      if(iprt.gt.2) then
       write(6,'(''=========== orbital pair list =========='')')
      endif

      call p2gtid('optwf:no_active',no_active,0,1)
      do 60 io=ncore+1,ndetorb
c Omit empty orbitals
       if(iflag(2,io).eq.0) goto 60
       do 50 jo=ncore+1,ndetorb+nadorb
c Omit if io and jo are the same
        if(io.eq.jo) goto 50
c Omit if io and jo have different symmetry
        if(irrep(io).ne.irrep(jo)) goto 50
c Omit if io and jo are both doubly occupied in all determinants
        if((iflag(1,io).eq.0).and.(iflag(1,jo).eq.0)) goto 50
c Omit if io and jo are both active orbitals
        if(no_active.ne.0.and.iflag(1,io).ne.0.and.iflag(2,jo).ne.0) goto 50
c Omit if we only want to mix according to the table mixvirt
        if(iwmix_virt(io,jo).eq.0) goto 50
c Include: io is occupied in some determinant and jo not
        do 40 iab=1,2
          n0=0
          n1=nup
          if(iab.eq.2) then
            n0=nup
            n1=ndn
          endif
          m(iab)=0
          do 30 k=1,ndet
            do 15 ie=1,n1
              if(iworbd(ie+n0,k).eq.io) then
                iesave=ie
                goto 20
              endif
 15         continue
            goto 30
 20         continue
            do 25 ie=1,n1
 25           if(iworbd(ie+n0,k).eq.jo) goto 30
            m(iab)=m(iab)+1
            iodet(iab,m(iab))=k
            iopos(iab,m(iab))=iesave
 30       continue
 40     continue
        if(m(1)+m(2).eq.0) then
          if(iprt.gt.3) write(6,'(''no appropriate determinant for '',2i4)') io,jo
          goto 50
        endif

c Define new operator (new variation) and its terms
        noporb=noporb+1
        if(noporb.gt.MXORBOP) then
          write(6,'(''noporb,max_orb'',2i5)') noporb,MXORBOP
          call fatal_error('ORB_DEFINE: too many terms, increase MXORBOP')
        endif

        ideriv(1,noporb)=io
        ideriv(2,noporb)=jo
        ideriv_iab(noporb)=0
        if(m(1).gt.0) ideriv_iab(noporb)=1
        if(m(2).gt.0) ideriv_iab(noporb)=ideriv_iab(noporb)+2

        do iab=1,2
          n0=0
          n1=nup
          if(iab.eq.2) then
            n0=nup
            n1=ndn
          endif
          ideriv_ref(noporb,iab)=0
          do i=1,n1
            if(iworbd(i+n0,kref).eq.io) then
              ideriv_ref(noporb,iab)=1
              irepcol_ref(noporb,iab)=i
            endif
          enddo
        enddo
        if(iprt.gt.2) write(6,'(a16,i4,a8,i4,i5,a15,i4)') 'new variation: ',noporb,' pair ',io,jo,' spin ',ideriv_iab(noporb)

 50    continue
 60   continue

      norbterm=noporb
      write(6,'(''number of orbital variations: '',2i8)') norbterm

c if mix_n, optorb_define called mutiple times with method=sr_n or lin_d
      if(icount_orbdef.eq.1) call p2gtad('optwf:method',method,'linear',1)
      if(method.eq.'linear') then

        if(MXREDUCED.ne.MXORBOP) call fatal_error('READ_INPUT: MXREDUCED.ne.MXORBOP')
        nreduced=norbterm
       elseif(method.eq.'sr_n'.or.method.eq.'lin_d'.or.method.eq.'mix_n') then
        nreduced=1
      endif

      icount_orbdef=icount_orbdef+1

      return
      end
c-----------------------------------------------------------------------
      subroutine check_orbitals

c Do not compute virtual orbitals during single-electron move
      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optjas.h'
      include 'optci.h'
      include 'optorb.h'

      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb

      save nadorb_save

      nadorb_save=nadorb
      nadorb=0

      return

      entry check_orbitals_reset

      nadorb=nadorb_save

      return
      end
