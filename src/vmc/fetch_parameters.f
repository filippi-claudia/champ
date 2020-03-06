c-----------------------------------------------------------------------
      subroutine fetch_parameters(p)
      implicit real*8(a-h,o-z)

      dimension p(*)

      n=0

      ip=1
      call fetch_jastrow(p(1),n)
      nparm=n

      call fetch_lcao(p(nparm+1),n)
      nparm=nparm+n

      call fetch_ci(p(nparm+1),n)
      nparm=nparm+n

      return
      end
c-----------------------------------------------------------------------
      subroutine fetch_jastrow(p,n)
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use jaspar, only: nspin1, nspin2, sspin, sspinn, is
      use jaspar3, only: a, b, c, fck, nord, scalek

      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'

      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc
      common /bparm/ nspin2b,nocuspb

      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj
      common /optwf_wjas/ iwjasa(83,MCTYP3X),iwjasb(83,3),iwjasc(83,MCTYPE),iwjasf(15,MCTYPE)
      common /optwf_nparmj/ nparma(MCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE)

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      dimension p(*)

      if(ioptjas.eq.0) return

      iparm=0
      do 50 ict=1,nctype
        do 50 i=1,nparma(ict)
          iparm=iparm+1
   50     p(iparm)=a4(iwjasa(i,ict),ict,1)
      do 60 i=1,nparmb(1)
        iparm=iparm+1
   60   p(iparm)=b(iwjasb(i,1),1,1)
      do 70 ict=1,nctype
        do 70 i=1,nparmc(ict)
          iparm=iparm+1
   70     p(iparm)=c(iwjasc(i,ict),ict,1)

      n=iparm

      return
      end
c-----------------------------------------------------------------------
      subroutine fetch_lcao(p,n)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'
      include 'optorb_cblk.h'

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj

      dimension p(*)

      if(ioptorb.eq.0) return

      do 30 i=1,norbterm
 30    p(i)=0.d0
      n=norbterm

      return
      end
c-----------------------------------------------------------------------
      subroutine fetch_ci(p,n)
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use dets, only: cdet, ndet
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'


      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj

      dimension p(*)

      if(ioptci.eq.0) return

      if(ncsf.eq.0) then
        do 35 idet=2,ndet
          p(idet-1)=cdet(idet,1,1)
 35     continue
        n=ndet-1
       else
        do 50 icsf=2,ncsf
          p(icsf-1)=ccsf(icsf,1,1)
 50     continue
        n=ncsf-1
      endif

c     do 90 j=1,nstates
c90     write(6,'(''csf ='',1000f20.15)') (ccsf(i,j,iadiag),i=1,ncsf)

      return
      end
c-----------------------------------------------------------------------
