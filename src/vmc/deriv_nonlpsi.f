      function deriv_psinl(u,rshifti,rshiftj,rri,rrj,gn,it)
c Written by Claudia Filippi, modified by Cyrus Umrigar

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use jaspar3, only: a, c
      use jaspar4, only: nordc
      use jaspar6, only: asymp_r
      use optwf_wjas, only: iwjasc
      use wfsec, only: iwf
      use contr2, only: ijas

      use vardep, only: cdep, iwdepend, nvdepend

      use cuspmat4, only: d, iwc4, nterms
      implicit real*8(a-h,o-z)





      parameter (zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,eps=1.d-12)

      dimension rshifti(3),rshiftj(3),gn(*)
      dimension uu(0:MORDJ),ss(0:MORDJ),tt(0:MORDJ)

      if(ijas.ge.4.and.ijas.le.6) then

        deriv_psinl=0
        if(nordc.eq.0) return

        if(rri.eq.asymp_r .or. rrj.eq.asymp_r) return
        do 37 k=1,3
   37     if(abs(rshifti(k)-rshiftj(k)).gt.eps) return

        uu(1)=u
        rrri=rri
        rrrj=rrj
        call switch_scale(uu(1))
        call switch_scale(rrri)
        call switch_scale(rrrj)

        uu(0)=1
        ss(0)=2
        tt(0)=1
        do 40 jp=1,nordc
          uu(jp)=uu(1)*uu(jp-1)
          ss(jp)=rrri**jp+rrrj**jp
   40     tt(jp)=(rrri*rrrj)**jp

        ll=0
        jj=1
        jparm=1
        do 50 n=2,nordc
          do 50 k=n-1,0,-1
            if(k.eq.0) then
              l_hi=n-k-2
             else
              l_hi=n-k
            endif
            do 50 l=l_hi,0,-1
              m=(n-k-l)/2
              if(2*m.eq.n-k-l) then
                ll=ll+1
                p=uu(k)*ss(l)*tt(m)
                deriv_psinl=deriv_psinl+c(ll,it,iwf)*p

                ideriv=0
                if(ll.eq.iwjasc(jparm,it)) then
                  ideriv=2
                 else
                  do 31 id=1,2*(nordc-1)
                    if(ll.eq.iwc4(id)) then
                      jj=id
                      if(nvdepend(jj,it).gt.0) ideriv=1
                    endif
   31             continue
                endif

                if(ideriv.eq.1) then
                  do 43 id=1,nvdepend(jj,it)
                    iparm=iwdepend(jj,id,it)
   43               gn(iparm)=gn(iparm)+cdep(jj,id,it)*p
c                 jj=jj+1
                 elseif(ideriv.eq.2) then
                  gn(jparm)=gn(jparm)+p
                  jparm=jparm+1
                endif
              endif
   50   continue

      endif

      return
      end

c-----------------------------------------------------------------------
      function deriv_psianl(rri,gn,it)

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use jaspar3, only: a, c

      use jaspar4, only: a4, norda
      use jaspar6, only: asymp_jasa, asymp_r
      use optwf_nparmj, only: nparma
      use optwf_wjas, only: iwjasa
      use wfsec, only: iwf
      use contr2, only: ijas
      implicit real*8(a-h,o-z)


      parameter(one=1.d0)




      dimension gn(*)

c Note: This routine is only called with iwf=1, but parts of it are
c written for general iwf, whereas others (asymp_r) assume iwf=1.

      if(rri.eq.asymp_r) then
        deriv_psianl=0
        return
      endif

      if(ijas.ge.4.and.ijas.le.6) then

        deriv_psianl=a4(1,it,iwf)*rri/(one+a4(2,it,iwf)*rri)-asymp_jasa(it)
        do 20 i=2,norda
   20     deriv_psianl=deriv_psianl+a4(i+1,it,iwf)*rri**i
        do 30 jparm=1,nparma(it)
            if(iwjasa(jparm,it).eq.1) then
              top=rri
              bot=one+a4(2,it,iwf)*rri
              gen=top/bot-asymp_r/(1+a4(2,it,iwf)*asymp_r)
             elseif(iwjasa(jparm,it).eq.2) then
              top=-a4(1,it,iwf)*rri*rri
              bot=one+a4(2,it,iwf)*rri
              bot=bot*bot
              gen=top/bot+a4(1,it,iwf)*asymp_r**2/(1+a4(2,it,iwf)*asymp_r)**2
             else
              iord=iwjasa(jparm,it)-1
              gen=rri**iord-asymp_r**iord
            endif
            gn(jparm)=gn(jparm)+gen
   30    continue
      endif

      return
      end

c-----------------------------------------------------------------------
      function deriv_psibnl(u,gn,isb,ipar)

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use jaspar, only: sspinn, is
      use jaspar3, only: a, b, c

      use jaspar4, only: nordb
      use jaspar6, only: asymp_jasb, asymp_r
      use optwf_nparmj, only: nparmb
      use optwf_wjas, only: iwjasb
      use wfsec, only: iwf
      use contr2, only: ijas
      implicit real*8(a-h,o-z)










      parameter(one=1.d0)




      dimension gn(*)

c Note: This routine is only called with iwf=1, but parts of it are
c written for general iwf, whereas others (asymp_r) assume iwf=1.

c Not updated for ijas=5,6 because we will probably stay with ijas=4
c If we want to use ijas=5,6 update this routine similarly to psi.f
      if(ijas.ge.5) stop 'ijas >= 5 not implemented in psibnl'

      if(u.eq.asymp_r) then
        deriv_psibnl=0
        return
      endif

      fee=b(1,isb,iwf)*u/(one+b(2,isb,iwf)*u)

      deriv_psibnl=sspinn*fee-asymp_jasb(ipar+1)
      if(ijas.ge.4.and.ijas.le.6) then
        do 10 i=2,nordb
   10     deriv_psibnl=deriv_psibnl+b(i+1,isb,iwf)*u**i
      endif

      do 20 jparm=1,nparmb(isb)
        if(iwjasb(jparm,isb).eq.1) then
          top=u
          bot=one+b(2,isb,iwf)*u
          gee=sspinn*(top/bot-asymp_r/(1+b(2,isb,iwf)*asymp_r))
         elseif(iwjasb(jparm,isb).eq.2) then
          top=-b(1,isb,iwf)*u*u
          bot=one+b(2,isb,iwf)*u
          bot=bot*bot
          gee=sspinn*(top/bot+b(1,isb,iwf)*asymp_r**2/(1+b(2,isb,iwf)*asymp_r)**2)
         else
          iord=iwjasb(jparm,isb)-1
          gee=u**iord-asymp_r**iord
        endif
        gn(jparm)=gn(jparm)+gee
  20  continue

      return
      end
