      subroutine optx_jas_ci_sum(p,q,enew,eold)

      use derivjas, only: gvalue
      use gradhessjo, only: denergy_old, gvalue_old
      use mix_jas_ci, only: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci
      use optwf_contrl, only: ioptci, ioptjas
      use optwf_parms, only: nparmj
      use deloc_dj_m, only: denergy

      use ci000, only: iciprt, nciprim, nciterm

      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'mstates.h'
      include 'optjas.h'
      include 'optci.h'
      include 'optci_cblk.h'


      if(ioptjas.eq.0.or.ioptci.eq.0) return

      do 10 i=1,nparmj
        do 10 j=1,nciterm
        dj_o_ci(i,j)=dj_o_ci(i,j)  +p*gvalue(i)*ci_o(j)+q*gvalue_old(i)*ci_o_old(j)
        dj_oe_ci(i,j)=dj_oe_ci(i,j)+p*gvalue(i)*ci_o(j)*enew+q*gvalue_old(i)*ci_o_old(j)*eold
        de_o_ci(i,j)=de_o_ci(i,j)  +p*denergy(i,1)*ci_o(j)+q*denergy_old(i,1)*ci_o_old(j)
  10    dj_de_ci(i,j)=dj_de_ci(i,j)+p*gvalue(i)*ci_de(j)+q*gvalue_old(i)*ci_de_old(j)

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_jas_ci_init

      use mix_jas_ci, only: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci
      use optwf_contrl, only: ioptci, ioptjas
      use optwf_parms, only: nparmj
      use ci000, only: iciprt, nciprim, nciterm

      implicit real*8(a-h,o-z)




      include 'vmc.h'
      include 'mstates.h'
      include 'optjas.h'
      include 'optci.h'
      include 'optci_cblk.h'



      if(ioptjas.eq.0.or.ioptci.eq.0) return

      do 10 i=1,nparmj
        do 10 j=1,nciterm
          dj_o_ci(i,j)=0
          dj_oe_ci(i,j)=0
          de_o_ci(i,j)=0
  10      dj_de_ci(i,j)=0

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_jas_ci_dump(iu)

      use mix_jas_ci, only: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci

      use optwf_contrl, only: ioptci, ioptjas
      use optwf_parms, only: nparmj
      use ci000, only: iciprt, nciprim, nciterm

      implicit real*8(a-h,o-z)




      include 'vmc.h'
      include 'mstates.h'
      include 'optjas.h'
      include 'optci.h'
      include 'optci_cblk.h'



      if(ioptjas.eq.0.or.ioptci.eq.0) return
      write(iu) ((dj_o_ci(i,j),dj_oe_ci(i,j),dj_de_ci(i,j),de_o_ci(i,j),i=1,nparmj),j=1,nciterm)

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_jas_ci_rstrt(iu)

      use mix_jas_ci, only: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci

      use optwf_contrl, only: ioptci, ioptjas
      use optwf_parms, only: nparmj
      use ci000, only: iciprt, nciprim, nciterm

      implicit real*8(a-h,o-z)




      include 'vmc.h'
      include 'mstates.h'
      include 'optjas.h'
      include 'optci.h'
      include 'optci_cblk.h'



      if(ioptjas.eq.0.or.ioptci.eq.0) return
      read(iu) ((dj_o_ci(i,j),dj_oe_ci(i,j),dj_de_ci(i,j),de_o_ci(i,j),i=1,nparmj),j=1,nciterm)

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_jas_ci_fin(passes,eave)
      use csfs, only: ccsf, ncsf

      use dets, only: cdet
      use gradhess_ci, only: grad_ci
      use gradhess_jas, only: grad_jas
      use gradhess_mix_jas_ci, only: h_mix_jas_ci, s_mix_jas_ci

      use mix_jas_ci, only: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci

      use optwf_contrl, only: ioptci, ioptjas
      use optwf_parms, only: nparmj
      use gradhessj, only: de, dj, dj_e

      use ci000, only: iciprt, nciprim, nciterm

      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'mstates.h'
      include 'optjas.h'
      include 'force.h'
      include 'optci.h'
      include 'optci_cblk.h'


      dimension oelocav(MXCITERM),eav(MXCITERM)

      if(ioptjas.eq.0.or.ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      if(method.eq.'hessian') then

c Compute mix Hessian
      do 10 i=1,nparmj
        do 10 j=1,nciterm
          h1=2*(2*(dj_oe_ci(i,j)-eave*dj_o_ci(i,j))-dj(i,1)*grad_ci(j)-grad_jas(i)*ci_o_cum(j))
          h2=de_o_ci(i,j)-de(i,1)*ci_o_cum(j)/passes
     &         +dj_de_ci(i,j)-dj(i,1)*ci_de_cum(j)/passes
  10      h_mix_jas_ci(i,j)=(h1+h2)/passes

      write(21,*) nciterm
      write(21,*) ((h_mix_jas_ci(i,j),j=1,nciterm),i=1,nparmj)

      elseif(method.eq.'linear') then

      if(ncsf.eq.0) then
        do 20 i=1,nciterm
          oelocav(i)=0
          eav(i)=0
          do 20 j=1,nciterm
            oelocav(i)=oelocav(i)+ci_oe_cum(i,j)*cdet(j,1,1)/passes
  20        eav(i)=eav(i)+ci_oe_cum(j,i)*cdet(j,1,1)/passes
       else
        do 25 i=1,ncsf
          oelocav(i)=0
          eav(i)=0
          do 25 j=1,ncsf
            oelocav(i)=oelocav(i)+ci_oe_cum(i,j)*ccsf(j,1,1)/passes
  25        eav(i)=eav(i)+ci_oe_cum(j,i)*ccsf(j,1,1)/passes
      endif

      do 30 i=1,nparmj
        do 30 j=1,nciterm
c Overlap s_jas_ci
          s_mix_jas_ci(i,j)=(dj_o_ci(i,j)-dj(i,1)*ci_o_cum(j)/passes)/passes
c H matrix h_jas_ci
          h_mix_jas_ci(i,j)=(dj_de_ci(i,j)+dj_oe_ci(i,j)
     &    +eave*dj(i,1)*ci_o_cum(j)/passes-dj(i,1)*eav(j)-ci_o_cum(j)*dj_e(i,1)/passes)/passes
c H matrix h_ci_jas
   30     h_mix_jas_ci(i+nparmj,j)=(de_o_ci(i,j)+dj_oe_ci(i,j)
     &    +eave*dj(i,1)*ci_o_cum(j)/passes-dj(i,1)*oelocav(j)-ci_o_cum(j)*(de(i,1)+dj_e(i,1))/passes)/passes
         
      endif

      return
      end
