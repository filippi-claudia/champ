      subroutine optx_orb_ci_sum(p,q)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'mstates.h'
      include 'optci.h'
      include 'optci_cblk.h'
      include 'optorb.h'
      include 'optorb_cblk.h'

      common /mix_orb_ci/ ci_o_o(MXCITERM,MXREDUCED),ci_o_oe(MXCITERM,MXREDUCED),
     &ci_de_o(MXCITERM,MXREDUCED),ci_o_ho(MXCITERM,MXREDUCED)

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      if(ioptorb.eq.0.or.ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do 10 j=1,nreduced
       do 10 i=1,nciterm
        ci_o_o(i,j)=ci_o_o(i,j)  +p*ci_o(i)*orb_o(j,1)+q*ci_o_old(i)*orb_o_old(j,1)
        ci_de_o(i,j)=ci_de_o(i,j)+p*ci_de(i)*orb_o(j,1)+q*ci_de_old(i)*orb_o_old(j,1)
        ci_o_ho(i,j)=ci_o_ho(i,j)+p*ci_o(i)*orb_ho(j,1)+q*ci_o_old(i)*orb_ho_old(j,1)
  10    ci_o_oe(i,j)=ci_o_oe(i,j)+p*ci_o(i)*orb_oe(j,1)+q*ci_o_old(i)*orb_oe_old(j,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_orb_ci_init

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'optci.h'
      include 'optorb.h'

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      common /mix_orb_ci/ ci_o_o(MXCITERM,MXREDUCED),ci_o_oe(MXCITERM,MXREDUCED),
     &ci_de_o(MXCITERM,MXREDUCED),ci_o_ho(MXCITERM,MXREDUCED)

      if(ioptorb.eq.0.or.ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do 10 j=1,nreduced
        do 10 i=1,nciterm
          ci_o_o(i,j)=0
          ci_o_oe(i,j)=0
          ci_o_ho(i,j)=0
  10      ci_de_o(i,j)=0

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_orb_ci_dump(iu)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'optci.h'
      include 'optorb.h'

      common /mix_orb_ci/ ci_o_o(MXCITERM,MXREDUCED),ci_o_oe(MXCITERM,MXREDUCED),
     &ci_de_o(MXCITERM,MXREDUCED),ci_o_ho(MXCITERM,MXREDUCED)

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      if(ioptorb.eq.0.or.ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return
      write(iu) ((ci_o_o(i,j),ci_o_oe(i,j),ci_o_ho(i,j),ci_de_o(i,j),i=1,nciterm),j=1,nreduced)

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_orb_ci_rstrt(iu)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'optci.h'
      include 'optorb.h'

      common /mix_orb_ci/ ci_o_o(MXCITERM,MXREDUCED),ci_o_oe(MXCITERM,MXREDUCED),
     &ci_de_o(MXCITERM,MXREDUCED),ci_o_ho(MXCITERM,MXREDUCED)

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      if(ioptorb.eq.0.or.ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return
      write(iu) ((ci_o_o(i,j),ci_o_oe(i,j),ci_o_ho(i,j),ci_de_o(i,j),i=1,nciterm),j=1,nreduced)

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_orb_ci_fin(passes,eave)
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use dets, only: cdet, ndet
      use gradhess_ci, only: grad_ci, h_ci, s_ci
      use gradhess_mix_orb_ci, only: h_mix_ci_orb, s_mix_ci_orb
      implicit real*8(a-h,o-z)




      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optci.h'
      include 'optci_cblk.h'
      include 'optorb.h'
      include 'optorb_cblk.h'
      include 'optjas.h'

      parameter(MPARMALL=MPARMJ+MXCIREDUCED+MXREDUCED)

      common /mix_orb_ci/ ci_o_o(MXCITERM,MXREDUCED),ci_o_oe(MXCITERM,MXREDUCED),
     &ci_de_o(MXCITERM,MXREDUCED),ci_o_ho(MXCITERM,MXREDUCED)

c     common /gradhess_orb/ grad_orb(MXORBOP),h_orb(MXMATDIM),s_orb(MXMATDIM)

      common /gradhess_all/ grad(MPARMALL),h(MPARMALL,MPARMALL),s(MPARMALL,MPARMALL)



      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj

      dimension oelocav(MXCITERM),eav(MXCITERM)

      if(ioptorb.eq.0.or.ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      if(method.eq.'hessian') then

      ishift=nparmj+nciterm-1
      do 30 i=1,nciterm
        do 30 j=1,norbprim
          h1=2*(ci_o_oe(i,j)-eave*ci_o_o(i,j)-ci_o_cum(i)*grad(j+ishift)-grad_ci(i)*orb_o_cum(j,1))
          h2=2*(ci_de_o(i,j)-ci_de_cum(i)*orb_o_cum(j,1)/passes)
   30     h_mix_ci_orb(i,j)=(h1+h2)/passes

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

      do 70 i=1,nciterm
        do 70 j=1,nreduced
c Overlap s_ij
          s_mix_ci_orb(i,j)=(ci_o_o(i,j)-ci_o_cum(i)*orb_o_cum(j,1)/passes)/passes
c Hamiltonian ci_orb
          h_mix_ci_orb(i,j)=(ci_o_ho(i,j)+eave*ci_o_cum(i)*orb_o_cum(j,1)/passes
     &    -ci_o_cum(i)*orb_ho_cum(j,1)/passes-orb_o_cum(j,1)*oelocav(i))/passes
c Hamiltonian orb_ci
   70     h_mix_ci_orb(i+nciterm,j)=(ci_de_o(i,j)+ci_o_oe(i,j)+eave*ci_o_cum(i)*orb_o_cum(j,1)/passes
     &    -ci_o_cum(i)*orb_oe_cum(j,1)/passes-orb_o_cum(j,1)*eav(i))/passes

      endif

      return
      end
