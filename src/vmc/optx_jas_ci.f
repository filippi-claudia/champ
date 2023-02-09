      module optx_jas_ci
      contains
      subroutine optx_jas_ci_sum(p,q,enew,eold)

      use derivjas, only: gvalue
      use gradhessjo, only: denergy_old, gvalue_old
      use mix_jas_ci, only: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci
      use optwf_control, only: ioptci, ioptjas
      use optwf_parms, only: nparmj
      use deloc_dj_m, only: denergy
      use ci000, only: nciterm
      use ci001_blk, only: ci_o
      use ci002_blk, only: ci_o_old
      use ci004_blk, only: ci_de, ci_de_old
      use precision_kinds, only: dp
      use csfs, only: nstates
      use vmc_mod, only: stoj
      use contrl_file, only: ounit

      implicit none

      integer :: i, j, k
      real(dp) :: enew, eold, p, q

      if(ioptjas.eq.0.or.ioptci.eq.0) return

      !STU in principle shouldn't the denergy (which have a nparmj index), (found in de_o_ci) use a stoj(k), and not just a '1' as shown. Then these will all be nbjx (mixed) sized quantities.
      do k=1,nstates !STU add jastrow/orb mapping here 
        do j=1,nciterm !STU ci variables need an index added eventually
          do i=1,nparmj !STU i think you should use your stobjx mapping IF ci_o is indeed a istate object (because psidi(istate)) then all of  these are istate objects!!!
            dj_o_ci(i,j,k)=dj_o_ci(i,j,k)+p*gvalue(i,stoj(k))*ci_o(j,k)
     &              +q*gvalue_old(i,stoj(k))*ci_o_old(j,k)
            dj_oe_ci(i,j,k)=dj_oe_ci(i,j,k)+p*gvalue(i,stoj(k))*ci_o(j,k)*enew
     &              +q*gvalue_old(i,stoj(k))*ci_o_old(j,k)*eold
            de_o_ci(i,j,k)=de_o_ci(i,j,k)+p*denergy(i,k)*ci_o(j,k)
     &              +q*denergy_old(i,k)*ci_o_old(j,k)
            dj_de_ci(i,j,k)=dj_de_ci(i,j,k)+p*gvalue(i,stoj(k))*ci_de(j,k)
     &              +q*gvalue_old(i,stoj(k))*ci_de_old(j,k)

c            dj_o_ci(i,j,stobjx(k))=dj_o_ci(i,j,stobjx(k))+p*gvalue(i,stoj(k))*ci_o(j,stoo(k))
c     &              +q*gvalue_old(i,stoj(k))*ci_o_old(j,stoo(k))
c            dj_oe_ci(i,j,stobjx(k))=dj_oe_ci(i,j,stobjx(k))+p*gvalue(i,stoj(k))*ci_o(j,stoo(k))*enew
c     &              +q*gvalue_old(i,stoj(k))*ci_o_old(j,stoo(k))*eold
c            de_o_ci(i,j,stoo(k))=de_o_ci(i,j,stoo(k))+p*denergy(i,1)*ci_o(j,stoo(k))
c     &              +q*denergy_old(i,1)*ci_o_old(j,stoo(k))
c            dj_de_ci(i,j,stobjx(k))=dj_de_ci(i,j,stobjx(k))+p*gvalue(i,stoj(k))*ci_de(j,stoo(k))
c     &              +q*gvalue_old(i,stoj(k))*ci_de_old(j,stoo(k))

c            write(ounit,*) "state,stobjx(k),stoj(k),stoo(k),ci,jas,gn,cn,go,co",
c     &                      k, stobjx(k), stoj(k),stoo(k), j, i, gvalue(i,stoj(k)), ci_o(j,stoo(k)),
c     &                      gvalue_old(i,stoj(k)), ci_o_old(j,stoo(k))
c            write(ounit,*) "dj_o_ci(jas,ci,stobjx(state)), dj_oe_ci, de_o_ci(i,j,stoo(k)), dj_de_ci",
c     &                      dj_o_ci(i,j,stobjx(k)), dj_oe_ci(i,j,stobjx(k)),
c     &                      de_o_ci(i,j,stoo(k)), dj_de_ci(i,j,stobjx(k))

          enddo
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_jas_ci_init

      use mix_jas_ci, only: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci
      use optwf_control, only: ioptci, ioptjas
      use optwf_parms, only: nparmj
      use ci000, only: nciterm
      use vmc_mod, only: nwftypemax, stoo, stobjx
      use csfs, only: nstates

      implicit none

      integer :: i, j, k

      if(ioptjas.eq.0.or.ioptci.eq.0) return

      do k=1,nstates
        do i=1,nparmj
          do j=1,nciterm
            dj_o_ci(i,j,k)=0
            dj_oe_ci(i,j,k)=0
            de_o_ci(i,j,k)=0 !STU change to stobjx if we use jastrow dependent denergy, was stoo
            dj_de_ci(i,j,k)=0
          enddo
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_jas_ci_dump(iu)

      use mix_jas_ci, only: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci
      use optwf_control, only: ioptci, ioptjas
      use optwf_parms, only: nparmj
      use ci000, only: nciterm

      implicit none

      integer :: i, iu, j, k


      if(ioptjas.eq.0.or.ioptci.eq.0) return

      k=1 !STU need to set up for nwftypemax, nbjx.

      write(iu) ((dj_o_ci(i,j,k),dj_oe_ci(i,j,k),dj_de_ci(i,j,k),de_o_ci(i,j,k),i=1,nparmj),j=1,nciterm)

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_jas_ci_rstrt(iu)

      use mix_jas_ci, only: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci
      use optwf_control, only: ioptci, ioptjas
      use optwf_parms, only: nparmj
      use ci000, only: nciterm

      implicit none

      integer :: i, iu, j, k

      if(ioptjas.eq.0.or.ioptci.eq.0) return

      k=1 !STU need to set up for nwftypemax, nbjx

      read(iu) ((dj_o_ci(i,j,k),dj_oe_ci(i,j,k),dj_de_ci(i,j,k),de_o_ci(i,j,k),i=1,nparmj),j=1,nciterm)

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_jas_ci_fin(passes,eave)

      use optci, only: mxciterm
      use csfs, only: ccsf, ncsf
      use slater, only: cdet
      use gradhess_ci, only: grad_ci
      use gradhess_jas, only: grad_jas
      use gradhess_mix_jas_ci, only: h_mix_jas_ci, s_mix_jas_ci
      use mix_jas_ci, only: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci
      use optwf_control, only: ioptci, ioptjas
      use optwf_parms, only: nparmj
      use gradhessj, only: de, dj, dj_e
      use ci000, only: nciterm
      use ci005_blk, only: ci_o_cum
      use ci006_blk, only: ci_de_cum
      use ci008_blk, only: ci_oe_cum
      use optwf_control, only: method
      use precision_kinds, only: dp

      implicit none

      integer :: i, j, k
      real(dp) :: eave, h1, h2, passes
      real(dp), dimension(mxciterm) :: oelocav
      real(dp), dimension(mxciterm) :: eav

      if(ioptjas.eq.0.or.ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      k=1 !STU need to set up with nbjx, stobjx etc, not needed for sr_n

      if(method.eq.'hessian') then

c Compute mix Hessian
      do i=1,nparmj
        do j=1,nciterm
          h1=2*(2*(dj_oe_ci(i,j,k)-eave*dj_o_ci(i,j,k))-dj(i,1)*grad_ci(j)-grad_jas(i)*ci_o_cum(j,k))
          h2=de_o_ci(i,j,k)-de(i,1)*ci_o_cum(j,k)/passes
     &         +dj_de_ci(i,j,k)-dj(i,1)*ci_de_cum(j,k)/passes
          h_mix_jas_ci(i,j)=(h1+h2)/passes
        enddo
      enddo

      write(21,*) nciterm
      write(21,*) ((h_mix_jas_ci(i,j),j=1,nciterm),i=1,nparmj)

      elseif(method.eq.'linear') then

      if(ncsf.eq.0) then
        do i=1,nciterm
          oelocav(i)=0
          eav(i)=0
          do j=1,nciterm
            oelocav(i)=oelocav(i)+ci_oe_cum(i,j,k)*cdet(j,1,1)/passes
            eav(i)=eav(i)+ci_oe_cum(j,i,k)*cdet(j,1,1)/passes
          enddo
        enddo
       else
        do i=1,ncsf
          oelocav(i)=0
          eav(i)=0
          do j=1,ncsf
            oelocav(i)=oelocav(i)+ci_oe_cum(i,j,k)*ccsf(j,1,1)/passes
            eav(i)=eav(i)+ci_oe_cum(j,i,k)*ccsf(j,1,1)/passes
          enddo
        enddo
      endif

      do i=1,nparmj
        do j=1,nciterm
c Overlap s_jas_ci
          s_mix_jas_ci(i,j)=(dj_o_ci(i,j,k)-dj(i,1)*ci_o_cum(j,k)/passes)/passes
c H matrix h_jas_ci
          h_mix_jas_ci(i,j)=(dj_de_ci(i,j,k)+dj_oe_ci(i,j,k)
     &    +eave*dj(i,1)*ci_o_cum(j,k)/passes-dj(i,1)*eav(j)-ci_o_cum(j,k)*dj_e(i,1)/passes)/passes
c H matrix h_ci_jas
          h_mix_jas_ci(i+nparmj,j)=(de_o_ci(i,j,k)+dj_oe_ci(i,j,k)
     &    +eave*dj(i,1)*ci_o_cum(j,k)/passes-dj(i,1)*oelocav(j)-ci_o_cum(j,k)*(de(i,1)+dj_e(i,1))/passes)/passes
        enddo
      enddo

      endif

      return
      end
      end module
