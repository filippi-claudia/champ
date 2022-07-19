      module optx_orb_ci
      contains
      subroutine optx_orb_ci_sum(p,q)

      use optwf_control, only: ioptci, ioptorb
      use mix_orb_ci, only: ci_de_o, ci_o_ho, ci_o_o, ci_o_oe
      use orb_mat_001, only: orb_ho, orb_o, orb_oe
      use orb_mat_002, only: orb_ho_old, orb_o_old, orb_oe_old
      use ci000, only: nciterm
      use ci001_blk, only: ci_o
      use ci002_blk, only: ci_o_old
      use ci004_blk, only: ci_de, ci_de_old
      use optorb_cblock, only: nreduced
      use precision_kinds, only: dp
      use optwf_control, only: method

      implicit none

      integer :: i, j
      real(dp) :: p, q

      if(ioptorb.eq.0.or.ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do j=1,nreduced
       do i=1,nciterm
        ci_o_o(i,j)=ci_o_o(i,j)  +p*ci_o(i)*orb_o(j,1)+q*ci_o_old(i)*orb_o_old(j,1)
        ci_de_o(i,j)=ci_de_o(i,j)+p*ci_de(i)*orb_o(j,1)+q*ci_de_old(i)*orb_o_old(j,1)
        ci_o_ho(i,j)=ci_o_ho(i,j)+p*ci_o(i)*orb_ho(j,1)+q*ci_o_old(i)*orb_ho_old(j,1)
        ci_o_oe(i,j)=ci_o_oe(i,j)+p*ci_o(i)*orb_oe(j,1)+q*ci_o_old(i)*orb_oe_old(j,1)
       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_orb_ci_init

      use optwf_control, only: ioptci, ioptorb
      use mix_orb_ci, only: ci_de_o, ci_o_ho, ci_o_o, ci_o_oe
      use ci000, only: nciterm
      use optorb_cblock, only: nreduced
      use optwf_control, only: method

      implicit none

      integer :: i, j

      if(ioptorb.eq.0.or.ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do j=1,nreduced
        do i=1,nciterm
          ci_o_o(i,j)=0
          ci_o_oe(i,j)=0
          ci_o_ho(i,j)=0
          ci_de_o(i,j)=0
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_orb_ci_dump(iu)

      use optwf_control, only: ioptci, ioptorb
      use mix_orb_ci, only: ci_de_o, ci_o_ho, ci_o_o, ci_o_oe
      use ci000, only: nciterm
      use optorb_cblock, only: nreduced
      use optwf_control, only: method

      implicit none

      integer :: i, iu, j

      if(ioptorb.eq.0.or.ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return
      write(iu) ((ci_o_o(i,j),ci_o_oe(i,j),ci_o_ho(i,j),ci_de_o(i,j),i=1,nciterm),j=1,nreduced)

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_orb_ci_rstrt(iu)

      use optwf_control, only: ioptci, ioptorb
      use mix_orb_ci, only: ci_de_o, ci_o_ho, ci_o_o, ci_o_oe
      use ci000, only: nciterm
      use optorb_cblock, only: nreduced
      use optwf_control, only: method

      implicit none

      integer :: i, iu, j

      if(ioptorb.eq.0.or.ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return
      write(iu) ((ci_o_o(i,j),ci_o_oe(i,j),ci_o_ho(i,j),ci_de_o(i,j),i=1,nciterm),j=1,nreduced)

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_orb_ci_fin(passes,eave)

      use optci, only: mxciterm
      use csfs, only: ccsf, ncsf
      use gradhess_ci, only: grad_ci
      use gradhess_mix_orb_ci, only: h_mix_ci_orb, s_mix_ci_orb
      use optwf_control, only: ioptci, ioptorb
      use optwf_parms, only: nparmj
      use optorb_cblock, only: norbprim
      use mix_orb_ci, only: ci_de_o, ci_o_ho, ci_o_o, ci_o_oe
      use orb_mat_003, only: orb_o_cum
      use orb_mat_004, only: orb_oe_cum
      use orb_mat_005, only: orb_ho_cum
      use gradhess_all, only: grad
      use ci000, only: nciterm
      use ci005_blk, only: ci_o_cum
      use ci006_blk, only: ci_de_cum
      use ci008_blk, only: ci_oe_cum
      use optorb_cblock, only: nreduced
      use precision_kinds, only: dp
      use optwf_control, only: method
      use slater, only: cdet

      implicit none

      integer :: i, ishift, j
      real(dp) :: eave, h1, h2, passes
      real(dp), dimension(mxciterm) :: oelocav
      real(dp), dimension(mxciterm) :: eav


c     common /gradhess_orb/ grad_orb(norbterm),h_orb(MXMATDIM),s_orb(MXMATDIM)


      if(ioptorb.eq.0.or.ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      if(method.eq.'hessian') then

      ishift=nparmj+nciterm-1
      do i=1,nciterm
        do j=1,norbprim
          h1=2*(ci_o_oe(i,j)-eave*ci_o_o(i,j)-ci_o_cum(i)*grad(j+ishift)-grad_ci(i)*orb_o_cum(j,1))
          h2=2*(ci_de_o(i,j)-ci_de_cum(i)*orb_o_cum(j,1)/passes)
          h_mix_ci_orb(i,j)=(h1+h2)/passes
        enddo
      enddo

      elseif(method.eq.'linear') then

      if(ncsf.eq.0) then
        do i=1,nciterm
          oelocav(i)=0
          eav(i)=0
          do j=1,nciterm
            oelocav(i)=oelocav(i)+ci_oe_cum(i,j)*cdet(j,1,1)/passes
            eav(i)=eav(i)+ci_oe_cum(j,i)*cdet(j,1,1)/passes
          enddo
        enddo
       else
        do i=1,ncsf
          oelocav(i)=0
          eav(i)=0
          do j=1,ncsf
            oelocav(i)=oelocav(i)+ci_oe_cum(i,j)*ccsf(j,1,1)/passes
            eav(i)=eav(i)+ci_oe_cum(j,i)*ccsf(j,1,1)/passes
          enddo
        enddo
      endif

      do i=1,nciterm
        do j=1,nreduced
c Overlap s_ij
          s_mix_ci_orb(i,j)=(ci_o_o(i,j)-ci_o_cum(i)*orb_o_cum(j,1)/passes)/passes
c Hamiltonian ci_orb
          h_mix_ci_orb(i,j)=(ci_o_ho(i,j)+eave*ci_o_cum(i)*orb_o_cum(j,1)/passes
     &    -ci_o_cum(i)*orb_ho_cum(j,1)/passes-orb_o_cum(j,1)*oelocav(i))/passes
c Hamiltonian orb_ci
          h_mix_ci_orb(i+nciterm,j)=(ci_de_o(i,j)+ci_o_oe(i,j)+eave*ci_o_cum(i)*orb_o_cum(j,1)/passes
     &    -ci_o_cum(i)*orb_oe_cum(j,1)/passes-orb_o_cum(j,1)*eav(i))/passes
        enddo
      enddo

      endif

      return
      end
      end module
