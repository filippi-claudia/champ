      module optx_jas_orb
      contains
      subroutine optx_jas_orb_sum(wtg_new,wtg_old,iflag)

      use csfs, only: nstates
      use derivjas, only: gvalue
      use gradhessjo, only: denergy_old, gvalue_old
      use optwf_control, only: ioptjas, ioptorb
      use optwf_parms, only: nparmj
      use deloc_dj_m, only: denergy
      use mix_jas_orb, only: de_o, dj_ho, dj_o, dj_oe
      use orb_mat_001, only: orb_ho, orb_o, orb_oe
      use orb_mat_002, only: orb_ho_old, orb_o_old, orb_oe_old
      use optwf_control, only: method
      use optorb_cblock, only: nreduced
      use precision_kinds, only: dp
      use vmc_mod, only: stoj

      implicit none

      integer :: i, iflag, istate, j, js
      real(dp) :: p, q
      real(dp), dimension(*) :: wtg_new
      real(dp), dimension(*) :: wtg_old

      if(ioptjas.eq.0.or.ioptorb.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do istate=1,nstates
      !STU must add jastrow state mapping here for gvalue and gvalue_old
      ! although this is not called during sr_n apparently. maybe
      ! istate=1 in gvalue.
      js=stoj(istate)
      p=wtg_new(istate)
      do i=1,nparmj
        do j=1,nreduced
        dj_o(i,j,istate)=dj_o(i,j,istate)  +p*gvalue(i,js)*orb_o(j,istate)
        dj_oe(i,j,istate)=dj_oe(i,j,istate)+p*gvalue(i,js)*orb_oe(j,istate)
        dj_ho(i,j,istate)=dj_ho(i,j,istate)+p*gvalue(i,js)*orb_ho(j,istate)
        de_o(i,j,istate)=de_o(i,j,istate)  +p*denergy(i,istate)*orb_o(j,istate)
        enddo
      enddo

      enddo

      if(iflag.eq.0) return

      do istate=1,nstates
      js=stoj(istate)

      q=wtg_old(istate)
      do i=1,nparmj
        do j=1,nreduced
        dj_o(i,j,istate)=dj_o(i,j,istate)  +q*gvalue_old(i,js)*orb_o_old(j,istate)
        dj_oe(i,j,istate)=dj_oe(i,j,istate)+q*gvalue_old(i,js)*orb_oe_old(j,istate)
        dj_ho(i,j,istate)=dj_ho(i,j,istate)+q*gvalue_old(i,js)*orb_ho_old(j,istate)
        de_o(i,j,istate)=de_o(i,j,istate)  +q*denergy_old(i,istate)*orb_o_old(j,istate)
        enddo
      enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_jas_orb_init

      use csfs, only: nstates
      use optwf_control, only: ioptjas, ioptorb
      use optwf_parms, only: nparmj
      use mix_jas_orb, only: de_o, dj_ho, dj_o, dj_oe
      use optwf_control, only: method
      use optorb_cblock, only: nreduced

      implicit none

      integer :: i, istate, j




      if(ioptjas.eq.0.or.ioptorb.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do istate=1,nstates

      do i=1,nparmj
        do j=1,nreduced
          dj_o(i,j,istate)=0
          dj_oe(i,j,istate)=0
          dj_ho(i,j,istate)=0
          de_o(i,j,istate)=0
        enddo
      enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_jas_orb_dump(iu)

      use csfs, only: nstates
      use optwf_control, only: ioptjas, ioptorb
      use optwf_parms, only: nparmj
      use mix_jas_orb, only: de_o, dj_ho, dj_o, dj_oe
      use optwf_control, only: method
      use optorb_cblock, only: nreduced

      implicit none

      integer :: i, istate, iu, j


      if(ioptjas.eq.0.or.ioptorb.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do istate=1,nstates
      write(iu) ((dj_o(i,j,istate),dj_oe(i,j,istate),dj_ho(i,j,istate),de_o(i,j,istate),i=1,nparmj),j=1,nreduced)
      enddo


      return
      end
c-----------------------------------------------------------------------
      subroutine optx_jas_orb_rstrt(iu)

      use csfs, only: nstates
      use optwf_control, only: ioptjas, ioptorb
      use optwf_parms, only: nparmj
      use mix_jas_orb, only: de_o, dj_ho, dj_o, dj_oe
      use optwf_control, only: method
      use optorb_cblock, only: nreduced

      implicit none

      integer :: i, istate, iu, j


      if(ioptjas.eq.0.or.ioptorb.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do istate=1,nstates
      read(iu) ((dj_o(i,j,istate),dj_oe(i,j,istate),dj_ho(i,j,istate),de_o(i,j,istate),i=1,nparmj),j=1,nreduced)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine optx_jas_orb_fin(wcum,ecum)

      use optorb_mod, only: mxreduced
      use csfs, only: nstates
      use gradhess_mix_jas_orb, only: h_mix_jas_orb, s_mix_jas_orb
      use optwf_control, only: ioptjas, ioptorb, iuse_orbeigv, iapprox
      use optwf_parms, only: nparmj
      use sa_weights, only: weights
      use gradhessj, only: de, dj, dj_e
      use mix_jas_orb, only: de_o, dj_ho, dj_o, dj_oe
      use orb_mat_003, only: orb_o_cum
      use orb_mat_004, only: orb_oe_cum
      use orb_mat_005, only: orb_ho_cum
      use optwf_control, only: method
      use optorb_cblock, only: nreduced
      ! I think this one is not needed ...
      ! use gradhess_jas, only: grad_jas
      use precision_kinds, only: dp
      use error, only: fatal_error
      implicit none

      integer :: i, istate, j
      real(dp) :: eave, grad_jas, h1, h2, passes
      real(dp) :: passesi, wts
      real(dp), dimension(*) :: wcum
      real(dp), dimension(*) :: ecum
      real(dp), dimension(mxreduced) :: grad_orb


      if(ioptjas.eq.0.or.ioptorb.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do i=1,nparmj
        do j=1,nreduced
          s_mix_jas_orb(i,j)=0
          h_mix_jas_orb(i,j)=0
          h_mix_jas_orb(i+nparmj,j)=0
        enddo
      enddo

c Hessian method
      if(method.eq.'hessian') then

c Exact mixed Hessian terms not implemented
      if(iuse_orbeigv.ne.0) then
        call fatal_error('OPTX_JAS_ORB_FIN: exact mix hessian formulas not implemented')
       else
c Approximate mixed Hessian terms using orbital eigenvalues
        do istate=1,nstates

        wts=weights(istate)

        passes=wcum(istate)
        passesi=1/passes
        eave=ecum(istate)*passesi

        do j=1,nreduced
          grad_orb(j)=2*(orb_oe_cum(j,istate)-eave*orb_o_cum(j,istate))*passesi
        enddo

        do i=1,nparmj
          grad_jas=2*(dj_e(i,istate)-eave*dj(i,istate))*passesi
          do j=1,nreduced
            h1=2*(2*(dj_oe(i,j,istate)-eave*dj_o(i,j,istate))-dj(i,istate)*grad_orb(j)-grad_jas*orb_o_cum(j,istate))
            h2=2*(de_o(i,j,istate)-de(i,istate)*orb_o_cum(j,istate)*passesi)
            h_mix_jas_orb(i,j)=h_mix_jas_orb(i,j)+wts*(h1+h2)*passesi
          enddo
        enddo

        enddo

      endif

c Linear method
      elseif(method.eq.'linear') then

c Exact Hamiltonian mixed terms on semi-orthogonal basis
      if(iuse_orbeigv.eq.0) then

      do istate=1,nstates

      wts=weights(istate)

      passes=wcum(istate)
      passesi=1/passes
      eave=ecum(istate)*passesi

      do i=1,nparmj
        do j=1,nreduced
c Overlap jas_orb
          s_mix_jas_orb(i,j)=s_mix_jas_orb(i,j)+wts*(dj_o(i,j,istate)-dj(i,istate)*orb_o_cum(j,istate)*passesi)*passesi
c Hamiltonian jas_orb
          h_mix_jas_orb(i,j)=h_mix_jas_orb(i,j)+wts*(dj_ho(i,j,istate)+
     &    (eave*dj(i,istate)*orb_o_cum(j,istate)-dj(i,istate)*orb_ho_cum(j,istate)-
     &    orb_o_cum(j,istate)*dj_e(i,istate))*passesi)*passesi
c Hamiltonian orb_jas
          h_mix_jas_orb(i+nparmj,j)=h_mix_jas_orb(i+nparmj,j)+wts*(de_o(i,j,istate)+dj_oe(i,j,istate)+
     &    (eave*dj(i,istate)*orb_o_cum(j,istate)-dj(i,istate)*orb_oe_cum(j,istate)-
     &    orb_o_cum(j,istate)*(de(i,istate)+dj_e(i,istate)))*passesi)*passesi

        enddo
      enddo

      enddo

      endif

      endif

      if(method.eq.'linear') then
      if(iapprox.eq.1) then
        do i=1,nparmj
          do j=1,nreduced
            h_mix_jas_orb(i,j)=h_mix_jas_orb(i+nparmj,j)
          enddo
        enddo
       elseif(iapprox.eq.2) then
        do i=1,nparmj
          do j=1,nreduced
            s_mix_jas_orb(i,j)=0
            h_mix_jas_orb(i,j)=0
            h_mix_jas_orb(i+nparmj,j)=0
          enddo
        enddo
       elseif(iapprox.eq.3) then
        do i=1,nparmj
          do j=1,nreduced
            h_mix_jas_orb(i,j)=0.5*(h_mix_jas_orb(i,j)+h_mix_jas_orb(i+nparmj,j))
            h_mix_jas_orb(i+nparmj,j)=h_mix_jas_orb(i,j)
          enddo
        enddo
      endif

      endif

      return
      end
      end module
