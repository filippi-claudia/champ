module nonloc_grid_mod
contains
      subroutine nonloc_grid(iel,iw,x,psid,imove,t_norm,ipropose)

      use casula,  only: icasula,t_vpsp
      use config,  only: xold_dmc
      use contrl_file, only: ounit
      use contrl_per, only: iperiodic
      use contrldmc, only: tau
      use distance_mod, only: r_en,rvec_en
      use distances_mod, only: distances
      use multideterminant_tmove_mod, only: multideterminant_tmove
      use multiple_geo, only: iwf,iwftype
      use nonloc_pot_mod, only: nonloc_pot
      use optwf_control, only: ioptci,ioptjas,ioptorb
      use optwf_parms, only: nparmj
      use precision_kinds, only: dp
      use qua,     only: nquad,xq,yq,zq
      use random_mod, only: random_dp
      use system,  only: cent,ncent,nelec
      use vmc_mod, only: nbjx
      implicit none

      integer :: i, i1, i2, ic, ic_good
      integer :: iel, iel_good, ii, imove
      integer :: ioptci_sav, ioptjas_sav, ioptorb_sav, iq
      integer :: iq_good, iw, ipropose
      real(dp) :: costh, p, pe, psid
      real(dp) :: psidi, ri, t_cum
      real(dp) :: t_norm, t_normi, tauprim
      real(dp), dimension(2,nbjx) :: vpsp_det
      real(dp), dimension(nparmj,nbjx) :: dvpsp_dj
      real(dp), dimension(*) :: x
      real(dp), parameter :: one = 1.0_dp !was uninit

! here vpsp_det and dvpsp_det are dummy

      iwf=iwftype(1)

      tauprim=tau
      if(icasula.gt.0)then
        call distances(iel,xold_dmc(1,1,iw,1))

        ioptjas_sav=ioptjas
        ioptorb_sav=ioptorb
        ioptci_sav=ioptci
        ioptjas=0
        ioptorb=0
        ioptci=0

        call nonloc_pot(xold_dmc(1,1,iw,1),rvec_en,r_en,pe,vpsp_det,dvpsp_dj,t_vpsp,iel,1)

        call multideterminant_tmove(psid,iel)

        ioptjas=ioptjas_sav
        ioptorb=ioptorb_sav
        ioptci=ioptci_sav

        i1=iel
        i2=iel
       else

        i1=1
        i2=nelec
      endif
      imove=0

!     do i=i1,i2
!     write(ounit,'(''t_vpsp before = '',100e10.2)') ((t_vpsp(ic,iq,i),ic=1,ncent),iq=1,nquad)
!     enddo

      t_norm=0.d0
      psidi=1.d0/psid
      do ii=i1,i2
        i=ii
        if(i.gt.nelec) i=i-nelec
        if(i.lt.1) i=i+nelec
        do iq=1,nquad
          do ic=1,ncent
            t_vpsp(ic,iq,i)=t_vpsp(ic,iq,i)
            if(t_vpsp(ic,iq,i).gt.0.d0) t_vpsp(ic,iq,i)=0.d0
            t_norm=t_norm-t_vpsp(ic,iq,i)
          enddo
        enddo
      enddo
      t_norm=1.d0+t_norm*tauprim
      t_normi=1.d0/t_norm
!     write(ounit,*) 'tnormi=',t_normi
!     do i=i1,i2
!     write(ounit,'(''t_vpsp after = '',100f14.6)') ((-tauprim*t_vpsp(ic,iq,i)*t_normi,ic=1,ncent),iq=1,nquad)
!     enddo
      if (ipropose.eq.0) return
      if(t_norm.eq.1.d0) return

      t_cum=0.d0
      p=random_dp()
      do ii=i1,i2
        i=ii
        if(i.gt.nelec) i=i-nelec
        if(i.lt.1) i=i+nelec
        do iq=1,nquad
          do ic=1,ncent
            t_cum=t_cum-tauprim*t_vpsp(ic,iq,i)*t_normi
            if(t_cum.gt.p)then
              ic_good=ic
              iq_good=iq
              iel_good=i
              imove=1
              go to 30
            endif
          enddo
        enddo
      enddo

      30 if(imove.eq.1)then
        iq=iq_good
        ic=ic_good
        iel=iel_good
        if(icasula.lt.0) call distances(iel,xold_dmc(1,1,iw,1))
        ri=one/r_en(iel,ic)
        costh=rvec_en(1,iel,ic)*xq(iq) &
             +rvec_en(2,iel,ic)*yq(iq) &
             +rvec_en(3,iel,ic)*zq(iq)
        costh=costh*ri

        x(1)=r_en(iel,ic)*xq(iq)+cent(1,ic)
        x(2)=r_en(iel,ic)*yq(iq)+cent(2,ic)
        x(3)=r_en(iel,ic)*zq(iq)+cent(3,ic)

!       write(ounit,*) 'moved B',iw,iel,(xold_dmc(kk,iel,iw,1),kk=1,3)
!       write(ounit,*) 'moved A',iw,iel,(x(kk),kk=1,3)
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine t_vpsp_sav


      use casula,  only: t_vpsp
      use precision_kinds, only: dp
      use pseudo_mod, only: MPS_QUAD
      use qua,     only: nquad
      use system,  only: ncent,ncent_tot,nelec
      use save_mod, only: t_vpsp_save

      implicit none

      integer :: i, ic, iq


      if (.not.allocated(t_vpsp_save)) allocate(t_vpsp_save(ncent_tot, MPS_QUAD, nelec))

! real(dp), dimension(ncent_tot, MPS_QUAD, nelec) :: t_vpsp_save
! save t_vpsp_save

      do i=1,nelec
        do iq=1,nquad
          do ic=1,ncent
            t_vpsp_save(ic,iq,i)=t_vpsp(ic,iq,i)
          enddo
        enddo
      enddo

      end subroutine

      subroutine t_vpsp_get
      use casula,  only: t_vpsp
      use qua,     only: nquad
      use system,  only: ncent,nelec
      use save_mod, only: t_vpsp_save

      implicit none

      integer :: i, ic, iq

      do i=1,nelec
        do iq=1,nquad
          do ic=1,ncent
            t_vpsp(ic,iq,i)=t_vpsp_save(ic,iq,i)
          enddo
        enddo
      enddo

      return
      end
end module
