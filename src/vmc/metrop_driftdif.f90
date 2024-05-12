module metrop_driftdif
contains
      subroutine metrop1_moveall(ipass,irun)

      use acuest_mod, only: acues1,acusig
      use config,  only: eold,psido,psijo,psi2n,psi2o
      use config,  only: vnew,vold,xnew,xold
      use contrl_file, only: ounit
      use control, only: ipr
      use determinante_mod, only: compute_determinante_grad
      use estsum,  only: acc,esum,esum1,pesum,tpbsum
      use error,   only: fatal_error
      use forcewt, only: wsum
      use gauss_mod, only: gauss
      use hpsi_mod, only: hpsi
      use metropolis, only: vmc_tau
      use mstates_ctrl, only: iguiding
      use multiple_geo, only: nforce, MFORCE
      use optci_mod, only: optci_sum, optci_save
      use optjas_mod, only: optjas_sum, optjas_save
      use optorb_f_mod, only: check_orbitals,check_orbitals_reset
      use optorb_f_mod, only: optorb_sum, optorb_save
      use optx_jas_ci, only: optx_jas_ci_sum
      use optx_jas_orb, only: optx_jas_orb_sum
      use optx_orb_ci, only: optx_orb_ci_sum
      use pcm_vmc, only: pcm_sum, pcm_save
      use precision_kinds, only: dp
      use prop_vmc, only: prop_sum, prop_save
      use random_mod, only: random_dp
      use strech_mod, only: strech
      use system,  only: nelec

      implicit none

      integer :: i, ii, ifr, ipass, irun, istate, k 
      real(dp) :: ajacob, dfus, dfus2o, dfus2n, drift, dx
      real(dp) :: p, q, rttau, tau, wstrn, wstro
      real(dp), dimension(1) :: ekinn, ekino
      real(dp), dimension(1) :: psidn, psijn
      real(dp), dimension(1) :: p_dum, q_dum
      real(dp), dimension(MFORCE)  :: enew
      real(dp), dimension(3) :: xbac
      real(dp), dimension(3,nelec) :: xstrech
      real(dp), dimension(1) :: zero_array = 0.0_dp
      real(dp), dimension(1) :: one_array = 1.0_dp
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0

      if(iguiding.gt.0) call fatal_error ('Guiding wave function only in mov1')

      tau=vmc_tau
      rttau=sqrt(tau)

      call check_orbitals

! Start loop over electrons
      dfus2o=zero
      do i=1,nelec
        !call compute_determinante_grad(i,psido(1),psido,psijo,vold(1,i),1)

        do k=1,3
          drift=vold(k,i)*tau
          dfus=gauss()*rttau
          dfus2o=dfus2o+dfus**2
          xnew(k,i)=xold(k,i)+drift+dfus
        enddo

        if(ipr.ge.1) then
          write(ounit,'(''xold='',9d12.4)') (xold(k,i),k=1,3)
          write(ounit,'(''vold='',9d12.4)') (vold(k,i),k=1,3)
          write(ounit,'(''xnew='',9d12.4)') (xnew(k,i),k=1,3)
          write(ounit,*)
        endif
      enddo

! calculate psi etc. at new configuration

! loop over secondary configurations
      do ifr=2,nforce
        call strech(xnew,xstrech,ajacob,ifr,1)
        call hpsi(xstrech,psidn,psijn,ekinn,enew(ifr),ipass,ifr)
        psi2n(ifr)=2*(dlog(dabs(psidn(1)))+psijn(1))+dlog(ajacob)
      enddo

      call check_orbitals_reset

! primary configuration
      if(nforce.gt.1) call strech(xnew,xstrech,ajacob,1,0)
      call hpsi(xnew,psidn,psijn,ekinn,enew(1),ipass,1)
      psi2n(1)=2*(dlog(dabs(psidn(1)))+psijn(1))

      do i=1,nelec
        call compute_determinante_grad(i,psidn(1),psidn,psijn,vnew(1,i),1)
      enddo

! calculate probability for reverse transition
      dfus2n=0.d0
      do i=1,nelec
        do k=1,3
          drift=vnew(k,i)*tau
          xbac(k)=xnew(k,i)+drift
          dfus=xbac(k)-xold(k,i)
          dfus2n=dfus2n+dfus**2
        enddo
      enddo

! p is the probability of accepting new move
      p=exp(psi2n(1)-psi2o(1,1))*exp(-0.5*(dfus2n-dfus2o)/tau)

      if(ipr.ge.1) write(ounit,'('',p'',d12.4)') p

      p=dmin1(one,p)
      q=one-p

! form expected values of e, pe, etc.
  208 esum1(1)=           p*enew(1)+q*eold(1,1)
      esum(1,1)=esum(1,1)+p*enew(1)+q*eold(1,1)
      wsum(1,1)=wsum(1,1)+1
      pesum(1)=pesum(1)+p*(enew(1)-ekinn(1))+q*(eold(1,1)-ekino(1))
      tpbsum(1)=tpbsum(1)+p*ekinn(1)+q*ekino(1)

      call pcm_sum(p,q)
      call prop_sum(p,q)

      p_dum(1)=p
      q_dum(1)=q
      call optjas_sum(p_dum,q_dum,enew,eold,0)
      call optorb_sum(p_dum,q_dum,enew,eold,0)
      call optci_sum(p,q,enew(1),eold(1,1))

      call acues1(one_array)

      do ifr=2,nforce
        wstro=exp(psi2o(1,ifr)-psi2o(1,1))
        wstrn=exp(psi2n(ifr)-psi2n(1))
        esum(1,ifr)=esum(1,ifr)+p*enew(ifr)*wstrn+q*eold(1,ifr)*wstro
        wsum(1,ifr)=wsum(1,ifr)+p*wstrn+q*wstro
      enddo

! accept new move with probability p
      acc=acc+p
      if (random_dp().lt.p) then
! move is accepted so update positions etc.
        do i=1,nelec
          do k=1,3
            xold(k,i)=xnew(k,i)
            vold(k,i)=vnew(k,i)
          enddo
        enddo
        do ifr=1,nforce
          eold(1,ifr)=enew(ifr)
          psi2o(1,ifr)=psi2n(ifr)
        enddo
        ekino(1)=ekinn(1)

        call optjas_save
        call optorb_save
        call optci_save
        call pcm_save
        call prop_save

        if(ipr.ge.1) write(ounit,*)'METROP ACCEPT'
       else
        if(ipr.ge.1) write(ounit,*)'METROP REJECT'
      endif

      esum1(1)=eold(1,1)
      call acusig(one_array)

      return
      end
end module metrop_driftdif
