module metrop_mov1_slat
contains
      subroutine metrop6(ipass,irun)
! Written by Cyrus Umrigar
! Uses the accelerated Metropolis method described in:
! 1) Accelerated Metropolis Method, C.J. Umrigar, PRL 71, 408 (1993).
! 2) Variational Monte Carlo Basics and Applications to Atoms and Molecules,
!    C.J. Umrigar, in "Quantum Monte Carlo Methods in Physics and Chemistry",
!    edited by M.P. Nightingale and C.J. Umrigar. NATO ASI Series, Series C,
!    Mathematical and Physical Sciences, Vol. C-525,
!    (Kluwer Academic Publishers, Boston, 1999)
      use acuest_mod, only: acues1,acusig
      use config, only: delttn,eold,nearestn,nearesto,psi2n,psi2o
      use config, only: psido,psijo,rminn,rminno,rmino,rminon,rvminn
      use config, only: rvminno,rvmino,rvminon,vnew,vold,xnew
      use config, only: xold
      use constants, only: pi
      use control, only: ipr
      use contrl_file, only: ounit
      use csfs,    only: nstates
      use determinant_psig_mod, only: determinant_psig
      use determinante_mod, only: compute_determinante_grad
      use detsav_mod, only: detsav
      use distances_mod, only: distancese_restore
      use estsum,  only: acc,esum,esum1,pesum,tpbsum
      use force_analytic, only: force_analy_sum
      use forcewt, only: wsum
      use gammai_mod, only: gammai
      use hpsi_mod, only: hpsi
      use hpsie,   only: psie
      use inputflags, only: eps_node_cutoff,node_cutoff
      use jassav_mod, only: jassav
      use metropolis, only: deltar,deltat,fbias
      use mmpol,   only: mmpol_efield
      use mmpol_cntrl, only: ich_mmpol
      use mmpol_vmc, only: mmpol_sum
      use mstates_ctrl, only: iguiding
      use mstates_mod, only: MSTATES
      use multideterminant_mod, only: update_ymat
      use multiple_geo, only: nforce
      use multiple_states, only: efficiency_sample
      use nodes_distance_mod, only: nodes_distance,rnorm_nodes_num
      use optci_mod, only: optci_sum
      use optjas_mod, only: optjas_sum
      use optorb_f_mod, only: check_orbitals,check_orbitals_reset
      use optorb_f_mod, only: optorb_sum
      use optwf_handle_wf, only: optwf_store
      use optx_jas_ci, only: optx_jas_ci_sum
      use optx_jas_orb, only: optx_jas_orb_sum
      use optx_orb_ci, only: optx_orb_ci_sum
      use pcm_cntrl, only: ichpol
      use pcm_mod, only: qpcm_efield
      use pcm_vmc, only: pcm_sum
      use precision_kinds, only: dp
      use prop_vmc, only: prop_sum
      use pseudo,  only: nloc
      use random_mod, only: random_dp
      use stats,   only: rejmax
      use strech_mod, only: strech
      use system,  only: cent,iwctype,ncent,nelec,nup,znuc
      use tmpnode, only: distance_node_sum
      use vmc_mod, only: nwftypejas, stoj
      use mpitimer, only: elapsed_time
      use mpi
      use optwf_control, only: ioptjas

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
      use qmckl_data
#endif
      implicit none 

      integer :: i, iab, ic, iel, iflag_dn, rc
      integer :: iflag_up, iflagb, iflagt, iflagz
      integer :: ifr, igeometrical, ii, ipass
      integer :: irun, istate, ierr
      integer :: j, jel, k, nearn, nearo
      real(dp) :: ajacob, arean, areao, bot
      real(dp) :: clim, co, cosphi, costht
      real(dp) :: deltri, deltt
      real(dp) :: dist, distance_node, dmin1, dot
      real(dp) :: fmax, fmax2, fxnp, fxop
      real(dp) :: g32dif, g32dif1, g32dif2, g52bot
      real(dp) :: g52dif, g52dif1, g52dif2, g52top
      real(dp) :: g52zer, p, phitry, phizer
      real(dp) :: psidg, psig, q
      real(dp) :: r, ratio, raver, ravern
      real(dp) :: rbot, rmax1, rmax2 = 0d0, rnew
      real(dp) :: rnorm, rnorm_nodes, rold, root
      real(dp) :: rratio, rtest, rtest2, rtop
      real(dp) :: rtry, rzero, sintht, term
      real(dp) :: term2, top, vnewp
      real(dp) :: vnewr, voldp, voldr, wstro
      real(dp) :: xprime, yprime, z, zcusp
      real(dp) :: zebot, zeta, zetop, zezer
      real(dp) :: zprime, zrbot, zrtop, zrzer
      real(dp), dimension(3,nelec) :: xstrech
      real(dp), dimension(3) :: xaxis
      real(dp), dimension(3) :: yaxis
      real(dp), dimension(3) :: zaxis
      real(dp), dimension(3) :: ddx_ref
      real(dp), dimension(MSTATES) :: ekino
      real(dp), dimension(MSTATES) :: psidn
      real(dp), dimension(nwftypejas) :: psijn
      real(dp), dimension(MSTATES) :: wtg
      real(dp), dimension(MSTATES) :: wtg_sqrt
      real(dp), dimension(MSTATES) :: zero_array = 0.0_dp
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: four = 4.d0
      real(dp), parameter :: half = 0.5d0
      real(dp), parameter :: d3b2 = 1.5d0
      real(dp), parameter :: d5b2 = 2.5d0
      real(dp), parameter :: d2b3 = .666666666666667d0
      real(dp), parameter :: eps = 1.d-10
      real(dp), parameter :: g5b2 = 1.329340388179137d0

! g3b2, g5b2 are gamma3/2), gamma(5/2)

! The moves are now being made in local r,theta phi coordinates.

! The foll. additions have been made:
! 1) Slater form of Tij.
! 2) Make theta_max a function of r
! 3) Generalize to molecules. This requires geometric rejections.

      deltri=one/deltar

      call check_orbitals

      do i=1,nelec

        if(i.le.nup) then
          iab=1
          iflag_up=2
          iflag_dn=3
         else
          iab=2
          iflag_up=3
          iflag_dn=2
        endif

        if(iguiding.eq.0) then
          psidg=psido(1)
          psig=psido(1)*exp(psijo(1))
        else
          call determinant_psig(psido,psijo,psig)
        endif
        call compute_determinante_grad(i,psig,psido,psijo,vold(1,i),1)
        fxop=one
        nearo=nearesto(i)
        if(nloc.eq.0) then
          zcusp=znuc(iwctype(nearo))
        else
          zcusp=zero
        endif
! Although rmino is saved, recalculate it, otherwise there may be a
! build up of errors via zaxis. I don't think so.
        rmino(i)=dsqrt(rvmino(1,i)**2+rvmino(2,i)**2+rvmino(3,i)**2)
! Choose lower and upper values of r sampling
        rbot=rmino(i)*deltri
        rtop=rmino(i)*deltar
! Calculate magnitude of the velocity in the radial direction
        voldr=zero
        do ic=1,3
          voldr=voldr+vold(ic,i)*rvmino(ic,i)
        enddo

        voldr=voldr/rmino(i)

! Place x-axis along direction of angular change and
! Calculate the velocity in the phi direction
        voldp=zero
        do ic=1,3
          zaxis(ic)=rvmino(ic,i)/rmino(i)
          xaxis(ic)=vold(ic,i)-voldr*zaxis(ic)
          voldp=voldp+xaxis(ic)**2
        enddo
        voldp=dsqrt(voldp)
        if(voldp.lt.eps) then
          xaxis(1)=eps*(one-zaxis(1)**2)
          xaxis(2)=eps*(-zaxis(1)*zaxis(2))
          xaxis(3)=eps*(-zaxis(1)*zaxis(3))
          voldp=eps*sqrt(one-zaxis(1)**2)
        endif
        do ic=1,3
          xaxis(ic)=xaxis(ic)/voldp
        enddo

! Limit radial component of velocity.
! It may be a good idea to limit it if it is positive too.
        voldr=max(voldr,-2*znuc(iwctype(nearo)))
!       voldr=min(voldr,2*znuc(iwctype(nearo)))

! y-axis is cross-product of z and x axes
        yaxis(1)=zaxis(2)*xaxis(3)-zaxis(3)*xaxis(2)
        yaxis(2)=zaxis(3)*xaxis(1)-zaxis(1)*xaxis(3)
        yaxis(3)=zaxis(1)*xaxis(2)-zaxis(2)*xaxis(1)

! Temporary test of fbias
!       voldr=voldr*fbias
        voldp=voldp*fbias

        root=(zcusp+voldr)*(zcusp+voldr-four/rmino(i))
        if(root.ge.zero) then
          root=sqrt(root)
          zeta=half*(zcusp-voldr-root)
          if(zeta.le.zero) zeta=half*(zcusp-voldr+root)
          if(zeta.le.zero) zeta=eps
         else
          if(voldr.lt.zero) then
            zeta=-voldr-eps
           else
            zeta=one
          endif
        endif
        co=(zeta+voldr)/(one-(zeta+voldr)*rmino(i))

! Use Slater approx for radial fn
! Determine the maximum value of radial function for rejection sampling
        root=dsqrt((d3b2*co-zeta)**2+two*zeta*co)
        rmax1=(d3b2*co-zeta+root)/(two*zeta*co)
        if(rmax1.lt.rbot) rmax1=rbot
        if(rmax1.gt.rtop) rmax1=rtop
        fmax=sqrt(rmax1)*abs(one+co*rmax1)*dexp(-zeta*rmax1)
        if(co.lt.zero) then
          rmax2=(d3b2*co-zeta-root)/(two*zeta*co)
          if(rmax2.lt.rbot) rmax2=rbot
          if(rmax2.gt.rtop) rmax2=rtop
          fmax2=sqrt(rmax2)*abs(one+co*rmax2)*dexp(-zeta*rmax2)
          fmax=max(fmax,fmax2)
        endif

!   Sample sqrt(r_f)*abs(1+co*r_f)*exp(-zeta*r_f) by rejection
        bot=sqrt(rmino(i))*abs(one+co*rmino(i))*dexp(-zeta*rmino(i))
      40   rtry=((deltar-deltri)*random_dp()+deltri)*rmino(i)
          top=sqrt(rtry)*abs(one+co*rtry)*dexp(-zeta*rtry)
          ratio=top/fmax
          rejmax=max(rejmax,ratio)
          if(ratio.gt.random_dp()) goto 50
        goto 40
      50   fxop=fxop*top/bot

!   Calculate the integral of T
        rzero=-one/co
        zrbot=zeta*rbot
        zrtop=zeta*rtop
        zebot=dsqrt(zrbot**3)*dexp(-zrbot)
        zetop=dsqrt(zrtop**3)*dexp(-zrtop)
        g52bot=gammai(d5b2,zrbot,zrbot*zebot,iflagb)
        g52top=gammai(d5b2,zrtop,zrtop*zetop,iflagt)
        if(rzero.lt.rbot .or. rzero.gt.rtop) then
          if(iflagb*iflagt.eq.1) then
            g52dif=g52top-g52bot
           else
            g52dif=g52top-g52bot+g5b2
          endif
          g32dif=d2b3*(g52dif+zetop-zebot)
          areao=dabs(g32dif+co*g52dif/zeta)/(bot*dsqrt(zeta**3))
         else
          zrzer=zeta*rzero
          zezer=dsqrt(zrzer**3)*dexp(-zrzer)
          g52zer=gammai(d5b2,zrzer,zrzer*zezer,iflagz)
          if(iflagb*iflagz.eq.1) then
            g52dif1=g52zer-g52bot
           else
            g52dif1=g52zer-g52bot+g5b2
          endif
          g32dif1=d2b3*(g52dif1+zezer-zebot)
          if(iflagt*iflagz.eq.1) then
            g52dif2=g52top-g52zer
           else
            g52dif2=g52top-g52zer+g5b2
          endif
          g32dif2=d2b3*(g52dif2+zetop-zezer)
          areao=(dabs(g32dif1+co*g52dif1/zeta) &
                +dabs(g32dif2+co*g52dif2/zeta))/ &
                (bot*dsqrt(zeta**3))
        endif

! Sample cos(theta)
        raver=half*(rmino(i)+rtry)
        deltt=thetamx(raver,znuc(iwctype(nearo)))
        costht=one-deltt*random_dp()
        zprime=rtry*costht
        sintht=dsqrt(one-costht*costht)

! For molecules deltt may not be the same for forward and reverse
! moves, so it is necessary to include this in areao and arean
        areao=areao*deltt

! Truncate phi variation if it goes through zero
! Sample phi by rejection. Note it is OK to have a) theta or sin(theta)
! and b) rtry/rold(i) or raver, as long as the forward and reverse probs.
! are consistent.
! If we do not limit term to be <=1 then use commented out lines.
        term=dmin1(voldp*raver*sintht,one)
!lim    term=voldp*raver*sintht
        fmax=one+term
      60   phitry=pi*random_dp()
          cosphi=dcos(phitry)
          top=one+term*cosphi
!lim      top=abs(one+term*cosphi)
          if(top.gt.random_dp()*fmax) goto 70
        goto 60
      70   fxop=fxop*top

!lim    if(term.gt.one) then
!lim      phizer=dacos(-one/term)
!lim      areao=areao*((two/pi)*(phizer+term*dsin(phizer))-one)
!lim    endif

! Calculate x and y coordinates in local coordinate system
        xprime=rtry*sintht*dcos(phitry)
        yprime=dsqrt(max(zero,rtry*rtry-xprime**2-zprime**2))
        if(random_dp().lt.half) yprime=-yprime

! Convert back to original coordinate system
        do ic=1,3
          rvminno(ic,i)=xaxis(ic)*xprime+yaxis(ic)*yprime+zaxis(ic) &
          *zprime
          xnew(ic,i)=rvminno(ic,i)+cent(ic,nearo)
        enddo
        rminno(i)=rtry

! Do geometrical rejections for molecules
        rminn(i)=99.d9
        do j=1,ncent
          dist=zero
          do ic=1,3
            dist=dist+(xnew(ic,i)-cent(ic,j))**2
          enddo
          if(dist.lt.rminn(i)) then
            rminn(i)=dist
            nearestn(i)=j
          endif
        enddo
        nearn=nearestn(i)
        rminn(i)=dsqrt(rminn(i))
        rminon(i)=zero
        dot=zero
        do  ic=1,3
          rvminn(ic,i)=xnew(ic,i)-cent(ic,nearestn(i))
          rvminon(ic,i)=xold(ic,i)-cent(ic,nearestn(i))
          rminon(i)=rminon(i)+rvminon(ic,i)**2
          dot=dot+rvminn(ic,i)*rvminon(ic,i)
        enddo
        rminon(i)=dsqrt(rminon(i))
        dot=dot/(rminn(i)*rminon(i))
        costht=dot
        sintht=dsqrt(one-costht*costht)
        ravern=half*(rminn(i)+rminon(i))
        delttn(i)=thetamx(ravern,znuc(iwctype(nearestn(i))))
        igeometrical=0
        if(rminon(i).gt.rminn(i)*deltar .or. dot.lt.one-delttn(i))then
          igeometrical=1
          if(ipr.gt.3) write(ounit,*) 'igeo',i,igeometrical
          p=zero
          q=one
          goto 208
        endif

! rratio^2 is needed for the density of the angular moves
        rratio=rminno(i)/rminon(i)

        if(ipr.ge.1) then
          rtest=dsqrt(rvminn(1,i)**2+rvminn(2,i)**2+rvminn(3,i)**2)
          rtest2=dsqrt(xprime**2+yprime**2+zprime**2)
          write(ounit,'(''ELECTRON'',i6)') i
          write(ounit,'(''rtest,rtest2,rtry'',9d14.6)')rtest,rtest2,rtry, &
          rtest-rtry,rtest2-rtry
          write(ounit,'(''psido='',9d12.4)') psido(1)
          write(ounit,'(''vold='',9d12.4)') (vold(ic,i),ic=1,3)
          write(ounit,'(''voldr,voldp='',9d12.4)') voldr,voldp
          write(ounit,'(''axes='',(3f8.4,3x))') xaxis,yaxis,zaxis
          write(ounit,'(''rmino(i),rmax1,rmax2,rzero'',9d12.4)') &
          rmino(i),rmax1,rmax2,rzero
          write(ounit,'(''rtry,costht,sintht,phitry'',9d12.4)') rtry,costht, &
          sintht,phitry
          write(ounit,'(''fxop'',9d12.4)') fxop
        endif

! calculate psi at new configuration
      iel=i

      call psie(iel,xnew,psidn,psijn,ipass,0)

      if(iguiding.eq.0) then
        psidg=psidn(1)
        psig=psidn(1)*exp(psijn(1))
      else
        call determinant_psig(psidn,psijn,psig)
      endif

      if(psig.eq.0.d0) then
        p=zero
        q=one
        goto 208
      endif

      call compute_determinante_grad(iel,psig,psidn,psijn,vnew(1,iel),0)

      if(ipr.gt.1) then
        write(ounit,'(''psidn,psig ='',2d12.4)') psidn(1),psig
      endif

      psi2n(1)=2*(dlog(dabs(psig)))

      if(node_cutoff.ne.0) then
        do jel=1,nup
          if(jel.ne.iel) call compute_determinante_grad(jel,psig,psidn,psijn,vnew(1,jel),iflag_up)
        enddo

        do jel=nup+1,nelec
          if(jel.ne.iel) call compute_determinante_grad(jel,psig,psidn,psijn,vnew(1,jel),iflag_dn)
        enddo

        call nodes_distance(vnew,distance_node,0)
        rnorm_nodes=rnorm_nodes_num(distance_node,eps_node_cutoff)/distance_node

        psi2n(1)=psi2n(1)+2*dlog(rnorm_nodes)

        if(ipr.gt.1) then
          write(ounit,'(''distance_node='',d12.4)') distance_node
          write(ounit,'(''rnorm_nodes='',d12.4)') rnorm_nodes
          write(ounit,'(''psid2n_ncut='',f9.4)') psi2n(1)
        endif
      endif

! calculate probability for reverse transition
      fxnp=one
! Choose lower and upper values of r sampling
        rbot=rminn(i)*deltri
        rtop=rminn(i)*deltar
! Calculate magnitude of the velocity in the radial direction
        vnewr=zero
        do ic=1,3
          vnewr=vnewr+vnew(ic,i)*rvminn(ic,i)
        enddo
        vnewr=vnewr/rminn(i)

! Place x-axis along direction of angular change and
! Calculate the velocity in the phi direction
        vnewp=zero
        do ic=1,3
          xaxis(ic)=vnew(ic,i)-vnewr*rvminn(ic,i)/rminn(i)
          vnewp=vnewp+xaxis(ic)**2
        enddo
        vnewp=dsqrt(vnewp)
        if(vnewp.lt.eps) then
          xaxis(1)=eps*(one-rvminn(1,i)**2)
          xaxis(2)=eps*(-rvminn(1,i)*rvminn(2,i))
          xaxis(3)=eps*(-rvminn(1,i)*rvminn(3,i))
          vnewp=eps*sqrt(one-rvminn(1,i)**2*(two-rminn(i)**2))
        endif
        do ic=1,3
          xaxis(ic)=xaxis(ic)/vnewp
        enddo

! Limit radial component of velocity.
! It may be a good idea to limit it if it is positive too.
        nearn=nearestn(i)
        if(nloc.eq.0) then
          zcusp=znuc(iwctype(nearn))
         else
          zcusp=zero
        endif
        vnewr=max(vnewr,-2*znuc(iwctype(nearn)))
!       vnewr=min(vnewr,2*znuc(iwctype(nearn)))

! Temporary test of fbias
!       vnewr=vnewr*fbias
        vnewp=vnewp*fbias

        root=(zcusp+vnewr)*(zcusp+vnewr-four/rminn(i))
        if(root.ge.zero) then
          root=sqrt(root)
          zeta=half*(zcusp-vnewr-root)
          if(zeta.le.zero) zeta=half*(zcusp-vnewr+root)
          if(zeta.le.zero) zeta=eps
         else
          if(vnewr.lt.zero) then
            zeta=-vnewr-eps
           else
            zeta=one
          endif
        endif
        co=(zeta+vnewr)/(one-(zeta+vnewr)*rminn(i))

        bot=sqrt(rminn(i))*abs(one+co*rminn(i))*dexp(-zeta*rminn(i))
        top=sqrt(rminon(i))*abs(one+co*rminon(i))*dexp(-zeta*rminon(i))
        fxnp=fxnp*top/bot

! Calculate the integral of T
        rzero=-one/co
        zrbot=zeta*rbot
        zrtop=zeta*rtop
        zebot=dsqrt(zrbot**3)*dexp(-zrbot)
        zetop=dsqrt(zrtop**3)*dexp(-zrtop)
        g52bot=gammai(d5b2,zrbot,zrbot*zebot,iflagb)
        g52top=gammai(d5b2,zrtop,zrtop*zetop,iflagt)
        if(rzero.lt.rbot .or. rzero.gt.rtop) then
          if(iflagb*iflagt.eq.1) then
            g52dif=g52top-g52bot
           else
            g52dif=g52top-g52bot+g5b2
          endif
          g32dif=d2b3*(g52dif+zetop-zebot)
          arean=dabs(g32dif+co*g52dif/zeta)/(bot*dsqrt(zeta**3))
         else
          zrzer=zeta*rzero
          zezer=dsqrt(zrzer**3)*dexp(-zrzer)
          g52zer=gammai(d5b2,zrzer,zrzer*zezer,iflagz)
          if(iflagb*iflagz.eq.1) then
            g52dif1=g52zer-g52bot
           else
            g52dif1=g52zer-g52bot+g5b2
          endif
          g32dif1=d2b3*(g52dif1+zezer-zebot)
          if(iflagt*iflagz.eq.1) then
            g52dif2=g52top-g52zer
           else
            g52dif2=g52top-g52zer+g5b2
          endif
          g32dif2=d2b3*(g52dif2+zetop-zezer)
          arean=(dabs(g32dif1+co*g52dif1/zeta) &
                +dabs(g32dif2+co*g52dif2/zeta))/ &
                (bot*dsqrt(zeta**3))
        endif

! For molecules deltt may not be the same for forward and reverse
! moves, so it is necessary to include this in areao and arean
        arean=arean*delttn(i)

! Truncate phi variation if it goes through zero
! Note it is OK to have a) theta or sin(theta)
! and b) rtry/rold(i) or raver, as long as the forward and reverse probs.
! are consistent.
! If we do not limit term to be <=1 then use commented out lines.
        term=dmin1(vnewp*ravern*sintht,one)
!lim    term=vnewp*raver*sintht
! Determine cos(phi)
        cosphi=zero
        rnorm=zero
        do ic=1,3
          term2=rvminon(ic,i)/rminon(i)-costht*rvminn(ic,i)/rminn(i)
          rnorm=rnorm+term2*term2
          cosphi=cosphi+term2*xaxis(ic)
        enddo
        cosphi=cosphi/dsqrt(rnorm)
        fxnp=fxnp*(one+term*cosphi)
!lim    fxnp=fxnp*abs(one+term*cosphi)

!lim    if(term.gt.one) then
!lim      phizer=dacos(-one/term)
!lim      arean=arean*((two/pi)*(phizer+term*dsin(phizer))-one)
!lim    endif

! p is the probability of accepting new move
      p=rratio**2*exp(psi2n(1)-psi2o(1,1))*dabs((fxnp*areao)/(fxop*arean))

        if(ipr.ge.1) then
          write(ounit,'(''rminn,rvminn,vnew,vnewr'',9d12.4)') &
          rminn(i),(rvminn(ic,i),ic=1,3),(vnew(ic,i),ic=1,3),vnewr
          write(ounit,'(''vnew='',9d12.4)') (vnew(ic,i),ic=1,3)
          write(ounit,'(''vnewr,vnewp='',9d12.4)') vnewr,vnewp
          write(ounit,'(''axes='',(3f8.4,3x))') xaxis,yaxis,zaxis
          write(ounit,'(''rminn(i),rmax1,rmax2,rzero'',9d12.4)') &
          rminn(i),rmax1,rmax2,rzero
          write(ounit,'(''rtry,costht,sintht,phitry,cosphi'',9d12.4)') rtry, &
          costht,sintht,phitry,cosphi
          write(ounit,'(''fxop,fxnp,areao,arean,psi2n,psi2o,p'',9d12.4)') &
                        fxop,fxnp,areao,arean,psi2n(1),psi2o(1,1),p
          if(dabs(vnew(1,i))+dabs(vnew(1,i))+dabs(vnew(1,i)).gt.10d+8) then
            do ii=1,i-1
              write(ounit,*) (xold(k,ii),k=1,3)
            enddo
            write(ounit,*) (xnew(k,i),k=1,3)
            do ii=i+1,nelec
              write(ounit,*) (xold(k,ii),k=1,3)
            enddo
          endif
        endif

      p=dmin1(one,p)
      q=one-p

      208 continue

! accept new move with probability p
! Note when one electron moves the velocity on all electrons change.
      if (random_dp().lt.p) then
!UNDO
#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
        !if (ioptjas.eq.0) then
        rc = qmckl_get_jastrow_champ_single_accept(qmckl_ctx(qmckl_no_ctx))
        !endif
#endif
        rmino(i)=rminn(i)
        nearesto(i)=nearestn(i)
        psi2o(1,1)=psi2n(1)
        do ic=1,3
          xold(ic,i)=xnew(ic,i)
          rvmino(ic,i)=rvminn(ic,i)
        enddo
        if(node_cutoff.gt.0) then
          do ic=1,3
            do ii=1,nelec
              vold(ic,ii)=vnew(ic,ii)
            enddo
          enddo
        endif
        do istate=1,nstates
          psido(istate)=psidn(istate)
          psijo(stoj(istate))=psijn(stoj(istate))
        enddo
        acc=acc+one
        call jassav(i,0)
        call detsav(i,0)
        if(ipr.ge.1) write(ounit,*)'METROP ACCEPT'
       else
        if(ipr.ge.1) write(ounit,*)'METROP REJECT'
        do ic=1,3
          xnew(ic,i)=xold(ic,i)
        enddo
        if(igeometrical.eq.0) call distancese_restore(i)
      endif

      call update_ymat(i)

      enddo

! loop over secondary configurations
      do ifr=2,nforce
        call strech(xold,xstrech,ajacob,ifr,1)
        call hpsi(xstrech,psido,psijo,ekino,eold(1,ifr),ipass,ifr)
        do istate=1,nstates
          j=stoj(istate)
          psi2o(istate,ifr)=2*(dlog(dabs(psido(istate)))+psijo(j))+dlog(ajacob)
        enddo
      enddo

      call check_orbitals_reset

! primary configuration
      if(nforce.gt.1) call strech(xold,xstrech,ajacob,1,0)
      call hpsi(xold,psido,psijo,ekino,eold(1,1),ipass,1)
      do istate=1,nstates
         j=stoj(istate)
         psi2o(istate,1)=2*(dlog(dabs(psido(istate)))+psijo(j))
      enddo

      if(iguiding.eq.0) then
        psidg=psido(1)
        psig=psido(1)*exp(psijo(1))
       else
        call determinant_psig(psido,psijo,psig)
      endif

      if(ipr.gt.1) then
        write(ounit,'(''psid,psig ='',2d12.4)') psido(1),psig
      endif


      rnorm_nodes=1.d0
      if(node_cutoff.gt.0) then
        do jel=1,nelec
          call compute_determinante_grad(jel,psig,psido(1),psijo,vold(1,jel),1)
        enddo
        call nodes_distance(vold,distance_node,1)
        rnorm_nodes=rnorm_nodes_num(distance_node,eps_node_cutoff)/distance_node
        psig=psig*rnorm_nodes
        if(ipr.gt.1) then
          write(ounit,'(''distance_node='',d12.4)') distance_node
          write(ounit,'(''rnorm_nodes='',d12.4)') rnorm_nodes
          write(ounit,'(''psig_ncut='',d12.4)') psidg
        endif
        distance_node_sum=distance_node_sum+distance_node
      endif

      do istate=1,nstates
        j=stoj(istate)
        wtg_sqrt(istate)=psido(istate)*exp(psijo(j))/psig
        wtg(istate)=wtg_sqrt(istate)*wtg_sqrt(istate)

! form expected values of e, pe, etc.
        esum1(istate)=eold(istate,1)
        wsum(istate,1)=wsum(istate,1)+wtg(istate)
        esum(istate,1)=esum(istate,1)+eold(istate,1)*wtg(istate)
        pesum(istate)=pesum(istate)+(eold(istate,1)-ekino(istate))*wtg(istate)
        tpbsum(istate)=tpbsum(istate)+ekino(istate)*wtg(istate)
      enddo

      if(ipr.gt.1) write(ounit,'(''energy reweighted '',d12.4)') eold(1,1)*wtg(1)

! normal component efield on cavity surface to compute a new set of polarization charges
      if(ichpol.eq.1) call qpcm_efield(nelec,xold)
! efield dovuto agli elettroni sui siti dei dipoli
      if(ich_mmpol.eq.1) call mmpol_efield(nelec,xold)

! use 'new' not 'old' value
      call pcm_sum(wtg(1),0.d0)
      call mmpol_sum(wtg(1),0.d0)
      call prop_sum(wtg(1),0.d0)
      call force_analy_sum(wtg(1),0.d0,eold(1,1),0.0d0)

      call optjas_sum(wtg,zero_array,eold(1,1),eold(1,1),0)
      call optorb_sum(wtg,zero_array,eold(1,1),eold(1,1),0)
      call optci_sum(wtg(1),0.d0,eold(1,1),eold(1,1))

      call optx_jas_orb_sum(wtg(1),zero_array,0)
      call optx_jas_ci_sum(wtg(1),0.d0,eold(1,1),eold(1,1))
      call optx_orb_ci_sum(wtg(1),0.d0)

      if(irun.eq.1) call optwf_store(ipass,wtg,wtg_sqrt,psido,eold(1,1))

      call efficiency_sample(ipass,psido,psijo,psig)

      call acues1(wtg)

      do istate=1,nstates
        esum1(istate)=eold(istate,1)
      enddo
      call acusig(wtg)

      do ifr=2,nforce
        do istate=1,nstates
          wstro=exp(psi2o(istate,ifr)-psi2o(istate,1))*wtg(istate)
          esum(istate,ifr)=esum(istate,ifr)+eold(istate,ifr)*wstro
          wsum(istate,ifr)=wsum(istate,ifr)+wstro
        enddo
      enddo

! rewrite psi2o for next metropolis step if you are sampling guiding
      if(iguiding.gt.0) psi2o(1,1)=2*(dlog(dabs(psig)))

      if(node_cutoff.gt.0) then
        psi2o(1,1)=psi2o(1,1)+2*dlog(rnorm_nodes)
      endif
      return
      end

!contains
      elemental pure function thetamx(r,z)
      use metropolis, only: deltat
      use precision_kinds, only: dp

      implicit none
      real(dp), intent(in) :: r, z
      real(dp)             :: thetamx
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0

      thetamx=deltat+(two-deltat)/(one+(z*r)**2)

      end function
end module
