module metrop_slat
contains
      subroutine metrop6_moveall(ipass,irun)
! Written by Cyrus Umrigar
! Uses the accelerated Metropolis method described in:
! 1) Accelerated Metropolis Method, C.J. Umrigar, PRL 71, 408 (1993).
! 2) Variational Monte Carlo Basics and Applications to Atoms and Molecules,
!    C.J. Umrigar, in "Quantum Monte Carlo Methods in Physics and Chemistry",
!    edited by M.P. Nightingale and C.J. Umrigar. NATO ASI Series, Series C,
!    Mathematical and Physical Sciences, Vol. C-525,
!    (Kluwer Academic Publishers, Boston, 1999)

      use acuest_mod, only: acues1,acusig
      use config,  only: delttn,eold,nearestn,nearesto,psi2n,psi2o
      use config,  only: psido,psijo,rminn,rminno,rmino,rminon,rvminn
      use config,  only: rvminno,rvmino,rvminon,vnew,vold,xnew,xold
      use constants, only: pi
      use contrl_file, only: ounit
      use control, only: ipr
      use determinante_mod, only: compute_determinante_grad
      use estsum,  only: acc,esum,esum1,pesum,tpbsum
      use error,   only: fatal_error
      use forcewt, only: wsum
      use gammai_mod, only: gammai
      use hpsi_mod, only: hpsi
      use hpsie,   only: psie
      use metropolis, only: deltar,deltat,fbias
      use metrop_mov1_slat, only: thetamx
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
      use pseudo,  only: nloc
      use random_mod, only: random_dp
      use stats,   only: rejmax
      use strech_mod, only: strech
      use system,  only: cent,iwctype,ncent,nelec,znuc

      implicit none

      integer :: i, iab, ic, iflag_dn
      integer :: iflag_up, iflagb, iflagt, iflagz
      integer :: ifr, igeometrical, ii, ipass
      integer :: irun, istate, j, k, nearn, nearo
      real(dp) :: ajacob, arean, areao, bot
      real(dp) :: co, cosphi, deltri, deltt
      real(dp) :: dist, dmin1, dot
      real(dp) :: fmax, fmax2, fxnp, fxop
      real(dp) :: g32dif, g32dif1, g32dif2, g52bot
      real(dp) :: g52dif, g52dif1, g52dif2, g52top
      real(dp) :: g52zer, p, phitry, phizer
      real(dp) :: q, r, ratio
      real(dp) :: rbot, rmax1, rmax2 = 0d0, rnew
      real(dp) :: rnorm, rold, root
      real(dp) :: rratio, rtest, rtest2, rtop
      real(dp) :: rtry, rzero, term
      real(dp) :: term2, top, vnewp
      real(dp) :: vnewr, voldp, voldr, wstro, wstrn
      real(dp) :: xprime, yprime, z, zcusp
      real(dp) :: zebot, zeta, zetop, zezer
      real(dp) :: zprime, zrbot, zrtop, zrzer
      real(dp), dimension(nelec) :: costht, sintht
      real(dp), dimension(nelec) :: raver, ravern
      real(dp), dimension(3,nelec) :: xstrech
      real(dp), dimension(3) :: xaxis, yaxis, zaxis
      real(dp), dimension(1) :: ekinn, ekino
      real(dp), dimension(1) :: psidn, psijn
      real(dp), dimension(1) :: p_dum, q_dum
      real(dp), dimension(MFORCE)  :: enew
      real(dp), dimension(1) :: zero_array = 0.0_dp
      real(dp), dimension(1) :: one_array = 1.0_dp
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

      if(iguiding.gt.0) call fatal_error ('Guiding wave function only in mov1')

      deltri=one/deltar
      fxop=one
      areao=one
      rratio=one

      call check_orbitals

! Start loop over electrons
      do i=1,nelec
        !call compute_determinante_grad(i,psido(1),psido,psijo,vold(1,i),1)

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

        if(ipr.ge.1) then
          write(ounit,'(''voldr,voldp'',f7.1,f7.2)') voldr,voldp
          write(ounit,'(''rmino(i),voldr,zeta,co='',9f10.5)') &
          rmino(i),voldr,zeta,co,(co-zeta-co*zeta*rmino(i))/ &
          (one+co*rmino(i))
        endif

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
          areao=areao*dabs(g32dif+co*g52dif/zeta)/(bot*dsqrt(zeta**3))
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
          areao=areao* &
          (dabs(g32dif1+co*g52dif1/zeta) &
          +dabs(g32dif2+co*g52dif2/zeta))/ &
          (bot*dsqrt(zeta**3))
        endif

! Sample cos(theta)
        raver(i)=half*(rmino(i)+rtry)
        deltt=thetamx(raver(i),znuc(iwctype(nearo)))
        costht(i)=one-deltt*random_dp()
        zprime=rtry*costht(i)
        sintht(i)=dsqrt(one-costht(i)*costht(i))

! For molecules deltt may not be the same for forward and reverse
! moves, so it is necessary to include this in areao and arean
        areao=areao*deltt

! Truncate phi variation if it goes through zero
! Sample phi by rejection. Note it is OK to have a) theta or sin(theta)
! and b) rtry or raver(i), as long as the forward and reverse probs.
! are consistent.
! If we do not limit term to be <=1 then use commented out lines.
        term=dmin1(voldp*raver(i)*sintht(i),one)
!lim    term=voldp*raver(i)*sintht(i)
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
        xprime=rtry*sintht(i)*cosphi
        yprime=dsqrt(max(zero,rtry*rtry-xprime**2-zprime**2))
        if(random_dp().lt.half) yprime=-yprime

! Convert back to original coordinate system
        do ic=1,3
          rvminno(ic,i)=xaxis(ic)*xprime+yaxis(ic)*yprime+zaxis(ic)*zprime
          xnew(ic,i)=rvminno(ic,i)+cent(ic,nearo)
        enddo
        rminno(i)=rtry

! Do geometrical rejections for molecules
        rminn(i)=99.d9
        do 85 j=1,ncent
          dist=zero
          do 84 ic=1,3
   84       dist=dist+(xnew(ic,i)-cent(ic,j))**2
          if(dist.lt.rminn(i)) then
            rminn(i)=dist
            nearestn(i)=j
          endif
   85     continue
        nearn=nearestn(i)
        rminn(i)=dsqrt(rminn(i))

        rminon(i)=zero
        dot=zero
        do 86  ic=1,3
          rvminn(ic,i)=xnew(ic,i)-cent(ic,nearestn(i))
          rvminon(ic,i)=xold(ic,i)-cent(ic,nearestn(i))
          rminon(i)=rminon(i)+rvminon(ic,i)**2
   86     dot=dot+rvminn(ic,i)*rvminon(ic,i)
        rminon(i)=dsqrt(rminon(i))
        dot=dot/(rminn(i)*rminon(i))
        costht(i)=dot
        sintht(i)=dsqrt(one-costht(i)*costht(i))
        ravern(i)=half*(rminn(i)+rminon(i))
        delttn(i)=thetamx(ravern(i),znuc(iwctype(nearestn(i))))
        if(rminon(i).gt.rminn(i)*deltar .or. dot.lt.one-delttn(i))then
          p=zero
          q=one
! Set psi2n(1) to any non-zero number to prevent wstrn from being infinite
! on the 1st step of a restart.  Value of psi2n(1) is irrelevant since p=0.
          psi2n(1)=one
          goto 208
        endif

! rratio^2 is needed for the density of the angular moves
        rratio=rratio*rminno(i)/rminon(i)

        if(ipr.ge.1) then
          rtest=dsqrt(rvminn(1,i)**2+rvminn(2,i)**2+rvminn(3,i)**2)
          rtest2=dsqrt(xprime**2+yprime**2+zprime**2)
          write(ounit,'(''rtest,rtest2,rtry'',9d14.6)')rtest,rtest2,rtry, &
          rtest-rtry,rtest2-rtry
          write(ounit,'(''vold='',9d12.4)') (vold(ic,i),ic=1,3)
          write(ounit,'(''voldr,voldp='',9d12.4)') voldr,voldp
          write(ounit,'(''axes='',(3f8.4,3x))') xaxis,yaxis,zaxis
          write(ounit,'(''rmino(i),rmax1,rmax2,rzero'',9d12.4)')rmino(i), &
          rmax1,rmax2,rzero
          write(ounit,'(''rtry,costht(i),sintht(i),phitry'',9f9.4)') rtry, &
          costht(i),sintht(i),phitry
          write(ounit,'(''fxop'',9d12.4)') fxop,areao
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
      fxnp=one
      arean=one
      do i=1,nelec
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
        if(ipr.ge.1) then
          write(ounit,'(''rminn(i),vnewr,zeta,co='',9f10.5)') &
          rminn(i),vnewr,zeta,co,(co-zeta-co*zeta*rminn(i))/ &
          (one+co*rminn(i))
        endif

        bot=sqrt(rminn(i))*abs(one+co*rminn(i))*dexp(-zeta*rminn(i))
        top=sqrt(rminon(i))*abs(one+co*rminon(i))*dexp(-zeta*rminon(i))
        fxnp=fxnp*top/bot

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
          arean=arean*dabs(g32dif+co*g52dif/zeta)/(bot*dsqrt(zeta**3))
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
          arean=arean* &
          (dabs(g32dif1+co*g52dif1/zeta) &
          +dabs(g32dif2+co*g52dif2/zeta))/ &
          (bot*dsqrt(zeta**3))
        endif


! For molecules deltt may not be the same for forward and reverse
! moves, so it is necessary to include this in areao and arean
        arean=arean*delttn(i)

! Truncate phi variation if it goes through zero
! Note it is OK to have a) theta or sin(theta)
! and b) rminon(i) or ravern(i), as long as the forward and reverse
! probs are consistent.
! If we do not limit term to be <=1 then use commented out lines.
        term=dmin1(vnewp*ravern(i)*sintht(i),one)
!lim    term=vnewp*ravern(i)*sintht(i)
! Determine cos(phi)
        cosphi=zero
        rnorm=zero
        do ic=1,3
          term2=rvminon(ic,i)/rminon(i)-costht(i)*rvminn(ic,i)/rminn(i)
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

        if(ipr.ge.1) then
          write(ounit,'(''vnewr,vnewp'',f7.1,f7.2)') vnewr,vnewp
          write(ounit,'(''rvminn'',9f10.4)') (rvminn(ic,i),ic=1,3)
          write(ounit,'(''vnew='',9d12.4)') (vnew(ic,i),ic=1,3)
          write(ounit,'(''vnewr,vnewp='',9d12.4)') vnewr,vnewp
          write(ounit,'(''axes='',(3f8.4,3x))') xaxis,yaxis,zaxis
          write(ounit,'(''costht,sintht,phitry,cosphi'',9f9.4)') &
          costht(i),sintht(i),phitry,cosphi
          write(ounit,'(9d12.4)') fxnp
          write(ounit,'(''fxnp'',9d12.4)') fxnp,arean
          write(ounit,*)
        endif
      enddo

! p is the probability of accepting new move
      p=rratio**2*exp(psi2n(1)-psi2o(1,1))*dabs((fxnp*areao)/(fxop*arean))

      if(ipr.ge.1) then
        write(ounit,'(''fxop,fxnp,areao,arean,psi2n,psi2o,p'',9d12.4)') &
                    fxop,fxnp,areao,arean,psi2n(1),psi2o(1,1),p
        write(ounit,*)
      endif

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
          rmino(i)=rminn(i)
          nearesto(i)=nearestn(i)
          do ic=1,3
            xold(ic,i)=xnew(ic,i)
            rvmino(ic,i)=rvminn(ic,i)
            vold(ic,i)=vnew(ic,i)
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
end module metrop_slat
