      subroutine metrop6(ipass,irun)
c Written by Cyrus Umrigar
c Uses the accelerated Metropolis method described in:
c 1) Accelerated Metropolis Method, C.J. Umrigar, PRL 71, 408 (1993).
c 2) Variational Monte Carlo Basics and Applications to Atoms and Molecules,
c    C.J. Umrigar, in "Quantum Monte Carlo Methods in Physics and Chemistry",
c    edited by M.P. Nightingale and C.J. Umrigar. NATO ASI Series, Series C,
c    Mathematical and Physical Sciences, Vol. C-525,
c    (Kluwer Academic Publishers, Boston, 1999)
      use vmc_mod, only: MELEC, MORB, MDET
      use vmc_mod, only: nrad
      use vmc_mod, only: delri
      use atom, only: znuc, cent, iwctype, ncent
      use mstates_mod, only: MSTATES

      use const, only: pi, fbias, nelec, ipr
      use config, only: delttn, eold, nearestn, nearesto, peo, psi2n, psi2o
      use config, only: psido, psijo, rminn, rminno, rmino, rminon, rvminn, rvminno, rvmino, rvminon
      use config, only: tjfoo, vnew, vold, xnew, xold
      use csfs, only: nstates
      use elec, only: nup
      use estsum, only: acc, esum, esum1, pesum, r2sum, tjfsum, tpbsum
      use forcepar, only: nforce
      use forcewt, only: wsum
      use kinet, only: dtdx2n, dtdx2o
      use stats, only: rejmax
      use step, only: ekin, ekin2, rprob, suc, trunfb, try
      use tmpnode, only: distance_node_sum
      use const2, only: deltar, deltat
      use contr3, only: mode
      use pseudo, only: nloc
      use mmpol_cntrl, only: ich_mmpol
      use mstates_ctrl, only: iguiding
      use pcm_cntrl, only: ichpol
      use method_opt, only: method
      use multislatern, only: ddorbn, detn, dorbn, orbn
      use inputflags, only: node_cutoff, eps_node_cutoff

      implicit real*8(a-h,o-z)


      parameter (zero=0.d0,one=1.d0,two=2.d0,four=4.d0,half=0.5d0)
      parameter (d3b2=1.5d0,d5b2=2.5d0,d2b3=.666666666666667d0)
      parameter (eps=1.d-10)
c     parameter (g3b2=.886226925452758d0)
      parameter (g5b2=1.329340388179137d0)
c g3b2, g5b2 are gamma3/2), gamma(5/2)


c The moves are now being made in local r,theta phi coordinates.

c The foll. additions have been made:
c 1) Slater form of Tij.
c 2) Make theta_max a function of r
c 3) Generalize to molecules. This requires geometric rejections.

c The foll. still need to be tried:
c 1) Quadratic, gaussian, Morse and Exp(-zeta*r)+co*Exp(-r) forms of Tij
c    Last 2 are prob. best

 
c TMP
      dimension xstrech(3,MELEC)
      dimension xaxis(3),yaxis(3),zaxis(3),idist(MELEC)
      dimension ddx_ref(3)
      dimension psidn(MSTATES) ,wtg(MSTATES)
      
c     area(ri,r1,r2,v)=dabs((one/sqrt(ri))*
c    &(r2**d3b2*(two*(one-v*ri)/3+.4d0*v*r2)
c    &-r1**d3b2*(two*(one-v*ri)/3+.4d0*v*r1)))

      thetamx(r,z)=deltat+(two-deltat)/(one+(z*r)**2)

      mode='vmc_mov1    '

      deltri=one/deltar

      call check_orbitals
      
      do 300 i=1,nelec

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
          psig=psido(1)
         else
          call determinant_psig(psido,psig)
        endif
        call compute_determinante_grad(i,psig,psido,vold(1,i),1)

        fxop=one
        nearo=nearesto(i)
        if(nloc.eq.0) then
          zcusp=znuc(iwctype(nearo))
         else
          zcusp=zero
        endif
c Although rmino is saved, recalculate it, otherwise there may be a
c build up of errors via zaxis. I don't think so.
        rmino(i)=dsqrt(rvmino(1,i)**2+rvmino(2,i)**2+rvmino(3,i)**2)
c Choose lower and upper values of r sampling
        rbot=rmino(i)*deltri
        rtop=rmino(i)*deltar
c Calculate magnitude of the velocity in the radial direction
        voldr=zero
        do 10 ic=1,3
   10     voldr=voldr+vold(ic,i)*rvmino(ic,i)
        
        voldr=voldr/rmino(i)

c Place x-axis along direction of angular change and
c Calculate the velocity in the phi direction
        voldp=zero
        do 20 ic=1,3
          zaxis(ic)=rvmino(ic,i)/rmino(i)
          xaxis(ic)=vold(ic,i)-voldr*zaxis(ic)
   20     voldp=voldp+xaxis(ic)**2
        voldp=dsqrt(voldp)
        if(voldp.lt.eps) then
          xaxis(1)=eps*(one-zaxis(1)**2)
          xaxis(2)=eps*(-zaxis(1)*zaxis(2))
          xaxis(3)=eps*(-zaxis(1)*zaxis(3))
          voldp=eps*sqrt(one-zaxis(1)**2)
        endif
        do 30 ic=1,3
   30     xaxis(ic)=xaxis(ic)/voldp

c Limit radial component of velocity.
c It may be a good idea to limit it if it is positive too.
        voldr=max(voldr,-2*znuc(iwctype(nearo)))
c       voldr=min(voldr,2*znuc(iwctype(nearo)))

c y-axis is cross-product of z and x axes
        yaxis(1)=zaxis(2)*xaxis(3)-zaxis(3)*xaxis(2)
        yaxis(2)=zaxis(3)*xaxis(1)-zaxis(1)*xaxis(3)
        yaxis(3)=zaxis(1)*xaxis(2)-zaxis(2)*xaxis(1)

c Temporary test of fbias
c       voldr=voldr*fbias
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

c       write(6,'(''rmino(i),voldr,zeta,co='',9f10.5)')
c    &  rmino(i),voldr,zeta,co,(co-zeta-co*zeta*rmino(i))/
c    &  (one+co*rmino(i))

c Use Slater approx for radial fn
c Determine the maximum value of radial function for rejection sampling
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

c   Sample sqrt(r_f)*abs(1+co*r_f)*exp(-zeta*r_f) by rejection
        bot=sqrt(rmino(i))*abs(one+co*rmino(i))*dexp(-zeta*rmino(i))
   40   rtry=((deltar-deltri)*rannyu(0)+deltri)*rmino(i)
          top=sqrt(rtry)*abs(one+co*rtry)*dexp(-zeta*rtry)
          ratio=top/fmax
          rejmax=max(rejmax,ratio)
          if(ratio.gt.rannyu(0)) goto 50
        goto 40
   50   fxop=fxop*top/bot

c   Calculate the integral of T
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
          areao=(dabs(g32dif1+co*g52dif1/zeta)
     &          +dabs(g32dif2+co*g52dif2/zeta))/
     &          (bot*dsqrt(zeta**3))
        endif

c Sample cos(theta)
        raver=half*(rmino(i)+rtry)
        deltt=thetamx(raver,znuc(iwctype(nearo)))
        costht=one-deltt*rannyu(0)
        zprime=rtry*costht
        sintht=dsqrt(one-costht*costht)

c For molecules deltt may not be the same for forward and reverse
c moves, so it is necessary to include this in areao and arean
        areao=areao*deltt

c Truncate phi variation if it goes through zero
c Sample phi by rejection. Note it is OK to have a) theta or sin(theta)
c and b) rtry/rold(i) or raver, as long as the forward and reverse probs.
c are consistent.
c If we do not limit term to be <=1 then use commented out lines.
        term=dmin1(voldp*raver*sintht,one)
clim    term=voldp*raver*sintht
        fmax=one+term
   60   phitry=pi*rannyu(0)
          cosphi=dcos(phitry)
          top=one+term*cosphi
clim      top=abs(one+term*cosphi)
          if(top.gt.rannyu(0)*fmax) goto 70
        goto 60
   70   fxop=fxop*top

clim    if(term.gt.one) then
clim      phizer=dacos(-one/term)
clim      areao=areao*((two/pi)*(phizer+term*dsin(phizer))-one)
clim    endif

c Calculate x and y coordinates in local coordinate system
        xprime=rtry*sintht*dcos(phitry)
        yprime=dsqrt(max(zero,rtry*rtry-xprime**2-zprime**2))
        if(rannyu(0).lt.half) yprime=-yprime

c Convert back to original coordinate system
        do 80 ic=1,3
          rvminno(ic,i)=xaxis(ic)*xprime+yaxis(ic)*yprime+zaxis(ic)
     &    *zprime
   80     xnew(ic,i)=rvminno(ic,i)+cent(ic,nearo)
        rminno(i)=rtry

c Do geometrical rejections for molecules
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
        costht=dot
        sintht=dsqrt(one-costht*costht)
        ravern=half*(rminn(i)+rminon(i))
        delttn(i)=thetamx(ravern,znuc(iwctype(nearestn(i))))
        igeometrical=0
        if(rminon(i).gt.rminn(i)*deltar .or. dot.lt.one-delttn(i))then
          igeometrical=1
          if(ipr.gt.3) write(6,*) 'igeo',i,igeometrical
          p=zero
          q=one
          goto 208
        endif

c rratio^2 is needed for the density of the angular moves
        rratio=rminno(i)/rminon(i)

        if(ipr.ge.1) then
          rtest=dsqrt(rvminn(1,i)**2+rvminn(2,i)**2+rvminn(3,i)**2)
          rtest2=dsqrt(xprime**2+yprime**2+zprime**2)
          write(6,'(''ELECTRON'',i6)') i
          write(6,'(''rtest,rtest2,rtry'',9d14.6)')rtest,rtest2,rtry,
     &    rtest-rtry,rtest2-rtry
          write(6,'(''psido='',9d12.4)') psido(1)
          write(6,'(''vold='',9d12.4)') (vold(ic,i),ic=1,3)
          write(6,'(''voldr,voldp='',9d12.4)') voldr,voldp
          write(6,'(''axes='',(3f8.4,3x))') xaxis,yaxis,zaxis
          write(6,'(''rmino(i),rmax1,rmax2,rzero'',9f9.4)')
     &    rmino(i),rmax1,rmax2,rzero
          write(6,'(''rtry,costht,sintht,phitry'',9f9.4)') rtry,costht,
     &    sintht,phitry
          write(6,'(''fxop'',9f12.4)') fxop
        endif

c calculate psi at new configuration
      iel=i
      call psie(iel,xnew,psidn,psijn,ipass,0)
      if(iguiding.eq.0) then
        psig=psidn(1)
       else
        call determinant_psig(psidn,psig)
      endif

      if(psig.eq.0.d0) then
        p=zero
        q=one
        goto 208
      endif

      call compute_determinante_grad(iel,psig,psidn,vnew(1,iel),0)

      if(ipr.gt.1) then
        write(6,'(''psidn,psig ='',2d12.4)') psidn(1),psig
      endif

      psi2n(1)=2*(dlog(dabs(psig))+psijn)

      if(node_cutoff.ne.0) then
        do 100 jel=1,nup
  100     if(jel.ne.iel) call compute_determinante_grad(jel,psig,psidn,vnew(1,jel),iflag_up)

        do 105 jel=nup+1,nelec
  105     if(jel.ne.iel) call compute_determinante_grad(jel,psig,psidn,vnew(1,jel),iflag_dn)

        call nodes_distance(vnew,distance_node,0)
        rnorm_nodes=rnorm_nodes_num(distance_node,eps_node_cutoff)/distance_node

        psi2n(1)=psi2n(1)+2*dlog(rnorm_nodes)

        if(ipr.gt.1) then
          write(6,'(''distance_node='',d12.4)') distance_node
          write(6,'(''rnorm_nodes='',d12.4)') rnorm_nodes
          write(6,'(''psid2n_ncut='',f9.4)') psi2n(1)
        endif
      endif

c calculate probability for reverse transition
      fxnp=one
c Choose lower and upper values of r sampling
        rbot=rminn(i)*deltri
        rtop=rminn(i)*deltar
c Calculate magnitude of the velocity in the radial direction
        vnewr=zero
        do 110 ic=1,3
  110     vnewr=vnewr+vnew(ic,i)*rvminn(ic,i)
        vnewr=vnewr/rminn(i)

c Place x-axis along direction of angular change and
c Calculate the velocity in the phi direction
        vnewp=zero
        do 120 ic=1,3
          xaxis(ic)=vnew(ic,i)-vnewr*rvminn(ic,i)/rminn(i)
  120     vnewp=vnewp+xaxis(ic)**2
        vnewp=dsqrt(vnewp)
        if(vnewp.lt.eps) then
          xaxis(1)=eps*(one-rvminn(1,i)**2)
          xaxis(2)=eps*(-rvminn(1,i)*rvminn(2,i))
          xaxis(3)=eps*(-rvminn(1,i)*rvminn(3,i))
          vnewp=eps*sqrt(one-rvminn(1,i)**2*(two-rminn(i)**2))
        endif
        do 130 ic=1,3
  130     xaxis(ic)=xaxis(ic)/vnewp

c Limit radial component of velocity.
c It may be a good idea to limit it if it is positive too.
        nearn=nearestn(i)
        if(nloc.eq.0) then
          zcusp=znuc(iwctype(nearn))
         else
          zcusp=zero
        endif
        vnewr=max(vnewr,-2*znuc(iwctype(nearn)))
c       vnewr=min(vnewr,2*znuc(iwctype(nearn)))

c Temporary test of fbias
c       vnewr=vnewr*fbias
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

c Calculate the integral of T
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
          arean=(dabs(g32dif1+co*g52dif1/zeta)
     &          +dabs(g32dif2+co*g52dif2/zeta))/
     &          (bot*dsqrt(zeta**3))
        endif

c For molecules deltt may not be the same for forward and reverse
c moves, so it is necessary to include this in areao and arean
        arean=arean*delttn(i)

c Truncate phi variation if it goes through zero
c Note it is OK to have a) theta or sin(theta)
c and b) rtry/rold(i) or raver, as long as the forward and reverse probs.
c are consistent.
c If we do not limit term to be <=1 then use commented out lines.
        term=dmin1(vnewp*ravern*sintht,one)
clim    term=vnewp*raver*sintht
c Determine cos(phi)
        cosphi=zero
        rnorm=zero
        do 160 ic=1,3
          term2=rvminon(ic,i)/rminon(i)-costht*rvminn(ic,i)/rminn(i)
          rnorm=rnorm+term2*term2
  160     cosphi=cosphi+term2*xaxis(ic)
        cosphi=cosphi/dsqrt(rnorm)
        fxnp=fxnp*(one+term*cosphi)
clim    fxnp=fxnp*abs(one+term*cosphi)

clim    if(term.gt.one) then
clim      phizer=dacos(-one/term)
clim      arean=arean*((two/pi)*(phizer+term*dsin(phizer))-one)
clim    endif

c p is the probability of accepting new move
      p=rratio**2*exp(psi2n(1)-psi2o(1,1))*dabs((fxnp*areao)/(fxop*arean))

        if(ipr.ge.1) then
          write(6,'(''rminn,rvminn,vnew,vnewr'',9f10.4)')
     &    rminn(i),(rvminn(ic,i),ic=1,3),(vnew(ic,i),ic=1,3),vnewr
          write(6,'(''vnew='',9d12.4)') (vnew(ic,i),ic=1,3)
          write(6,'(''vnewr,vnewp='',9d12.4)') vnewr,vnewp
          write(6,'(''axes='',(3f8.4,3x))') xaxis,yaxis,zaxis
          write(6,'(''rminn(i),rmax1,rmax2,rzero'',9f9.4)')
     &    rminn(i),rmax1,rmax2,rzero
          write(6,'(''rtry,costht,sintht,phitry,cosphi'',9f9.4)') rtry,
     &    costht,sintht,phitry,cosphi
          write(6,'(''fxop,fxnp,areao,arean,psi2n,psi2o,p'',9f9.4)')
     &                  fxop,fxnp,areao,arean,psi2n(1),psi2o(1,1),p
          if(dabs(vnew(1,i))+dabs(vnew(1,i))+dabs(vnew(1,i)).gt.10d+8) then
            do ii=1,i-1
              write(6,*) (xold(k,ii),k=1,3)
            enddo
            write(6,*) (xnew(k,i),k=1,3)
            do ii=i+1,nelec
              write(6,*) (xold(k,ii),k=1,3)
            enddo
          endif
        endif

      p=dmin1(one,p)
      q=one-p

  208 continue
c Calculate as a function of the distance to the nucleus
c 1) acceptance,  2) force-bias truncation probability,
c 3) kinetic energy and it's fluctuation
c The K.E. is not quite correct, since we should use p times new
c and q times old, and keep track of which bin the old was in
      rold=dsqrt(xold(1,i)**2+xold(2,i)**2+xold(3,i)**2)
      rnew=dsqrt(xnew(1,i)**2+xnew(2,i)**2+xnew(3,i)**2)
      itryo=min(int(delri*rold)+1,nrad)
      itryn=min(int(delri*rnew)+1,nrad)
      try(itryo)=try(itryo)+1
      suc(itryo)=suc(itryo)+p
      if(try(itryo).lt.0.) write(6,'(''itryo,try'',i5,d13.5)')itryo,
     &try(itryo)
      if(suc(itryo).lt.0.) write(6,'(''itryo,suc'',i5,d13.5)')itryo,
     &suc(itryo)
      if(voldp*raver*sintht.gt.one) trunfb(itryo)=trunfb(itryo)+1
      rprob(itryo)=rprob(itryo)+q
      rprob(itryn)=rprob(itryn)+p
      do 210 ic=1,3
  210   r2sum=r2sum+p*xnew(ic,i)**2+q*xold(ic,i)**2

c accept new move with probability p
c Note when one electron moves the velocity on all electrons change.
      if (rannyu(0).lt.p) then
        idist(i)=itryn
        rmino(i)=rminn(i)
        nearesto(i)=nearestn(i)
        psi2o(1,1)=psi2n(1)
        do 240 ic=1,3
          xold(ic,i)=xnew(ic,i)
  240     rvmino(ic,i)=rvminn(ic,i)
        if(node_cutoff.gt.0) then
          do 245 ic=1,3
            do 245 ii=1,nelec
  245         vold(ic,ii)=vnew(ic,ii)
        endif
        do istate=1,nstates
          psido(istate)=psidn(istate)
        enddo
        psijo=psijn
        acc=acc+one
        call jassav(i,0)
        call detsav(i,0)
        if(ipr.ge.1) write(6,*)'METROP ACCEPT'
       else
        if(ipr.ge.1) write(6,*)'METROP REJECT'
        idist(i)=itryo
        do 250 ic=1,3
  250     xnew(ic,i)=xold(ic,i)
        if(igeometrical.eq.0) call distancese_restore(i)
      endif

      call update_ymat(i)

  300 continue

c loop over secondary configurations
      do 350 ifr=2,nforce
        call strech(xold,xstrech,ajacob,ifr,1)
        call hpsi(xstrech,psido(1),psijo,eold(1,ifr),ipass,ifr)
        do 350 istate=1,nstates
  350     psi2o(istate,ifr)=2*(dlog(dabs(psido(istate)))+psijo)+dlog(ajacob)

      call check_orbitals_reset
c primary configuration
      if(nforce.gt.1) call strech(xold,xstrech,ajacob,1,0)
      call hpsi(xold,psido(1),psijo,eold(1,1),ipass,1)
      do 355 istate=1,nstates
  355    psi2o(istate,1)=2*(dlog(dabs(psido(istate)))+psijo)

      if(iguiding.eq.0) then
        psidg=psido(1)
       else
        call determinant_psig(psido,psidg)
      endif

      if(ipr.gt.1) then
        write(6,'(''psid,psig ='',2d12.4)') psido(1),psidg
      endif

      rnorm_nodes=1.d0
      if(node_cutoff.gt.0) then
        do 356 jel=1,nelec
  356     call compute_determinante_grad(jel,psidg,psido(1),vold(1,jel),1)
        call nodes_distance(vold,distance_node,1)
        rnorm_nodes=rnorm_nodes_num(distance_node,eps_node_cutoff)/distance_node
        psidg=psido(1)*rnorm_nodes
        if(ipr.gt.1) then
          write(6,'(''distance_node='',d12.4)') distance_node
          write(6,'(''rnorm_nodes='',d12.4)') rnorm_nodes
          write(6,'(''psig_ncut='',d12.4)') psidg
        endif
        distance_node_sum=distance_node_sum+distance_node
      endif

      do 360 istate=1,nstates
        wtg(istate)=psido(istate)/psidg
        wtg(istate)=wtg(istate)*wtg(istate)
     
c form expected values of e, pe, etc.
        esum1(istate)=eold(istate,1)
        wsum(istate,1)=wsum(istate,1)+wtg(istate)
        esum(istate,1)=esum(istate,1)+eold(istate,1)*wtg(istate)
        pesum(istate)=pesum(istate)+peo(istate)*wtg(istate)
        tpbsum(istate)=tpbsum(istate)+(eold(istate,1)-peo(istate))*wtg(istate)
  360   tjfsum(istate)=tjfsum(istate)+tjfoo*wtg(istate)

      if(ipr.gt.1) write(6,'(''energy reweighted '',d12.4)') eold(1,1)*wtg(1)

c normal component efield on cavity surface to compute a new set of polarization charges
      if(ichpol.eq.1) call qpcm_efield(nelec,xold)
c efield dovuto agli elettroni sui siti dei dipoli
      if(ich_mmpol.eq.1) call mmpol_efield(nelec,xold)

c use 'new' not 'old' value
      call pcm_sum(wtg,0.d0)
      call mmpol_sum(wtg,0.d0)
      call prop_sum(wtg,0.d0)
      call force_analy_sum(wtg,0.d0,eold(1,1),0)

      call optjas_sum(wtg(1),0.d0,eold(1,1),eold(1,1),0)
      call optorb_sum(wtg(1),0.d0,eold(1,1),eold(1,1),0)
      call optci_sum(wtg(1),0.d0,eold(1,1),eold(1,1))

      call optx_jas_orb_sum(wtg,0.d0,0)
      call optx_jas_ci_sum(wtg,0.d0,eold(1,1),eold(1,1))
      call optx_orb_ci_sum(wtg,0.d0)

      if(irun.eq.1) call optwf_store(ipass,wtg,psido,eold(1,1))

      call efficiency_sample(ipass,psido,psidg)

      call acues1(wtg)

      do 370 istate=1,nstates
  370   esum1(istate)=eold(istate,1)
      call acusig(wtg)

      do 380 ifr=2,nforce
        do 380 istate=1,nstates
          wstro=exp(psi2o(istate,ifr)-psi2o(istate,1))*wtg(istate)
          esum(istate,ifr)=esum(istate,ifr)+eold(istate,ifr)*wstro
  380     wsum(istate,ifr)=wsum(istate,ifr)+wstro
      do 390 i=1,nelec
        dtdx2o(i)=dtdx2n(i)
        ekin(idist(i))=ekin(idist(i))+dtdx2o(i)*wtg(1)
  390   ekin2(idist(i))=ekin2(idist(i))+dtdx2o(i)**2*wtg(1)

c rewrite psi2o for next metropolis step if you are sampling guiding 
      if(iguiding.gt.0) psi2o(1,1)=2*(dlog(dabs(psidg))+psijo)

      if(node_cutoff.gt.0) then
        psi2o(1,1)=psi2o(1,1)+2*dlog(rnorm_nodes)
      endif

      return
      end
