      subroutine dmc_ps
c Written by Cyrus Umrigar and Claudia Filippi
c Uses the diffusion Monte Carlo algorithm described in:
c 1) A Diffusion Monte Carlo Algorithm with Very Small Time-Step Errors,
c    C.J. Umrigar, M.P. Nightingale and K.J. Runge, J. Chem. Phys., 99, 2865 (1993)
c modified to do accept/reject after single-electron moves and to
c remove portions related to nuclear cusps.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Control variables are:
c idmc         < 0     VMC
c              > 0     DMC
c abs(idmc)    = 1 **  simple kernel using dmc.brock.f
c              = 2     good kernel using dmc_good or dmc_good_inhom
c ipq         <= 0 *   do not use expected averages
c             >= 1     use expected averages (mostly in all-electron move algorithm)
c itau_eff    <=-1 *   always use tau in branching (not implemented)
c              = 0     use 0 / tau for acc /nacc moves in branching
c             >= 1     use tau_eff (calcul in equilibration runs) for all moves
c iacc_rej    <=-1 **  accept all moves (except possibly node crossings)
c              = 0 **  use weights rather than accept/reject
c             >= 1     use accept/reject
c icross      <=-1 **  kill walkers that cross nodes (not implemented)
c              = 0     reject walkers that cross nodes
c             >= 1     allow walkers to cross nodes
c                      (OK since crossing prob. goes as tau^(3/2))
c icuspg      <= 0     approximate cusp in Green function
c             >= 1     impose correct cusp condition on Green function
c icut_br     <= 0     do not limit branching
c             >= 1 *   use smooth formulae to limit branching to (1/2,2)
c                      (bad because it makes energies depend on E_trial)
c icut_e      <= 0     do not limit energy
c             >= 1 *   use smooth formulae to limit energy (not implemented)

c *  => bad option, modest deterioration in efficiency or time-step error
c ** => very bad option, big deterioration in efficiency or time-step error
c So, idmc=6,66 correspond to the foll. two:
c 2 1 1 1 0 0 0 0 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
c 2 1 0 1 1 0 0 0 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
c Another reasonable choice is:
c 2 1 0 1 1 1 1 0 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      use vmc_mod, only: nrad
      use vmc_mod, only: delri
      use const, only: etrial, hb, ipr, nelec
      use forcepar, only: istrech, nforce
      use age, only: iage, ioldest, ioldestmx
      use contrldmc, only: iacc_rej, icross, icut_br, icut_e, idmc, ipq, nfprod, rttau, tau
      use atom, only: cent
      use estcum, only: ipass
      use config, only: d2o, peo_dmc, psido_dmc, psijo_dmc, vold_dmc, xold_dmc
      use stats, only: acc, dfus2ac, dfus2un, dr2ac, dr2un, nacc, nodecr, trymove
      use estsum, only: efsum1, egsum1, esum1_dmc
      use estsum, only: pesum_dmc, r2sum, risum, tausum, tjfsum_dmc, tpbsum_dmc
      use estsum, only: wfsum1, wgsum1, wsum1
      use force_dmc, only: itausec, nwprod
      use derivest, only: derivsum
      use step, only: rprob
      use branch, only: eest, eigv, eold, ff, fprod, nwalk, pwt, wdsumo, wgdsumo, wt
      use branch, only: wthist
      use casula, only: i_vpsp, icasula
      use jacobsave, only: ajacob, ajacold
      use elec, only: nup
      use velratio, only: fratio, xdrifted
      use control_dmc, only: dmc_irstar, dmc_nconf
      use inputflags, only: node_cutoff, eps_node_cutoff, icircular, idrifdifgfunc
      use precision_kinds, only: dp
      use contrl_file,    only: ounit
      implicit none

      interface
         function rannyu(idum)
          use precision_kinds, only: dp
         implicit none
         integer,intent(in) :: idum
         real(dp) :: rannyu
         end function rannyu
      end interface

      interface
         function gauss()
          use precision_kinds, only: dp
         implicit none
         real(dp) :: gauss
         end function gauss
      end interface

      integer :: i, iaccept, iel
      integer :: iflag_dn, iflag_up, ifr, ii
      integer :: imove, ipmod, ipmod2, iw
      integer :: iwmod, j, jel, k
      integer :: ncall, ncount_casula, nmove_casula
      integer, dimension(nelec) :: itryo
      integer, dimension(nelec) :: itryn
      integer, dimension(nelec) :: iacc_elec
      real(dp) :: d2n, den, deo, dfus
      real(dp) :: dfus2n, dfus2o, distance_node, distance_node_ratio2
      real(dp) :: dmin1, dr2, drifdif, drifdifgfunc
      real(dp) :: drifdifr, drifdifs, drift, dwt
      real(dp) :: dx, e_cutoff, enew
      real(dp) :: ewtn, ewto, expon, ffi
      real(dp) :: ffn, fration, ginv
      real(dp) :: p, pen, pp, psi2savo
      real(dp) :: psidn, psijn, q, r2n
      real(dp) :: r2o, r2sume, risume
      real(dp) :: rminn, rmino, rnorm_nodes, rnorm_nodes_new
      real(dp) :: rnorm_nodes_num, rnorm_nodes_old, ro, taunow
      real(dp) :: tauprim, tratio, v2new, v2old
      real(dp) :: v2sumn, v2sumo, vav2sumn, vav2sumo
      real(dp) :: vavvn, vavvo, vavvt, wtg
      real(dp) :: wtg_derivsum1, wtnow
      real(dp), dimension(3, nelec) :: xstrech
      real(dp), dimension(3) :: xnew
      real(dp), dimension(3, nelec) :: vnew
      real(dp), dimension(3) :: xbac
      real(dp), dimension(3, nelec) :: xdriftedn
      real(dp), dimension(nelec) :: unacp
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: half = .5d0
      real(dp), parameter :: adrift = 0.5d0



      data ncall /0/

c     term=(sqrt(two*pi*tau))**3/pi

      eps_node_cutoff=eps_node_cutoff*sqrt(tau)
      e_cutoff=0.2d0*sqrt(nelec/tau)

      if(idmc.lt.0) then
        expon=1
        dwt=1
      endif
c Undo products
      ipmod=mod(ipass,nfprod)
      ipmod2=mod(ipass+1,nfprod)
      ginv=min(1.d0,tau)
      ffn=eigv*(wdsumo/dmc_nconf)**ginv
      ffi=one/ffn
      fprod=fprod*ffn/ff(ipmod)
      ff(ipmod)=ffn

c Undo weights
      iwmod=mod(ipass,nwprod)

c Store (well behaved velocity/velocity)
      if(ncall.eq.0.and.dmc_irstar.eq.0) then
        do iw=1,nwalk
          do ifr=1,nforce

            vav2sumo=zero
            v2sumo=zero
            do i=1,nelec

c Tau secondary in drift is one (first time around)
              tratio=one

              v2old=vold_dmc(1,i,iw,ifr)**2+vold_dmc(2,i,iw,ifr)**2
     &        +vold_dmc(3,i,iw,ifr)**2
              vavvt=(dsqrt(one+two*adrift*v2old*tau*tratio)-one)/
     &              (adrift*v2old)
              vavvo=vavvt/(tau*tratio)
              vav2sumo=vav2sumo+vavvo*vavvo*v2old
              v2sumo=v2sumo+v2old

              do k=1,3
                xdrifted(k,i,iw,ifr)=xold_dmc(k,i,iw,ifr)+vold_dmc(k,i,iw,ifr)*vavvt
              enddo
            enddo
            fratio(iw,ifr)=dsqrt(vav2sumo/v2sumo)
          enddo
        enddo
        ncall=ncall+1
      endif

      imove=0
      ioldest=0
      ncount_casula=0
      nmove_casula=0
      do iw=1,nwalk
c Loop over primary walker

        call distances(0,xold_dmc(1,1,iw,1))
c Set nuclear coordinates and n-n potential (0 flag = no strech e-coord)
        if(nforce.gt.1)
     &  call strech(xold_dmc(1,1,iw,1),xold_dmc(1,1,iw,1),ajacob,1,0)

        call walkstrdet(iw)
        call walkstrjas(iw)

c Sample Green function for forward move
        r2sume=zero
        risume=zero
        dfus2ac=zero
        dfus2un=zero
        drifdif=zero
        iaccept=0

        if(icasula.eq.3) then
          do i=1,nelec
            imove=0
            call nonloc_grid(i,iw,xnew,psido_dmc(iw,1),imove)
            ncount_casula=ncount_casula+1
            if(imove.gt.0) then
              call psiedmc(i,iw,xnew,psidn,psijn,0)
              nmove_casula=nmove_casula+1

              call compute_determinante_grad(i,psidn,psidn,vnew(1,i),0)
              iaccept=1
              iage(iw)=0
              do k=1,3
                xold_dmc(k,i,iw,1)=xnew(k)
                vold_dmc(k,i,iw,1)=vnew(k,i)
              enddo
              psido_dmc(iw,1)=psidn
              psijo_dmc(iw,1)=psijn
              call jassav(i,0)
              call detsav(i,0)

              call update_ymat(i)
             else
              call distancese_restore(i)
            endif
          enddo
          if(nforce.gt.1.and.istrech.gt.0) then
            do ifr=1,nforce
              call strech(xold_dmc(1,1,iw,1),xstrech,ajacob,ifr,1)
              do k=1,3
                do j=1,nelec
                  xold_dmc(k,j,iw,ifr)=xstrech(k,j)
                enddo
              enddo
            enddo
          endif
        endif

        dwt=1
        pp = 1.
        do i=1,nelec

          if(i.le.nup) then
            iflag_up=2
            iflag_dn=3
           else
            iflag_up=3
            iflag_dn=2
          endif

          call compute_determinante_grad(i,psido_dmc(iw,1),psido_dmc(iw,1),vold_dmc(1,i,iw,1),1)

c Use more accurate formula for the drift
          v2old=vold_dmc(1,i,iw,1)**2+vold_dmc(2,i,iw,1)**2+vold_dmc(3,i,iw,1)**2
c Tau primary -> tratio=one
          vavvt=(dsqrt(one+two*adrift*v2old*tau)-one)/(adrift*v2old)

          dr2=zero
          dfus2o=zero
          do k=1,3
            drift=vavvt*vold_dmc(k,i,iw,1)
            dfus=gauss()*rttau
            dx=drift+dfus
            dr2=dr2+dx**2
            dfus2o=dfus2o+dfus**2
            xnew(k)=xold_dmc(k,i,iw,1)+dx
          enddo

          if(ipr.ge.1) then
            write(ounit,'(''xold_dmc'',2i4,9f8.5)') iw,i,(xold_dmc(k,i,iw,1),k=1,3)
            write(ounit,'(''vold_dmc'',2i4,9f8.5)') iw,i,(vold_dmc(k,i,iw,1),k=1,3)
            write(ounit,'(''psido_dmc'',2i4,9f8.5)') iw,i,psido_dmc(iw,1)
            write(ounit,'(''xnewdr'',2i4,9f8.5)') iw,i,(xnew(k),k=1,3)
          endif

c calculate psi and velocity at new configuration
          call psiedmc(i,iw,xnew,psidn,psijn,0)

          call compute_determinante_grad(i,psidn,psidn,vnew(1,i),0)

          distance_node_ratio2=1.d0
          if(node_cutoff.gt.0) then
            do jel=1,nup
              if(jel.ne.i) call compute_determinante_grad(jel,psidn,psidn,vnew(1,jel),iflag_up)
            enddo

            do jel=nup+1,nelec
              if(jel.ne.i) call compute_determinante_grad(jel,psidn,psidn,vnew(1,jel),iflag_dn)
            enddo

            call nodes_distance(vold_dmc(1,1,iw,1),distance_node,1)
            rnorm_nodes_old=rnorm_nodes_num(distance_node,eps_node_cutoff)/distance_node

            call nodes_distance(vnew,distance_node,0)
            rnorm_nodes_new=rnorm_nodes_num(distance_node,eps_node_cutoff)/distance_node
            distance_node_ratio2=(rnorm_nodes_new/rnorm_nodes_old)**2
          endif

c Check for node crossings
          if(psidn*psido_dmc(iw,1).le.zero) then
            nodecr=nodecr+1
            if(icross.le.0) then
              p=zero
              goto 160
            endif
          endif

c Calculate Green function for the reverse move

          v2new=vnew(1,i)**2+vnew(2,i)**2+vnew(3,i)**2
          vavvt=(dsqrt(one+two*adrift*v2new*tau)-one)/(adrift*v2new)

          dfus2n=zero
          do k=1,3
            drift=vavvt*vnew(k,i)
            xbac(k)=xnew(k)+drift
            dfus=xbac(k)-xold_dmc(k,i,iw,1)
            dfus2n=dfus2n+dfus**2
          enddo

          if(ipr.ge.1) then
            write(ounit,'(''xold_dmc'',9f10.6)')(xold_dmc(k,i,iw,1),k=1,3),
     &      (xnew(k),k=1,3), (xbac(k),k=1,3)
            write(ounit,'(''dfus2o'',9f10.6)')dfus2o,dfus2n,
     &      psido_dmc(iw,1),psidn,psijo_dmc(iw,1),psijn
          endif

          p=(psidn/psido_dmc(iw,1))**2*exp(2*(psijn-psijo_dmc(iw,1)))*
     &    exp((dfus2o-dfus2n)/(two*tau))*distance_node_ratio2

          if(ipr.ge.1) write(ounit,'(''p'',11f10.6)')
     &    p,(psidn/psido_dmc(iw,1))**2*exp(2*(psijn-psijo_dmc(iw,1))),
     &    exp((dfus2o-dfus2n)/(two*tau)),psidn,psido_dmc(iw,1),
     &    psijn,psijo_dmc(iw,1),dfus2o,dfus2n

c The following is one reasonable way to cure persistent configurations
c Not needed if itau_eff <=0 and in practice we have never needed it even
c otherwise
          if(iage(iw).gt.50) p=p*1.1d0**(iage(iw)-50)

          pp=pp*p
          ! pp = p*p ! <= double check that with Claudia
          p=dmin1(one,p)
  160     q=one-p

          acc=acc+p
          trymove=trymove+1
          dfus2ac=dfus2ac+p*dfus2o
          dfus2un=dfus2un+dfus2o
          dr2ac=dr2ac+p*dr2
          dr2un=dr2un+dr2

c Calculate density and moments of r for primary walk
          r2o=zero
          r2n=zero
          rmino=zero
          rminn=zero
          do k=1,3
            r2o=r2o+xold_dmc(k,i,iw,1)**2
            r2n=r2n+xnew(k)**2
            rmino=rmino+(xold_dmc(k,i,iw,1)-cent(k,1))**2
            rminn=rminn+(xnew(k)-cent(k,1))**2
          enddo
          rmino=sqrt(rmino)
          rminn=sqrt(rminn)
          itryo(i)=min(int(delri*rmino)+1,nrad)
          itryn(i)=min(int(delri*rminn)+1,nrad)

c If we are using weights rather than accept/reject
          if(iacc_rej.le.0) then
            p=one
            q=zero
          endif

          iacc_elec(i)=0
          if(rannyu(0).lt.p) then
            iaccept=1
            nacc=nacc+1
            iacc_elec(i)=1
            if(ipq.le.0) p=one

            iage(iw)=0
            do k=1,3
              drifdif=drifdif+(xold_dmc(k,i,iw,1)-xnew(k))**2
              xold_dmc(k,i,iw,1)=xnew(k)
            enddo
            psido_dmc(iw,1)=psidn
            psijo_dmc(iw,1)=psijn
            call jassav(i,0)
            call detsav(i,0)

           else
            if(ipq.le.0) p=zero
            call distancese_restore(i)
          endif
          q=one-p

c Calculate moments of r and save rejection probability for primary walk
          r2sume=r2sume+(q*r2o+p*r2n)
          risume=risume+(q/dsqrt(r2o)+p/dsqrt(r2n))
          unacp(i)=q

          call update_ymat(i)

        enddo

c Effective tau for branching
        tauprim=tau*dfus2ac/dfus2un

        do ifr=1,nforce

          if(ifr.eq.1) then
c Primary configuration
            if(icasula.lt.0) i_vpsp=icasula
            drifdifr=one
            if(nforce.gt.1)
     &      call strech(xold_dmc(1,1,iw,1),xold_dmc(1,1,iw,1),ajacob,1,0)
            call hpsi(xold_dmc(1,1,iw,1),psidn,psijn,enew,ipass,1)
            call walksav_det(iw)
            call walksav_jas(iw)
            if(icasula.lt.0) call multideterminant_tmove(psidn,0)
c           call t_vpsp_sav(iw)
            call t_vpsp_sav
            i_vpsp=0
            rnorm_nodes=1.d0
            if(node_cutoff.gt.0) then
              call nodes_distance(vold_dmc(1,1,iw,1),distance_node,1)
              rnorm_nodes=rnorm_nodes_num(distance_node,eps_node_cutoff)/distance_node
            endif
           else
c Secondary configuration
            if(istrech.eq.0) then
              call strech(xold_dmc(1,1,iw,ifr),xold_dmc(1,1,iw,ifr),ajacob,ifr,0)
              drifdifr=one
c No streched positions for electrons
              do i=1,nelec
                do k=1,3
                  xold_dmc(k,i,iw,ifr)=xold_dmc(k,i,iw,1)
                enddo
              enddo
              ajacold(iw,ifr)=one
             else
c Compute streched electronic positions for all nucleus displacement
              call strech(xold_dmc(1,1,iw,1),xstrech,ajacob,ifr,1)
              drifdifs=zero
              do i=1,nelec
                do k=1,3
                  drifdifs=drifdifs+(xstrech(k,i)-xold_dmc(k,i,iw,ifr))**2
                  xold_dmc(k,i,iw,ifr)=xstrech(k,i)
                enddo
              enddo
              ajacold(iw,ifr)=ajacob
              if(drifdif.eq.0.d0) then
                drifdifr=one
               else
                drifdifr=drifdifs/drifdif
              endif
            endif
            if(icasula.lt.0) i_vpsp=icasula
            call hpsi(xold_dmc(1,1,iw,ifr),psidn,psijn,enew,ipass,ifr)
            i_vpsp=0
          endif

          do i=1,nelec
              call compute_determinante_grad(i,psidn,psidn,vold_dmc(1,i,iw,ifr),1)
          enddo

          vav2sumn=zero
          v2sumn=zero
          do i=1,nelec

c Use more accurate formula for the drift and tau secondary in drift
            tratio=one
            if(ifr.gt.1.and.itausec.eq.1) tratio=drifdifr

            v2old=vold_dmc(1,i,iw,ifr)**2+vold_dmc(2,i,iw,ifr)**2
     &      +vold_dmc(3,i,iw,ifr)**2
            vavvt=(dsqrt(one+two*adrift*v2old*tau*tratio)-one)/
     &          (adrift*v2old)
            vavvn=vavvt/(tau*tratio)

            vav2sumn=vav2sumn+vavvn**2*v2old
            v2sumn=v2sumn+v2old

            do k=1,3
              xdriftedn(k,i)=xold_dmc(k,i,iw,ifr)+vold_dmc(k,i,iw,ifr)*vavvt
            enddo
          enddo
          fration=dsqrt(vav2sumn/v2sumn)

          taunow=tauprim*drifdifr

          if(ipr.ge.1)write(ounit,'(''wt'',9f10.5)') wt(iw),etrial,eest

          if(icut_e.eq.0) then
            ewto=eest-(eest-eold(iw,ifr))*fratio(iw,ifr)
            ewtn=eest-(eest-enew)*fration
           else
            deo=eest-eold(iw,ifr)
            den=eest-enew
            ewto=eest-sign(1.d0,deo)*min(e_cutoff,dabs(deo))
            ewtn=eest-sign(1.d0,den)*min(e_cutoff,dabs(den))
          endif

          if(idmc.gt.0) then
            expon=(etrial-half*(ewto+ewtn))*taunow
            if(icut_br.le.0) then
              dwt=dexp(expon)
             else
              dwt=0.5d0+1/(1+exp(-4*expon))
            endif
          endif

c If we are using weights rather than accept/reject
          if(iacc_rej.eq.0) dwt=dwt*pp

c Exercise population control if dmc or vmc with weights
          if(idmc.gt.0.or.iacc_rej.eq.0) dwt=dwt*ffi

          drifdifgfunc=1.d0
          if(nforce.gt.1.and.idrifdifgfunc.gt.0) then
            dfus=0
            do i=1,nelec
              if(iacc_elec(i).gt.0) then
                do k=1,3
                  dfus=dfus+(xold_dmc(k,i,iw,ifr)-xdrifted(k,i,iw,ifr))**2
                enddo
              endif
            enddo
            dfus=0.5d0*dfus/tau
            drifdifgfunc=-dfus
c           drifdifgfunc=drifdifgfunc-2*log(rnorm_nodes)
          endif
c Set weights and product of weights over last nwprod steps
          if(ifr.eq.1) then

            wt(iw)=wt(iw)*dwt
            wtnow=wt(iw)
            pwt(iw,ifr)=pwt(iw,ifr)+log(dwt)+drifdifgfunc-wthist(iw,iwmod,ifr)
            wthist(iw,iwmod,ifr)=dlog(dwt)+drifdifgfunc

           elseif(ifr.gt.1) then

            ro=log(ajacold(iw,ifr))
            if(idrifdifgfunc.eq.0) ro=0.d0
            pwt(iw,ifr)=pwt(iw,ifr)+dlog(dwt)+drifdifgfunc+ro-wthist(iw,iwmod,ifr)
            wthist(iw,iwmod,ifr)=dlog(dwt)+drifdifgfunc+ro
            wtnow=wt(iw)*dexp(pwt(iw,ifr)-pwt(iw,1))

          endif

          wtnow=wtnow/rnorm_nodes**2
c         if(idrifdifgfunc.eq.0)wtnow=wtnow/rnorm_nodes**2

          if(ipr.ge.1)write(ounit,'(''eold,enew,wt'',9f10.5)')
     &    eold(iw,ifr),enew,wtnow

          if(idmc.gt.0) then
            wtg=wtnow*fprod
           else
            wtg=wtnow
          endif

          if(ifr.eq.1) then
            r2sum=r2sum+wtg*r2sume
            risum=risum+wtg*risume
            do i=1,nelec
              rprob(itryo(i))=rprob(itryo(i))+wtg*unacp(i)
              rprob(itryn(i))=rprob(itryn(i))+wtg*(one-unacp(i))
            enddo
          endif
          tausum(ifr)=tausum(ifr)+wtg*taunow

          if(dabs((enew-etrial)/etrial).gt.0.2d+0) then
           write(18,'(i6,f8.2,2d10.2,(8f8.4))') ipass,
     &     enew-etrial,psidn,psijn,(xnew(ii),ii=1,3)
          endif

          if(wt(iw).gt.3) write(18,'(i6,i4,3f8.2,30f8.4)') ipass,iw,
     &    wt(iw),enew-etrial,eold(iw,ifr)-etrial,(xnew(ii),ii=1,3)

          eold(iw,ifr)=enew
          ! peo_dmc(iw,ifr)=pen <= pen is undefined, could be den mispelled ?
          if(icut_e .ne. 0) then  ! <= check that with claudia
            peo_dmc(iw,ifr)=den
          else
            peo_dmc(iw,ifr)=0.0_dp
          end if
          ! d2o(iw,ifr)=d2n <= dn2 is not initialized I don't think it' s mispelled
          d2o(iw,ifr)=0.0_dp ! TODO : Old stuff to be fixed in the future
          psido_dmc(iw,ifr)=psidn
          psijo_dmc(iw,ifr)=psijn
          fratio(iw,ifr)=fration
          do i=1,nelec
            do k=1,3
              xdrifted(k,i,iw,ifr)=xdriftedn(k,i)
            enddo
          enddo
          call prop_save_dmc(iw)
          call pcm_save(iw)
          call mmpol_save(iw)

          if(ifr.eq.1) then
            if(iaccept.eq.0) then
              iage(iw)=iage(iw)+1
              ioldest=max(ioldest,iage(iw))
              ioldestmx=max(ioldestmx,iage(iw))
            endif

            psi2savo=2*(dlog(dabs(psido_dmc(iw,1)))+psijo_dmc(iw,1))

            wsum1(ifr)=wsum1(ifr)+wtnow
            esum1_dmc(ifr)=esum1_dmc(ifr)+wtnow*eold(iw,ifr)
            pesum_dmc(ifr)=pesum_dmc(ifr)+wtg*peo_dmc(iw,ifr)
            tpbsum_dmc(ifr)=tpbsum_dmc(ifr)+wtg*(eold(iw,ifr)-peo_dmc(iw,ifr))
            tjfsum_dmc(ifr)=tjfsum_dmc(ifr)-wtg*half*hb*d2o(iw,ifr)

            derivsum(1,ifr)=derivsum(1,ifr)+wtg*eold(iw,ifr)

            if(idrifdifgfunc.gt.0) then
              derivsum(2,ifr)=derivsum(2,ifr)+wtg*eold(iw,ifr)*pwt(iw,ifr)
              derivsum(3,ifr)=derivsum(3,ifr)+wtg*pwt(iw,ifr)
             else
              derivsum(2,ifr)=derivsum(2,ifr)+wtg*eold(iw,ifr)*(pwt(iw,ifr)+psi2savo)
              derivsum(3,ifr)=derivsum(3,ifr)+wtg*(pwt(iw,ifr)+psi2savo)
            endif

            call prop_sum_dmc(0.d0,wtg,iw)
            call pcm_sum(0.d0,wtg,iw)
            call mmpol_sum(0.d0,wtg,iw)

            call optjas_sum(wtg,0.d0,eold(iw,1),eold(iw,1),0)
            call optorb_sum(wtg,0.d0,eold(iw,1),eold(iw,1),0)
            call optci_sum(wtg,0.d0,eold(iw,1),eold(iw,1))

            call optx_jas_orb_sum(wtg,0.d0,0)
            call optx_jas_ci_sum(wtg,0.d0,eold(iw,1),eold(iw,1))
            call optx_orb_ci_sum(wtg,0.d0)

           else
c           write(ounit,*) 'IN DMC',ajacold(iw,ifr)
            ro=1.d0
            if(idrifdifgfunc.eq.0) ro=ajacold(iw,ifr)*psido_dmc(iw,ifr)**2*exp(2*psijo_dmc(iw,ifr)-psi2savo)

            wsum1(ifr)=wsum1(ifr)+wtnow*ro
            esum1_dmc(ifr)=esum1_dmc(ifr)+wtnow*eold(iw,ifr)*ro
            pesum_dmc(ifr)=pesum_dmc(ifr)+wtg*peo_dmc(iw,ifr)*ro
            tpbsum_dmc(ifr)=tpbsum_dmc(ifr)+wtg*(eold(iw,ifr)-peo_dmc(iw,ifr))*ro
            tjfsum_dmc(ifr)=tjfsum_dmc(ifr)-wtg*half*hb*d2o(iw,ifr)*ro

            wtg=wt(iw)*fprod/rnorm_nodes**2
            wtg_derivsum1=wtg
c           if(idrifdifgfunc.eq.0)then
c             wtg=wt(iw)*fprod/rnorm_nodes**2
c             wtg_derivsum1=wtg
c            else
c             wtg=wt(iw)*fprod
c             wtg_derivsum1=wtg/rnorm_nodes**2
c           endif

            derivsum(1,ifr)=derivsum(1,ifr)+wtg_derivsum1*eold(iw,ifr)

            if(idrifdifgfunc.gt.0) then
              derivsum(2,ifr)=derivsum(2,ifr)+wtg*eold(iw,1)*pwt(iw,ifr)
              derivsum(3,ifr)=derivsum(3,ifr)+wtg*pwt(iw,ifr)
            else
              ro=log(ajacold(iw,ifr))+2*(log(abs(psido_dmc(iw,ifr)))+psijo_dmc(iw,ifr))
              derivsum(2,ifr)=derivsum(2,ifr)+wtg*eold(iw,1)*(pwt(iw,ifr)+ro)
              derivsum(3,ifr)=derivsum(3,ifr)+wtg*(pwt(iw,ifr)+ro)
            endif
          endif

        enddo
c       write(ounit,*) 'IN DMC',ajacold(iw,2)

c       wtg=wt(iw)*fprod/rnorm_nodes**2
c       write(*,*)'prima ',wtg,eold(iw,2),pwt(iw,2),ajacold(iw,2),psido_dmc(iw,2),psijo_dmc(iw,2),idrifdifgfunc
c       call deriv(wtg,eold(iw,1),pwt(iw,1),ajacold(iw,1),psido_dmc(iw,1),psijo_dmc(iw,1),idrifdifgfunc)
c       call deriv(wtg,eold,pwt,ajacold,psido_dmc,psijo_dmc,idrifdifgfunc,iw,mwalk)

        if(icasula.eq.-1) then

c Set nuclear coordinates (0 flag = no strech e-coord)
          if(nforce.gt.1)
     &    call strech(xold_dmc(1,1,iw,1),xold_dmc(1,1,iw,1),ajacob,1,0)

          call walkstrdet(iw)
          call walkstrjas(iw)
c         call t_vpsp_get(iw)
          call t_vpsp_get

          imove=0
          call nonloc_grid(iel,iw,xnew,psido_dmc(iw,1),imove)

          ncount_casula=ncount_casula+1
          if(imove.gt.0) then
            call psiedmc(iel,iw,xnew,psidn,psijn,0)
            nmove_casula=nmove_casula+1

c           call compute_determinante_grad(iel,psidn,psidn,vnew(1,iel),0)

            iage(iw)=0
            do k=1,3
              xold_dmc(k,iel,iw,1)=xnew(k)
            enddo
c 290         vold_dmc(k,iel,iw,1)=vnew(k,iel)
            psido_dmc(iw,1)=psidn
            psijo_dmc(iw,1)=psijn
            call jassav(iel,0)
            call detsav(iel,0)

            if(iel.le.nup) call update_ymat(nup)
            if(iel.gt.nup) call update_ymat(nelec)

            call walksav_det(iw)
            call walksav_jas(iw)
            if(nforce.gt.1.and.istrech.gt.0) then
              do ifr=1,nforce
                call strech(xold_dmc(1,1,iw,1),xstrech,ajacob,ifr,1)
                do k=1,3
                  do i=1,nelec
                     xold_dmc(k,i,iw,ifr)=xstrech(k,i)
                  enddo
                enddo
              enddo
            endif


          endif
        endif

      call average(1)
      enddo

      if(wsum1(1).gt.1.1d0*dmc_nconf) write(18,'(i6,9d12.4)') ipass,ffn,fprod,
     &fprod/ff(ipmod2),wsum1(1),wgdsumo

      if(idmc.gt.0.or.iacc_rej.eq.0) then
        wfsum1=wsum1(1)*ffn
        efsum1=esum1_dmc(1)*ffn
      endif
      do ifr=1,nforce
        if(idmc.gt.0.or.iacc_rej.eq.0) then
          wgsum1(ifr)=wsum1(ifr)*fprod
          egsum1(ifr)=esum1_dmc(ifr)*fprod
         else
          wgsum1(ifr)=wsum1(ifr)
          egsum1(ifr)=esum1_dmc(ifr)
        endif
      enddo

      call splitj
      if(icasula.eq.0) ncount_casula=1
      if(ipr.gt.-2) write(11,'(i8,f9.6,f12.5,f11.6,i5,f11.5)') ipass,ffn,
     &wsum1(1),esum1_dmc(1)/wsum1(1),nwalk
     &,float(nmove_casula)/float(ncount_casula)

      return
      end
