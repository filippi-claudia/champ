      subroutine hpsi(coord,psid,psij,energy,ipass,ifr)

c Written by Cyrus Umrigar, modified by Claudia Filippi and A. Scemama
c modified by Claudio Amovilli and Franca Floris for PCM and QM-MMPOl

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pseudo.h'
      include 'force.h'
      include 'pcm.h'
      include 'mmpol.h'
      include 'efield.h'
      include 'optjas.h'
      include 'mstates.h'

c Calculates energy

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contr2/ ijas,icusp,icusp2,isc,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
      common /elec/ nup,ndn

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb

      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT)
     &,r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)

      common /slater/ slmui(MMAT_DIM),slmdi(MMAT_DIM)
     &,fpu(3,MMAT_DIM),fpd(3,MMAT_DIM)
     &,fppu(MMAT_DIM),fppd(MMAT_DIM)
     &,ddx(3,MELEC),d2dx2(MELEC)
      common /multislater/ detu(MDET),detd(MDET)

      common /jaso/ fso(MELEC,MELEC),fijo(3,MELEC,MELEC)
     &,d2ijo(MELEC,MELEC),d2o,fsumo,fjo(3,MELEC)

      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,lpot(MCTYPE),nloc

      common /casula/ t_vpsp(MCENT,MPS_QUAD,MELEC),icasula,i_vpsp
      common /pcm_hpsi/ pepcms,pepcmv,qopcm,enfpcm(MCHS)
      common /mmpol_hpsi/ peQMdp,peQMq,eek_pol(3,MCHMM)

      common /velocity_jastrow/vj(3,MELEC),vjn(3,MELEC)

      common /multidet/ kref,numrep_det(MDET,2),irepcol_det(MELEC,MDET,2),ireporb_det(MELEC,MDET,2)
     & ,iwundet(MDET,2),iactv(2),ivirt(2)

      common /Bloc/ b(MORB,MELEC),xmatu(MELEC**2),xmatd(MELEC**2)
     & ,tildem(MELEC,MORB,2)

      common /ycompact/ ymat(MORB,MELEC,2,MSTATES),dymat(MORB,MELEC,2,MSTATES)

      common /force_analy/ iforce_analy

      dimension coord(3,*),psid(*),energy(*)
      dimension denergy(MSTATES),eloc_det(MDET,2),vpsp_det(2),dvpsp_dj(MPARMJ)

      iwf=iwftype(ifr)

c pe_ee computed in distances (called from distances if iperiodic != 0)
c pe_en(loc) computed in distances if nloc=0 or called from distances if iperiodic!=0
c pe_en(loc) computed in nonloc_pot if nloc!=0 and iperiodic=0
c pe_en(nonloc) computed in nonloc_pot if nloc !=0

c distances needed for Jastrow, determinants, and potential energy
      call distances(0,coord)

c local potential contributions
      call pot_local(pe_local)

c external potential on a grid (e.g. MM from CPMD)
      call p2gtid('qmmm:iqmmm',iqmmm,0,1)
      if(iqmmm.eq.1) then
        ext_pot=0
        call qmmm_extpot_ene(coord,nelec,ext_pot)
        pe_local=pe_local+ext_pot
      endif

c external charges
      if(iefield.eq.1) then
        ext_pot=0
        call efield_extpot_ene(coord,nelec,ext_pot)
        pe_local=pe_local+ext_pot
      endif

c PCM polarization charges 
      if(ipcm.gt.1) then
        pepcms=0
        pepcmv=0
        call pcm_extpot_ene(coord,nelec,pepcms,pepcmv)
        pepcm=pepcms+pepcmv
        pe_local=pe_local+pepcm
      endif

c QM-MMPOL (charges+induced dipoles) 
      if(immpol.gt.1) then
        peQMdp=0
        peQMq=0
        call mmpol_extpot_ene(coord,nelec,peQMdp,peQMq)
        peQM=peQMdp+peQMq
        pe_local=pe_local+peQM
      endif

      if(ipr.ge.3) write(6,'(''pe_loc before nonloc_pot'',9f12.5)') pe_local

c get contribution from jastrow (also compute derivatives wrt parameters and nuclei)
      if(ianalyt_lap.eq.1) then
        call jastrow(coord,vj,d2j,psij,ifr)
       else
        call jastrow_num(coord,vj,d2j,psij)
      endif
      if(ipr.ge.3) write(6,'(''d2j,psij'',9f12.5)') d2j,psij

c compute reference determinant, its derivatives, and kinetic contribution to B_eloc and its derivatives
      icheck=0
  2   call determinant(coord,rvec_en,r_en)

      icheck=icheck+1
      call check_detref(ipass,icheck,newref)
      if(newref.gt.0) then
        if (ioptorb.ne.0) call optorb_define
        goto 2
      endif

      call compute_bmatrices_kin

c compute pseudo-potential contribution
c nonloc_pot must be called after determinant because slater matrices are needed
      if(nloc.gt.0) 
     &  call nonloc_pot(coord,rshift,rvec_en,r_en,pe_local,vpsp_det,dvpsp_dj,t_vpsp,i_vpsp,ifr)

      if(ipr.ge.3) then 
        write(6,'(''pe_loc after nonloc_pot'',9f12.5)') pe_local
        write(6,'(''pe_ref after nonloc_pot'',9f12.5)') (vpsp_det(ii),ii=1,2)
      endif

      call multideterminant_hpsi(vj,vpsp_det,eloc_det)

      e_other=pe_local-hb*d2j
      do 10 i=1,nelec
   10   e_other=e_other-hb*(vj(1,i)**2+vj(2,i)**2+vj(3,i)**2)

      do 30 istate=1,nstates

c combine determinantal quantities to obtain trial wave function
        call determinant_psit(psid(istate),istate)

c compute energy using Ymat
        denergy(istate)=0
        do 20 iab=1,2
          nel=nup
          if(iab.eq.2) nel=ndn
          do 20 jrep=ivirt(iab),norb
            do 20 irep=iactv(iab),nel
   20         denergy(istate)=denergy(istate)+ymat(jrep,irep,iab,istate)*tildem(irep,jrep,iab)
        denergy(istate)=denergy(istate)*detu(kref)*detd(kref)/psid(istate)

        energy(istate)=denergy(istate)+eloc_det(kref,1)+eloc_det(kref,2)+e_other

        if(ipr.ge.2) then
          write(6,'(''state'',i4)') istate
          write(6,'(''psid,psij'',9d12.5)') psid(istate),psij
          write(6,'(''psitot   '',e18.11)') psid(istate)*exp(psij)
          do k=1,ndet
c           write(6,'(''psitot_k '',i6,3e18.8)') k, detu(k),detd(k),detu(k)*detd(k)*exp(psij)
            write(6,'(''psitot_k '',i6,3e18.8)') k, detu(k),detd(k),cdet(k,1,1)*detu(k)*detd(k)*exp(psij)
          enddo
c         do 25 i=1,nelec
c           do 25 k=1,3
c  25         write(6,'(''vj'',2e18.11)') vj(k,i)
          if(ipr.ge.3) write(6,'(''energy'',9f16.10)') energy(istate)
        endif

   30 continue

      if(ifr.eq.1) then
        if(iforce_analy.eq.1) call compute_force(psid(1),denergy(1))

        call optorb_compute(psid,energy,denergy)

        call optjas_deloc(psid,energy,dvpsp_dj,vj)

        call optci_deloc(eloc_det,e_other,psid,energy)

        call prop_compute(coord)
      endif

      return
      end
