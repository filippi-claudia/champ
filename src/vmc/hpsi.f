      subroutine hpsi(coord,psid,psij,energy,ipass,ifr)
c     Written by Cyrus Umrigar, modified by Claudia Filippi and A. Scemama
c     modified by Claudio Amovilli and Franca Floris for PCM and QM-MMPOl
      use optjas, only: MPARMJ
      use vmc_mod, only: MELEC, MDET, MCENT
      use vmc_mod, only: MMAT_DIM, MMAT_DIM2
      use const, only: hb, nelec, ipr
      use mstates_mod, only: MSTATES
      use csfs, only: nstates
      use dets, only: cdet, ndet
      use elec, only: ndn, nup
      use mmpol_hpsi, only: peQMdp, peQMq
      use multidet, only: iactv, ivirt, kref
      use pcm_hpsi, only: pepcms, pepcmv
      use wfsec, only: iwf, iwftype
      use ycompact, only: ymat
      use casula, only: i_vpsp, t_vpsp
      use coefs, only: norb
      use contr2, only: ianalyt_lap
      use Bloc, only: tildem
      use force_analy, only: iforce_analy
      use pseudo, only: nloc
      use velocity_jastrow, only: vj
      use mmpol_cntrl, only: immpol
      use efield, only: iefield
      use distance_mod, only: rshift, r_en, rvec_en
      use pcm_cntrl, only: ipcm
      use slater, only: d2dx2, ddx, fp, fpp, slmi
      use distance_mod, only: rshift, r_en, rvec_en, r_ee, rvec_ee
      use multislater, only: detiab
      use inputflags, only: iqmmm

      implicit real*8(a-h,o-z)

c     Calculates energy

      dimension coord(3,*),psid(*),energy(*),psij(*),d2j(MSTATES)
      dimension denergy(MSTATES),eloc_det(MDET,2,MSTATES),
     &     vpsp_det(2,MSTATES),dvpsp_dj(MPARMJ)

      iwf=iwftype(ifr)

c     pe_ee computed in distances (called from distances if iperiodic != 0)
c     pe_en(loc) computed in distances if nloc=0 or called from distances if iperiodic!=0
c     pe_en(loc) computed in nonloc_pot if nloc!=0 and iperiodic=0
c     pe_en(nonloc) computed in nonloc_pot if nloc !=0

c     distances needed for Jastrow, determinants, and potential energy
      call distances(0,coord)

c     local potential contributions
      call pot_local(pe_local)

c     external potential on a grid (e.g. MM from CPMD)
      if(iqmmm.eq.1) then
         ext_pot=0.0d0
         call qmmm_extpot_ene(coord,nelec,ext_pot)
         pe_local=pe_local+ext_pot
      endif

c     external charges
      if(iefield.eq.1) then
         ext_pot=0.0d0
         call efield_extpot_ene(coord,nelec,ext_pot)
         pe_local=pe_local+ext_pot
      endif

c     PCM polarization charges 
      if(ipcm.gt.1) then
         pepcms=0.0d0
         pepcmv=0.0d0
         call pcm_extpot_ene(coord,nelec,pepcms,pepcmv)
         pepcm=pepcms+pepcmv
         pe_local=pe_local+pepcm
      endif

c     QM-MMPOL (charges+induced dipoles) 
      if(immpol.gt.1) then
         peQMdp=0.0d0
         peQMq=0.0d0
         call mmpol_extpot_ene(coord,nelec,peQMdp,peQMq)
         peQM=peQMdp+peQMq
         pe_local=pe_local+peQM
      endif

      if(ipr.ge.3) write(6,'(''pe_loc before nonloc_pot'',9f12.5)') pe_local

c     get contribution from jastrow (also compute derivatives wrt parameters and nuclei)
      do istate=1,nstates
         if(ianalyt_lap.eq.1) then
            call jastrow(coord,vj(:,:,istate),d2j(istate),psij(istate),ifr,istate)
         else
            call jastrow_num(coord,vj(:,:,istate),d2j(istate),psij(istate))
         endif
         if(ipr.ge.3) then
            write(6, *) "STATE", istate
            write(6,'(''d2j,psij'',9f12.5)') d2j(istate),psij(istate)
         endif
      enddo

c     compute reference determinant, its derivatives,
c     and kinetic contribution to B_eloc and its derivatives

      call determinant(ipass,coord,rvec_en,r_en)
      call compute_bmatrices_kin

c     compute pseudo-potential contribution
c     nonloc_pot must be called after determinant because slater matrices are needed

      if(nloc.gt.0) then
         call nonloc_pot(coord,rshift,rvec_en,r_en,pe_local,
     &        vpsp_det,dvpsp_dj,t_vpsp,i_vpsp,ifr)
      endif

      if(ipr.ge.3) then 
         write(6,'(''pe_loc after nonloc_pot'',9f12.5)') pe_local
         do istate=1,nstates
            write(6,*) "STATE", istate
            write(6,'(''pe_ref after nonloc_pot'',9f12.5)') (vpsp_det(ii,istate),ii=1,2)
         enddo
      endif

      call multideterminant_hpsi(vj,vpsp_det,eloc_det)

      do istate=1,nstates
         e_other=pe_local-hb*d2j(istate)
         do i=1,nelec
            e_other=e_other-hb*(vj(1,i,istate)**2+vj(2,i,istate)**2+vj(3,i,istate)**2)
         enddo

c     combine determinantal quantities to obtain trial wave function
         call determinant_psit(psid(istate),istate)
c     compute energy using Ymat
         denergy(istate)=0
         do iab=1,2
            nel=nup
            if(iab.eq.2) nel=ndn
            do jrep=ivirt(iab),norb
               do irep=iactv(iab),nel
                  denergy(istate)=denergy(istate)+ymat(jrep,irep,iab,istate)
     &                 *tildem(irep,jrep,iab,istate)
               enddo 
            enddo 
         enddo

         denergy(istate)=denergy(istate)*detiab(kref,1,istate)
     &        *detiab(kref,2,istate)/psid(istate)

         energy(istate)=denergy(istate)+eloc_det(kref,1,istate)
     &        +eloc_det(kref,2,istate)+e_other

         if(ipr.ge.2) then
            write(6,'(''state'',i4)') istate
            write(6,'(''psid,psij'',9d12.5)') psid(istate),psij(istate)
            write(6,'(''psitot   '',e18.11)') psid(istate)*exp(psij(istate))
            if(ipr.ge.3) write(6,'(''energy'',9f16.10)') energy(istate)
         endif
      enddo

      if(ifr.eq.1) then
         if(iforce_analy.eq.1) call compute_force(psid(1),denergy(1))
         call optorb_compute(psid,energy,denergy)
         call optjas_deloc(psid,energy,dvpsp_dj,vj)
         call optci_deloc(eloc_det,e_other,psid,energy)
         call prop_compute(coord)
      endif

      end subroutine
