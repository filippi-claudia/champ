      module hpsi_mod
      contains
      subroutine hpsi(coord,psid,psij,ekin,energy,ipass,ifr)

c Written by Cyrus Umrigar, modified by Claudia Filippi and A. Scemama
c modified by Claudio Amovilli and Franca Floris for PCM and QM-MMPOl

      use Bloc,    only: tildem,tildemkin
      use casula,  only: i_vpsp,t_vpsp
      use constants, only: hb
      use contrl_file, only: ounit
      use control, only: ipr
      use csfs,    only: nstates
      use determinant_mod, only: compute_bmatrices_kin,determinant
      use determinant_psit_mod, only: determinant_psit
      use distance_mod, only: r_en,rshift,rvec_en
      use distances_mod, only: distances
      use efield,  only: iefield
      use efield_f_mod, only: efield_extpot_ene
      use force_analytic, only: compute_force
      use inputflags, only: iqmmm
      use jastrow, only: ianalyt_lap
      use jastrow_mod, only: jastrow_f => jastrow
      use jastrow_num_mod, only: jastrow_num
      use m_force_analytic, only: iforce_analy
      use mmpol,   only: mmpol_extpot_ene
      use mmpol_cntrl, only: immpol
      use mmpol_hpsi, only: peQMdp,peQMq
      use mstates_mod, only: MSTATES
      use multidet, only: iactv,ivirt
      use multideterminant_mod, only: multideterminant_hpsi
      use multiple_geo, only: iwf,iwftype
      use multislater, only: detiab
      use nonloc_pot_mod, only: nonloc_pot
      use optci_mod, only: optci_deloc
      use optjas_mod, only: optjas_deloc
      use optorb_f_mod, only: optorb_compute
      use optwf_parms, only: nparmj
      use pcm_cntrl, only: ipcm
      use pcm_hpsi, only: pcms,pcmv
      use pcm_mod, only: pcm_extpot_ene
      use pot_local_mod, only: pot_local
      use precision_kinds, only: dp
      use properties_mod, only: prop_compute
      use pseudo,  only: nloc
      use qmmm_pot, only: qmmm_extpot_ene
      use slater,  only: kref,ndet, norb
      use system,  only: ndn,nelec,nup
      use velocity_jastrow, only: vj
      use ycompact, only: ymat

      implicit none

      integer :: i, iab, ifr, ii, ipass
      integer :: irep, istate, jrep, nel
      real(dp) :: d2j, e_other, ekin_other, ext_pot, peQM, pe_local
      real(dp) :: pepcm, psij
      real(dp), dimension(3, *) :: coord
      real(dp), dimension(*) :: psid
      real(dp), dimension(*) :: energy
      real(dp), dimension(*) :: ekin
      real(dp), dimension(MSTATES) :: dekin
      real(dp), dimension(MSTATES) :: denergy
      real(dp), dimension(ndet, 2) :: eloc_det
      real(dp), dimension(2) :: ekin_det
      real(dp), dimension(2) :: vpsp_det
      real(dp), dimension(nparmj) :: dvpsp_dj

c Calculates energy


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
        pcms=0
        pcmv=0
        call pcm_extpot_ene(coord,nelec,pcms,pcmv)
        pepcm=pcms+pcmv
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

      if(ipr.ge.3) write(ounit,'(''pe_loc before nonloc_pot'',9f12.5)') pe_local

c get contribution from jastrow (also compute derivatives wrt parameters and nuclei)
      if(ianalyt_lap.eq.1) then
        call jastrow_f(coord,vj,d2j,psij,ifr)
       else
        call jastrow_num(coord,vj,d2j,psij)
      endif
      if(ipr.ge.3) write(ounit,'(''d2j,psij'',9f12.5)') d2j,psij

c compute reference determinant, its derivatives, and kinetic contribution to B_eloc and its derivatives
      call determinant(ipass,coord,rvec_en,r_en)

      call compute_bmatrices_kin

c compute pseudo-potential contribution
c nonloc_pot must be called after determinant because slater matrices are needed

      if(nloc.gt.0)
     &  call nonloc_pot(coord,rshift,rvec_en,r_en,pe_local,vpsp_det,dvpsp_dj,t_vpsp,i_vpsp,ifr)

      if(ipr.ge.3) then
        write(ounit,'(''pe_loc after nonloc_pot'',9f12.5)') pe_local
        write(ounit,'(''pe_ref after nonloc_pot'',9f12.5)') (vpsp_det(ii),ii=1,2)
      endif

      call multideterminant_hpsi(vj,ekin_det,vpsp_det,eloc_det)

      ekin_other=-hb*d2j
      do i=1,nelec
        ekin_other=ekin_other-hb*(vj(1,i)**2+vj(2,i)**2+vj(3,i)**2)
      enddo
      e_other=ekin_other+pe_local

      do istate=1,nstates
c combine determinantal quantities to obtain trial wave function
        call determinant_psit(psid(istate),istate)
c compute energy using Ymat
        denergy(istate)=0
        dekin(istate)=0
        do iab=1,2
          nel=nup
          if(iab.eq.2) nel=ndn
          do jrep=ivirt(iab),norb
            do irep=iactv(iab),nel
              denergy(istate)=denergy(istate)+ymat(jrep,irep,iab,istate)*tildem(irep,jrep,iab)
              dekin(istate)=dekin(istate)+ymat(jrep,irep,iab,istate)*tildemkin(irep,jrep,iab)
            enddo
          enddo
        enddo
        denergy(istate)=denergy(istate)*detiab(kref,1)*detiab(kref,2)/psid(istate)
        dekin(istate)=dekin(istate)*detiab(kref,1)*detiab(kref,2)/psid(istate)

        energy(istate)=denergy(istate)+eloc_det(kref,1)+eloc_det(kref,2)+e_other
        ekin(istate)=dekin(istate)+ekin_det(1)+ekin_det(2)+ekin_other

        if(ipr.ge.2) then
          write(ounit,'(''state'',i4)') istate
          write(ounit,'(''psid,psij'',9d12.5)') psid(istate),psij
          write(ounit,'(''psitot   '',e18.11)') psid(istate)*exp(psij)
c         do k=1,ndet
c           write(ounit,'(''psitot_k '',i6,3e18.8)') k, detiab(k,1),detiab(k,2),detiab(k,1)*detiab(k,2)*exp(psij)
c           write(ounit,'(''psitot_k '',i6,3e18.8)') k, detiab(k,1),detiab(k,2),cdet(k,1,1)*detiab(k,1)*detiab(k,2)*exp(psij)
c         enddo
c         do 25 i=1,nelec
c           do 25 k=1,3
c  25         write(ounit,'(''vj'',2e18.11)') vj(k,i)
          if(ipr.ge.3) write(ounit,'(''energy'',9f16.10)') energy(istate)
        endif

      enddo
      if(ifr.eq.1) then
        if(iforce_analy.eq.1) call compute_force(psid(1),denergy(1))

        call optorb_compute(psid,energy,denergy)

        call optjas_deloc(psid,energy,dvpsp_dj,vj)

        call optci_deloc(eloc_det,e_other,psid(1),energy(1))

        call prop_compute(coord)
      endif

      return
      end
      end module
