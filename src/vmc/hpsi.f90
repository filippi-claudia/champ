module hpsi_mod
contains
      subroutine hpsi(coord,psid,psij,ekin,energy,ipass,ifr)
!> @author Cyrus Umrigar
!> @author Claudia Filippi
!> @author Anthony Scemama
!> @author Claudio Amovilli
!> @author Franca Floris
!>
!> Written by Cyrus Umrigar, modified by Claudia Filippi and A. Scemama
!> modified by Claudio Amovilli and Franca Floris for PCM and QM-MMPOl

      use Bloc,    only: tildem,tildemkin,b,bkin
      use casula,  only: i_vpsp,t_vpsp
      use constants, only: hb
      use contrl_file, only: ounit
      use contrldmc, only: icut_e
      use control, only: ipr
      use csfs,    only: nstates
      use determinant_mod, only: compute_bmatrices_kin,determinant
      use determinant_psit_mod, only: determinant_psit
      use distance_mod, only: r_en,rvec_en
      use distances_mod, only: distances
      use efield,  only: iefield
      use efield_f_mod, only: efield_extpot_ene
      use force_analytic, only: compute_force
      use fragments, only: eloc_i, elocfrag, ifragelec, potnnfrag, nfrag
      use inputflags, only: iqmmm
      use jastrow_mod, only: jastrow_factor
      use m_force_analytic, only: iforce_analy
      use mmpol,   only: mmpol_extpot_ene
      use mmpol_cntrl, only: immpol
      use mmpol_hpsi, only: peQMdp,peQMq
      use mstates_mod, only: MSTATES
      use multidet, only: iactv,ivirt
      use multideterminant_mod, only: multideterminant_hpsi
      use multiple_geo, only: iwf,iwftype, nforce
      use multislater, only: detiab
      use nonloc_pot_mod, only: nonloc_pot
      use vmc_mod,        only: nwftypejas, nbjx, stoo, stoj, stobjx
      use optci_mod, only: optci_deloc
      use optjas_mod, only: optjas_deloc
      use optorb_f_mod, only: optorb_compute
      use orbval, only: nadorb
      use optwf_parms, only: nparmj
      use pcm_cntrl, only: ipcm
      use pcm_hpsi, only: pcms,pcmv
      use pcm_mod, only: pcm_extpot_ene
      use pot_local_mod, only: pot_local
      use precision_kinds, only: dp
      use properties_mod, only: prop_compute
      use pseudo,  only: nloc
      use qmmm_pot, only: qmmm_extpot_ene
      use slater,  only: kref,ndet,norb
      use system,  only: ndn,nelec,nup
      use velocity_jastrow, only: vj
      use ycompact, only: ymat

      implicit none

      integer :: i, iab, ifr, ii, ipass, j, o, x
      integer :: irep, istate, jrep, nel, iorb, i1, i2, iparm
      real(dp) :: e_other, ekin_other, ext_pot, peQM, pe_local, pepcm
      real(dp) :: tmp
      real(dp), dimension(3, *) :: coord
      real(dp), dimension(*) :: psid
      real(dp), dimension(*) :: psij
      real(dp), dimension(*) :: energy
      real(dp), dimension(nwftypejas) :: d2j
      real(dp), dimension(nstates) :: denergy
      real(dp), dimension(ndet, 2, nbjx) :: eloc_det
      real(dp), dimension(2, nbjx) :: vpsp_det, ekin_det
      real(dp), dimension(nparmj,nbjx) :: dvpsp_dj
      real(dp), dimension(*) :: ekin
      real(dp), dimension(MSTATES) :: dekin

      eloc_i = 0.d0 
      elocfrag = potnnfrag
! Calculates energy


!      iwf=iwftype(ifr)

! pe_ee computed in distances (called from distances if iperiodic != 0)
! pe_en(loc) computed in distances if nloc=0 or called from distances if iperiodic!=0
! pe_en(loc) computed in nonloc_pot if nloc!=0 and iperiodic=0
! pe_en(nonloc) computed in nonloc_pot if nloc !=0

! distances needed for Jastrow, determinants, and potential energy
      call distances(0,coord)
! local potential contributions
      call pot_local(coord,pe_local)

! external potential on a grid (e.g. MM from CPMD)
      if(iqmmm.eq.1) then
        ext_pot=0
        call qmmm_extpot_ene(coord,nelec,ext_pot)
        pe_local=pe_local+ext_pot
      endif

! external charges
      if(iefield.eq.1) then
        ext_pot=0
        call efield_extpot_ene(coord,nelec,ext_pot)
        pe_local=pe_local+ext_pot
      endif

! PCM polarization charges
      if(ipcm.gt.1) then
        pcms=0
        pcmv=0
        call pcm_extpot_ene(coord,nelec,pcms,pcmv)
        pepcm=pcms+pcmv
        pe_local=pe_local+pepcm
      endif

! QM-MMPOL (charges+induced dipoles)
      if(immpol.gt.1) then
        peQMdp=0
        peQMq=0
        call mmpol_extpot_ene(coord,nelec,peQMdp,peQMq)
        peQM=peQMdp+peQMq
        pe_local=pe_local+peQM
      endif

      if(ipr.ge.3) write(ounit,'(''pe_loc before nonloc_pot'',9f12.5)') pe_local

! get contribution from jastrow (also compute derivatives wrt parameters and nuclei)
      if(nforce.gt.1) iwf=iwftype(ifr)

      call jastrow_factor(coord,vj,d2j,psij,ifr)

! reset iwf to 1, why?
      if(nwftypejas.gt.1) iwf = 1

      if(ipr.ge.3) then
        do ii=1,nwftypejas
          write(ounit,'(''jastrow,d2j,psij'',i4,9f12.5)') ii, d2j(ii), psij(ii)
        enddo
      endif

! compute reference determinant, its derivatives, and kinetic contribution to B_eloc and its derivatives
      call determinant(ipass,coord,rvec_en,r_en)
      call compute_bmatrices_kin
! compute pseudo-potential contribution
! nonloc_pot must be called after determinant because slater matrices are needed

      if(nloc.gt.0) then
        call nonloc_pot(coord,rvec_en,r_en,pe_local,vpsp_det,dvpsp_dj,t_vpsp,i_vpsp,ifr)
      else
        do x=1,nbjx
          vpsp_det(1,x)=0.d0
          vpsp_det(2,x)=0.d0
          do iparm=1,nparmj
            dvpsp_dj(iparm,x)=0.d0
          enddo
        enddo
        do i=1,nelec
          do x=1,nbjx
            do iorb=1,norb+nadorb
              b(iorb,i,x)=bkin(iorb,i,x)
            enddo
          enddo
        enddo
      endif

      if(ipr.ge.3) then
        write(ounit,'(''pe_loc after nonloc_pot'',9f12.5)') pe_local
        do i=1,nstates
          write(ounit,'(''pe_ref after nonloc_pot, state ,up/down: '',i4,9f12.5)') i,(vpsp_det(ii,stobjx(i)),ii=1,2)
        enddo
      endif

      
      call multideterminant_hpsi(vj,ekin_det,vpsp_det,eloc_det)

      do istate=1,nstates
        j=stoj(istate)
        o=stoo(istate)
        x=stobjx(istate)

        ekin_other=-hb*d2j(j)
        do i=1,nelec
          tmp = -hb*(vj(1,i,j)**2+vj(2,i,j)**2+vj(3,i,j)**2)
          ekin_other=ekin_other+tmp
          if (icut_e.lt.0) then
            eloc_i(i)=eloc_i(i)+tmp
          end if
          if (nfrag.gt.1) then
            elocfrag(ifragelec(i)) = elocfrag(ifragelec(i)) + tmp
          endif
        enddo
        e_other=ekin_other+pe_local
! combine determinantal quantities to obtain trial wave function
        call determinant_psit(psid(istate),istate)
! compute energy using Ymat
        denergy(istate)=0
        dekin(istate)=0
        do iab=1,2
          nel=nup
          if(iab.eq.2) nel=ndn
          do jrep=ivirt(iab),norb
            do irep=iactv(iab),nel
              denergy(istate)=denergy(istate)+ymat(jrep,irep,iab,istate)*tildem(irep,jrep,iab,x)
              dekin(istate)=dekin(istate)+ymat(jrep,irep,iab,istate)*tildemkin(irep,jrep,iab,x)
            enddo
          enddo
        enddo
        denergy(istate)=denergy(istate)*detiab(kref,1,o)*detiab(kref,2,o)/psid(istate)
        dekin(istate)=dekin(istate)*detiab(kref,1,o)*detiab(kref,2,o)/psid(istate)

        energy(istate)=denergy(istate)+eloc_det(kref,1,x)+eloc_det(kref,2,x)+e_other
        ekin(istate)=dekin(istate)+ekin_det(1,x)+ekin_det(2,x)+ekin_other

        if(ipr.ge.2) then
          write(ounit,'(''state'',i4)') istate
          write(ounit,'(''psid,psij'',9d12.5)') psid(istate),psij(stoj(istate))
          write(ounit,'(''psitot   '',e18.11)') psid(istate)*exp(psij(stoj(istate)))
!         do k=1,ndet
!           write(ounit,'(''psitot_k '',i6,3e18.8)') k, detiab(k,1),detiab(k,2),detiab(k,1)*detiab(k,2)*exp(psij)
!           write(ounit,'(''psitot_k '',i6,3e18.8)') k, detiab(k,1),detiab(k,2),cdet(k,1,1)*detiab(k,1)*detiab(k,2)*exp(psij)
!         enddo
!         do 25 i=1,nelec
!           do 25 k=1,3
!  25         write(ounit,'(''vj'',2e18.11)') vj(k,i)
          if(ipr.ge.3) write(ounit,'(''energy'',9f16.10)') energy(istate)
        endif

      enddo
      if(ifr.eq.1) then
        if(iforce_analy.eq.1) call compute_force(psid(1),denergy(1))

        call optorb_compute(psid,energy,denergy)
        call optjas_deloc(psid,energy,dvpsp_dj,vj)
        call optci_deloc(eloc_det,e_other,psid,energy)
        call prop_compute(coord)
      endif

      return
      end
end module
