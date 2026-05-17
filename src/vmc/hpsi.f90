module hpsi_mod
contains

  subroutine hpsi(coord, psid, psij, ekin, energy, ipass, ifr)
    !> @author Cyrus Umrigar
    !> @author Claudia Filippi
    !> @author Anthony Scemama
    !> @author Claudio Amovilli
    !> @author Franca Floris
    !>
    !> Written by Cyrus Umrigar, modified by Claudia Filippi and A. Scemama.
    !> Modified by Claudio Amovilli and Franca Floris for PCM and QM-MMPOL.

    use Bloc, only: b, bkin, tildem, tildemkin, xmat, xmatkin
    use bxmatrices, only: bxmatrix
    use casula, only: i_vpsp, t_vpsp
    use constants, only: hb
    use contrl_file, only: ounit
    use contrldmc, only: icut_e
    use control, only: ipr
    use csfs, only: nstates
    use determinant_mod, only: compute_bmatrices_kin, determinant
    use determinant_psit_mod, only: determinant_psit
    use distance_mod, only: r_en, rvec_en
    use distances_mod, only: distances
    use force_analytic, only: compute_force
    use fragments, only: eloc_i, elocfrag, ifragelec, nfrag, potnnfrag
    use jastrow_mod, only: jastrow_factor
    use m_force_analytic, only: iforce_analy
    use mstates_mod, only: MSTATES
    use multidet, only: iactv, ivirt
    use multideterminant_mod, only: multideterminant_hpsi
    use multiple_geo, only: iwf, iwftype, nforce
    use multislater, only: detiab
    use nonloc_pot_mod, only: nonloc_pot
    use optci_mod, only: optci_deloc
    use optjas_mod, only: optjas_deloc
    use optorb_f_mod, only: optorb_compute
    use optwf_control, only: ioptorb
    use optwf_parms, only: nparmj
    use orbval, only: nadorb
    use pot_local_mod, only: pot_local
    use precision_kinds, only: dp
    use properties_mod, only: prop_compute
    use pseudo, only: nloc
    use slater, only: d2dx2, ddx, kref, ndet, norb
    use system, only: ndn, nelec, nup
    use velocity_jastrow, only: vj
    use vmc_mod, only: nbjx, nwftypejas, stoj, stobjx, stoo
    use ycompact, only: ymat

    implicit none

    integer :: i, iab, ifr, ii, iparm, ipass, iorb, irep, ish
    integer :: istate, j, jrep, nel, o, x
    real(dp) :: e_other, ekin_other, pe_local
    real(dp) :: tmp
    real(dp), dimension(3, *) :: coord
    real(dp), dimension(*) :: psid
    real(dp), dimension(*) :: psij
    real(dp), dimension(*) :: energy
    real(dp), dimension(*) :: ekin
    real(dp), dimension(nwftypejas) :: d2j
    real(dp), dimension(nstates) :: denergy
    real(dp), dimension(ndet, 2, nbjx) :: eloc_det
    real(dp), dimension(2, nbjx) :: vpsp_det, ekin_det
    real(dp), dimension(nparmj, nbjx) :: dvpsp_dj
    real(dp), dimension(MSTATES) :: dekin

    eloc_i = 0.d0
    elocfrag = potnnfrag

    ! pe_ee is computed in distances.
    ! pe_en(loc) is computed in distances when nloc=0 or iperiodic/=0.
    ! pe_en(loc) is computed in nonloc_pot when nloc/=0 and iperiodic=0.
    ! pe_en(nonloc) is computed in nonloc_pot when nloc/=0.

    ! Distances needed for Jastrow, determinants, and potential energy.
    call distances(0, coord)

    ! Local potential contributions.
    call pot_local(coord, pe_local)
    call add_external_potential_contributions(coord, pe_local)

    if (ipr .ge. 3) write(ounit,'(''pe_loc before nonloc_pot'',9f12.5)') pe_local

    ! Get contribution from Jastrow and derivatives with respect to parameters and nuclei.
    if (nforce .gt. 1) iwf = iwftype(ifr)

    call jastrow_factor(coord, vj, d2j, psij, ifr)

    ! TODO: clarify why multi-Jastrow calculations reset iwf here.
    if (nwftypejas .gt. 1) iwf = 1

    if (ipr .ge. 3) then
      do ii = 1, nwftypejas
        write(ounit,'(''jastrow,d2j,psij'',i4,9f12.5)') ii, d2j(ii), psij(ii)
      enddo
    endif

    ! Compute reference determinant, its derivatives, and kinetic contribution.
    call determinant(ipass, coord, rvec_en, r_en)
    call compute_bmatrices_kin

    ! Compute pseudo-potential contribution.
    ! nonloc_pot must be called after determinant because slater matrices are needed.
    if (nloc .gt. 0) then
      call nonloc_pot(coord, rvec_en, r_en, pe_local, vpsp_det, dvpsp_dj, &
          t_vpsp, i_vpsp, ifr)
    else
      do x = 1, nbjx
        vpsp_det(1, x) = 0.d0
        vpsp_det(2, x) = 0.d0
        do iparm = 1, nparmj
          dvpsp_dj(iparm, x) = 0.d0
        enddo
      enddo
      do i = 1, nelec
        do x = 1, nbjx
          do iorb = 1, norb + nadorb
            b(iorb, i, x) = bkin(iorb, i, x)
          enddo
        enddo
      enddo
    endif

    if (ipr .ge. 3) then
      write(ounit,'(''pe_loc after nonloc_pot'',9f12.5)') pe_local
      do i = 1, nstates
        write(ounit,'(''pe_ref after nonloc_pot, state ,up/down: '',i4,9f12.5)') &
            i, (vpsp_det(ii, stobjx(i)), ii = 1, 2)
      enddo
    endif

    ! Single determinant local energy calculations.
    do istate = 1, nstates
      j = stoj(istate)
      o = stoo(istate)
      x = stobjx(istate)
      nel = nup
      ish = 0
      do iab = 1, 2
        if (iab .eq. 2) then
          nel = ndn
          ish = nup
        endif
        ekin_det(iab, x) = 0.d0
        do i = 1, nel
          tmp = -hb * (d2dx2(i + ish, o) + &
              2.d0 * (vj(1, i + ish, j) * ddx(1, i + ish, o) + &
                      vj(2, i + ish, j) * ddx(2, i + ish, o) + &
                      vj(3, i + ish, j) * ddx(3, i + ish, o)))
          if (icut_e .lt. 0) then
            eloc_i(i + ish) = eloc_i(i + ish) + tmp
          endif
          if (nfrag .gt. 1) then
            elocfrag(ifragelec(i + ish)) = elocfrag(ifragelec(i + ish)) + tmp
          endif
          ekin_det(iab, x) = ekin_det(iab, x) + tmp
        enddo
        eloc_det(kref, iab, x) = ekin_det(iab, x) + vpsp_det(iab, x)
      enddo

      if (ndet .ne. 1 .or. iforce_analy .ne. 0 .or. ioptorb .ne. 0) then
        call bxmatrix(kref, xmat(1, 1, x), xmat(1, 2, x), b(1, 1, x), x)
        call bxmatrix(kref, xmatkin(1, 1, x), xmatkin(1, 2, x), bkin(1, 1, x), x)
      endif
    enddo

    if (.not. (ndet .eq. 1 .and. ioptorb .eq. 0)) then
      call multideterminant_hpsi(vj, ekin_det, vpsp_det, eloc_det)
    endif

    do istate = 1, nstates
      j = stoj(istate)
      o = stoo(istate)
      x = stobjx(istate)

      ekin_other = -hb * d2j(j)
      do i = 1, nelec
        tmp = -hb * (vj(1, i, j)**2 + vj(2, i, j)**2 + vj(3, i, j)**2)
        ekin_other = ekin_other + tmp
        if (icut_e .lt. 0) then
          eloc_i(i) = eloc_i(i) + tmp
        endif
        if (nfrag .gt. 1) then
          elocfrag(ifragelec(i)) = elocfrag(ifragelec(i)) + tmp
        endif
      enddo
      e_other = ekin_other + pe_local

      ! Combine determinantal quantities to obtain the trial wave function.
      call determinant_psit(psid(istate), istate)

      ! Compute energy using Ymat.
      denergy(istate) = 0.d0
      dekin(istate) = 0.d0
      do iab = 1, 2
        nel = nup
        if (iab .eq. 2) nel = ndn
        do jrep = ivirt(iab), norb
          do irep = iactv(iab), nel
            denergy(istate) = denergy(istate) + &
                ymat(jrep, irep, iab, istate) * tildem(irep, jrep, iab, x)
            dekin(istate) = dekin(istate) + &
                ymat(jrep, irep, iab, istate) * tildemkin(irep, jrep, iab, x)
          enddo
        enddo
      enddo
      denergy(istate) = denergy(istate) * detiab(kref, 1, o) * &
          detiab(kref, 2, o) / psid(istate)
      dekin(istate) = dekin(istate) * detiab(kref, 1, o) * &
          detiab(kref, 2, o) / psid(istate)

      energy(istate) = denergy(istate) + eloc_det(kref, 1, x) + &
          eloc_det(kref, 2, x) + e_other
      ekin(istate) = dekin(istate) + ekin_det(1, x) + ekin_det(2, x) + ekin_other

      if (ipr .ge. 2) then
        write(ounit,'(''state'',i4)') istate
        write(ounit,'(''psid,psij'',9d12.5)') psid(istate), psij(stoj(istate))
        write(ounit,'(''psitot   '',e18.11)') psid(istate) * exp(psij(stoj(istate)))
!       do k=1,ndet
!         write(ounit,'(''psitot_k '',i6,3e18.8)') k, detiab(k,1),detiab(k,2),detiab(k,1)*detiab(k,2)*exp(psij)
!         write(ounit,'(''psitot_k '',i6,3e18.8)') k, detiab(k,1),detiab(k,2),cdet(k,1,1)*detiab(k,1)*detiab(k,2)*exp(psij)
!       enddo
!       do 25 i=1,nelec
!         do 25 k=1,3
!25       write(ounit,'(''vj'',2e18.11)') vj(k,i)
        if (ipr .ge. 3) write(ounit,'(''energy'',9f16.10)') energy(istate)
      endif

    enddo

    if (ifr .eq. 1) then
      if (iforce_analy .eq. 1) call compute_force(psid(1), denergy(1))

      call optorb_compute(psid, energy, denergy)
      call optjas_deloc(psid, energy, dvpsp_dj, vj)
      call optci_deloc(eloc_det, e_other, psid, energy)
      call prop_compute(coord)
    endif

  end subroutine hpsi

  subroutine add_external_potential_contributions(coord, pe_local)
    use efield, only: iefield
    use efield_f_mod, only: efield_extpot_ene
    use inputflags, only: iqmmm
    use mmpol, only: mmpol_extpot_ene
    use mmpol_cntrl, only: immpol
    use mmpol_hpsi, only: peQMdp, peQMq
    use pcm_cntrl, only: ipcm
    use pcm_hpsi, only: pcms, pcmv
    use pcm_mod, only: pcm_extpot_ene
    use precision_kinds, only: dp
    use qmmm_pot, only: qmmm_extpot_ene
    use system, only: nelec

    implicit none

    real(dp), dimension(3, *) :: coord
    real(dp), intent(inout) :: pe_local
    real(dp) :: ext_pot, pepcm, peQM

    ! External potential on a grid, e.g. MM from CPMD.
    if (iqmmm .eq. 1) then
      ext_pot = 0.d0
      call qmmm_extpot_ene(coord, nelec, ext_pot)
      pe_local = pe_local + ext_pot
    endif

    ! External charges.
    if (iefield .eq. 1) then
      ext_pot = 0.d0
      call efield_extpot_ene(coord, nelec, ext_pot)
      pe_local = pe_local + ext_pot
    endif

    ! PCM polarization charges. pcms and pcmv are accumulated later by pcm_vmc.
    if (ipcm .gt. 1) then
      pcms = 0.d0
      pcmv = 0.d0
      call pcm_extpot_ene(coord, nelec, pcms, pcmv)
      pepcm = pcms + pcmv
      pe_local = pe_local + pepcm
    endif

    ! QM-MMPOL charges and induced dipoles. peQMdp and peQMq are accumulated later by mmpol_vmc.
    if (immpol .gt. 1) then
      peQMdp = 0.d0
      peQMq = 0.d0
      call mmpol_extpot_ene(coord, nelec, peQMdp, peQMq)
      peQM = peQMdp + peQMq
      pe_local = pe_local + peQM
    endif

  end subroutine add_external_potential_contributions

end module hpsi_mod
