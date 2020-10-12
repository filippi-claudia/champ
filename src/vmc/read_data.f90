subroutine read_znuc(iu)
!znuc inp
!KEYDOC nuclear charge for each atom type and ghost type

    use vmc_mod, only: MCTYPE
    use atom, only: znuc, nctype
    use ghostatom, only: newghostype
    use inputflags, only: iznuc

    implicit real*8(a - h, o - z)

    call p2gti('atoms:nctype', nctype, 1)
    call p2gtid('atoms:addghostype', newghostype, 0, 1)
    if (nctype + newghostype .gt. MCTYPE) call fatal_error('INPUT: nctype+newghostype > MCTYPE')

    allocate (znuc(nctype + newghostype))

    call incpos(iu, itmp, 1)
    read (iu, *) (znuc(i), i=1, nctype + newghostype)
    iznuc = 1
    call p2chkend(iu, 'znuc')
end

subroutine read_lcao(norb_tmp, nbasis_tmp, iwft, filename)
    ! lcao i i i=1 a=<input>
    !KEYDOC Orbital coefficients wrt complete basis.
    !KEYDOC Usage:  {        t lcao  norb,nbasis,filename,norbv}
    !KEYDOC norb: number of orbitals for trial wave function
    !KEYDOC nbasis: number of basis functiobns
    !KEYDOC iwft: wave function type (used when nforce>1 and wftype>1)
    !KEYDOC filename: file containing orbitals coefficients

    use vmc_mod, only: MORB, MBASIS
    use coefs, only: coef, nbasis, norb
    use inputflags, only: ilcao

    use pcm_fdc, only: fs

    ! was not in master but is needed
    use wfsec, only: nwftype

    implicit real*8(a - h, o - z)

    ! fs NOTE: additional variable norbv for efp orbitals removed

    character filename*(*)

    call file(iu, filename, 'old', 1, 0)
    nbasis = nbasis_tmp
    norb = norb_tmp
    nototal = norb
    if (nbasis .gt. MBASIS) call fatal_error('LCAO: nbasis > MBASIS')
    if (nototal .gt. MORB) call fatal_error('LCAO: number of orbitals > MORB')
    call p2gtid('general:nwftype', nwftype, 1, 1)

    if (iwft .gt. nwftype) call fatal_error('LCAO: wave function type > nwftype')

    allocate (coef(nbasis, norb, nwftype))

    do i = 1, nototal
        call incpos(iu, itmp, 1)
        read (iu, *) (coef(j, i, iwft), j=1, nbasis)
    enddo

    ilcao = ilcao + 1
    if (filename .eq. '<input>') then
        call p2chkend(iu, 'lcao')
    endif
end subroutine read_lcao

subroutine read_geometry(iu)
    ! geometry inp
    !KEYDOC position and type for each atom and ghost atom

    use vmc_mod, only: MCENT
    use atom, only: cent, iwctype, ncent
    use ghostatom, only: nghostcent
    use inputflags, only: igeometry

    implicit real*8(a - h, o - z)

    call p2gti('atoms:natom', ncent, 1)
    call p2gtid('atoms:nghostcent', nghostcent, 0, 1)
    if (ncent + nghostcent .gt. MCENT) call fatal_error('INPUT: ncent+nghostcent > MCENT')

    allocate (cent(3, ncent + nghostcent))

    do i = 1, ncent + nghostcent
        call incpos(iu, itmp, 1)
        read (iu, *) (cent(k, i), k=1, 3), iwctype(i)
    enddo
    igeometry = 1
    call p2chkend(iu, 'geometry')
end subroutine read_geometry

subroutine read_exponents(iu, iwft)
    ! exponents inp i=1
    !KEYDOC Basis function exponents (only if no numerical basis)

    use coefs, only: nbasis
    use basis, only: zex
    use inputflags, only: iexponents
    use wfsec, only: nwftype

    implicit real*8(a - h, o - z)

    call p2gtid('general:nwftype', nwftype, 1, 1)
    write (6, *) 'nbasis', nbasis

    allocate (zex(nbasis, nwftype))

    call incpos(iu, itmp, 1)
    read (iu, *) (zex(i, iwft), i=1, nbasis)
    iexponents = iexponents + 1
    call p2chkend(iu, 'exponents')
end

subroutine read_determinants(iu, nd, iwft)
! determinants inp i i=1
!KEYDOC CI coefficients and occupation of determinants in wf

    use vmc_mod, only: MELEC, MDET
    use dets, only: cdet, ndet
    use dorb_m, only: iworbd
    use inputflags, only: ideterminants
    use wfsec, only: nwftype
    use csfs, only: nstates

    ! not sure if needed but it's called
    use const, only: nelec
    implicit real*8(a - h, o - z)

    call p2gtid('general:nwftype', nwftype, 1, 1)

    ndet = nd

    if (ndet .gt. MDET) then
        write (6, *) "ndet=", ndet
        write (6, *) "MDET=", MDET
        call fatal_error('DET: ndet > MDET')
    endif

    call p2gti('electrons:nelec', nelec, 1)
    if (nelec .gt. MELEC) call fatal_error('INPUT: nelec exceeds MELEC')
    call incpos(iu, itmp, 1)

    allocate (cdet(ndet, nstates, nwftype))

    read (iu, *) (cdet(i, 1, iwft), i=1, ndet)

    allocate (iworbd(nelec, ndet))
    do i = 1, ndet
        call incpos(iu, itmp, 1)
        read (iu, *) (iworbd(j, i), j=1, nelec)
    enddo

    ideterminants = ideterminants + 1
    call p2chkend(iu, 'determinants')
    write (6, *) 'done det'
end

subroutine read_multideterminants(iu, nd)
! multideterminants inp i
!KEYDOC CI coefficients and occupation of determinants in wf
    use dets, only: ndet
    use const, only: nelec
    use multidet, only: irepcol_det, ireporb_det, numrep_det, iwundet
    use inputflags, only: imultideterminants

    implicit real*8(a - h, o - z)

    if (nd .ne. ndet - 1) call fatal_error('INPUT: problem in multidet')

    write (6, *) 'nelec', nelec
    write (6, *) 'ndet', ndet

    allocate (iwundet(ndet, 2))
    allocate (numrep_det(ndet, 2))
    allocate (irepcol_det(nelec, ndet, 2))
    allocate (ireporb_det(nelec, ndet, 2))

    call incpos(iu, itmp, 1)
    do k = 2, nd + 1
        read (iu, *) (numrep_det(k, iab), iab=1, 2)
        do iab = 1, 2
            do irep = 1, numrep_det(k, iab)
                read (iu, *) irepcol_det(irep, k, iab), ireporb_det(irep, k, iab)
            enddo
        enddo
    enddo

    imultideterminants = imultideterminants + 1
    call p2chkend(iu, 'multideterminants')
end subroutine read_multideterminants

subroutine read_jastrow_parameter(iu, iwft)
! jastrow_parameter inp i=1
!KEYDOC Parameters of Jastrow factor (depends on value of ijas!)

    use jaspar, only: nspin1, nspin2
    use elec, only: ndn
    use jaspar3, only: a, b, c, scalek
    use jaspar4, only: a4, norda, nordb, nordc
    use jaspar6, only: cutjas
    use bparm, only: nocuspb, nspin2b
    use contr2, only: ifock, ijas
    use contr2, only: isc
    use inputflags, only: ijastrow_parameter
    use wfsec, only: nwftype
    use atom, only: ncent, nctype

    implicit real*8(a - h, o - z)

    call p2gti('jastrow:ijas', ijas, 1)
    call p2gti('jastrow:isc', isc, 1)
    call p2gtid('jastrow:nspin1', nspin1, 1, 1)
    call p2gtid('jastrow:nspin2', nspin2, 1, 1)
    call p2gtid('jastrow:ifock', ifock, 0, 1)

    call p2gti('atoms:natom', ncent, 1)
    call p2gti('atoms:nctype', nctype, 1)

    call p2gtid('general:nwftype', nwftype, 1, 1)
    if (ijas .lt. 4 .or. ijas .gt. 6) call fatal_error('JASTROW: only ijas=4,5,6 implemented')
    if (ndn .eq. 1 .and. nspin2 .eq. 3) call fatal_error('JASTROW: 1 spin down and nspin2=3')

    if ((ijas .eq. 4 .or. ijas .eq. 5) .and. &
        (isc .ne. 2 .and. isc .ne. 4 .and. isc .ne. 6 .and. isc .ne. 7 .and. &
         isc .ne. 12 .and. isc .ne. 14 .and. isc .ne. 16 .and. isc .ne. 17)) &
        call fatal_error('JASTROW: if ijas=4 or 5, isc must be one of 2,4,6,7,12,14,16,17')

    if ((ijas .eq. 6) .and. (isc .ne. 6 .and. isc .ne. 7)) &
        call fatal_error('JASTROW: if ijas=6, isc must be 6 or 7')

    nspin2b = iabs(nspin2)
    nocuspb = 0
    if (nspin2 .lt. 0) then
        if (nspin2 .eq. -1) nocuspb = 1
        nspin2 = 1
    endif

    allocate (scalek(nwftype))

    if (ijas .ge. 4 .and. ijas .le. 6) then
        if (ifock .gt. 0) call fatal_error('JASTROW: fock not yet implemented for ijas=4,5,6')
        read (iu, *) norda, nordb, nordc
        if (isc .ge. 2) read (iu, *) scalek(iwft), a21
        mparmja = 2 + max(0, norda - 1)
        mparmjb = 2 + max(0, nordb - 1)
        mparmjc = nterms4(nordc)

        allocate (a4(mparmja, nctype, nwftype))

        do it = 1, nctype
            read (iu, *) (a4(iparm, it, iwft), iparm=1, mparmja)
            call incpos(iu, itmp, 1)
        enddo

        allocate (b(mparmjb, 2, nwftype))

        do isp = nspin1, nspin2b
            read (iu, *) (b(iparm, isp, iwft), iparm=1, mparmjb)
            call incpos(iu, itmp, 1)
        enddo

        allocate (c(mparmjc, nctype, nwftype))

        do it = 1, nctype
            read (iu, *) (c(iparm, it, iwft), iparm=1, mparmjc)
            call incpos(iu, itmp, 1)
        enddo

    endif
    !Read cutoff for Jastrow4, 5, 6
    if (isc .eq. 6 .or. isc .eq. 7) read (iu, *) cutjas

    ijastrow_parameter = ijastrow_parameter + 1
    call p2chkend(iu, 'jastrow_parameter')
    write (6, *) 'done jastrow'
end subroutine read_jastrow_parameter

subroutine read_bas_num_info(iu, numeric)
    ! basis inp i
    !KEYDOC Basis function types and pointers to radial parts tables
    ! qmc_bf_info inp i
    !KEYDOC alternative name for keyword basis because of GAMBLE input
    use numbas_mod, only: MRWF
    use vmc_mod, only: MCTYPE
    use numbas, only: iwrwf, numr
    use numbas1, only: iwlbas, nbastyp
    use basis, only: n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz
    use basis, only: n4s, n4p, n4fxxx, n4fyyy, n4fzzz, n4fxxy, n4fxxz, n4fyyx, n4fyyz
    use basis, only: n4fzzx, n4fzzy, n4fxyz, nsa, npa, ndzra, ndxya, ndxza, ndyza, ndx2a, ndz2a
    use inputflags, only: ibasis_num

    use atom, only: nctype
    use ghostatom, only: newghostype

    implicit real*8(a - h, o - z)

    call p2gti('atoms:nctype', nctype, 1)
    call p2gtid('atoms:addghostype', newghostype, 0, 1)
    if (nctype + newghostype .gt. MCTYPE) call fatal_error('ATOMS: nctype+newghostype > MCTYPE')

    nctot = nctype + newghostype
    allocate (nbastyp(nctot))
    allocate (n1s(nctot))
    allocate (n2s(nctot))
    allocate (n2p(3, nctot))
    allocate (n3s(nctot))
    allocate (n3p(3, nctot))
    allocate (n3dzr(nctot))
    allocate (n3dx2(nctot))
    allocate (n3dxy(nctot))
    allocate (n3dxz(nctot))
    allocate (n3dyz(nctot))
    allocate (n4s(nctot))
    allocate (n4p(3, nctot))
    allocate (n4fxxx(nctot))
    allocate (n4fyyy(nctot))
    allocate (n4fzzz(nctot))
    allocate (n4fxxy(nctot))
    allocate (n4fxxz(nctot))
    allocate (n4fyyx(nctot))
    allocate (n4fyyz(nctot))
    allocate (n4fzzx(nctot))
    allocate (n4fzzy(nctot))
    allocate (n4fxyz(nctot))
    allocate (nsa(nctot))
    allocate (npa(3, nctot))
    allocate (ndzra(nctot))
    allocate (ndz2a(nctot))
    allocate (ndxya(nctot))
    allocate (ndxza(nctot))
    allocate (ndx2a(nctot))
    allocate (ndyza(nctot))

    numr = numeric
    do i = 1, nctype + newghostype
        read (iu, *) n1s(i), n2s(i), (n2p(j, i), j=1, 3), &
            n3s(i), (n3p(j, i), j=1, 3), &
            n3dzr(i), n3dx2(i), n3dxy(i), n3dxz(i), n3dyz(i), &
            n4s(i), (n4p(j, i), j=1, 3), &
            n4fxxx(i), n4fyyy(i), n4fzzz(i), n4fxxy(i), n4fxxz(i), &
            n4fyyx(i), n4fyyz(i), n4fzzx(i), n4fzzy(i), n4fxyz(i), &
            nsa(i), (npa(j, i), j=1, 3), &
            ndzra(i), ndx2a(i), ndxya(i), ndxza(i), ndyza(i)
        call incpos(iu, itmp, 1)
        if (numr .gt. 0) then
            if (n2s(i) .ne. 0 .or. n3s(i) .ne. 0 .or. n4s(i) .ne. 0 .or. &
                n3p(1, i) .ne. 0 .or. n3p(2, i) .ne. 0 .or. n3p(3, i) .ne. 0 .or. &
                n4p(1, i) .ne. 0 .or. n4p(2, i) .ne. 0 .or. n4p(3, i) .ne. 0 .or. &
                nsa(i) .ne. 0 .or. npa(1, i) .ne. 0 .or. npa(2, i) .ne. 0 .or. &
                npa(3, i) .ne. 0 .or. ndzra(i) .ne. 0 .or. ndx2a(i) .ne. 0 .or. &
                ndxya(i) .ne. 0 .or. ndxza(i) .ne. 0 .or. ndyza(i) .ne. 0) &
                call fatal_error('BASIS: n1s,n2p,n3d only for numerical basis')

            nbastyp(i) = iabs(n1s(i)) &
                         + iabs(n2p(1, i)) + iabs(n2p(2, i)) + iabs(n2p(3, i)) &
                         + iabs(n3dzr(i)) + iabs(n3dx2(i)) &
                         + iabs(n3dxy(i)) + iabs(n3dxz(i)) + iabs(n3dyz(i)) &
                         + iabs(n4fxxx(i)) + iabs(n4fyyy(i)) + iabs(n4fzzz(i)) + iabs(n4fxxy(i)) + iabs(n4fxxz(i)) &
                         + iabs(n4fyyx(i)) + iabs(n4fyyz(i)) + iabs(n4fzzx(i)) + iabs(n4fzzy(i)) + iabs(n4fxyz(i))

            if (nbastyp(i) .gt. MRWF) call fatal_error('BASIS: nbastyp > MRWF')

            read (iu, *) (iwrwf(ib, i), ib=1, nbastyp(i))
            call incpos(iu, itmp, 1)
        else
            if (n4fxxx(i) .ne. 0 .or. n4fyyy(i) .ne. 0 .or. n4fzzz(i) .ne. 0 .or. &
                n4fxxy(i) .ne. 0 .or. n4fxxz(i) .ne. 0 .or. n4fyyx(i) .ne. 0 .or. &
                n4fyyz(i) .ne. 0 .or. n4fzzx(i) .ne. 0 .or. n4fzzy(i) .ne. 0 .or. &
                n4fxyz(i) .ne. 0) call fatal_error('BASIS: n4f only for numerical basis')
        endif
    enddo

    if (numr .gt. 0) then

        do i = 1, nctype + newghostype
            jj = 0
            do j = 1, iabs(n1s(i))
                jj = jj + 1
                iwlbas(jj, i) = 1
            enddo
            do j = 1, iabs(n2p(1, i))
                jj = jj + 1
                iwlbas(jj, i) = 2
            enddo
            do j = 1, iabs(n2p(2, i))
                jj = jj + 1
                iwlbas(jj, i) = 3
            enddo
            do j = 1, iabs(n2p(3, i))
                jj = jj + 1
                iwlbas(jj, i) = 4
            enddo
            do j = 1, iabs(n3dzr(i))
                jj = jj + 1
                iwlbas(jj, i) = 5
            enddo
            do j = 1, iabs(n3dx2(i))
                jj = jj + 1
                iwlbas(jj, i) = 6
            enddo
            do j = 1, iabs(n3dxy(i))
                jj = jj + 1
                iwlbas(jj, i) = 7
            enddo
            do j = 1, iabs(n3dxz(i))
                jj = jj + 1
                iwlbas(jj, i) = 8
            enddo
            do j = 1, iabs(n3dyz(i))
                jj = jj + 1
                iwlbas(jj, i) = 9
            enddo
            do j = 1, iabs(n4fxxx(i))
                jj = jj + 1
                iwlbas(jj, i) = 10
            enddo
            do j = 1, iabs(n4fyyy(i))
                jj = jj + 1
                iwlbas(jj, i) = 11
            enddo
            do j = 1, iabs(n4fzzz(i))
                jj = jj + 1
                iwlbas(jj, i) = 12
            enddo
            do j = 1, iabs(n4fxxy(i))
                jj = jj + 1
                iwlbas(jj, i) = 13
            enddo
            do j = 1, iabs(n4fxxz(i))
                jj = jj + 1
                iwlbas(jj, i) = 14
            enddo
            do j = 1, iabs(n4fyyx(i))
                jj = jj + 1
                iwlbas(jj, i) = 15
            enddo
            do j = 1, iabs(n4fyyz(i))
                jj = jj + 1
                iwlbas(jj, i) = 16
            enddo
            do j = 1, iabs(n4fzzx(i))
                jj = jj + 1
                iwlbas(jj, i) = 17
            enddo
            do j = 1, iabs(n4fzzy(i))
                jj = jj + 1
                iwlbas(jj, i) = 18
            enddo
            do j = 1, iabs(n4fxyz(i))
                jj = jj + 1
                iwlbas(jj, i) = 19
            enddo

        enddo
    endif
    ibasis_num = 1
    call p2chkend(iu, 'basis')
end subroutine read_bas_num_info

subroutine read_lattice(iu)
    ! lattice inp
    !KEYDOC Lattice vectors of primitive and simulation cell
    implicit real*8(a - h, o - z)
    call do_read_lattice(iu)
end subroutine read_lattice

subroutine read_forces(iu)
    ! forces_displace inp
    !KEYDOC Displacement parameters and wave function types

    use force_mod, only: MFORCE
    use vmc_mod, only: MCENT
    use forcepar, only: nforce
    use forcestr, only: delc
    use wfsec, only: iwftype
    use inputflags, only: iforces

    use atom, only: ncent

    implicit real*8(a - h, o - z)

    call p2gti('atoms:natom', ncent, 1)
    if (ncent .gt. MCENT) call fatal_error('FORCES: ncent > MCENT')

    call p2gtid('general:nforce', nforce, 1, 1)
    if (nforce .gt. MFORCE) call fatal_error('FORCES: nforce > MFORCE')

    allocate (delc(3, MCENT, MFORCE))
    allocate (iwftype(nforce))

    do i = 1, nforce
        do ic = 1, ncent
            call incpos(iu, itmp, 1)
            read (iu, *) (delc(k, ic, i), k=1, 3)
        enddo
    enddo

    call incpos(iu, itmp, 1)
    read (iu, *) (iwftype(i), i=1, nforce)
    if (iwftype(1) .ne. 1) call fatal_error('INPUT: iwftype(1) ne 1')

    iforces = 1
    call p2chkend(iu, 'forces')

end subroutine read_forces

subroutine read_csf(ncsf_read, nstates_read, fn)
    ! csf i i=1 a=<input>

    use vmc_mod, only: MDET
    use csfs, only: ccsf, ncsf, nstates
    use mstates_mod, only: MSTATES
    use inputflags, only: icsfs
    use wfsec, only: nwftype
    implicit real*8(a - h, o - z)

    character fn*(*)
    call p2gtid('general:nwftype', nwftype, 1, 1)
    call ptfile(iu, fn, 'old')

    ncsf = ncsf_read
    if (ncsf .gt. MDET) call fatal_error('CSF: too many csf')

    nstates = nstates_read
    if (nstates .gt. MSTATES) call fatal_error('CSF: too many states')

    allocate (ccsf(ndet, nstates, nwftype))

    do i = 1, nstates
        read (iu, *) (ccsf(j, i, 1), j=1, ncsf)
    enddo

    icsfs = 1

    if (fn .eq. '<input>') then
        call p2chkend(iu, 'csf')
    endif

end subroutine read_csf

subroutine read_csfmap(fn)
    ! csfmap a=<input>
    !KEYDOC Read mapping between csf and determinants.
    use vmc_mod, only: MDET
    use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates
    use mstates_mod, only: MDETCSFX
    use dets, only: cdet, ndet
    use wfsec, only: nwftype

    implicit real*8(a - h, o - z)

    character fn*(*)

    call ptfile(iu, fn, 'old')
    call p2gtid('general:nwftype', nwftype, 1, 1)

    read (iu, *) ncsf_check, ndet_check, nmap_check
    write (6, '(''csfmap'',3i4)') ncsf_check, ndet_check, nmap_check
    if (ndet_check .ne. ndet) call fatal_error('CSFMAP: wrong number of determinants')
    if (ncsf_check .ne. ncsf) call fatal_error('CSFMAP: wrong number of csf')
    if (nmap_check .gt. float(MDET)*MDET) call fatal_error('CSFMAP: too many determinants in map list')

    nptr = 1
    do i = 1, ncsf
        read (iu, *) nterm
        iadet(i) = nptr
        ibdet(i) = nptr + nterm - 1
        do j = 1, nterm
            read (iu, *) id, c
            icxdet(nptr) = id
            cxdet(nptr) = c
            nptr = nptr + 1
            if (nptr .gt. MDET*MDETCSFX) call fatal_error('CSFMAP: problem with nmap')
        enddo
    enddo
    if (nmap_check .ne. nptr - 1) call fatal_error('CSFMAP: problem with nmap')
    nmap = nptr

    if (.not. allocated(cdet)) allocate (cdet(ndet, nstates, nwftype))

    write (6, '(''Warning: det coef overwritten with csf'')')
    do k = 1, nstates
        do j = 1, ndet
            cdet(j, k, 1) = 0
        enddo
        do icsf = 1, ncsf
            do j = iadet(icsf), ibdet(icsf)
                jx = icxdet(j)
                cdet(jx, k, 1) = cdet(jx, k, 1) + ccsf(icsf, k, 1)*cxdet(j)
            enddo
        enddo
    enddo

    if (fn .eq. '<input>') then
        call p2chkend(iu, 'csfmap')
    endif

    return
end subroutine read_csfmap

subroutine read_jasderiv(iu)
!    jasderiv inp
    use optjas, only: MPARMJ
    use atom, only: nctype
    use jaspar, only: nspin1, is
    use jaspar4, only: norda, nordb, nordc
    use jaspointer, only: npoint, npointa
    use numbas, only: numr

    use optwf_nparmj, only: nparma, nparmb, nparmc, nparmf
    use optwf_parms, only: nparmj
    use optwf_wjas, only: iwjasa, iwjasb, iwjasc
    use bparm, only: nspin2b
    use contr2, only: ijas
    use contr2, only: isc
    implicit real*8(a - h, o - z)

    na1 = 1
    na2 = nctype

    if (.not. allocated(nparma)) allocate (nparma(MCTYP3X))
    if (.not. allocated(nparmb)) allocate (nparmb(3))
    if (.not. allocated(nparmc)) allocate (nparmc(nctype))
    if (.not. allocated(nparmf)) allocate (nparmf(nctype))

    read (iu, *) (nparma(ia), ia=na1, na2), &
        (nparmb(isp), isp=nspin1, nspin2b), &
        (nparmc(it), it=1, nctype), &
        (nparmf(it), it=1, nctype)

    if (ijas .ge. 4 .and. ijas .le. 6) then
        do it = 1, nctype
            if (numr .eq. 0) then
                ! All-electron with analytic slater basis
                if ((nparma(it) .gt. 0 .and. norda .eq. 0) .or. (nparma(it) .gt. norda + 1)) then
                    write (6, '(''it,norda,nparma(it)'',3i5)') it, norda, nparma(it)
                    stop 'nparma too large for norda'
                endif
            else
                ! Pseudopotential with numerical basis: cannot vary a(1) or a(2)
                if (norda .eq. 1) stop 'makes no sense to have norda=1 for numr>0'
                if ((norda .eq. 0 .and. nparma(it) .gt. 0) .or. (norda .gt. 0 .and. nparma(it) .gt. norda - 1)) then
                    write (6, '(''it,norda,nparma(it)'',3i5)') it, norda, nparma(it)
                    stop 'nparma too large for norda'
                endif
            endif

            if (isc .le. 7 .and. &
                ((nordc .le. 2 .and. nparmc(it) .gt. 0) &
                 .or. (nordc .eq. 3 .and. nparmc(it) .gt. 2) &
                 .or. (nordc .eq. 4 .and. nparmc(it) .gt. 7) &
                 .or. (nordc .eq. 5 .and. nparmc(it) .gt. 15) &
                 .or. (nordc .eq. 6 .and. nparmc(it) .gt. 27) &
                 .or. (nordc .eq. 7 .and. nparmc(it) .gt. 43))) then
                write (6, '(''it,nordc,nparmc(it)'',3i5)') it, nordc, nparmc(it)
                stop 'nparmc too large for nordc in J_een with cusp conds'
            endif

            if (isc .gt. 7 .and. &
                ((nordc .le. 1 .and. nparmc(it) .gt. 0) &
                 .or. (nordc .eq. 2 .and. nparmc(it) .gt. 2) &
                 .or. (nordc .eq. 3 .and. nparmc(it) .gt. 6) &
                 .or. (nordc .eq. 4 .and. nparmc(it) .gt. 13) &
                 .or. (nordc .eq. 5 .and. nparmc(it) .gt. 23) &
                 .or. (nordc .eq. 6 .and. nparmc(it) .gt. 37) &
                 .or. (nordc .eq. 7 .and. nparmc(it) .gt. 55))) then
                write (6, '(''it,nordc,nparmc(it)'',3i5)') it, nordc, nparmc(it)
                stop 'nparmc too large for nordc without cusp conds'
            endif

        enddo

        ! For the b coefs. we assume that b(1) is fixed by the cusp-cond.
        do isp = 1, nspin1, nspin2b
            if (nparmb(isp) .gt. nordb) then
                write (6, '(''isp,nordb,nparmb(isp)'',3i5)') isp, nordb, nparmb(isp)
                stop 'nparmb too large for nordb'
            endif
        enddo
    endif

    ! compute nparmj
    nparmj = 0
    npointa(1) = 0
    do ia = na1, na2
        if (ia .gt. 1) npointa(ia) = npointa(ia - 1) + nparma(ia - 1)
        nparmj = nparmj + nparma(ia)
    enddo
    do isp = nspin1, nspin2b
        nparmj = nparmj + nparmb(isp)
    enddo
    npoint(1) = nparmj
    do it = 1, nctype
        if (it .gt. 1) npoint(it) = npoint(it - 1) + nparmc(it - 1)
        nparmj = nparmj + nparmc(it) + nparmf(it)
    enddo

    if (nparmj .gt. MPARMJ) call fatal_error('JASDERIV: MPARMJ too small')

    do it = 1, nctype
        read (iu, *) (iwjasa(iparm, it), iparm=1, nparma(it))
    enddo
    do isp = nspin1, nspin2b
        read (iu, *) (iwjasb(iparm, isp), iparm=1, nparmb(isp))
    enddo
    do it = 1, nctype
        read (iu, *) (iwjasc(iparm, it), iparm=1, nparmc(it))
    enddo

    call p2chkend(iu, 'jasderiv')
end subroutine read_jasderiv

subroutine read_sym(nsym, mo, fn)
    ! sym_labels i i a=<input>
    ! KEYDOC Read symmetry information
    use coefs, only: norb
    use optorb, only: irrep
    implicit real*8(a - h, o - z)

    character fn*(*)
    character atmp*80

    call ptfile(iu, fn, 'old')
    nirrep = nsym
    if (norb .ne. 0 .and. norb .ne. mo) then
        write (6, '(2i5)') norb, mo
        call fatal_error('READSYM: wrong number of orbitals')
    else
        norb = mo
    endif

    ! Ignore irrep text labels
    read (iu, '(a80)') atmp

    ! allocate
    allocate (irrep(norb))

    ! read data
    read (iu, *) (irrep(io), io=1, norb)

    if (fn .eq. '<input>') then
        call p2chkend(iu, 'sym_labels')
    endif
end subroutine read_sym