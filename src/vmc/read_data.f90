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
end
