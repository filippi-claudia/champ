subroutine inputzex
    ! Set the exponents to one when using a numerical basis
    use force_mod, only: MWF
    use numbas, only: numr
    use coefs, only: nbasis
    use basis, only: zex

    ! are they needed ??!!
    use contrl_per, only: iperiodic
    use wfsec, only: nwftype
    use method_opt, only: method
      implicit none

      integer :: i, iwft



    if( (method(1:3) == 'lin')) then
        if (.not. allocated(zex)) allocate (zex(nbasis, 3))
    else
        if (.not. allocated(zex)) allocate (zex(nbasis, nwftype))
    endif

    if (numr .eq. 0 .and. iperiodic .eq. 0) &
        call fatal_error('ZEX: numr=0 and iperiodic=0 but no zex are inputed')

    zex = 1
    return
end subroutine inputzex

subroutine inputcsf
    ! Check that the required blocks are there in the input


    use csfs, only: ncsf, nstates
    use inputflags, only: ici_def
    use ci000, only: nciprim, nciterm

    ! are they needed ??!!
    use optwf_contrl, only: ioptci
    implicit none


    nstates = 1
    ncsf = 0

    if (ioptci .ne. 0 .and. ici_def .eq. 1) nciterm = nciprim
    return
end subroutine inputcsf

subroutine multideterminants_define(iflag, icheck)

    use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
    use vmc_mod, only: MORB, MCENT, MCTYPE, MCTYP3X
    use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
    use vmc_mod, only: radmax, delri
    use vmc_mod, only: NEQSX, MTERMS
    use vmc_mod, only: MCENT3, NCOEF, MEXCIT
    use const, only: nelec
    use csfs, only: cxdet, iadet, ibdet, icxdet, ncsf, nstates
    use dets, only: cdet, ndet
    use elec, only: ndn, nup
    use multidet, only: iactv, irepcol_det, ireporb_det, ivirt, iwundet, kref, numrep_det, allocate_multidet
    use coefs, only: norb
    use dorb_m, only: iworbd

    use contrl_file,    	only: ounit, errunit

    ! not sure about that one either ....
    use wfsec, only: nwftype

      use precision_kinds, only: dp
      implicit none


    interface
    function idiff(j, i, iab)
        integer, intent(in) :: j
        integer, intent(in) :: i
        integer, intent(in) :: iab
        integer :: idiff
    endfunction
    end interface

      integer :: i, iab, icheck, icsf, idist
      integer :: iflag, in, iphase, iref
      integer :: irep, isav, ish, istate
      integer :: isub, iw, iwf, iwref
      integer :: j, k, kref_old, l
      integer :: ndet_dist, nel
      integer, dimension(nelec) :: iswapped
      integer, dimension(ndet) :: itotphase



    save kref_old

    if (nup .gt. nelec/2) call fatal_error('INPUT: nup exceeds nelec/2')
    ndn = nelec - nup

    !!call p2gtid('general:nwftype', nwftype, 1, 1)
    if (nwftype .gt. MWF) call fatal_error('INPUT: nwftype exceeds MWF')

    if (iflag .eq. 0) then
        kref = 1
    else
        if (kref .gt. 1 .and. icheck .eq. 1) then
            kref = 1
            goto 2
        endif
1       kref = kref + 1
        if (kref .gt. ndet) call fatal_error('MULTIDET_DEFINE: kref > ndet')

2       if (idiff(kref_old, kref, iflag) .eq. 0) goto 1
        write (ounit, *) 'kref change', iflag, kref_old, kref
    endif
    kref_old = kref

    if (.not. allocated(iwundet)) allocate (iwundet(ndet, 2))
    if (.not. allocated(numrep_det)) allocate (numrep_det(ndet, 2))
    if (.not. allocated(irepcol_det)) allocate (irepcol_det(nelec, ndet, 2))
    if (.not. allocated(ireporb_det)) allocate (ireporb_det(nelec, ndet, 2))

    do iab = 1, 2
        numrep_det(kref, iab) = 0
    enddo

    do k = 1, ndet
        itotphase(k) = 0
        if (k .eq. kref) goto 5
        do iab = 1, 2
            nel = nup
            ish = 0
            if (iab .eq. 2) then
                nel = ndn
                ish = nup
            endif
            numrep_det(k, iab) = 0
            do iref = 1, nel
                iwref = iworbd(iref + ish, kref)
                in = 0
                do i = 1, nel
                    iw = iworbd(i + ish, k)
                    if (iw .eq. iwref) in = 1
                enddo
                if (in .eq. 0) then
                    numrep_det(k, iab) = numrep_det(k, iab) + 1
                    irepcol_det(numrep_det(k, iab), k, iab) = iref
                endif
            enddo
            isub = 0
            do i = 1, nel
                iw = iworbd(i + ish, k)
                in = 0
                do iref = 1, nel
                    iwref = iworbd(iref + ish, kref)
                    if (iw .eq. iwref) in = 1
                enddo
                if (in .eq. 0) then
                    isub = isub + 1
                    ireporb_det(isub, k, iab) = iw
                endif
            enddo
            if (isub .ne. numrep_det(k, iab)) then
                write (ounit, *) isub, numrep_det(k, iab)
                stop 'silly error'
            endif
            do irep = 1, nel
                iswapped(irep) = iworbd(irep + ish, kref)
            enddo
            do irep = 1, numrep_det(k, iab)
                iswapped(irepcol_det(irep, k, iab)) = ireporb_det(irep, k, iab)
            enddo
            iphase = 0
            do i = 1, nel
                if (iworbd(i + ish, k) .ne. iswapped(i)) then
                    do l = i + 1, nel
                        if (iswapped(l) .eq. iworbd(i + ish, k)) then
                            isav = iswapped(i)
                            iswapped(i) = iswapped(l)
                            iswapped(l) = isav
                            iphase = iphase + 1
                        endif
                    enddo
                endif
            enddo

            itotphase(k) = itotphase(k) + iphase
        enddo
        do iwf = 1, nwftype
            do istate = 1, nstates
                cdet(k, istate, iwf) = cdet(k, istate, iwf)*(-1)**itotphase(k)
            enddo
        enddo
5       continue
    enddo

    do k = 1, ndet
        if (k .eq. kref) goto 6
        do i = 1, nelec
            iworbd(i, k) = iworbd(i, kref)
        enddo
        do iab = 1, 2
            ish = 0
            if (iab .eq. 2) ish = nup
            do irep = 1, numrep_det(k, iab)
                iworbd(irepcol_det(irep, k, iab) + ish, k) = ireporb_det(irep, k, iab)
            enddo
        enddo
6       continue
    enddo

    call allocate_multidet()

    iactv(1) = nup + 1
    iactv(2) = ndn + 1
    ivirt(1) = nup + 1
    ivirt(2) = ndn + 1
    do k = 1, ndet
        if (k .eq. kref) go to 8
        do iab = 1, 2
            do irep = 1, numrep_det(k, iab)
        if (irepcol_det(irep, k, iab) .ne. 0 .and. irepcol_det(irep, k, iab) .lt. iactv(iab)) iactv(iab) = irepcol_det(irep, k, iab)
                if (ireporb_det(irep, k, iab) .lt. ivirt(iab)) ivirt(iab) = ireporb_det(irep, k, iab)
            enddo
        enddo
8       continue
    enddo

    write (ounit, *) ' Multideterminants :: '
    write (ounit, *) 'norb  =', norb
    write (ounit, *) 'iactv =', (iactv(iab), iab=1, 2)
    write (ounit, *) 'ivirt =', (ivirt(iab), iab=1, 2)

    idist = 1
    if (idist .eq. 0) then
        do iab = 1, 2
            do i = 1, ndet
                iwundet(i, iab) = i
            enddo
        enddo
    else
        do iab = 1, 2
            do i = 1, ndet
                iwundet(i, iab) = i
                if (i .eq. kref) goto 10
                if (idiff(kref, i, iab) .eq. 0) then
                    iwundet(i, iab) = kref
                    goto 10
                endif
                do j = 1, i - 1
                    if (idiff(j, i, iab) .eq. 0) then
                        iwundet(i, iab) = j
                        go to 10
                    endif
                enddo
10              continue
            enddo
        enddo
        do iab = 1, 2
            ndet_dist = 0
            do i = 1, ndet
                if (iwundet(i, iab) .eq. i) then
                    ndet_dist = ndet_dist + 1
                endif
            enddo
            write (ounit, *) iab, ndet_dist, ' distinct out of ', ndet
        enddo
    endif

    do icsf = 1, ncsf
        do j = iadet(icsf), ibdet(icsf)
            k = icxdet(j)
            cxdet(j) = cxdet(j)*(-1)**itotphase(k)
        enddo
    enddo

    return
end subroutine multideterminants_define

subroutine inputforces
! Set all force displacements to zero
!    use force_mod, only: MWF
!    use vmc_mod, only: MCENT
!    use force_mod, only: MFORCE
    use forcepar, only: nforce
    use forcestr, only: delc
    use wfsec, only: iwftype, nwftype
    use contrl_file, only: errunit
    use atom, only: ncent

    implicit none
    integer             :: i

    if (.not. allocated(delc)) allocate (delc(3, ncent, nforce))
    if (.not. allocated(iwftype)) allocate (iwftype(nforce))

    call set_displace_zero(nforce)


    if (nwftype .eq. 1) then
        iwftype = 1             ! set array to 1 for all nforce elements
    elseif (nwftype .eq. nforce) then
        do i = 1, nforce
            iwftype(i) = i
        enddo
    else
        call fatal_error('FORCES: need to specify iwftype')
    endif

end subroutine inputforces

subroutine inputdet()
    ! Set the cdet to be equal
    use dets, only: cdet, ndet
    use csfs, only: nstates
!    use vmc_mod, only: MORB
!    use mstates_mod, only: MSTATES
    use wfsec, only: nwftype
    use method_opt, only: method

    implicit none
    integer             :: iwft, k

    if( (method(1:3) == 'lin')) then
        if (.not. allocated(cdet)) allocate(cdet(ndet,nstates,3))
    else
        if (.not. allocated(cdet)) allocate(cdet(ndet,nstates,nwftype))
    endif

    do iwft = 2, nwftype
        do k = 1, ndet
            cdet(k, 1, iwft) = cdet(k, 1, 1)
        enddo
    enddo

end subroutine inputdet

subroutine inputlcao()
    ! Set the lcao to be equal
    use vmc_mod, only: MORB
    use coefs, only: coef, nbasis, norb
    use wfsec, only: nwftype
    use method_opt, only: method

    implicit none
    integer             :: iwft, i,j


    if( (method(1:3) == 'lin')) then
        if (.not. allocated(coef)) allocate (coef(nbasis, MORB, 3))
    else
        if (.not. allocated(coef)) allocate (coef(nbasis, MORB, nwftype))
    endif

    do iwft = 2, nwftype
        do i = 1, norb
            do j = 1, nbasis
                coef(j, i, iwft) = coef(j, i, 1)
            enddo
        enddo
    enddo

end subroutine inputlcao

subroutine inputjastrow()
    ! Set the jastrow to be equal

    use jaspar, only: nspin1, nspin2
    use jaspar3, only: a, b, c, scalek
    use jaspar4, only: a4, norda, nordb, nordc
    use bparm, only: nspin2b
    use contr2, only: ifock, ijas
    use contr2, only: isc
    use wfsec, only: nwftype
    use atom, only: ncent, nctype

      implicit none

      interface
      function nterms4(nord)
          integer, intent(in) :: nord
          integer :: nterms4
      end function nterms4
      end interface

      integer :: iparm, isp, it, iwft, mparmja
      integer :: mparmjb, mparmjc


    if (.not. allocated(scalek)) allocate (scalek(nwftype))

    if (ijas .ge. 4 .and. ijas .le. 6) then
        mparmja = 2 + max(0, norda - 1)
        mparmjb = 2 + max(0, nordb - 1)
        mparmjc = nterms4(nordc)

        if (.not. allocated(a4)) allocate (a4(mparmja, nctype, nwftype))
        if (.not. allocated(b))  allocate (b(mparmjb, 2, nwftype))
        if (.not. allocated(c))  allocate (c(mparmjc, nctype, nwftype))

        do iwft = 2, nwftype
            scalek(iwft) = scalek(1)
            do it = 1, nctype
                do iparm = 1, mparmja
                    a4(iparm, it, iwft) = a4(iparm, it, 1)
                enddo
            enddo

            do isp = nspin1, nspin2b
                do iparm = 1, mparmjb
                    b(iparm, isp, iwft) = b(iparm, isp, 1)
                enddo
            enddo

            do it = 1, nctype
                do iparm = 1, mparmjc
                    c(iparm, it, iwft) = c(iparm, it, 1)
                enddo
            enddo
        enddo
    endif

end subroutine inputjastrow

subroutine set_displace_zero(nforce_tmp)
!    use vmc_mod, only: MCENT
    use pcm, only: MCHS
    use forcestr, only: delc
    use pcm_force, only: sch_s
    use pcm_cntrl, only: ipcm
    use pcm_parms, only: ch, nchs

    use atom, only: ncent

    implicit none
    integer         :: i, j, nforce_tmp

    if (.not. allocated(delc)) allocate (delc(3, ncent, nforce_tmp))

    delc = 0.d0

    if (.not. allocated(sch_s)) allocate (sch_s(MCHS, nforce_tmp))

    if (ipcm .eq. 3) then
        do i = 1, nforce_tmp
            do j = 1, nchs
                sch_s(j, i) = ch(j)
            enddo
        enddo
    endif

    return
end subroutine set_displace_zero

subroutine modify_zmat_define

    use vmc_mod, only: MCENT
    use grdntsmv, only: igrdmv
    use atom, only: ncent
      implicit none

      integer :: ic, k


    !call p2gti('atoms:natom', ncent, 1)
    ! if (ncent .gt. MCENT) call fatal_error('MODIFY_ZMATRIX: ncent > MCENT')

    if (.not. allocated(igrdmv)) allocate (igrdmv(3, ncent))

    do ic = 1, ncent
        do k = 1, 3
            igrdmv(k, ic) = 1
        enddo
    enddo

    return
end subroutine modify_zmat_define

subroutine hessian_zmat_define

    use vmc_mod, only: MCENT
    use grdnthes, only: hessian_zmat
    use atom, only: ncent

      implicit none

      integer :: ic, k


    !call p2gti('atoms:natom', ncent, 1)
    ! if (ncent .gt. MCENT) call fatal_error('HESSIAN_ZMATRIX: ncent > MCENT')

    if (.not. allocated(hessian_zmat)) allocate (hessian_zmat(3, ncent))

    do ic = 1, ncent
        do k = 1, 3
            hessian_zmat(k, ic) = 1.d0
        enddo
    enddo

    return
end subroutine hessian_zmat_define
