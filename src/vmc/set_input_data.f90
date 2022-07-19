module set_input_data
use error, only: fatal_error
contains
subroutine inputzex
    ! Set the exponents to one when using a numerical basis
    use multiple_geo, only: MWF
    use numbas, only: numr
    use coefs, only: nbasis
    use basis, only: zex

    ! are they needed ??!!
    use contrl_per, only: iperiodic
    use wfsec, only: nwftype
    use method_opt, only: method
    use precision_kinds, only: dp
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

    use multiple_geo, only: MFORCE, MFORCE_WT_PRD, MWF
    use vmc_mod, only: nrad, nordj, nordj1, nmat_dim, nmat_dim2
    use vmc_mod, only: radmax, delri
    use vmc_mod, only: neqsx
    use csfs, only: cxdet, iadet, ibdet, icxdet, ncsf, nstates
    use dets, only: cdet, ndet
    use multidet, only: iactv, irepcol_det, ireporb_det, ivirt, iwundet, kref, numrep_det, allocate_multidet
    use multidet, only: k_det, ndetiab, ndet_req, k_det2, k_aux, ndetiab2, ndetsingle, kref_old
    use coefs, only: norb
    use dorb_m, only: iworbd

    use contrl_file,    	only: ounit, errunit

    ! not sure about that one either ....
    use wfsec, only: nwftype
    use multideterminant_mod, only: idiff
      use system, only: nelec
      use system, only: nup
      use system, only: ndn

    implicit none

    integer :: i, iab, icheck, icsf, idist
    integer :: iflag, in, iphase, iref
    integer :: irep, isav, ish, istate
    integer :: isub, iw, iwf, iwref
    integer :: j, k, l
    integer :: ndet_dist, nel, kk, kun, kaux, naux, kkn, kn
    integer, dimension(nelec) :: iswapped
    integer, dimension(ndet) :: itotphase
    integer, dimension(nelec) :: auxdet


    if (nup .lt. nelec/2) call fatal_error('INPUT: nelec/2 exceeds nup')
    ndn = nelec - nup

    !!call p2gtid('general:nwftype', nwftype, 1, 1)
    if (nwftype .gt. MWF) call fatal_error('INPUT: nwftype exceeds MWF')
    

! This part remains just for the initialization calls set kref and kref_old first time    
    if (iflag .eq. 0) then
       kref = 1
       kref_old = kref

    endif
!    write (ounit, *) 'I am in multi_def',kref,kref_old
       
    if (.not. allocated(iwundet)) allocate (iwundet(ndet, 2))
    if (.not. allocated(numrep_det)) allocate (numrep_det(ndet, 2))
    if (.not. allocated(irepcol_det)) allocate (irepcol_det(nelec, ndet, 2))
    if (.not. allocated(ireporb_det)) allocate (ireporb_det(nelec, ndet, 2))

    do iab = 1, 2
        numrep_det(kref, iab) = 0
    enddo

    do k = 1, ndet
        itotphase(k) = 0
        if (k .ne. kref) then
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
        endif
    enddo

    do k = 1, ndet
       if (k .ne. kref) then
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
        endif
    enddo

    call allocate_multidet()

    iactv(1) = nup + 1
    iactv(2) = ndn + 1
    ivirt(1) = nup + 1
    ivirt(2) = ndn + 1
    do k = 1, ndet
        if (k .ne. kref) then
        do iab = 1, 2
            do irep = 1, numrep_det(k, iab)
        if (irepcol_det(irep, k, iab) .ne. 0 .and. irepcol_det(irep, k, iab) .lt. iactv(iab)) iactv(iab)=irepcol_det(irep, k, iab)
                if (ireporb_det(irep, k, iab) .lt. ivirt(iab)) ivirt(iab) = ireporb_det(irep, k, iab)
            enddo
        enddo
        endif
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
                if (i .ne. kref) then
                if (idiff(kref, i, iab) .eq. 0) then
                   iwundet(i, iab) = kref
                 else
                    j=1
                    do while (idiff(j, i, iab) .ne. 0  .and. j.ne.i)
                       j= j+1
                    enddo
                    iwundet(i, iab) = j
                 endif
              endif
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

    
    !reshufling arrays to avoid redundancy of unequivalent determinats
    k_det=0
    ndetiab=0
    do iab = 1, 2
       kk=0
       kaux=0
       do k = 1, ndet
          if(iwundet(k,iab).eq.k.and.k.ne.kref) then
             if(numrep_det(k, iab).gt.0) then
                
                kk=kk+1
                k_det(k,iab)=kk
    

                auxdet=irepcol_det(:, kk, iab)                
                irepcol_det(:, kk, iab) = irepcol_det(:, k, iab)
                irepcol_det(:, k, iab) = auxdet
                
                auxdet=ireporb_det(:, kk, iab)
                ireporb_det(:, kk, iab) = ireporb_det(:, k, iab)
                ireporb_det(:, k, iab) = auxdet
                
                naux=numrep_det(kk, iab)
                numrep_det(kk, iab)=numrep_det(k, iab)
                numrep_det(k, iab)=naux
                
                
             endif
          endif
       enddo
       
       
       ndetiab(iab)=kk

!       print*,"multiple",iab,kk,kaux
!       print*,"ndetiab",iab,ndetiab(iab)
!    ordering first excitations at the begining for specialization
       kkn=0
       do kk=1,ndetiab(iab)

          if (numrep_det(kk, iab).eq.1) then
             kkn=kkn+1
             
             k=0
             do while (k_det(k,iab).ne.kk)
                k=k+1
             enddo
             kaux=k_det(k,iab)

             kn=0
             do while (k_det(kn,iab).ne.kkn)
                kn=kn+1
             enddo
             
!             print*,"kk",kk,"k_det(kk,",iab,")",kaux
             k_det(k,iab)=kkn
             k_det(kn,iab)=kk


             auxdet=irepcol_det(:, kkn, iab)                
             irepcol_det(:, kkn, iab) = irepcol_det(:, kk, iab)
             irepcol_det(:, kk, iab) = auxdet
                
             auxdet=ireporb_det(:, kkn, iab)
             ireporb_det(:, kkn, iab) = ireporb_det(:, kk, iab)
             ireporb_det(:, kk, iab) = auxdet
                
             naux=numrep_det(kkn, iab)
             numrep_det(kkn, iab)=numrep_det(kk, iab)
             numrep_det(kk, iab)=naux
             
             
          
          endif
          
       enddo

       ndetsingle(iab)=kkn
       
       !       print*,"kkn singles",iab,kkn
!       print*,"kkn singles",iab,ndetsingle(iab)
       
       
       
    enddo

    !setting larger number of required determinants for unequivalent or unique determinants 
    if(ndetiab(1).le.ndetiab(2)) then
       ndet_req=ndetiab(1)
    else
       ndet_req=ndetiab(2)
    endif

    !arrays for all not equivalent to kref
    do iab = 1, 2
       kk=0
       do k = 1, ndet
          kun=iwundet(k,iab)
          if(kun.ne.kref.and.k.ne.kref) then
!          if(numrep_det(k,iab).gt.1) then
             kk=kk+1
             k_det2(kk,iab)=k
             k_aux(kk,iab)=k_det(kun,iab)
!          endif
       endif
       enddo
       ndetiab2(iab)=kk
    enddo
    
    
    return
end subroutine multideterminants_define

subroutine inputforces
! Set all force displacements to zero
!    use multiple_geo, only: MWF
!    use multiple_geo, only: MFORCE
    use forcepar, only: nforce
    use forcestr, only: delc
    use wfsec, only: iwftype, nwftype
    use contrl_file, only: errunit
    use system, only: ncent
    use precision_kinds, only: dp

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
!    use mstates_mod, only: MSTATES
    use wfsec, only: nwftype
    use method_opt, only: method
    use precision_kinds, only: dp

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
    use vmc_mod, only: norb_tot
    use coefs, only: coef, nbasis, norb
    use wfsec, only: nwftype
    use method_opt, only: method
    use precision_kinds, only: dp

    implicit none
    integer             :: iwft, i,j


    if( (method(1:3) == 'lin')) then
        if (.not. allocated(coef)) allocate (coef(nbasis, norb_tot, 3))
    else
        if (.not. allocated(coef)) allocate (coef(nbasis, norb_tot, nwftype))
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
    use jaspar3, only: b, c, scalek
    use jaspar4, only: a4, norda, nordb, nordc
    use bparm, only: nspin2b
    use contr2, only: ijas
    use contr2, only: isc
    use wfsec, only: nwftype
    use system, only: ncent, nctype
    use precision_kinds, only: dp
    use jastrow4_mod, only: nterms4

      implicit none

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
    use pcm, only: MCHS
    use forcestr, only: delc
    use pcm_force, only: sch_s
    use pcm_cntrl, only: ipcm
    use pcm_parms, only: ch, nchs

    use system, only: ncent
    use precision_kinds, only: dp

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

    use grdntsmv, only: igrdmv
    use system, only: ncent
    implicit none

    integer :: ic, k


    if (.not. allocated(igrdmv)) allocate (igrdmv(3, ncent))

    do ic = 1, ncent
        do k = 1, 3
            igrdmv(k, ic) = 1
        enddo
    enddo

    return
end subroutine modify_zmat_define

subroutine hessian_zmat_define

    use grdnthes, only: hessian_zmat
    use system, only: ncent
    use precision_kinds, only: dp

    implicit none

    integer :: ic, k

    if (.not. allocated(hessian_zmat)) allocate (hessian_zmat(3, ncent))

    do ic = 1, ncent
        do k = 1, 3
            hessian_zmat(k, ic) = 1.d0
        enddo
    enddo

    return
end subroutine hessian_zmat_define
end module 
