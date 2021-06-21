
subroutine header_printing()
    !> This subroutine prints the header in each output file. It contains some
    !! useful information about the compilers, version of the code, input and output file names.
    !! @author Ravindra Shinde (r.l.shinde@utwente.nl)

    use mpi
    use mpiconf, only: idtask, nproc
    use, intrinsic :: iso_fortran_env, only: iostat_end
    use contrl_file,    only: file_input, file_output, file_error
    use contrl_file,    only: ounit, errunit

    implicit none

    integer                             :: status, i
    character(len=8)                    :: date
    character(len=10)                   :: time
    character(len=40)                   :: env_variable
    character(len=100)                  :: input_filename, output



    write(ounit,*) "____________________________________________________________________"
    write(ounit,*)
    write(ounit,*)
    write(ounit,*) ' .d8888b.   888    888         d8888  888b     d888  8888888b. '
    write(ounit,*) 'd88P  Y88b  888    888        d88888  8888b   d8888  888   Y88b'
    write(ounit,*) '888    888  888    888       d88P888  88888b.d88888  888    888'
    write(ounit,*) '888         8888888888      d88P 888  888Y88888P888  888   d88P'
    write(ounit,*) '888         888    888     d88P  888  888 Y888P 888  8888888P" '
    write(ounit,*) '888    888  888    888    d88P   888  888  Y8P  888  888       '
    write(ounit,*) 'Y88b  d88P  888    888   d8888888888  888   "   888  888       '
    write(ounit,*) ' "Y8888P"   888    888  d88P     888  888       888  888       '
    write(ounit,*)
    write(ounit,*) "____________________________________________________________________"
    write(ounit,*)
    write(ounit,*) ' Cornell Holland Ab-initio Materials Package'
    write(ounit,*)
    write(ounit,*)

    write(ounit,*) " information about the contributors goes here"
    write(ounit,*)
    write(ounit,*)
    write(ounit,*) " https://github.com/filippi-claudia/champ"
    write(ounit,*)
    write(ounit,*)

    write(ounit,*) " paper to cite for this code goes here"
    write(ounit,*)
    write(ounit,*)
    write(ounit,*)
    write(ounit,*)

    write(ounit,*) " license information goes here"

    write(ounit,*) "____________________________________________________________________"
    write(ounit,*)
    write(ounit,*)
    write(ounit,*)
    write(ounit,*)

    call date_and_time(date=date,time=time)
    write(ounit, '(12a)') " Calculation started on     :: ",  &
                            date(1:4), "-", date(5:6), "-", date(7:8), " at ",  time(1:2), ":", time(3:4), ":", time(5:6)
    call get_command_argument(number=0, value=output)
    write(ounit, '(2a)') " Executable                 :: ",   output

#if defined(GIT_BRANCH)
    write(ounit,'(2a)')  " Git branch                 :: ", GIT_BRANCH
#endif

#if defined(GIT_HASH)
    write(ounit,'(2a)')  " Git commit hash            :: ", GIT_HASH
#endif

#if defined(COMPILER)
    write(ounit,'(2a)')  " Compiler                   :: ", COMPILER
#endif

#if defined(COMPILER_VERSION)
    write(ounit,'(2a)')  " Compiler version           :: ", COMPILER_VERSION
#endif

    call hostnm(output)
    write(ounit, '(2a)') " Hostname                   :: ",   output
    call get_environment_variable ("PWD", output)
    write(ounit, '(2a)') " Current directory          :: ",   output
    call get_environment_variable ("USER", output)
    write(ounit, '(2a)') " Username                   :: ",   output
    write(ounit, '(2a)') " Input file                 :: ",   file_input
    write(ounit, '(2a)') " Output file                :: ",   file_output
    write(ounit, '(2a)') " Error file                 :: ",   file_error
    write(ounit, '(4a)') " Code compiled on           :: ",__DATE__, " at ", __TIME__
    write(ounit, '(a,i0)') " Number of processors       :: ", nproc
    write(ounit,*)



end subroutine header_printing


subroutine read_molecule_file(file_molecule)
    !> This subroutine reads the .xyz molecule file. It then computes the
    !! number of types of atoms, nuclear charges (from the symbol), and
    !! number of valence electrons if pseudopotential is provided.
    !! @author Ravindra Shinde (r.l.shinde@utwente.nl)
    !! @date
    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid
    use atom,               only: znuc, cent, pecent, iwctype, nctype, ncent, ncent_tot, nctype_tot, symbol, atomtyp
    use ghostatom, 		    only: newghostype, nghostcent
    use inputflags,         only: igeometry
    use periodic_table,     only: atom_t, element
    use contrl_file,        only: ounit, errunit

    implicit none

    !   local use
    character(len=72), intent(in)   :: file_molecule
    character(len=40)               :: temp1, temp2, temp3, temp4
    character(len=80)               :: comment
    integer                         :: iostat, i, j, k, iunit
    logical                         :: exist
    type(atom_t)                    :: atoms
    character(len=2), allocatable   :: unique(:)

    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: float_format   = '(A, T60, f12.8)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading
    write(ounit,*) '-----------------------------------------------------------------------'
    write(ounit,string_format)  " Reading molecular coordinates from the file :: ",  trim(file_molecule)
    write(ounit,*) '-----------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_molecule, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_molecule, iostat=iostat, action='read' )
            if (iostat .ne. 0) stop "Problem in opening the molecule file"
        else
            call fatal_error (" molecule file "// trim(file_molecule) // " does not exist.")
        endif

        read(iunit,*) ncent
    endif
    call bcast(ncent)

    write(ounit,fmt=int_format) " Number of atoms ::  ", ncent
    write(ounit,*)

    if (.not. allocated(cent)) allocate(cent(3,ncent))
    if (.not. allocated(symbol)) allocate(symbol(ncent))
    if (.not. allocated(iwctype)) allocate(iwctype(ncent))
    if (.not. allocated(unique)) allocate(unique(ncent))

    if (wid) read(iunit,'(A)')  comment
    call bcast(comment)

    write(ounit,*) "Comment from the molecule file :: ", trim(comment)
    write(ounit,*)

    if (wid) then
        do i = 1, ncent
            read(iunit,*) symbol(i), cent(1,i), cent(2,i), cent(3,i)
        enddo
    endif
    call bcast(symbol)
    call bcast(cent)

    if (wid) close(iunit)


    ! Count unique type of elements
    nctype = 1
    unique(1) = symbol(1)
    do j= 2, ncent
        if (any(unique == symbol(j) ))  cycle
        nctype = nctype + 1
        unique(nctype) = symbol(j)
    enddo

    write(ounit,fmt=int_format) " Number of distinct types of elements (nctype) :: ", nctype
    write(ounit,*)

    if (.not. allocated(atomtyp)) allocate(atomtyp(nctype))
    if (.not. allocated(znuc)) allocate(znuc(nctype))

    ! get the correspondence for each atom according to the rule defined for atomtypes
    do j = 1, ncent
        do k = 1, nctype
            if (symbol(j) == unique(k))   iwctype(j) = k
        enddo
    enddo

    ! Get the correspondence rule
    do k = 1, nctype
        atomtyp(k) = unique(k)
    enddo

    if (allocated(unique)) deallocate(unique)

    ! Get the znuc for each unique atom
    do j = 1, nctype
        atoms = element(atomtyp(j))
        znuc(j) = atoms%nvalence
    enddo

    ncent_tot = ncent + nghostcent
    nctype_tot = nctype + newghostype

    write(ounit,*) '-----------------------------------------------------------------------'
    write(ounit,'(a, t15, a, t27, a, t39, a, t45, a)') 'Symbol', 'x', 'y', 'z', 'Type'
    write(ounit,'(t14, a, t26, a, t38, a )') '(A)', '(A)', '(A)'
    write(ounit,*) '-----------------------------------------------------------------------'

    do j= 1, ncent
        write(ounit,'(A4, 2x, 3F12.6, 2x, i3)') symbol(j), (cent(i,j),i=1,3), iwctype(j)
    enddo

    write(ounit,*) '-----------------------------------------------------------------------'
    write(ounit,*) " Values of znuc (number of valence electrons) "
    write(ounit,'(10F12.6)') (znuc(j), j = 1, nctype)
    write(ounit,*) '-----------------------------------------------------------------------'
    write(ounit,*)
end subroutine read_molecule_file


subroutine read_determinants_file(file_determinants)
    !> This subroutine reads the single state determinant file.
    !! @author Ravindra Shinde

    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid
    use, intrinsic :: iso_fortran_env, only: iostat_eor
    use contrl_file,    only: ounit, errunit
    use dets,           only: cdet, ndet
    use vmc_mod,        only: MDET
    use dorb_m,         only: iworbd
    use coefs,          only: norb
    use inputflags,     only: ideterminants
    use wfsec,          only: nwftype
    use csfs,           only: nstates

    use elec,           only: ndn, nup
    use const,          only: nelec

    implicit none

    !   local use
    character(len=72), intent(in)   :: file_determinants
    character(len=80)               :: temp1, temp2, temp3
    integer                         :: iostat, i, j, iunit, counter
    logical                         :: exist, skip = .true.

    !   Formatting
    character(len=100)               :: int_format     = '(A, T40, I8)'
    character(len=100)               :: string_format  = '(A, T40, A)'

    !   External file reading
    write(ounit,*) '------------------------------------------------------'
    write(ounit,string_format)  " Reading determinants from the file :: ",  trim(file_determinants)
    write(ounit,*) '------------------------------------------------------'

    if (wid) then
        inquire(file=file_determinants, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_determinants, iostat=iostat, action='read' )
            if (iostat .ne. 0) stop "Problem in opening the determinant file"
        else
            call fatal_error (" determinant file "// trim(file_determinants) // " does not exist.")
        endif
    endif

    ndn  = nelec - nup

    write(ounit,*)
    write(ounit,int_format) " Number of total electrons ", nelec
    write(ounit,int_format) " Number of alpha electrons ", nup
    write(ounit,int_format) " Number of beta  electrons ", ndn
    write(ounit,*)


    ! to escape the comments before the "determinants" line
    if (wid) then
        do while (skip)
            read(iunit,*, iostat=iostat) temp1
            temp1 = trim(temp1)
            if (temp1 == "determinants") then
                backspace(iunit)
                skip = .false.
            endif
        enddo
    endif

!   Read the first main line
    if (wid) then
        read(iunit, *, iostat=iostat)  temp2, ndet, nwftype
        if (iostat == 0) then
            if (trim(temp2) == "determinants") write(ounit,int_format) " Number of determinants ", ndet
        else
            call fatal_error ("Error in reading number of determinants / number of wavefunction types")
        endif
    endif
    call bcast(ndet)
    call bcast(nwftype)

    ! Note the hack here about capitalized variables. DEBUG
    MDET = ndet

    if (.not. allocated(cdet)) allocate(cdet(ndet,1,nwftype))

    if (wid) then
        read(iunit,*, iostat=iostat) (cdet(i,1,1), i=1,ndet)
        if (iostat /= 0) call fatal_error( "Error in determinant coefficients ")
    endif
    call bcast(cdet)

    write(ounit,*)
    write(ounit,*) " Determinant coefficients "
    write(ounit,'(10(1x, f11.8, 1x))') (cdet(i,1,1), i=1,ndet)

!       allocate the orbital mapping array
    if (.not. allocated(iworbd)) allocate(iworbd(nelec, ndet))

    if (wid) then
        do i = 1, ndet
            read(iunit,*, iostat=iostat) (iworbd(j,i), j=1,nelec)
            if (iostat /= 0) call fatal_error("Error in reading orbital -- determinants mapping ")
        enddo
    endif
    call bcast(iworbd)
    ! This part replaces a call to verify_orbitals
    !if(any(iworbd .gt. norb))  call fatal_error('INPUT: iworbd > norb')


    write(ounit,*)
    write(ounit,*) " Orbitals <--> Determinants mapping :: which orbitals enter in which dets"
    do i = 1, ndet
        write(ounit,'(<nelec>(i4, 1x))') (iworbd(j,i), j=1,nelec)
    enddo

    if (wid) then
        read(iunit,*) temp1
        if (temp1 == "end" ) write(ounit,*) " Single state determinant file read successfully "
        ideterminants = ideterminants + 1
    endif
    call bcast(ideterminants)

    if (wid) close(iunit)

end subroutine read_determinants_file

subroutine read_multideterminants_file(file_multideterminants)
    !> This subroutine reads the multideterminants file. The first line appears as 'multideterminants' ndet_local
    !! CI coefficients and occupation of determinants in wf
    !! @author Ravindra Shinde
    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid
    use contrl_file,    only: ounit, errunit
    use dets, only: ndet
    use const, only: nelec
    use multidet, only: irepcol_det, ireporb_det, numrep_det, iwundet
    use inputflags, only: imultideterminants

    implicit none

    !   local use
    character(len=72), intent(in)   :: file_multideterminants
    character(len=80)               :: temp1, temp2, temp3
    integer                         :: iostat, k, iunit, ndet_local, iab, irep
    logical                         :: exist, skip = .true.

    !   Formatting
    character(len=100)               :: int_format     = '(A, T40, I8)'
    character(len=100)               :: string_format  = '(A, T40, A)'

    !   External file reading
    write(ounit,*) '------------------------------------------------------'
    write(ounit,string_format)  " Reading multideterminants from the file :: ",  trim(file_multideterminants)
    write(ounit,*) '------------------------------------------------------'

    if (wid) then
        inquire(file=file_multideterminants, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_multideterminants, iostat=iostat, action='read' )
            if (iostat .ne. 0) stop "Problem in opening the multideterminant file"
        else
            call fatal_error (" Multideterminant file "// trim(file_multideterminants) // " does not exist.")
        endif

        read (iunit, *, iostat=iostat) temp1, ndet_local
        if (iostat /= 0) call fatal_error("Error in reading multideterminant file :: expecting 'multideterminants', ndet")

        if (trim(temp1) /= "multideterminants") then
            call fatal_error ("Error in reading multideterminant file :: expecting 'multideterminants'")
        elseif (ndet_local .ne. ndet ) then
            call fatal_error ('Error: ndet not matching with previous records')
        endif
    endif

    if (.not. allocated(iwundet)) allocate(iwundet(ndet, 2))
    if (.not. allocated(numrep_det)) allocate(numrep_det(ndet, 2))
    if (.not. allocated(irepcol_det)) allocate(irepcol_det(nelec, ndet, 2))
    if (.not. allocated(ireporb_det)) allocate(ireporb_det(nelec, ndet, 2))

    if (wid) then
        do k = 2, ndet_local
            read (iunit, *, iostat=iostat) (numrep_det(k, iab), iab=1, 2)
            do iab = 1, 2
                do irep = 1, numrep_det(k, iab)
                    read (iunit, *, iostat=iostat) irepcol_det(irep, k, iab), ireporb_det(irep, k, iab)
                enddo
            enddo
        enddo
        imultideterminants = imultideterminants + 1
    endif
    call bcast(numrep_det)
    call bcast(irepcol_det)
    call bcast(ireporb_det)
    call bcast(imultideterminants)

    if (wid) close(iunit)

end subroutine read_multideterminants_file




subroutine read_jastrow_file(file_jastrow)
    ! This subroutine reads jastrow parameters from a file.
    ! Ravindra

    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use, intrinsic :: iso_fortran_env, only: iostat_eor !, iostat_eof
    use contrl_file,    only: ounit, errunit

    use force_mod,          only: MWF
    use jaspar,             only: nspin1, nspin2
    use elec,               only: ndn
    use jaspar3,            only: a, b, c, scalek
    use jaspar4,            only: a4, norda, nordb, nordc
    use jaspar6,            only: cutjas, cutjasi, allocate_jaspar6
    use bparm,              only: nocuspb, nspin2b
    use contr2,             only: ifock, ijas
    use contr2,             only: isc
    use inputflags,         only: ijastrow_parameter
    use wfsec,              only: nwftype
    use atom,               only: ncent, nctype
    use precision_kinds,    only: dp
    use contrl_per, 		only: iperiodic
    use jaspar6, 			only: asymp_jasa, asymp_jasb, asymp_r, c1_jas6, c1_jas6i, c2_jas6
    implicit none

    !   local use
    character(len=72), intent(in)   :: file_jastrow
    character(len=40)               :: temp1, temp2, temp3, temp4, temp5
    integer                         :: iunit, iostat, it, isp, iparm, iwft
    integer                         :: mparmja, mparmjb, mparmjc, nterms4
    logical                         :: exist
    real(dp)                        :: a21, cutjas_tmp
    integer                         :: i, j

    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading
    write(ounit,*) '---------------------------------------------------------------------------'
    write(ounit,string_format)  " Reading jastrow parameters from the file :: ",  trim(file_jastrow)
    write(ounit,*) '---------------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_jastrow, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_jastrow, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error("Problem in opening the jastrow file")
        else
            call fatal_error (" Jastrow file "// trim(file_jastrow) // " does not exist.")
        endif
    endif


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

    ! read the first word of the file
    if (wid) then
        read (iunit, *, iostat=iostat)  temp2, iwft
        if (iostat == 0) then
            if (trim(temp2) == "jastrow_parameter") &
            write(ounit,int_format) " Jastrow parameters being read : type of wavefunctions :: ", iwft
        else
            call fatal_error ("Error in reading jastrow parameters / number of wavefunction types")
        endif
    endif
    call bcast(iwft)

    allocate (scalek(nwftype))

    if (ijas .ge. 4 .and. ijas .le. 6) then
        if (ifock .gt. 0) call fatal_error('JASTROW: fock not yet implemented for ijas=4,5,6')
        if (wid) read (iunit, *) norda, nordb, nordc
        call bcast(norda)
        call bcast(nordb)
        call bcast(nordc)

        write(ounit, '(3(A,i4))') " norda = ", norda, "; nordb = ", nordb, "; nordc = ", nordc

        if (isc .ge. 2) then
            if (wid) read (iunit, *) scalek(iwft), a21
        endif
        call bcast(scalek)
        call bcast(a21)
        write(ounit, '(2(A,f12.6))') " scalek = ", scalek(iwft), "; a21 = ", a21

        mparmja = 2 + max(0, norda - 1)
        mparmjb = 2 + max(0, nordb - 1)
        mparmjc = nterms4(nordc)

        allocate (a4(mparmja, nctype, nwftype))

        write(ounit, '(A)') "Jastrow parameters :: "
        write(ounit, '(A)') "mparmja : "
        do it = 1, nctype
            if (wid) read (iunit, *) (a4(iparm, it, iwft), iparm=1, mparmja)
            write(ounit, '(<mparmja>(2X,f12.8))') (a4(iparm, it, iwft), iparm=1, mparmja)
        enddo
        call bcast(a4)

        allocate (b(mparmjb, 2, nwftype))

        write(ounit, '(A)') "mparmjb : "
        do isp = nspin1, nspin2b
            if (wid) read (iunit, *) (b(iparm, isp, iwft), iparm=1, mparmjb)
            write(ounit, '(<mparmjb>(2X,f12.8))') (b(iparm, isp, iwft), iparm=1, mparmjb)
        enddo
        call bcast(b)

        allocate (c(mparmjc, nctype, nwftype))

        write(ounit, '(A)') "mparmjc : "
        do it = 1, nctype
            if (wid) read (iunit, *) (c(iparm, it, iwft), iparm=1, mparmjc)
            write(ounit, '(<mparmjc>(2X,f12.8))') (c(iparm, it, iwft), iparm=1, mparmjc)
        enddo
        call bcast(c)

    endif

    !Read cutoff for Jastrow4, 5, 6
    if (isc .eq. 6 .or. isc .eq. 7) then
        if (wid) read (iunit, *) cutjas
        write(iunit, '(A,2X,f12.8)') " cutjas = ", cutjas
    endif
    call bcast(cutjas)

    ijastrow_parameter = ijastrow_parameter + 1
    call bcast(ijastrow_parameter)

    if (wid) close(iunit)

end subroutine read_jastrow_file


subroutine read_orbitals_file(file_orbitals)
    ! Ravindra

    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use contrl_file,    only: ounit, errunit
    use coefs, only: coef, nbasis, norb
    use inputflags, only: ilcao
    use orbval, only: nadorb
    use pcm_fdc, only: fs
    use vmc_mod, only: MORB, MBASIS
    ! was not in master but is needed
    use wfsec, only: nwftype

    implicit none

!   local use
    character(len=72), intent(in)   :: file_orbitals
    character(len=40)               :: temp1, temp2
    character(len=120)              :: temp3
    integer                         :: iunit, iostat, iwft
    integer                         :: iorb, ibasis, i, k, counter
    logical                         :: exist
    logical                         :: skip = .true.

    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: string_format  = '(A, T60, A)'
    character(len=100)               :: float_format   = '(A, T60, f12.8)'

    !   External file reading
    write(ounit,*) '---------------------------------------------------------------------------'
    write(ounit,string_format)  " Reading LCAO orbitals from the file :: ",  trim(file_orbitals)
    write(ounit,*) '---------------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_orbitals, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_orbitals, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error("Problem in opening the LCAO orbitals file")
        else
            call fatal_error (" LCAO file "// trim(file_orbitals) // " does not exist.")
        endif
    endif

    ! to escape the comments before the "lcao nbasis norb" line
    if (wid) then
        do while (skip)
            read(iunit,*, iostat=iostat) temp1
            temp1 = trim(temp1)
            if (temp1 == "lcao") then
                backspace(iunit)
                skip = .false.
            endif
        enddo
    endif
    ! read the first line
    if (wid) read(iunit, *, iostat=iostat)  temp1, nbasis, norb, iwft
    call bcast(nbasis)
    call bcast(norb)
    call bcast(iwft)

    if (wid) then
        if (iostat == 0) then
            if (trim(temp1) == "lcao") then
                write(ounit,int_format) " Number of basis functions ", nbasis
                write(ounit,int_format) " Number of lcao orbitals ", norb
                write(ounit,int_format) " Type of wave functions ", iwft
                ! Note the hack with capitalized variables DEBUG
                MBASIS = nbasis
                MORB = norb
            endif
        else
            write(ounit, *) " Check ", temp1, nbasis, norb, iwft
            call fatal_error ("Error in reading number of lcao orbitals / basis / number of wavefunction types")
        endif
    endif
    call bcast(MBASIS)
    call bcast(MORB)

    if (iwft .gt. nwftype) call fatal_error('LCAO: wave function type > nwftype')

    if (.not. allocated(coef)) allocate (coef(nbasis, norb, nwftype))

    do iorb = 1, norb
        if (wid) then
            read (iunit, *, iostat=iostat) (coef(ibasis, iorb, iwft), ibasis=1, nbasis)
            if (iostat /= 0) call fatal_error( "Error in reading lcao orbitals ")
        endif
    enddo
    call bcast(coef)

    write(ounit,*)
    write(ounit,*) " LCAO orbitals "

    temp3 = '(T8, T14, i3, T28, i3, T42, i3, T56, i3, T70, i3, T84, i3, T98, i3, T112, i3, T126, i3, T140, i3)'
    ! print orbs in blocks of 10
    counter = 0
    do k = 10, nbasis, 10
!        write(ounit,*) " Orbitals  ", k-9 , "  to ", k
        write(ounit, fmt=temp3 )  (i, i = k-9, k)
        do iorb = 1, norb
            write(ounit, '(A,i5,A, 10(1x, f12.8, 1x))') "[", iorb, "] ", (coef(ibasis, iorb, iwft), ibasis=k-9, k)
        enddo
        counter = counter + 10
    enddo


    ! Remaining block
    write(ounit, fmt=temp3 )  (i, i = counter, nbasis)
    do iorb = 1, norb
        write(ounit, '(A,i5,A, 10(1x, f12.8, 1x))') "[", iorb, "] ", (coef(ibasis, iorb, iwft), ibasis=counter, nbasis)
    enddo
    ilcao = ilcao + 1

    call bcast(ilcao)
    if (wid) close(iunit)

    write(ounit,*) "----------------------------------------------------------"

end subroutine read_orbitals_file

subroutine read_csf_file(file_determinants)
    ! This subroutine reads the csf coefficients from the determinant file.
    ! Ravindra

    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use, intrinsic :: iso_fortran_env!, only: is_iostat_end
    use contrl_file,    only: ounit, errunit
    use vmc_mod, only: MDET
    use csfs, only: ccsf, ncsf, nstates
    use mstates_mod, only: MSTATES
    use inputflags, only: icsfs
    use wfsec, only: nwftype
    use dets, only: ndet, cdet
!   Not sure about the following two lines
    use ci000, only: nciprim, nciterm
    use optwf_contrl, only: ioptci
    implicit none

    !   local use
    character(len=72), intent(in)   :: file_determinants
    character(len=40)               :: temp1, temp2, temp3, temp4, temp5
    integer                         :: iostat, i, j, iunit
    logical                         :: exist, printed
    logical                         :: found = .false.

    !   Formatting
    character(len=100)              :: int_format     = '(A, T40, I8)'
    character(len=100)              :: string_format  = '(A, T40, A)'

    !   External file reading
    write(ounit,*) '------------------------------------------------------'
    write(ounit,string_format)  " Reading csf from the file :: ",  trim(file_determinants)
    write(ounit,*) '------------------------------------------------------'

    if (wid) then
        inquire(file=file_determinants, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_determinants, iostat=iostat, action='read' )
            if (iostat .ne. 0) stop "Problem in opening the determinant file for reading csfs"
        else
            call fatal_error (" determinant file "// trim(file_determinants) // " does not exist.")
        endif
    endif

    if (wid) then
        do while (.not. found)
            read(iunit,*, iostat=iostat) temp1
            if (is_iostat_end(iostat)) exit
            temp1 = trim(temp1)
            if (temp1 == "csf") then
                backspace(iunit)
                found = .true.
            endif
        enddo
    endif
    call bcast(found)

    ! if there is no mention of "csf" in the file
    if (.not. found) then
        ! No csf information present. One to one mapping cdet == ccsf
        nstates = 1
        ncsf = ndet
        if (ioptci .ne. 0) nciterm = nciprim
        printed = .true.
        ! if there is no mention of "csf" in the file:: allocate and assign
        if (.not. allocated(ccsf)) allocate(ccsf(ndet, nstates, nwftype))
        do i = 1, nstates
            do j = 1, ndet
                ccsf(j,i,nwftype) = cdet(j,i,nwftype)
            enddo
        enddo
        ! printing
        write(ounit,int_format) " Number of configuration state functions (csf) ", ncsf
        write(ounit,int_format) " Number of states (nstates) ", nstates

        write(ounit,*)
        write(ounit,*) " CSF coefficients not present in the determinant file "
        write(ounit,*) " CSF coefficients map one-to-one to determinant coefficients "
        write(ounit,*)

    else
        ! read the same line again to get the ncsf and nstates
        if (wid) read(iunit, *, iostat=iostat)  temp2, ncsf, nstates
        call bcast(ncsf)
        call bcast(nstates)

        if (.not. allocated(ccsf)) allocate(ccsf(ndet, nstates, nwftype))
        if (wid) then
            do i = 1, nstates
                read(iunit,*) (ccsf(j,i,1), j=1,ncsf)
            enddo
        endif
        call bcast(ccsf)

        write(ounit,int_format) " Number of configuration state functions (csf) ", ncsf
        write(ounit,int_format) " Number of states (nstates) ", nstates

        write(ounit,*)
        write(ounit,*) " CSF coefficients from an external file "

        write(ounit,'(10(1x, a9, i3, 1x))') ((" State: ", i), i =1, nstates)
        do j = 1, ncsf
            write(ounit,'(10(1x, f12.8, 1x))') (ccsf(j,i,1), i=1,nstates)
        enddo
    endif

    icsfs = 1
    call bcast(icsfs)

    if (wid) close(iunit)

end subroutine read_csf_file

subroutine read_csfmap_file(file_determinants)
    ! This subroutine reads the csf coefficients from the determinant file.
    ! Ravindra

    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use, intrinsic :: iso_fortran_env
    use contrl_file,    only: ounit, errunit
    use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates
    use dets, only: cdet, ndet
    use wfsec, only: nwftype
    use mstates_mod, only: MDETCSFX
    use precision_kinds,    only: dp

    implicit none

    !   local use
    character(len=72), intent(in)   :: file_determinants
    character(len=40)               :: temp1, temp2, temp3, temp4, temp5
    integer                         :: iostat, i, j, k, iunit
    integer                         :: icsf, jx
    integer                         :: nptr, nterm, id, nmap
    integer                         :: ncsf_check, ndet_check, nmap_check
    real(dp)                        :: c
    logical                         :: exist, printed
    logical                         :: found = .false.

    !   Formatting
    character(len=100)              :: int_format     = '(A, T40, I8)'
    character(len=100)              :: string_format  = '(A, T40, A)'

    !   External file reading
    write(ounit,*) '------------------------------------------------------'
    write(ounit,string_format)  " Reading csfmap from the file :: ",  trim(file_determinants)
    write(ounit,*) '------------------------------------------------------'

    if (wid) then
        inquire(file=file_determinants, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_determinants, iostat=iostat, action='read' )
            if (iostat .ne. 0) stop "Problem in opening the determinant file for reading csfmap"
        else
            call fatal_error (" determinant file "// trim(file_determinants) // " does not exist.")
        endif
    endif

    if (wid) then
        do while (.not. found)
            read(iunit,*, iostat=iostat) temp1
            if (is_iostat_end(iostat)) exit
            temp1 = trim(temp1)
            if (temp1 == "csfmap") then
                backspace(iunit)
                found = .true.
            endif
        enddo
    endif
    call bcast(found)

    ! if there is no mention of "csfmap" in the file
    if (.not. found) then
        ! No csfmap information present. One to one mapping cdet == ccsf
        ! Check this part carefully
        if (.not. allocated(cxdet)) allocate (cxdet(ndet*MDETCSFX))     ! why MDETCSFX
        if (.not. allocated(iadet)) allocate (iadet(ndet))
        if (.not. allocated(ibdet)) allocate (ibdet(ndet))
        if (.not. allocated(icxdet)) allocate (icxdet(ndet*MDETCSFX))   ! why MDETCSFX

        do i = 1, ncsf
            iadet(i) = i
            ibdet(i) = i
            icxdet(i) = i
            cxdet(i) = 1.0d0
        enddo

        write(ounit,*) " Determinant - CSF has one-to-one mapping  "
    else
        ! read from the file.
        if (wid) then
            read(iunit, *, iostat=iostat)  temp2, ncsf_check, ndet_check, nmap_check
            if (iostat == 0) then
                if (ndet_check .ne. ndet) call fatal_error('CSFMAP: wrong number of determinants')
                if (ncsf_check .ne. ncsf) call fatal_error('CSFMAP: wrong number of csf')
                if (nmap_check .gt. float(ndet)*ndet) call fatal_error('CSFMAP: too many determinants in map list')
            endif
        endif
        call bcast(ncsf_check)
        call bcast(ndet_check)
        call bcast(nmap_check)

        if (.not. allocated(cxdet)) allocate (cxdet(ndet*MDETCSFX))     ! why MDETCSFX
        if (.not. allocated(iadet)) allocate (iadet(ndet))
        if (.not. allocated(ibdet)) allocate (ibdet(ndet))
        if (.not. allocated(icxdet)) allocate (icxdet(ndet*MDETCSFX))   ! why MDETCSFX

        nptr = 1
        do i = 1, ncsf
            if (wid) read (iunit, *) nterm
            call bcast(nterm)
            iadet(i) = nptr
            ibdet(i) = nptr + nterm - 1
            do j = 1, nterm
                if (wid) read (iunit, *) id, c
                call bcast(id)
                call bcast(c)
                icxdet(nptr) = id
                cxdet(nptr) = c
                nptr = nptr + 1
                if (nptr .gt. ndet*MDETCSFX) call fatal_error ('Error in CSFMAP:: problem with nmap')
            enddo
        enddo
        if (nmap_check .ne. nptr - 1) call fatal_error ('Error in CSFMAP:: not enough nmaps / file is corrupt')
        nmap = nptr

        if (.not. allocated(cdet)) allocate (cdet(ndet, nstates, nwftype))

        write(ounit, '(''Warning: det coef overwritten with csf'')')

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

        write(ounit,int_format) " Number of configuration state functions (csf) ", ncsf
        write(ounit,int_format) " Number of determinants (ndet) ", ndet
        write(ounit,int_format) " Number of mappings (nmap) ", nmap
        write(ounit,*)

        write(ounit,*) " Determinant - CSF mapping  "

        do icsf = 1, ncsf
            write(ounit,'(i4)') icsf
            do j = iadet(icsf), ibdet(icsf)
                jx = icxdet(j)
                write(ounit,'(1(5x, i4, t12, f12.8, 1x))') icxdet(j), cxdet(j)
            enddo
        enddo
        write(ounit,*)
    endif

    if (wid) close(iunit)

    write(ounit,*) '------------------------------------------------------'

end subroutine read_csfmap_file




subroutine read_exponents_file(file_exponents)
    ! Read basis function exponents (only if no numerical basis)
    ! Ravindra

    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use contrl_file,        only: ounit, errunit
    use coefs,              only: nbasis
    use basis,              only: zex
    use inputflags,         only: iexponents
    use wfsec,              only: nwftype

    implicit none

    !   local use
    character(len=72), intent(in)   :: file_exponents
    character(len=40)               :: temp1, temp2
    integer                         :: iostat, i, iwft, iunit
    logical                         :: exist

    !   Formatting
    character(len=100)              :: int_format     = '(A, T40, I8)'
    character(len=100)              :: string_format  = '(A, T40, A)'

    !   External file reading
    write(ounit,*) '------------------------------------------------------'
    write(ounit,string_format)  " Reading exponents from the file :: ",  trim(file_exponents)
    write(ounit,*) '------------------------------------------------------'

    if (wid) then
        inquire(file=file_exponents, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_exponents, iostat=iostat, action='read' )
            if (iostat .ne. 0) stop "Problem in opening the exponents file for reading csfs"
        else
            call fatal_error (" exponents file "// trim(file_exponents) // " does not exist.")
        endif
    endif

    write(ounit, *) 'nbasis', nbasis
    write(ounit, *) 'nwftype', nwftype

    if (.not. allocated(zex)) allocate (zex(nbasis, nwftype))

    do iwft = 1, nwftype
        if (wid) then
            read(iunit,*, iostat=iostat)  (zex(i, iwft), i=1, nbasis)

            if (iostat /= 0) call fatal_error( "Error in reading exponents from the exponent file ")
        endif

        write(ounit,*)
        write(ounit,*) " Basis set exponents "

        write(ounit,'(10(1x, f11.8, 1x))') (zex(i, iwft), i=1, nbasis)
        write(ounit,*)
    enddo
    call bcast(zex)

    if (wid) close(iunit)

end subroutine read_exponents_file


subroutine read_jasderiv_file(file_jastrow_der)
    ! Read jastrow derivatives
    ! Ravindra

    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use contrl_file,        only: ounit, errunit
    use optjas,             only: MPARMJ
    use atom,               only: nctype
    use jaspar,             only: nspin1, is
    use jaspar4,            only: norda, nordb, nordc
    use jaspointer,         only: npoint, npointa
    use numbas,             only: numr

    use optwf_nparmj,       only: nparma, nparmb, nparmc, nparmf
    use optwf_parms,        only: nparmj
    use optwf_wjas,         only: iwjasa, iwjasb, iwjasc, iwjasf
    use bparm,              only: nspin2b
    use contr2,             only: ijas
    use contr2,             only: isc
    use vmc_mod,            only: MCTYP3X
    use atom,               only: nctype_tot

    implicit none

    !   local use
    character(len=72), intent(in)   :: file_jastrow_der
    character(len=40)               :: temp1, temp2, temp3, temp4, temp5
    integer                         :: iunit, iostat
    integer                         :: na1, na2, it, isp, iparm, ia
    logical                         :: exist, found = .false.
    !real(dp)                        ::

    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading
    write(ounit,*) '---------------------------------------------------------------------------'
    write(ounit,string_format)  " Reading jastrow derivative parameters from the file :: ",  trim(file_jastrow_der)
    write(ounit,*) '---------------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_jastrow_der, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_jastrow_der, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error( "Problem in opening the jastrow derivative file")
        else
            call fatal_error (" Jastrow derivative file "// trim(file_jastrow_der) // " does not exist.")
        endif
    endif

    na1 = 1
    na2 = nctype
    MCTYP3X = max(3, nctype) !BUG DEBUG nctype_tot

    if (.not. allocated(nparma)) allocate (nparma(MCTYP3X))
    if (.not. allocated(nparmb)) allocate (nparmb(3))
    if (.not. allocated(nparmc)) allocate (nparmc(nctype))
    if (.not. allocated(nparmf)) allocate (nparmf(nctype))

    if (.not. allocated(iwjasa)) allocate (iwjasa(83, MCTYP3X))
    if (.not. allocated(iwjasb)) allocate (iwjasb(83, 3))
    if (.not. allocated(iwjasc)) allocate (iwjasc(83, nctype))
    if (.not. allocated(iwjasf)) allocate (iwjasf(15, nctype))

    if (.not. allocated(npoint)) allocate (npoint(MCTYP3X))
    if (.not. allocated(npointa)) allocate (npointa(3*MCTYP3X))

    ! to escape the comments before the "jasderiv" line
    if (wid) then
        do while (.not. found)
            read(iunit,*, iostat=iostat) temp1
            if (is_iostat_end(iostat)) exit
            temp1 = trim(temp1)
            if (temp1 == "jasderiv") then
                found = .true.
            endif
        enddo
    endif
    call bcast(found)

    ! Read the file if jasderiv keyword is found
    if (.not. found) then
        call fatal_error ("Error in reading the first line of jastrow derivative file.")
    else
        if (wid) then
            read (iunit, *, iostat=iostat) (nparma(ia), ia=na1, na2), &
                            (nparmb(isp), isp=nspin1, nspin2b), &
                            (nparmc(it), it=1, nctype), &
                            (nparmf(it), it=1, nctype)
            if (iostat /= 0 ) error stop "Problem reading jasderiv file"
        endif

        call bcast(nparma)
        call bcast(nparmb)
        call bcast(nparmc)
        call bcast(nparmf)
        write(ounit, '(A,10i4)') " nparma = ", (nparma(ia), ia=na1, na2)
        write(ounit, '(A,10i4)') " nparmb = ", (nparmb(isp), isp=nspin1, nspin2b)
        write(ounit, '(A,10i4)') " nparmc = ", (nparmc(it), it=1, nctype)
        write(ounit, '(A,10i4)') " nparmf = ", (nparmf(it), it=1, nctype)

        if (ijas .ge. 4 .and. ijas .le. 6) then
            do it = 1, nctype
                if (numr .eq. 0) then
                    ! All-electron with analytic slater basis
                    if ((nparma(it) .gt. 0 .and. norda .eq. 0) .or. (nparma(it) .gt. norda + 1)) then
                        write(ounit, '(''it,norda,nparma(it)'',3i5)') it, norda, nparma(it)
                        call fatal_error( 'nparma too large for norda')
                    endif
                else
                    ! Pseudopotential with numerical basis: cannot vary a(1) or a(2)
                    if (norda .eq. 1) call fatal_error( 'makes no sense to have norda=1 for numr>0')
                    if ((norda .eq. 0 .and. nparma(it) .gt. 0) .or. (norda .gt. 0 .and. nparma(it) .gt. norda - 1)) then
                        write(ounit, '(''it,norda,nparma(it)'',3i5)') it, norda, nparma(it)
                        call fatal_error( 'nparma too large for norda')
                    endif
                endif

                if (isc .le. 7 .and. &
                    ((nordc .le. 2 .and. nparmc(it) .gt. 0) &
                        .or. (nordc .eq. 3 .and. nparmc(it) .gt. 2) &
                        .or. (nordc .eq. 4 .and. nparmc(it) .gt. 7) &
                        .or. (nordc .eq. 5 .and. nparmc(it) .gt. 15) &
                        .or. (nordc .eq. 6 .and. nparmc(it) .gt. 27) &
                        .or. (nordc .eq. 7 .and. nparmc(it) .gt. 43))) then
                    write(ounit, '(''it,nordc,nparmc(it)'',3i5)') it, nordc, nparmc(it)
                    call fatal_error( 'nparmc too large for nordc in J_een with cusp conds')
                endif

                if (isc .gt. 7 .and. &
                    ((nordc .le. 1 .and. nparmc(it) .gt. 0) &
                        .or. (nordc .eq. 2 .and. nparmc(it) .gt. 2) &
                        .or. (nordc .eq. 3 .and. nparmc(it) .gt. 6) &
                        .or. (nordc .eq. 4 .and. nparmc(it) .gt. 13) &
                        .or. (nordc .eq. 5 .and. nparmc(it) .gt. 23) &
                        .or. (nordc .eq. 6 .and. nparmc(it) .gt. 37) &
                        .or. (nordc .eq. 7 .and. nparmc(it) .gt. 55))) then
                    write(ounit, '(''it,nordc,nparmc(it)'',3i5)') it, nordc, nparmc(it)
                    call fatal_error( 'nparmc too large for nordc without cusp conds')
                endif

            enddo

            ! For the b coefs. we assume that b(1) is fixed by the cusp-cond.
            do isp = 1, nspin1, nspin2b
                if (nparmb(isp) .gt. nordb) then
                    write(ounit, '(''isp,nordb,nparmb(isp)'',3i5)') isp, nordb, nparmb(isp)
                    call fatal_error( 'nparmb too large for nordb')
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

        call bcast(nparmj)
        if (nparmj .gt. MPARMJ) call fatal_error('JASDERIV: MPARMJ too small')

        do it = 1, nctype
            if (wid) read (iunit, *) (iwjasa(iparm, it), iparm=1, nparma(it))
            write(ounit, '(A,10i4)') " iwjasa = ", (iwjasa(iparm, it), iparm=1, nparma(it))
        enddo
        call bcast(iwjasa)
        do isp = nspin1, nspin2b
            if (wid) read (iunit, *) (iwjasb(iparm, isp), iparm=1, nparmb(isp))
            write(ounit, '(A,10i4)') " iwjasb = ", (iwjasb(iparm, isp), iparm=1, nparmb(isp))
        enddo
        call bcast(iwjasb)
        do it = 1, nctype
            if (wid) read (iunit, *) (iwjasc(iparm, it), iparm=1, nparmc(it))
            write(ounit, '(A,10i4)') " iwjasc = ", (iwjasc(iparm, it), iparm=1, nparmc(it))
        enddo
        call bcast(iwjasc)
            ! end of reading the jasderiv file block

    endif ! if the file containing jasderiv exists

    if (wid) close(iunit)

end subroutine read_jasderiv_file


subroutine read_forces_file(file_forces)
    !
    ! Ravindra
    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use atom,               only: symbol
    use contrl_file,        only: ounit, errunit
    use forcepar,           only: nforce
    use forcestr,           only: delc
    use wfsec,              only: iwftype
    use inputflags,         only: iforces

    use atom,               only: ncent

    implicit none

    !   local use
    character(len=72), intent(in)   :: file_forces
    character(len=40)               :: temp1, temp2, temp3, temp4, temp5
    integer                         :: iunit, iostat
    integer                         :: i,ic,j, k
    logical                         :: exist, skip = .true.
    !real(dp)                        ::

    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading
    write(ounit,*) '-----------------------------------------------------------------------'
    write(ounit,string_format)  " Reading force displacements from the file :: ",  trim(file_forces)
    write(ounit,*) '-----------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_forces, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_forces, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error( "Problem in opening the forces file")
        else
            call fatal_error (" Forces file "// trim(file_forces) // " does not exist.")
        endif
    endif

    if (.not. allocated(delc)) allocate (delc(3, ncent, nforce))
    if (.not. allocated(iwftype)) allocate (iwftype(nforce))

    if (wid) then
        read (iunit, *, iostat=iostat) (iwftype(i), i=1, nforce)
        if (iostat /= 0) call fatal_error( "Error in reading iwftype")
        if (iwftype(1) .ne. 1) call fatal_error( 'INPUT: iwftype(1) ne 1')
    endif
    call bcast(iwftype)

    do i = 1, nforce
        do ic = 1, ncent
            if (wid) then
                read (iunit, *, iostat=iostat) delc(1, ic, i), delc(2, ic, i), delc(3, ic, i)
                if (iostat /= 0) call fatal_error( "Error in reading delc")
            endif
        enddo
    enddo
    call bcast(delc)


    do i = 1, nforce
      write(ounit,'(a,i4)') 'Number ::',i
      write(ounit,*) '-----------------------------------------------------------------------'
      write(ounit,'(a, t15, a, t27, a, t39, a, t45)') 'Symbol', 'x', 'y', 'z'
      write(ounit,'(t14, a, t26, a, t38, a )') '(A)', '(A)', '(A)'
      write(ounit,*) '-----------------------------------------------------------------------'
      do j= 1, ncent
          write(ounit,'(A4, 2x, 3F12.6)') symbol(j), (delc(k, j, i),k=1,3)
      enddo
    enddo

    iforces = 1
    call bcast(iforces)

    if (wid) close(iunit)
end subroutine read_forces_file

subroutine read_symmetry_file(file_symmetry)
    ! Ravindra

    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid, idtask

    use contrl_file,        only: ounit, errunit
    use coefs,              only: norb
    use optorb,             only: irrep
    use vmc_mod,            only: MORB

    implicit none

    !   local use
    character(len=72), intent(in)   :: file_symmetry
    character(len=40)               :: temp1, temp2, label
    integer                         :: iunit, iostat
    integer                         :: io, nsym, mo
    logical                         :: exist, skip = .true.


    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading
    write(ounit,*) '---------------------------------------------------------------------------'
    write(ounit,string_format)  " Reading orbital symmetries from the file :: ",  trim(file_symmetry)
    write(ounit,*) '---------------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_symmetry, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_symmetry, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error( "Problem in opening the symmetry file")
        else
            call fatal_error (" Orbital symmetries file "// trim(file_symmetry) // " does not exist.")
        endif
    endif

    if (wid) then
        read (iunit, *, iostat=iostat) label, nsym, mo
        if (iostat /= 0) call fatal_error( "Error in reading symmetry file :: expecting 'sym_labels', nsym, norb")
    endif
    call bcast(label)
    call bcast(nsym)
    call bcast(mo)

    if (trim(label) == "sym_labels") then
        if (norb /= mo) call fatal_error( "Number of orbitals not consistent with previous records")
    else
        call fatal_error (" Orbital symmetries file "// trim(file_symmetry) // " is corrupt.")
    endif


    ! Ignore irrep text labels
    if (wid) read (iunit, '(a80)') temp2
    write(ounit, *) "Irreducible representation correspondence for all norb orbitals"
    write(ounit, *) temp2

    ! safe allocate
    if (.not. allocated(irrep)) allocate (irrep(norb))

    ! read data
    if (wid) then
        read (iunit, *, iostat=iostat) (irrep(io), io=1, norb)
        if (iostat /= 0) call fatal_error( "Error in reading symmetry file :: expecting irrep correspondence for all norb orbitals")
    endif
    call bcast(irrep)

    write(ounit, '(10(1x, i3))') (irrep(io), io=1, norb)

    if (wid) close(iunit)
end subroutine read_symmetry_file


subroutine read_optorb_mixvirt_file(file_optorb_mixvirt)
    !
    ! Ravindra
    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use contrl_file,        only: ounit, errunit
    use optorb_mix,         only: iwmix_virt, norbopt, norbvirt
    use coefs,              only: norb
    use inputflags,         only: ioptorb_mixvirt

    implicit none

    !   local use
    character(len=72), intent(in)   :: file_optorb_mixvirt
    character(len=40)               :: temp1, temp2
    integer                         :: iunit, iostat, io, jo
    integer                         :: moopt, movirt
    logical                         :: exist, skip = .true.


    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading
    write(ounit,*) '---------------------------------------------------------------------------'
    write(ounit,string_format)  " Reading optorb_mixvirt from the file :: ",  trim(file_optorb_mixvirt)
    write(ounit,*) '---------------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_optorb_mixvirt, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_optorb_mixvirt, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error( "Problem in opening the optorb_mixvirt file")
        else
            call fatal_error (" optorb_mixvirt file "// trim(file_optorb_mixvirt) // " does not exist.")
        endif

        read (iunit, *, iostat=iostat) temp1, moopt, movirt
        if (iostat /= 0) call fatal_error( "Error in reading optorb_mixvirt file :: expecting 'optorb_mixvirt', norbopt, norbvirt")
    endif
    call bcast(temp1)

    if (trim(temp1) == "optorb_mixvirt") then
        if (moopt .gt. norb) call fatal_error( "Number of orbitals for optimization are greater than the total orbitals")
    else
        call fatal_error (" optorb_mixvirt file "// trim(file_optorb_mixvirt) // " is corrupt.")
    endif

    norbopt     =   moopt
    norbvirt    =   movirt

    call bcast(norbopt)
    call bcast(norbvirt)


    if (.not. allocated(iwmix_virt)) allocate (iwmix_virt(norbopt, norbvirt))

    do io = 1, norbopt
        if (wid) then
            read (iunit, *, iostat=iostat) (iwmix_virt(io, jo), jo=1, norbvirt)
            if (iostat /= 0) call fatal_error( "Error in reading optorb_mixvirt file :: incomplete data")
        endif
    enddo
    call bcast(iwmix_virt)

    write(ounit, *) "Printing which virtual orbitals are mixed with the occupied ones "
    do io = 1, norbopt
        write(ounit, '(10(1x, i5))') (iwmix_virt(io, jo), jo=1, norbvirt)
    enddo

    ioptorb_mixvirt = 1
    call bcast(ioptorb_mixvirt)

    if (wid) close(iunit)
end subroutine read_optorb_mixvirt_file



subroutine read_eigenvalues_file(file_eigenvalues)

    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use contrl_file,        only: ounit, errunit
    use coefs,              only: norb
    use vmc_mod,            only: MORB
    use optorb,             only: orb_energy

    implicit none

    !   local use
    character(len=72), intent(in)   :: file_eigenvalues
    character(len=40)               :: temp1, temp2
    integer                         :: iunit, iostat
    integer                         :: io, mo
    logical                         :: exist, skip = .true.


    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading
    write(ounit,*) '---------------------------------------------------------------------------'
    write(ounit,string_format)  " Reading orbital eigenvalues from the file :: ",  trim(file_eigenvalues)
    write(ounit,*) '---------------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_eigenvalues, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_eigenvalues, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error( "Problem in opening the eigenvalues file")
        else
            call fatal_error (" Orbital eigenvalues file "// trim(file_eigenvalues) // " does not exist.")
        endif

        read (iunit, *, iostat=iostat) temp1, mo
        if (iostat /= 0) call fatal_error( "Error in reading eigenvalues file :: expecting 'eigenvalues / energies', norb")
    endif
    call bcast(temp1)
    call bcast(mo)

    if ((trim(temp1) == "eigenvalues")  .or. (trim(temp1) == "energies")) then
        if (norb /= mo) call fatal_error( "Number of orbitals not consistent with previous records")
    else
        call fatal_error (" Orbital eigenvalues file "// trim(file_eigenvalues) // " is corrupt.")
    endif


    ! safe allocate
    if (.not. allocated(orb_energy)) allocate (orb_energy(norb))

    ! read data
    if (wid) then
        read (iunit, *, iostat=iostat) (orb_energy(io), io=1, norb)
        if (iostat /= 0) call fatal_error( "Error in reading eigenvalues file :: expecting eigenvalues of all norb orbitals")
    endif
    call bcast(orb_energy)

    write(ounit, *) "Eigenvalues of all orbitals"
    write(ounit, '(10(1x, i3))') (orb_energy(io), io=1, norb)

    if (wid) close(iunit)
end subroutine read_eigenvalues_file



subroutine read_basis_num_info_file(file_basis_num_info)
    ! Ravindra
    ! Basis function types and pointers to radial parts tables
    ! alternative name for keyword basis because of GAMBLE inputword basis because of GAMBLE input
    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use contrl_file,    only: ounit, errunit
    use numbas_mod, only: MRWF
    use vmc_mod, only: MCTYPE, MBASIS
    use numbas, only: iwrwf, numr
    use numbas1, only: iwlbas, nbastyp
    use basis, only: n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz
    use basis, only: n4s, n4p, n4fxxx, n4fyyy, n4fzzz, n4fxxy, n4fxxz, n4fyyx, n4fyyz
    use basis, only: n4fzzx, n4fzzy, n4fxyz, nsa, npa, ndzra, ndxya, ndxza, ndyza, ndx2a, ndz2a
    use inputflags, only: ibasis_num
    use coefs, only: nbasis
    use general, only: pooldir

    use atom, only: nctype
    use ghostatom, only: newghostype

    implicit none

    !   local use
    character(len=72), intent(in)   :: file_basis_num_info
    character(len=40)               :: temp1, temp2
    integer                         :: iunit, iostat
    integer                         :: i,j, jj, ib, nctot
    logical                         :: exist, skip = .true.


    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading
    write(ounit,*) '---------------------------------------------------------------------------'
    write(ounit,'(a)')  " Reading Basis function types and pointers to radial parts tables from the file :: ", &
                          trim(file_basis_num_info)
    write(ounit,*) '---------------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_basis_num_info, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_basis_num_info, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error( "Problem in opening the basis num info file")
        else
            call fatal_error (" Basis num info file "// trim(file_basis_num_info) // " does not exist.")
        endif

        read (iunit, *, iostat=iostat) temp1, numr
        if (iostat /= 0) call fatal_error( "Error in reading basis num info file :: expecting 'qmc_bf_info / basis', numr")

        if (.not. ((trim(temp1) == "qmc_bf_info")  .or. (trim(temp1) == "basis"))) then
            call fatal_error( "Error in reading basis num info file :: expecting 'qmc_bf_info / basis'")
        endif
    endif
    call bcast(numr)

    nctot = nctype + newghostype    ! DEBUG:: this statement might go. ghosttypes built-in

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

    if (nbasis .eq. 0) then
        call fatal_error('Please Load LCAO before basis info in the input file')
    endif

    allocate (iwlbas(nbasis, nctot))
    allocate (iwrwf(nbasis, nctot))

    if (wid) then
        do i = 1, nctype + newghostype
            read (iunit, *, iostat=iostat) n1s(i), n2s(i), (n2p(j, i), j=1, 3), &
                n3s(i), (n3p(j, i), j=1, 3), &
                n3dzr(i), n3dx2(i), n3dxy(i), n3dxz(i), n3dyz(i), &
                n4s(i), (n4p(j, i), j=1, 3), &
                n4fxxx(i), n4fyyy(i), n4fzzz(i), n4fxxy(i), n4fxxz(i), &
                n4fyyx(i), n4fyyz(i), n4fzzx(i), n4fzzy(i), n4fxyz(i), &
                nsa(i), (npa(j, i), j=1, 3), &
                ndzra(i), ndx2a(i), ndxya(i), ndxza(i), ndyza(i)
            if (iostat /= 0) call fatal_error( "Error in reading basis num info file")
            write (ounit, '(100i3)') n1s(i), n2s(i), (n2p(j, i), j=1, 3), &
                n3s(i), (n3p(j, i), j=1, 3), &
                n3dzr(i), n3dx2(i), n3dxy(i), n3dxz(i), n3dyz(i), &
                n4s(i), (n4p(j, i), j=1, 3), &
                n4fxxx(i), n4fyyy(i), n4fzzz(i), n4fxxy(i), n4fxxz(i), &
                n4fyyx(i), n4fyyz(i), n4fzzx(i), n4fzzy(i), n4fxyz(i), &
                nsa(i), (npa(j, i), j=1, 3), &
                ndzra(i), ndx2a(i), ndxya(i), ndxza(i), ndyza(i)


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

                read (iunit, *, iostat=iostat) (iwrwf(ib, i), ib=1, nbastyp(i))
                if (iostat /= 0) call fatal_error( "Error in reading basis num info file")
                write(ounit, '(100i3)') (iwrwf(ib, i), ib=1, nbastyp(i))
                write(ounit, *)

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

    endif

    call bcast(numr)
    call bcast(iwlbas)
    call bcast(iwrwf)
    call bcast(nbastyp)
    call bcast(n1s)
    call bcast(n2s)
    call bcast(n2p)
    call bcast(n3s)
    call bcast(n3p)
    call bcast(n3dzr)
    call bcast(n3dx2)
    call bcast(n3dxy)
    call bcast(n3dxz)
    call bcast(n3dyz)
    call bcast(n4s)
    call bcast(n4p)
    call bcast(n4fxxx)
    call bcast(n4fyyy)
    call bcast(n4fzzz)
    call bcast(n4fxxy)
    call bcast(n4fxxz)
    call bcast(n4fyyx)
    call bcast(n4fyyz)
    call bcast(n4fzzx)
    call bcast(n4fzzy)
    call bcast(n4fxyz)
    call bcast(nsa)
    call bcast(npa)
    call bcast(ndzra)
    call bcast(ndz2a)
    call bcast(ndxya)
    call bcast(ndxza)
    call bcast(ndx2a)
    call bcast(ndyza)

    ibasis_num = 1
    call bcast(ibasis_num)
    if (wid) close(iunit)
end subroutine read_basis_num_info_file


subroutine read_dmatrix_file(file_dmatrix)
    ! Ravindra (no=ndetorb, ns=nweight)
    !INPUT dmatrix i i a=<input>
    !KEYDOC Read diagonal density matrix information.
    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use contrl_file,    only: ounit, errunit
    use precision_kinds, only: dp
    use vmc_mod, only: MORB
    use csfs, only: nstates
    use sa_weights, only: iweight, nweight, weights
    use mstates_mod, only: MSTATES
    use coefs, only: norb
    use optorb, only: dmat_diag

    implicit none

    !   local use
    character(len=72), intent(in)   :: file_dmatrix
    character(len=40)               :: temp1, temp2
    integer                         :: iunit, iostat
    integer                         :: i,j, iw, ndetorb, ipr
    logical                         :: exist, skip = .true.

    real(dp), dimension(:), allocatable :: dmat
    integer, dimension(:), allocatable  :: iwdmat

    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading
    write(ounit,*) '---------------------------------------------------------------------------'
    write(ounit,string_format)  " Reading dmatrix the file :: ",  trim(file_dmatrix)
    write(ounit,*) '---------------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_dmatrix, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_dmatrix, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error( "Problem in opening the dmatrix file")
        else
            call fatal_error (" dmatrix file "// trim(file_dmatrix) // " does not exist.")
        endif

        read (iunit, *, iostat=iostat) temp1, ndetorb, nweight
        if (iostat /= 0) call fatal_error( "Error in reading dmatrix file :: expecting 'dmatrix', ndetorb, nweight")

        if (.not. (trim(temp1) == "dmatrix") ) then
            call fatal_error( "Error in reading dmatrix file :: expecting 'dmatrix'")
        endif
    endif

    call bcast(ndetorb)
    call bcast(nweight)


    allocate (dmat(norb))
    allocate (iwdmat(nstates))

    if (ndetorb .gt. norb) call fatal_error( 'READ_DMATRIX: wrong number of orbitals')

    allocate (weights(nstates))
    allocate (iweight(nstates))


    if (wid) read (iunit, *) (iwdmat(i), i=1, nweight)
    call bcast(iwdmat)

    do iw = 1, nweight
        if (iwdmat(iw) .ne. iweight(iw)) call fatal_error('READ_DMATRIX: iwdmat')
    enddo

    allocate (dmat_diag(norb))
    dmat_diag = 0.0d0

    do iw = 1, nweight
        if (wid) read (iunit, *, iostat=iostat) (dmat(j), j=1, ndetorb)
        call bcast(dmat)
        do j = 1, ndetorb
            dmat_diag(j) = dmat_diag(j) + weights(iw)*dmat(j)
        enddo
    enddo
    call bcast(dmat_diag)

    !DEBUG why
    do i = 1, ndetorb
        if (dabs(dmat_diag(i) - 1.d0) .lt. 1.d-6) dmat_diag(i) = 1.d0
    enddo

    if (ipr .gt. 2) then
        write(ounit, '(''diagonal elements of the density matrix'')')
        write(ounit, '(100f10.6)') (dmat_diag(i), i=1, ndetorb)
    endif

    deallocate (dmat)
    deallocate (iwdmat)
    if (wid) close(iunit)

end subroutine read_dmatrix_file


subroutine read_cavity_spheres_file(file_cavity_spheres)
    ! Ravindra
    ! Read centers of cavity spheres and radii
    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use contrl_file,        only: ounit, errunit
    use pcm_parms,          only: nesph, re, re2
    use pcm_parms,          only: xe, ye, ze

    implicit none

    !   local use
    character(len=72), intent(in)   :: file_cavity_spheres
    character(len=40)               :: key
    integer                         :: iunit, iostat
    integer                         :: i,j
    logical                         :: exist, skip = .true.

    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading
    write(ounit,*) '---------------------------------------------------------------------------'
    write(ounit,string_format)  " Reading cavity spheres from the file :: ",  trim(file_cavity_spheres)
    write(ounit,*) '---------------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_cavity_spheres, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_cavity_spheres, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error( "Problem in opening the cavity spheres file")
        else
            call fatal_error (" cavity spheres file "// trim(file_cavity_spheres) // " does not exist.")
        endif


        read (iunit, *, iostat=iostat) key, nesph
        if (iostat /= 0) call fatal_error( "Error in reading cavity spheres file :: expecting 'cavity_spheres', nspheres")

        if (.not. (trim(key) == "cavity_spheres") ) then
            error stop "Error in reading cavity_spheres file :: expecting 'cavity_spheres'"
        endif

    endif
    call bcast(nesph)

    if (.not. allocated(re)) allocate (re(nesph))
    if (.not. allocated(re2)) allocate (re2(nesph))
    if (.not. allocated(xe)) allocate (xe(nesph))
    if (.not. allocated(ye)) allocate (ye(nesph))
    if (.not. allocated(ze)) allocate (ze(nesph))

    do i = 1, nesph
        if (wid) read (iunit, *) xe(i), ye(i), ze(i), re(i)
        re2(i) = re(i)*re(i)
    enddo
    call bcast(xe)
    call bcast(ye)
    call bcast(ze)
    call bcast(re)
    call bcast(re2)

    if (wid) close(iunit)
end subroutine read_cavity_spheres_file


subroutine read_gradients_cartesian_file(file_gradients_cartesian)
    !Ravindra
    !INPUT gradients_cartesian inp
    !KEYDOC Read for which x,y,z cartesian coordiantes of
    !KEYDOC atoms energy gradients are to be calculated for.

    !     Originally written by Omar Valsson
    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use contrl_file,        only: ounit, errunit
    use vmc_mod,            only: MCENT
    use forcepar,           only: nforce
    use force_mod,          only: MFORCE
    use forcestr,           only: delc
    use grdntsmv,           only: igrdaidx, igrdcidx, igrdmv
    use grdntspar,          only: delgrdxyz, igrdtype, ngradnts
    use wfsec,              only: iwftype
    use inputflags,         only: igradients

    use atom,               only: ncent

    implicit none

!   local use
    character(len=72), intent(in)   :: file_gradients_cartesian
    character(len=40)               :: key
    integer                         :: iunit, iostat
    integer                         :: i,ia, ic, k
    logical                         :: exist, skip = .true.

    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading
    write(ounit,*) '---------------------------------------------------------------------------'
    write(ounit,string_format)  " Reading gradients cartesian from the file :: ",  trim(file_gradients_cartesian)
    write(ounit,*) '---------------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_gradients_cartesian, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_gradients_cartesian, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error( "Problem in opening the gradients_cartesian file")
        else
            call fatal_error ( " Gradients cartesian file "// trim(file_gradients_cartesian) // " does not exist.")
        endif

        read (iunit, *, iostat=iostat) key
        if (iostat /= 0) call fatal_error( "Error in reading gradients cartesian file :: expecting 'gradients_cartesian'")

        if (.not. (trim(key) == "gradients_cartesian") ) then
            error stop "Error in reading gradients cartesian file :: expecting 'gradients_cartesian'"
        endif

    endif

    if (igrdtype .ne. 1) call fatal_error('GRADIENTS_CARTESIAN: igrdtype /= 1')
    if ((2*ngradnts + 1) .ne. nforce) call fatal_error('GRADIENTS_CARTESIAN: (2*ngradnts+1)  /=  nforce')

    if (.not. allocated(delc)) allocate (delc(3, ncent, MFORCE))
    if (.not. allocated(igrdaidx)) allocate (igrdaidx(MFORCE))
    if (.not. allocated(igrdcidx)) allocate (igrdcidx(MFORCE))
    if (.not. allocated(igrdmv)) allocate (igrdmv(3, ncent))

    ! initialize the values to zero

    iwftype = 1
    igrdmv  = 0
    delc    = 0.0d0

    ia = 2
    if (wid) then
        do ic = 1, ncent
            read (iunit, *, iostat=iostat) (igrdmv(k, ic), k=1, 3)
            do k = 1, 3
                if (igrdmv(k, ic) .lt. 0 .or. igrdmv(k, ic) .gt. 1) then
                    call fatal_error('GRADIENTS_CARTESIAN: igrdmv \= 0,1')
                endif
                if (igrdmv(k, ic) .eq. 1) then
                    igrdaidx(ia/2) = ic
                    igrdcidx(ia/2) = k
                    delc(k, ic, ia) = delgrdxyz
                    delc(k, ic, ia + 1) = -delgrdxyz
                    ia = ia + 2
                endif
            enddo
        enddo
    endif

    call bcast(igrdmv)
    call bcast(igrdaidx)
    call bcast(igrdcidx)
    call bcast(delc)

    igradients = 1
    call bcast(igradients)

    if (wid) close(iunit)
end subroutine read_gradients_cartesian_file

subroutine read_gradients_zmatrix_file(file_gradients_zmatrix)
    ! Ravindra
    ! Read for which Z matrix (internal) coordiantes of
    ! atoms energy gradients are to be calculated for.

    ! Originally written by Omar Valsson.
    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use contrl_file,    only: ounit, errunit
    use vmc_mod, only: MCENT
    use forcepar, only: nforce
    use force_mod, only: MFORCE
    use forcestr, only: delc
    use grdntsmv, only: igrdaidx, igrdcidx, igrdmv
    use grdntspar, only: delgrdba, delgrdbl, delgrdda, igrdtype, ngradnts
    use zmatrix, only: izmatrix
    use wfsec, only: iwftype
    use inputflags, only: igradients

    use atom, only: ncent

    implicit none

!   local use
    character(len=72), intent(in)   :: file_gradients_zmatrix
    character(len=40)               :: key
    integer                         :: iunit, iostat
    integer                         :: i,ia, ic, k
    logical                         :: exist, skip = .true.

    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading
    write(ounit,*) '---------------------------------------------------------------------------'
    write(ounit,string_format)  " Reading gradients zmatrix from the file :: ",  trim(file_gradients_zmatrix)
    write(ounit,*) '---------------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_gradients_zmatrix, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_gradients_zmatrix, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error( "Problem in opening the gradients_zmatrix file")
        else
            call fatal_error ( " Gradients zmatrix file "// trim(file_gradients_zmatrix) // " does not exist.")
        endif

        read (iunit, *, iostat=iostat) key
        if (iostat /= 0) call fatal_error( "Error in reading gradients zmatrix file :: expecting 'gradients_zmatrix'")

        if (.not. (trim(key) == "gradients_zmatrix") ) then
            call fatal_error( "Error in reading gradients zmatrix file :: expecting 'gradients_zmatrix'")
        endif
    endif

    if (igrdtype .ne. 2) call fatal_error('GRADIENTS_ZMATRIX: igrdtype /= 2')
    if (izmatrix .ne. 1) call fatal_error('GRADIENTS_ZMATRIX: No Z matrix connection matrix')
    if ((2*ngradnts + 1) .ne. nforce) call fatal_error('GRADIENTS_ZMATRIX: (2*ngradnts+1)  /=  nforce')

    if (.not. allocated(delc)) allocate (delc(3, ncent, MFORCE))
    if (.not. allocated(igrdaidx)) allocate (igrdaidx(MFORCE))
    if (.not. allocated(igrdcidx)) allocate (igrdcidx(MFORCE))
    if (.not. allocated(igrdmv)) allocate (igrdmv(3, ncent))

    ! initialize the values to zero

    iwftype = 1
    igrdmv  = 0
    delc    = 0.0d0


    ia = 2
    if (wid) then
        do ic = 1, ncent
            read (iunit, *, iostat=iostat) (igrdmv(k, ic), k=1, 3)
            do k = 1, 3
                if (igrdmv(k, ic) .lt. 0 .or. igrdmv(k, ic) .gt. 1) call fatal_error('GRADIENTS_ZMATRIX: igrdmv \= 0,1')
                if (igrdmv(k, ic) .eq. 1) then
                    igrdaidx(ia/2) = ic
                    igrdcidx(ia/2) = k
                    call grdzmat_displ(k, ic, ia, +1.0d0)
                    call grdzmat_displ(k, ic, ia + 1, -1.0d0)
                    ia = ia + 2
                endif
            enddo
        enddo
    endif
    call bcast(igrdmv)
    call bcast(igrdaidx)
    call bcast(igrdcidx)
    call bcast(delc)

    igradients = 1
    call bcast(igradients)

    if (wid) close(iunit)

end subroutine read_gradients_zmatrix_file

subroutine read_modify_zmatrix_file(file_modify_zmatrix)
    ! Ravindra
    !
    ! Read for which Z matrix (internal) coordiantes of
    ! atoms energy gradients are to be calculated for.
    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use contrl_file,        only: ounit, errunit
    use grdntsmv,           only: igrdmv
    use inputflags,         only: imodify_zmat

    use atom,               only: ncent

    implicit none

!   local use
    character(len=72), intent(in)   :: file_modify_zmatrix
    character(len=40)               :: key
    integer                         :: iunit, iostat
    integer                         :: ic,k
    logical                         :: exist, skip = .true.

    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading
    write(ounit,*) '---------------------------------------------------------------------------'
    write(ounit,string_format)  " Reading modify zmatrix from the file :: ",  trim(file_modify_zmatrix)
    write(ounit,*) '---------------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_modify_zmatrix, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_modify_zmatrix, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error( "Problem in opening the modify_zmatrix file")
        else
            call fatal_error (" modify zmatrix file "// trim(file_modify_zmatrix) // " does not exist.")
        endif


        read (iunit, *, iostat=iostat) key
        if (iostat /= 0) call fatal_error( "Error in reading modify zmatrix file")

        if (.not. (trim(key) == "modify_zmatrix") ) then
            call fatal_error( "Error in reading modify zmatrix file :: expecting 'modify_zmatrix'")
        endif

    endif


    if (.not. allocated(igrdmv)) allocate (igrdmv(3, ncent))

    if (wid) then
        do ic = 1, ncent
            read (iunit, *, iostat=iostat) (igrdmv(k, ic), k=1, 3)
            do k = 1, 3
                if (igrdmv(k, ic) .lt. 0 .or. igrdmv(k, ic) .gt. 1) then
                    call fatal_error('MODIFY_ZMATRIX: igrdmv \= 0,1')
                endif
            enddo
        enddo
    endif
    call bcast(igrdmv)

    imodify_zmat = 1
    call bcast(imodify_zmat)

    if (wid) close(iunit)
end subroutine read_modify_zmatrix_file

subroutine read_hessian_zmatrix_file(file_hessian_zmatrix)
    ! Ravindra
    !
    ! Read for which Z matrix (internal) coordiantes of
    ! atoms energy gradients are to be calculated for.
    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use contrl_file,        only: ounit, errunit
    use grdnthes,           only: hessian_zmat
    use inputflags,         only: ihessian_zmat
    use atom,               only: ncent

    implicit none

!   local use
    character(len=72), intent(in)   :: file_hessian_zmatrix
    character(len=40)               :: key
    integer                         :: iunit, iostat
    integer                         :: ic,k
    logical                         :: exist, skip = .true.

    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading
    write(ounit,*) '---------------------------------------------------------------------------'
    write(ounit,string_format)  " Reading hessian zmatrix from the file :: ",  trim(file_hessian_zmatrix)
    write(ounit,*) '---------------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_hessian_zmatrix, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_hessian_zmatrix, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error( "Problem in opening the hessian_zmatrix file")
        else
            call fatal_error (" hessian zmatrix file "// trim(file_hessian_zmatrix) // " does not exist.")
        endif


        read (iunit, *, iostat=iostat) key
        if (iostat /= 0) call fatal_error( "Error in reading hessian zmatrix file")

        if (.not. (trim(key) == "hessian_zmatrix") ) then
            call fatal_error( "Error in reading hessian zmatrix file :: expecting 'hessian_zmatrix'")
        endif

    endif

    if (.not. allocated(hessian_zmat)) allocate (hessian_zmat(3, ncent))

    if (wid) then
        do ic = 1, ncent
            read (iunit, *) (hessian_zmat(k, ic), k=1, 3)
            do k = 1, 3
                if (hessian_zmat(k, ic) .le. 0) then
                    call fatal_error('HESSIAN_ZMATRIX: hess <=  0')
                endif
            enddo
        enddo
    endif

    call bcast(hessian_zmat)

    ihessian_zmat = 1
    call bcast(ihessian_zmat)

    if (wid) close(iunit)
end subroutine read_hessian_zmatrix_file


subroutine read_zmatrix_connection_file(file_zmatrix_connection)
    ! Ravindra
    !
    ! Read the atom connection matrix for the Z matrix.
    ! It is need when calculating forces in Z matrix
    ! coordinates.

    ! Originally written by Omar Valsson
    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use contrl_file,    only: ounit, errunit
    use atom, only: cent, ncent
    use zmatrix, only: czcart, czint, czcart_ref, izcmat, izmatrix
    use inputflags, only: izmatrix_check

    implicit none

!   local use
    character(len=72), intent(in)   :: file_zmatrix_connection
    character(len=40)               :: key
    integer                         :: iunit, iostat
    integer                         :: k, ic
    logical                         :: exist, skip = .true.

    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading
    write(ounit,*) '---------------------------------------------------------------------------'
    write(ounit,string_format)  " Reading zmatrix connection matrix from the file :: ",  trim(file_zmatrix_connection)
    write(ounit,*) '---------------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_zmatrix_connection, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_zmatrix_connection, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error( "Problem in opening the zmatrix connection matrix file")
        else
            call fatal_error (" zmatrix connection matrix file "// trim(file_zmatrix_connection) // " does not exist.")
        endif


        read (iunit, *, iostat=iostat) key
        if (iostat /= 0) call fatal_error( "Error in reading zmatrix connection matrix file")

        if (.not. (trim(key) == "zmatrix_connectionmatrix") ) then
            call fatal_error( "Error in reading zmatrix connection matrix file :: expecting 'zmatrix_connectionmatrix'")
        endif

    endif


    if (.not. allocated(czcart)) allocate (czcart(3, ncent))
    if (.not. allocated(czint)) allocate (czint(3, ncent))
    if (.not. allocated(czcart_ref)) allocate (czcart_ref(3, 3))
    if (.not. allocated(izcmat)) allocate (izcmat(3, ncent))

    czcart_ref = cent

    izcmat     = 0
    czint      = 0.0d0
    czcart     = cent

    if (wid) then
        do ic = 1, ncent
            read (iunit, *, iostat=iostat) (izcmat(k, ic), k=1, 3)
            do k = 1, 3
                if (izcmat(k, ic) .ge. ic) call fatal_error('ZMATRIX: Error in connection matrix')
            enddo
        enddo
    endif
    call bcast(izcmat)

    call cart2zmat(ncent, czcart, izcmat, czint)
    call zmat2cart_rc(ncent, izcmat, czint, czcart, czcart_ref)

    call bcast(czcart)
    call bcast(czint)
    call bcast(czcart_ref)

    izmatrix = 1
    izmatrix_check = 1
    call bcast(izmatrix)
    call bcast(izmatrix_check)

    if (wid) close(iunit)
end subroutine read_zmatrix_connection_file

subroutine read_efield_file(file_efield) !ncharges_tmp, iscreen_tmp
    ! Ravindra
    use custom_broadcast,   only: bcast
    use mpiconf,            only: wid

    use contrl_file,    only: ounit, errunit
    use efield_mod, only: MCHARGES
    use efield_blk, only: ascreen, bscreen, qcharge, xcharge, ycharge, zcharge
    use efield, only: iscreen, ncharges
    use inputflags, only: icharge_efield

    implicit none

!   local use
    character(len=72), intent(in)   :: file_efield
    character(len=40)               :: key
    integer                         :: iunit, iostat
    integer                         :: ncharges_tmp, iscreen_tmp, i
    logical                         :: exist, skip = .true.

    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading
    write(ounit,*) '---------------------------------------------------------------------------'
    write(ounit,string_format)  " Reading efield from the file :: ",  trim(file_efield)
    write(ounit,*) '---------------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_efield, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_efield, iostat=iostat, action='read' )
            if (iostat .ne. 0) error stop "Problem in opening the efield file"
        else
            call fatal_error ( " efield file "// trim(file_efield) // " does not exist.")
        endif


        read (iunit, *, iostat=iostat) key, ncharges_tmp, iscreen_tmp
        if (iostat /= 0) error stop "Error in reading efield file"

        if (.not. (trim(key) == "efield") ) then
            error stop "Error in reading efield file :: expecting 'efield'"
        endif

    endif

!    call file(iu, filename, 'old', 1, 0)  <-- whats is this?
    ncharges = ncharges_tmp
    iscreen  = iscreen_tmp
    write(ounit, *) 'reading in', ncharges, ' charges!'

    if (ncharges .gt. MCHARGES) call fatal_error('EFIELD: ncharges > MCHARGES')

    if (.not. allocated(ascreen)) allocate (ascreen(ncharges))
    if (.not. allocated(bscreen)) allocate (bscreen(ncharges))
    if (.not. allocated(qcharge)) allocate (qcharge(ncharges))
    if (.not. allocated(xcharge)) allocate (xcharge(ncharges))
    if (.not. allocated(ycharge)) allocate (ycharge(ncharges))
    if (.not. allocated(zcharge)) allocate (zcharge(ncharges))

    if (wid) then
        do i = 1, ncharges
            read (iunit, *, iostat=iostat) xcharge(i), ycharge(i), zcharge(i), qcharge(i), ascreen(i), bscreen(i)
        enddo
    endif
    call bcast(xcharge)
    call bcast(ycharge)
    call bcast(zcharge)
    call bcast(qcharge)
    call bcast(ascreen)
    call bcast(bscreen)

    icharge_efield = icharge_efield + 1
    write(ounit, *) 'icharge_efield=', icharge_efield

    call bcast(icharge_efield)

    if (wid) close(iunit)
end subroutine read_efield_file