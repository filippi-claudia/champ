module trexio_read_data
    use error, only : fatal_error
    contains

    subroutine read_trexio_molecule_file(file_trexio)
        !> This subroutine reads the .hdf5 trexio generated file/folder. It then computes the
        !! number of types of atoms, nuclear charges (from the symbol), and
        !! number of valence electrons if pseudopotential is provided.
        !! @author Ravindra Shinde (r.l.shinde@utwente.nl)
        !! @date 07 October 2021
        use custom_broadcast,   only: bcast
        use mpiconf,            only: wid
        use atom,               only: znuc, cent, pecent, iwctype, nctype, ncent, ncent_tot, nctype_tot, symbol, atomtyp
        use ghostatom, 		    only: newghostype, nghostcent
        use inputflags,         only: igeometry
        use periodic_table,     only: atom_t, element
        use elec,           	only: ndn, nup
        use const,          	only: nelec
        use contrl_file,        only: ounit, errunit
        use general,            only: pooldir
        use precision_kinds, only: dp
#if defined(TREXIO_FOUND)
        use trexio
        use contrl_file,        only: backend
#endif

        implicit none

        !   local use
        character(len=72), intent(in)   :: file_trexio
        character(len=40)               :: temp1, temp2, temp3, temp4
        character(len=80)               :: comment, file_trexio_path
        integer                         :: iostat, i, j, k, iunit
        logical                         :: exist
        type(atom_t)                    :: atoms
        character(len=2), allocatable   :: unique(:)

        ! trexio
        integer(8)                      :: trex_molecule_file
        integer                         :: rc = 1

        !   Formatting
        character(len=100)              :: int_format     = '(A, T60, I0)'
        character(len=100)              :: float_format   = '(A, T60, f12.8)'
        character(len=100)              :: string_format  = '(A, T60, A)'

        trex_molecule_file = 0

        !   External file reading

        if((file_trexio(1:6) == '$pool/') .or. (file_trexio(1:6) == '$POOL/')) then
            file_trexio_path = pooldir // file_trexio(7:)
        else
            file_trexio_path = file_trexio
        endif

        write(ounit,*) '-----------------------------------------------------------------------'
        write(ounit,string_format)  " Reading molecular coordinates from the trexio file :: ",  file_trexio_path
        write(ounit,*) '-----------------------------------------------------------------------'

        ! Check if the file exists
        if (wid) then
#if defined(TREXIO_FOUND)
            trex_molecule_file = trexio_open(file_trexio_path, 'r', backend, rc)
            rc = trexio_read_nucleus_num(trex_molecule_file, ncent)
            rc = trexio_read_electron_up_num(trex_molecule_file, nup)
            rc = trexio_read_electron_dn_num(trex_molecule_file, ndn)
#endif
        endif
        call bcast(ncent)
        call bcast(nup)
        call bcast(ndn)

        nelec = nup + ndn

        ! Do the allocations based on the ncent
        if (.not. allocated(cent))    allocate(cent(3,ncent))
        if (.not. allocated(symbol))  allocate(symbol(ncent))
        if (.not. allocated(iwctype)) allocate(iwctype(ncent))
        if (.not. allocated(unique))  allocate(unique(ncent))

        if (wid) then
#if defined(TREXIO_FOUND)
        rc = trexio_read_nucleus_coord(trex_molecule_file, cent)
        rc = trexio_read_nucleus_label(trex_molecule_file, symbol, 3)
        rc = trexio_close(trex_molecule_file)
#endif
        endif
        call bcast(cent)
        call bcast(symbol)


        write(ounit,fmt=int_format) " Number of atoms ::  ", ncent
        write(ounit,*)

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
        write(ounit,'(t14, a, t26, a, t38, a )') '(bohr)', '(bohr)', '(bohr)'
        write(ounit,*) '-----------------------------------------------------------------------'

        do j= 1, ncent
            write(ounit,'(A4, 2x, 3F12.6, 2x, i3)') symbol(j), (cent(i,j),i=1,3), iwctype(j)
        enddo

        write(ounit,*) '-----------------------------------------------------------------------'
        write(ounit,*) " Values of znuc (number of valence electrons) "
        write(ounit,'(10F12.6)') (znuc(j), j = 1, nctype)
        write(ounit,*) '-----------------------------------------------------------------------'
        write(ounit,*)
    end subroutine read_trexio_molecule_file


    subroutine read_trexio_orbitals_file(file_trexio)
        !> This subroutine reads the .hdf5 trexio generated file/folder. It then reads the
        !! number of molecular and atomic orbitals and their corresponding coefficients.
        !! @author Ravindra Shinde (r.l.shinde@utwente.nl)
        !! @date 12 October 2021
        use custom_broadcast,   only: bcast
        use mpiconf,            only: wid
        use contrl_file,        only: ounit, errunit
        use coefs,              only: coef, nbasis, norb
        use inputflags,         only: ilcao
        use orbval,             only: nadorb
        use pcm_fdc,            only: fs
        use vmc_mod,            only: norb_tot
        use wfsec,              only: nwftype
        use general,            only: pooldir
        use method_opt,         only: method
        use precision_kinds, only: dp
#if defined(TREXIO_FOUND)
        use trexio
        use contrl_file,        only: backend
#endif
        implicit none

    !   local use
        character(len=72), intent(in)   :: file_trexio
        character(len=40)               :: temp1, temp2
        character(len=120)              :: temp3, file_trexio_path
        integer                         :: iunit, iostat, iwft
        integer                         :: iorb, ibasis, i, k, counter
        logical                         :: exist
        logical                         :: skip = .true.

        !   Formatting
        character(len=100)               :: int_format     = '(A, T60, I0)'
        character(len=100)               :: string_format  = '(A, T60, A)'
        character(len=100)               :: float_format   = '(A, T60, f12.8)'

        ! trexio
        integer(8)                      :: trex_orbitals_file
        integer                         :: rc = 1

        iwft = 0
        trex_orbitals_file = 0
        !   External file reading

        if((file_trexio(1:6) == '$pool/') .or. (file_trexio(1:6) == '$POOL/')) then
            file_trexio_path = pooldir // file_trexio(7:)
        else
            file_trexio_path = file_trexio
        endif

        write(ounit,*) '---------------------------------------------------------------------------'
        write(ounit,string_format)  " Reading LCAO orbitals from the file :: ",  trim(file_trexio_path)
        write(ounit,*) '---------------------------------------------------------------------------'
        ! Check if the file exists

        if (wid) then
#if defined(TREXIO_FOUND)
            trex_orbitals_file = trexio_open(file_trexio_path, 'r', backend, rc)
            rc = trexio_read_mo_num(trex_orbitals_file, norb)
            rc = trexio_read_ao_num(trex_orbitals_file, nbasis)
#endif
        endif
        call bcast(norb)
        call bcast(nbasis)

        ! Do the array allocations
        if( (method(1:3) == 'lin')) then
            if (.not. allocated(coef)) allocate (coef(nbasis, norb, 3))
        else
            if (.not. allocated(coef)) allocate (coef(nbasis, norb, nwftype))
        endif

        ! Read the orbitals
        if (wid) then
#if defined(TREXIO_FOUND)
            rc = trexio_read_mo_coefficient(trex_orbitals_file, coef(:,:,1))
#endif
        endif
        call bcast(coef)

        ! IMPORTANT
        ! The orbital ordering should be made consistent with the CHAMP code.
        ! call a function to transform the ordering.
        ! DEBUG #148
        ! Close the trexio file
#if defined(TREXIO_FOUND)
        if (wid) rc = trexio_close(trex_orbitals_file)
#endif


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


        write(ounit,*) "----------------------------------------------------------"

    end subroutine read_trexio_orbitals_file



end module