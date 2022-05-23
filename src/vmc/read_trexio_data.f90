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


    subroutine read_trexio_basis_file(file_trexio)
        !> This subroutine reads the .hdf5 trexio generated file/folder.
        !! It reads the exponents, coefficients, number of basis functions,
        !! shell angular momentum, and number of shells.
        !! @author Ravindra Shinde (r.l.shinde@utwente.nl)
        !! @date 23 May 2022
        use custom_broadcast,   only: bcast
        use mpiconf,            only: wid
        use periodic_table,     only: atom_t, element
        use contrl_file,        only: ounit, errunit
        use general,            only: pooldir
        use precision_kinds, only: dp

#if defined(TREXIO_FOUND)
        use trexio
        use contrl_file,        only: backend
#endif

        implicit none

        !for local use.  To be read from trexio file
        integer                         :: basis_num_shell
        integer                         :: basis_num_prim
        integer, allocatable            :: basis_nucleus_index(:)
        integer, allocatable            :: basis_shell_index(:)
        integer, allocatable            :: basis_shell_ang_mom(:)
        real(dp), allocatable           :: basis_shell_factor(:)
        real(dp), allocatable           :: basis_exponent(:)
        real(dp), allocatable           :: basis_coefficient(:)
        real(dp), allocatable           :: basis_prim_factor(:)

        integer                         :: ao_num
        integer,allocatable             :: ao_shell(:)
        real(dp),allocatable            :: ao_normalization(:)

        logical                         :: ao_cartesian

        ! for local use.
        character(len=72), intent(in)   :: file_trexio
        character(len=128)              :: file_trexio_path
        integer                         :: iostat, i, j, k, iunit
        logical                         :: exist
        type(atom_t)                    :: atoms
        character(len=2), allocatable   :: unique(:)

        ! trexio
        integer(8)                      :: trex_basis_file
        integer                         :: rc = 1

        !   Formatting
        character(len=128)              :: int_format     = '(A, T60, I0)'
        character(len=128)              :: float_format   = '(A, T60, f12.8)'
        character(len=128)              :: string_format  = '(A, T60, A128)'


        trex_basis_file = 0

        !   External file reading

        if((file_trexio(1:6) == '$pool/') .or. (file_trexio(1:6) == '$POOL/')) then
            file_trexio_path = pooldir // file_trexio(7:)
        else
            file_trexio_path = file_trexio
        endif

        write(ounit,*) '-----------------------------------------------------------------------'
        write(ounit,*) " Reading Basis Set information from the trexio file :: ", len(file_trexio_path), trim(adjustl(file_trexio_path))
        write(ounit,*) '-----------------------------------------------------------------------'

        ! Check if the file exists
        if (wid) then
#if defined(TREXIO_FOUND)
            trex_basis_file = trexio_open(file_trexio_path, 'r', backend, rc)
            ! write(*,*) "trexio_open :: ", rc
            rc = trexio_read_basis_prim_num(trex_basis_file, basis_num_prim)
            ! write(*,*) "trexio_read_basis_prim_num :: ", rc
            rc = trexio_read_basis_shell_num(trex_basis_file, basis_num_shell)
            ! write(*,*) "trexio_read_basis_shell_num :: ", rc
            rc = trexio_read_ao_num(trex_basis_file, ao_num)
            ! write(*,*) "trexio_read_ao_num :: ", rc
#endif
        endif
        call bcast(basis_num_prim)
        call bcast(basis_num_shell)
        call bcast(ao_num)

        ! Do the allocations based on the number of shells and primitives
        if (.not. allocated(basis_nucleus_index))    allocate(basis_nucleus_index(basis_num_shell))
        if (.not. allocated(basis_shell_index))      allocate(basis_shell_index(basis_num_prim))
        if (.not. allocated(basis_shell_ang_mom))    allocate(basis_shell_ang_mom(basis_num_shell))
        if (.not. allocated(basis_shell_factor))     allocate(basis_shell_factor(basis_num_shell))
        if (.not. allocated(basis_exponent))         allocate(basis_exponent(basis_num_prim))
        if (.not. allocated(basis_coefficient))      allocate(basis_coefficient(basis_num_prim))
        if (.not. allocated(basis_prim_factor))      allocate(basis_prim_factor(basis_num_prim))
        if (.not. allocated(ao_shell))               allocate(ao_shell(ao_num))
        if (.not. allocated(ao_normalization))       allocate(ao_normalization(ao_num))

        if (wid) then
#if defined(TREXIO_FOUND)
            trex_basis_file = trexio_open(file_trexio_path, 'r', backend, rc)
            ! write(*,*) "trexio_open :: ", rc
            rc = trexio_read_basis_nucleus_index(trex_basis_file, basis_nucleus_index)
            ! write(*,*) "trexio_read_basis_nucleus_index :: ", rc
            rc = trexio_read_basis_shell_index(trex_basis_file, basis_shell_index)
            ! write(*,*) "trexio_read_basis_shell_index :: ", rc
            rc = trexio_read_basis_shell_ang_mom(trex_basis_file, basis_shell_ang_mom)
            ! write(*,*) "trexio_read_basis_shell_ang_mom :: ", rc
            rc = trexio_read_basis_shell_factor(trex_basis_file, basis_shell_factor)
            ! write(*,*) "trexio_read_basis_shell_factor :: ", rc
            rc = trexio_read_basis_exponent(trex_basis_file, basis_exponent)
            ! write(*,*) "trexio_read_basis_exponent :: ", rc
            rc = trexio_read_basis_coefficient(trex_basis_file, basis_coefficient)
            ! write(*,*) "trexio_read_basis_coefficient :: ", rc
            rc = trexio_read_basis_prim_factor(trex_basis_file, basis_prim_factor)
            ! write(*,*) "trexio_read_basis_prim_factor :: ", rc
            rc = trexio_read_ao_shell(trex_basis_file, ao_shell)
            ! write(*,*) "trexio_read_ao_shell :: ", rc
            rc = trexio_read_ao_normalization(trex_basis_file, ao_normalization)
            ! write(*,*) "trexio_read_ao_normalization :: ", rc
#endif
        endif
        call bcast(basis_nucleus_index)
        call bcast(basis_shell_index)
        call bcast(basis_shell_ang_mom)
        call bcast(basis_shell_factor)
        call bcast(basis_exponent)
        call bcast(basis_coefficient)
        call bcast(basis_prim_factor)
        call bcast(ao_shell)
        call bcast(ao_normalization)



        write(ounit,fmt=int_format) " Number of primitives ::  ", basis_num_prim
        write(ounit,fmt=int_format) " Number of shells     ::  ", basis_num_shell
        write(ounit,fmt=int_format) " Number of AO         ::  ", ao_num
        write(ounit,*) " Nucleus Index      ::  ", (basis_nucleus_index(i), i=1, basis_num_shell)
        write(ounit,*) " Shell Index        ::  ", (basis_shell_index(i), i=1, basis_num_prim)
        write(ounit,*) " Shell Ang. Mom.    ::  ", (basis_shell_ang_mom(i), i=1, basis_num_shell)
        write(ounit,*) " Shell Factor       ::  ", (basis_shell_factor(i), i=1, basis_num_shell)
        write(ounit,*) " Exponent           ::  ", (basis_exponent(i), i=1, basis_num_prim)
        write(ounit,*) " Coefficient        ::  ", (basis_coefficient(i), i=1, basis_num_prim)
        write(ounit,*) " Prim. Factor       ::  ", (basis_prim_factor(i), i=1, basis_num_prim)
        write(ounit,*) " AO Shell           ::  ", (ao_shell(i), i=1, ao_num)
        write(ounit,*) " AO Normalization   ::  ", (ao_normalization(i), i=1, ao_num)




        write(ounit,*)


        write(ounit,*) '-----------------------------------------------------------------------'

        ! do j= 1, ncent
        !     write(ounit,'(A4, 2x, 3F12.6, 2x, i3)') symbol(j), (cent(i,j),i=1,3), iwctype(j)
        ! enddo


        write(ounit,*)
    end subroutine read_trexio_basis_file




end module