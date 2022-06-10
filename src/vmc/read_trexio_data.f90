module trexio_read_data
    use error, only : fatal_error
    use precision_kinds,        only: dp
    use array_utils,            only: unique_elements

    private
    public :: dp
    public :: read_trexio_molecule_file
    public :: read_trexio_symmetry_file
    public :: read_trexio_orbitals_file
    public :: read_trexio_basis_file
    public :: read_trexio_determinant_file
    public :: read_trexio_ecp_file
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
        character(len=*), intent(in)   :: file_trexio
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
        write(ounit,*) " Reading molecular coordinates from the trexio file :: ",  file_trexio_path
        write(ounit,*) '-----------------------------------------------------------------------'

        ! Check if the file exists
        if (wid) then
#if defined(TREXIO_FOUND)
            trex_molecule_file = trexio_open(file_trexio_path, 'r', backend, rc)
            call trexio_assert(rc, TREXIO_SUCCESS)
            rc = trexio_read_nucleus_num(trex_molecule_file, ncent)
            call trexio_assert(rc, TREXIO_SUCCESS)
            rc = trexio_read_electron_up_num(trex_molecule_file, nup)
            call trexio_assert(rc, TREXIO_SUCCESS)
            rc = trexio_read_electron_dn_num(trex_molecule_file, ndn)
            call trexio_assert(rc, TREXIO_SUCCESS)
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
        call trexio_assert(rc, TREXIO_SUCCESS)
        rc = trexio_read_nucleus_label(trex_molecule_file, symbol, 3)
        call trexio_assert(rc, TREXIO_SUCCESS)
        rc = trexio_close(trex_molecule_file)
        call trexio_assert(rc, TREXIO_SUCCESS)
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
        write(ounit,*) " Reading LCAO orbitals from the file :: ",  trim(file_trexio_path)
        write(ounit,*) '---------------------------------------------------------------------------'
        ! Check if the file exists

        if (wid) then
#if defined(TREXIO_FOUND)
            trex_orbitals_file = trexio_open(file_trexio_path, 'r', backend, rc)
            call trexio_assert(rc, TREXIO_SUCCESS)
            rc = trexio_read_mo_num(trex_orbitals_file, norb)
            call trexio_assert(rc, TREXIO_SUCCESS)
            rc = trexio_read_ao_num(trex_orbitals_file, nbasis)
            call trexio_assert(rc, TREXIO_SUCCESS)
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
            call trexio_assert(rc, TREXIO_SUCCESS)
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

        ! The following to be used to store the information
        use numbas_mod,         only: MRWF, MRWF_PTS
        use atom,               only: znuc, nctype, nctype_tot, ncent_tot
        use atom,               only: symbol, atomtyp
        use vmc_mod,            only: NCOEF
        use ghostatom,          only: newghostype
        use const,              only: ipr
        use numbas,             only: arg, d2rwf, igrid, nr, nrbas, r0, rwf
        use numbas,             only: allocate_numbas
        use coefs,              only: nbasis
        use numexp,             only: ae, ce, ab, allocate_numexp
        use pseudo,             only: nloc
        use general,            only: filename, filenames_bas_num

        ! For processing the stored information
        use atom, 			    only: atomtyp
        use general, 			only: pooldir, bas_id
        use contrl_file,        only: ounit, errunit
        use spline2_mod,        only: spline2
        use fitting_methods,    only: exp_fit

#if defined(TREXIO_FOUND)
        use trexio
        use contrl_file,        only: backend
        use error,              only: trexio_error
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
        integer                         :: iostat, ic, ir, i, j, k, l, iunit, tcount1, tcount2, tcount3, tcount4, tcount5
        logical                         :: exist
        type(atom_t)                    :: atoms

        ! trexio
        integer(8)                      :: trex_basis_file
        integer                         :: rc = 1

        !   Formatting
        character(len=128)              :: int_format     = '(A, T60, I0)'
        character(len=128)              :: float_format   = '(A, T60, f12.8)'
        character(len=128)              :: string_format  = '(A, T60, A128)'

        ! Grid related
        integer                         :: gridtype=3
        integer                         :: gridpoints=2000
        real(dp)                        :: gridarg=1.003
        real(dp)                        :: gridr0=20.0
        real(dp)                        :: gridr0_save = 20.0
        real(dp)                        :: rgrid(2000)  ! Grid points
        integer, dimension(nctype_tot)  :: icusp
        integer                         :: cartesian_shells(5) = (/1, 3, 6, 10, 15/)
        real(dp)                        :: r, r2, r3, val   ! local values

        integer, dimension(:), allocatable :: atom_index(:), shell_index_atom(:), nshells_per_atom(:)
        integer, dimension(:), allocatable :: prim_index_atom(:), nprims_per_atom(:)
        integer, dimension(:), allocatable :: unique_atom_index(:), shell_prim_correspondence(:)
        integer                         :: count
        character(len=2), allocatable   :: unique(:) ! unique symbols of atoms

        trex_basis_file = 0

        !   External file reading

        if((file_trexio(1:6) == '$pool/') .or. (file_trexio(1:6) == '$POOL/')) then
            file_trexio_path = pooldir // file_trexio(7:)
        else
            file_trexio_path = file_trexio
        endif

        write(ounit,*) '-----------------------------------------------------------------------'
        write(ounit,*) " Reading Basis Set information from the trexio file :: ", trim(adjustl(file_trexio_path))
        write(ounit,*) '-----------------------------------------------------------------------'

        ! Check if the file exists
        if (wid) then
#if defined(TREXIO_FOUND)
            trex_basis_file = trexio_open(file_trexio_path, 'r', backend, rc)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio file open error', __FILE__, __LINE__)
            rc = trexio_read_basis_prim_num(trex_basis_file, basis_num_prim)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_basis_prim_num', __FILE__, __LINE__)
            rc = trexio_read_basis_shell_num(trex_basis_file, basis_num_shell)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_basis_shell_num', __FILE__, __LINE__)
            rc = trexio_read_ao_num(trex_basis_file, ao_num)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_ao_num', __FILE__, __LINE__)
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
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio file open error', __FILE__, __LINE__)
            rc = trexio_read_basis_nucleus_index(trex_basis_file, basis_nucleus_index)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_basis_nucleus_index', __FILE__, __LINE__)
            rc = trexio_read_basis_shell_index(trex_basis_file, basis_shell_index)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_basis_shell_index', __FILE__, __LINE__)
            rc = trexio_read_basis_shell_ang_mom(trex_basis_file, basis_shell_ang_mom)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_basis_shell_ang_mom', __FILE__, __LINE__)
            print *, "shell ang mom ", basis_shell_ang_mom
            rc = trexio_read_basis_shell_factor(trex_basis_file, basis_shell_factor)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_basis_shell_factor', __FILE__, __LINE__)
            rc = trexio_read_basis_exponent(trex_basis_file, basis_exponent)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_basis_exponent', __FILE__, __LINE__)
            rc = trexio_read_basis_coefficient(trex_basis_file, basis_coefficient)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_basis_coefficient', __FILE__, __LINE__)
            rc = trexio_read_basis_prim_factor(trex_basis_file, basis_prim_factor)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_basis_prim_factor', __FILE__, __LINE__)
            rc = trexio_read_ao_shell(trex_basis_file, ao_shell)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_ao_shell', __FILE__, __LINE__)
            rc = trexio_read_ao_normalization(trex_basis_file, ao_normalization)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_ao_normalization', __FILE__, __LINE__)
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



        write(ounit,fmt=int_format) " Number of primitives  ::  ", basis_num_prim
        write(ounit,fmt=int_format) " Number of shells      ::  ", basis_num_shell
        write(ounit,fmt=int_format) " Number of AO          ::  ", ao_num
        write(ounit,*)
        write(ounit,*) '-----------------------------------------------------------------------'
        write(ounit,*)

        ! Processing the basis set information to get the numerican grid
        gridr0_save = gridr0


        ! Get the number of shells per atom (information needed to reshuffle AOs)

        allocate(atom_index(basis_num_shell))
        allocate(nshells_per_atom(basis_num_shell))
        allocate(shell_index_atom(basis_num_shell))




        call unique_elements(basis_num_shell, basis_nucleus_index, atom_index, count, nshells_per_atom, shell_index_atom)

        print*, "Number of unique elements :: ", count
        print*, "Unique elements index :: ", atom_index(1:count)
        print*, "frequency :: ", nshells_per_atom(1:count)
        print*, "result", shell_index_atom(1:count)

        ! count                     :: "Number of unique elements"
        ! atom_index(1:count)       :: "Unique elements index (not used here)"
        ! nshells_per_atom(1:count) :: "frequency or count of shells per atom"
        ! shell_index_atom(1:count) :: "index number of the shell for each atom"
        ! The shells per atom can be obtained by accessing the shell_index_atom
        ! for a given atom index by the slice of size frequency.


        ! Now get the number of primitives per atom and their indices
        allocate(nprims_per_atom(basis_num_prim))
        allocate(prim_index_atom(basis_num_prim))

        prim_index_atom(1) = 1
        nprims_per_atom(1) = 0
        tcount3 = 0; tcount4 = 0
        do i = 1, ncent_tot
            do j = 1, nshells_per_atom(i)   ! frequency
                tcount4 = tcount4 + 1
                nprims_per_atom(i) = nprims_per_atom(i) + cartesian_shells(basis_shell_ang_mom(tcount4)+1)
            enddo
            if (i .ne. ncent_tot) prim_index_atom(i+1) = prim_index_atom(i) + nprims_per_atom(i)
        enddo

        print*, "prim_index_atom(1:ncent_tot) :: ", prim_index_atom(1:ncent_tot)
        print*, "nprims_per_atom(1:ncent_tot) :: ", nprims_per_atom(1:ncent_tot)
        print*, "tcount4 :: ", tcount4


        ! Obtain the number of unique types of atoms stored in the hdf5 file.
        print*, "number of types of atoms :: ", nctype_tot
        print*, "nucleus shell index ", basis_nucleus_index
        if (.not. allocated(unique)) allocate(unique(nctype_tot))
        if (.not. allocated(unique_atom_index)) allocate(unique_atom_index(nctype_tot))


        print*, "symbol ", symbol
        print*, "atom_type ", atomtyp


        tcount1 = 1; tcount2 = 1
        unique_atom_index(1) = 1
        unique(1) = symbol(1)
        do j= 2, ncent_tot
            if (any(unique == symbol(j) ))  then
                cycle
            endif
            print*, "j ", j, "symbol ", symbol(j)
            tcount1 = tcount1 + 1
            unique_atom_index(tcount1) = j
            unique(tcount1) = symbol(j)
        enddo
        print *, "tcount1 ", tcount1, "unique ", unique(1:tcount1)
        print*, "unique ", unique
        print*, "unique atom index ", unique_atom_index


        ! start putting in the information in the arrays and variables
        gridtype=3
        gridpoints=2000
        gridarg=1.003
        gridr0=20.0
        gridr0_save = gridr0

        ! Do the necessary allocation for the numerical basis set
        call allocate_numbas()
        call allocate_numexp()

        ! count                     :: "Number of unique elements"
        ! atom_index(1:count)       :: "Unique elements index (not used here)"
        ! nshells_per_atom(1:count) :: "frequency or count of shells per atom"
        ! shell_index_atom(1:count) :: "index number of the shell for each atom"
        ! The shells per atom can be obtained by accessing the shell_index_atom
        ! for a given atom index by the slice of size frequency.

        ! Populate the rgrid array.
        do i = 1, gridpoints
            if (gridtype .eq. 1) then
                rgrid(i) = gridr0 + i*gridarg
            else if (gridtype .eq. 2) then
                rgrid(i) = gridr0 * gridarg**i
            else if (gridtype .eq. 3) then
                rgrid(i) = gridr0*gridarg**i - gridr0
            endif
        enddo

        tcount5 = 0
        allocate(shell_prim_correspondence(basis_num_prim))
        do j = 1, basis_num_shell
            if (basis_shell_ang_mom(j) .eq. 0) then
                tcount5 = tcount5 + 1
                shell_prim_correspondence(tcount5) = 0
            elseif (basis_shell_ang_mom(j) .eq. 1) then
                tcount5 = tcount5 + 1
                shell_prim_correspondence(tcount5:tcount5+3) = 1
            elseif (basis_shell_ang_mom(j) .eq. 2) then
                shell_prim_correspondence(j:j+6) = 2
            elseif (basis_shell_ang_mom(j) .eq. 3) then
                shell_prim_correspondence(j:j+10) = 3
            elseif (basis_shell_ang_mom(j) .eq. 4) then
                shell_prim_correspondence(j:j+15) = 4
            elseif (basis_shell_ang_mom(j) .eq. 5) then
                shell_prim_correspondence(j:j+21) = 5
            endif
        enddo

        print *, "shell prim correspondence ", shell_prim_correspondence




        do ic = 1, nctype_tot            ! loop over all the unique atoms
            nrbas(ic)   = nshells_per_atom(shell_index_atom(ic))
            igrid(ic)   = gridtype       ! grid type default is 3
            nr(ic)      = gridpoints     ! number of grid points default is 2000
            arg(ic)     = gridarg        ! grid spacing default is 1.003
            r0(ic)      = gridr0         ! grid origin default is 20.0
            icusp(ic)   = 0              ! default is 0

            write(ounit,*)
            write(ounit,'(A, T60,  A)')     " For Nucleus           ::  ", unique(ic)
            write(ounit,'(A, T60, I0)')     " Number of Shells      ::  ", nrbas(ic)
            write(ounit,'(A, T60, I0)')     " Grid type             ::  ", igrid(ic)
            write(ounit,'(A, T60, I0)')     " Number of Grid Points ::  ", nr(ic)
            write(ounit,'(A, T56, F10.6)')  " Grid spacing          ::  ", arg(ic)
            write(ounit,'(A, T56, F10.6)')  " Grid origin           ::  ", r0(ic)
            write(ounit,'(A, T60, I0)')     " Icusp                 ::  ", icusp(ic)
            write(ounit,*)

            ! Make space for the special case when nloc == 0


            if (gridtype .eq. 3) gridr0 = gridr0/(gridarg**(gridpoints-1)-1)

            print*, "all primitive exponents ", basis_exponent(1:35)
            ! loop over all the primitives for the unique atom
            val = 0.0d0
            do k = prim_index_atom(unique_atom_index(ic)), prim_index_atom(unique_atom_index(ic)) + nprims_per_atom(unique_atom_index(ic)) - 1
                ! k is index of primitives that needs to used for adding to the grid.
                ! gnorm(exponents[j], shell_ang_mom) * coefficients[j] * np.exp(-exponents[j]*r2)
                print*, "the primi list k ", k
            enddo



            do j = 1, basis_num_shell   ! loop over all the shells
                ! select the shells corresponding to the unique atoms only
                if (unique_atom_index(ic) == basis_nucleus_index(j)) then
                    ! j is the running shell index for the unique atom i
                    print *, "j ", j, "basis_nucleus_index ", basis_nucleus_index(j), "basis_shell ang mom  ", basis_shell_ang_mom(j)
                    ! loop over all the gridpoints to add contracted Gaussians over the grid
                    ! list of exponents


                    do ir = 1, gridpoints
                        ! Generate the grid here for the unique atom
                        r = rgrid(ir)
                        r2 = r*r
                        r3 = r2*r
                        val = 0.0d0
                        ! rwf(ir,j,ic,iwf) = compute_grid_value(r, basis_shell_ang_mom(j), basis_shell_coeff(j), basis_shell_exp(j), val)
                    enddo
                endif
            enddo
        enddo


        ! Put in the read information in the x(ir) and rwf(ir,j,ic,iwf) arrays
        !     do ir=1,nr(ic)
        !       read(iunit,*,iostat=iostat) x(ir),(rwf(ir,irb,ic,iwf),irb=1,nrbas(ic))
        !     enddo




    contains

    ! Functions needed only for the radial basis set generation
    double precision function shell_to_grid(l, exponents, coefficients) result (val)
    use m_trexio_basis,             only: gnorm
    implicit None
    integer, intent(in)             :: l
    double precision, intent(in)    :: exponents(:)
    double precision, intent(in)    :: coefficients(:)
    integer                         :: i,j
    double precision                :: r, r2

    do i = 1, gridpoints
        r = rgrid(i)
        r2 = r*r
        val = 0.0d0
        do j = 1, size(exponents)
            val = val + gnorm(exponents(j), l) * coefficients(j) * dexp(-exponents(j)*r2)
        enddo
    enddo
    end function shell_to_grid

    end subroutine read_trexio_basis_file



    subroutine read_trexio_symmetry_file(file_trexio)
        ! Ravindra

        use custom_broadcast,   only: bcast
        use mpiconf,            only: wid, idtask

        use contrl_file,        only: ounit, errunit
        use coefs,              only: norb
        use optorb,             only: irrep
        use vmc_mod,            only: norb_tot
        use general,            only: pooldir
        use precision_kinds,    only: dp
        use array_utils,        only: unique_string_elements
         


#if defined(TREXIO_FOUND)
        use trexio
        use contrl_file,        only: backend
        use error,              only: trexio_error
#endif

        implicit none

        !   local use
        character(len=72), intent(in)   :: file_trexio
        character(len=128)              :: file_trexio_path
        integer                         :: iostat, i, j, k, iunit
        logical                         :: exist, skip = .true.
        character(len=40)               :: label
        integer                         :: io, nsym, mo_num
        character(len=3), allocatable   :: mo_symmetry(:)


        ! trexio
        integer(8)                      :: trex_symmetry_file
        integer                         :: rc = 1

        character(len=3), dimension(:), allocatable :: unique_irrep       ! The output
        integer                                     :: num_irrep          ! The number of unique elements




        !   Formatting
        character(len=100)               :: int_format     = '(A, T60, I0)'
        character(len=100)               :: string_format  = '(A, T60, A)'

        !   External file reading

        if((file_trexio(1:6) == '$pool/') .or. (file_trexio(1:6) == '$POOL/')) then
            file_trexio_path = pooldir // file_trexio(7:)
        else
            file_trexio_path = file_trexio
        endif

        write(ounit,*) '-----------------------------------------------------------------------'
        write(ounit,*) " Reading orbital symmetries information from the trexio file :: ", trim(adjustl(file_trexio_path))
        write(ounit,*) '-----------------------------------------------------------------------'

        ! Check if the file exists
        if (wid) then
#if defined(TREXIO_FOUND)
            trex_symmetry_file = trexio_open(file_trexio_path, 'r', backend, rc)
            rc = trexio_read_mo_num(trex_symmetry_file, mo_num)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_mo_num failed', __FILE__, __LINE__)
#endif
        endif
        call bcast(mo_num)
        ! safe allocate
        if (.not. allocated(irrep)) allocate (irrep(mo_num))
        if (.not. allocated(mo_symmetry)) allocate (mo_symmetry(mo_num))
        if (.not. allocated(unique_irrep)) allocate (unique_irrep(mo_num))

        if (wid) then
#if defined(TREXIO_FOUND)
            rc = trexio_read_mo_symmetry(trex_symmetry_file, mo_symmetry, 2)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_mo_symmetry failed', __FILE__, __LINE__)
#endif
        endif
        call bcast(mo_symmetry)

        write(ounit,fmt=int_format) " Number of molecular orbital symmetries read ::  ", mo_num

        call unique_string_elements(mo_num, mo_symmetry, unique_irrep, num_irrep)


        write(ounit,fmt=int_format) " Number of irreducible representations       ::  ", num_irrep
        write(ounit,*) 
        write(ounit,'(a)')          " Irreducible representations correspondence  ::  "
        write(ounit,'(1x,10(a2,a,i2,a,x))') (unique_irrep(i), "=", i ,";", i=1, num_irrep)
        write(ounit,*) 

        ! get the correspondence for each atom according to the rule defined for atomtypes
        do j = 1, mo_num
            do k = 1, num_irrep
                if (mo_symmetry(j) == unique_irrep(k))   irrep(j) = k
            enddo
        enddo


        write(ounit,*)  "Irreducible representation correspondence for all molecular orbitals"
        write(ounit, '(10(1x, i3))') (irrep(i), i=1, mo_num)


    end subroutine read_trexio_symmetry_file


    subroutine read_trexio_determinant_file(file_trexio)
        !> This subroutine reads the .hdf5 trexio generated file/folder. It then reads the
        !> determinant coefficients and orbital occupations .
        !! @author Ravindra Shinde (r.l.shinde@utwente.nl)
        !! @date 25 May 2022
        use custom_broadcast,   only: bcast
        use mpiconf,            only: wid
        use, intrinsic :: iso_fortran_env, only: iostat_eor
        use contrl_file,    only: ounit, errunit
        use general,            only: pooldir
        use dets,           only: cdet, ndet
        use dorb_m,         only: iworbd
        use coefs,          only: norb
        use inputflags,     only: ideterminants
        use wfsec,          only: nwftype
        use csfs,           only: nstates
        use mstates_mod,    only: MSTATES
        use general,        only: pooldir
        use elec,           only: ndn, nup
        use const,          only: nelec
        use method_opt,     only: method
        use precision_kinds, only: dp

#if defined(TREXIO_FOUND)
        use trexio
        use contrl_file,        only: backend
        use error,              only: trexio_error
#endif

        implicit none

        !   local use
        character(len=72), intent(in)   :: file_trexio
        character(len=40)               :: temp1, temp2, temp3, temp4
        character(len=80)               :: comment, file_trexio_path
        integer                         :: iostat, i, j, k, iunit
        logical                         :: exist
        character(len=2), allocatable   :: unique(:)

        ! trexio
        integer(8)                      :: trex_determinant_file
        integer                         :: rc = 1

        !   Formatting
        character(len=100)              :: int_format     = '(A, T60, I0)'
        character(len=100)              :: float_format   = '(A, T60, f12.8)'
        character(len=100)              :: string_format  = '(A, T60, A)'

        ! determinant data (debugging)
        integer*8, allocatable :: det_list(:)
        integer*8 :: read_buf_det_size      ! how many do you want
        integer*8 :: jj, offset_det_read = 0    ! How many first you want to skip
        integer*8 :: chunk_det_read = 1
        integer*8 :: determinant_num
        integer   :: int64_num           ! Number of intergers required per spin component
        ! orbital lists (debugging)
        integer*4, allocatable :: orb_list_up(:), orb_list_dn(:)
        integer*8, allocatable :: orb_list(:)
        integer*4 :: occ_num_up, occ_num_dn, occupied_num


        trex_determinant_file = 0

        !   External file reading

        if((file_trexio(1:6) == '$pool/') .or. (file_trexio(1:6) == '$POOL/')) then
            file_trexio_path = pooldir // file_trexio(7:)
        else
            file_trexio_path = file_trexio
        endif

        file_trexio_path = 'hcn.trexio.hdf5'

        write(ounit,*) '-----------------------------------------------------------------------'
        write(ounit,*) " Reading determinants from the trexio file :: ",  file_trexio_path
        write(ounit,*) '-----------------------------------------------------------------------'

        ! Check if the file exists
        if (wid) then
#if defined(TREXIO_FOUND)
            trex_determinant_file = trexio_open(file_trexio_path, 'r', backend, rc)
            call trexio_assert(rc, TREXIO_SUCCESS)
            rc = trexio_has_determinant_num (trex_determinant_file)
            if (rc == TREXIO_SUCCESS) then
                rc = trexio_read_determinant_num(trex_determinant_file, ndet)
                call trexio_assert(rc, TREXIO_SUCCESS)
                rc = trexio_get_int64_num(trex_determinant_file, int64_num)
                call trexio_assert(rc, TREXIO_SUCCESS)
            else
                write(errunit,*) "trexio file does not have number of determinant  stored :: ", rc
                call trexio_error(rc, TREXIO_SUCCESS, 'trexio_has_determinant_num failed', __FILE__, __LINE__)
            endif
#endif
        endif
        call bcast(ndet)
        call bcast(int64_num)

        determinant_num = ndet
        write(ounit,int_format) " Number of determinants (read from trexio) :: ", ndet

!       Do the allocations based on the number of determinants and the method
        if( (method(1:3) == 'lin')) then
            if (.not. allocated(cdet)) allocate(cdet(ndet,MSTATES,3))
        else
            if (.not. allocated(cdet)) allocate(cdet(ndet,MSTATES,nwftype))
        endif

        read_buf_det_size = determinant_num

        allocate(det_list(ndet))

        if (wid) then
#if defined(TREXIO_FOUND)
        rc = trexio_read_determinant_coefficient(trex_determinant_file, offset_det_read, read_buf_det_size, cdet(:,1,nwftype))
        call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_determinant_coeff failed', __FILE__, __LINE__)
#endif
        endif
        call bcast(cdet)


        write(ounit,*)
        write(ounit,*) " Determinant coefficients "
        write(ounit,'(10(1x, f11.8, 1x))') (cdet(i,1,nwftype), i=1, read_buf_det_size)

!       allocate the orbital mapping array
        if (.not. allocated(iworbd)) allocate(iworbd(nelec, ndet))
        allocate(orb_list(ndet))
        allocate(orb_list_up(int64_num))
        allocate(orb_list_dn(int64_num))

        write(ounit, *)
        write(ounit, *) "Orbitals <--> Determinants mapping read from a trexio file :: "
        write(ounit, *) "Serial numbers of orbitals that are occupied               :: "
        write(ounit, *) "'alpha (spin up)'  <---------------------->  'beta (spin down)' "
        write(ounit, *)
        ! convert one given determinant into lists of orbitals

        read_buf_det_size = int64_num
        offset_det_read = 0
        do jj = 1, determinant_num
            rc = trexio_read_determinant_list(trex_determinant_file, offset_det_read, read_buf_det_size, orb_list(jj))
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_determinant_list failed', __FILE__, __LINE__)
            rc = trexio_to_orbital_list_up_dn(int64_num, orb_list(jj), orb_list_up, orb_list_dn, occ_num_up, occ_num_dn)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_to_orbital_list_up_dn filed', __FILE__, __LINE__)
            write(ounit,'(<occ_num_up>(i4,1x), 2x, <occ_num_dn>(i4,1x))') (orb_list_up(i), i = 1, occ_num_up), (orb_list_dn(i), i = 1, occ_num_dn)

            do i = 1, occ_num_up
                iworbd(i, jj) = orb_list_up(i)
            enddo
            do i = 1, occ_num_dn
                iworbd(occ_num_up + i, jj) = orb_list_dn(i)
            enddo

            offset_det_read = jj
        enddo


        write(ounit,*) '-----------------------------------------------------------------------'
        write(ounit,*)
    end subroutine read_trexio_determinant_file


    subroutine read_trexio_ecp_file(file_trexio)
        !> This subroutine reads the .hdf5 trexio generated file/folder. It then reads the
        !> ECP information for all the unique atoms.
        !! @author Ravindra Shinde (r.l.shinde@utwente.nl)
        !! @date 01 June 2022

        use custom_broadcast,   only: bcast
        use mpiconf,            only: wid

#if defined(TREXIO_FOUND)
        use trexio
        use contrl_file,        only: backend
        use error,              only: trexio_error
#endif

        use pseudo_mod,         only: MPS_L, MGAUSS, MPS_QUAD
        use atom,               only: nctype, atomtyp
        use gauss_ecp,          only: ecp_coef, ecp_exponent, necp_power, necp_term
        use gauss_ecp,          only: allocate_gauss_ecp
        use pseudo,             only: lpot
        use qua,                only: nquad, wq, xq0, yq0, zq0
        use general,            only: pooldir, filename, pp_id, filenames_ps_gauss
        use contrl_file,        only: ounit, errunit
        use rotqua_mod,         only: gesqua

        use precision_kinds,    only: dp

        implicit none

        !   local use
        character(len=72), intent(in)   :: file_trexio
        character(len=40)               :: temp1, temp2, temp3, temp4
        character(len=80)               :: comment, file_trexio_path
        logical                         :: exist, skip = .true.

        ! trexio
        integer(8)                      :: trex_ecp_file
        integer                         :: rc = 1

        ! local variables
        integer                         :: ecp_num
        integer, allocatable            :: flat_ecp_ang_mom(:)
        integer, allocatable            :: flat_ecp_nucleus_index(:)
        integer, allocatable            :: flat_ecp_max_ang_mom_plus_1(:)
        integer, allocatable            :: flat_ecp_power(:)
        integer, allocatable            :: flat_ecp_z_core(:)
        real(dp), allocatable           :: flat_ecp_coefficient(:)
        real(dp), allocatable           :: flat_ecp_exponent(:)


        !   Formatting
        character(len=100)              :: int_format     = '(A, T60, I0)'
        character(len=100)              :: float_format   = '(A, T60, f12.8)'
        character(len=100)              :: string_format  = '(A, T60, A)'

        integer         :: i, ic, idx, l
        integer         :: iunit, iostat, counter = 0

        character*80 label

        trex_ecp_file = 0

        !   External file reading

        if((file_trexio(1:6) == '$pool/') .or. (file_trexio(1:6) == '$POOL/')) then
            file_trexio_path = pooldir // file_trexio(7:)
        else
            file_trexio_path = file_trexio
        endif


        write(ounit,*) '-----------------------------------------------------------------------'
        write(ounit,*) " Reading ECP data from the trexio file :: ",  file_trexio_path
        write(ounit,*) '-----------------------------------------------------------------------'



        ! Check if the file exists
        if (wid) then
#if defined(TREXIO_FOUND)
            trex_ecp_file = trexio_open(file_trexio_path, 'r', backend, rc)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio file open', __FILE__, __LINE__)
            rc = trexio_read_ecp_num(trex_ecp_file, ecp_num)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_ecp_num', __FILE__, __LINE__)
#endif
        endif
        call bcast(ecp_num)
        allocate (flat_ecp_ang_mom(ecp_num))
        allocate (flat_ecp_nucleus_index(ecp_num))
        allocate (flat_ecp_max_ang_mom_plus_1(ecp_num))
        allocate (flat_ecp_power(ecp_num))
        allocate (flat_ecp_z_core(ecp_num))
        allocate (flat_ecp_coefficient(ecp_num))
        allocate (flat_ecp_exponent(ecp_num))

        if (wid) then
#if defined(TREXIO_FOUND)
            rc = trexio_read_ecp_ang_mom(trex_ecp_file, flat_ecp_ang_mom)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_ecp_ang_mom', __FILE__, __LINE__)
            rc = trexio_read_ecp_nucleus_index(trex_ecp_file, flat_ecp_nucleus_index)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_ecp_nucleus_index', __FILE__, __LINE__)
            rc = trexio_read_ecp_max_ang_mom_plus_1(trex_ecp_file, flat_ecp_max_ang_mom_plus_1)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_ecp_max_ang_mom_plus_1', __FILE__, __LINE__)
            rc = trexio_read_ecp_power(trex_ecp_file, flat_ecp_power)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_ecp_power', __FILE__, __LINE__)
            rc = trexio_read_ecp_z_core(trex_ecp_file, flat_ecp_z_core)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_ecp_z_core', __FILE__, __LINE__)
            rc = trexio_read_ecp_coefficient(trex_ecp_file, flat_ecp_coefficient)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_ecp_coefficient', __FILE__, __LINE__)
            rc = trexio_read_ecp_exponent(trex_ecp_file, flat_ecp_exponent)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_ecp_exponent', __FILE__, __LINE__)
#endif
        endif
        call bcast(flat_ecp_ang_mom)
        call bcast(flat_ecp_nucleus_index)
        call bcast(flat_ecp_max_ang_mom_plus_1)
        call bcast(flat_ecp_power)
        call bcast(flat_ecp_z_core)
        call bcast(flat_ecp_coefficient)
        call bcast(flat_ecp_exponent)


!         do ic=1,nctype
!           if (wid) then
!             if (nctype.gt.100) call fatal_error('READPS_GAUSS: nctype>100')
!             filename =  trim(pooldir) // trim(pp_id) // ".gauss_ecp.dat." // atomtyp(ic)

!             inquire(file=filename, exist=exist)
!             if (exist) then
!               open (newunit=iunit,file=filename, iostat=iostat, action='read', status='old')
!               if (iostat .ne. 0) error stop "Problem in opening the pseudopotential file (Gaussian)"
!             else
!               call fatal_error( " Pseudopotential file (Gaussian) "// filename // " does not exist.")
!             endif

!           !   External file reading
!             write(ounit,*) '-----------------------------------------------------------------------'
!             write(ounit,'(4a)')  " Reading ECP pseudopotential for ", trim(atomtyp(ic))," from the file :: ", trim(filename)
!             write(ounit,*) '-----------------------------------------------------------------------'

!         ! label

!             read(iunit,'(a80)',iostat=iostat) label
!             if (iostat .ne. 0) then
!               write(errunit,'(a)') "Error:: Problem in reading the pseudopotential file: label"
!               write(errunit,'(2a)') "Stats for nerds :: in file ",__FILE__
!               write(errunit,'(a,i6)') "at line ", __LINE__
!             endif
!             write(ounit,'(a,i4,a,a80)') 'ECP for atom type ', ic, ' label = ', adjustl(label)
!           endif
!           call bcast(label)
!           ! max projector
!           if (.not. allocated(lpot)) allocate (lpot(nctype))

!           if (wid) then
!             read(iunit,*,iostat=iostat) lpot(ic)
!             if (iostat .ne. 0) then
!               write(errunit,'(a)') "Error:: Problem in reading the pseudopotential file: lpot"
!               write(errunit,'(2a)') "Stats for nerds :: in file ",__FILE__
!               write(errunit,'(a,i6)') "at line ", __LINE__
!             endif
!             write(ounit,'(a,i4,a,i4)') 'ECP for atom type ', ic, ' lpot = ', lpot(ic)

!             if(lpot(ic).gt.MPS_L) call fatal_error('READPS_GAUSS: increase MPS_L')
!           endif
!           call bcast(lpot)
!         ! read terms of local part and all non-local parts
!         ! local part first in file, but stored at index lpot
!         ! non-local l=0 at index 1 etc, up to lpot-1

!           call allocate_gauss_ecp()
!           do l=1,lpot(ic)
!               if(l.eq.1)then
!                 idx=lpot(ic)
!                 else
!                 idx=l-1
!               endif
!               if (wid) then
!                 read(iunit,*,iostat=iostat) necp_term(idx,ic)
!                 if (iostat .ne. 0) then
!                     write(errunit,'(a)') "Error:: Problem in reading the pseudopotential file: necp_term"
!                     write(errunit,'(2a)') "Stats for nerds :: in file ",__FILE__
!                     write(errunit,'(a,i6)') "at line ", __LINE__
!                 endif
!               endif
!               call bcast(necp_term)

!               if(necp_term(idx,ic).gt.MGAUSS) call fatal_error('READPS_GAUSS: increase MGAUSS')

!               write(ounit,'(a,2i6)') '    component, #terms ', l,necp_term(idx,ic)

!               do i=1,necp_term(idx,ic)
!                 if (wid) then
!                   read(iunit,*,iostat=iostat) ecp_coef(i,idx,ic), necp_power(i,idx,ic),ecp_exponent(i,idx,ic)

!                   if (iostat .ne. 0) then
!                     write(errunit,'(a)') "Error:: Problem in reading the pseudopotential file: ecp_coeff, power, ecp_exponents"
!                     write(errunit,'(2a)') "Stats for nerds :: in file ",__FILE__
!                     write(errunit,'(a,i6)') "at line ", __LINE__
!                   endif
!                   write(ounit,'(a,f16.8,i2,f16.8)') '    coef, power, expo ', ecp_coef(i,idx,ic), &
!                                                           necp_power(i,idx,ic), ecp_exponent(i,idx,ic)
!                 endif
!               enddo
!               call bcast(ecp_coef)
!               call bcast(necp_power)
!               call bcast(ecp_exponent)
!           enddo

!           if (wid) close(iunit)
!         enddo

!         if (.not. allocated(wq)) allocate (wq(MPS_QUAD))
!         if (.not. allocated(xq0)) allocate (xq0(MPS_QUAD))
!         if (.not. allocated(yq0)) allocate (yq0(MPS_QUAD))
!         if (.not. allocated(zq0)) allocate (zq0(MPS_QUAD))

!         call gesqua(nquad,xq0,yq0,zq0,wq)
!         call bcast(wq)
!         call bcast(xq0)
!         call bcast(yq0)
!         call bcast(zq0)
        return
      end subroutine read_trexio_ecp_file


end module
