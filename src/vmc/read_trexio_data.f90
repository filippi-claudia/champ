module trexio_read_data
    use error, only : fatal_error
    use precision_kinds,        only: dp
    use m_trexio_basis,         only: gnorm
    use array_utils,            only: unique_elements

    private
    public :: dp
    public :: read_trexio_molecule_file
    public :: read_trexio_symmetry_file
    public :: read_trexio_orbitals_file
    public :: read_trexio_basis_file
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
        write(ounit,*) " Reading molecular coordinates from the trexio file :: ",  file_trexio_path
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
        write(ounit,*) " Reading LCAO orbitals from the file :: ",  trim(file_trexio_path)
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
        integer                         :: iostat, ic, i, j, k, l, iunit, tcount1, tcount2, tcount3, tcount4
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
        real(dp)                        :: rgrid(2000)
        integer, dimension(nctype_tot)  :: icusp
        integer                         :: cartesian_shells(5) = (/1, 3, 6, 10, 15/)
        real(dp)                        :: r, r2, r3, val   ! local values

        integer, dimension(:), allocatable :: atom_index(:), shell_index_atom(:), nshells_per_atom(:)
        integer, dimension(:), allocatable :: prim_index_atom(:), nprims_per_atom(:)
        integer, dimension(:), allocatable :: unique_atom_index(:)
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
            ! write(*,*) "trexio_read_basis_nucleus_index :: ", basis_nucleus_index
            rc = trexio_read_basis_shell_index(trex_basis_file, basis_shell_index)
            ! write(*,*) "trexio_read_basis_shell_index :: ", basis_shell_index
            rc = trexio_read_basis_shell_ang_mom(trex_basis_file, basis_shell_ang_mom)
            ! write(*,*) "trexio_read_basis_shell_ang_mom :: ", basis_shell_ang_mom
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

        ! print*, "prim_index_atom(1:ncent_tot) :: ", prim_index_atom(1:ncent_tot)
        ! print*, "nprims_per_atom(1:ncent_tot) :: ", nprims_per_atom(1:ncent_tot)


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





            do j = 1, basis_num_shell   ! loop over all the shells
                ! select the shells corresponding to the unique atoms only
                if (unique_atom_index(ic) == basis_nucleus_index(j)) then
                    ! j is the running shell index for the unique atom i
                    print *, "j ", j, "basis_nucleus_index ", basis_nucleus_index(j), "basis_shell ang mom  ", basis_shell_ang_mom(j)
                    ! loop over all the gridpoints to add contracted Gaussians over the grid
                    do i = 1, gridpoints
                        ! Generate the grid here for the unique atom
                        r = rgrid(i)
                        r2 = r*r
                        r3 = r2*r
                        val = 0.0d0

                    enddo
                endif
            enddo
        enddo




        ! Extract the shell angular momentum information only for unique type
        ! of atoms.



    !     def compute_grid():
    !     # Compute the radial grid r for a given number of grid points
    !     # and grid type
    !     for i in range(gridpoints):
    !         if gridtype == 1:
    !             r = gridr0 + i*gridarg
    !         elif gridtype == 2:
    !             r = gridr0 * gridarg**i
    !         elif gridtype == 3:
    !             r = gridr0 * gridarg**i - gridr0
    !         bgrid[:,i] = r
    !     return bgrid

    ! def add_function(shell_ang_mom, exponents, coefficients, shell, bgrid):
    !     # put a new function on the grid
    !     # The function is defined by the exponent, coefficient and type
    !     for i in range(gridpoints):
    !         r = bgrid[shell+1, i]
    !         r2 = r*r
    !         r3 = r2*r
    !         value = 0.0
    !         for j in range(len(exponents)):
    !             value += gnorm(exponents[j], shell_ang_mom) * coefficients[j] * np.exp(-exponents[j]*r2)
    !             # print ("each value k, ib,", i ,j , value)

    !         bgrid[shell+1,i] = value

    !     return

    ! if filename is not None:
    !     if isinstance(filename, str):
    !         unique_elements, indices = np.unique(nucleus_label, return_index=True)

    !         for i in range(len(unique_elements)):
    !             # Write down an radial basis grid file in the new champ v2.0 format for each unique atom type
    !             filename_basis_grid = "BASISGRID." + 'basis.' + unique_elements[i]
    !             with open(filename_basis_grid, 'w') as file:

    !                 # Common numbers
    !                 gridtype=3
    !                 gridpoints=2000
    !                 gridarg=1.003
    !                 gridr0=20.0

    !                 number_of_shells_per_atom = list_nshells[indices[i]]

    !                 shell_ang_mom_per_atom_list = []
    !                 for ind, val in enumerate(dict_basis["nucleus_index"]):
    !                     if val == indices[i]:
    !                         shell_ang_mom_per_atom_list.append(dict_basis["shell_ang_mom"][ind])

    !                 shell_ang_mom_per_atom_count = Counter(shell_ang_mom_per_atom_list)

    !                 total_shells = sum(shell_ang_mom_per_atom_count.values())

    !                 shells_per_atom = {}
    !                 for count in shell_ang_mom_per_atom_count:
    !                     shells_per_atom[count] = shell_ang_mom_per_atom_count[count]

    !                 bgrid = np.zeros((number_of_shells_per_atom+1, gridpoints))


    !                 ## The main part of the file starts here
    !                 gridr0_save = gridr0
    !                 if gridtype == 3:
    !                     gridr0 = gridr0/(gridarg**(gridpoints-1)-1)


    !                 bgrid = compute_grid()  # Compute the grid, store the results in bgrid

    !                 ### Note temp index should point to shell index of unique atoms
    !                 # get the exponents and coefficients of unique atom types
    !                 counter = 0
    !                 for ind, val in enumerate(dict_basis["nucleus_index"]):
    !                     if val == indices[i]:
    !                         shell_index_unique_atom = index_radial[indices[i]][counter]
    !                         list_contracted_exponents =  contr[shell_index_unique_atom]["exponent"]
    !                         list_contracted_coefficients =  contr[shell_index_unique_atom]["coefficient"]
    !                         add_function(dict_basis["shell_ang_mom"][ind], list_contracted_exponents, list_contracted_coefficients, counter, bgrid)
    !                         counter += 1


    !                 # file writing part
    !                 file.write(f"{number_of_shells_per_atom} {gridtype} {gridpoints} {gridarg:0.6f} {gridr0_save:0.6f} {0}\n")
    !                 np.savetxt(file, np.transpose(bgrid), fmt=' %.12e')

    !             file.close()
    !     else:
    !         raise ValueError
    ! # If filename is None, return a string representation of the output.
    ! else:
    !     return None



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


#if defined(TREXIO_FOUND)
        use trexio
        use contrl_file,        only: backend
#endif

        implicit none

        !   local use
        character(len=72), intent(in)   :: file_trexio
        character(len=128)              :: file_trexio_path
        integer                         :: iostat, i, j, k, iunit
        logical                         :: exist, skip = .true.
        character(len=40)               :: label
        integer                         :: io, nsym, mo_num
        character(len=10), allocatable  :: mo_symmetry(:)


        ! trexio
        integer(8)                      :: trex_symmetry_file
        integer                         :: rc = 1

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
            write(*,*) "trexio_open from symmetry :: ", rc
            rc = trexio_read_mo_num(trex_symmetry_file, mo_num)
            write(*,*) "trexio_read_mo_num from symmetry:: ", rc
#endif
        endif
        call bcast(mo_num)
        write(*,*) "trexio symmetry file :: ", trex_symmetry_file
        ! safe allocate
        if (.not. allocated(irrep)) allocate (irrep(mo_num))
        if (.not. allocated(mo_symmetry)) allocate (mo_symmetry(mo_num))

!         if (wid) then
! #if defined(TREXIO_FOUND)
!             rc = trexio_read_mo_symmetry(trex_symmetry_file, mo_symmetry)
!             ! write(*,*) "trexio_read_mo_symmetry :: ", rc
! #endif
!         endif
!         call bcast(mo_symmetry)

        write(ounit,fmt=int_format) " Number of MOs ::  ", mo_num
        ! write(ounit,*) " MO Symmetry ::  ", (mo_symmetry(i), i=1, mo_num)


        ! ! read data
        ! if (wid) then
        !     read (iunit, *, iostat=iostat) (irrep(io), io=1, norb)
        !     if (iostat/=0) call fatal_error("Error in reading symmetry file :: expecting irrep correspondence for all norb orbitals")
        ! endif
        ! call bcast(irrep)

        ! write(ounit, '(10(1x, i3))') (irrep(io), io=1, norb)

        ! if (wid) close(iunit)
    end subroutine read_trexio_symmetry_file




end module