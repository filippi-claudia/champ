module trexio_read_data
    use error, only: fatal_error
    use mpi
    use precision_kinds, only: dp
    use array_utils, only: unique_elements

    logical :: trexio_has_group_molecule      = .false.
    logical :: trexio_has_group_symmetry      = .false.
    logical :: trexio_has_group_orbitals      = .false.
    logical :: trexio_has_group_basis         = .false.
    logical :: trexio_has_group_determinant   = .false.
    logical :: trexio_has_group_ecp           = .false.

    character(:), allocatable  :: file_trexio_new, file_trexio_path

    private
    public :: dp

    public :: trexio_has_group_molecule
    public :: trexio_has_group_symmetry
    public :: trexio_has_group_orbitals
    public :: trexio_has_group_basis
    public :: trexio_has_group_determinant
    public :: trexio_has_group_ecp
    public :: file_trexio_new
    public :: file_trexio_path


#if defined(TREXIO_FOUND)
    public :: read_trexio_molecule_file
    public :: read_trexio_symmetry_file
    public :: read_trexio_orbitals_file
    public :: read_trexio_basis_file
    public :: update_trexio_orbitals
    public :: read_trexio_determinant_file
    public :: read_trexio_ecp_file
    public :: write_trexio_basis_num_info_file
#endif
    contains

#if defined(TREXIO_FOUND)
    subroutine read_trexio_molecule_file(file_trexio)
        !> This subroutine reads the .hdf5 trexio generated file/folder. It then computes the
        !! number of types of atoms, nuclear charges (from the symbol), and
        !! number of valence electrons if pseudopotential is provided.
        !! @author Ravindra Shinde (r.l.shinde@utwente.nl)
        !! @date 07 October 2021
        use custom_broadcast, only: bcast
        use mpiconf, only: wid
        use system, only: znuc, cent, iwctype, nctype, ncent, ncent_tot, nctype_tot, symbol, atomtyp
        use system, 		    only: newghostype, nghostcent
        use inputflags, only: igeometry
        use periodic_table, only: atom_t, element
        use contrl_file, only: ounit, errunit
        use general, only: pooldir
        use precision_kinds, only: dp
        use pseudo, only: nloc
        use system, only: nelec
        use system, only: nup
        use system, only: ndn
        use trexio
        use contrl_file, only: backend
        use multiple_geo, only: pecent

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
            trex_molecule_file = trexio_open(file_trexio_path, 'r', backend, rc)
            call trexio_assert(rc, TREXIO_SUCCESS)
            rc = trexio_read_nucleus_num(trex_molecule_file, ncent)
            if (trexio_has_nucleus(trex_molecule_file) == 0) trexio_has_group_molecule = .true.
            call trexio_assert(rc, TREXIO_SUCCESS)
            rc = trexio_read_electron_up_num(trex_molecule_file, nup)
            call trexio_assert(rc, TREXIO_SUCCESS)
            rc = trexio_read_electron_dn_num(trex_molecule_file, ndn)
            call trexio_assert(rc, TREXIO_SUCCESS)
        endif
        call bcast(trexio_has_group_molecule)
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
        rc = trexio_read_nucleus_coord(trex_molecule_file, cent)
        call trexio_assert(rc, TREXIO_SUCCESS)
        rc = trexio_read_nucleus_label(trex_molecule_file, symbol, 3)
        call trexio_assert(rc, TREXIO_SUCCESS)
        rc = trexio_close(trex_molecule_file)
        call trexio_assert(rc, TREXIO_SUCCESS)
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
            if (nloc == 0) then
                znuc(j) = atoms%znuclear
            else
                znuc(j) = atoms%znuclear - atoms%core
            endif
        enddo

        ncent_tot = ncent + nghostcent
        nctype_tot = nctype + newghostype

        write(ounit,*) '-----------------------------------------------------------------------'
        write(ounit,'(a, t15, a, t27, a, t39, a, t45, a)') 'Symbol', 'x', 'y', 'z', 'Type'
        write(ounit,'(t14, a, t26, a, t38, a )') '(bohr)', '(bohr)', '(bohr)'
        write(ounit,*) '-----------------------------------------------------------------------'

        do j= 1, ncent
            write(ounit,'(A4, 2x, 3F12.8, 2x, i3)') symbol(j), (cent(i,j),i=1,3), iwctype(j)
        enddo

        write(ounit,*) '-----------------------------------------------------------------------'
        write(ounit,*) " Values of znuc (number of valence electrons) "
        write(ounit,'(10F12.6)') (znuc(j), j = 1, nctype)
        write(ounit,*) '-----------------------------------------------------------------------'
        write(ounit,*)
    end subroutine read_trexio_molecule_file


    subroutine read_trexio_orbitals_file(file_trexio, build_only_basis)
        !> This subroutine reads the .hdf5 trexio generated file/folder. It then reads the
        !! number of molecular and atomic orbitals and their corresponding coefficients.
        !! @author Ravindra Shinde (r.l.shinde@utwente.nl)
        !! @date 12 October 2021
        use basis, only: ns, np, nd, nf, ng
        use custom_broadcast, only: bcast
        use mpiconf, only: wid
        use contrl_file, only: ounit, errunit
        use system, only: ncent, ncent_tot, iwctype, nctype_tot
        use system, only: newghostype, nghostcent
        use inputflags, only: ilcao
        use numbas, only: iwrwf, numr
        use custom_broadcast,   only: bcast
        use mpiconf,            only: wid
        use contrl_file,        only: ounit, errunit
        use coefs,              only: nbasis
        use csfs,               only: nstates
        use inputflags,         only: ilcao
        use numbas,             only: nrbas
        use numbas1,            only: iwlbas, nbastyp
        use orbval,             only: nadorb
        use pcm_fdc,            only: fs
        use vmc_mod,            only: norb_tot, nwftypeorb
        use multiple_geo,       only: nwftype
        use general,            only: pooldir
        use optwf_control,      only: method
        use precision_kinds, only: dp
        use slater,             only: coef
        use vmc_mod, only: nwftypeorb, nstoo, nstoomax, otos, extrao, nstoo_tot

        use error, only: trexio_error
        use trexio
        use contrl_file, only: backend
        use m_trexio_basis,     only: slm_per_l, index_slm, num_rad_per_cent
        use m_trexio_basis,     only: basis_num_shell, basis_shell_ang_mom
        use m_trexio_basis,     only: num_ao_per_cent, ao_radial_index
      use slater, only: norb

        implicit none

    !   local use
        character(len=72), intent(in)   :: file_trexio
        logical, intent(in), optional   :: build_only_basis
        character(len=40)               :: temp1, temp2
        character(len=120)              :: temp3, file_trexio_path
        integer                         :: iunit, iostat, iwft
        integer                         :: iorb, ibasis, i, j, k, l, ic, it
        integer                         :: counter, count1, count2
        integer                         :: index_ao, index_ao_copy
        integer                         :: lower_range, upper_range
        integer                         :: count
        integer                         :: cum_rad_per_cent, cum_ao_per_cent
        integer                         :: lower_rad_range, upper_rad_range
        logical                         :: exist
        logical                         :: skip = .true.

!       trexio
        integer, allocatable            :: basis_nucleus_index(:), ao_frequency(:), unique_index(:)
        integer, allocatable            :: res(:)

        !   Formatting
        character(len=100)               :: int_format     = '(A, T60, I0)'
        character(len=100)               :: string_format  = '(A, T60, A)'
        character(len=100)               :: float_format   = '(A, T60, f12.8)'

        ! trexio
        integer(8)                      :: trex_orbitals_file
        integer                         :: rc = 1, ii, jj, ind, index

        iwft = 1
        trex_orbitals_file = 0
        !   External file reading

        if((file_trexio(1:6) == '$pool/') .or. (file_trexio(1:6) == '$POOL/')) then
            file_trexio_path = pooldir // file_trexio(7:)
        else
            file_trexio_path = file_trexio
        endif

        if (.not. build_only_basis) then
        write(ounit,*) '---------------------------------------------------------------------------'
        write(ounit,*) " Reading LCAO orbitals from the file :: ",  trim(file_trexio_path)
        write(ounit,*) '---------------------------------------------------------------------------'
        endif
        ! Check if the file exists

        if (wid) then
            trex_orbitals_file = trexio_open(file_trexio_path, 'r', backend, rc)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio file open error', __FILE__, __LINE__)
            rc = trexio_read_mo_num(trex_orbitals_file, norb_tot)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_mo_num', __FILE__, __LINE__)
            rc = trexio_read_ao_num(trex_orbitals_file, nbasis)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_ao_num', __FILE__, __LINE__)
            rc = trexio_read_basis_shell_num(trex_orbitals_file, basis_num_shell)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_basis_shell_num', __FILE__, __LINE__)
        endif
        call bcast(norb_tot)
        call bcast(nbasis)
        call bcast(basis_num_shell)

        norb = norb_tot        ! norb will get updated later. norb_tot is fixed


        ! Do the array allocations
        if( (method(1:3) == 'lin')) then
            if (.not. allocated(coef)) allocate (coef(nbasis, norb_tot, 3))
        elseif(method == 'sr_n') then
            if (.not. allocated(coef)) allocate (coef(nbasis, norb_tot, nstates))
            nwftypeorb=nstates
            nstoo_tot=nstates
            nstoomax=1
            extrao=1
            if (.not. allocated(nstoo)) allocate (nstoo(nwftypeorb))
            if (.not. allocated(otos)) allocate (otos(nwftypeorb,1))
            do i=1,nstates
               nstoo(i)=1
               otos(i,1)=i
            enddo
        else
            if (.not. allocated(coef)) allocate (coef(nbasis, norb_tot, nwftype))
        endif

        ! Do the allocations based on the number of shells and primitives
        if (.not. allocated(basis_nucleus_index))    allocate(basis_nucleus_index(basis_num_shell))
        if (.not. allocated(basis_shell_ang_mom))    allocate(basis_shell_ang_mom(basis_num_shell))
        if (.not. allocated(ao_radial_index))        allocate(ao_radial_index(nbasis))
        if (.not. allocated(index_slm))              allocate(index_slm(nbasis))
        if (.not. allocated(num_rad_per_cent))       allocate(num_rad_per_cent(ncent_tot))
        if (.not. allocated(num_ao_per_cent))        allocate(num_ao_per_cent(ncent_tot))

        if (.not. build_only_basis) then
        ! Read the orbitals
        if (wid) then
            if (trexio_has_mo(trex_orbitals_file) == 0) trexio_has_group_orbitals = .true.
            rc = trexio_read_mo_coefficient(trex_orbitals_file, coef(:,:,1))
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_mo_coeffs', __FILE__, __LINE__)
        endif
        call bcast(trexio_has_group_orbitals)
        call bcast(coef(:,:,1))
        endif

        ! Make a copy of orbital coeffs for multiple states
        if( (method == 'sr_n') .and. (nstates .gt. 1)) then
            do i=2,nstates
              coef(:,:,i)=coef(:,:,1)
            enddo
        endif

!   Generate the basis information (which radial to be read for which Slm)
        if (wid) then
            rc = trexio_read_basis_shell_ang_mom(trex_orbitals_file, basis_shell_ang_mom)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_basis_shell_ang_mom', __FILE__, __LINE__)
            rc = trexio_read_basis_nucleus_index(trex_orbitals_file, basis_nucleus_index)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_basis_nucleus_index', __FILE__, __LINE__)
        endif
        call bcast(basis_shell_ang_mom)
        call bcast(basis_nucleus_index)

!         if (wid) then
! #if defined(TREXIO_FOUND)
!             rc = trexio_close(trex_orbitals_file)
!             call trexio_error(rc, TREXIO_SUCCESS, 'trexio_close trex_orbital_file', __FILE__, __LINE__)
! #endif
!         endif


        numr = 1            ! Debug Check this statement. Not sure how to store multiple bfinfo files in single trexio

        ! Generate the index of the slm for each AO
        ! i.e. index_slm(i) = the slm index of the i-th AO
        ! The list follows the following order:

        !   l   |   1  2  3  4    5   6   7   8   9   10
        ! ------+-----------------------------------------
        !   y   |   s  x  y  z    xx  xy  xz  yy  yz  zz

        !           11  12  13  14  15  16  17  18  19  20
        !       ------------------------------------------
        !           xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
        !
        !          21   22   23   24   25   26   27   28   29   30   31   32   33   34   35
        !       +-----------------------------------------------------------------------------
        !          xxxx xxxy xxxz xxyy xxyz xxzz xyyy xyyz xyzz xzzz yyyy yyyz yyzz yzzz zzzz

        counter = 0; count1 = 1; count2 = 0
        cum_rad_per_cent = 0
        cum_ao_per_cent  = 0
        index_ao = 0; jj = 1; index_ao_copy = 0
        ! The following loop will generate the index_slm array which tells
        ! which AO is of which type (from the above list)
        do l = 1, basis_num_shell
            k = basis_shell_ang_mom(l)

            counter = counter + slm_per_l(k+1)
            count2 = 0;
            do ii = 1, slm_per_l(k+1)

                index_ao = index_ao + 1
                index_slm(index_ao) = sum(slm_per_l(1:k)) + count2 + 1

                count2 = count2 + 1

                cum_ao_per_cent = cum_ao_per_cent + 1
            end do

            jj = jj + 1

            cum_rad_per_cent = cum_rad_per_cent + 1

            ! The following if loop is for counting the number of radial functions
            ! and number of AOs per center
            if (count1 == basis_nucleus_index(l)) then
                num_rad_per_cent(count1) = cum_rad_per_cent
                do ii = 1, slm_per_l(k+1)
                    index_ao_copy = index_ao_copy + 1
                    ao_radial_index(index_ao_copy) = cum_rad_per_cent
                enddo
                num_ao_per_cent(count1) = cum_ao_per_cent
            else
                cum_rad_per_cent = 1
                index_ao_copy = index_ao_copy + 1
                ao_radial_index(index_ao_copy) = cum_rad_per_cent
                cum_ao_per_cent  = 1
                count1 = count1 + 1
            end if
        enddo ! loop on shells

        if (.not. allocated(nbastyp)) allocate (nbastyp(nctype_tot))

        if (.not. allocated(ns)) allocate (ns(nctype_tot))
        if (.not. allocated(np)) allocate (np(nctype_tot))
        if (.not. allocated(nd)) allocate (nd(nctype_tot))
        if (.not. allocated(nf)) allocate (nf(nctype_tot))
        if (.not. allocated(ng)) allocate (ng(nctype_tot))


        ! Obtain the index of radials for each unique center (iwrwf)
        if (.not. allocated(iwlbas)) allocate (iwlbas(nbasis, nctype_tot))
        if (.not. allocated(nrbas)) allocate (nrbas(nctype_tot), source=0)
        if (.not. allocated(iwrwf))  allocate (iwrwf(nbasis, nctype_tot))
        if (.not. allocated(ao_frequency)) allocate (ao_frequency(nbasis), source=0)    ! ao upto g orbitals
        if (.not. allocated(unique_index)) allocate (unique_index(nbasis), source=0)    ! ao upto g orbitals
        if (.not. allocated(res)) allocate (res(basis_num_shell), source=0)                       ! shells upto g orbitals

        do i = 1, ncent_tot
            if (numr .gt. 0) then
                nbastyp(iwctype(i)) = num_ao_per_cent(i)
                nrbas(iwctype(i)) = num_rad_per_cent(i)
            endif
        enddo

        ! Generate iwlbas and iwrwf
        lower_range = 1; lower_rad_range = 1
        do ic=1,ncent+nghostcent
            it=iwctype(ic)
            upper_range = lower_range + num_ao_per_cent(ic) - 1
            upper_rad_range = lower_rad_range + num_rad_per_cent(ic) - 1
            iwlbas(1:num_ao_per_cent(ic), it) = index_slm(lower_range:upper_range)
            iwrwf(1:num_ao_per_cent(ic), it) = ao_radial_index(lower_range:upper_range)

            ! The following block of code is for generating the ns,np,nf,nd, and ng arrays
            call unique_elements(num_rad_per_cent(ic), basis_shell_ang_mom(lower_rad_range:upper_rad_range), &
                                 res, count, ao_frequency(1:num_rad_per_cent(ic)), &
                                 unique_index(1:num_rad_per_cent(ic)))
            ns(iwctype(ic)) = ao_frequency(1)
            np(iwctype(ic)) = ao_frequency(2)
            nd(iwctype(ic)) = ao_frequency(3)
            nf(iwctype(ic)) = ao_frequency(4)
            ng(iwctype(ic)) = ao_frequency(5)
            ! End of block of code
            lower_range = upper_range + 1
            lower_rad_range = upper_rad_range + 1
        enddo

        if (allocated(basis_nucleus_index)) deallocate(basis_nucleus_index)
        if (allocated(ao_frequency)) deallocate(ao_frequency)
        if (allocated(unique_index)) deallocate(unique_index)
        if (allocated(res)) deallocate(res)

        if (.not. build_only_basis) then
        write(ounit,int_format) " Number of basis functions ", nbasis
        write(ounit,int_format) " Number of lcao orbitals ", norb
        write(ounit,int_format) " Type of wave functions ", iwft
        write(ounit,*) "Orbital coefficients are written to the output.log file"




        write(ounit,*)
        ilcao = ilcao + 1
        write(ounit,*) "----------------------------------------------------------"
        endif

    end subroutine read_trexio_orbitals_file

    subroutine write_trexio_basis_num_info_file(file_trexio)
        !> This subroutine prints the generated information about the numerical basis.
        !> The data is already generated in the read_trexio_orbitals_file.
        !! @author Ravindra Shinde (r.l.shinde@utwente.nl)
        !! @date 06 July 2022
        use custom_broadcast,   only: bcast
        use mpiconf,            only: wid
        use contrl_file,        only: ounit, errunit
        use system,               only: nctype_tot
        use inputflags,         only: ibasis_num
        use numbas,             only: iwrwf, numr
        use numbas1,            only: iwlbas, nbastyp

        use basis,              only: ns, np, nd, nf, ng

        use multiple_geo,              only: nwftype

        implicit none

    !   local use
        character(len=72), intent(in)   :: file_trexio
        integer                         :: i, ib
        logical                         :: build_only_basis

        write(ounit,*) '---------------------------------------------------------------------------'
        write(ounit,'(a)')  " Basis function types and pointers to radial parts tables generated from :: " &
                            // trim(file_trexio)
        write(ounit,*) '---------------------------------------------------------------------------'

        call read_trexio_orbitals_file(file_trexio, build_only_basis=.true.)

        do i = 1, nctype_tot
            write (ounit, '(100i3)') ns(i), np(i), nd(i), nf(i), ng(i)

            if (numr .gt. 0) then
                write(ounit, '(100i3)') (iwlbas(ib, i), ib=1, nbastyp(i))
                write(ounit, '(100i3)') (iwrwf(ib, i), ib=1, nbastyp(i))
            endif
            write(ounit,*)
        enddo

        write(ounit,*)
        ibasis_num =  1
        write(ounit,*) "----------------------------------------------------------"

    end subroutine write_trexio_basis_num_info_file


    subroutine update_trexio_orbitals
        !> This subroutine updates the .hdf5 trexio generated file/folder. It then reads the
        !! number of molecular and atomic orbitals and their corresponding coefficients.
        !! @author Ravindra Shinde (r.l.shinde@utwente.nl)
        !! @date 08 November 2023
        use mpiconf,            only: wid
        use custom_broadcast,   only: bcast
        use contrl_file,        only: ounit, errunit
        use mpiconf,            only: wid
        use contrl_file,        only: ounit, errunit
        use coefs,              only: nbasis
        use csfs,               only: nstates
        use inputflags,         only: scalecoef
        use vmc_mod,            only: norb_tot, nwftypeorb
        use mpiconf,            only: wid
        use multiple_geo,       only: nwftype
        use optwf_control,      only: ioptorb
        use slater,             only: coef
        use general,            only: pooldir
        use precision_kinds,    only: dp
        use error,              only: trexio_error
        use trexio
        use contrl_file,        only: backend
        use slater,             only: norb


#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND)
        use qmckl_data
#endif

        implicit none

    !   local use
        character(len=120)              :: file_trexio_path

    ! trexio
        integer(8)                      :: trex_orbitals_file
        integer                         :: k, rc, ierr

        if(ioptorb.eq.0) return

        trex_orbitals_file = 0

        !   External file reading
#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND)
        if((file_trexio_new(1:6) == '$pool/') .or. (file_trexio_new(1:6) == '$POOL/')) then
            file_trexio_path = pooldir // file_trexio_new(7:)
        else
            file_trexio_path = file_trexio_new
        endif

        write(ounit,*) '---------------------------------------------------------------------------'
        write(ounit,*) " Updating LCAO orbitals from the file :: ",  trim(file_trexio_path)
        write(ounit,*) '---------------------------------------------------------------------------'

        ! Destroy the existing QMCkl context first
        write(ounit, *) " QMCkl destroying the old context " , qmckl_ctx , " successfully "
        rc = qmckl_context_destroy(qmckl_ctx)
        call trexio_error(rc, TREXIO_SUCCESS, 'trexio_close trex_update_mo', __FILE__, __LINE__)
        write(ounit, '(a)') " QMCkl old context destroyed successfully "

        if(wid) then

        ! Check if the file exists
            trex_orbitals_file = trexio_open(file_trexio_path, 'u', backend, rc)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio file open error', __FILE__, __LINE__)
            write(ounit, '(a)') "File opened successfully "

        ! Update the orbitals
            do k=1,nwftypeorb
              rc = trexio_write_mo_coefficient(trex_orbitals_file, coef(:,:,k)/scalecoef)
              call trexio_error(rc, TREXIO_SUCCESS, 'trexio_write_mo_coeffs', __FILE__, __LINE__)
              write(ounit, '(a)') " MO coeffs updated successfully "
            enddo

        ! Close the file
            rc = trexio_close(trex_orbitals_file)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_close trex_update_mo', __FILE__, __LINE__)
            write(ounit, '(a)') "File " // trim(file_trexio_path) // " closed successfully "
        endif
        call MPI_Barrier( MPI_COMM_WORLD, ierr )

        ! Create a new QMCkl context with the new trexio file
        qmckl_ctx = qmckl_context_create()
        rc = qmckl_trexio_read(qmckl_ctx, file_trexio_path, 1_8*len(trim(file_trexio_path)))
        call trexio_error(rc, TREXIO_SUCCESS, 'INPUT: QMCkl error: Unable to read TREXIO file', __FILE__, __LINE__)
        write(ounit, *) " QMCkl new context created  ", qmckl_ctx,  " successfully "

#endif

        write(ounit,*) "----------------------------------------------------------"

    end subroutine update_trexio_orbitals


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
        use system,               only: znuc, nctype, nctype_tot, ncent_tot
        use system,               only: symbol, atomtyp
        use vmc_mod,            only: NCOEF
        use system,          only: newghostype
        use control,            only: ipr
        use numbas,             only: arg, d2rwf, igrid, nr, nrbas, r0, rwf, rmaxwf
        use numbas,             only: allocate_numbas
        use coefs,              only: nbasis
        use numexp,             only: ae, ce, allocate_numexp
        use pseudo,             only: nloc
        use general,            only: filename, filenames_bas_num

        ! For processing the stored information
        use system,             only: atomtyp
        use general,            only: pooldir, bas_id
        use contrl_file,        only: ounit, errunit
        use spline2_mod,        only: spline2
        use fitting_methods,    only: exp_fit

        use trexio
        use contrl_file,        only: backend
        use error,              only: trexio_error
        use m_trexio_basis,     only: gnorm

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

        ! for local use.
        character(len=72), intent(in)   :: file_trexio
        character(len=128)              :: file_trexio_path
        integer                         :: i, j, k, tcount1, tcount2, tcount3, tcount4
        integer                         :: counter, counter_shell, lower_shell, upper_shell, lower_prim, upper_prim

        ! trexio
        integer(8)                      :: trex_basis_file
        integer                         :: rc = 1

        !   Formatting
        character(len=128)              :: int_format     = '(A, T60, I0)'

        ! Grid related
        integer                         :: gridtype=3
        integer                         :: gridpoints=2000
        real(dp)                        :: gridarg=1.003d0
        real(dp)                        :: gridr0=20.0d0
        real(dp)                        :: gridr0_save = 20.0d0
        real(kind=dp), dimension(2000)  :: rgrid  ! Grid points
        integer, dimension(nctype_tot)  :: icusp
        real(dp)                        :: r, r2

        integer, dimension(:), allocatable :: atom_index(:), shell_index_atom(:), nshells_per_atom(:)
        integer, dimension(:), allocatable :: prim_index_atom(:), nprims_per_atom(:)
        integer, dimension(:), allocatable :: unique_atom_index(:), shell_prim_correspondence(:)
        integer                         :: count
        character(len=2), allocatable   :: unique(:) ! unique symbols of atoms

        ! needed for spline
        real(dp), dimension(MRWF_PTS)       ::  x, work
        real(dp), dimension(ncoef)          ::  y
        real(dp), dimension(ncoef*ncoef)    ::  dmatr
        real(dp), dimension(nbasis)         ::  l
        integer, dimension(ncoef)           :: ipiv
        integer         :: ic, ir, irb, ii, jj, ll, icoef, iff
        integer         :: iwf = 1
        integer         :: info
        real(dp)        :: val, dwf1, wfm, dwfn, dwfm, temp

        ! rmax cutoff
        real(dp)                            :: cutoff_rmax = 1.0d-12

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
            trex_basis_file = trexio_open(file_trexio_path, 'r', backend, rc)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio file open error', __FILE__, __LINE__)
            rc = trexio_read_basis_shell_num(trex_basis_file, basis_num_shell)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_basis_shell_num', __FILE__, __LINE__)
            rc = trexio_read_ao_num(trex_basis_file, ao_num)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_ao_num', __FILE__, __LINE__)
            rc = trexio_read_basis_prim_num(trex_basis_file, basis_num_prim)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_basis_prim_num', __FILE__, __LINE__)
            if (trexio_has_basis(trex_basis_file) == 0) trexio_has_group_basis = .true.
        endif
        call bcast(trexio_has_group_basis)
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
            trex_basis_file = trexio_open(file_trexio_path, 'r', backend, rc)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio file open error', __FILE__, __LINE__)
            rc = trexio_read_basis_nucleus_index(trex_basis_file, basis_nucleus_index)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_basis_nucleus_index', __FILE__, __LINE__)
            rc = trexio_read_basis_shell_index(trex_basis_file, basis_shell_index)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_basis_shell_index', __FILE__, __LINE__)
            rc = trexio_read_basis_shell_ang_mom(trex_basis_file, basis_shell_ang_mom)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_basis_shell_ang_mom', __FILE__, __LINE__)
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

        if (wid) then
            rc = trexio_close(trex_basis_file)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_close trex_basis_file', __FILE__, __LINE__)
        endif

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
        allocate(nshells_per_atom(basis_num_shell), source=0)
        allocate(shell_index_atom(basis_num_shell))

        call unique_elements(basis_num_shell, basis_nucleus_index, atom_index, count, nshells_per_atom, shell_index_atom)

        ! count                     :: "Number of unique elements"
        ! atom_index(1:count)       :: "Unique elements index (not used here)"
        ! nshells_per_atom(1:count) :: "frequency or count of shells per atom"
        ! shell_index_atom(1:count) :: "index number of the shell for each atom"
        ! The shells per atom can be obtained by accessing the shell_index_atom
        ! for a given atom index by the slice of size frequency.


        ! Now get the number of primitives per atom and their indices
        ! also obtain the number of primitives per shell for all the atoms
        allocate(nprims_per_atom(ncent_tot), source=0)
        allocate(prim_index_atom(ncent_tot))
        allocate(shell_prim_correspondence(basis_num_shell))

        prim_index_atom(1) = 1
        nprims_per_atom(1) = 0
        tcount4 = 0; tcount2 = 0
        do i = 1, ncent_tot
            tcount3 = 0
            do j = 1, nshells_per_atom(i)   ! frequency
                tcount4 = tcount4 + 1
                tcount3 = tcount3 + 1
                tcount1 = 0
                do k = 1, basis_num_prim
                    if (tcount4 == basis_shell_index(k)) then
                        tcount1 = tcount1 + 1
                        shell_prim_correspondence(tcount4) = tcount1
                    endif
                enddo

                nprims_per_atom(i) = nprims_per_atom(i) + shell_prim_correspondence(tcount4)
            enddo
            if (i .ne. ncent_tot) prim_index_atom(i+1) = prim_index_atom(i) + nprims_per_atom(i)
        enddo

        ! Obtain the number of unique types of atoms stored in the hdf5 file.
        if (.not. allocated(unique)) allocate(unique(nctype_tot))
        if (.not. allocated(unique_atom_index)) allocate(unique_atom_index(nctype_tot))


        tcount1 = 1; tcount2 = 1
        unique_atom_index(1) = 1
        unique(1) = symbol(1)
        do j= 2, ncent_tot
            if (any(unique == symbol(j) ))  then
                cycle
            endif
            tcount1 = tcount1 + 1
            unique_atom_index(tcount1) = j
            unique(tcount1) = symbol(j)
        enddo

        ! start putting in the information in the arrays and variables
        gridtype=3
        gridpoints=2000
        gridarg=1.003d0
        gridr0=20.0d0
        gridr0_save = gridr0

        ! Do the necessary allocation for the numerical basis set
        call allocate_numbas()
        call allocate_numexp()

        if (gridtype .eq. 3) gridr0 = gridr0/(gridarg**(gridpoints-1)-1)

        ! Populate the rgrid array.
        do i = 1, gridpoints
            if (gridtype .eq. 1) then
                rgrid(i) = gridr0 + (i-1)*gridarg
            else if (gridtype .eq. 2) then
                rgrid(i) = gridr0 * gridarg**(i-1)
            else if (gridtype .eq. 3) then
                rgrid(i) = gridr0*gridarg**(i-1) - gridr0
            endif
            x(i) = rgrid(i)
        enddo


        do ic = 1, nctype + newghostype            ! loop over all the unique atoms

            nrbas(ic)   = nshells_per_atom(unique_atom_index(ic))
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

            ! All electron case when nloc equals zero.
            ! Make sure that the input file has nloc=0.
            if(nloc.eq.0) then
                do irb = 1, nrbas(ic)
                    l(irb) = 0
                enddo
            endif


            ! loop over all the primitives for the unique atom
            ! The lower and upper indices of primitive indices
            lower_shell = shell_index_atom(unique_atom_index(ic))
            upper_shell = shell_index_atom(unique_atom_index(ic)) + nshells_per_atom(unique_atom_index(ic)) - 1

            lower_prim = prim_index_atom(unique_atom_index(ic))
            upper_prim = prim_index_atom(unique_atom_index(ic)) + nprims_per_atom(unique_atom_index(ic)) - 1


            ! select the shells corresponding to the unique atoms only
            ! j is the running shell index for the unique atom i
            do i = 1, nr(ic)
                r = rgrid(i)
                r2 = r*r

                counter = lower_prim
                counter_shell = 1
                do j = lower_shell, upper_shell
                    ! loop on primitives in the given shell
                    val = 0.0d0
                    do k = counter, counter + shell_prim_correspondence(j) -1
                        val = val + gnorm(basis_exponent(k), basis_shell_ang_mom(j)) &
                                  * basis_coefficient(k) * dexp(-basis_exponent(k)*r2)
                    enddo
                    counter = counter + shell_prim_correspondence(j)
                    rwf(i,counter_shell,ic,1) = val
                    counter_shell = counter_shell + 1
                enddo

            enddo

            do irb=1,nrbas(ic)

                if(nloc.eq.0.and.l(irb).eq.0.and.icusp(ic).eq.1) then

        ! c small radii wf(r)=ce1-znuc*ce1*r+ce3*r**2+ce4*r**3+ce5*r**4
                do ii=1,NCOEF-1
                    dmatr(ii)=1.d0-znuc(ic)*x(ii)
                enddo
                y(1)=rwf(1,irb,ic,iwf)
                ll=NCOEF-1
                do jj=2,NCOEF-1
                    y(jj)=rwf(jj,irb,ic,iwf)
                    do ii=2,NCOEF-1
                    ll=ll+1
                    dmatr(ll)=x(ii)**jj
                    enddo
                enddo

                call dgesv(NCOEF-1,1,dmatr,NCOEF-1,ipiv,y,NCOEF,info)
                ce(1,irb,ic,iwf)=y(1)
                ce(2,irb,ic,iwf)=-znuc(ic)*ce(1,irb,ic,iwf)
                ce(3,irb,ic,iwf)=y(2)
                ce(4,irb,ic,iwf)=y(3)
                ce(5,irb,ic,iwf)=y(4)
                else

        ! c small radii wf(r)=ce1+ce2*r+ce3*r**2+ce4*r**3+ce5*r**4
                ll=0
                do jj=1,NCOEF
                    y(jj)=rwf(jj,irb,ic,iwf)
                    do ii=1,NCOEF
                    ll=ll+1
                    dmatr(ll)=x(ii)**(jj-1)
                    enddo
                enddo
                call dgesv(NCOEF,1,dmatr,NCOEF,ipiv,y,NCOEF,info)

                do icoef=1,NCOEF
                    ce(icoef,irb,ic,iwf)=y(icoef)
                enddo
                endif



        ! c       if(ipr.gt.1) then
                if (wid) then
                    write(45,'(''basis = '',i4)') irb
                    write(45,'(''check the small radius expansion'')')
                    write(45,'(''coefficients'',1p10e22.10)') &
                                (ce(iff,irb,ic,iwf),iff=1,NCOEF)
                    write(45,'(''check the small radius expansion'')')
                    write(45,'(''irad, rad, extrapolated value, correct value'')')
                    do ir=1,10
                        val=ce(1,irb,ic,iwf)
                        do icoef=2,NCOEF
                        val=val+ce(icoef,irb,ic,iwf)*x(ir)**(icoef-1)
                        enddo
                        write(45,'(i2,1p3e22.14)')ir,x(ir),val,rwf(ir,irb,ic,iwf)
                    enddo
                endif
        ! c       endif

                dwf1=0.d0
                do icoef=2,NCOEF
                dwf1=dwf1+(icoef-1)*ce(icoef,irb,ic,iwf)*x(1)**(icoef-2)
                enddo

        ! Update the rmax at the point where rwf goes below cutoff (scanning from right to left)
                rmaxwf(irb, ic) = x(nr(ic))
                rloop: do ir=nr(ic),1,-1
                  if (dabs(rwf(ir,irb,ic,iwf)) .gt. cutoff_rmax ) then
                    rmaxwf(irb, ic) = x(ir)
                    exit rloop
                  endif
                enddo rloop

                if (wid) write(45,'(a,i0,a,i0,a,g0)') "Initial rmax for center = ",ic, " basis = ",irb, " is ", rmaxwf(irb, ic)

! Nonzero basis at the boundary : Do exponential fitting.
                if(dabs(rmaxwf(irb,ic)-x(nr(ic))).lt.1.0d-10) then
                    call exp_fit(x(nr(ic)-9:nr(ic)),rwf(nr(ic)-9:nr(ic),irb,ic,iwf), 10, ae(1,irb,ic,iwf), ae(2,irb,ic,iwf))
                    if(ae(2,irb,ic,iwf).lt.0) call fatal_error ('BASIS_READ_NUM: ak<0')

                    rmaxwf(irb,ic)=-dlog(cutoff_rmax/dabs(ae(1,irb,ic,iwf)))/ae(2,irb,ic,iwf)

                    if (wid) then
                        write(45,'(a)') 'check the large radius expansion'
                        write(45,'(a,g0,2x,g0)') 'Exponential fitting parameters : ', ae(1,irb,ic,iwf), ae(2,irb,ic,iwf)

                        write(45,'(a,i0,a,i0,a,g0)') "Final rmax (fit) for center = ",ic, " basis = ",irb, " is ", rmaxwf(irb, ic)
                        write(45, '(a)') 'irad,         rad                  rwf value            expo fit'
                        do ir=1,10
                        temp = ae(1,irb,ic,iwf)*dexp(-ae(2,irb,ic,iwf)*x(nr(ic)-ir))
                        write(45,'(i3,2x,1p4e22.14)') ir,x(nr(ic)-ir),rwf(nr(ic)-ir,irb,ic,iwf), temp
                        enddo
                    endif

                    dwfn=-ae(2,irb,ic,iwf)*ae(1,irb,ic,iwf)*dexp(-ae(2,irb,ic,iwf)*x(nr(ic)))

                else
                    dwfn=0.d0
                endif
                if (wid) write(45,*) 'dwf1,dwfn',dwf1,dwfn

                call spline2(x,rwf(1,irb,ic,iwf),nr(ic),dwf1,dwfn, d2rwf(1,irb,ic,iwf), work)

            enddo ! loop on irb : number of radial shells
        enddo ! loop on ic : the unique atom types

        ! Do the deallocations of local arrays
        if (allocated(unique)) deallocate(unique)
        if (allocated(unique_atom_index)) deallocate(unique_atom_index)

        if (allocated(nprims_per_atom)) deallocate(nprims_per_atom)
        if (allocated(prim_index_atom)) deallocate(prim_index_atom)
        if (allocated(shell_prim_correspondence)) deallocate(shell_prim_correspondence)

        if (allocated(atom_index)) deallocate(atom_index)
        if (allocated(nshells_per_atom)) deallocate(nshells_per_atom)
        if (allocated(shell_index_atom)) deallocate(shell_index_atom)

        if (allocated(basis_nucleus_index))    deallocate(basis_nucleus_index)
        if (allocated(basis_shell_index))      deallocate(basis_shell_index)
        if (allocated(basis_shell_ang_mom))    deallocate(basis_shell_ang_mom)
        if (allocated(basis_shell_factor))     deallocate(basis_shell_factor)
        if (allocated(basis_exponent))         deallocate(basis_exponent)
        if (allocated(basis_coefficient))      deallocate(basis_coefficient)
        if (allocated(basis_prim_factor))      deallocate(basis_prim_factor)
        if (allocated(ao_shell))               deallocate(ao_shell)
        if (allocated(ao_normalization))       deallocate(ao_normalization)

    end subroutine read_trexio_basis_file



    subroutine read_trexio_symmetry_file(file_trexio)
        ! Ravindra

        use custom_broadcast,   only: bcast
        use mpiconf,            only: wid, idtask

        use contrl_file,        only: ounit, errunit
        use slater,             only: norb
        use optorb,             only: irrep
        use vmc_mod,            only: norb_tot
        use general,            only: pooldir
        use precision_kinds,    only: dp
        use array_utils,        only: unique_string_elements


        use trexio
        use contrl_file,        only: backend
        use error,              only: trexio_error

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
            trex_symmetry_file = trexio_open(file_trexio_path, 'r', backend, rc)
            rc = trexio_read_mo_num(trex_symmetry_file, mo_num)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_mo_num failed', __FILE__, __LINE__)
            if (trexio_has_mo_symmetry(trex_symmetry_file) == 0) trexio_has_group_symmetry = .true.
        endif
        call bcast(mo_num)
        call bcast(trexio_has_group_symmetry)
        ! safe allocate
        if (.not. allocated(irrep)) allocate (irrep(mo_num))
        if (.not. allocated(mo_symmetry)) allocate (mo_symmetry(mo_num))
        if (.not. allocated(unique_irrep)) allocate (unique_irrep(mo_num))


        if (trexio_has_group_symmetry) then
            if (wid) then
                rc = trexio_read_mo_symmetry(trex_symmetry_file, mo_symmetry, 2)
                call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_mo_symmetry failed', __FILE__, __LINE__)
            endif
            call bcast(mo_symmetry)
        else ! set default symmetry if information not present in the file
            mo_symmetry(:) = 'A'
            irrep(:) = 1
        endif


        write(ounit,fmt=int_format) " Number of molecular orbital symmetries read ::  ", mo_num

        call unique_string_elements(mo_num, mo_symmetry, unique_irrep, num_irrep)


        write(ounit,fmt=int_format) " Number of irreducible representations       ::  ", num_irrep
        write(ounit,*)
        write(ounit,'(a)')          " Irreducible representations correspondence  ::  "
        write(ounit,'(1x,10(a2,a,i2,a,1x))') (unique_irrep(i), "=", i ,";", i=1, num_irrep)
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
        use contrl_file,        only: ounit, errunit
        use general,            only: pooldir
        use slater,             only: cdet, ndet
        use dorb_m,             only: iworbd
        use slater,             only: norb
        use inputflags,         only: ideterminants
        use multiple_geo,              only: nwftype
        use csfs,               only: nstates
        use mstates_mod,        only: MSTATES
        use general,            only: pooldir
        use system,             only: ndn, nup
        use system,             only: nelec
        use optwf_control,      only: method
        use precision_kinds,    only: dp

        use trexio
        use contrl_file,        only: backend
        use error,              only: trexio_error

        implicit none

        !   local use
        character(len=72), intent(in)   :: file_trexio
        character(len=40)               :: temp
        character(len=80)               :: comment, file_trexio_path
        integer                         :: iostat, i, j, k, iunit, jj
        logical                         :: exist
        character(len=2), allocatable   :: unique(:)

        ! trexio
        integer(8)                      :: trex_determinant_file
        integer                         :: rc = 1

        !   Formatting
        character(len=100)              :: int_format     = '(A, T60, I0)'
        character(len=100)              :: float_format   = '(A, T60, f12.8)'
        character(len=100)              :: string_format  = '(A, T60, A)'

        ! determinant data
        integer*8, allocatable          :: buffer(:,:,:)
        integer(8)                      :: offset, icount, BUFSIZE
        integer                         :: int64_num, m           ! Number of intergers required per spin component
        integer*8                       :: determinant_num
        integer*4, allocatable          :: orb_list_up(:), orb_list_dn(:)
        integer*4                       :: occ_num_up, occ_num_dn, occupied_num


        trex_determinant_file = 0

        !   External file reading

        if((file_trexio(1:6) == '$pool/') .or. (file_trexio(1:6) == '$POOL/')) then
            file_trexio_path = pooldir // file_trexio(7:)
        else
            file_trexio_path = file_trexio
        endif

        write(ounit,*) '-----------------------------------------------------------------------'
        write(ounit,*) " Reading determinants from the trexio file :: ",  file_trexio_path
        write(ounit,*) '-----------------------------------------------------------------------'

        ! Check if the file exists
        if (wid) then
            trex_determinant_file = trexio_open(file_trexio_path, 'r', backend, rc)
            call trexio_assert(rc, TREXIO_SUCCESS)
            rc = trexio_has_determinant (trex_determinant_file)
            if (rc == TREXIO_SUCCESS) then
                rc = trexio_read_determinant_num(trex_determinant_file, ndet)
                call trexio_assert(rc, TREXIO_SUCCESS)
                rc = trexio_get_int64_num(trex_determinant_file, int64_num)
                call trexio_assert(rc, TREXIO_SUCCESS)
            else
                write(errunit,*) "trexio file does not have number of determinant  stored :: ", rc
                call trexio_error(rc, TREXIO_SUCCESS, 'trexio_has_group_determinant failed', __FILE__, __LINE__)
            endif
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

        BUFSIZE = determinant_num
        offset = 0_8

        allocate(buffer(int64_num, 2, BUFSIZE))
        allocate(orb_list_up(int64_num*64), orb_list_dn(int64_num*64))


        if (wid) then
        rc = trexio_read_determinant_coefficient(trex_determinant_file, offset, BUFSIZE, cdet(:,1,nwftype))
        if (trexio_has_determinant_coefficient(trex_determinant_file) == 0) trexio_has_group_determinant = .true.
        call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_determinant_coeff failed', __FILE__, __LINE__)
        endif
        call bcast(trexio_has_group_determinant)
        call bcast(cdet)


        write(ounit,*)
        write(ounit,*) " Determinant coefficients "
        write(ounit,'(10(1x, f11.8, 1x))') (cdet(i,1,nwftype), i=1, BUFSIZE)

!       allocate the orbital mapping array
        if (.not. allocated(iworbd)) allocate(iworbd(nelec, determinant_num))

        write(ounit, *)
        write(ounit, *) "Orbitals <--> Determinants mapping read from a trexio file :: "
        write(ounit, *) "Serial numbers of orbitals that are occupied               :: "
        write(ounit, *) "'alpha (spin up)'  <---------------------->  'beta (spin down)' "
        write(ounit, *)

        offset = 0_8
        icount = BUFSIZE

        if (wid) then
            do while (icount == BUFSIZE)
                if (offset < ndet) then
                    rc = trexio_read_determinant_list(trex_determinant_file, offset, icount, buffer)
                    call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_determinant_list failed', __FILE__, __LINE__)
                    offset = offset + icount
                else
                    icount = 0
                end if

                do m=1,icount
                    rc = trexio_to_orbital_list_up_dn(int64_num, buffer(1,1,m), orb_list_up, orb_list_dn, occ_num_up, occ_num_dn)
                    call trexio_error(rc, TREXIO_SUCCESS, 'trexio_to_orbital_list_up_dn failed', __FILE__, __LINE__)
                    write(temp, '(1x,a,i0,a,i0,a)') '(', occ_num_up, '(i4,1x),', occ_num_dn, '(i4,1x))'
                    write(ounit, temp) (orb_list_up(i), i = 1, occ_num_up), (orb_list_dn(i), i = 1, occ_num_dn)

                    do i = 1, occ_num_up
                        iworbd(i, m) = orb_list_up(i)
                    enddo
                    do i = 1, occ_num_dn
                        iworbd(occ_num_up + i, m) = orb_list_dn(i)
                    enddo

                end do
            end do

            deallocate(buffer)
            deallocate(orb_list_up)
            deallocate(orb_list_dn)
        endif
        call bcast(iworbd)

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

        use trexio
        use contrl_file,        only: backend
        use error,              only: trexio_error

        use pseudo_mod,         only: MPS_L, MGAUSS, MPS_QUAD
        use system,               only: symbol, nctype_tot, ncent_tot
        use gauss_ecp,          only: ecp_coef, ecp_exponent, necp_power, necp_term
        use gauss_ecp,          only: allocate_gauss_ecp
        use pseudo,             only: lpot
        use qua,                only: nquad, wq, xq0, yq0, zq0
        use general,            only: pooldir
        use contrl_file,        only: ounit
        use rotqua_mod,         only: gesqua

        use precision_kinds,    only: dp

        implicit none

        !   local use
        character(len=72), intent(in)   :: file_trexio
        character(len=80)               :: file_trexio_path


        ! trexio
        integer(8)                      :: trex_ecp_file
        integer                         :: rc = 1

        ! local variables
        integer                         :: ecp_num
        integer, allocatable            :: flat_ecp_ang_mom(:)
        integer, allocatable            :: idx_array(:)
        integer, allocatable            :: flat_ecp_nucleus_index(:)
        integer, allocatable            :: flat_ecp_max_ang_mom_plus_1(:)
        integer, allocatable            :: flat_ecp_power(:)
        integer, allocatable            :: flat_ecp_z_core(:)
        real(dp), allocatable           :: flat_ecp_coefficient(:)
        real(dp), allocatable           :: flat_ecp_exponent(:)

        integer, allocatable            :: unique_atom_index(:)
        character(len=2), allocatable   :: unique(:) ! unique symbols of atoms

        integer, dimension(:), allocatable :: atom_index(:), component_index_atom(:), components_per_atom(:)
        integer, dimension(:), allocatable :: nterms_per_component(:), term_index_component(:)
        integer                         :: count, lower_comp, upper_comp, counter_comp


        integer         :: i, ic, idx, l, tcount1, j, count_term

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
            trex_ecp_file = trexio_open(file_trexio_path, 'r', backend, rc)
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio file open', __FILE__, __LINE__)
            rc = trexio_read_ecp_num(trex_ecp_file, ecp_num)
            if (trexio_has_ecp(trex_ecp_file) == 0) trexio_has_group_ecp = .true.
            call trexio_error(rc, TREXIO_SUCCESS, 'trexio_read_ecp_num', __FILE__, __LINE__)
        endif
        call bcast(trexio_has_group_ecp)
        call bcast(ecp_num)
        allocate (flat_ecp_ang_mom(ecp_num))
        allocate (idx_array(ecp_num))
        allocate (flat_ecp_nucleus_index(ecp_num))
        allocate (flat_ecp_max_ang_mom_plus_1(ncent_tot))
        allocate (flat_ecp_power(ecp_num))
        allocate (flat_ecp_z_core(ncent_tot))
        allocate (flat_ecp_coefficient(ecp_num))
        allocate (flat_ecp_exponent(ecp_num))

        if (wid) then
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
        endif
        call bcast(flat_ecp_ang_mom)
        call bcast(flat_ecp_nucleus_index)
        call bcast(flat_ecp_max_ang_mom_plus_1)
        call bcast(flat_ecp_power)
        call bcast(flat_ecp_z_core)
        call bcast(flat_ecp_coefficient)
        call bcast(flat_ecp_exponent)

        ! Obtain the number of unique types of atoms stored in the hdf5 file.
        if (.not. allocated(unique)) allocate(unique(nctype_tot))
        if (.not. allocated(unique_atom_index)) allocate(unique_atom_index(nctype_tot))


        tcount1 = 1
        unique_atom_index(1) = 1
        unique(1) = symbol(1)
        do j= 2, ncent_tot
            if (any(unique == symbol(j) ))  then
                cycle
            endif
            tcount1 = tcount1 + 1
            unique_atom_index(tcount1) = j
            unique(tcount1) = symbol(j)
        enddo

        if (.not. allocated(lpot)) allocate (lpot(nctype_tot))
        call allocate_gauss_ecp()

        allocate(atom_index(ecp_num), source=0)
        allocate(components_per_atom(ncent_tot*10), source=0)
        allocate(component_index_atom(ncent_tot*10), source=0)

        call unique_elements(ecp_num, flat_ecp_nucleus_index, atom_index, count, components_per_atom, component_index_atom)

        ! count                     :: "Number of unique elements". Not used
        ! atom_index(1:count)       :: "Unique elements index (not used here)"
        ! components_per_atom(1:count) :: "frequency or count of shells per atom"
        ! component_index_atom(1:count) :: "index number of the shell for each atom"
        ! The components per atom can be obtained by accessing the component_index_atom
        ! for a given atom index by the slice of size frequency.


        do ic = 1, nctype_tot

            ! loop over all the primitives for the unique atom
            ! The lower and upper indices of primitive indices
            lower_comp = component_index_atom(unique_atom_index(ic))
            upper_comp = component_index_atom(unique_atom_index(ic)) + components_per_atom(unique_atom_index(ic)) - 1

            lpot(ic) = flat_ecp_max_ang_mom_plus_1(unique_atom_index(ic)) + 1

            write(ounit,'(a,i4,a,a)') 'ECP for atom type ', ic, ' Element = ', unique(ic)
            write(ounit,*) '-----------------------------------------------------------------------'
            write(ounit,*)
            write(ounit,'(a,i4,a,i4)') 'ECP for atom type ', ic, ' lpot = ', lpot(ic)

            if(lpot(ic).gt.MPS_L) call fatal_error('READPS_GAUSS: increase MPS_L')

            if (.not. allocated(nterms_per_component)) allocate(nterms_per_component(lpot(ic)))
            if (.not. allocated(term_index_component)) allocate(term_index_component(lpot(ic)))


            count_term = 0
            do i = lower_comp, upper_comp
                if(flat_ecp_ang_mom(i) .eq. lpot(ic)-1) then
                    idx_array(i) = lpot(ic)
                    else
                    count_term = count_term + 1
                    idx_array(i) = count_term
                endif
            enddo


            counter_comp = 0;
            do l = 1, lpot(ic)
                atom_index = 0
                call unique_elements(components_per_atom(unique_atom_index(ic)), flat_ecp_ang_mom(lower_comp:upper_comp), atom_index, count, nterms_per_component, term_index_component)

                write(ounit,*) "nterms_per_component", nterms_per_component

                write(ounit,*)
                write(ounit,'(a,2i6)') '    component, #terms ', l, nterms_per_component(l)

                do i=1,nterms_per_component(l)
                    idx=idx_array(lower_comp+counter_comp)
                    necp_term(idx,ic) = nterms_per_component(l)

                    ecp_coef(i,idx,ic) = flat_ecp_coefficient(lower_comp + counter_comp)
                    necp_power(i,idx,ic) = flat_ecp_power(lower_comp + counter_comp) + 2
                    ecp_exponent(i,idx,ic) = flat_ecp_exponent(lower_comp + counter_comp)
                    counter_comp = counter_comp + 1

                    write(ounit,'(a,f16.8,i2,f16.8)') '    coef, power, expo ', ecp_coef(i,idx,ic), &
                    necp_power(i,idx,ic), ecp_exponent(i,idx,ic)

                enddo
            enddo  ! loop on l upto lpot(ic)
            if (allocated(nterms_per_component)) deallocate(nterms_per_component)
            if (allocated(term_index_component)) deallocate(term_index_component)
            write(ounit,*) '-----------------------------------------------------------------------'
            write(ounit,*)
        enddo

        deallocate(atom_index)
        deallocate(components_per_atom)
        deallocate(component_index_atom)

        deallocate(flat_ecp_ang_mom)
        deallocate(flat_ecp_nucleus_index)
        deallocate(flat_ecp_max_ang_mom_plus_1)
        deallocate(flat_ecp_power)
        deallocate(flat_ecp_z_core)
        deallocate(flat_ecp_coefficient)
        deallocate(flat_ecp_exponent)

        if (.not. allocated(wq)) allocate (wq(MPS_QUAD))
        if (.not. allocated(xq0)) allocate (xq0(MPS_QUAD))
        if (.not. allocated(yq0)) allocate (yq0(MPS_QUAD))
        if (.not. allocated(zq0)) allocate (zq0(MPS_QUAD))

        call gesqua(nquad,xq0,yq0,zq0,wq)
    end subroutine read_trexio_ecp_file
#endif


end module
