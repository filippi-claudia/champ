module parser_read_data
      use error,   only: fatal_error
      implicit none
contains
subroutine header_printing()
    !> This subroutine prints the header in each output file. It contains some
    !! useful information about the compilers, version of the code, input and output file names.
    !! @author Ravindra Shinde (r.l.shinde@utwente.nl)

      use contrl_file, only: errunit,file_error,file_input,file_output
      use contrl_file, only: ounit
      use mpi_f08
      use mpiconf, only: idtask,nproc
#if defined(TREXIO_FOUND)
      use trexio
#endif
    use, intrinsic :: iso_fortran_env, only: iostat_end

    implicit none

    integer                             :: status, i
    character(len=8)                    :: date
    character(len=10)                   :: time
    character(len=40)                   :: env_variable
    character(len=100)                  :: input_filename, output
    character(len=MPI_MAX_PROCESSOR_NAME):: hostname


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

#if defined(GIT_HEAD_BRANCH)
    write(ounit,'(2a)')  " Git branch                 :: ", GIT_HEAD_BRANCH
#endif

#if defined(GIT_REVISION_HASH)
    write(ounit,'(2a)')  " Git commit hash            :: ", GIT_REVISION_HASH
#endif

#if defined(CMAKE_Fortran_COMPILER)
    write(ounit,'(2a)')  " Compiler                   :: ", CMAKE_Fortran_COMPILER
#endif

#if defined(CMAKE_Fortran_COMPILER_VERSION)
    write(ounit,'(2a)')  " Compiler version           :: ", CMAKE_Fortran_COMPILER_VERSION
#endif

#if defined(TARGET_ARCHITECTURE)
    write(ounit,'(2a)')  " Vectorization Instructions :: ", TARGET_ARCHITECTURE
#endif

#if defined(HDF5_VERSION)
    write(ounit,'(2a)')  " HDF5 library version       :: ", HDF5_VERSION
#endif

    call MPI_Get_processor_name(hostname,i)
    write(ounit, '(2a)') " Hostname                   :: ",   hostname
    call get_environment_variable ("PWD", output)
    write(ounit, '(2a)') " Current directory          :: ",   output
    call get_environment_variable ("USER", output)
    write(ounit, '(2a)') " Username                   :: ",   output
    write(ounit, '(2a)') " Input file                 :: ",   file_input
    write(ounit, '(2a)') " Output file                :: ",   file_output
    write(ounit, '(2a)') " Error file                 :: ",   file_error
    write(ounit, '(4a)') " Code compiled on           :: ",__DATE__, " at ", __TIME__
    write(ounit, '(a,i0)') " Number of processors       :: ", nproc
#if defined(TREXIO_FOUND)
    if (TREXIO_SUCCESS == 0) write(ounit,*) "TREXIO library version     :: ", TREXIO_PACKAGE_VERSION
#endif
    write(ounit,*)



end subroutine header_printing


subroutine read_molecule_file(file_molecule)
    !> This subroutine reads the .xyz molecule file. It then computes the
    !! number of types of atoms, nuclear charges (from the symbol), and
    !! number of valence electrons if pseudopotential is provided.
    !! @author Ravindra Shinde (r.l.shinde@utwente.nl)
    !! @date
      use contrl_file, only: errunit,ounit
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use inputflags, only: igeometry
      use m_string_operations, only: wordcount
      use mpiconf, only: wid
      use multiple_geo, only: pecent
      use periodic_table, only: atom_t,element
      use precision_kinds, only: dp
      use system,  only: atomtyp,cent,iwctype,ncent,ncent_tot,nctype
      use system,  only: nctype_tot,newghostype,nghostcent,symbol,znuc

    implicit none

    !   local use
    character(len=72), intent(in)   :: file_molecule
    character(len=40)               :: temp1, temp2, temp3, temp4
    character(len=80)               :: comment, file_molecule_path, line
    integer                         :: iostat, i, j, k, iunit, count
    logical                         :: exist
    type(atom_t)                    :: atoms
    character(len=2), allocatable   :: unique(:)
    double precision, allocatable   :: nval(:)

    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: float_format   = '(A, T60, f12.8)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading

    if((file_molecule(1:6) == '$pool/') .or. (file_molecule(1:6) == '$POOL/')) then
            file_molecule_path = pooldir // file_molecule(7:)
    else
            file_molecule_path = file_molecule
    endif

    write(ounit,*) '-----------------------------------------------------------------------'
    write(ounit,string_format)  " Reading molecular coordinates from the file :: ",  file_molecule_path
    write(ounit,*) '-----------------------------------------------------------------------'

    if (wid) then
            inquire(file=file_molecule_path, exist=exist)
            if (exist) then
                    open (newunit=iunit,file=file_molecule_path, iostat=iostat, action='read' )
                    if (iostat .ne. 0) stop "Problem in opening the molecule file"
            else
                    call fatal_error (" molecule file "// pooldir // trim(file_molecule) // " does not exist.")
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
    if (.not. allocated(nval)) allocate(nval(ncent))
    unique = ''
    symbol = ''

    if (wid) read(iunit,'(A)')  comment
    call bcast(comment)

    write(ounit,*) "Comment from the molecule file :: ", trim(comment)
    write(ounit,*)

    if (wid) then
            read(iunit,'(A)')  line
            backspace(iunit)
    endif
    call bcast(line)
    count = wordcount(line)

    if (count == 4) then
            ! Read the symbol and coords only
            if (wid) then
                    do i = 1, ncent
                    read(iunit,*) symbol(i), cent(1,i), cent(2,i), cent(3,i)
                    enddo
            endif
            call bcast(symbol)
            call bcast(cent)

            if (wid) close(iunit)
    else
            ! Read the symbol, coords, and nvalence (or znuc)
            if (wid) then
                    do i = 1, ncent
                    read(iunit,*) symbol(i), cent(1,i), cent(2,i), cent(3,i), nval(i)
                    enddo
            endif
            call bcast(symbol)
            call bcast(cent)
            call bcast(nval)
            if (wid) close(iunit)
    endif

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
    if (symbol(j) == unique(k))  then
            iwctype(j) = k
            if (count .gt. 4) znuc(k) = nval(j)
    endif
    enddo
    enddo

    ! Get the correspondence rule
    do k = 1, nctype
    atomtyp(k) = unique(k)
    enddo

    if (allocated(unique)) deallocate(unique)

    if (count == 4) then
            ! Get the znuc for each unique atom
            do j = 1, nctype
            atoms = element(atomtyp(j))
            znuc(j) = atoms%nvalence
            enddo
    endif

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
end subroutine read_molecule_file





subroutine read_determinants_file(file_determinants)
    !> This subroutine reads the single state determinant file.
    !! @author Ravindra Shinde

      use contrl_file, only: errunit,ounit
      use csfs,    only: nstates
      use custom_broadcast, only: bcast
      use dorb_m,  only: iworbd
      use general, only: pooldir
      use inputflags, only: ideterminants
      use mpiconf, only: wid
      use mstates_mod, only: MSTATES
      use multiple_geo, only: nwftype
      use optwf_control, only: method
      use precision_kinds, only: dp
      use slater,  only: cdet,ndet,norb
      use system,  only: ndn,nelec,nup
    use, intrinsic :: iso_fortran_env, only: iostat_eor

    implicit none

    !   local use
    character(len=72), intent(in)   :: file_determinants
    character(len=80)               :: temp1, temp2, temp3
    integer                         :: iostat, i, j, iunit, counter, istate, itmp
    logical                         :: exist, skip = .true., found = .false.

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

    nwftype = 0
!   nstates = 1  ! read from the input file
    if (wid) then
        do while (.not. found)
            read(iunit,*, iostat=iostat) temp1
            if (is_iostat_end(iostat)) exit
            temp1 = trim(temp1)
            if (temp1 == "determinants") then
                nwftype = nwftype + 1
            endif
        enddo
        rewind(iunit)
    endif

    if (wid) then ! so you can put whatever you want before this line
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
        read(iunit, *, iostat=iostat)  temp2, ndet, itmp
        if (iostat == 0) then
            if (trim(temp2) == "determinants") write(ounit,int_format) " Number of determinants ", ndet
        else
            call fatal_error ("Error in reading number of determinants / number of wavefunction types")
        endif
    endif
    call bcast(ndet)
    call bcast(nwftype)


    if( (method(1:3) == 'lin')) then
        if (.not. allocated(cdet)) allocate(cdet(ndet,MSTATES,3))
    else
        if (.not. allocated(cdet)) allocate(cdet(ndet,MSTATES,nwftype))
    endif
    !       allocate the orbital mapping array
    if (.not. allocated(iworbd)) allocate(iworbd(nelec, ndet))

    write(ounit,int_format) " # of sets of dets read from the file ", nwftype

    if (wid) then
        do itmp = 1, nwftype
            read(iunit,*, iostat=iostat) (cdet(i,1,itmp), i=1,ndet)
            if (iostat /= 0) call fatal_error( "Error in determinant coefficients ")
            write(ounit,*)
            write(ounit,*) " Determinant coefficients "
            write(ounit,'(10(1x, f11.8, 1x))') (cdet(i,1,itmp), i=1,ndet)

            do i = 1, ndet
                read(iunit,*, iostat=iostat) (iworbd(j,i), j=1,nelec)
                if (iostat /= 0) call fatal_error("Error in reading orbital -- determinants mapping ")
            enddo
            read(iunit,*, iostat=iostat) temp1
        enddo
    endif
    call bcast(cdet)
    call bcast(iworbd)
    call bcast(temp1)
    ! This part replaces a call to verify_orbitals
    !if(any(iworbd .gt. norb))  call fatal_error('INPUT: iworbd > norb')


    write(ounit,*)
    write(ounit,*) " Orbitals <--> Determinants mapping :: which orbitals enter in which dets"
    write(temp3, '(a,i0,a)') '(', nelec, '(i4, 1x))'
    do i = 1, ndet
        !write(ounit,'(<nelec>(i4, 1x))') (iworbd(j,i), j=1,nelec)  !Intel Version
        write(ounit,temp3) (iworbd(j,i), j=1,nelec)                 !GNU version
    enddo

    if (temp1 == "end" ) then
        write(ounit,*) " Determinant file read successfully "
        ideterminants = ideterminants + 1
    endif

    if (wid) close(iunit)

end subroutine read_determinants_file

subroutine read_multideterminants_file(file_multideterminants)
    !> This subroutine reads the multideterminants file. The first line appears as 'multideterminants' ndet_local
    !! CI coefficients and occupation of determinants in wf
    !! @author Ravindra Shinde
      use contrl_file, only: errunit,ounit
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use inputflags, only: imultideterminants
      use mpiconf, only: wid
      use multidet, only: irepcol_det,ireporb_det,numrep_det
      use slater,  only: iwundet,ndet
      use system,  only: nelec

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

      use bparm,   only: nocuspb,nspin2b
      use contrl_file, only: errunit,ounit
      use contrl_per, only: iperiodic
      use control, only: mode
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use inputflags, only: ijastrow_parameter
      use jastrow, only: norda,nordb,nordc
      use jaspar6, only: asymp_r,c1_jas6,c1_jas6i,c2_jas6,cutjas,cutjasi
      use jastrow, only: a4,asymp_jasa,asymp_jasb,b,c,ijas,isc,neqsx
      use jastrow, only: nordj,nordj1,nspin1,nspin2,scalek
      use jastrow4_mod, only: nterms4
      use mpiconf, only: wid
      use multiple_geo, only: MWF,nwftype
      use optwf_control, only: method, ioptwf
      use vmc_mod, only: nwftypeorb, nwftypejas, nstoj, jtos, stoj, nstojmax, extraj, nstoj_tot
      use precision_kinds, only: dp
      use system,  only: ncent,nctype,ndn
      use, intrinsic :: iso_fortran_env, only: iostat_eor !, iostat_eof
    implicit none

    !   local use
    character(len=72), intent(in)   :: file_jastrow
    character(len=40)               :: temp1, temp2, temp3, temp4, temp5
    integer                         :: iunit, iostat, it, isp, iparm, iwft
    integer                         :: mparmja, mparmjb, mparmjc, int1, int2
    logical                         :: exist, found = .false.
    real(dp)                        :: cutjas_tmp
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

    !check how many sets of jastrows are listed. make only for sr_n, multistate?
    extraj = 0
    if ((method .eq. 'sr_n')) then
        nwftypejas = 0
        extraj = 1
        if (wid) then
            do while (.not. found)
                read(iunit,*, iostat=iostat) temp2
                if (is_iostat_end(iostat)) exit
                temp1 = trim(temp2)
                if (temp2 == "jastrows_to_states") then
                    nwftypejas = nwftypejas + 1
                endif
            enddo
            rewind(iunit)
            if ((nwftypejas .eq. 0)) then
                nwftypejas = 1
                extraj = 0
            endif
            if(mode(1:3) == 'vmc') write(ounit,'(A,i4,A)') " Found ", nwftypejas, " types of jastrows in the input file. "
            if (extraj .eq. 0) then
               if(mode(1:3) == 'vmc') write(ounit,int_format) " For a multistate sr_n optimization , a 'jastrows_to_states' line is required before each jastrow's parameters are specified "
               if(mode(1:3) == 'vmc') write(ounit,'(A,i4,A)') " Assuming 1-state sr_n optimization for the moment. "
            endif
            !if (nwftypejas .ne. nwftypeorb) then
                !write(ounit,int_format) " Using sr_n with multiple states requires a 'jastrow' line before each state's parameters are specified "
                !call fatal_error ( " Number of states specified in jastrow file do not match lcao file ")
            !endif
        endif
        call bcast(nwftypejas)
        call bcast(extraj)
        allocate (nstoj(nwftypejas))
        i = 0
        if (extraj.eq.1) then
          if (wid) then
            do while (.not. found)
              read(iunit,*, iostat=iostat) temp1, int1
              if (is_iostat_end(iostat)) exit
              temp1 = trim(temp1)
              if (temp1 == "jastrows_to_states") then
                i = i + 1
                nstoj(i) = int1
              endif
            enddo
            rewind(iunit)
          endif
        else
          nstoj(1) = 1 ! will redefine later, no, we don't need it.
          if (ioptwf.ge.1) then
            if(mode(1:3) == 'vmc') write(ounit,int_format) " This is an sr_n optimization. Assuming 1 jastrow type and 1-state for the moment."
          else
            if(mode(1:3) == 'vmc') write(ounit,int_format) " No optimization to be performed. Will decide jastrow mapping once we know how many states. "
          endif
        endif
        call bcast(nstoj)
    else ! assuming 1-jastrow listing with sr_n does not need 'jastrows line', annoying because could imply use for all states.
        allocate (nstoj(1))
        nstoj(1)=1
        call bcast(nstoj)
    endif

    nstoj_tot=0
    do i=1,nwftypejas
      nstoj_tot=nstoj_tot+nstoj(i)
    enddo

    if(mode(1:3) == 'vmc') write(ounit,'(A,i4,A,i4,A)') " Found ", nwftypejas, " jastrow types to be assigned to ", nstoj_tot, " states. If only sampling, ignore."

    ! read the first word of the file
    if (wid) then
        read (iunit, *, iostat=iostat)  temp2, iwft
        if (iostat == 0) then
            if (trim(temp2) == "jastrow_parameter") &
            write(ounit,int_format) " Jastrow parameters being read : type of wavefunctions :: ", iwft
            if((method .eq. 'sr_n')) then 
                !if (nwftypejas .ne. nwftypeorb) then
                    !call fatal_error ( " Number of states specified in jastrow file do not match lcao file ")
                !else
              nwftype = nwftypejas
                !endif
            endif
        else
            call fatal_error ("Error in reading jastrow parameters / number of wavefunction types")
        endif
    endif
    call bcast(nwftype)
    call bcast(iwft)

    if((method(1:3) == 'lin')) then
        allocate (scalek(3))
    else
        allocate (scalek(nwftype)) 
    endif
          

    if (ijas .ge. 4 .and. ijas .le. 6) then
        if (wid) read (iunit, *) norda, nordb, nordc
        call bcast(norda)
        call bcast(nordb)
        call bcast(nordc)

        nordj = max(norda, nordb, nordc)
        nordj1 = nordj + 1
        neqsx = 6*nordj
        
        write(ounit, '(3(A,i4))') " norda = ", norda, "; nordb = ", nordb, "; nordc = ", nordc

        if (isc .ge. 2) then
            if (wid) read (iunit, *) scalek(iwft) ! we set iwft = 1 for 'sr_n' so only the (1) will have the saclek value
        endif

        if (method.eq.'sr_n') then 
          do i=1,nwftype
            scalek(i)=scalek(1)
          enddo
        endif

        call bcast(scalek)
        write(ounit, '(A,f12.6)') " scalek = ", scalek(iwft)

        mparmja = 2 + max(0, norda - 1)
        mparmjb = 2 + max(0, nordb - 1)
        mparmjc = nterms4(nordc)

        if( (method(1:3) == 'lin')) then
            allocate (a4(mparmja, nctype, 3))
            allocate (b(mparmjb, 2, 3))
            allocate (c(mparmjc, nctype, 3))
        else
            allocate (a4(mparmja, nctype, nwftype))
            allocate (b(mparmjb, 2, nwftype))
            allocate (c(mparmjc, nctype, nwftype))
        endif
       
        nstojmax=maxval(nstoj)
        if(mode(1:3) == 'vmc') write(ounit,*) "nstoj,nstojmax", nstoj(1), nstojmax
        if(mode(1:3) == 'vmc') write(ounit,*) "Ignore if only sampling."
        allocate(jtos(nwftypejas,nstojmax))

        if (wid) then 
        do iwft = 1, nwftype
            write(ounit, '(A)') "Jastrow parameters :: "
            if( (method .eq. 'sr_n')) write(ounit, '(A,i0)') "Jastrow type :: ", iwft
            write(ounit, '(A)') "mparmja : "
            write(temp3, '(a,i0,a)') '(', mparmja, '(2X, f12.8))'
            if (method .eq. 'sr_n') then
                if (extraj .eq. 1) then
                    read (iunit, *, iostat=iostat) temp1, int1, (jtos(iwft,i),i=1,nstoj(iwft))
                    write(ounit,*) 'jastrows_to_states'
                    write(ounit,*) '# of states assigned: ', int1 
                    write(ounit,*) 'states assigned: ', (jtos(iwft,i),i=1,nstoj(iwft))
                    int2=int1+1
                    do i=int2,nstojmax ! setting unused spaces to zero.
                       jtos(iwft,i) = 0
                    enddo
                    if (trim(temp1) .ne. "jastrows_to_states") then
                        write(ounit, *) " Expected 'jastrows_to_states' on this line ", temp1
                        call fatal_error ("Error in reading jastrow parameters")
                    endif
                    if (iostat .ne. 0) call fatal_error ("Error in reading jastrows_to_states")
                else
                  jtos(1,1)=1
                endif
            else
                jtos(1,1)=1
            endif
            do it = 1, nctype
                read (iunit, *) (a4(iparm, it, iwft), iparm=1, mparmja)
                !write(ounit, '(<mparmja>(2X,f12.8))') (a4(iparm, it, iwft), iparm=1, mparmja)  !Intel Version
                if (mparmja .ne. 0) write(ounit, temp3) (a4(iparm, it, iwft), iparm=1, mparmja)                     !GNU version
            enddo
            write(ounit, '(A)') "mparmjb : "
            write(temp3, '(a,i0,a)') '(', mparmjb, '(2X, f12.8))'
            do isp = nspin1, nspin2b
                read (iunit, *) (b(iparm, isp, iwft), iparm=1, mparmjb)
                !write(ounit, '(<mparmjb>(2X,f12.8))') (b(iparm, isp, iwft), iparm=1, mparmjb)  !Intel Version
                if (mparmjb .ne. 0) write(ounit, temp3) (b(iparm, isp, iwft), iparm=1, mparmjb)                     !GNU version
            enddo
            write(ounit, '(A)') "mparmjc : "
            write(temp3, '(a,i0,a)') '(', mparmjc, '(2X, f12.8))'
            do it = 1, nctype
                read (iunit, *) (c(iparm, it, iwft), iparm=1, mparmjc)
                !write(ounit, '(<mparmjc>(2X,f12.8))') (c(iparm, it, iwft), iparm=1, mparmjc)   !Intel Version
                if (mparmjc .ne. 0) write(ounit, temp3) (c(iparm, it, iwft), iparm=1, mparmjc)                      !GNU version
            enddo
        enddo
        endif
        call bcast(a4)
        call bcast(b)
        call bcast(c)
        !call bcast(extraj)
        if (extraj.gt.0) then
           call bcast(jtos)
           call bcast(nstoj)
        endif
    endif

    !Read cutoff for Jastrow4, 5, 6
    if (isc .eq. 6 .or. isc .eq. 7) then
        if (wid) then
          read (iunit, *) cutjas
          write(iunit, '(A,2X,f12.8)') " cutjas = ", cutjas
        endif
    endif
    call bcast(cutjas)

    ijastrow_parameter = ijastrow_parameter + 1
    call bcast(ijastrow_parameter)

    if (wid) close(iunit)

    ! for now do this:
    if (method .eq. 'sr_n') nwftype=1

end subroutine read_jastrow_file


subroutine read_orbitals_file(file_orbitals)
    ! Ravindra

      use coefs,   only: nbasis
      use contrl_file, only: errunit,ounit
      use control, only: mode
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use inputflags, only: ilcao
      use mpiconf, only: wid
      use multiple_geo, only: nwftype
      use optwf_control, only: method, ioptwf
      use orbval,  only: nadorb
      use pcm_fdc, only: fs
      use precision_kinds, only: dp
      use slater,  only: coef,norb
      use vmc_mod, only: norb_tot, nwftypeorb, nstoo, nstoomax, otos, extrao, nstoo_tot      
      use write_orb_loc_mod, only: write_orb_loc

    ! was not in master but is needed

    implicit none

!   local use
    character(len=72), intent(in)   :: file_orbitals
    character(len=40)               :: temp1, temp2
    character(len=120)              :: temp3
    integer                         :: iunit, iostat, iwft
    integer                         :: iorb, ibasis, i, k, counter, int1, int2
    logical                         :: exist
    logical                         :: skip = .true., found = .false.

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

    !check how many sets of orbitals are listed
    extrao = 0
    if (method .eq. 'sr_n') then
        nwftypeorb = 0
        extrao = 1
        if (wid) then
            do while (.not. found)
                read(iunit,*, iostat=iostat) temp1
                if (is_iostat_end(iostat)) exit
                temp1 = trim(temp1)
                if (temp1 == "orbitals_to_states") then
                    nwftypeorb = nwftypeorb + 1
                endif
            enddo
            rewind(iunit)
            if (nwftypeorb .eq. 0) then
                nwftypeorb = 1
                extrao = 0
            endif
            if (mode(1:3) == 'vmc') write(ounit,'(A,i4,A)') " Found ", nwftypeorb, " sets of orbitals in the input file. "
            if (extrao .eq. 0) then
               if (mode(1:3) == 'vmc') write(ounit,int_format) " For a multistate sr_n optimization , a 'orbitals_to_states' line is required before each orbital set's parameters are specified "
               if (mode(1:3) == 'vmc') write(ounit,'(A,i4,A)') " Assuming 1-state sr_n optimization for the moment. "
            endif
            nwftype = nwftypeorb !nwftype is used by force stuff outside this subroutine. don't really know how nwftype is used
        endif
        call bcast(nwftypeorb)
        call bcast(nwftype)
        call bcast(extrao)
        allocate (nstoo(nwftypeorb))
        i = 0
        if (extrao.eq.1) then
          if (wid) then
            do while (.not. found)
               read(iunit,*, iostat=iostat) temp1, int1
               if (is_iostat_end(iostat)) exit
               temp1 = trim(temp1)
               if (temp1 == "orbitals_to_states") then
                 i = i + 1
                 nstoo(i) = int1
               endif
            enddo
            rewind(iunit)
          endif
        else
          nstoo(1)=1
          if(ioptwf.eq.1) then
            if (mode(1:3) == 'vmc') write(ounit,int_format) " This is an sr_n optimization. Assuming 1 orbital set and 1-state for the moment. "
          else
            if (mode(1:3) == 'vmc') write(ounit,int_format) " No optimization to be performed. Will decide orbital mapping once we know how many states. "
          endif
        endif
        call bcast(nstoo)
    else
        allocate (nstoo(1))
        nstoo(1)=1
        call bcast(nstoo)
    endif

    nstoo_tot=0
    do i=1,nwftypeorb
      nstoo_tot=nstoo_tot+nstoo(i)
    enddo

    ! to escape the comments before the "lcao nbasis norb" line
    if (wid) then
        if (mode(1:3) == 'vmc') write(ounit,'(A,i4,A,i4,A)') " Found ", nwftypeorb, " orbital sets to be assigned to ", nstoo_tot, " states. Ignore if only sampling"
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
    if (wid) read(iunit, *, iostat=iostat)  temp1, norb_tot, nbasis, iwft
    call bcast(nbasis)
    call bcast(norb_tot)
    call bcast(iwft)

    if (wid) then
        if (iostat == 0) then
            if (trim(temp1) == "lcao") then
                write(ounit,int_format) " Number of basis functions ", nbasis
                write(ounit,int_format) " Number of lcao orbitals ", norb_tot
                write(ounit,int_format) " Type of wave functions ", iwft
                norb = norb_tot     ! norb will get updated later. norb_tot is fixed
            endif
        else
            write(ounit, *) " Check ", temp1, norb, nbasis, iwft
            call fatal_error ("Error in reading number of lcao orbitals / basis / number of wavefunction types")
        endif
    endif
    call bcast(norb)

    if (iwft .gt. nwftype) call fatal_error('LCAO: wave function type > nwftype') ! I don't think nwftype is calculated anywhere, besides default

    if( (method(1:3) == 'lin')) then
        if (.not. allocated(coef)) allocate (coef(nbasis, norb_tot, 3))
    else
        if (.not. allocated(coef)) allocate (coef(nbasis, norb_tot, nwftype))
    endif

    nstoomax=maxval(nstoo)
    if (mode(1:3) == 'vmc') write(ounit,*) "nstoo,nstoomax", nstoo(1), nstoomax
    if (mode(1:3) == 'vmc') write(ounit,*) "Ignore, if only sampling"
    allocate(otos(nwftypeorb,nstoomax))

    if (wid) then ! Moved if to outside, all the important quantities are bcast later.
    do iwft = 1, nwftype ! make sure nwftype = 1 at this point for all other methods.
      if (method .eq. 'sr_n') then
        if (extrao .eq. 1) then
          read (iunit, *, iostat=iostat) temp1, int1, (otos(iwft,i),i=1,nstoo(iwft))
          int2=int1+1
          do i=int2,nstoomax ! setting unused to 0
            otos(iwft,i) = 0
          enddo
          if (trim(temp1) .ne. "orbitals_to_states") then
            write(ounit, *) " Expected 'orbitals_to_states' on this line ", temp1
            call fatal_error ("Error in reading number of lcao orbitals / basis / number of wavefunction types")
            if (iostat .ne. 0) call fatal_error ("Error in reading orbitals_to_states.")
          endif
        else
          otos(1,1)=1
        endif
      else
        otos(1,1)=1
      endif
      do iorb = 1, norb
        !if (wid) then
          read (iunit, *, iostat=iostat) (coef(ibasis, iorb, iwft), ibasis=1, nbasis)
          if (iostat /= 0) call fatal_error( "Error in reading lcao orbitals ")
        !endif
      enddo
    enddo
    endif
    ! printing of the lcao orbitals coefficients will be done by write_orb_loc subroutine.
    write(ounit,*) "Orbital coefficients are written to the output.log file"
    call bcast(coef)
    if (extrao.gt.0) then
       call bcast(otos)
       call bcast(nstoo)
    endif

    ilcao = ilcao + 1

    call bcast(ilcao)
    if (wid) close(iunit)

    write(ounit,*) "----------------------------------------------------------"

    ! for now do this:
    if (method .eq. 'sr_n') nwftype=1

end subroutine read_orbitals_file


subroutine read_csf_file(file_determinants)
    ! This subroutine reads the csf coefficients from the determinant file.
    ! Ravindra

      use ci000,   only: nciprim,nciterm
      use contrl_file, only: errunit,ounit
      use csfs,    only: ccsf,ncsf,nstates
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use inputflags, only: icsfs
      use mpiconf, only: wid
      use mstates_mod, only: MSTATES
      use multiple_geo, only: nwftype
      use optwf_control, only: ioptci,method
      use precision_kinds, only: dp
      use slater,  only: cdet,ndet
      use vmc_mod, only: nwftypeorb, nstoj_tot, nstoo_tot      
      use, intrinsic :: iso_fortran_env!, only: is_iostat_end
!   Not sure about the following two lines
    implicit none

    !   local use
    character(len=72), intent(in)   :: file_determinants
    character(len=40)               :: temp1, temp2, temp3, temp4, temp5
    integer                         :: iostat, i, j, iunit
    integer                         :: nstates_local
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
        
        if (nstates .gt. 1) then
            call fatal_error (" Provide CSF information for multiple states. Exiting!")
        endif

        ncsf = ndet
        if((method.eq.'sr_n'.and.ioptci.eq.1)) then
            if (nstates.ne.nstoj_tot.or.nstates.ne.nstoo_tot) then
                call fatal_error (" determinant file: number of states found do not match those in jastrow and lcao. ")
            endif
        endif
        if (ioptci .ne. 0) nciterm = nciprim
        printed = .true.
        ! if there is no mention of "csf" in the file:: allocate and assign
        if( (method(1:3) == 'lin')) then
            if (.not. allocated(ccsf)) allocate(ccsf(ndet, nstates, 3))
        else
            if (.not. allocated(ccsf)) allocate(ccsf(ndet, nstates, nwftype))
        endif


        do j = 1, ndet
            ccsf(j,1,nwftype) = cdet(j,1,nwftype)
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
        if (wid) read(iunit, *, iostat=iostat)  temp2, ncsf, nstates_local
        call bcast(ncsf)
        call bcast(nstates_local)

        if (nstates_local .ne. nstates) then
           call fatal_error ("Number of states read from the determinant file does not match with the input file. Exiting!")
        endif

        if((method.eq.'sr_n'.and.ioptci.eq.1)) then
            if (nstates.ne.nstoj_tot.or.nstates.ne.nstoo_tot) then
                write(ounit,*) 'nstates, nstoj_tot, nstoo_tot', nstates, nstoj_tot, nstoo_tot
                call fatal_error (" determinant file: number of states found for csf's do not match those in jastrow and lcao. ")
            endif
        endif

        if( (method(1:3) == 'lin')) then
            if (.not. allocated(ccsf)) allocate(ccsf(ncsf, nstates, 3)) 
        else
            if (.not. allocated(ccsf)) allocate(ccsf(ncsf, nstates, nwftype))
        endif

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

        write(ounit,'(8x,10(1x, a9, i3, 1x))') (" State: ", i, i =1, nstates)
        do j = 1, ncsf
            write(ounit,'(a,i5,a,10(1x, f12.8, 1x))') "[", j, "] ", (ccsf(j,i,1), i=1,nstates)
        enddo
    endif

    icsfs = 1
    call bcast(icsfs)

    if (wid) close(iunit)

end subroutine read_csf_file

subroutine read_csfmap_file(file_determinants)
    ! This subroutine reads the csf coefficients from the determinant file.
    ! Ravindra

      use contrl_file, only: errunit,ounit
      use csfs,    only: ccsf,cxdet,iadet,ibdet,icxdet,ncsf,nstates
      use custom_broadcast, only: bcast
      use dets,    only: nmap
      use general, only: pooldir
      use mpiconf, only: wid
      use multiple_geo, only: nwftype
      use precision_kinds, only: dp
      use slater,  only: cdet,ndet

    use, intrinsic :: iso_fortran_env

    implicit none

    !   local use
    character(len=72), intent(in)   :: file_determinants
    character(len=40)               :: temp1, temp2, temp3, temp4, temp5
    integer                         :: iostat, i, j, k, iunit
    integer                         :: icsf, jx
    integer                         :: nptr, nterm, id
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
        nmap = ndet
        if (.not. allocated(cxdet)) allocate (cxdet(nmap))
        if (.not. allocated(iadet)) allocate (iadet(ndet))
        if (.not. allocated(ibdet)) allocate (ibdet(ndet))
        if (.not. allocated(icxdet)) allocate (icxdet(nmap))

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
            else
              call fatal_error('CSFMAP: error reading in file')
            endif
        endif
        call bcast(ncsf_check)
        call bcast(ndet_check)
        call bcast(nmap_check)

        if (.not. allocated(cxdet)) allocate (cxdet(nmap_check))
        if (.not. allocated(iadet)) allocate (iadet(ndet))
        if (.not. allocated(ibdet)) allocate (ibdet(ndet))
        if (.not. allocated(icxdet)) allocate (icxdet(nmap_check))

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
                if (id > ndet) call fatal_error ('id in CSFMAP larger than ndet_check')
                icxdet(nptr) = id
                cxdet(nptr) = c
                nptr = nptr + 1
                if (nptr - 1 .gt. nmap_check) call fatal_error ('Error in CSFMAP:: problem with nmap')
            enddo
        enddo
        if (nmap_check .ne. nptr - 1) call fatal_error ('Error in CSFMAP:: not enough nmaps / file is corrupt')
        nmap = nptr - 1
!        if (allocated(cdet)) deallocate (cdet)
!        if (.not. allocated(cdet)) allocate (cdet(ndet, nstates, nwftype)) 
        write(ounit, *) "ndet, nstates, nwftype", ndet, nstates, nwftype
        write(ounit, *) "ncsf", ncsf
        write(ounit, '(''Warning: det coef overwritten with csf'')')

        do k = 1, nstates
            do j = 1, ndet
                cdet(j, k, 1) = 0.0d0
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
                write(ounit,'(1(5x, i0, t15, f12.8, 1x))') icxdet(j), cxdet(j)
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

      use basis,   only: zex
      use coefs,   only: nbasis
      use contrl_file, only: errunit,ounit
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use inputflags, only: iexponents
      use mpiconf, only: wid
      use multiple_geo, only: nwftype
      use optwf_control, only: method
      use precision_kinds, only: dp

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

    if( (method(1:3) == 'lin')) then
        if (.not. allocated(zex)) allocate (zex(nbasis, 3))
    else
        if (.not. allocated(zex)) allocate (zex(nbasis, nwftype))
    endif

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

      use bparm,   only: nspin2b
      use contrl_file, only: errunit,ounit
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use jastrow, only: norda,nordb,nordc
      use jaspointer, only: npoint,npointa
      use jastrow, only: ijas,is,isc,nspin1
      use mpiconf, only: wid
      use numbas,  only: numr
      use optwf_nparmj, only: nparma,nparmb,nparmc,nparmf
      use optwf_parms, only: nparmj
      use optwf_wjas, only: iwjasa,iwjasb,iwjasc,iwjasf
      use system,  only: nctype,nctype_tot
      use vmc_mod, only: nctyp3x

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
    write(ounit,string_format)  " Reading jastrow derivative parameters from the file :: ", trim(file_jastrow_der)
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
    nctyp3x = max(3, nctype_tot)

    if (.not. allocated(nparma)) allocate (nparma(nctyp3x))
    if (.not. allocated(nparmb)) allocate (nparmb(3))
    if (.not. allocated(nparmc)) allocate (nparmc(nctype))
    if (.not. allocated(nparmf)) allocate (nparmf(nctype))

    if (.not. allocated(iwjasa)) allocate (iwjasa(83, nctyp3x))
    if (.not. allocated(iwjasb)) allocate (iwjasb(83, 3))
    if (.not. allocated(iwjasc)) allocate (iwjasc(83, nctype))
    if (.not. allocated(iwjasf)) allocate (iwjasf(15, nctype))

    if (.not. allocated(npoint)) allocate (npoint(nctyp3x))
    if (.not. allocated(npointa)) allocate (npointa(3*nctyp3x))

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

        if (wid) then
        do it = 1, nctype
            read (iunit, *) (iwjasa(iparm, it), iparm=1, nparma(it))
            write(ounit, '(A,30i4)') " iwjasa = ", (iwjasa(iparm, it), iparm=1, nparma(it))
        enddo
        do isp = nspin1, nspin2b
            read (iunit, *) (iwjasb(iparm, isp), iparm=1, nparmb(isp))
            write(ounit, '(A,30i4)') " iwjasb = ", (iwjasb(iparm, isp), iparm=1, nparmb(isp))
        enddo
        do it = 1, nctype
            read (iunit, *) (iwjasc(iparm, it), iparm=1, nparmc(it))
            write(ounit, '(A,30i4)') " iwjasc = ", (iwjasc(iparm, it), iparm=1, nparmc(it))
        enddo
            ! end of reading the jasderiv file block
        endif
        call bcast(iwjasa)
        call bcast(iwjasb)
        call bcast(iwjasc)

    endif ! if the file containing jasderiv exists

    if (wid) close(iunit)

end subroutine read_jasderiv_file


subroutine read_forces_file(file_forces)
    !
    ! Ravindra
      use contrl_file, only: errunit,ounit
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use inputflags, only: iforces
      use mpiconf, only: wid
      use multiple_geo, only: delc,iwftype,nforce
      use precision_kinds, only: dp
      use system,  only: ncent,symbol


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
    write(ounit,string_format)  " Reading force displacements from the file :: ", trim(file_forces)
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

      use contrl_file, only: errunit,ounit
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use mpiconf, only: idtask,wid
      use optorb,  only: irrep
      use precision_kinds, only: dp
      use slater,  only: norb
      use vmc_mod, only: norb_tot


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
    write(ounit,string_format)  " Reading orbital symmetries from the file :: ", trim(file_symmetry)
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
        if (iostat /= 0) call fatal_error( "Error in reading symmetry file :: expecting 'sym_labels', nsym, norb_tot")
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
    call bcast(temp2)
    write(ounit, *) "Irreducible representation correspondence for all norb orbitals"
    write(ounit, *) temp2

    ! safe allocate
    if (.not. allocated(irrep)) allocate (irrep(norb_tot))

    ! read data
    if (wid) then
        read (iunit, *, iostat=iostat) (irrep(io), io=1, norb)
        if (iostat/=0) call fatal_error("Error in reading symmetry file :: expecting irrep correspondence for all norb orbitals")
    endif
    call bcast(irrep)

    write(ounit, '(10(1x, i3))') (irrep(io), io=1, norb)

    if (wid) close(iunit)
end subroutine read_symmetry_file


subroutine read_optorb_mixvirt_file(file_optorb_mixvirt)
    !
    ! Ravindra
      use contrl_file, only: errunit,ounit
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use inputflags, only: ioptorb_mixvirt
      use mpiconf, only: wid
      use optorb_mix, only: iwmix_virt,norbopt,norbvirt
      use precision_kinds, only: dp
      use slater,  only: norb


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
    write(ounit,string_format)  " Reading optorb_mixvirt from the file :: ", trim(file_optorb_mixvirt)
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

      use contrl_file, only: errunit,ounit
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use mpiconf, only: wid
      use optorb,  only: orb_energy
      use precision_kinds, only: dp
      use slater,  only: norb
      use vmc_mod, only: norb_tot


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
    write(ounit,string_format)  " Reading orbital eigenvalues from the file :: ", trim(file_eigenvalues)
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
    if (.not. allocated(orb_energy)) allocate (orb_energy(norb_tot))

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


      use basis, only: ns, np, nd, nf, ng
      use coefs,   only: nbasis
      use contrl_file, only: errunit,ounit
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use fdf, only: fdf_defined
      use inputflags, only: ibasis_num
      use mpiconf, only: wid
      use numbas,  only: iwrwf,numr,nrbas
      use numbas1, only: iwlbas,nbastyp
      use numbas_mod, only: MRWF
      use precision_kinds, only: dp
      use system,  only: nctype,newghostype
      use write_orb_loc_mod, only: write_orb_loc


    implicit none

    !   local use
    character(len=72), intent(in)   :: file_basis_num_info
    character(len=40)               :: temp1, temp2
    character(len=72)               :: file_path
    integer                         :: iunit, iostat
    integer                         :: i,j, jj, ib, nctot
    logical                         :: exist, skip = .true.


    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I0)'
    character(len=100)               :: string_format  = '(A, T60, A)'

    !   External file reading
    write(ounit,*) '---------------------------------------------------------------------------'
    write(ounit,'(a)')  " Reading Basis function types and pointers to radial parts tables from the file :: ", &
                        pooldir // trim(file_basis_num_info)
    write(ounit,*) '---------------------------------------------------------------------------'

    if((file_basis_num_info(1:6) == '$pool/') .or. (file_basis_num_info(1:6) == '$POOL/')) then
        file_path = pooldir // file_basis_num_info(7:)
    else
        file_path = file_basis_num_info
    endif

    if (wid) then
        inquire(file=file_path, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_path, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error( "Problem in opening the basis num info file")
        else
            call fatal_error (" Basis num info file "// pooldir // trim(file_basis_num_info) // " does not exist.")
        endif

        ! To skip the comments in the file
        do while (skip)
            read(iunit,*, iostat=iostat) temp1
            temp1 = trim(temp1)
            if (temp1 == "qmc_bf_info") then
                backspace(iunit)
                skip = .false.
            endif
        enddo


        read (iunit, *, iostat=iostat) temp1, numr
        if (iostat /= 0) call fatal_error( "Error in reading basis num info file :: expecting 'qmc_bf_info / basis', numr")

        if (.not. ((trim(temp1) == "qmc_bf_info")  .or. (trim(temp1) == "basis"))) then
            call fatal_error( "Error in reading basis num info file :: expecting 'qmc_bf_info / basis'")
        endif
    endif
    call bcast(numr)

    nctot = nctype + newghostype    ! DEBUG:: this statement might go. ghosttypes built-in

    if (fdf_defined("trexio")) then 
    allocate (nbastyp(nctot))
    allocate (ns(nctot))
    allocate (np(nctot))
    allocate (nd(nctot))
    allocate (nf(nctot))
    allocate (ng(nctot))
    else
    if (.not. allocated(nbastyp)) allocate (nbastyp(nctot))
    if (.not. allocated(ns)) allocate (ns(nctot))
    if (.not. allocated(np)) allocate (np(nctot))
    if (.not. allocated(nd)) allocate (nd(nctot))
    if (.not. allocated(nf)) allocate (nf(nctot))
    if (.not. allocated(ng)) allocate (ng(nctot))
    endif

    if (.not. allocated(iwlbas)) allocate (iwlbas(nbasis, nctot))
    if (.not. allocated(iwrwf))  allocate (iwrwf(nbasis, nctot))
    if (.not. allocated(nrbas)) allocate (nrbas(nctype + newghostype), source=0)

    if (wid) then
        do i = 1, nctype + newghostype
            read (iunit, *, iostat=iostat) nbastyp(i),  ns(i), np(i), nd(i), nf(i), ng(i)

            if (iostat /= 0) call fatal_error( "Error in reading basis num info file")
            write (ounit, '(a,i0)')  "Number of AOs for this center :: ", nbastyp(i)
            write (ounit, '(5(a,i0))') "ns= ", ns(i), " np= ", np(i), " nd= ", nd(i), " nf= ", nf(i), " ng= ", ng(i)

            if (numr .gt. 0) then
                nrbas(i) = ns(i) + np(i) + nd(i) + nf(i) + ng(i)

                if (nbastyp(i) .gt. MRWF) call fatal_error('BASIS: nbastyp > MRWF')

                read (iunit, *) (iwlbas(ib, i), ib=1, nbastyp(i))

                if (iostat /= 0) call fatal_error( "Error in reading basis num info file :: iwlbas")
                write(ounit, '(100i3)') (iwlbas(ib, i), ib=1, nbastyp(i))

                read (iunit, *, iostat=iostat) (iwrwf(ib, i), ib=1, nbastyp(i))
                if (iostat /= 0) call fatal_error( "Error in reading basis num info file :: iwrwf")
                write(ounit, '(100i3)') (iwrwf(ib, i), ib=1, nbastyp(i))
                write(ounit, *)
            else
                call fatal_error('BASIS: numerical basis functions not supported')
            endif
        enddo

    endif

    call bcast(numr)
    call bcast(iwlbas)
    call bcast(iwrwf)
    call bcast(nbastyp)
    call bcast(ns)
    call bcast(np)
    call bcast(nd)
    call bcast(nf)
    call bcast(ng)

    ibasis_num = 1
    call bcast(ibasis_num)
    if (wid) close(iunit)

    write(ounit,*) "Orbital coefficients are written to the output.log file"
    call write_orb_loc()

end subroutine read_basis_num_info_file


subroutine read_dmatrix_file(file_dmatrix)
    ! Ravindra (no=ndetorb, ns=nweight)
    !INPUT dmatrix i i a=<input>
    !KEYDOC Read diagonal density matrix information.
      use contrl_file, only: errunit,ounit
      use csfs,    only: nstates
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use mpiconf, only: wid
      use mstates_mod, only: MSTATES
      use optorb,  only: dmat_diag
      use precision_kinds, only: dp
      use sa_weights, only: iweight,nweight,weights
      use slater,  only: norb
      use vmc_mod, only: norb_tot


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

    ipr = 0

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


    allocate (dmat(norb_tot))
    allocate (iwdmat(nstates))

    if (ndetorb .gt. norb) call fatal_error( 'READ_DMATRIX: wrong number of orbitals')

    allocate (weights(nstates))
    allocate (iweight(nstates))


    if (wid) read (iunit, *) (iwdmat(i), i=1, nweight)
    call bcast(iwdmat)

    do iw = 1, nweight
        if (iwdmat(iw) .ne. iweight(iw)) call fatal_error('READ_DMATRIX: iwdmat')
    enddo

    allocate (dmat_diag(norb_tot))
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
      use contrl_file, only: errunit,ounit
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use mpiconf, only: wid
      use pcm_parms, only: nesph,re,re2,xe,ye,ze
      use precision_kinds, only: dp


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
    write(ounit,string_format)  " Reading cavity spheres from the file :: ", trim(file_cavity_spheres)
    write(ounit,*) '---------------------------------------------------------------------------'

    if (wid) then
        inquire(file=file_cavity_spheres, exist=exist)
        if (exist) then
            open (newunit=iunit,file=file_cavity_spheres, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error( "Problem in opening the cavity spheres file")
        else
            call fatal_error (" cavity spheres file "//trim(file_cavity_spheres) // " does not exist.")
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
      use contrl_file, only: errunit,ounit
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use grdntsmv, only: igrdaidx,igrdcidx,igrdmv
      use grdntspar, only: delgrdxyz,igrdtype,ngradnts
      use inputflags, only: igradients
      use mpiconf, only: wid
      use multiple_geo, only: MFORCE,delc,iwftype,nforce
      use precision_kinds, only: dp
      use system,  only: ncent


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
    write(ounit,string_format)  " Reading gradients cartesian from the file :: ", trim(file_gradients_cartesian)
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
      use contrl_file, only: errunit,ounit
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use grdntsmv, only: igrdaidx,igrdcidx,igrdmv
      use grdntspar, only: delgrdba,delgrdbl,delgrdda,igrdtype,ngradnts
      use inputflags, only: igradients
      use misc_grdnts, only: grdzmat_displ
      use mpiconf, only: wid
      use multiple_geo, only: MFORCE,delc,iwftype,nforce
      use precision_kinds, only: dp
      use system,  only: ncent
      use zmatrix, only: izmatrix

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
    write(ounit,string_format)  " Reading gradients zmatrix from the file :: ", trim(file_gradients_zmatrix)
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
      use contrl_file, only: errunit,ounit
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use grdntsmv, only: igrdmv
      use inputflags, only: imodify_zmat
      use mpiconf, only: wid
      use system,  only: ncent


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
    write(ounit,string_format)  " Reading modify zmatrix from the file :: ", trim(file_modify_zmatrix)
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
      use contrl_file, only: errunit,ounit
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use grdnthes, only: hessian_zmat
      use inputflags, only: ihessian_zmat
      use mpiconf, only: wid
      use precision_kinds, only: dp
      use system,  only: ncent

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
    write(ounit,string_format)  " Reading hessian zmatrix from the file :: ", trim(file_hessian_zmatrix)
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
      use contrl_file, only: errunit,ounit
      use custom_broadcast, only: bcast
      use general, only: pooldir
      use inputflags, only: izmatrix_check
      use m_zmat_tools, only: cart2zmat,zmat2cart_rc
      use mpiconf, only: wid
      use precision_kinds, only: dp
      use system,  only: cent,ncent
      use zmatrix, only: czcart,czcart_ref,czint,izcmat,izmatrix

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
    write(ounit,string_format)  " Reading zmatrix connection matrix from the file :: ", trim(file_zmatrix_connection)
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
      use contrl_file, only: errunit,ounit
      use custom_broadcast, only: bcast
      use efield,  only: iscreen,ncharges
      use efield_blk, only: ascreen,bscreen,qcharge,xcharge,ycharge
      use efield_blk, only: zcharge
      use efield_mod, only: MCHARGES
      use general, only: pooldir
      use inputflags, only: icharge_efield
      use mpiconf, only: wid
      use precision_kinds, only: dp

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
    write(ounit,string_format)  " Reading efield from the file :: ", trim(file_efield)
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
end module
