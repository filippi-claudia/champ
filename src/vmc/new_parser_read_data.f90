
subroutine header_printing()
    ! Ravindra
    use mpi_f08
    use, intrinsic :: iso_fortran_env, only: iostat_end
    implicit none 
    
    integer                             :: status, i
    character(len=8)                    :: date
    character(len=10)                   :: time
    character(len=40)                   :: env_variable
    character(len=100)                  :: input_filename, output
    
    

    write(*,*) "____________________________________________________________________"
    write(*,*)                   
    write(*,*)
    write(*,*) ' .d8888b.   888    888         d8888  888b     d888  8888888b. '
    write(*,*) 'd88P  Y88b  888    888        d88888  8888b   d8888  888   Y88b'
    write(*,*) '888    888  888    888       d88P888  88888b.d88888  888    888'
    write(*,*) '888         8888888888      d88P 888  888Y88888P888  888   d88P'
    write(*,*) '888         888    888     d88P  888  888 Y888P 888  8888888P" '
    write(*,*) '888    888  888    888    d88P   888  888  Y8P  888  888       '
    write(*,*) 'Y88b  d88P  888    888   d8888888888  888   "   888  888       '
    write(*,*) ' "Y8888P"   888    888  d88P     888  888       888  888       '
    write(*,*)
    write(*,*) "____________________________________________________________________"
    write(*,*)
    write(*,*) ' Cornell Holland Ab-initio Materials Package'
    write(*,*)
    write(*,*)

    write(*,*) " information about the contributors goes here"
    write(*,*)
    write(*,*)
    write(*,*)
    write(*,*)

    write(*,*) " paper to cite for this code goes here"
    write(*,*)
    write(*,*)
    write(*,*)
    write(*,*)

    write(*,*) " license information goes here"

    write(*,*) "____________________________________________________________________"
    write(*,*)
    write(*,*)
    write(*,*)
    write(*,*)

    call date_and_time(date=date,time=time)
    write(*, '(12a)') " Calculation started on     :: ",   date(1:4), "-", date(5:6), "-", date(7:8), " at ",  time(1:2), ":", time(3:4), ":", time(5:6)                                                               
    call get_environment_variable ("PWD", output)
    write(*, '(2a)') " Current directory          :: ",   output
    call get_command_argument(number=0, value=output)
    write(*, '(2a)') " Executable                 :: ",   output
    call hostnm(output)
    write(*, '(2a)') " Hostname                   :: ",   output
    write(*,*)
    write(*,*)
    write(*,*)
    write(*,*)

    ! Get Command line arguments
    ! i = 0
    ! do
    !     call get_command_argument(number=i, value=output)
    !     if (len_trim(output) == 0) exit
    !     print*, "number ", i, " output", trim(output)
    !     i = i+1
    ! end do


    
    
    ! compiled date and time 
    ! source directory 
    ! branch   
    ! git commit id 

    ! input filename  
    ! nproc           
   
   
end subroutine header_printing


subroutine read_molecule_file(file_molecule)
    ! This subroutine reads the .xyz molecule file.
    ! Ravindra

    use atom, only: znuc, cent, pecent, iwctype, nctype, ncent, ncent_tot, nctype_tot, symbol, atomtyp    
    use ghostatom, only: nghostcent
    use inputflags, only: igeometry
    use periodic_table, only: atom_t, element

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
    character(len=100)               :: int_format     = '(A, T30, I8)'
    character(len=100)               :: float_format   = '(A, T60, f12.8)'    
    character(len=100)               :: string_format  = '(A, T60, A)'  
  
    !   External file reading
    write(6,*) '-----------------------------------------------------------------------'      
    write(6,string_format)  " Reading molecular coordinates from the file :: ",  trim(file_molecule)
    write(6,*) '-----------------------------------------------------------------------'      

    inquire(file=file_molecule, exist=exist)
    if (exist) then
        open (newunit=iunit,file=file_molecule, iostat=iostat, action='read' )
        if (iostat .ne. 0) stop "Problem in opening the molecule file"
    else
        error stop " molecule file "// trim(file_molecule) // " does not exist."
    endif

    read(iunit,*) ncent
    write(*,fmt=int_format) " Number of atoms ::  ", ncent
    write(*,*)

    if (.not. allocated(cent)) allocate(cent(3,ncent))
    if (.not. allocated(symbol)) allocate(symbol(ncent)) 
    if (.not. allocated(iwctype)) allocate(iwctype(ncent))              
    if (.not. allocated(unique)) allocate(unique(ncent))  
    
    read(iunit,'(A)')  comment
    write(*,*) "Comment from the molecule file :: ", trim(comment)
    write(*,*)

    do i = 1, ncent
        read(iunit,*) symbol(i), cent(1,i), cent(2,i), cent(3,i)
    enddo
    close(iunit)


    ! Count unique type of elements
    nctype = 1 
    unique(1) = symbol(1)
    do j= 2, ncent  
        if (any(unique == symbol(j) ))  cycle
        nctype = nctype + 1 
        unique(nctype) = symbol(j)
    enddo

    write(*,*) " Number of distinct types of elements (nctype) :: ", nctype 
    write(*,*)

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

    write(6,*) 'Atomic symbol, coordinates, and iwctype from the molecule coordinates file '
    write(*,*)
    do j= 1, ncent
        write(6,'(A4,3F10.6, i3)') symbol(j), (cent(i,j),i=1,3), iwctype(j)
    enddo

    write(*,*)
    write(*,*) " Values of znuc (number of valence electrons) "
    write(*,'(10F10.6)') (znuc(j), j = 1, nctype)
    write(*,*)
end subroutine read_molecule_file


subroutine read_determinants_file(file_determinants)
    ! This subroutine reads the single state determinant file.
    ! Ravindra

    use, intrinsic :: iso_fortran_env, only: iostat_eor   
    use dets,           only: cdet, ndet
    use dorb_m,         only: iworbd
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
    write(6,*) '------------------------------------------------------'      
    write(6,string_format)  " Reading determinants from the file :: ",  trim(file_determinants)
    write(6,*) '------------------------------------------------------'      

    inquire(file=file_determinants, exist=exist)
    if (exist) then
        open (newunit=iunit,file=file_determinants, iostat=iostat, action='read' )
        if (iostat .ne. 0) stop "Problem in opening the determinant file"
    else
        error stop " determinant file "// trim(file_determinants) // " does not exist."
    endif

    ndn  = nelec - nup        

    write(*,*)     
    write(*,int_format) " Number of total electrons ", nelec
    write(*,int_format) " Number of alpha electrons ", nup        
    write(*,int_format) " Number of beta  electrons ", ndn
    write(*,*) 


    ! to escape the comments before the "lcao nbasis norb" line
    do while (skip)
        read(iunit,*, iostat=iostat) temp1
        temp1 = trim(temp1)
        if (temp1 == "determinants") then
            backspace(iunit)
            skip = .false. 
        endif
    enddo

!   Read the first main line
    read(iunit, *, iostat=iostat)  temp2, ndet, nwftype
    if (iostat == 0) then 
        if (trim(temp2) == "determinants") write(*,int_format) " Number of determinants ", ndet 
    else
        error stop "Error in reading number of determinants / number of wavefunction types"
    endif


    if (.not. allocated(cdet)) allocate(cdet(ndet,1,nwftype))           

    read(iunit,*, iostat=iostat) (cdet(i,1,1), i=1,ndet)
    if (iostat /= 0) error stop "Error in determinant coefficients "

    write(*,*)         
    write(*,*) " Determinant coefficients "
    write(*,'(10(1x, f11.8, 1x))') (cdet(i,1,1), i=1,ndet)   
    
!       allocate the orbital mapping array        
    if (.not. allocated(iworbd)) allocate(iworbd(nelec, ndet))
    
    do i = 1, ndet
        read(iunit,*, iostat=iostat) (iworbd(j,i), j=1,nelec)
        if (iostat /= 0) error stop "Error in reading orbital -- determinants mapping "
    enddo
    
    write(*,*)     
    write(*,*) " Orbitals <--> Determinants mapping :: which orbitals enter in which dets"
    do i = 1, ndet
        write(*,'(<nelec>(i4, 1x))') (iworbd(j,i), j=1,nelec)
    enddo
    
    read(iunit,*) temp1
    if (temp1 == "end" ) write(*,*) " Single state determinant file read successfully "

    close(iunit)
end subroutine read_determinants_file


subroutine read_jastrow_file(file_jastrow)
    ! This subroutine reads jastrow parameters from a file.
    ! Ravindra

    use, intrinsic :: iso_fortran_env, only: iostat_eor !, iostat_eof   

    use force_mod,          only: MWF
    use jaspar,             only: nspin1, nspin2
    use elec,               only: ndn
    use jaspar3,            only: a, b, c, scalek
    use jaspar4,            only: a4, norda, nordb, nordc
    use jaspar6,            only: cutjas
    use bparm,              only: nocuspb, nspin2b
    use contr2,             only: ifock, ijas
    use contr2,             only: isc
    use inputflags,         only: ijastrow_parameter
    use wfsec,              only: nwftype
    use atom,               only: ncent, nctype
    use precision_kinds,    only: dp

    implicit none

    !   local use  
    character(len=72), intent(in)   :: file_jastrow
    character(len=40)               :: temp1, temp2, temp3, temp4, temp5   
    integer                         :: iunit, iostat, it, isp, iparm, iwft
    integer                         :: mparmja, mparmjb, mparmjc, nterms4
    logical                         :: exist
    real(dp)                        :: a21

    !   Formatting
    character(len=100)               :: int_format     = '(A, T60, I8)'
    character(len=100)               :: string_format  = '(A, T60, A)'  
  
    !   External file reading
    write(6,*) '---------------------------------------------------------------------------'      
    write(6,string_format)  " Reading jastrow parameters from the file :: ",  trim(file_jastrow)
    write(6,*) '---------------------------------------------------------------------------'      
    
    inquire(file=file_jastrow, exist=exist)
    if (exist) then
        open (newunit=iunit,file=file_jastrow, iostat=iostat, action='read' )
        if (iostat .ne. 0) error stop "Problem in opening the jastrow file"
    else
        error stop " Jastrow file "// trim(file_jastrow) // " does not exist."
    endif



    if (ijas .lt. 4 .or. ijas .gt. 6) error stop 'JASTROW: only ijas=4,5,6 implemented'
    if (ndn .eq. 1 .and. nspin2 .eq. 3) error stop 'JASTROW: 1 spin down and nspin2=3'

    if ((ijas .eq. 4 .or. ijas .eq. 5) .and. &
        (isc .ne. 2 .and. isc .ne. 4 .and. isc .ne. 6 .and. isc .ne. 7 .and. &
         isc .ne. 12 .and. isc .ne. 14 .and. isc .ne. 16 .and. isc .ne. 17)) &
         error stop 'JASTROW: if ijas=4 or 5, isc must be one of 2,4,6,7,12,14,16,17'

    if ((ijas .eq. 6) .and. (isc .ne. 6 .and. isc .ne. 7)) &
        error stop 'JASTROW: if ijas=6, isc must be 6 or 7'

    nspin2b = iabs(nspin2)
    nocuspb = 0
    if (nspin2 .lt. 0) then
        if (nspin2 .eq. -1) nocuspb = 1
        nspin2 = 1
    endif

    ! read the first word of the file
    read(iunit, *, iostat=iostat)  temp2, iwft
    if (iostat == 0) then 
        if (trim(temp2) == "jastrow_parameter") write(*,int_format) " Jastrow parameters being read : type of wavefunctions :: ", iwft
    else
        error stop "Error in reading jastrow parameters / number of wavefunction types"
    endif

    allocate (scalek(nwftype))

    if (ijas .ge. 4 .and. ijas .le. 6) then
        if (ifock .gt. 0) error stop 'JASTROW: fock not yet implemented for ijas=4,5,6'
        read (iunit, *) norda, nordb, nordc
        write (*, '(3(A,i4))') " norda = ", norda, "; nordb = ", nordb, "; nordc = ", nordc 

        if (isc .ge. 2) read (iunit, *) scalek(iwft), a21
        write (*, '(2(A,f12.6))') " scalek = ", scalek(iwft), "; a21 = ", a21

        mparmja = 2 + max(0, norda - 1)
        mparmjb = 2 + max(0, nordb - 1)
        mparmjc = nterms4(nordc)

        allocate (a4(mparmja, nctype, nwftype))

        write (*, '(A)') "Jastrow parameters :: "
        do it = 1, nctype
            read (iunit, *) (a4(iparm, it, iwft), iparm=1, mparmja)
            write (*, '(<mparmja>(2X,f12.8))') (a4(iparm, it, iwft), iparm=1, mparmja)
        enddo

        allocate (b(mparmjb, 2, nwftype))

        do isp = nspin1, nspin2b
            read (iunit, *) (b(iparm, isp, iwft), iparm=1, mparmjb)
            write (*, '(<mparmjb>(2X,f12.8))') (b(iparm, isp, iwft), iparm=1, mparmjb)
        enddo

        allocate (c(mparmjc, nctype, nwftype))

        do it = 1, nctype
            read (iunit, *) (c(iparm, it, iwft), iparm=1, mparmjc)
            write (*, '(<mparmjc>(2X,f12.8))') (c(iparm, it, iwft), iparm=1, mparmjc)
        enddo

    endif
    !Read cutoff for Jastrow4, 5, 6
    if (isc .eq. 6 .or. isc .eq. 7) then 
        read (iunit, *) cutjas
        write(iunit, '(A,2X,f12.8)') " cutjas = ", cutjas
    endif

    ijastrow_parameter = ijastrow_parameter + 1
    
    close(iunit)

end subroutine read_jastrow_file


subroutine read_orbitals_file(file_orbitals)
    
    use coefs, only: coef, nbasis, norb
    use inputflags, only: ilcao
    use orbval, only: nadorb
    use pcm_fdc, only: fs

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
    character(len=100)               :: int_format     = '(A, T60, I8)'
    character(len=100)               :: string_format  = '(A, T60, A)'  
    character(len=100)               :: float_format   = '(A, T60, f12.8)'    

    !   External file reading
    write(6,*) '---------------------------------------------------------------------------'      
    write(6,string_format)  " Reading LCAO orbitals from the file :: ",  trim(file_orbitals)
    write(6,*) '---------------------------------------------------------------------------'      
    
    inquire(file=file_orbitals, exist=exist)
    if (exist) then
        open (newunit=iunit,file=file_orbitals, iostat=iostat, action='read' )
        if (iostat .ne. 0) error stop "Problem in opening the LCAO orbitals file"
    else
        error stop " Jastrow file "// trim(file_orbitals) // " does not exist."
    endif

    ! to escape the comments before the "lcao nbasis norb" line
    do while (skip)
        read(iunit,*, iostat=iostat) temp1
        temp1 = trim(temp1)
        if (temp1 == "lcao") then
            backspace(iunit)
            skip = .false. 
        endif
    enddo

    ! read the first line 
    read(iunit, *, iostat=iostat)  temp1, nbasis, norb, iwft

    if (iostat == 0) then 
        if (trim(temp2) == "lcao") then
            write(*,int_format) " Number of basis functions ", nbasis
            write(*,int_format) " Number of lcao orbitals ", norb            
            write(*,int_format) " Type of wave functions ", iwft
        endif
    else
        write(*, *) " Check ", temp1, nbasis, norb, iwft
        error stop "Error in reading number of lcao orbitals / basis / number of wavefunction types"
    endif

    ! Fix the maximum size of all array relative
    ! to MOs with the maximum number of MOs
    ! this may not be needed later
    !MORB = norb

    if (iwft .gt. nwftype) error stop 'LCAO: wave function type > nwftype'

    if (.not. allocated(coef)) allocate (coef(nbasis, norb, nwftype))           

    do iorb = 1, norb
        read (iunit, *, iostat=iostat) (coef(ibasis, iorb, iwft), ibasis=1, nbasis)
    enddo
    if (iostat /= 0) error stop "Error in reading lcao orbitals "

    write(*,*)         
    write(*,*) " LCAO orbitals "

    temp3 = '(T8, T14, i3, T28, i3, T42, i3, T56, i3, T70, i3, T84, i3, T98, i3, T112, i3, T126, i3, T140, i3)'
    ! print orbs in blocks of 10
    counter = 0
    do k = 10, nbasis, 10
!        write(*,*) " Orbitals  ", k-9 , "  to ", k       
        write (*, fmt=temp3 )  (i, i = k-9, k)
        do iorb = 1, norb
            write (*, '(A,i5,A, 10(1x, f12.8, 1x))') "[", iorb, "] ", (coef(ibasis, iorb, iwft), ibasis=k-9, k)
        enddo
        counter = counter + 10
    enddo


    ! Remaining block
    write (*, fmt=temp3 )  (i, i = counter, nbasis)        
    do k = counter, nbasis
!        write(*,*) " Orbitals  ", counter , "  to ", nbasis
        do iorb = 1, norb
            write (*, '(A,i5,A, 10(1x, f12.8, 1x))') "[", iorb, "] ", (coef(ibasis, iorb, iwft), ibasis=counter, nbasis)
        enddo
    enddo

    close(iunit)
    write(*,*) "----------------------------------------------------------"        

end subroutine read_orbitals_file

subroutine read_csf_file(file_determinants)
    ! This subroutine reads the csf coefficients from the determinant file.
    ! Ravindra

    use, intrinsic :: iso_fortran_env!, only: is_iostat_end
    use vmc_mod, only: MDET
    use csfs, only: ccsf, ncsf, nstates
    use mstates_mod, only: MSTATES
    use inputflags, only: icsfs
    use wfsec, only: nwftype
    use dets, only: ndet

    implicit none

    !   local use  
    character(len=72), intent(in)   :: file_determinants
    character(len=40)               :: temp1, temp2, temp3, temp4, temp5   
    integer                         :: iostat, i, j, iunit
    logical                         :: exist

    !   Formatting
    character(len=100)              :: int_format     = '(A, T40, I8)'
    character(len=100)              :: string_format  = '(A, T40, A)'  
    
    !   External file reading
    write(6,*) '------------------------------------------------------'      
    write(6,string_format)  " Reading csf from the file :: ",  trim(file_determinants)
    write(6,*) '------------------------------------------------------'      

    inquire(file=file_determinants, exist=exist)
    if (exist) then
        open (newunit=iunit,file=file_determinants, iostat=iostat, action='read' )
        if (iostat .ne. 0) stop "Problem in opening the determinant file for reading csfs"
    else
        error stop " determinant file "// trim(file_determinants) // " does not exist."
    endif        

    do 
        read(iunit,*, iostat=iostat) temp1
        temp1 = trim(temp1)
        if (is_iostat_end(iostat)) exit


        if (temp1 == "csf") then
            backspace(iunit)   ! go a line back
            read(iunit, *, iostat=iostat)  temp2, ncsf, nstates
            write(*,*) " Number of csf and nstates ", ncsf, nstates
            if (iostat == 0) then 
                if (.not. allocated(ccsf)) allocate(ccsf(ncsf, nstates, nwftype))    
                do i = 1, nstates
                    read(iunit,*, iostat=iostat) (ccsf(j,i,1), j=1,ncsf)
                enddo
                if (iostat /= 0) error stop "Error in reading csf coefficients "
            else
                error stop "Error in reading number of csfs / number of states"
            endif
        endif 
        
        
    enddo

    write(*,*)         
    write(*,*) " CSF coefficients "
    
    do i = 1, nstates
        write(*,*) " State :: ", i , " out of ", nstates
        write(*,'(10(1x, f11.8, 1x))') (ccsf(j,i,1), j=1,ncsf)   
        write(*,*) 
    enddo   
    close(iunit)

end subroutine read_csf_file

subroutine read_csfmap_file(file_determinants)
    ! This subroutine reads the csf coefficients from the determinant file.
    ! Ravindra

    use, intrinsic :: iso_fortran_env
    use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates
    use dets, only: cdet, ndet
    use wfsec, only: nwftype
    use mstates_mod, only: MDETCSFX
    use precision_kinds,    only: dp

    implicit none

    !   local use  
    character(len=72), intent(in)   :: file_determinants
    character(len=40)               :: temp1, temp2, temp3, temp4, temp5   
    integer                         :: iostat, i, j, iunit
    integer                         :: nptr, nterm, id, nmap
    real(dp)                        :: c
    logical                         :: exist

    !   Formatting
    character(len=100)              :: int_format     = '(A, T40, I8)'
    character(len=100)              :: string_format  = '(A, T40, A)'  
    
    !   External file reading
    write(6,*) '------------------------------------------------------'      
    write(6,string_format)  " Reading csfmap from the file :: ",  trim(file_determinants)
    write(6,*) '------------------------------------------------------'      

    inquire(file=file_determinants, exist=exist)
    if (exist) then
        open (newunit=iunit,file=file_determinants, iostat=iostat, action='read' )
        if (iostat .ne. 0) stop "Problem in opening the determinant file for reading csfmap"
    else
        error stop " determinant file "// trim(file_determinants) // " does not exist."
    endif        

    
    
    do 
        read(iunit,*, iostat=iostat) temp1
        temp1 = trim(temp1)
        if (is_iostat_end(iostat)) exit


        if (temp1 == "csfmap") then
            backspace(iunit)   ! go a line back
            read(iunit, *, iostat=iostat)  temp2, ncsf, ndet, nmap
            write(*,*) " Number of csf, number of determinants, and number of mappings ", ncsf, ndet, nmap
            if (iostat == 0) then 
                if (.not. allocated(cxdet)) allocate (cxdet(ndet*MDETCSFX))
                if (.not. allocated(iadet)) allocate (iadet(ndet))
                if (.not. allocated(ibdet)) allocate (ibdet(ndet))
                if (.not. allocated(icxdet)) allocate (icxdet(ndet*MDETCSFX))                
                
                nptr = 1
                do i = 1, ncsf
                    read (iunit, *) nterm
                    iadet(i) = nptr
                    ibdet(i) = nptr + nterm - 1
                    do j = 1, nterm
                        read (iunit, *) id, c
                        icxdet(nptr) = id
                        cxdet(nptr) = c
                        nptr = nptr + 1
                        if (nptr .gt. ndet*MDETCSFX) error stop 'Error in CSFMAP:: problem with nmap'
                    enddo
                enddo

                if (nmap .ne. nptr - 1) error stop 'Error in CSFMAP:: not enough nmaps / file is corrupt'
                nmap = nptr
            
                if (.not. allocated(cdet)) allocate (cdet(ndet, nstates, nwftype))
        
                ! write (6, '(''Warning: det coef overwritten with csf'')')
                ! do k = 1, nstates
                !     do j = 1, ndet
                !         cdet(j, k, 1) = 0
                !     enddo
                !     do icsf = 1, ncsf
                !         do j = iadet(icsf), ibdet(icsf)
                !             jx = icxdet(j)
                !             cdet(jx, k, 1) = cdet(jx, k, 1) + ccsf(icsf, k, 1)*cxdet(j)
                !         enddo
                !     enddo
                ! enddo
                

            else
                error stop "Error in reading number of csfs, number of determinants, or number of mappings"
            endif
        endif         
    enddo

    write(*,*)         
    close(iunit)
    
    ! do i = 1, nstates
    !     write(*,*) " State :: ", i , " out of ", nstates
    !     write(*,'(10(1x, f11.8, 1x))') (ccsf(j,i,1), j=1,ncsf)   
    !     write(*,*) 
    ! enddo   

end subroutine read_csfmap_file




subroutine read_exponents_file(file_exponents)
    ! Read basis function exponents (only if no numerical basis)
    ! Ravindra
    
    use coefs, only: nbasis
    use basis, only: zex
    use inputflags, only: iexponents
    use wfsec, only: nwftype

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
    write(6,*) '------------------------------------------------------'      
    write(6,string_format)  " Reading exponents from the file :: ",  trim(file_exponents)
    write(6,*) '------------------------------------------------------'      

    inquire(file=file_exponents, exist=exist)
    if (exist) then
        open (newunit=iunit,file=file_exponents, iostat=iostat, action='read' )
        if (iostat .ne. 0) stop "Problem in opening the exponents file for reading csfs"
    else
        error stop " exponents file "// trim(file_exponents) // " does not exist."
    endif        


    write (6, *) 'nbasis', nbasis
    write (6, *) 'nwftype', nwftype

    if (.not. allocated(zex)) allocate (zex(nbasis, nwftype))    
    
    do iwft = 1, nwftype
        read(iunit,*, iostat=iostat)  (zex(i, iwft), i=1, nbasis)

        if (iostat /= 0) error stop "Error in reading exponents from the exponent file "

        write(*,*)         
        write(*,*) " Basis set exponents "
        
        write(*,'(10(1x, f11.8, 1x))') (zex(i, iwft), i=1, nbasis)
        write(*,*) 
    enddo
    close(iunit)

end subroutine read_exponents_file


