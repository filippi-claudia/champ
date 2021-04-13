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
    character(len=40)               :: temp1, temp2, temp3, temp4, temp5   
    integer                         :: iostat, i, j, iunit
    logical                         :: exist

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

    ! read the first word of the file
    read(iunit,*) temp1
    temp1 = trim(temp1)

    if (temp1(1:1) == "#") then
        write(*,*) "Comment from the file ", temp1
        read(iunit, *, iostat=iostat)  temp2, ndet, nwftype
        if (iostat == 0) then 
            if (trim(temp2) == "determinants") write(*,int_format) " Number of determinants ", ndet 
        else
            error stop "Error in reading number of determinants / number of wavefunction types"
        endif
    else
        backspace(iunit)   ! go a line back
        read(iunit, *, iostat=iostat)  temp2, ndet, nwftype
        if (iostat == 0) then 
            if (trim(temp2) == "determinants") write(*,int_format) " Number of determinants ", ndet 
        else
            error stop "Error in reading number of determinants / number of wavefunction types"
        endif
    endif 

    !BUG :: the number 1 should be replaced by nstates. undefined at this point
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
        print *, " unit of file opened ", iunit
        if (iostat .ne. 0) stop "Problem in opening the jastrow file"
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
    read(iunit, *, iostat=iostat)  temp2, nwftype
    if (iostat == 0) then 
        if (trim(temp2) == "jastrow_parameter") write(*,int_format) " Jastrow parameters being read : types ", nwftype
    else
        error stop "Error in reading jastrow parameters / number of wavefunction types"
    endif

    allocate (scalek(nwftype))

    if (ijas .ge. 4 .and. ijas .le. 6) then
        if (ifock .gt. 0) error stop 'JASTROW: fock not yet implemented for ijas=4,5,6'
        read (iunit, *) norda, nordb, nordc
        write (*, *) norda, nordb, nordc
        if (isc .ge. 2) read (iunit, *) scalek(iwft), a21
        write (*, *) scalek(iwft), a21
        mparmja = 2 + max(0, norda - 1)
        mparmjb = 2 + max(0, nordb - 1)
        mparmjc = nterms4(nordc)

        allocate (a4(mparmja, nctype, nwftype))

        do it = 1, nctype
            read (iunit, *) (a4(iparm, it, iwft), iparm=1, mparmja)
            write (*, *) (a4(iparm, it, iwft), iparm=1, mparmja)
        enddo

        allocate (b(mparmjb, 2, nwftype))

        do isp = nspin1, nspin2b
            read (iunit, *) (b(iparm, isp, iwft), iparm=1, mparmjb)
            write (*, *) (b(iparm, isp, iwft), iparm=1, mparmjb)
        enddo

        allocate (c(mparmjc, nctype, nwftype))

        do it = 1, nctype
            read (iunit, *) (c(iparm, it, iwft), iparm=1, mparmjc)
            write (*, *) (c(iparm, it, iwft), iparm=1, mparmjc)
        enddo

    endif
    !Read cutoff for Jastrow4, 5, 6
    if (isc .eq. 6 .or. isc .eq. 7) read (iunit, *) cutjas

    ijastrow_parameter = ijastrow_parameter + 1
    call p2chkend(iunit, 'jastrow_parameter')
    
    close(iunit)

end subroutine read_jastrow_file

