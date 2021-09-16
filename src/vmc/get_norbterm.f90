subroutine get_norbterm
  !> This subroutine computes the number of orbital parameters needed for
  !> the allocations of the arrays.
  !> @author Ravindra Shinde
  !> @email  r.l.shinde@utwente.nl
  !> @date   14-09-2021

    use optorb_mod, only: mxreduced
    use vmc_mod, only: norb_tot
    use const, only: nelec
    use dets, only: ndet
    use elec, only: ndn, nup
    use multidet, only: kref
    use optorb_mix, only: norbopt, norbvirt , iwmix_virt
    use coefs, only: norb, next_max
    use dorb_m, only: iworbd
    use optorb, only: irrep
    use optorb_cblock, only: norbterm
    use orb_mat_022, only: ideriv
    use orb_mat_033, only: ideriv_iab, ideriv_ref, irepcol_ref
    use method_opt, only: method
    use optorb_cblock, only: nreduced
    use orbval, only: nadorb, ndetorb, orb
    use optwf_contrl, only: ncore, no_active
    use contrl_file,    only: ounit, errunit

    implicit none

    integer :: i, iab, icount_orbdef, ie, iesave
    integer :: io, iocc, iprt, iterm
    integer :: j, jo, k, n0
    integer :: n1, noporb
    integer :: local_norb, local_norb_tot, local_nadorb, local_ndetorb, local_ncore, local_noporb, local_next_max
    integer :: local_norbopt, local_norbvirt, local_norbterm, local_nreduced
    integer, dimension(2, ndet) :: iodet
    integer, dimension(2, ndet) :: iopos
    integer, dimension(2, norb_tot) :: iflag
    integer, dimension(2) :: ne
    integer, dimension(2) :: m
    integer, dimension(norb_tot, norb_tot) :: local_iwmix_virt


    data icount_orbdef /1/

    save icount_orbdef


    local_norb        = norb
    local_norb_tot    = norb_tot
    local_nadorb      = nadorb
    local_ndetorb     = ndetorb
    local_ncore       = ncore
    local_ndetorb     = ndetorb
    local_next_max    = next_max
    local_norbopt     = norbopt
    local_norbvirt    = norbvirt
    local_norbterm    = norbterm
    local_nreduced    = nreduced

    if (.not. allocated(iwmix_virt)) allocate (iwmix_virt(norb_tot, norb_tot))
    if (.not. allocated(irrep)) allocate (irrep(norb_tot))

    ! check the following condition against known cases
    if (method .eq. "linear") irrep = 0   !initialize irrep array to 0

    iprt=3

    ndn=nelec-nup

    ne(1)=nup
    ne(2)=nelec
    local_ndetorb=0

    do i=1,ndet
     do j=1,nelec
      if(iworbd(j,i).gt.norb) then
!       write(ounit,1) i,j,iworbd(j,i),norb
       call fatal_error('VERIFY: orbital index out of range')
      endif
      if(iworbd(j,i).gt.local_ndetorb) then
        local_ndetorb=iworbd(j,i)
      endif
     enddo
    enddo
1   format('Det ',i4,' column ',i4,' orb index ',i4,' norb ',i4)

!   Number of external orbitals for orbital optimization
    local_next_max  = local_norb  - local_ndetorb
    if(nadorb.gt.next_max) local_nadorb = local_next_max
    ! write(ounit, '(a, t40, i0)' ) 'norb', local_norb
    ! write(ounit, '(a, t40, i0)') 'nadorb', local_nadorb
    ! write(ounit, '(a, t40, i0)') 'ndet_orb', local_ndetorb
    ! write(ounit, '(a, t40, i0)') 'next_max', local_next_max

    ! if(iprt.gt.0) then
    !  write(ounit,'(''Determinantal orbitals in orbital optimization: '',i0)') local_ndetorb
    !  write(ounit,'(''External orbitals in orbital optimization:      '',i0)') local_nadorb
    !  write(ounit,'(''Total orbitals in orbital optimization:         '',i0)') local_nadorb + local_ndetorb - local_ncore
    ! endif

!   norb becomes number of occupied orbitals in wf; total needed are norb+adorb
    local_norb=local_ndetorb

!   Omit doubly occupied in all input determinants
    do 5 i=1,local_ndetorb
      iflag(1,i)=0
      do 3 k=1,ndet
        iocc=0
        do 2 j=1,nelec
 2        if(iworbd(j,k).eq.i) iocc=iocc+1
        if(iocc.ne.2) then
          iflag(1,i)=1
          goto 5
        endif
 3    continue
 5  continue

!   Omit empty orbitals

    do 6 i=1,local_ndetorb
     iflag(2,i)=0
     do 6 k=1,ndet
      do 6 j=1,nelec
 6      if(iworbd(j,k).eq.i) iflag(2,i)=1
    do 8 i=local_ndetorb+1,local_ndetorb+local_nadorb
     iflag(1,i)=1
 8   iflag(2,i)=0

    if(local_norbopt.eq.0.or.local_norbvirt.eq.0) then
      do 9 io=1,local_ndetorb
       do 9 jo=local_ncore+1,local_ndetorb+local_nadorb
 9      iwmix_virt(io,jo)=jo
    elseif(local_norbopt.ne.local_ndetorb.or.local_norbvirt.lt.local_nadorb) then
!     write(ounit,'(''get_norbterm: norbopt,ndetorb'',2i6)') local_norbopt,local_ndetorb
!     write(ounit,'(''get_norbterm: noptvirt,nadorb'',2i6)') local_norbvirt,local_nadorb
     call fatal_error('get_norbterm: Mixvirt block, inconsistent')
    endif

!  Orbital variation io -> io+a*jo
!  io: occupied orbitals in twf
!  jo: all orbitals
!  omitted if not same symmetry, or io empty, or both doubly occupied

    local_noporb=0

    ! if(iprt.gt.2) then
    !     write(ounit,*) '(''=========== establish no. orb variations =========='')'
    ! endif

    do 60 io=local_ncore+1,local_ndetorb
!   Omit empty orbitals
     if(iflag(2,io).eq.0) goto 60
     do 50 jo=local_ncore+1,local_ndetorb+local_nadorb
!   Omit if io and jo are the same
      if(io.eq.jo) goto 50
!   Omit if io and jo have different symmetry
      if(irrep(io).ne.irrep(jo)) goto 50
!   Omit if io and jo are both doubly occupied in all determinants
      if((iflag(1,io).eq.0).and.(iflag(1,jo).eq.0)) goto 50
!   Omit if io and jo are both active orbitals
      if(no_active.ne.0.and.iflag(1,io).ne.0.and.iflag(2,jo).ne.0) goto 50
!   Omit if we only want to mix according to the table mixvirt
      if(iwmix_virt(io,jo).eq.0) goto 50
!   Include: io is occupied in some determinant and jo not
      do 40 iab=1,2
        n0=0
        n1=nup
        if(iab.eq.2) then
          n0=nup
          n1=ndn
        endif
        m(iab)=0
        do 30 k=1,ndet
          do 15 ie=1,n1
            if(iworbd(ie+n0,k).eq.io) then
              iesave=ie
              goto 20
            endif
15         continue
          goto 30
20         continue
          do 25 ie=1,n1
25           if(iworbd(ie+n0,k).eq.jo) goto 30
          m(iab)=m(iab)+1
30       continue
40     continue
      if(m(1)+m(2).eq.0) then
!        if(iprt.gt.3) write(ounit,'(''no appropriate determinant for '',2i4)') io,jo
        goto 50
      endif

!   Define new operator (new variation) and its terms
      local_noporb=local_noporb+1


50    continue
60   continue

    local_norbterm=local_noporb
!    write(ounit,'(''number of orbital variations: '',i0)') local_norbterm

!   if mix_n, optorb_define called mutiple times with method=sr_n or lin_d
    if(method.eq.'linear') then
      nreduced=local_norbterm
    elseif(method.eq.'sr_n'.or.method.eq.'lin_d'.or.method.eq.'mix_n') then
      nreduced=1
    endif

    norbterm = local_norbterm
    return
end subroutine get_norbterm
