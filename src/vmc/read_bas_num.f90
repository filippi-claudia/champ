module read_bas_num_mod
      use error,   only: fatal_error
interface !LAPACK interface
  SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    INTEGER            INFO, LDA, LDB, N, NRHS
    INTEGER            IPIV( * )
    DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
  END SUBROUTINE
end interface
contains
      subroutine read_bas_num(iwf)

      use coefs,   only: nbasis
      use contrl_file, only: errunit,ounit
      use control, only: ipr
      use custom_broadcast, only: bcast
      use fitting_methods, only: exp_fit
      use general, only: bas_id,filename,filenames_bas_num,pooldir
     use mpiconf, only: wid
      use numbas,  only: allocate_numbas,arg,d2rwf,igrid,nr,nrbas,r0
      use numbas,  only: rwf, rmaxwf
      use numbas_mod, only: MRWF,MRWF_PTS
      use numexp,  only: ae,allocate_numexp,ce
      use precision_kinds, only: dp
      use pseudo,  only: nloc
      use spline2_mod, only: spline2
      use system,  only: atomtyp,nctype,nctype_tot,newghostype,znuc
      use vmc_mod, only: NCOEF

! c Written by Claudia Filippi
! c Modified by F. Schautz to use fancy file names
! c Reads in localized orbitals on a radial grid

      implicit none

      integer         :: ic, ir, irb, ii, jj, ll, icoef, iff
      integer         :: iwf, info
      real(dp)        :: val, dwf1, wfm, dwfn, dwfm, temp
      integer         :: iunit, iostat = 0, counter = 0
      real(dp)        :: cutoff_rmax = 1.0d-12
      logical         :: exist, skip = .true.

      real(dp), dimension(MRWF_PTS)       ::  x, work
      real(dp), dimension(ncoef)          ::  y
      real(dp), dimension(ncoef*ncoef)    ::  dmatr
      real(dp), dimension(nbasis)         ::  l
      integer, dimension(ncoef)           :: ipiv
      integer, dimension(nctype)          :: icusp

! c nrbas = number of numerical orbitals for each center
! c igrid = 1 linear r(i+1)=arg+r(i), r(1)=r0
! c         2 exponential r(i+1)=arg*r(i), r(1)=r0
! c         3 shifted exponential r(i+1)=r0*(arg**(i-1)-1)

      if (.not. allocated(filenames_bas_num)) allocate(filenames_bas_num(nctype))

      call allocate_numbas()
      call allocate_numexp()

      do ic=1,nctype+newghostype
        if (ic .gt. 999) call fatal_error('READ_BAS_NUM: atomtyp > 999')
        filename =  trim(pooldir) //  trim(bas_id) // ".basis." // atomtyp(ic)

        if (wid) then
          inquire(file=filename, exist=exist)
          if (exist) then
            open (newunit=iunit,file=filename, iostat=iostat, action='read', status='old')
            if (iostat .ne. 0) error stop "Problem in opening the numerical basis file"
          else
            call fatal_error( " Numerical Basis file "// filename // " does not exist.")
          endif
        endif

        write(ounit,*) '-----------------------------------------------------------------------'
        write(ounit,'(4a)')  " Reading numerical basis for ", trim(atomtyp(ic))," from the file :: ", trim(filename)
        write(ounit,*) '-----------------------------------------------------------------------'

        if (wid) read(iunit,*, iostat=iostat) nrbas(ic),igrid(ic),nr(ic),arg(ic),r0(ic),icusp(ic)
        call bcast(nrbas)
        call bcast(igrid)
        call bcast(nr)
        call bcast(arg)
        call bcast(r0)
        call bcast(icusp)

        write(45,*) "Reading the content ", nrbas(ic),igrid(ic),nr(ic),arg(ic),r0(ic),icusp(ic)
        write(45,*) "(Reading basis grid file = [ ",  trim(filename), " ] )"
        write(45,'(''center type '',i4,'' nrbas,igrid,nr,arg,r0 ='',2i4,i5,2f10.5)') &
        ic,nrbas(ic),igrid(ic),nr(ic),arg(ic),r0(ic)

        if(nrbas(ic).gt.MRWF) call fatal_error('READ_BAS_NUM: nrbas gt MRWF')
        if(nr(ic).gt.MRWF_PTS) call fatal_error('READ_BAS_NUM: nr gt MRWF_PTS')
        if(igrid(ic).ne.1.and.igrid(ic).ne.2.and.igrid(ic).ne.3) &
        call fatal_error('READ_BAS_NUM: grid not implemented')

        if(nloc.eq.0) then
          if (wid) read(iunit,*,iostat=iostat) (l(irb),irb=1,nrbas(ic))
          if (iostat .ne. 0) then
            write(errunit,'(a)') "Error:: Problem in reading the numerical basis file"
            write(errunit,'(2a)') "Stats for nerds :: in file ",__FILE__
            write(errunit,'(a,i6)') "at line ", __LINE__
          endif
          call bcast(l)
        endif
        if(nloc.eq.0) write(ounit,*) "nloc = 0 :: ", (l(irb),irb=1,nrbas(ic))

        if (wid) then
          do ir=1,nr(ic)
            read(iunit,*,iostat=iostat) x(ir),(rwf(ir,irb,ic,iwf),irb=1,nrbas(ic))
          enddo
          if (iostat .ne. 0) then
            write(errunit,'(a)') "Error:: Problem in reading the numerical basis file"
            write(errunit,'(2a)') "Stats for nerds :: in file ",__FILE__
            write(errunit,'(a,i6)') "at line ", __LINE__
          endif
        endif
        call bcast(x)
        call bcast(rwf)

!        Get the rmax value for each center. Set the cutoff to 10^-12

        if (.not. allocated(rmaxwf)) allocate (rmaxwf(nrbas(ic), nctype_tot))



        if(igrid(ic).eq.2.and.arg(ic).le.1.d0) arg(ic)=x(2)/x(1)
        if(igrid(ic).eq.3) r0(ic)=r0(ic)/(arg(ic)**(nr(ic)-1)-1.d0)
        call bcast(arg)
        call bcast(r0)

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
! c       endif

        dwf1=0.d0
        do icoef=2,NCOEF
          dwf1=dwf1+(icoef-1)*ce(icoef,irb,ic,iwf)*x(1)**(icoef-2)
        enddo

        ! Update the rmax at the point where rwf goes below cutoff (scanning from right to left)
        rmaxwf(irb, ic) = 20.0d0
        rloop: do ir=nr(ic),1,-1
          if (dabs(rwf(ir,irb,ic,iwf)) .gt. cutoff_rmax ) then
            rmaxwf(irb, ic) = x(ir)
            exit rloop
          endif
        enddo rloop

        write(45,'(a,i0,a,i0,a,g0)') "Initial rmax for center = ",ic, " basis = ",irb, " is ", rmaxwf(irb, ic)

! Nonzero basis at the boundary : Do exponential fitting.
        if(dabs(rmaxwf(irb,ic)-x(nr(ic))).lt.1.0d-10) then
          call exp_fit(x(nr(ic)-9:nr(ic)),rwf(nr(ic)-9:nr(ic),irb,ic,iwf), 10, ae(1,irb,ic,iwf), ae(2,irb,ic,iwf))
          if(ae(2,irb,ic,iwf).lt.0) call fatal_error ('BASIS_READ_NUM: ak<0')

          rmaxwf(irb,ic)=-dlog(cutoff_rmax/ae(1,irb,ic,iwf))/ae(2,irb,ic,iwf)

          write(45,'(a)') 'check the large radius expansion'
          write(45,'(a,g0,2x,g0)') 'Exponential fitting parameters : ', ae(1,irb,ic,iwf), ae(2,irb,ic,iwf)

          write(45,'(a,i0,a,i0,a,g0)') "Final rmax (fit) for center = ",ic, " basis = ",irb, " is ", rmaxwf(irb, ic)
          write(45, '(a)') 'irad,         rad                  rwf value            expo fit'
          do ir=1,10
            temp = ae(1,irb,ic,iwf)*dexp(-ae(2,irb,ic,iwf)*x(nr(ic)-ir))
            write(45,'(i3,2x,1p4e22.14)') ir,x(nr(ic)-ir),rwf(nr(ic)-ir,irb,ic,iwf), temp
          enddo

          dwfn=-ae(2,irb,ic,iwf)*ae(1,irb,ic,iwf)*dexp(-ae(2,irb,ic,iwf)*x(nr(ic)))

         else
          dwfn=0.d0
        endif
        write(45,*) 'dwf1,dwfn',dwf1,dwfn

        call spline2(x(:),rwf(:,irb,ic,iwf),nr(ic),dwf1,dwfn,d2rwf(:,irb,ic,iwf),work)

        enddo
      enddo

      return
      end subroutine read_bas_num

!  This is a temporary shift of the following subroutine

subroutine readps_gauss
  ! read 'Quantum-chemist' gauss pseudopotentials
  ! file format: one text file with basename gauss_ecp.dat
  !              for each atom type
  ! first line : arbitrary label (written to log-file)
  ! second line: number of projectors + 1 (i.e. total number of components)
  ! remaining lines: components in the order (local,L=0,L=1 ...)
  !     repeated for each component
  !        number terms
  !        repeated for each term in this component
  !          coefficient power exponent
  !
  ! NOTE: as usual power n means r**(n-2)
  !
      use contrl_file, only: errunit,ounit
      use custom_broadcast, only: bcast
      use gauss_ecp, only: allocate_gauss_ecp,ecp_coef,ecp_exponent
      use gauss_ecp, only: necp_power,necp_term
      use general, only: filename,filenames_ps_gauss,pooldir,pp_id
      use mpiconf, only: wid
      use precision_kinds, only: dp
      use pseudo,  only: lpot
      use pseudo_mod, only: MGAUSS,MPS_L,MPS_QUAD
      use qua,     only: nquad,wq,xq0,yq0,zq0
      use rotqua_mod, only: gesqua
      use system,  only: atomtyp,nctype


  implicit none

  integer         :: i, ic, idx, l
  integer         :: iunit, iostat, counter = 0
  logical         :: exist, skip = .true.
  real(dp) ::  necp_power_tmp
  
  character*80 label


  do ic=1,nctype
    if (wid) then
      if (nctype.gt.100) call fatal_error('READPS_GAUSS: nctype>100')
      filename =  trim(pooldir) // trim(pp_id) // ".gauss_ecp.dat." // atomtyp(ic)

      inquire(file=filename, exist=exist)
      if (exist) then
        open (newunit=iunit,file=filename, iostat=iostat, action='read', status='old')
        if (iostat .ne. 0) error stop "Problem in opening the pseudopotential file (Gaussian)"
      else
        call fatal_error( " Pseudopotential file (Gaussian) "// filename // " does not exist.")
      endif

    !   External file reading
      write(ounit,*) '-----------------------------------------------------------------------'
      write(ounit,'(4a)')  " Reading ECP pseudopotential for ", trim(atomtyp(ic))," from the file :: ", trim(filename)
      write(ounit,*) '-----------------------------------------------------------------------'

  ! label

      read(iunit,'(a80)',iostat=iostat) label
      if (iostat .ne. 0) then
        write(errunit,'(a)') "Error:: Problem in reading the pseudopotential file: label"
        write(errunit,'(2a)') "Stats for nerds :: in file ",__FILE__
        write(errunit,'(a,i6)') "at line ", __LINE__
      endif
      write(ounit,'(a,i4,a,a80)') 'ECP for atom type ', ic, ' label = ', adjustl(label)
    endif
    call bcast(label)
    ! max projector
    if (.not. allocated(lpot)) allocate (lpot(nctype))

    if (wid) then
      read(iunit,*,iostat=iostat) lpot(ic)
      if (iostat .ne. 0) then
        write(errunit,'(a)') "Error:: Problem in reading the pseudopotential file: lpot"
        write(errunit,'(2a)') "Stats for nerds :: in file ",__FILE__
        write(errunit,'(a,i6)') "at line ", __LINE__
      endif
      write(ounit,'(a,i4,a,i4)') 'ECP for atom type ', ic, ' lpot = ', lpot(ic)

      if(lpot(ic).gt.MPS_L) call fatal_error('READPS_GAUSS: increase MPS_L')
    endif
    call bcast(lpot)
  ! read terms of local part and all non-local parts
  ! local part first in file, but stored at index lpot
  ! non-local l=0 at index 1 etc, up to lpot-1

    call allocate_gauss_ecp()
    do l=1,lpot(ic)
        if(l.eq.1)then
          idx=lpot(ic)
          else
          idx=l-1
        endif
        if (wid) then
          read(iunit,*,iostat=iostat) necp_term(idx,ic)
          if (iostat .ne. 0) then
              write(errunit,'(a)') "Error:: Problem in reading the pseudopotential file: necp_term"
              write(errunit,'(2a)') "Stats for nerds :: in file ",__FILE__
              write(errunit,'(a,i6)') "at line ", __LINE__
          endif
        endif
        call bcast(necp_term)

        if(necp_term(idx,ic).gt.MGAUSS) call fatal_error('READPS_GAUSS: increase MGAUSS')

        write(ounit,'(a,2i6)') '    component, #terms ', l,necp_term(idx,ic)

        do i=1,necp_term(idx,ic)
          if (wid) then
            read(iunit,*,iostat=iostat) ecp_coef(i,idx,ic), necp_power_tmp, ecp_exponent(i,idx,ic)

            if (iostat .ne. 0) then
              write(errunit,'(a)') "Error:: Problem in reading the pseudopotential file: ecp_coeff, power, ecp_exponents"
              write(errunit,'(2a)') "Stats for nerds :: in file ",__FILE__
              write(errunit,'(a,i6)') "at line ", __LINE__
           endif
           necp_power(i,idx,ic)=int(necp_power_tmp)
            write(ounit,'(a,f16.8,i2,f16.8)') '    coef, power, expo ', ecp_coef(i,idx,ic), &
                                                    necp_power(i,idx,ic), ecp_exponent(i,idx,ic)
          endif
        enddo
        call bcast(ecp_coef)
        call bcast(necp_power)
        call bcast(ecp_exponent)
    enddo

    if (wid) close(iunit)
  enddo

  if (.not. allocated(wq)) allocate (wq(MPS_QUAD))
  if (.not. allocated(xq0)) allocate (xq0(MPS_QUAD))
  if (.not. allocated(yq0)) allocate (yq0(MPS_QUAD))
  if (.not. allocated(zq0)) allocate (zq0(MPS_QUAD))

  call gesqua(nquad,xq0,yq0,zq0,wq)
  call bcast(wq)
  call bcast(xq0)
  call bcast(yq0)
  call bcast(zq0)
  return
end subroutine readps_gauss

! subroutine trexio_read_ecp_gaussian
!   !> This subroutine reads the ECP pseudopotential in Gaussian format
!   !! from the trexio hdf5/text file.
!   !! first line : arbitrary label (written to log-file)
!   !! second line: number of projectors + 1 (i.e. total number of components)
!   !! remaining lines: components in the order (local,L=0,L=1 ...)
!   !!     repeated for each component
!   !!        number terms
!   !!        repeated for each term in this component
!   !!          coefficient power exponent
!   !!
!   !! NOTE: as usual power n means r**(n-2)
!   !!
!   !! @author Ravindra Shinde (r.l.shinde@utwente.nl)
!   use custom_broadcast,   only: bcast
!   use mpiconf,            only: wid

!   use pseudo_mod, only: MPS_L, MGAUSS, MPS_QUAD
!   use system, only: nctype, atomtyp
!   use gauss_ecp, only: ecp_coef, ecp_exponent, necp_power, necp_term
!   use gauss_ecp, only: allocate_gauss_ecp
!   use pseudo, only: lpot
!   use qua, only: nquad, wq, xq0, yq0, zq0
!   use general, only: pooldir, filename, pp_id, filenames_ps_gauss
!   use contrl_file,        only: ounit, errunit

!   use precision_kinds, only: dp
!   use trexio
!   implicit none

!   integer                   :: i, ic, idx, l
!   integer                   :: iunit, iostat, counter = 0
!   logical                   :: exist, skip = .true.
!   character(:), allocatable :: label

!   ! trexio
!   character(:), allocatable :: file_trexio, file_trexio_path
!   integer(8)                :: trex_ecp_file
!   integer                   :: rc = 1


!   if((file_trexio(1:6) == '$pool/') .or. (file_trexio(1:6) == '$POOL/')) then
!     file_trexio_path = pooldir // file_trexio(7:)
!   else
!     file_trexio_path = file_trexio
!   endif



!   do ic=1,nctype
!     if (wid) then
!       if (nctype.gt.100) call fatal_error('READPS_GAUSS: nctype>100')
!       filename =  trim(pooldir) // trim(pp_id) // ".gauss_ecp.dat." // atomtyp(ic)

!       inquire(file=filename, exist=exist)
!       if (exist) then
!         open (newunit=iunit,file=filename, iostat=iostat, action='read', status='old')
!         if (iostat .ne. 0) error stop "Problem in opening the pseudopotential file (Gaussian)"
!       else
!         call fatal_error( " Pseudopotential file (Gaussian) "// filename // " does not exist.")
!       endif

!     !   External file reading
!       write(ounit,*) '-----------------------------------------------------------------------'
!       write(ounit,'(4a)')  " Reading ECP pseudopotential for ", trim(atomtyp(ic))," from the file :: ", trim(filename)
!       write(ounit,*) '-----------------------------------------------------------------------'

!   ! label

!       read(iunit,'(a80)',iostat=iostat) label
!       if (iostat .ne. 0) then
!         write(errunit,'(a)') "Error:: Problem in reading the pseudopotential file: label"
!         write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!       endif
!       write(ounit,'(a,i4,a,a80)') 'ECP for atom type ', ic, ' label = ', adjustl(label)
!     endif
!     call bcast(label)
!     ! max projector
!     if (.not. allocated(lpot)) allocate (lpot(nctype))

!     if (wid) then
!       read(iunit,*,iostat=iostat) lpot(ic)
!       if (iostat .ne. 0) then
!         write(errunit,'(a)') "Error:: Problem in reading the pseudopotential file: lpot"
!         write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!       endif
!       write(ounit,'(a,i4,a,i4)') 'ECP for atom type ', ic, ' lpot = ', lpot(ic)

!       if(lpot(ic).gt.MPS_L) call fatal_error('READPS_GAUSS: increase MPS_L')
!     endif
!     call bcast(lpot)
!   ! read terms of local part and all non-local parts
!   ! local part first in file, but stored at index lpot
!   ! non-local l=0 at index 1 etc, up to lpot-1

!     call allocate_gauss_ecp()
!     do l=1,lpot(ic)
!         if(l.eq.1)then
!           idx=lpot(ic)
!           else
!           idx=l-1
!         endif
!         if (wid) then
!           read(iunit,*,iostat=iostat) necp_term(idx,ic)
!           if (iostat .ne. 0) then
!               write(errunit,'(a)') "Error:: Problem in reading the pseudopotential file: necp_term"
!               write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!           endif
!         endif
!         call bcast(necp_term)

!         if(necp_term(idx,ic).gt.MGAUSS) call fatal_error('READPS_GAUSS: increase MGAUSS')

!         write(ounit,'(a,2i6)') '    component, #terms ', l,necp_term(idx,ic)

!         do i=1,necp_term(idx,ic)
!           if (wid) then
!             read(iunit,*,iostat=iostat) ecp_coef(i,idx,ic), necp_power(i,idx,ic),ecp_exponent(i,idx,ic)

!             if (iostat .ne. 0) then
!               write(errunit,'(a)') "Error:: Problem in reading the pseudopotential file: ecp_coeff, power, ecp_exponents"
!               write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!             endif
!             write(ounit,'(a,f16.8,i2,f16.8)') '    coef, power, expo ', ecp_coef(i,idx,ic),necp_power(i,idx,ic), ecp_exponent(i,idx,ic)
!           endif
!         enddo
!         call bcast(ecp_coef)
!         call bcast(necp_power)
!         call bcast(ecp_exponent)
!     enddo

!     if (wid) close(iunit)
!   enddo

!   if (.not. allocated(wq)) allocate (wq(MPS_QUAD))
!   if (.not. allocated(xq0)) allocate (xq0(MPS_QUAD))
!   if (.not. allocated(yq0)) allocate (yq0(MPS_QUAD))
!   if (.not. allocated(zq0)) allocate (zq0(MPS_QUAD))

!   call gesqua(nquad,xq0,yq0,zq0,wq)
!   call bcast(wq)
!   call bcast(xq0)
!   call bcast(yq0)
!   call bcast(zq0)
!   return
! end subroutine readps_gauss

end module
