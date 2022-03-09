      module optwf_handle_wf
      use error, only: fatal_error
      use jastrow4_mod,       only: nterms4
      interface ! LAPACK interface
        SUBROUTINE dcopy(N,DX,INCX,DY,INCY)
!*  -- Reference BLAS level1 routine --
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
          INTEGER INCX,INCY,N
          DOUBLE PRECISION DX(*),DY(*)
        end subroutine
      end interface
      contains
c-----------------------------------------------------------------------
      subroutine write_wf(iwf_fit,iter)

      use mpiconf, only: idtask
      use mpi

      implicit none

      integer :: index, iter, iwf_fit
      character*40 filetype,wf,itn

      if(idtask.ne.0) return

      if(iter.lt.0) then
        filetype='_optimal.'//wf(1:index(wf,' ')-1)
       else
        write(wf,'(i1)') iwf_fit
        if(iter.lt.10) then
          write(itn,'(i1)') iter
         elseif(iter.lt.100) then
          write(itn,'(i2)') iter
         elseif(iter.lt.1000) then
          write(itn,'(i3)') iter
        endif
        filetype='_optimal.'//wf(1:index(wf,' ')-1)//'.iter'//itn(1:index(itn,' ')-1)
      endif

      call write_jastrow(iwf_fit,filetype)
      call write_lcao(iwf_fit,filetype)
      call write_ci(iwf_fit,filetype)

      return
      end
c-----------------------------------------------------------------------
      subroutine write_wf_best
      implicit none

      call restore_jastrow_best
      call restore_lcao_best
      call restore_ci_best
      call write_wf(1,-1)

      return
      end
c-----------------------------------------------------------------------
      subroutine write_jastrow(iwf_fit,filetype)

      use atom, only: nctype
      use jaspar, only: nspin1, nspin2
      use jaspar3, only: b, c, scalek
      use jaspar4, only: a4, norda, nordb, nordc
      use optwf_contrl, only: ioptjas
      use optwf_nparmj, only: nparma, nparmb, nparmc
      use contr2, only: ianalyt_lap, ijas, isc
      use precision_kinds, only: dp

      implicit none

      integer :: i, ict, index, iwf_fit, mparmja
      integer :: mparmjb, mparmjc
      character*50 fmt
      character*40 filename,filetype

      if(ioptjas.eq.0) return

      filename='jastrow'//filetype(1:index(filetype,' ')-1)

      open(2,file=filename,status='unknown')

      write(2,'(''&jastrow ianalyt_lap'',i2,'' ijas'',i2,'' isc'',i2,
     &'' nspin1'',i2,'' nspin2'',i2)') ianalyt_lap,ijas,isc,nspin1,nspin2
      write(2,*)
      write(2,'(''jastrow_parameter'',i4)') iwf_fit
      write(2,'(3i3,a28)') norda,nordb,nordc,' norda,nordb,nordc'
c tmp
      write(2,'(f13.8,a15)') scalek(1),' scalek'
      mparmja=2+max(0,norda-1)
      mparmjb=2+max(0,nordb-1)
      mparmjc=nterms4(nordc)
      if(mparmja.gt.0) then
        write(fmt,'(''(''i2,''f13.8,a28)'')') mparmja
       else
        write(fmt,'(''(a28)'')')
      endif
      do ict=1,nctype
        write(2,fmt) (a4(i,ict,1),i=1,mparmja),' (a(iparmj),iparmj=1,nparma)'
      enddo

      if(mparmjb.gt.0) then
        write(fmt,'(''(''i2,''f13.8,a28)'')') mparmjb
       else
        write(fmt,'(''(a28)'')')
      endif
      write(2,fmt) (b(i,1,1),i=1,mparmjb),' (b(iparmj),iparmj=1,nparmb)'

      if(mparmjc.gt.0) then
        write(fmt,'(''(''i2,''f13.8,a28)'')') mparmjc
       else
        write(fmt,'(''(a28)'')')
      endif
      do ict=1,nctype
        write(2,fmt) (c(i,ict,1),i=1,mparmjc),' (c(iparmj),iparmj=1,nparmc)'
      enddo
      write(2,'(''end'')')
      close(2)

      return
      end
c-----------------------------------------------------------------------
      subroutine write_lcao(iwf_fit,filetype)

      use numbas, only: numr
      use optwf_contrl, only: ioptorb
      use coefs, only: coef, nbasis, norb
      use orbval, only: nadorb
      use inputflags, only: scalecoef
      use precision_kinds, only: dp

      implicit none

      integer :: i, index, iwf_fit, j
      character*40 filename,filetype

      ! call resize_tensor(coef, norb+nadorb, 2)

      if(ioptorb.eq.0) return

      filename='orbitals'//filetype(1:index(filetype,' ')-1)
      open(2,file=filename,status='unknown')
      write(2,'(''lcao '',3i4)') norb+nadorb,nbasis,iwf_fit

      do i=1,norb+nadorb
        write(2,'(1000e20.8)') (coef(j,i,1)/scalecoef,j=1,nbasis)
      enddo

      write(2,'(''end'')')
      close(2)

      return
      end
c-----------------------------------------------------------------------
      subroutine write_ci(iwf_fit,filetype)

      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates
      use dets, only: cdet, ndet
      use multidet, only: kref
      use optwf_contrl, only: ioptci
      use dorb_m, only: iworbd
      use elec, only: nup
      use const, only: nelec

      implicit none

      integer :: i, index, istate, iwf_fit, j
      integer :: k, nmap, nptr, nterm
      character*40 filename,filetype

      if(ioptci.eq.0) return

      filename='det'//filetype(1:index(filetype,' ')-1)
      open(2,file=filename,status='unknown')

      write(2,'(''&electrons nelec '',i4,'' nup '',i4)') nelec,nup
      write(2,'(''# kref'',i4)') kref
      do istate=1,nstates
      write(2,'(''# State '',i4)') istate
      write(2,'(''determinants'',i10,i4)') ndet,iwf_fit
      write(2,'(100f15.8)') (cdet(i,istate,1),i=1,ndet)
      do k=1,ndet
       write(2,'(100i4)') (iworbd(i,k),i=1,nelec)
      enddo
      enddo

      write(2,'(''end'')')

      if(ncsf.ne.0) then
        write(2,'(''csf '',i10,i4)') ncsf,nstates
        do i=1,nstates
          write(2,'(100f15.8)') (ccsf(j,i,1),j=1,ncsf)
        enddo
        write(2,'(''end'')')
c
        nmap=0
        do i=1,ncsf
          nmap=nmap+ibdet(i)-iadet(i)+1
        enddo
        write(2,'(''csfmap'')')
        write(2,'(3i10)') ncsf,ndet,nmap
        nptr=0
        do i=1,ncsf
         nterm=ibdet(i)-iadet(i)+1
         write(2,'(i10)') ibdet(i)-iadet(i)+1
         do j=1,nterm
           nptr=nptr+1
           write(2,'(i10,f10.6)') icxdet(nptr),cxdet(nptr)
         enddo
        enddo
        write(2,'(''end'')')
      endif

      close(2)

      return
      end
c-----------------------------------------------------------------------
      subroutine setup_wf

      implicit none

      integer :: k

      do k=2,3
        call copy_jastrow(k)
        call copy_lcao(k)
        call copy_ci(k)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine save_wf

      use optwf_contrl, only: ioptci, ioptjas, ioptorb

      implicit none

      if(ioptjas.ne.0) call save_jastrow
      if(ioptorb.ne.0) call save_lcao
      if(ioptci.ne.0) call save_ci

      return
      end
c-----------------------------------------------------------------------
      subroutine restore_wf(iadiag)

      use optwf_contrl, only: ioptci, ioptjas, ioptorb

      implicit none

      integer :: iadiag

      if(ioptjas.ne.0) call restore_jastrow(iadiag)
      if(ioptorb.ne.0) call restore_lcao(iadiag)
      if(ioptci.ne.0) call restore_ci(iadiag)

      return
      end
c-----------------------------------------------------------------------
      subroutine save_wf_best(ioptjas,ioptorb,ioptci)

      implicit none

      integer :: ioptci, ioptjas, ioptorb

      if(ioptjas.ne.0) call save_jastrow_best
      if(ioptorb.ne.0) call save_lcao_best
      if(ioptci.ne.0) call save_ci_best

      return
      end
c-----------------------------------------------------------------------
      subroutine save_jastrow

      use precision_kinds, only: dp
      use vmc_mod, only: nordj1
      use atom, only: nctype, nctype_tot
      use wfsec, only: nwftype
      use jaspar3, only: b, c
      use jaspar4, only: a4, norda, nordb, nordc

      implicit none

      integer :: i, iadiag, ict, mparmja, mparmjb
      integer :: mparmjc
      real(dp), allocatable, save :: a4_save(:,:,:)
      real(dp), allocatable, save :: b_save(:,:,:)
      real(dp), allocatable, save :: c_save(:,:,:)

      ! dimension a4_save(nordj1,nctype_tot,MWF),b_save(nordj1,2,MWF),
      ! dimension c_save(83,nctype_tot,MWF)
      ! save a4_save,b_save,c_save

      save mparmja,mparmjb,mparmjc

      if(.not.allocated(a4_save)) allocate(a4_save(nordj1,nctype_tot,nwftype))
      if(.not.allocated(b_save)) allocate(b_save(nordj1,2,nwftype))
      if(.not.allocated(c_save)) allocate(c_save(83,nctype_tot,nwftype))

c Save parameters corresponding to run generating hessian

      mparmja=2+max(0,norda-1)
      mparmjb=2+max(0,nordb-1)
      mparmjc=nterms4(nordc)

      do ict=1,nctype
        do i=1,mparmja
          a4_save(i,ict,1)=a4(i,ict,1)
        enddo
      enddo
      do i=1,mparmjb
        b_save(i,1,1)=b(i,1,1)
      enddo
      do ict=1,nctype
        do i=1,mparmjc
          c_save(i,ict,1)=c(i,ict,1)
        enddo
      enddo

      return

      entry restore_jastrow(iadiag)

      if(.not.allocated(a4_save)) allocate(a4_save(nordj1,nctype_tot,nwftype))
      if(.not.allocated(b_save)) allocate(b_save(nordj1,2,nwftype))
      if(.not.allocated(c_save)) allocate(c_save(83,nctype_tot,nwftype))

c Restore parameters corresponding to run generating hessian
      do ict=1,nctype
        do i=1,mparmja
          a4(i,ict,iadiag)=a4_save(i,ict,1)
        enddo
      enddo
      do i=1,mparmjb
        b(i,1,iadiag)=b_save(i,1,1)
      enddo
      do ict=1,nctype
        do i=1,mparmjc
          c(i,ict,iadiag)=c_save(i,ict,1)
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine save_lcao

      use precision_kinds, only: dp
      use vmc_mod, only: norb_tot
      use coefs, only: coef, nbasis, norb
      use wfsec, only: nwftype

      implicit none

      integer :: i, iadiag, j
      real(dp), allocatable, save :: coef_save(:,:,:)

      if (.not. allocated(coef_save)) allocate(coef_save(nbasis, norb_tot, nwftype))
      ! dimension coef_save(nbasis,norb,MWF)
      ! save coef_save

      do i=1,norb
       do j=1,nbasis
        coef_save(j,i,1)=coef(j,i,1)
       enddo
      enddo

      return

      entry restore_lcao(iadiag)
      if (.not. allocated(coef_save)) allocate(coef_save(nbasis, norb_tot, nwftype))

      do i=1,norb
       do j=1,nbasis
        coef(j,i,iadiag)=coef_save(j,i,1)
       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine save_ci

      use precision_kinds, only: dp
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates
      use mstates_mod, only: MSTATES
      use dets, only: cdet, ndet
      use set_input_data, only: multideterminants_define

      implicit none

      integer :: i, iadiag, icsf, j, k
      integer :: kx
      real(dp), ALLOCATABLE, save :: cdet_save(:,:)
      real(dp), ALLOCATABLE, save :: ccsf_save(:,:)

      if(.not. allocated(cdet_save)) allocate(cdet_save(ndet,MSTATES))
      if(.not. allocated(ccsf_save)) allocate(ccsf_save(ndet,MSTATES))

      ! dimension cdet_save(ndet,nstates),ccsf_save(ndet,nstates)
      ! save cdet_save,ccsf_save

      do j=1,nstates
        do i=1,ndet
          cdet_save(i,j)=cdet(i,j,1)
        enddo
      enddo

      do j=1,nstates
       do icsf=1,ncsf
        ccsf_save(icsf,j)=ccsf(icsf,j,1)
       enddo
      enddo

      return

      entry restore_ci(iadiag)
      if(.not. allocated(cdet_save)) allocate(cdet_save(ndet,MSTATES))
      if(.not. allocated(ccsf_save)) allocate(ccsf_save(ndet,MSTATES))

      do j=1,nstates
        do i=1,ndet
          cdet(i,j,iadiag)=cdet_save(i,j)
        enddo
      enddo

      do j=1,nstates
       do icsf=1,ncsf
        ccsf(icsf,j,iadiag)=ccsf_save(icsf,j)
       enddo
      enddo

c if kref (iwdetorb, cxdet) has changed
      if(ncsf.gt.0) then
        do j=1,nstates
          do k=1,ndet
            cdet(k,j,iadiag)=0
          enddo
          do icsf=1,ncsf
            do k=iadet(icsf),ibdet(icsf)
              kx=icxdet(k)
              cdet(kx,j,iadiag)=cdet(kx,j,iadiag)+ccsf(icsf,j,iadiag)*cxdet(k)
            enddo
          enddo
        enddo

c reset kref=1
      call multideterminants_define(0,0)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_jastrow(iadiag)

      use atom, only: nctype
      use jaspar3, only: b, c, scalek
      use jaspar4, only: a4, norda, nordb, nordc

      implicit none

      integer :: i, iadiag, ict, mparmja, mparmjb
      integer :: mparmjc

      mparmja=2+max(0,norda-1)
      mparmjb=2+max(0,nordb-1)
      mparmjc=nterms4(nordc)

      scalek(iadiag)=scalek(1)
      do ict=1,nctype
        do i=1,mparmja
          a4(i,ict,iadiag)=a4(i,ict,1)
        enddo
      enddo
      do i=1,mparmjb
        b(i,1,iadiag)=b(i,1,1)
      enddo
      do ict=1,nctype
        do i=1,mparmjc
          c(i,ict,iadiag)=c(i,ict,1)
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine copy_lcao(iadiag)

      use coefs, only: coef, nbasis, norb
      use orbval, only: nadorb

      implicit none

      integer :: i, iadiag, j

      do i=1,norb+nadorb
       do j=1,nbasis
        coef(j,i,iadiag)=coef(j,i,1)
       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_ci(iadiag)

      use csfs, only: ccsf, ncsf, nstates
      use dets, only: cdet, ndet

      implicit none

      integer :: i, iadiag, icsf, j

      do j=1,nstates
        do i=1,ndet
          cdet(i,j,iadiag)=cdet(i,j,1)
        enddo
      enddo

      do j=1,nstates
       do icsf=1,ncsf
        ccsf(icsf,j,iadiag)=ccsf(icsf,j,1)
       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_zex(iadiag)

      use coefs, only: nbasis
      use basis, only: zex

      implicit none

      integer :: i, iadiag

      do i=1,nbasis
        zex(i,iadiag)=zex(i,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine save_jastrow_best

      use precision_kinds, only: dp
      use vmc_mod, only: nordj1
      use atom, only: nctype, nctype_tot
      use wfsec, only: nwftype
      use jaspar3, only: b, c
      use jaspar4, only: a4, norda, nordb, nordc

      implicit none

      integer :: i, ict, mparmja, mparmjb, mparmjc
      real(dp), allocatable, save :: a4_best(:,:,:)
      real(dp), allocatable, save :: b_best(:,:,:)
      real(dp), allocatable, save :: c_best(:,:,:)

      ! dimension a4_best(nordj1,nctype_tot,MWF),b_best(nordj1,2,MWF),
      ! dimension c_best(83,nctype_tot,MWF)
      ! save a4_best,b_best,c_best

      save mparmja,mparmjb,mparmjc

      if(.not.allocated(a4_best)) allocate(a4_best(nordj1,nctype_tot,nwftype))
      if(.not.allocated(b_best)) allocate(b_best(nordj1,2,nwftype))
      if(.not.allocated(c_best)) allocate(c_best(83,nctype_tot,nwftype))

c Save parameters corresponding to run generating hessian

      mparmja=2+max(0,norda-1)
      mparmjb=2+max(0,nordb-1)
      mparmjc=nterms4(nordc)

      do ict=1,nctype
        do i=1,mparmja
          a4_best(i,ict,1)=a4(i,ict,1)
        enddo
      enddo
      do i=1,mparmjb
        b_best(i,1,1)=b(i,1,1)
      enddo
      do ict=1,nctype
        do i=1,mparmjc
          c_best(i,ict,1)=c(i,ict,1)
        enddo
      enddo

      return

      entry restore_jastrow_best
      if(.not.allocated(a4_best)) allocate(a4_best(nordj1,nctype_tot,nwftype))
      if(.not.allocated(b_best)) allocate(b_best(nordj1,2,nwftype))
      if(.not.allocated(c_best)) allocate(c_best(83,nctype_tot,nwftype))

c Restore parameters corresponding to run generating hessian
      do ict=1,nctype
        do i=1,mparmja
          a4(i,ict,1)=a4_best(i,ict,1)
        enddo
      enddo
      do i=1,mparmjb
        b(i,1,1)=b_best(i,1,1)
      enddo
      do ict=1,nctype
        do i=1,mparmjc
          c(i,ict,1)=c_best(i,ict,1)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine save_lcao_best

      use precision_kinds, only: dp
      use vmc_mod, only: norb_tot
      use coefs, only: coef, nbasis, norb
      use wfsec, only: nwftype

      implicit none

      integer :: i, j
      real(dp), allocatable, save :: coef_best(:,:,:)

      if (.not. allocated(coef_best)) allocate(coef_best(nbasis, norb_tot, nwftype))
      ! dimension coef_best(nbasis,norb,MWF)
      ! save coef_best

      do i=1,norb
       do j=1,nbasis
        coef_best(j,i,1)=coef(j,i,1)
       enddo
      enddo

      return

      entry restore_lcao_best

c     if(ioptorb.eq.0) return
      if (.not. allocated(coef_best)) allocate(coef_best(nbasis, norb_tot, nwftype))
      do i=1,norb
       do j=1,nbasis
        coef(j,i,1)=coef_best(j,i,1)
       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine save_ci_best

      use precision_kinds, only: dp
      use csfs, only: ccsf, ncsf, nstates
      use csfs, only: cxdet, iadet, ibdet, icxdet
      use mstates_mod, only: MSTATES
      use dets, only: cdet, ndet
      use set_input_data, only: multideterminants_define

      implicit none

      integer :: i, icsf, j, k, kx
      real(dp), ALLOCATABLE, save :: cdet_best(:,:)
      real(dp), ALLOCATABLE, save :: ccsf_best(:,:)

      if(.not. allocated(cdet_best)) allocate(cdet_best(ndet,MSTATES))
      if(.not. allocated(ccsf_best)) allocate(ccsf_best(ndet,MSTATES))

      ! dimension cdet_best(ndet,nstates),ccsf_best(ndet,nstates)
      ! save cdet_best,ccsf_best

      do j=1,nstates
        do i=1,ndet
          cdet_best(i,j)=cdet(i,j,1)
        enddo
      enddo

      do j=1,nstates
       do icsf=1,ncsf
        ccsf_best(icsf,j)=ccsf(icsf,j,1)
       enddo
      enddo

      return

      entry restore_ci_best
      if(.not. allocated(cdet_best)) allocate(cdet_best(ndet,MSTATES))
      if(.not. allocated(ccsf_best)) allocate(ccsf_best(ndet,MSTATES))

c     if(ioptci.eq.0) return

      do j=1,nstates
        do i=1,ndet
          cdet(i,j,1)=cdet_best(i,j)
        enddo
      enddo

      do j=1,nstates
       do icsf=1,ncsf
        ccsf(icsf,j,1)=ccsf_best(icsf,j)
       enddo
      enddo

c if kref (iwdetorb, cxdet) has changed
      if(ncsf.gt.0) then
        do j=1,nstates
          do k=1,ndet
            cdet(k,j,1)=0
          enddo
          do icsf=1,ncsf
            do k=iadet(icsf),ibdet(icsf)
              kx=icxdet(k)
              cdet(kx,j,1)=cdet(kx,j,1)+ccsf(icsf,j,1)*cxdet(k)
            enddo
          enddo
        enddo

c reset kref=1
      call multideterminants_define(0,0)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_parameters(dparm,iflag,iadiag)

      use precision_kinds, only: dp

      implicit none

      integer :: iadiag, iflag
      real(dp), dimension(*) :: dparm

      iflag=0
      call compute_jastrow(dparm,iflag,iadiag)

      if(iflag.ne.0) return

      call compute_lcao(dparm,iadiag)
      call compute_ci(dparm,iadiag)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_jastrow(dparm,iflag,iadiag)

      use atom, only: nctype
      use jaspar3, only: b, c
      use jaspar4, only: a4
      use optwf_contrl, only: ioptjas
      use optwf_nparmj, only: nparma, nparmb, nparmc
      use optwf_wjas, only: iwjasa, iwjasb, iwjasc
      use precision_kinds, only: dp
      use cuspinit4_mod, only: cuspinit4
      use cuspexact4_mod, only: cuspexact4

      implicit none

      integer :: i, iadiag, ict, iflag, iparm
      real(dp), dimension(*) :: dparm

      if(ioptjas.eq.0) return

c Set up cusp conditions
      call cuspinit4(0)

c Add change to old parameters
      iparm=0
      do ict=1,nctype
        do i=1,nparma(ict)
          iparm=iparm+1
          a4(iwjasa(i,ict),ict,iadiag)=a4(iwjasa(i,ict),ict,iadiag)-dparm(iparm)
        enddo
      enddo
      do i=1,nparmb(1)
        iparm=iparm+1
        b(iwjasb(i,1),1,iadiag)=b(iwjasb(i,1),1,iadiag)-dparm(iparm)
      enddo
      do ict=1,nctype
        do i=1,nparmc(ict)
          iparm=iparm+1
          c(iwjasc(i,ict),ict,iadiag)=c(iwjasc(i,ict),ict,iadiag)-dparm(iparm)
        enddo
      enddo
      call cuspexact4(0,iadiag)

c Check parameters a2 and b2 > -scalek
      call check_parms_jas(iflag)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_lcao(dparm,iadiag)

      use vmc_mod, only: norb_tot
      use optwf_contrl, only: ioptorb
      use optwf_parms, only: nparmd, nparmj
      use coefs, only: coef, nbasis, norb
      use optorb_cblock, only: norbterm
      use orb_mat_022, only: ideriv
      use precision_kinds, only: dp

      implicit none

      integer :: i, iadiag, io, j, jo
      real(dp), dimension(nbasis, norb_tot) :: acoef
      real(dp), dimension(*) :: dparm

      if(ioptorb.eq.0) return

      do i=1,norb
       do j=1,nbasis
        acoef(j,i)=coef(j,i,iadiag)
       enddo
      enddo

c Update the orbitals
      do i=1,norbterm
       io=ideriv(1,i)
       jo=ideriv(2,i)
       do j=1,nbasis
        acoef(j,io)=acoef(j,io)-dparm(i+nparmj+nparmd)*coef(j,jo,iadiag)
       enddo
      enddo

      do i=1,norb
       do j=1,nbasis
        coef(j,i,iadiag)=acoef(j,i)
       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_ci(dparm,iadiag)

      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates
      use dets, only: cdet, ndet
      use optwf_contrl, only: ioptci, ioptjas, ioptorb
      use optwf_parms, only: nparmj
      use method_opt, only: method
      use precision_kinds, only: dp

      implicit none

      integer :: i, iadiag, icsf, idet, j
      integer :: jx, k
      real(dp) :: c90
      real(dp), dimension(*) :: dparm

      if(ioptci.eq.0) return

c Update the ci coef
      if((method.eq.'linear'.or.method.eq.'lin_d').and.ioptjas+ioptorb.eq.0) then
        do k=1,nstates

          if(ncsf.eq.0) then
            do idet=1,ndet
              cdet(idet,k,1)=dparm(idet+ndet*(k-1))
            enddo
           else
            do j=1,ndet
              cdet(j,k,1)=0
            enddo
            do icsf=1,ncsf
              do j=iadet(icsf),ibdet(icsf)
                jx=icxdet(j)
                cdet(jx,k,1)=cdet(jx,k,1)+dparm(icsf+ncsf*(k-1))*cxdet(j)
              enddo
              ccsf(icsf,k,1)=dparm(icsf+ncsf*(k-1))
            enddo
          endif

        enddo
       else
         if(ncsf.eq.0) then
           do idet=2,ndet
             cdet(idet,1,iadiag)=cdet(idet,1,iadiag)-dparm(idet-1+nparmj)
           enddo
          else
           do icsf=2,ncsf
             do j=iadet(icsf),ibdet(icsf)
               jx=icxdet(j)
               cdet(jx,1,iadiag)=cdet(jx,1,iadiag)-dparm(icsf-1+nparmj)*cxdet(j)
             enddo
             ccsf(icsf,1,iadiag)=ccsf(icsf,1,iadiag)-dparm(icsf-1+nparmj)
           enddo
         endif
      endif

c     do 90 j=1,nstates
c90     write(ounit,'(''csf ='',1000f20.15)') (ccsf(i,j,iadiag),i=1,ncsf)

      return
      end
c-----------------------------------------------------------------------
      subroutine check_parms_jas(iflag)

      use atom, only: nctype
      use jaspar3, only: b, scalek
      use jaspar4, only: a4
      use optwf_nparmj, only: nparma, nparmb
      use optwf_wjas, only: iwjasa, iwjasb
      use precision_kinds, only: dp
      use contrl_file,    only: ounit

      implicit none

      integer :: i, ict, iflag, iflaga, iflagb
      real(dp) :: scalem

      iflag=0
      iflaga=0
      iflagb=0

      scalem=-scalek(1)
      do ict=1,nctype
        do i=1,nparma(ict)
          if(iwjasa(i,ict).eq.2.and.a4(2,ict,1).le.scalem) iflaga=1
        enddo
      enddo
      if(iflaga.eq.1) then
        do ict=1,nctype
          write(ounit,'(''a2 < -scalek'',f10.5)') a4(2,ict,1)
        enddo
      endif
      do i=1,nparmb(1)
        if(iwjasb(i,1).eq.2.and.b(2,1,1).le.scalem) iflagb=1
      enddo
      if(iflagb.eq.1) write(ounit,'(''b2 < -scalek'',f10.5)') b(2,1,1)

      if(iflaga.eq.1.or.iflagb.eq.1) iflag=1

      return
      end
c-----------------------------------------------------------------------
      subroutine test_solution_parm(nparm,dparm,
     &              dparm_norm,dparm_norm_min,add_diag,iflag)
      use contrl_file,    only: ounit
      use precision_kinds, only: dp
      use contrl_file,    only: ounit

      implicit none

      integer :: i, iflag, nparm
      real(dp) :: add_diag, dparm_norm, dparm_norm_min
      real(dp), dimension(*) :: dparm

      iflag=0
      if(add_diag.le.0.d0) return

c Calculate rms change in parameters
      dparm_norm=0
      do i=1,nparm
        dparm_norm=dparm_norm+dparm(i)**2
      enddo
      dparm_norm=sqrt(dparm_norm/nparm)

      write(ounit,'(''dparm_norm,adiag ='',3g12.5)')
     &dparm_norm,add_diag

      if(dparm_norm.gt.dparm_norm_min) iflag=1

      return
      end
c-----------------------------------------------------------------------
      subroutine save_nparms

      use optwf_contrl, only: ioptci, ioptjas, ioptorb
      use optwf_parms, only: nparmd, nparmj
      use optorb_cblock, only: norbterm, nreduced
      use ci000, only: nciterm
      use contrl_file,    only: ounit
      use orbval, only: nadorb
      implicit none

      integer :: nciterm_sav, norbterm_sav, nparmd_sav, nparmj_sav, nreduced_sav, nadorb_sav

      save nparmj_sav,norbterm_sav,nciterm_sav,nparmd_sav,nreduced_sav, nadorb_sav

      nparmj_sav=nparmj
      norbterm_sav=norbterm
      nreduced_sav=nreduced
      nciterm_sav=nciterm
      nparmd=max(nciterm-1,0)
      nparmd_sav=nparmd
      nadorb_sav=nadorb

      write(ounit,'(''Saved max number of parameters, nparmj,norb,nciterm,nciterm-1: '',5i5)') nparmj,norbterm,nciterm,nparmd
      return

      entry set_nparms

      nparmj=nparmj_sav
      nparmd=nparmd_sav
      norbterm=norbterm_sav
      nreduced=nreduced_sav
      nciterm=nciterm_sav
      nadorb=nadorb_sav

      if(ioptjas.eq.0) nparmj=0
      if(ioptorb.eq.0) then
        norbterm=0
        nreduced=0
        nadorb=0
      endif
      if(ioptci.eq.0) then
        nciterm=0
        nparmd=0
      endif

      write(ounit,'(''Max number of parameters set to nparmj,norb,nciterm,nciterm-1: '',5i5)') nparmj,norbterm,nciterm,nparmd
      call set_nparms_tot

      return
      end
c-----------------------------------------------------------------------
      subroutine set_nparms_tot

      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      use optwf_parms, only: nparmd, nparmj
      use optorb_cblock, only: norbterm
      use ci000, only: nciterm
      use method_opt, only: method
      use contrl_file,    only: ounit
      implicit none

      integer :: i0

c Note: we do not vary the first (i0) CI coefficient unless a run where we only optimize the CI coefs

      if(method.eq.'sr_n') then

        nparmd=max(nciterm-1,0)
        nparm=nparmj+nparmd+norbterm

      elseif(method.eq.'linear'.or.method.eq.'lin_d' .or. method.eq.'mix_n') then

       i0=0
       if(ioptci.ne.0) i0=1
       if(ioptjas.eq.0.and.ioptorb.eq.0) i0=0

       nparmd=max(nciterm-1,0)
       nparm=nparmj+norbterm+nciterm-i0

      elseif(method.eq.'mix_n') then

       nparm=nparmj+norbterm+nciterm

      endif

      write(ounit,'(/,''number of parms: total, Jastrow, CI, orbitals= '',4i5)')
     & nparm,nparmj,nciterm,norbterm

      return
      end
c-----------------------------------------------------------------------
      subroutine optwf_store(l,wt,psid,energy)
c store elocal and derivatives of psi for each configuration (call in vmc)

      use sr_mod, only: mparm, mconf
      use optwf_parms, only: nparmj
      use csfs, only: nstates
      use derivjas, only: gvalue
      use optwf_contrl, only: ioptci, ioptjas, ioptorb
      use optwf_func, only: ifunc_omega
      use optwf_parms, only: nparmj
      use sr_mat_n, only: elocal, nconf_n, sr_ho
      use sr_mat_n, only: sr_o, wtg
      use deloc_dj_m, only: denergy
      use force_analy, only: iforce_analy
      use optorb_cblock, only: norbterm
      use orb_mat_001, only: orb_ho, orb_o
      use ci000, only: nciterm
      use ci001_blk, only: ci_o
      use ci003_blk, only: ci_e
      use method_opt, only: method
      !use optwf_sr_mod, only: izvzb, i_sr_rescale
      use precision_kinds, only: dp
      use optgeo_lib, only: force_store

      implicit none

      integer :: i0, ii, ijasci, istate, j
      integer :: l, ntmp
      integer :: izvzb, i_sr_rescale
      real(dp), dimension(nparmj) :: tmp_ho
      real(dp), dimension(*) :: wt
      real(dp), dimension(*) :: psid
      real(dp), dimension(*) :: energy

      if(iforce_analy.gt.0.and.izvzb.eq.1) call force_store(l)

      if((method.ne.'sr_n'.and.method.ne.'lin_d').or.ioptjas+ioptorb+ioptci.eq.0)return

      i0=1
      if(method.eq.'lin_d'.and.ioptjas+ioptorb.eq.0) i0=0

      if(l.gt.mconf) call fatal_error('SR_STORE: l gt mconf')

      if(nparmj /= 0) call dcopy(nparmj,gvalue,1,sr_o(1,l),1)

      ntmp=max(nciterm-i0,0)
      if (ntmp /= 0) call dcopy(ntmp,ci_o(1+i0),1,sr_o(nparmj+1,l),1)

      ijasci=nparmj+ntmp
      if(ijasci+nstates*norbterm+nstates.gt.mparm) call fatal_error('SR_STORE: iparm gt mparm')

      do istate=1,nstates
        ii=ijasci+(istate-1)*norbterm
        if (norbterm /= 0) call dcopy(norbterm,orb_o(1,istate),1,sr_o(ii+1,l),1)
        elocal(l,istate)=energy(istate)
        wtg(l,istate)=wt(istate)
      enddo

      ii=ijasci+nstates*norbterm
      do istate=1,nstates
        sr_o(ii+istate,l)=psid(istate)
      enddo

      nconf_n=l

      if (.not.allocated(sr_ho)) return
c TO FIX: we are assuming optjas.ne.0 or optorb.ne.0 -> Otherwise, standard secular problem
      do j=1,nparmj
        tmp_ho(j)=denergy(j,1)+gvalue(j)*energy(1)
      enddo

      if(nparmj /= 0) call dcopy(nparmj,tmp_ho,1,sr_ho(1,l),1)

      if(ntmp /= 0) call dcopy(ntmp,ci_e(1+i0),1,sr_ho(nparmj+1,l),1)

      if(norbterm /= 0) call dcopy(norbterm,orb_ho(1,1),1,sr_ho(nparmj+ntmp+1,l),1)

      return
      end
c-----------------------------------------------------------------------
      end module
