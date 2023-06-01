      module optwf_handle_wf
      use error,   only: fatal_error
      use jastrow4_mod, only: nterms4
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

      use mpi
      use mpiconf, only: idtask

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

      use system, only: nctype
      use jastrow, only: nspin1, nspin2, b, c, scalek, a4, norda, nordb, nordc, ianalyt_lap, ijas, isc
      use bparm, only: nspin2b
      use optwf_control, only: ioptjas
      use optwf_nparmj, only: nparma, nparmb, nparmc
      use vmc_mod, only: nwftypejas, nstoj, extraj, jtos
      use precision_kinds, only: dp
      use system,  only: nctype

      implicit none

      integer :: i, isp, ict, index, iwf_fit, mparmja, k
      integer :: mparmjb, mparmjc
      character*50 fmt, temp
      character*40 filename,filetype

      if(ioptjas.eq.0) return

      filename='jastrow'//filetype(1:index(filetype,' ')-1)
      open(2,file=filename,status='unknown')
      write(2,'(''jastrow_parameter'',i4)') iwf_fit
      write(2,'(3i3,a28)') norda,nordb,nordc,' norda,nordb,nordc'
c tmp
      write(2,'(f13.8,a15)') scalek(1),' scalek'
      do k=1,nwftypejas
c        if (extraj.eq.1) write(2,'(''jastrows_to_states'',i6,<nstoj(k)>i4)')
c     &                          nstoj(k), (jtos(k,i),i=1,nstoj(k))             ! Intel version
        if (extraj.eq.1) write(temp,'(a,i0,a,a)') '(a,1x,i0,1x,',nstoj(k),'(1x,i0,1x)',')'
        if (extraj.eq.1) write(2,temp) "jastrows_to_states",nstoj(k),(jtos(k,i),i=1,nstoj(k)) ! GNU version

        mparmja=2+max(0,norda-1)
        mparmjb=2+max(0,nordb-1)
        mparmjc=nterms4(nordc)
        if(mparmja.gt.0) then
          write(fmt,'(''(''i2,''f13.8,a28)'')') mparmja
        else
          write(fmt,'(''(a28)'')')
        endif
        do ict=1,nctype
          write(2,fmt) (a4(i,ict,k),i=1,mparmja),' (a(iparmj),iparmj=1,nparma)'
        enddo

        if(mparmjb.gt.0) then
          write(fmt,'(''(''i2,''f13.8,a28)'')') mparmjb
        else
          write(fmt,'(''(a28)'')')
        endif
        do isp=1,nspin2b
          write(2,fmt) (b(i,isp,k),i=1,mparmjb),' (b(iparmj),iparmj=1,nparmb)'
        enddo

        if(mparmjc.gt.0) then
          write(fmt,'(''(''i2,''f13.8,a28)'')') mparmjc
        else
          write(fmt,'(''(a28)'')')
        endif
        do ict=1,nctype
          write(2,fmt) (c(i,ict,k),i=1,mparmjc),' (c(iparmj),iparmj=1,nparmc)'
        enddo
      enddo

      write(2,'(''end'')')
      close(2)
      return
      end
c-----------------------------------------------------------------------
      subroutine write_lcao(iwf_fit,filetype)

      use numbas, only: numr
      use optwf_control, only: ioptorb
      use coefs, only: nbasis
      use slater, only: norb, coef
      use orbval, only: nadorb
      use inputflags, only: scalecoef
      use vmc_mod, only: nwftypeorb, otos, nstoo, extrao
      use precision_kinds, only: dp
      use slater,  only: coef,norb

      implicit none

      integer :: i, index, iwf_fit, j, k
      character*40 filename,filetype, temp

      ! call resize_tensor(coef, norb+nadorb, 2)

      if(ioptorb.eq.0) return

      filename='orbitals'//filetype(1:index(filetype,' ')-1)
      open(2,file=filename,status='unknown')
      write(2,'(''lcao '',3i4)') norb+nadorb,nbasis,iwf_fit

      do k=1,nwftypeorb
c        if (extrao.eq.1) write(2,'(''orbitals_to_states'',i6,<nstoo(k)>i4)')
c     &                         nstoo(k), (otos(k,i),i=1,nstoo(k))           ! Intel version

        if (extrao.eq.1) write(temp,'(a,i0,a,a)') '(a,1x,i0,1x,',nstoo(k),'(1x,i0,1x)',')'
        if (extrao.eq.1) write(2,temp) "orbitals_to_states",nstoo(k),(otos(k,i),i=1,nstoo(k)) ! GNU version

        do i=1,norb+nadorb
          write(2,'(1000e20.8)') (coef(j,i,k)/scalecoef,j=1,nbasis)
        enddo
      enddo

      write(2,'(''end'')')
      close(2)

      return
      end
c-----------------------------------------------------------------------
      subroutine write_ci(iwf_fit,filetype)

      use csfs,    only: ccsf,cxdet,iadet,ibdet,icxdet,ncsf,nstates
      use dorb_m,  only: iworbd
      use optwf_control, only: ioptci
      use slater,  only: cdet,kref,ndet
      use system,  only: nelec,nup

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

      use optwf_control, only: ioptci,ioptjas,ioptorb

      implicit none

      if(ioptjas.ne.0) call save_jastrow
      if(ioptorb.ne.0) call save_lcao
      if(ioptci.ne.0) call save_ci

      return
      end
c-----------------------------------------------------------------------
      subroutine restore_wf(iadiag)

      use optwf_control, only: ioptci,ioptjas,ioptorb

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

      use bparm,   only: nspin2b
      use jastrow, only: norda,nordb,nordc
      use jastrow, only: a4,b,c,nordj1
      use multiple_geo, only: nwftype
      use precision_kinds, only: dp
      use system, only: nctype, nctype_tot
      use multiple_geo, only: nwftype
      use jastrow, only: b, c, a4, norda, nordb, nordc, nordj1
      use bparm, only: nspin2b
      use vmc_mod, only: nwftypejas
      use optwf_control, only: method

      implicit none

      integer :: i, isp, iadiag, ict, mparmja, mparmjb
      integer :: mparmjc, k
      real(dp), allocatable, save :: a4_save(:,:,:)
      real(dp), allocatable, save :: b_save(:,:,:)
      real(dp), allocatable, save :: c_save(:,:,:)

      ! dimension a4_save(nordj1,nctype_tot,MWF),b_save(nordj1,2,MWF),
      ! dimension c_save(83,nctype_tot,MWF)
      ! save a4_save,b_save,c_save

      save mparmja,mparmjb,mparmjc

      !if(.not.allocated(a4_save)) allocate(a4_save(nordj1,nctype_tot,nwftype))
      !if(.not.allocated(b_save)) allocate(b_save(nordj1,2,nwftype))
      !if(.not.allocated(c_save)) allocate(c_save(83,nctype_tot,nwftype))

c Save parameters corresponding to run generating hessian

      mparmja=2+max(0,norda-1)
      mparmjb=2+max(0,nordb-1)
      mparmjc=nterms4(nordc)

      if (method.eq.'sr_n'.and.nwftypejas.gt.1) then
        if(.not.allocated(a4_save)) allocate(a4_save(nordj1,nctype_tot,nwftypejas))
        if(.not.allocated(b_save)) allocate(b_save(nordj1,2,nwftypejas))
        if(.not.allocated(c_save)) allocate(c_save(83,nctype_tot,nwftypejas))
        do k=1,nwftypejas
          do ict=1,nctype
            do i=1,mparmja
              a4_save(i,ict,k)=a4(i,ict,k)
            enddo
          enddo
          do isp=1,nspin2b
            do i=1,mparmjb
              b_save(i,isp,k)=b(i,isp,k)
            enddo
          enddo
          do ict=1,nctype
            do i=1,mparmjc
              c_save(i,ict,k)=c(i,ict,k)
            enddo
          enddo
        enddo

      else
        if(.not.allocated(a4_save)) allocate(a4_save(nordj1,nctype_tot,nwftype))
        if(.not.allocated(b_save)) allocate(b_save(nordj1,2,nwftype))
        if(.not.allocated(c_save)) allocate(c_save(83,nctype_tot,nwftype))

        do ict=1,nctype
          do i=1,mparmja
            a4_save(i,ict,1)=a4(i,ict,1)
          enddo
        enddo
        do isp=1,nspin2b
          do i=1,mparmjb
            b_save(i,isp,1)=b(i,isp,1)
          enddo
        enddo
        do ict=1,nctype
          do i=1,mparmjc
            c_save(i,ict,1)=c(i,ict,1)
          enddo
        enddo

      endif

      return

      entry restore_jastrow(iadiag)

      if (method.eq.'sr_n'.and.nwftypejas.gt.1) then
        if(.not.allocated(a4_save)) allocate(a4_save(nordj1,nctype_tot,nwftypejas))
        if(.not.allocated(b_save)) allocate(b_save(nordj1,2,nwftypejas))
        if(.not.allocated(c_save)) allocate(c_save(83,nctype_tot,nwftypejas))

c Restore parameters corresponding to run generating hessian
        do k=1,nwftypejas
          do ict=1,nctype
            do i=1,mparmja
              a4(i,ict,k)=a4_save(i,ict,k)
            enddo
          enddo
          do isp=1,nspin2b
            do i=1,mparmjb
              b(i,isp,k)=b_save(i,isp,k)
            enddo
          enddo
          do ict=1,nctype
            do i=1,mparmjc
              c(i,ict,k)=c_save(i,ict,k)
            enddo
          enddo
        enddo

      else

        if(.not.allocated(a4_save)) allocate(a4_save(nordj1,nctype_tot,nwftype))
        if(.not.allocated(b_save)) allocate(b_save(nordj1,2,nwftype))
        if(.not.allocated(c_save)) allocate(c_save(83,nctype_tot,nwftype))

c Restore parameters corresponding to run generating hessian
        do ict=1,nctype
          do i=1,mparmja
            a4(i,ict,iadiag)=a4_save(i,ict,1)
          enddo
        enddo
        do isp=1,nspin2b
          do i=1,mparmjb
            b(i,isp,iadiag)=b_save(i,isp,1)
          enddo
        enddo
        do ict=1,nctype
          do i=1,mparmjc
            c(i,ict,iadiag)=c_save(i,ict,1)
          enddo
        enddo
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine save_lcao

      use coefs,   only: nbasis
      use multiple_geo, only: nwftype
      use precision_kinds, only: dp
      use vmc_mod, only: norb_tot, nwftypeorb
      use coefs, only: nbasis
      use slater, only: norb, coef
      use multiple_geo, only: nwftype
      use optwf_control, only: method

      implicit none

      integer :: i, iadiag, j, k
      real(dp), allocatable, save :: coef_save(:,:,:)

      if (method.eq.'sr_n'.and.nwftypeorb.gt.1) then
        if (.not. allocated(coef_save)) allocate(coef_save(nbasis, norb_tot, nwftypeorb))
        ! dimension coef_save(nbasis,norb,MWF)
        ! save coef_save

        do k=1,nwftypeorb
          do i=1,norb
            do j=1,nbasis
              coef_save(j,i,k)=coef(j,i,k)
            enddo
          enddo
        enddo

      else

        if (.not. allocated(coef_save)) allocate(coef_save(nbasis, norb_tot, nwftype))
        ! dimension coef_save(nbasis,norb,MWF)
        ! save coef_save

        do i=1,norb
          do j=1,nbasis
            coef_save(j,i,1)=coef(j,i,1)
          enddo
        enddo

      endif

      return

      entry restore_lcao(iadiag)

      if (method.eq.'sr_n'.and.nwftypeorb.gt.1) then
        if (.not. allocated(coef_save)) allocate(coef_save(nbasis, norb_tot, nwftypeorb))

        do k=1,nwftypeorb
          do i=1,norb
            do j=1,nbasis
              coef(j,i,k)=coef_save(j,i,k)
            enddo
          enddo
        enddo

      else

        if (.not. allocated(coef_save)) allocate(coef_save(nbasis, norb_tot, nwftype))

        do i=1,norb
          do j=1,nbasis
            coef(j,i,iadiag)=coef_save(j,i,1)
          enddo
        enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine save_ci

      use csfs,    only: ccsf,cxdet,iadet,ibdet,icxdet,ncsf,nstates
      use mstates_mod, only: MSTATES
      use precision_kinds, only: dp
      use set_input_data, only: multideterminants_define
      use slater,  only: cdet,ndet

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

      use system, only: nctype
      use jastrow, only: b, c, scalek, a4, norda, nordb, nordc
      use bparm, only: nspin2b
      use optwf_control, only: method
      use vmc_mod, only: nwftypejas

      implicit none

      integer :: i, isp, iadiag, ict, mparmja, mparmjb
      integer :: mparmjc, k

      mparmja=2+max(0,norda-1)
      mparmjb=2+max(0,nordb-1)
      mparmjc=nterms4(nordc)

      scalek(iadiag)=scalek(1)

      if (method.eq.'sr_n'.and.nwftypejas.gt.1) then
        do k=1,nwftypejas
          do ict=1,nctype
            do i=1,mparmja
              a4(i,ict,k)=a4(i,ict,k)
            enddo
          enddo
          do isp=1,nspin2b
            do i=1,mparmjb
              b(i,isp,k)=b(i,isp,k)
            enddo
          enddo
          do ict=1,nctype
            do i=1,mparmjc
              c(i,ict,k)=c(i,ict,k)
            enddo
          enddo
        enddo

      else

        do ict=1,nctype
          do i=1,mparmja
            a4(i,ict,iadiag)=a4(i,ict,1)
          enddo
        enddo
        do isp=1,nspin2b
          do i=1,mparmjb
            b(i,isp,iadiag)=b(i,isp,1)
          enddo
        enddo
        do ict=1,nctype
          do i=1,mparmjc
            c(i,ict,iadiag)=c(i,ict,1)
          enddo
        enddo

      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine copy_lcao(iadiag)

      use coefs, only: nbasis
      use slater, only: norb, coef
      use orbval, only: nadorb
      use optwf_control, only: method
      use vmc_mod, only: nwftypeorb

      implicit none

      integer :: i, iadiag, j, k

      if (method.eq.'sr_n'.and.nwftypeorb.gt.1) then
        do k=1,nwftypeorb
          do i=1,norb+nadorb
            do j=1,nbasis
              coef(j,i,k)=coef(j,i,k)
            enddo
          enddo
        enddo

      else

        do i=1,norb+nadorb
          do j=1,nbasis
            coef(j,i,iadiag)=coef(j,i,1)
          enddo
        enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_ci(iadiag)

      use csfs,    only: ccsf,ncsf,nstates
      use slater,  only: cdet,ndet

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

      use basis,   only: zex
      use coefs,   only: nbasis

      implicit none

      integer :: i, iadiag

      do i=1,nbasis
        zex(i,iadiag)=zex(i,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine save_jastrow_best

      use bparm,   only: nspin2b
      use jastrow, only: norda,nordb,nordc
      use jastrow, only: a4,b,c,nordj1
      use multiple_geo, only: nwftype
      use precision_kinds, only: dp
      use system,  only: nctype,nctype_tot

      implicit none

      integer :: i, isp, ict, mparmja, mparmjb, mparmjc
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
      do isp=1,nspin2b
        do i=1,mparmjb
          b_best(i,isp,1)=b(i,isp,1)
        enddo
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
      do isp=1,nspin2b
        do i=1,mparmjb
          b(i,isp,1)=b_best(i,isp,1)
        enddo
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

      use coefs,   only: nbasis
      use multiple_geo, only: nwftype
      use precision_kinds, only: dp
      use vmc_mod, only: norb_tot, nwftypeorb
      use optwf_control, only: method
      use coefs, only: nbasis
      use slater, only: norb, coef
      use multiple_geo, only: nwftype

      implicit none

      integer :: i, j, k
      real(dp), allocatable, save :: coef_best(:,:,:)

      !if (.not. allocated(coef_best)) allocate(coef_best(nbasis, norb_tot, nwftype))
      ! dimension coef_best(nbasis,norb,MWF)
      ! save coef_best

      !remove multistate, allocate nwftype, remove loops over nwftypeorb, replace, 3rd index with 1.
      if (method.eq.'sr_n'.and.nwftypeorb.gt.1) then
        if (.not. allocated(coef_best)) allocate(coef_best(nbasis, norb_tot, nwftypeorb))
        do k=1,nwftypeorb
          do i=1,norb
            do j=1,nbasis
              coef_best(j,i,k)=coef(j,i,k)
            enddo
          enddo
        enddo

      else

        if (.not. allocated(coef_best)) allocate(coef_best(nbasis, norb_tot, nwftype))
        do i=1,norb
          do j=1,nbasis
            coef_best(j,i,1)=coef(j,i,1)
          enddo
        enddo

      endif

      return

      entry restore_lcao_best

c     if(ioptorb.eq.0) return
      if (method.eq.'sr_n'.and.nwftypeorb.gt.1) then

        if (.not. allocated(coef_best)) allocate(coef_best(nbasis, norb_tot, nwftypeorb))
        do k=1, nwftypeorb
          do i=1,norb
            do j=1,nbasis
              coef(j,i,k)=coef_best(j,i,k)
            enddo
          enddo
        enddo

      else

        if (.not. allocated(coef_best)) allocate(coef_best(nbasis, norb_tot, nwftype))
        do i=1,norb
          do j=1,nbasis
            coef(j,i,1)=coef_best(j,i,1)
          enddo
        enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine save_ci_best

      use csfs,    only: ccsf,cxdet,iadet,ibdet,icxdet,ncsf,nstates
      use mstates_mod, only: MSTATES
      use precision_kinds, only: dp
      use set_input_data, only: multideterminants_define
      use slater,  only: cdet,ndet

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

      use system, only: nctype
      use jastrow, only: b, c, a4
      use bparm, only: nspin2b
      use optwf_control, only: ioptjas, method
      use optwf_nparmj, only: nparma, nparmb, nparmc
      use optwf_wjas, only: iwjasa, iwjasb, iwjasc
      use precision_kinds, only: dp
      use cuspinit4_mod, only: cuspinit4
      use cuspexact4_mod, only: cuspexact4
      use vmc_mod, only: nwftypejas, stoj
      use sr_mat_n, only: sr_state

      implicit none

      integer :: i, isp, iadiag, ict, iflag, iparm, j!, k
      real(dp), dimension(*) :: dparm

      if(ioptjas.eq.0) return

c Set up cusp conditions
      call cuspinit4(0)

c Add change to old parameters
      if (method.eq.'sr_n'.and.nwftypejas.gt.1) then
          j=stoj(sr_state)
          iparm=0
          do ict=1,nctype
            do i=1,nparma(ict)
              iparm=iparm+1
              a4(iwjasa(i,ict),ict,j)=a4(iwjasa(i,ict),ict,j)-dparm(iparm)
            enddo
          enddo
          do isp=1,nspin2b
            do i=1,nparmb(1)
              iparm=iparm+1
              b(iwjasb(i,isp),isp,j)=b(iwjasb(i,isp),isp,j)-dparm(iparm)
            enddo
          enddo
          do ict=1,nctype
            do i=1,nparmc(ict)
              iparm=iparm+1
              c(iwjasc(i,ict),ict,j)=c(iwjasc(i,ict),ict,j)-dparm(iparm)
            enddo
          enddo

      else

        iparm=0
        do ict=1,nctype
          do i=1,nparma(ict)
            iparm=iparm+1
            a4(iwjasa(i,ict),ict,iadiag)=a4(iwjasa(i,ict),ict,iadiag)-dparm(iparm)
          enddo
        enddo
        do isp=1,nspin2b
          do i=1,nparmb(1)
            iparm=iparm+1
            b(iwjasb(i,isp),isp,iadiag)=b(iwjasb(i,isp),isp,iadiag)-dparm(iparm)
          enddo
        enddo
        do ict=1,nctype
          do i=1,nparmc(ict)
            iparm=iparm+1
            c(iwjasc(i,ict),ict,iadiag)=c(iwjasc(i,ict),ict,iadiag)-dparm(iparm)
          enddo
        enddo

      endif

      call cuspexact4(0,iadiag)

c Check parameters a2 and b2 > -scalek
      call check_parms_jas(iflag)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_lcao(dparm,iadiag)

      use vmc_mod, only: norb_tot, nwftypeorb, stoo
      use optwf_control, only: ioptorb, method
      use optwf_parms, only: nparmd, nparmj
      use coefs, only: nbasis
      use slater, only: norb, coef
      use optorb_cblock, only: norbterm
      use optwf_control, only: ioptorb
      use optwf_parms, only: nparmd,nparmj
      use orb_mat_022, only: ideriv
      use precision_kinds, only: dp
      use sr_mat_n, only: sr_state

      implicit none

      integer :: i, iadiag, io, j, jo, o!, k
      real(dp), dimension(nbasis, norb_tot) :: acoef
      real(dp), dimension(*) :: dparm

      if(ioptorb.eq.0) return

      if (method.eq.'sr_n'.and.nwftypeorb.gt.1) then
        o=stoo(sr_state)
          do i=1,norb
            do j=1,nbasis
              acoef(j,i)=coef(j,i,o)
            enddo
          enddo

c Update the orbitals
          do i=1,norbterm
            io=ideriv(1,i)
            jo=ideriv(2,i)
            do j=1,nbasis
              acoef(j,io)=acoef(j,io)-dparm(i+nparmj+nparmd)*coef(j,jo,o)
            enddo
          enddo

          do i=1,norb
            do j=1,nbasis
              coef(j,i,o)=acoef(j,i)
            enddo
          enddo

      else

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
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_ci(dparm,iadiag)

      use csfs,    only: ccsf,cxdet,iadet,ibdet,icxdet,ncsf,nstates
      use optwf_control, only: ioptci,ioptjas,ioptorb,method
      use optwf_parms, only: nparmj
      use precision_kinds, only: dp
      use slater,  only: cdet,ndet
      use sr_mat_n, only: sr_state

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
              cdet(j,k,1)=0.0d0
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
            cdet(idet,sr_state,iadiag)=cdet(idet,sr_state,iadiag)-dparm(idet-1+nparmj)
          enddo
        else
          do icsf=2,ncsf
            do j=iadet(icsf),ibdet(icsf)
              jx=icxdet(j)
              cdet(jx,sr_state,iadiag)=cdet(jx,sr_state,iadiag)-dparm(icsf-1+nparmj)*cxdet(j)
            enddo
            ccsf(icsf,sr_state,iadiag)=ccsf(icsf,sr_state,iadiag)-dparm(icsf-1+nparmj)
          enddo
        endif
      endif

c     do 90 j=1,nstates
c90     write(ounit,'(''csf ='',1000f20.15)') (ccsf(i,j,iadiag),i=1,ncsf)

      return
      end
c-----------------------------------------------------------------------
      subroutine check_parms_jas(iflag)

      use system, only: nctype
      use jastrow, only: b, scalek, a4
      use bparm, only: nspin2b
      use optwf_nparmj, only: nparma, nparmb
      use optwf_wjas, only: iwjasa, iwjasb
      use optwf_control, only: method
      use precision_kinds, only: dp
      use contrl_file,    only: ounit
      use vmc_mod, only: nwftypejas

      implicit none

      integer :: i, isp, ict, iflag, iflaga, iflagb, k
      real(dp) :: scalem

      iflag=0
      iflaga=0
      iflagb=0

      scalem=-scalek(1)

      if (method.eq.'sr_n'.and.nwftypejas.gt.1) then

        do k=1,nwftypejas
          do ict=1,nctype
            do i=1,nparma(ict)
              if(iwjasa(i,ict).eq.2.and.a4(2,ict,k).le.scalem) iflaga=1
            enddo
          enddo

          if(iflaga.eq.1) then
            do ict=1,nctype
              write(ounit,*) "Jastrow type: ", k
              write(ounit,'(''a2 < -scalek'',f10.5)') a4(2,ict,k)
            enddo
          endif
          do isp=1,nspin2b
            do i=1,nparmb(1)
              if(iwjasb(i,isp).eq.2.and.b(2,isp,k).le.scalem) iflagb=1
            enddo
          enddo
          if(iflagb.eq.1) then
            write(ounit,*) "Jastrow type: ", k
            write(ounit,'(''b2 < -scalek'',f10.5)') b(2,1,k)
          endif
        enddo

      else

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
        do isp=1,nspin2b
          do i=1,nparmb(1)
            if(iwjasb(i,isp).eq.2.and.b(2,isp,1).le.scalem) iflagb=1
          enddo
        enddo
        if(iflagb.eq.1) then
          write(ounit,'(''b2 < -scalek'',f10.5)') b(2,1,1)
        endif

      endif

      if(iflaga.eq.1.or.iflagb.eq.1) iflag=1
      return
      end

c-----------------------------------------------------------------------
      subroutine test_solution_parm(nparm,dparm,
     &              dparm_norm,dparm_norm_min,add_diag,iflag)
      use contrl_file, only: ounit
      use precision_kinds, only: dp

      implicit none

      integer :: i, iflag, nparm
      real(dp) :: add_diag, dparm_norm, dparm_norm_min
      real(dp), dimension(*) :: dparm

      iflag=0
      if(add_diag.le.0.d0) return

c Calculate rms change in parameters
      dparm_norm=0.0d0
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

      use optwf_control, only: ioptci, ioptjas, ioptorb
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

      use optwf_control, only: ioptci, ioptjas, ioptorb, nparm
      use optwf_parms, only: nparmd, nparmj
      use optorb_cblock, only: norbterm
      use ci000, only: nciterm
      use optwf_control, only: method
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
      subroutine optwf_store(l,wt,wt_sqrt,psid,energy)
c store elocal and derivatives of psi for each configuration (call in vmc)

      use sr_mod, only: mparm, mconf
      use optwf_parms, only: nparmj
      use csfs, only: nstates
      use derivjas, only: gvalue
      use optwf_control, only: ioptci, ioptjas, ioptorb
      use optwf_func, only: ifunc_omega
      use optwf_parms, only: nparmj
      use sr_mat_n, only: elocal, nconf_n, sr_ho, ortho
      use sr_mat_n, only: sr_o, wtg
      use deloc_dj_m, only: denergy
      use m_force_analytic, only: iforce_analy
      use optorb_cblock, only: norbterm
      use orb_mat_001, only: orb_ho, orb_o
      use ci000, only: nciterm
      use ci001_blk, only: ci_o
      use ci003_blk, only: ci_e
      use optwf_control, only: method
      use sr_mod, only: izvzb, i_sr_rescale
      use precision_kinds, only: dp
      use optgeo_lib, only: force_store
      use vmc_mod, only: nwftypeorb, nwftypejas, stoj, stoo
      use contrl_file, only: ounit

      implicit none

      integer :: i0, ii, ijasci, istate, j
      integer :: l, ntmp, k
      real(dp), dimension(nparmj) :: tmp_ho
      real(dp), dimension(*) :: wt
      real(dp), dimension(*) :: wt_sqrt
      real(dp), dimension(*) :: psid
      real(dp), dimension(*) :: energy

      if(iforce_analy.gt.0.and.izvzb.eq.1) call force_store(l)

      if((method.ne.'sr_n').and.(method.ne.'lin_d').or.(ioptjas+ioptorb+ioptci.eq.0))return

      i0=1
      if(method.eq.'lin_d'.and.ioptjas+ioptorb.eq.0) i0=0

      if(l.gt.mconf) call fatal_error('SR_STORE: l gt mconf')

      if (method.eq.'sr_n'.and.ortho.eq.1.or.nstates.eq.1) then ! for sr_n w/ ortho, or sr_n 1-state
        do istate=1,nstates
          if(nparmj /= 0) call dcopy(nparmj,gvalue(1,stoj(istate)),1,sr_o(1,l,istate),1)
          ntmp=max(nciterm-i0,0)
          if (ntmp /= 0) call dcopy(ntmp,ci_o(1+i0,istate),1,sr_o(nparmj+1,l,istate),1)
          ijasci=nparmj+ntmp
          if((ijasci+norbterm)*nstates.gt.mparm) call fatal_error('SR_STORE: iparm gt mparm')
          ii=ijasci
          if (norbterm /= 0) call dcopy(norbterm,orb_o(1,istate),1,sr_o(ii+1,l,istate),1)
          elocal(l,istate)=energy(istate)
          wtg(l,istate)=wt(istate)

          ii=ijasci+norbterm
          sr_o(ii+1,l,istate)=psid(istate)
          sr_o(ii+2,l,istate)=wt_sqrt(istate)
        enddo
      else !for lin_d or sr_n (in mix_n)
        if(nparmj /= 0) call dcopy(nparmj,gvalue(1,1),1,sr_o(1,l,1),1)
        ntmp=max(nciterm-i0,0)
        if (ntmp /= 0) call dcopy(ntmp,ci_o(1+i0,1),1,sr_o(nparmj+1,l,1),1)
        ijasci=nparmj+ntmp
        if(ijasci+nstates*norbterm+nstates.gt.mparm) call fatal_error('SR_STORE: iparm gt mparm')
        do istate=1,nstates
          ii=ijasci+(istate-1)*norbterm
          if (norbterm /= 0) call dcopy(norbterm,orb_o(1,istate),1,sr_o(ii+1,l,1),1)
          elocal(l,istate)=energy(istate)
          wtg(l,istate)=wt(istate)
          ii=ijasci+nstates*norbterm
          sr_o(ii+istate,l,1)=psid(istate)
        enddo
      endif

      nconf_n=l

      if(method.eq.'sr_n'.and.i_sr_rescale.eq.0.and.izvzb.eq.0.and.ifunc_omega.eq.0) return

c TO FIX: we are assuming optjas.ne.0 or optorb.ne.0 -> Otherwise, standard secular problem
      do j=1,nparmj
        tmp_ho(j)=denergy(j,1)+gvalue(j,1)*energy(1)
      enddo

      if(nparmj /= 0) call dcopy(nparmj,tmp_ho,1,sr_ho(1,l),1)

      if(ntmp /= 0) call dcopy(ntmp,ci_e(1+i0,1),1,sr_ho(nparmj+1,l),1)

      if(norbterm /= 0) call dcopy(norbterm,orb_ho(1,1),1,sr_ho(nparmj+ntmp+1,l),1)

      return
      end
c-----------------------------------------------------------------------
      end module
