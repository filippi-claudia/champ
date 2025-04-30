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
!-----------------------------------------------------------------------
      subroutine write_wf(iwf_fit,iter)

      use mpi
      use mpiconf, only: idtask
#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
      use trexio_read_data, only: update_trexio_orbitals
#endif

      implicit none

      integer :: index, iter, iwf_fit
      character(len=40) filetype,wf,itn

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
      call update_trexio_orbitals
#endif

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

      call write_lcao(iwf_fit,filetype)
      call write_jastrow(iwf_fit,filetype)
      call write_ci(iwf_fit,filetype)

      return
      end
!-----------------------------------------------------------------------
      subroutine write_wf_best
      implicit none

      call restore_jastrow_best
      call restore_lcao_best
      call restore_ci_best
      call write_wf(1,-1)

      return
      end
!-----------------------------------------------------------------------
      subroutine write_jastrow(iwf_fit,filetype)

      use system, only: nctype
      use jastrow, only: nspin1, nspin2, b, c, scalek, a4, norda, nordb, nordc, ijas, isc
      use bparm, only: nspin2b
      use optwf_control, only: ioptjas
      use optwf_nparmj, only: nparma, nparmb, nparmc
      use vmc_mod, only: nwftypejas, nstoj, extraj, jtos
      use precision_kinds, only: dp
      use system,  only: nctype

      implicit none

      integer :: i, isp, ict, index, iwf_fit, mparmja, k
      integer :: mparmjb, mparmjc
      character(len=50) fmt, temp
      character(len=40) filename,filetype

      if(ioptjas.eq.0) return

      filename='jastrow'//filetype(1:index(filetype,' ')-1)
      open(2,file=filename,status='unknown')
      write(2,'(''jastrow_parameter'',i4)') iwf_fit
      write(2,'(3i3,a28)') norda,nordb,nordc,' norda,nordb,nordc'

      if(ijas.eq.4) then

!      tmp
         write(2,'(f13.8,a15)') scalek(1),' scalek'
         do k=1,nwftypejas
!     if (extraj.eq.1) write(2,'(''jastrows_to_states'',i6,<nstoj(k)>i4)')
!     &                          nstoj(k), (jtos(k,i),i=1,nstoj(k))             ! Intel version
            if (extraj.eq.1) write(temp,'(a,i0,a,a)') '(a,1x,i0,1x,',nstoj(k),'(1x,i0,1x)',')'
            if (extraj.eq.1) write(2,temp) "jastrows_to_states",nstoj(k),(jtos(k,i),i=1,nstoj(k)) ! GNU version

            mparmja=2+max(0,norda-1)
            mparmjb=2+max(0,nordb-1)
            mparmjc=nterms4(nordc)
            if(mparmja.gt.0) then
               write(fmt,'(''('',i2,''f13.8,a28)'')') mparmja
            else
               write(fmt,'(''(a28)'')')
            endif
            do ict=1,nctype
               write(2,fmt) (a4(i,ict,k),i=1,mparmja),' (a(iparmj),iparmj=1,nparma)'
            enddo

            if(mparmjb.gt.0) then
               write(fmt,'(''('',i2,''f13.8,a28)'')') mparmjb
            else
               write(fmt,'(''(a28)'')')
            endif
            do isp=1,nspin2b
               write(2,fmt) (b(i,isp,k),i=1,mparmjb),' (b(iparmj),iparmj=1,nparmb)'
            enddo

            if(mparmjc.gt.0) then
               write(fmt,'(''('',i2,''f13.8,a28)'')') mparmjc
            else
               write(fmt,'(''(a28)'')')
            endif
            do ict=1,nctype
               write(2,fmt) (c(i,ict,k),i=1,mparmjc),' (c(iparmj),iparmj=1,nparmc)'
            enddo
         enddo

      else

!     tmp

         do k=1,nwftypejas
!     if (extraj.eq.1) write(2,'(''jastrows_to_states'',i6,<nstoj(k)>i4)')
!     &                          nstoj(k), (jtos(k,i),i=1,nstoj(k))             ! Intel version
            if (extraj.eq.1) write(temp,'(a,i0,a,a)') '(a,1x,i0,1x,',nstoj(k),'(1x,i0,1x)',')'
            if (extraj.eq.1) write(2,temp) "jastrows_to_states",nstoj(k),(jtos(k,i),i=1,nstoj(k)) ! GNU version

            mparmja=norda
            mparmjb=nordb
            mparmjc=nterms4(nordc)
            if(mparmja.gt.0) then
               write(fmt,'(''('',i2,''f13.8,a28)'')') mparmja+1
               do ict=1,nctype
                  write(2,fmt) (a4(i,ict,k),i=1,mparmja+1),' (a(iparmj),iparmj=1,nparma)'
               enddo
            else
               write(fmt,'(''(a28)'')')
               do ict=1,nctype
                  write(2,fmt) ' (a(iparmj),iparmj=1,nparma)'
               enddo
            endif

            if(mparmjb.gt.0) then
               write(fmt,'(''('',i2,''f13.8,a28)'')') mparmjb+1
               do isp=1,nspin2b
                  write(2,fmt) (b(i,isp,k),i=1,mparmjb+1),' (b(iparmj),iparmj=1,nparmb)'
               enddo
            else
               write(fmt,'(''(a28)'')')
               do isp=1,nspin2b
                  write(2,fmt) ' (b(iparmj),iparmj=1,nparmb)'
               enddo
            endif

            if(mparmjc.gt.0) then
               write(fmt,'(''('',i2,''f13.8,a28)'')') mparmjc
            else
               write(fmt,'(''(a28)'')')
            endif
            do ict=1,nctype
               write(2,fmt) (c(i,ict,k),i=1,mparmjc),' (c(iparmj),iparmj=1,nparmc)'
            enddo
         enddo

      endif

      write(2,'(''end'')')
      close(2)
      return
      end
!-----------------------------------------------------------------------
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
      character(len=40) filename,filetype, temp

      if(ioptorb.eq.0) return

      filename='orbitals'//filetype(1:index(filetype,' ')-1)
      open(2,file=filename,status='unknown')
      write(2,'(''lcao '',3i4)') norb+nadorb,nbasis,iwf_fit

      do k=1,nwftypeorb
!        if (extrao.eq.1) write(2,'(''orbitals_to_states'',i6,<nstoo(k)>i4)')
!     &                         nstoo(k), (otos(k,i),i=1,nstoo(k))           ! Intel version

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
!-----------------------------------------------------------------------
      subroutine write_ci(iwf_fit,filetype)

      use csfs,    only: ccsf,cxdet,iadet,ibdet,icxdet,ncsf,nstates
      use dorb_m,  only: iworbd
      use optwf_control, only: ioptci
      use slater,  only: cdet,kref,ndet
      use system,  only: nelec,nup

      implicit none

      integer :: i, index, istate, iwf_fit, j
      integer :: k, nmap, nptr, nterm
      character(len=40) filename,filetype

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
!
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
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
      subroutine save_wf

      use optwf_control, only: ioptci,ioptjas,ioptorb

      implicit none

      if(ioptjas.ne.0) call save_jastrow
      if(ioptorb.ne.0) call save_lcao
      if(ioptci.ne.0) call save_ci

      return
      end
!-----------------------------------------------------------------------
      subroutine restore_wf(iadiag)

      use optwf_control, only: ioptci,ioptjas,ioptorb

      implicit none

      integer :: iadiag

      if(ioptjas.ne.0) call restore_jastrow(iadiag)
      if(ioptorb.ne.0) call restore_lcao(iadiag)
      if(ioptci.ne.0) call restore_ci(iadiag)

      return
      end
!-----------------------------------------------------------------------
      subroutine save_wf_best(ioptjas,ioptorb,ioptci)

      implicit none

      integer :: ioptci, ioptjas, ioptorb

      if(ioptjas.ne.0) call save_jastrow_best
      if(ioptorb.ne.0) call save_lcao_best
      if(ioptci.ne.0) call save_ci_best

      return
      end
!-----------------------------------------------------------------------
      subroutine save_jastrow

      use bparm,   only: nspin2b
      use jastrow, only: ijas, norda,nordb,nordc
      use jastrow, only: a4,b,c,nordj1
      use save_mod,only: mparmja, mparmjb, mparmjc
      use save_mod,only: a4_save, b_save, c_save
      use multiple_geo, only: MWF, nwftype
      use precision_kinds, only: dp
      use system, only: nctype, nctype_tot
      use jastrow, only: b, c, a4, norda, nordb, nordc, nordj1
      use bparm, only: nspin2b
      use vmc_mod, only: nwftypejas
      use optwf_control, only: method

      implicit none

      integer :: i, isp, iadiag, ict, k

! Save parameters corresponding to run generating hessian

      if(ijas.eq.1) then
         mparmja=norda
         mparmjb=nordb
      else
         mparmja=2+max(0,norda-1)
         mparmjb=2+max(0,nordb-1)
      endif

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
      end subroutine
!-----------------------------------------------------------------------
      subroutine restore_jastrow(iadiag)
      use bparm, only: nspin2b
      use jastrow, only: b, c, a4, norda, nordb, nordc, nordj1
      use scale_dist_mod, only: set_scale_dist
      use save_mod,only: mparmja, mparmjb, mparmjc
      use save_mod,only: a4_save, b_save, c_save
      use system, only: nctype

      implicit none

      integer :: i, isp, iadiag, ict

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

      call set_scale_dist(iadiag,0)

      return
      end
!-----------------------------------------------------------------------
      subroutine save_lcao

      use coefs,   only: nbasis
      use multiple_geo, only: nwftype
      use optwf_control, only: method
      use precision_kinds, only: dp
      use vmc_mod, only: norb_tot, nwftypeorb
      use save_mod, only: coef_save
      use slater, only: norb, coef

      implicit none

      integer :: i, iadiag, j, k

      if (method.eq.'sr_n'.and.nwftypeorb.gt.1) then
        if (.not. allocated(coef_save)) allocate(coef_save(nbasis, norb_tot, nwftypeorb))

        do k=1,nwftypeorb
          do i=1,norb
            do j=1,nbasis
              coef_save(j,i,k)=coef(j,i,k)
            enddo
          enddo
        enddo

      else
        if (.not. allocated(coef_save)) allocate(coef_save(nbasis, norb_tot, nwftype))

        do i=1,norb
          do j=1,nbasis
            coef_save(j,i,1)=coef(j,i,1)
          enddo
        enddo

      endif

      end subroutine
!-----------------------------------------------------------------------
      subroutine restore_lcao(iadiag)
      use coefs, only: nbasis
      use save_mod, only: coef_save
      use slater, only: norb, coef
      implicit none
      integer :: i, iadiag, j

      do i=1,norb
        do j=1,nbasis
          coef(j,i,iadiag)=coef_save(j,i,1)
        enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine save_ci

      use csfs,    only: ccsf,cxdet,iadet,ibdet,icxdet,ncsf,nstates
      use mstates_mod, only: MSTATES
      use precision_kinds, only: dp
      use save_mod, only: cdet_save, ccsf_save
      use set_input_data, only: multideterminants_define
      use slater,  only: cdet,ndet

      implicit none

      integer :: i, iadiag, icsf, j, k

      if(.not. allocated(cdet_save)) allocate(cdet_save(ndet,MSTATES))
      if(.not. allocated(ccsf_save)) allocate(ccsf_save(ndet,MSTATES))

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

      end subroutine
!-----------------------------------------------------------------------
      subroutine restore_ci(iadiag)
      use csfs,    only: ccsf,cxdet,iadet,ibdet,icxdet,ncsf,nstates
      use mstates_mod, only: MSTATES
      use save_mod, only: cdet_save, ccsf_save
      use set_input_data, only: multideterminants_define
      use slater,  only: cdet,ndet
      implicit none
      integer :: i, iadiag, icsf, j, k
      integer :: kx

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

! if kref (iwdorb, cxdet) has changed
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

      endif

! reset kref=1
      call multideterminants_define(0)

      return
      end
!-----------------------------------------------------------------------
      subroutine copy_jastrow(iadiag)

      use bparm, only: nspin2b
      use jastrow, only: b, c, scalek, a4, ijas, norda, nordb, nordc
      use scale_dist_mod, only: set_scale_dist
      use system, only: nctype

      implicit none

      integer :: i, isp, iadiag, ict, mparmja, mparmjb, mparmjc

      if(ijas.eq.1) then
        mparmja=norda
        mparmjb=nordb

        do ict=1,nctype
          a4(mparmja+1,ict,iadiag)=a4(mparmja+1,ict,1)
        enddo
        do isp=1,nspin2b
          b(mparmjb+1,isp,iadiag)=b(mparmjb+1,isp,1)
        enddo
      else
        mparmja=2+max(0,norda-1)
        mparmjb=2+max(0,nordb-1)
        scalek(iadiag)=scalek(1)
      endif

      mparmjc=nterms4(nordc)

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

      call set_scale_dist(iadiag,0)

      return
      end

!-----------------------------------------------------------------------
      subroutine copy_lcao(iadiag)

      use coefs, only: nbasis
      use slater, only: norb, coef
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
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
      subroutine copy_bas_num(iadiag)

      use system,  only: nctype,newghostype
      use numbas,  only: rwf, d2rwf,nr,nrbas
      use numexp,  only: ae,ce

      implicit none

      integer :: i, iadiag, ic, irb

      do ic=1,nctype+newghostype
        do irb=1,nrbas(ic)
          rwf(1:nr(ic),irb,ic,iadiag)=rwf(1:nr(ic),irb,ic,1)
          d2rwf(1:nr(ic),irb,ic,iadiag)=d2rwf(1:nr(ic),irb,ic,1)
          ce(1:5,irb,ic,iadiag)=ce(1:5,irb,ic,1)
          ae(1:2,irb,ic,iadiag)=ae(1:2,irb,ic,1)
        enddo
      enddo
          
      return
      end
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
      subroutine save_jastrow_best

      use bparm,   only: nspin2b
      use jastrow, only: ijas, norda,nordb,nordc
      use jastrow, only: a4,b,c,nordj1
      use multiple_geo, only: nwftype
      use precision_kinds, only: dp
      use save_mod, only: a4_best, b_best, c_best
      use save_mod, only:  mparmja_best, mparmjb_best, mparmjc_best
      use system,  only: nctype,nctype_tot

      implicit none

      integer :: i, isp, ict

      if(.not.allocated(a4_best)) allocate(a4_best(nordj1,nctype_tot,nwftype))
      if(.not.allocated(b_best))  allocate(b_best(nordj1,2,nwftype))
      if(.not.allocated(c_best))  allocate(c_best(83,nctype_tot,nwftype))

! Save parameters corresponding to run generating hessian

      if(ijas.eq.1) then
         mparmja_best=norda
         mparmjb_best=nordb
      else
         mparmja_best=2+max(0,norda-1)
         mparmjb_best=2+max(0,nordb-1)
      endif

      mparmjc_best=nterms4(nordc)

      do ict=1,nctype
        do i=1,mparmja_best
          a4_best(i,ict,1)=a4(i,ict,1)
        enddo
      enddo
      do isp=1,nspin2b
        do i=1,mparmjb_best
          b_best(i,isp,1)=b(i,isp,1)
        enddo
      enddo
      do ict=1,nctype
        do i=1,mparmjc_best
          c_best(i,ict,1)=c(i,ict,1)
        enddo
      enddo

      end subroutine
!-----------------------------------------------------------------------
      subroutine restore_jastrow_best

      use bparm,   only: nspin2b
      use jastrow, only: a4,b,c,nordj1
      use multiple_geo, only: nwftype
      use scale_dist_mod, only: set_scale_dist
      use save_mod, only: a4_best, b_best, c_best
      use save_mod, only:  mparmja_best, mparmjb_best, mparmjc_best
      use system,  only: nctype,nctype_tot

      implicit none

      integer :: i, isp, ict

      if(.not.allocated(a4_best)) allocate(a4_best(nordj1,nctype_tot,nwftype))
      if(.not.allocated(b_best))  allocate(b_best(nordj1,2,nwftype))
      if(.not.allocated(c_best))  allocate(c_best(83,nctype_tot,nwftype))

! Restore parameters corresponding to run generating hessian
      do ict=1,nctype
        do i=1,mparmja_best
          a4(i,ict,1)=a4_best(i,ict,1)
        enddo
      enddo
      do isp=1,nspin2b
        do i=1,mparmjb_best
          b(i,isp,1)=b_best(i,isp,1)
        enddo
      enddo
      do ict=1,nctype
        do i=1,mparmjc_best
          c(i,ict,1)=c_best(i,ict,1)
        enddo
      enddo

      call set_scale_dist(1,0)

      return
      end
!-----------------------------------------------------------------------
      subroutine save_lcao_best

      use coefs,   only: nbasis
      use multiple_geo, only: nwftype
      use optwf_control, only: method
      use precision_kinds, only: dp
      use save_mod, only: coef_best
      use slater, only: norb, coef
      use vmc_mod, only: norb_tot, nwftypeorb

      implicit none

      integer :: i, j, k

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

      end subroutine
!-----------------------------------------------------------------------
      subroutine restore_lcao_best

      use coefs, only: nbasis
      use multiple_geo, only: nwftype
      use optwf_control, only: method
      use save_mod, only: coef_best
      use slater, only: norb, coef
      use vmc_mod, only: norb_tot, nwftypeorb
      implicit none

      integer :: i, j, k

!     if(ioptorb.eq.0) return

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
!-----------------------------------------------------------------------
      subroutine save_ci_best

      use csfs,    only: ccsf,cxdet,iadet,ibdet,icxdet,ncsf,nstates
      use mstates_mod, only: MSTATES
      use precision_kinds, only: dp
      use save_mod, only: cdet_best, ccsf_best
      use set_input_data, only: multideterminants_define
      use slater,  only: cdet,ndet

      implicit none

      integer :: i, icsf, j

      if(.not. allocated(cdet_best)) allocate(cdet_best(ndet,nstates))
      if(.not. allocated(ccsf_best)) allocate(ccsf_best(ndet,nstates))

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

      end subroutine
!-----------------------------------------------------------------------
      subroutine restore_ci_best
      use csfs,    only: ccsf,cxdet,iadet,ibdet,icxdet,ncsf,nstates
      use mstates_mod, only: MSTATES
      use save_mod, only: cdet_best, ccsf_best
      use set_input_data, only: multideterminants_define
      use slater,  only: cdet,ndet
      implicit none

      integer :: i, icsf, j, k, kx

      if(.not. allocated(cdet_best)) allocate(cdet_best(ndet,nstates))
      if(.not. allocated(ccsf_best)) allocate(ccsf_best(ndet,nstates))

!     if(ioptci.eq.0) return

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

! if kref (iwdetorb, cxdet) has changed
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

! reset kref=1
      call multideterminants_define(0)
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine compute_parameters(dparm,iflag,iadiag)

      use precision_kinds, only: dp
      use contrl_file, only: ounit

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
!-----------------------------------------------------------------------
      subroutine compute_jastrow(dparm,iflag,iadiag)

      use control, only: ipr
      use cuspinit4_mod, only: cuspinit4
      use cuspexact4_mod, only: cuspexact4
      use jastrow, only: b, c, a4
      use bparm, only: nspin2b
      use optwf_control, only: ioptjas, method
      use optwf_nparmj, only: nparma, nparmb, nparmc
      use optwf_wjas, only: iwjasa, iwjasb, iwjasc
      use precision_kinds, only: dp
      use scale_dist_mod, only: set_scale_dist
      use sr_mat_n, only: sr_state
      use system, only: nctype
      use vmc_mod, only: nwftypejas, stoj
      use contrl_file, only: ounit

      implicit none

      integer :: i, iadiag, isp, ict, iflag, iparm, j
      real(dp), dimension(*) :: dparm

      if(ioptjas.eq.0) return

! Set up cusp conditions
      call cuspinit4(0)

      j=iadiag
      if (method.eq.'sr_n'.and.nwftypejas.gt.1) j=stoj(sr_state)

! Add change to old parameters
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

      call cuspexact4(0,j)

! Check parameters a2 and b2 > -scalek
      call check_parms_jas(j,iflag)

      call set_scale_dist(j,ipr)

      return
      end
!-----------------------------------------------------------------------
      subroutine compute_lcao(dparm,iadiag)

      use vmc_mod, only: norb_tot, nwftypeorb, stoo
      use optwf_control, only: ioptorb, method, orbitals_ortho
      use optwf_parms, only: nparmd, nparmj
      use orbval, only: nadorb
      use coefs, only: nbasis
      use slater, only: norb, coef
      use optorb_cblock, only: norbterm
      use optwf_control, only: ioptorb
      use optwf_parms, only: nparmd,nparmj
      use orb_mat_022, only: ideriv
      use precision_kinds, only: dp
      use sr_mat_n, only: sr_state
      use contrl_file, only: ounit

      implicit none

      integer :: i, iadiag, io, j, jo, k, o
      real(dp), dimension(nbasis, norb_tot) :: acoef
      real(dp), dimension(*) :: dparm
      real(dp), dimension(norb+nadorb,norb+nadorb) :: xmat, umat

      if(ioptorb.eq.0) return

      if (method.eq.'sr_n'.and.nwftypeorb.gt.1) then
        o=stoo(sr_state)
       else
        o=iadiag
      endif

      if(.not.orbitals_ortho) then

        do i=1,norb
          do j=1,nbasis
            acoef(j,i)=coef(j,i,o)
          enddo
        enddo

! Update the orbitals
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

        xmat(:,:)=0.d0
        do i=1,norbterm
          io=ideriv(1,i)
          jo=ideriv(2,i)
          xmat(io,jo)=-dparm(i+nparmj+nparmd)
        enddo

        do i=1,norb+nadorb
          do j=1,i-1
            if(xmat(i,j).ne.0.d0.and.xmat(j,i).ne.0.d0) then
              xmat(i,j)=0.5*(xmat(i,j)-xmat(j,i))
              xmat(j,i)=-xmat(i,j)
             elseif(xmat(i,j).ne.0.d0) then
              xmat(j,i)=-xmat(i,j)
             else
              xmat(i,j)=-xmat(j,i)
            endif
          enddo
        enddo

        call ortho_orbitals(norb+nadorb,xmat,umat)

        acoef(:,:)=0.d0
        do i=1,norb+nadorb
          do j=1,norb+nadorb
            do k=1,nbasis
              acoef(k,i)=acoef(k,i)+umat(i,j)*coef(k,j,o)
            enddo
          enddo
        enddo

        do i=1,norb+nadorb
          if(acoef(1,i)*coef(1,i,o).lt.0.d0) call fatal_error('COMPUTE_LCAO: orbitals have changed sign through ortho')
        enddo

        do i=1,norb+nadorb
          do j=1,nbasis
            coef(j,i,o)=acoef(j,i)
          enddo
        enddo

      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine compute_ci(dparm,iadiag)

      use csfs,    only: ccsf,cxdet,iadet,ibdet,icxdet,maxcsf,ncsf,nstates
      use optwf_control, only: ioptci,ioptjas,ioptorb,method
      use optwf_parms, only: nparmj
      use precision_kinds, only: dp
      use slater,  only: cdet,ndet
      use sr_mat_n, only: sr_state

      implicit none

      integer :: i, iadiag, icsf, idet, ish, j
      integer :: jx, k
      real(dp) :: c90
      real(dp), dimension(*) :: dparm

      if(ioptci.eq.0) return

! Update the ci coef
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
          do icsf=1,maxcsf(sr_state)-1
            do j=iadet(icsf),ibdet(icsf)
              jx=icxdet(j)
              cdet(jx,sr_state,iadiag)=cdet(jx,sr_state,iadiag)-dparm(icsf+nparmj)*cxdet(j)
            enddo
            ccsf(icsf,sr_state,iadiag)=ccsf(icsf,sr_state,iadiag)-dparm(icsf+nparmj)
          enddo

          ish=1
          do icsf=maxcsf(sr_state)+1,ncsf
            do j=iadet(icsf),ibdet(icsf)
              jx=icxdet(j)
              cdet(jx,sr_state,iadiag)=cdet(jx,sr_state,iadiag)-dparm(icsf-ish+nparmj)*cxdet(j)
            enddo
            ccsf(icsf,sr_state,iadiag)=ccsf(icsf,sr_state,iadiag)-dparm(icsf-ish+nparmj)
          enddo

        endif
      endif

!     do 90 j=1,nstates
!90     write(ounit,'(''csf ='',1000f20.15)') (ccsf(i,j,iadiag),i=1,ncsf)

      return
      end
!-----------------------------------------------------------------------
      subroutine check_parms_jas(k,iflag)

      use system, only: nctype
      use jastrow, only: ijas, b, scalek, a4
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

      if(ijas.eq.1) return
      iflag=0
      iflaga=0
      iflagb=0

      scalem=-scalek(1)

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

      if(iflaga.eq.1.or.iflagb.eq.1) iflag=1

      return
      end

!-----------------------------------------------------------------------
      subroutine test_solution_parm(nparm,dparm, &
                    dparm_norm,dparm_norm_min,add_diag,iflag)
      use contrl_file, only: ounit
      use precision_kinds, only: dp

      implicit none

      integer :: i, iflag, nparm
      real(dp) :: add_diag, dparm_norm, dparm_norm_min
      real(dp), dimension(*) :: dparm

      iflag=0
      if(add_diag.le.0.d0) return

! Calculate rms change in parameters
      dparm_norm=0.0d0
      do i=1,nparm
        dparm_norm=dparm_norm+dparm(i)**2
      enddo
      dparm_norm=sqrt(dparm_norm/nparm)

      write(ounit,'(''dparm_norm,adiag ='',3g12.5)') &
      dparm_norm,add_diag

      if(dparm_norm.gt.dparm_norm_min) iflag=1

      return
      end
!-----------------------------------------------------------------------
      subroutine save_nparms

      use optwf_control, only: ioptci, ioptjas, ioptorb
      use optwf_parms, only: nparmd, nparmj
      use optorb_cblock, only: norbterm, nreduced
      use ci000, only: nciterm
      use contrl_file,    only: ounit
      use orbval, only: nadorb
      use save_mod, only: nciterm_sav, norbterm_sav, nparmd_sav, nparmj_sav, nreduced_sav, nadorb_sav

      implicit none

      nparmj_sav=nparmj
      norbterm_sav=norbterm
      nreduced_sav=nreduced
      nciterm_sav=nciterm
      nparmd=max(nciterm-1,0)
      nparmd_sav=nparmd
      nadorb_sav=nadorb

      write(ounit,'(''Saved max number of parameters, nparmj,norb,nciterm,nciterm-1: '',5i5)') nparmj,norbterm,nciterm,nparmd
      end subroutine
!-----------------------------------------------------------------------
      subroutine set_nparms

      use ci000, only: nciterm
      use contrl_file,    only: ounit
      use optorb_cblock, only: norbterm, nreduced
      use optwf_control, only: ioptci, ioptjas, ioptorb
      use optwf_parms, only: nparmd, nparmj
      use orbval, only: nadorb
      use save_mod, only: nciterm_sav, norbterm_sav, nparmd_sav, nparmj_sav, nreduced_sav, nadorb_sav
      implicit none

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
!-----------------------------------------------------------------------
      subroutine set_nparms_tot

      use optwf_control, only: ioptci, ioptjas, ioptorb, nparm
      use optwf_parms, only: nparmd, nparmj
      use optorb_cblock, only: norbterm
      use ci000, only: nciterm
      use optwf_control, only: method
      use contrl_file,    only: ounit
      implicit none

      integer :: i0

! Note: we do not vary the first (i0) CI coefficient unless a run where we only optimize the CI coefs

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

      write(ounit,'(/,''number of parms: total, Jastrow, CI, orbitals= '',4i5)') &
       nparm,nparmj,nciterm,norbterm

      return
      end
!-----------------------------------------------------------------------
      subroutine optwf_store(l,wt,wt_sqrt,psid,energy)
! store elocal and derivatives of psi for each configuration (call in vmc)

      use sr_mod, only: mparm, mconf
      use optwf_parms, only: nparmj
      use csfs, only: maxcsf,nstates
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
      use control_vmc, only: vmc_nstep, vmc_nblk_max

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

      if(l.gt.mconf) then
        print*, "l",l, "mconf", mconf
        print*, "vmc_nstep", vmc_nstep, "vmc_nblk_max", vmc_nblk_max
        call fatal_error('SR_STORE: l gt mconf')
      endif

      if (method.eq.'sr_n'.and.(ortho.eq.1.or.nstates.eq.1)) then ! for sr_n w/ ortho, or sr_n 1-state
        do istate=1,nstates
          if(nparmj /= 0) call dcopy(nparmj,gvalue(1,stoj(istate)),1,sr_o(1,l,istate),1)
          ntmp=max(nciterm-i0,0)
!         if (ntmp /= 0) call dcopy(ntmp,ci_o(1+i0,istate),1,sr_o(nparmj+1,l,istate),1)
          if (ntmp /= 0) then
            if(maxcsf(istate).gt.1) call dcopy(maxcsf(istate)-1,ci_o(1,istate),1,sr_o(nparmj+1,l,istate),1)
            if(maxcsf(istate).lt.nciterm) &
            call dcopy(nciterm-maxcsf(istate),ci_o(maxcsf(istate)+1,istate),1,sr_o(nparmj+maxcsf(istate),l,istate),1)
          endif
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

! TO FIX: we are assuming optjas.ne.0 or optorb.ne.0 -> Otherwise, standard secular problem
      do j=1,nparmj
        tmp_ho(j)=denergy(j,1)+gvalue(j,1)*energy(1)
      enddo

      if(nparmj /= 0) call dcopy(nparmj,tmp_ho,1,sr_ho(1,l),1)

      if(ntmp /= 0) call dcopy(ntmp,ci_e(1+i0,1),1,sr_ho(nparmj+1,l),1)

      if(norbterm /= 0) call dcopy(norbterm,orb_ho(1,1),1,sr_ho(nparmj+ntmp+1,l),1)

      return
      end
!-----------------------------------------------------------------------
end module
