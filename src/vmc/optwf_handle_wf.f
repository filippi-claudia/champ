c-----------------------------------------------------------------------
      subroutine write_wf(iwf_fit,iter)
      use mpiconf, only: idtask
      implicit real*8(a-h,o-z)


      include 'mpif.h'
      include 'vmc.h'
      include 'force.h'


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
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

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
      use jaspar3, only: a, b, c, scalek

      use jaspar4, only: a4, norda, nordb, nordc
      use optwf_contrl, only: ioptjas
      use optwf_nparmj, only: nparma, nparmb, nparmc

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      character*50 fmt
      character*40 filename,filetype



      if(ioptjas.eq.0) return

      filename='jastrow'//filetype(1:index(filetype,' ')-1)

      open(2,file=filename,status='unknown')

      call p2gtid('jastrow:ianalyt_lap',ianalyt_lap,1,1)
      call p2gti('jastrow:ijas',ijas,1)
      call p2gti('jastrow:isc',isc,1)
      call p2gtid('jastrow:nspin1',nspin1,1,1)
      call p2gtid('jastrow:nspin2',nspin2,1,1)
      call p2gtid('jastrow:ifock',ifock,0,1)
      write(2,'(''&jastrow ianalyt_lap'',i2,'' ijas'',i2,'' isc'',i2,
     &'' nspin1'',i2,'' nspin2'',i2,'' ifock'',i2)') ianalyt_lap,ijas,isc,nspin1,nspin2,ifock
      write(2,*)
      write(2,'(''jastrow_parameter'',i4)') iwf_fit
      write(2,'(3i3,a28)') norda,nordb,nordc,' norda,nordb,nordc'
c tmp
      a21=0
      write(2,'(2f13.8,a15)') scalek(1),a21,' scalek,a21'
      mparmja=2+max(0,norda-1)
      mparmjb=2+max(0,nordb-1)
      mparmjc=nterms4(nordc)
      if(mparmja.gt.0) then
        write(fmt,'(''(''i2,''f13.8,a28)'')') mparmja
       else
        write(fmt,'(''(a28)'')')
      endif
      do 80 ict=1,nctype
   80   write(2,fmt) (a4(i,ict,1),i=1,mparmja),' (a(iparmj),iparmj=1,nparma)'

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
      do 90 ict=1,nctype
   90   write(2,fmt) (c(i,ict,1),i=1,mparmjc),' (c(iparmj),iparmj=1,nparmc)'
      write(2,'(''end'')')
      close(2)

      return
      end
c-----------------------------------------------------------------------
      subroutine write_lcao(iwf_fit,filetype)
      use numbas, only: numr

      use optwf_contrl, only: ioptorb
      use coefs, only: coef, nbasis, norb
      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'force.h'
      include 'numbas.h'


      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb



      character*40 filename,filetype

      dimension anorm(MBASIS)

      if(ioptorb.eq.0) return

      filename='orbitals'//filetype(1:index(filetype,' ')-1)
      open(2,file=filename,status='unknown')
      write(2,'(''lcao '',3i4)') norb+nadorb,nbasis,iwf_fit

      call p2gtfd('general:scalecoef',scalecoef,1.0d0,1)
      if(numr.gt.0) then
        do 20 i=1,norb+nadorb
   20     write(2,'(1000e20.8)') (coef(j,i,1)/scalecoef,j=1,nbasis)
      else
        call basis_norm(1,anorm,1)
        do 40 i=1,norb+nadorb
   40     write(2,'(1000e20.8)') (coef(j,i,1)/(anorm(j)*scalecoef),j=1,nbasis)
      endif

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
      implicit real*8(a-h,o-z)






      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      character*40 filename,filetype





      if(ioptci.eq.0) return

      filename='det'//filetype(1:index(filetype,' ')-1)
      open(2,file=filename,status='unknown')

      call p2gti('electrons:nelec',nelec,1)
      call p2gti('electrons:nup',nup,1)
      write(2,'(''&electrons nelec '',i4,'' nup '',i4)') nelec,nup
      write(2,'(''# kref'',i4)') kref
      do 1 istate=1,nstates
      write(2,'(''# State '',i4)') istate
      write(2,'(''determinants'',i10,i4)') ndet,iwf_fit
      write(2,'(100f15.8)') (cdet(i,istate,1),i=1,ndet)
      do 1 k=1,ndet
   1   write(2,'(100i4)') (iworbd(i,k),i=1,nelec)
 
      write(2,'(''end'')')

      if(ncsf.ne.0) then
        write(2,'(''csf '',i10,i4)') ncsf,nstates
        do i=1,nstates
          write(2,'(100f15.8)') (ccsf(j,i,1),j=1,ncsf)
        enddo
        write(2,'(''end'')')
c
        nmap=0
        do 5 i=1,ncsf
   5      nmap=nmap+ibdet(i)-iadet(i)+1
        write(2,'(''csfmap'')') 
        write(2,'(3i10)') ncsf,ndet,nmap
        nptr=0
        do 10 i=1,ncsf
         nterm=ibdet(i)-iadet(i)+1
         write(2,'(i10)') ibdet(i)-iadet(i)+1
         do 12 j=1,nterm
           nptr=nptr+1
           write(2,'(i10,f10.6)') icxdet(nptr),cxdet(nptr)
  12     enddo
  10    enddo
        write(2,'(''end'')')
      endif

      close(2)

      return
      end
c-----------------------------------------------------------------------
      subroutine setup_wf
      implicit real*8(a-h,o-z)
  
      do 10 k=2,3
        call copy_jastrow(k)
        call copy_lcao(k)
  10    call copy_ci(k)

      return
      end
c-----------------------------------------------------------------------
      subroutine save_wf
      use optwf_contrl, only: ioptci, ioptjas, ioptorb
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'


      if(ioptjas.ne.0) call save_jastrow
      if(ioptorb.ne.0) call save_lcao 
      if(ioptci.ne.0) call save_ci

      return
      end
c-----------------------------------------------------------------------
      subroutine restore_wf(iadiag)
      use optwf_contrl, only: ioptci, ioptjas, ioptorb
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'


      if(ioptjas.ne.0) call restore_jastrow(iadiag)
      if(ioptorb.ne.0) call restore_lcao(iadiag)
      if(ioptci.ne.0) call restore_ci(iadiag)

      return
      end
c-----------------------------------------------------------------------
      subroutine save_wf_best(ioptjas,ioptorb,ioptci)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      if(ioptjas.ne.0) call save_jastrow_best
      if(ioptorb.ne.0) call save_lcao_best
      if(ioptci.ne.0) call save_ci_best

      return
      end
c-----------------------------------------------------------------------
      subroutine save_jastrow

      use atom, only: nctype

      use jaspar3, only: a, b, c

      use jaspar4, only: a4, norda, nordb, nordc
      implicit real*8(a-h,o-z)




      include 'vmc.h'
      include 'force.h'


      dimension a4_save(MORDJ1,MCTYPE,MWF),b_save(MORDJ1,2,MWF),
     &c_save(83,MCTYPE,MWF)

      save a4_save,b_save,c_save
      save mparmja,mparmjb,mparmjc

c Save parameters corresponding to run generating hessian

      mparmja=2+max(0,norda-1)
      mparmjb=2+max(0,nordb-1)
      mparmjc=nterms4(nordc)

      do 50 ict=1,nctype
        do 50 i=1,mparmja
   50     a4_save(i,ict,1)=a4(i,ict,1)
      do 60 i=1,mparmjb
   60   b_save(i,1,1)=b(i,1,1)
      do 70 ict=1,nctype
        do 70 i=1,mparmjc
   70     c_save(i,ict,1)=c(i,ict,1)

      return

      entry restore_jastrow(iadiag)

c Restore parameters corresponding to run generating hessian
      do 80 ict=1,nctype
        do 80 i=1,mparmja
   80     a4(i,ict,iadiag)=a4_save(i,ict,1)
      do 90 i=1,mparmjb
   90   b(i,1,iadiag)=b_save(i,1,1)
      do 100 ict=1,nctype
        do 100 i=1,mparmjc
  100     c(i,ict,iadiag)=c_save(i,ict,1)

      return
      end

c-----------------------------------------------------------------------
      subroutine save_lcao
      use coefs, only: coef, nbasis, norb
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'


      dimension coef_save(MBASIS,MORB,MWF)

      save coef_save

      do 10 i=1,norb
       do 10 j=1,nbasis
   10   coef_save(j,i,1)=coef(j,i,1)

      return

      entry restore_lcao(iadiag)

      do 20 i=1,norb
       do 20 j=1,nbasis
   20   coef(j,i,iadiag)=coef_save(j,i,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine save_ci
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use dets, only: cdet, ndet
      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'



      dimension cdet_save(MDET,MSTATES),ccsf_save(MDET,MSTATES)
      save cdet_save,ccsf_save

      do 10 j=1,nstates
        do 10 i=1,ndet
   10     cdet_save(i,j)=cdet(i,j,1)

      do 20 j=1,nstates
       do 20 icsf=1,ncsf
   20   ccsf_save(icsf,j)=ccsf(icsf,j,1)

      return

      entry restore_ci(iadiag)

      do 30 j=1,nstates
        do 30 i=1,ndet
   30     cdet(i,j,iadiag)=cdet_save(i,j)

      do 40 j=1,nstates
       do 40 icsf=1,ncsf
   40   ccsf(icsf,j,iadiag)=ccsf_save(icsf,j)

c if kref (iwdetorb, cxdet) has changed
      if(ncsf.gt.0) then
        do 50 j=1,nstates
          do 45 k=1,ndet
   45       cdet(k,j,iadiag)=0
          do 50 icsf=1,ncsf
            do 50 k=iadet(icsf),ibdet(icsf)
              kx=icxdet(k)
              cdet(kx,j,iadiag)=cdet(kx,j,iadiag)+ccsf(icsf,j,iadiag)*cxdet(k)
   50  continue

c reset kref=1
      call multideterminants_define(0,0)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_jastrow(iadiag)
 
      use atom, only: nctype

      use jaspar3, only: a, b, c, scalek

      use jaspar4, only: a4, norda, nordb, nordc
      implicit real*8(a-h,o-z)




      include 'vmc.h'
      include 'force.h'


      mparmja=2+max(0,norda-1)
      mparmjb=2+max(0,nordb-1)
      mparmjc=nterms4(nordc)

      scalek(iadiag)=scalek(1)
      do 80 ict=1,nctype
        do 80 i=1,mparmja
   80     a4(i,ict,iadiag)=a4(i,ict,1)
      do 90 i=1,mparmjb
   90   b(i,1,iadiag)=b(i,1,1)
      do 100 ict=1,nctype
        do 100 i=1,mparmjc
  100     c(i,ict,iadiag)=c(i,ict,1)

      return
      end

c-----------------------------------------------------------------------
      subroutine copy_lcao(iadiag)
      use coefs, only: coef, nbasis, norb
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb

      do 20 i=1,norb+nadorb
       do 20 j=1,nbasis
   20   coef(j,i,iadiag)=coef(j,i,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_ci(iadiag)
      use csfs, only: ccsf, ncsf, nstates

      use dets, only: cdet, ndet
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'


      do 30 j=1,nstates
        do 30 i=1,ndet
   30     cdet(i,j,iadiag)=cdet(i,j,1)

      do 40 j=1,nstates
       do 40 icsf=1,ncsf
   40   ccsf(icsf,j,iadiag)=ccsf(icsf,j,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_zex(iadiag)

      use coefs, only: nbasis
      use basis, only: zex

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'


      do 20 i=1,nbasis
   20   zex(i,iadiag)=zex(i,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine save_jastrow_best

      use atom, only: nctype

      use jaspar3, only: a, b, c

      use jaspar4, only: a4, norda, nordb, nordc
      implicit real*8(a-h,o-z)




      include 'vmc.h'
      include 'force.h'


      dimension a4_best(MORDJ1,MCTYPE,MWF),b_best(MORDJ1,2,MWF),
     &c_best(83,MCTYPE,MWF)

      save a4_best,b_best,c_best
      save mparmja,mparmjb,mparmjc

c Save parameters corresponding to run generating hessian

      mparmja=2+max(0,norda-1)
      mparmjb=2+max(0,nordb-1)
      mparmjc=nterms4(nordc)

      do 50 ict=1,nctype
        do 50 i=1,mparmja
   50     a4_best(i,ict,1)=a4(i,ict,1)
      do 60 i=1,mparmjb
   60   b_best(i,1,1)=b(i,1,1)
      do 70 ict=1,nctype
        do 70 i=1,mparmjc
   70     c_best(i,ict,1)=c(i,ict,1)

      return

      entry restore_jastrow_best

c Restore parameters corresponding to run generating hessian
      do 80 ict=1,nctype
        do 80 i=1,mparmja
   80     a4(i,ict,1)=a4_best(i,ict,1)
      do 90 i=1,mparmjb
   90   b(i,1,1)=b_best(i,1,1)
      do 100 ict=1,nctype
        do 100 i=1,mparmjc
  100     c(i,ict,1)=c_best(i,ict,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine save_lcao_best
      use optwf_contrl, only: ioptorb
      use coefs, only: coef, nbasis, norb
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'



      dimension coef_best(MBASIS,MORB,MWF)

      save coef_best

      do 10 i=1,norb
       do 10 j=1,nbasis
   10   coef_best(j,i,1)=coef(j,i,1)

      return

      entry restore_lcao_best

      if(ioptorb.eq.0) return

      do 20 i=1,norb
       do 20 j=1,nbasis
   20   coef(j,i,1)=coef_best(j,i,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine save_ci_best
      use csfs, only: ccsf, ncsf, nstates

      use dets, only: cdet, ndet
      use optwf_contrl, only: ioptci
      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'



      dimension cdet_best(MDET,MSTATES),ccsf_best(MDET,MSTATES)
      save cdet_best,ccsf_best

      do 10 j=1,nstates
        do 10 i=1,ndet
   10     cdet_best(i,j)=cdet(i,j,1)

      do 20 j=1,nstates
       do 20 icsf=1,ncsf
   20   ccsf_best(icsf,j)=ccsf(icsf,j,1)

      return

      entry restore_ci_best

      if(ioptci.eq.0) return

      do 30 j=1,nstates
        do 30 i=1,ndet
   30     cdet(i,j,1)=cdet_best(i,j)

      do 40 j=1,nstates
       do 40 icsf=1,ncsf
   40   ccsf(icsf,j,1)=ccsf_best(icsf,j)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_parameters(dparm,iflag,iadiag)
      implicit real*8(a-h,o-z)

      dimension dparm(*)

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
      use jaspar3, only: a, b, c, scalek
      use jaspar4, only: a4
      use optwf_contrl, only: ioptjas
      use optwf_nparmj, only: nparma, nparmb, nparmc
      use optwf_wjas, only: iwjasa, iwjasb, iwjasc
      
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      dimension dparm(*)

      if(ioptjas.eq.0) return

c Set up cusp conditions
      call cuspinit4(0)

c Add change to old parameters
      iparm=0
      do 50 ict=1,nctype
        do 50 i=1,nparma(ict)
          iparm=iparm+1
   50     a4(iwjasa(i,ict),ict,iadiag)=a4(iwjasa(i,ict),ict,iadiag)-dparm(iparm)
      do 60 i=1,nparmb(1)
        iparm=iparm+1
   60   b(iwjasb(i,1),1,iadiag)=b(iwjasb(i,1),1,iadiag)-dparm(iparm)
      do 70 ict=1,nctype
        do 70 i=1,nparmc(ict)
          iparm=iparm+1
   70     c(iwjasc(i,ict),ict,iadiag)=c(iwjasc(i,ict),ict,iadiag)-dparm(iparm)
      call cuspexact4(0,iadiag)

c Check parameters a2 and b2 > -scalek
      call check_parms_jas(iflag)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_lcao(dparm,iadiag)
      use optwf_contrl, only: ioptorb
      use optwf_parms, only: nparmd, nparmj
      use coefs, only: coef, nbasis, norb
      use optorb_cblock, only: norbterm
      use orb_mat_022, only: ideriv

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'



      dimension acoef(MBASIS,MORB),dparm(*)

      if(ioptorb.eq.0) return

      do 10 i=1,norb
       do 10 j=1,nbasis
 10     acoef(j,i)=coef(j,i,iadiag)

c Update the orbitals
      do 30 i=1,norbterm
       io=ideriv(1,i)
       jo=ideriv(2,i)
       do 30 j=1,nbasis
 30     acoef(j,io)=acoef(j,io)-dparm(i+nparmj+nparmd)*coef(j,jo,iadiag)

      do 50 i=1,norb
       do 50 j=1,nbasis
 50     coef(j,i,iadiag)=acoef(j,i)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_ci(dparm,iadiag)
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use dets, only: cdet, ndet
      use optwf_contrl, only: ioptci, ioptjas, ioptorb
      use optwf_parms, only: nparmj
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      dimension dparm(*)

      if(ioptci.eq.0) return

c Update the ci coef
      if((method.eq.'linear'.or.method.eq.'lin_d').and.ioptjas+ioptorb.eq.0) then
        do 31 k=1,nstates

          if(ncsf.eq.0) then
            do 10 idet=1,ndet
              cdet(idet,k,1)=dparm(idet+ndet*(k-1))
 10         continue
           else
            do 15 j=1,ndet
 15           cdet(j,k,1)=0
            do 30 icsf=1,ncsf
              do 20 j=iadet(icsf),ibdet(icsf)
                jx=icxdet(j)
                cdet(jx,k,1)=cdet(jx,k,1)+dparm(icsf+ncsf*(k-1))*cxdet(j)
 20           continue
              ccsf(icsf,k,1)=dparm(icsf+ncsf*(k-1))
 30         continue
          endif

 31     continue
       else
         if(ncsf.eq.0) then
           do 35 idet=2,ndet
             cdet(idet,1,iadiag)=cdet(idet,1,iadiag)-dparm(idet-1+nparmj)
 35        continue
          else
           do 50 icsf=2,ncsf
             do 40 j=iadet(icsf),ibdet(icsf)
               jx=icxdet(j)
               cdet(jx,1,iadiag)=cdet(jx,1,iadiag)-dparm(icsf-1+nparmj)*cxdet(j)
 40          continue
             ccsf(icsf,1,iadiag)=ccsf(icsf,1,iadiag)-dparm(icsf-1+nparmj)
 50        continue
         endif
      endif

c     do 90 j=1,nstates
c90     write(6,'(''csf ='',1000f20.15)') (ccsf(i,j,iadiag),i=1,ncsf)

      return
      end
c-----------------------------------------------------------------------
      subroutine check_parms_jas(iflag)

      use atom, only: nctype

      use jaspar3, only: a, b, scalek

      use jaspar4, only: a4
      use optwf_nparmj, only: nparma, nparmb
      use optwf_wjas, only: iwjasa, iwjasb
      implicit real*8(a-h,o-z)







      include 'vmc.h'
      include 'force.h'



      iflag=0
      iflaga=0
      iflagb=0

      scalem=-scalek(1)
      do 50 ict=1,nctype
        do 50 i=1,nparma(ict)
   50     if(iwjasa(i,ict).eq.2.and.a4(2,ict,1).le.scalem) iflaga=1
      if(iflaga.eq.1) then
        do 55 ict=1,nctype
   55     write(6,'(''a2 < -scalek'',f10.5)') a4(2,ict,1)
      endif
      do 60 i=1,nparmb(1)
   60   if(iwjasb(i,1).eq.2.and.b(2,1,1).le.scalem) iflagb=1
      if(iflagb.eq.1) write(6,'(''b2 < -scalek'',f10.5)') b(2,1,1)
      
      if(iflaga.eq.1.or.iflagb.eq.1) iflag=1

      return
      end
c-----------------------------------------------------------------------
      subroutine test_solution_parm(nparm,dparm,
     &              dparm_norm,dparm_norm_min,add_diag,iflag)

      implicit real*8(a-h,o-z)

      include 'vmc.h'

      dimension dparm(*)

      iflag=0
      if(add_diag.le.0.d0) return

c Calculate rms change in parameters
      dparm_norm=0
      do 30 i=1,nparm
  30    dparm_norm=dparm_norm+dparm(i)**2
      dparm_norm=sqrt(dparm_norm/nparm)

      write(6,'(''dparm_norm,adiag ='',3g12.5)') 
     &dparm_norm,add_diag

      if(dparm_norm.gt.dparm_norm_min) iflag=1
      
      return
      end
c-----------------------------------------------------------------------
      subroutine save_nparms

      use optwf_contrl, only: ioptci, ioptjas, ioptorb
      use optwf_parms, only: nparmd, nparmj
      use optorb_cblock, only: norbterm
      use ci000, only: iciprt, nciprim, nciterm

      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'force.h'
      include 'optci.h'
      include 'optorb.h'


      save nparmj_sav,norbterm_sav,nciterm_sav,nparmd_sav,nreduced_sav

      nparmj_sav=nparmj
      norbterm_sav=norbterm
      nreduced_sav=nreduced
      nciterm_sav=nciterm
      nparmd=max(nciterm-1,0)
      nparmd_sav=nparmd

      write(6,'(''Saved max number of parameters, nparmj,norb,nciterm,nciterm-1: '',5i5)') nparmj,norbterm,nciterm,nparmd
      return

      entry set_nparms

      nparmj=nparmj_sav
      nparmd=nparmd_sav
      norbterm=norbterm_sav
      nreduced=nreduced_sav
      nciterm=nciterm_sav

      if(ioptjas.eq.0) nparmj=0
      if(ioptorb.eq.0) then
        norbterm=0
        nreduced=0
      endif
      if(ioptci.eq.0) then
        nciterm=0
        nparmd=0
      endif

      write(6,'(''Max number of parameters set to nparmj,norb,nciterm,nciterm-1: '',5i5)') nparmj,norbterm,nciterm,nparmd
      call set_nparms_tot

      return
      end
c-----------------------------------------------------------------------
      subroutine set_nparms_tot

      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      use optwf_parms, only: nparmd, nparmj
      use optorb_cblock, only: norbterm
      use ci000, only: iciprt, nciprim, nciterm

      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'force.h'
      include 'optjas.h'
      include 'optci.h'
      include 'optorb.h'
        

c Note: we do not vary the first (i0) CI coefficient unless a run where we only optimize the CI coefs

      if(method.eq.'sr_n') then

        nparmd=max(nciterm-1,0)
        nparm=nparmj+nparmd+norbterm

      elseif(method.eq.'linear'.or.method.eq.'lin_d') then
        
       i0=0
       if(ioptci.ne.0) i0=1
       if(ioptjas.eq.0.and.ioptorb.eq.0) i0=0

       nparmd=max(nciterm-1,0)
       nparm=nparmj+norbterm+nciterm-i0

      endif

      write(6,'(/,''number of parms: total, Jastrow, CI, orbitals= '',4i5)') 
     & nparm,nparmj,nciterm,norbterm

      return
      end
c-----------------------------------------------------------------------
      subroutine optwf_store(l,wt,psid,energy)
c store elocal and derivatives of psi for each configuration (call in vmc)

      use csfs, only: nstates
      use derivjas, only: gvalue
      use optwf_contrl, only: ioptci, ioptjas, ioptorb
      use optwf_func, only: ifunc_omega
      use optwf_parms, only: nparmj
      use sr_mat_n, only: elocal, nconf_n, sr_ho
      use sr_mat_n, only: sr_o, wtg
      use deloc_dj_m, only: denergy
      use force_analy, only: iforce_analy, iuse_zmat, alfgeo
      use optorb_cblock, only: norbterm
      use orb_mat_001, only: orb_ho, orb_o

      use ci000, only: iciprt, nciprim, nciterm

      use ci001_blk, only: ci_o, ci_oe
      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optjas.h'
      include 'optorb.h'
      include 'optci.h'
      include 'optci_cblk.h'
      include 'sr.h'

      dimension tmp_ho(MPARMJ),wt(*),psid(*),energy(*)

      call p2gtid('optgeo:izvzb',izvzb,0,1)
      call p2gtid('optwf:sr_rescale',i_sr_rescale,0,1)

      if(iforce_analy.gt.0.and.izvzb.eq.1) call force_store(l)

      if((method.ne.'sr_n'.and.method.ne.'lin_d').or.ioptjas+ioptorb+ioptci.eq.0)return

      i0=1
      if(method.eq.'lin_d'.and.ioptjas+ioptorb.eq.0) i0=0

      if(l.gt.MCONF) call fatal_error('SR_STORE: l gt MCONF')

      call dcopy(nparmj,gvalue,1,sr_o(1,l),1)

      ntmp=max(nciterm-i0,0)
      call dcopy(ntmp,ci_o(1+i0),1,sr_o(nparmj+1,l),1)

      ijasci=nparmj+ntmp
      if(ijasci+nstates*norbterm+nstates.gt.MPARM) call fatal_error('SR_STORE: iparm gt MPARM')

      do istate=1,nstates
        ii=ijasci+(istate-1)*norbterm
        call dcopy(norbterm,orb_o(1,istate),1,sr_o(ii+1,l),1)
        elocal(l,istate)=energy(istate)
        wtg(l,istate)=wt(istate)
      enddo

      ii=ijasci+nstates*norbterm
      do istate=1,nstates
        sr_o(ii+istate,l)=psid(istate)
      enddo
      
      nconf_n=l

      if(method.eq.'sr_n'.and.i_sr_rescale.eq.0.and.izvzb.eq.0.and.ifunc_omega.eq.0) return

c TO FIX: we are assuming optjas.ne.0 or optorb.ne.0 -> Otherwise, standard secular problem

      do 10 j=1,nparmj
  10    tmp_ho(j)=denergy(j,1)+gvalue(j)*energy(1)

      call dcopy(nparmj,tmp_ho,1,sr_ho(1,l),1)

      call dcopy(ntmp,ci_e(1+i0),1,sr_ho(nparmj+1,l),1)

      call dcopy(norbterm,orb_ho(1,1),1,sr_ho(nparmj+ntmp+1,l),1)
      
      return
      end
c-----------------------------------------------------------------------
