c-----------------------------------------------------------------------
      subroutine write_wf(iwf_fit,iter)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      character*40 filetype,wf,itn

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
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      character*50 fmt
      character*40 filename,filetype

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(83,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc
      common /bparm/ nspin2b,nocuspb

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj
      common /optwf_wjas/ iwjasa(83,MCTYP3X),iwjasb(83,3),iwjasc(83,MCTYPE),iwjasf(15,MCTYPE)
      common /optwf_nparmj/ nparma(MCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE)

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
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'numbas.h'

      common /numbas/ arg(MCTYPE),r0(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS,MCTYPE)

      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

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
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      character*40 filename,filetype

      common /dorb/ iworbd(MELEC,MDET)

      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      common /multidet/ kref,numrep_det(MDET,2),irepcol_det(MELEC,MDET,2),ireporb_det(MELEC,MDET,2)
     & ,iwundet(MDET,2),iactv(2),ivirt(2)

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj

      if(ioptci.eq.0) return

      filename='det'//filetype(1:index(filetype,' ')-1)
      open(2,file=filename,status='unknown')

      call p2gti('electrons:nelec',nelec,1)
      call p2gti('electrons:nup',nup,1)
      write(2,'(''&electrons nelec '',i4,'' nup '',i4)') nelec,nup
      write(2,'(''kref'',i4)') kref
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
