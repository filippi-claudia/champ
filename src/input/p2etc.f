c **********************************************************************
c p2etc.F
C     $Revision: $
c     routines for p2 that could be replaced by 
c     more adapted ones of the program that uses the p2 parser
c
c **********************************************************************
      subroutine noxxx(i)
      implicit double precision (a-h,o-z)
      common /noxx1/ nox
      nox=i
      end 
      subroutine fatal(msg)
      implicit double precision (a-h,o-z)
      common /noxx1/ nox
      character msg*(*)
      write (6,*) ' '
      print 10,msg
      write (6,*) ' '
 10   format('*** FATAL : ',a)
      if(nox.eq.0) then
       call theend(1)
      endif
      end
      subroutine warn(msg)
      common /warni/ iwrn
      character msg*(*)
      write (6,*) ' '
      print 10,msg
      write (6,*) ' '
      iwrn=iwrn+1
 10   format('*** WARNING : ',a)
      end

      subroutine infox(iprt0,itrce0)
C$INPUT info i i=-1
      implicit double precision (a-h,o-z)
      include 'inc/debug.inc'
      if(itrce0.ne.-1) itrce=itrce0
      iprt=iprt0
      print 10,iprt,itrce
 10   format ('info: print=',i3,' trace=',i3)
      end 

      subroutine intro(nm)
      implicit double precision (a-h,o-z)
      character nm*(*)
      include 'inc/debug.inc'
c
c      call atime(t1)
      t1=0
      ittop=ittop+1
      if(ittop.le.itrce) then
       print 10,ittop,nm
      endif
      if(ittop.gt.MXTSTCK) then
       call fatal('t-stack full')
      endif
      tim1(ittop)=t1
 10   format('[',i2,'] entering ',a,' ... ')
      end

      subroutine outro(nm)
      implicit double precision (a-h,o-z)
      character nm*(*)
      include 'inc/debug.inc'
c      call atime(tim2)
      tim2=0
      if(ittop.gt.0) then
       if(ittop.le.itrce) then
        print 10,ittop,nm,tim2-tim1(ittop)
       endif
       ittop=ittop-1
      else
       call fatal('t-stack empty')
      endif
 10   format('[',i2,']  leaving ',a,' (',f8.2,' s)')
      end

      subroutine fexit
      implicit double precision (a-h,o-z)
      call theend(1)
      end 

      subroutine theend(ier)
      implicit double precision (a-h,o-z)
      common /warni/ iwrn
c
      call flecle
c
      if(ier.eq.0) then
       if(iwrn.eq.0) then
        ixcode=0
       else
        print 100,iwrn
        ixcode=2
       endif
      else
       ixcode=1
      endif
 100  format ('** ',i5,' warning(s)')
      stop
      end 

c***  files ***
      subroutine nfile(iu,basen,num,ext,st,itmp,ity)
      implicit double precision (a-h,o-z)
      character basen*(*)
      character ext*(*)
      character st*(*)
      parameter(nlen=12,mlen=64)
      character tmp*(nlen)
      character fnam*(mlen)
      logical ex
c
      write(tmp,10) num
 10   format(i12)
      id1=idxfn(' ',tmp,nlen)
      id2=idxln(' ',tmp,nlen)
ctest>
c$$$      print 876,basen
c$$$      print 876,tmp(id1:id2)
c$$$      print 876,ext
c$$$ 876  format('--> [',a,']')
ctest<
      write(fnam,11) basen,tmp(id1:id2),ext
 11   format(a,a,'.',a)
      id1=idxfn(' ',fnam,mlen)
      id2=idxln(' ',fnam,mlen)
c      print 12,fnam(id1:id2)
c 12   format('FILENAME= ',a)
      if(st.eq.'inq') then
       inquire(file=fnam(id1:id2),exist=ex)
       if(ex) then
        call file(iu,fnam(id1:id2),'old',itmp,ity)
       else
        iu=0
       endif
      else
       call file(iu,fnam(id1:id2),st,itmp,ity)
      endif
      end 
      subroutine ptfile(iu,fn,st)
      implicit double precision (a-h,o-z)
      character fn*(*)
      character st*(*)
      call file(iu,fn,st,1,0)
      end 
      subroutine file(iu,fn,st,itmp,ity)
      implicit double precision (a-h,o-z)
      include 'inc/files.inc'
      character fn*(*)
      character st*(*)
      character s1*12
      logical ex
      common /inidon_file/inidon
      if( inidon.eq.0 ) then
       call  initf
       inidon = 1
      endif
c     status st: old,new,app,cls
c     itmp=1: permanent, itmp=2 temporary file
c     ity=0 text, ity=1 binary
c
      call intro('file')
      if(fn.eq.'stdin')then
       iu=5
       goto 999
      endif
      if(fn.eq.'stdout')then
       iu=6
       goto 999
      endif
      if(fn.eq.'<input>') then
       if(inptu.ne.0) then
        iu=inptu
        goto 999
       else
        call fatal('file: <input> not defined')
       endif
      endif
      k=0
      do i=1,nfle
       if((ifn(i).eq.fn).and.(istat(i).ne.0))then
        k=i
        goto 99
       endif
      enddo
 99   continue 
      if(st.eq.'old')then
       if(k.ne.0) then
        goto 200
       else
        s1='OLD'
        goto 100
       endif
      elseif(st.eq.'new')then
       if(k.ne.0)then
        print 11,fn
 11     format ('ERROR: file ',a,' already open')
        call fatal('file error') 
       else
        s1='NEW'
        goto 100
       endif
      elseif(st.eq.'app')then
       if(k.ne.0) then
        goto 200
       else
        s1='unknown'
        goto 100
       endif
      elseif((st.eq.'cls').or.(st.eq.'del'))then
       if(k.eq.0)then
        print 12,fn
 12     format ('WARNING: file ',a,' not open -- cannot close it')
        call warn('file error')  
       else
        if(istat(k).eq.1)then
         if(st.eq.'cls')then
          close(ifu(k),status='keep')
         else
          close(ifu(k),status='delete')
         endif
        elseif(istat(k).eq.2)then
         close(ifu(k),status='delete')
        else
         stop 'file:status?'
        endif
        istat(k)=0
       endif
      else
       print 20,st
 20    format ('sr file: unknown flag ',a)
       call fatal('file error') 
      endif
      goto 999
c     open new file...
 100  continue 
c     try to find empty slot ..
      do i=1,nfle
       if(istat(i).eq.0) then
        kfle=i
        iu=ifu(i)
        goto 101
       endif
      enddo
c     new slot ...
      nfle=nfle+1
      if(nfle.gt.MXFLS)then
       call fatal('sr file: too many files')
      endif
      iul=iul+1 
      if(iul.gt.MXUNT) then
       write (6,*) '> ',iul,MXUNT
       call fatal('file: unit number ot of range')
      endif
      kfle=nfle
      iu=iul
 101  continue 
      ifu(kfle)=iu
      ifn(kfle)=fn
      istat(kfle)=itmp
      itype(kfle)=ity
      icpos(iu)=0
      if((s1.eq.'OLD').or.(s1.eq.'old')) then
       inquire(file=fn,exist=ex)
       if(.not.ex) then
        print 31,fn
 31     format('file ',a,' does not exist')
        call fatal('file error') 
       endif
      endif
      if(ity.eq.0) then
       open(iu,file=fn,status=s1)
      elseif(ity.eq.1) then
       open(iu,file=fn,access='sequential',form='unformatted'
     $      ,status=s1)
      else
       write (6,*) 'bad file type  ',ity
       call fatal('file error') 
      endif
      goto 999
c     existing file...
 200  continue
      if(itype(k).eq.ity) then
       iu=ifu(k)
       if((irewnd.eq.1).and.(st.eq.'old')) then
        rewind(iu)
       endif
       goto 999
      else
       print 30,fn,itype(k),ity
 30    format ('sr file: wrong type. file=',a,' type=',i2,' req=',i2)
       call fatal('file error') 
      endif
 999  continue 
      call outro('file')
      end 
      subroutine setpos(iu,npos)
      implicit double precision (a-h,o-z)
      include 'inc/files.inc'
      icpos(iu)=npos
      end 
      subroutine incpos(iu,il,incr)
      implicit double precision (a-h,o-z)
C    increment and return file posision counter
      include 'inc/files.inc'
      icpos(iu)=icpos(iu)+incr
      il=icpos(iu)
      end 
      subroutine flecle
      implicit double precision (a-h,o-z)
c     close all files
      include 'inc/files.inc'
      if(nfle.gt.0) then
       do ifle=1,nfle
        if(istat(ifle).eq.1)then
         close(ifu(ifle),status='keep')
        elseif(istat(ifle).eq.2)then
         close(ifu(ifle),status='delete')
        endif
       enddo
      endif
      end 

      subroutine finfo(iflag)
C INPUT finfo i=0
      implicit double precision (a-h,o-z)
      include 'inc/files.inc'
      character s*10
      character ts*4
      write (6,*) '** FILE INFO **'
      print 9
      do  ifle=1,nfle
       if(istat(ifle).eq.0) then
        s='CLOSED   '
       elseif(istat(ifle).eq.1) then
        s='PERMANENT'
       elseif(istat(ifle).eq.2) then
        s='TEMPORARY'
       else
        write (6,*) ifle,istat(ifle)
        stop 'finfo:status?' 
       endif
       k1=idxfn(' ',ifn(ifle),MXFNL)
       k2=idxln(' ',ifn(ifle),MXFNL)
       if(itype(ifle).eq.0) then
        ts='TEXT'
       elseif(itype(ifle).eq.1) then
        ts='DATA'
       else
        write(ts,'(i4)') itype(ifle)
       endif
c     when iflag is set, dont list closed files
       if((iflag.eq.0).and.(istat(ifle).eq.0)) goto 99
       print 10,ifle,ifu(ifle),ifn(ifle)(k1:k2),ts,s,
     $      icpos(ifu(ifle))
 99    continue 
      enddo
      print 11,icpos(5)
 9    format ('num  f-unit  filename type state curr.line') 
 10   format (i4,' ',i4,' ',a,' ',a,' ',a,' ',i8)
 11   format('stdin : current pos = ',i12)
      end 

      subroutine rwf(fn)
C INPUT rewind a
      implicit double precision (a-h,o-z)
      include 'inc/files.inc'
      character fn*(*)
      if(fn.eq.'stdin') then
       rewind(5)
       call setpos(5,0)
       return 
      endif
      if(fn.eq.'stdout') then
       call warn('cannot rewind stdout')
       return
      endif
      k=0
      do i=1,nfle
       if((ifn(i).eq.fn).and.(istat(i).ne.0))then
        k=i
        goto 99
       endif
      enddo
 99   continue
      if(k.eq.0) then
       k1=idxfn(' ',ifn(i),MXFNL)
       k2=idxln(' ',ifn(i),MXFNL)
       print 11,ifn(i)(k1:k2)
 11    format ('rewind: file ',a,' not open')
       call fatal('file error') 
      else
       rewind(ifu(k))
       call setpos(ifu(k),0)
      endif
      end 


c     array ops
      subroutine cparri(ndim,isrc,itgt)
      implicit double precision (a-h,o-z)
      dimension isrc(ndim)
      dimension itgt(ndim)
      do i=1,ndim
       itgt(i)=isrc(i)
      enddo
      end 
      subroutine cparrf(ndim,fsrc,ftgt)
      implicit double precision (a-h,o-z)
      dimension fsrc(ndim)
      dimension ftgt(ndim)
      do i=1,ndim
       ftgt(i)=fsrc(i)
      enddo
      end 

CC********** conversion functions *******************
      function ifroms(s,ifl)
      implicit double precision (a-h,o-z)
      character s*(*)
      read(s,*,err=99) i
      ifroms=i
      return
 99   continue
      print 10,s
 10   format ('>> [',a,']')
      if(ifl.eq.0) then
       call warn('ifroms: conversion error')
      else
       call fatal('ifroms: conversion error')
      endif
      ifroms=0
      end 

      function dfroms(s,ifl)
      implicit double precision (a-h,o-z)
      character s*(*)
      read(s,*,err=99) f
      dfroms=f
      return
 99   continue
      print 10,'> ',s
 10   format (2a)
      if(ifl.eq.0) then
       call warn('dfroms: conversion error')
      else
       call fatal('dfroms: conversion error')
      endif
      dfroms=0.d0
      end 


C    clear vector
      subroutine iclear(k,n)
      implicit double precision (a-h,o-z)
      dimension k(n)
      do i=1,n
       k(i)=0
      enddo
      end 
      subroutine fclear(v,n)
      implicit double precision (a-h,o-z)
      dimension v(n)
      do i=1,n
       v(i)=0.d0
      enddo
      end 
C    read out vector elements
      subroutine  getvc1(src,np,tgt,np2,n1,n2)
      implicit double precision (a-h,o-z)
      dimension src(np)
      dimension tgt(np2)
      ix=1
      do i=n1,n2
       tgt(ix)=src(i)
       ix=ix+1
      enddo
      end 
      subroutine  getvc2(src,npa,npb,tgt,np2,idx2,n1,n2)
      implicit double precision (a-h,o-z)
      dimension src(npa,npb)
      dimension tgt(np2)
      ix=1
      do i=n1,n2
       tgt(ix)=src(i,idx2)
       ix=ix+1
      enddo
      end 
      subroutine arwnd(flag,i1)
C INPUT autorewind a=on i=0
      implicit double precision (a-h,o-z)
      include 'inc/files.inc'
      character flag*(*)
      call intro('arwnd')
      if(flag.eq.'on') then
       irewnd=1
      elseif(flag.eq.'off') then
       irewnd=0
      elseif(flag.eq.'arg') then
       itmp=irewnd
       irewnd=i1
       i1=itmp
      elseif(flag.eq.'?') then
       print 10,irewnd
 10    format ('autorewind= ',i2)
      else
       call fatal('arwnd: bad flag')
      endif
      call outro('arwnd')
      end 

      subroutine setinp(iu)
      implicit double precision (a-h,o-z)
      include 'inc/files.inc'
      inptu=iu
      end 

      subroutine initf
      implicit double precision (a-h,o-z)
      include 'inc/files.inc'
      nfle=0
      do i=1,MXFLS
       ifu(i)=0
       istat(i)=0
       itype(i)=0
      enddo
      do i=1,MXUNT
       icpos(i)=0
      enddo
      irewnd=1
      end 

      subroutine echcmd(cmd)
      implicit double precision (a-h,o-z)
      character cmd*(*)
      if(cmd.eq.'loop')  return
      if(cmd.eq.'end')   return
      if(cmd.eq.'gotol') return
      print 10
      print 11,cmd
      print 12
 10   format
     $('------------------------------------------------------------')
 11   format('Command ',a)
 12   format(' ')
      end

C **********************************************************************
      subroutine p2chkd(iu, delim, ierrflg, ierrtrn)
C **********************************************************************
c
c  Test whether next input line consists of delim.
c  
c Input: 
c   iu:      input unit
c   delim:   expected string
c   ierrflg: flag :  1 stop on error and print message
c
c Output (only if ierrflg not 1):
c ierrtrn: 0   delim found
c          1   delim not found
c          2   more than one field in line
C **********************************************************************
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      dimension idx1(MXF),idx2(MXF)
      character delim*(*)
      call gpnln2(iu,il,lne,MXLNE,idx1,idx2,mm,lbuf,'p2chkd')
      ierr=0
      if( mm .ne. 1 ) then
       ierr=2
      endif
      if( lne(idx1(1):idx2(1)) .ne. delim) then
       ierr=1
      endif
      if( (ierr.ne.0).and.(ierrflg.eq.1)) then
       print 10,iu,il
       print 11,delim
       print 12,lne(idx1(1):idx2(mm))
       stop
      endif
      if(ierrflg.ne.1) ierrtrn=ierr
      return
 10   format('p2chkd: Error : expected string not found.',
     $     ' unit ',i4,' line ',i6)
 11   format('expected : ',a)
 12   format('found : ',a)
      end  


      subroutine fn_open(fn,vn,flag)
C$INPUT open_file a a=f a=app
CKEYDOC {\tt open_file} {\sl filename} [ {\sl varname} [ {\sl flag}]] \\
CKEYDOC Open file and store assigned unit number in a variable.
      implicit double precision (a-h,o-z)
      character fn*(*)
      character vn*(*)
      character flag*(*)
      call ptfile(iu, fn, flag)
      call p2sti(iu, vn)
      print 10,fn,iu,vn
 10   format('fn_open: file= ',a,' unit= ',i4,' variable= ',a) 
      end 
C **********************************************************************
      subroutine p2chkend(iu, keyw)
C **********************************************************************
c
c  Test whether next input line consists the word end.
c  and issue an error if this is not the case
c  
c Input: 
c   iu:      input unit
c   keyw:    keyword for error message
C **********************************************************************
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      dimension idx1(MXF),idx2(MXF)
      character keyw*(*)
      call gpnln2(iu,il,lne,MXLNE,idx1,idx2,mm,lbuf,'p2chkend')
      ierr=0
      if( mm .ne. 1 ) then
       ierr=2
      endif
      if( .not.((lne(idx1(1):idx2(1)) .eq. 'end').or.
     $     (lne(idx1(1):idx2(1)) .eq. 'END'))) then
       ierr=1
      endif
      if(ierr.ne.0) then
       print 10,il,keyw
       print 12,lne(idx1(1):idx2(mm))
       stop
      endif
      return
 10   format('Error: <end> expected on line ',i4,' keyword ',a)
 12   format('( found : ',a,')')
      end  

      block data p2etc_data
      implicit double precision (a-h,o-z)

      common /warni/ iwrn
      common /in0/ itrce,iprt
      common /jt0/ ittop,idum003
      common /fle0/ nfle,iul
      common /fle2b/ inptu
      common /inidon_file/inidon

      data iwrn /0/
      data itrce /0/
      data iprt /0/
      data ittop /0/
      data nfle /0/
      data iul /50/
      data inptu /0/
      data inidon /0/
      end
