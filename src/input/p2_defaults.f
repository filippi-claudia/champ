c********************************************************************************
c     mini-parser p2
c     with support for macros, name lists and keywords with
c     default arguments
c     (c) Friedemann Schautz 1998-2002
c$Revision: 1.5 $
c********************************************************************************
      subroutine pver
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'  
      include 'inc/p2b.inc'  
      write (6,10)
      write (6,11) sfile
      write (6,12)
 10   format ('MINIPARSER P2 with support for macros and default args')
 11   format(' compiled command-list from file ',a) 
 12   format(' $Revision: 1.5 $')
      end 

      subroutine p2go(iu,isil)
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'  
      include 'inc/p2b.inc'  
      dimension idx1(MXF)
      dimension idx2(MXF)
ctest>  
      
      if(ip2deb.gt.0) then
       write (6,*) 'nesting level ',lnest
       write (6,*) 'parser dimensions'
       write (6,*) 'MAXKEY ',MXKEY
       write (6,*) 'MXVAR  ',MXVAR
       write (6,*) 'MXLNE  ',MXLNE
       write (6,*) 'MXIDL  ',MXIDL
       write (6,*) 'MXF    ',MXF
       write (6,*) 'MXDEF  ',MXDEF
      endif
ctest<
      iutmp=iu
      call setinp(iu)
cc
      lne=' '
cc
 1    continue 
      if(iimode.ne.0) then
       write (6,20) il,p2prmt
 20    format (i4,a)
      endif
      call gpnlne(iu,il,lne,nl,idx1,idx2,mm,lbuf,ierr)
      if(ip2deb.gt.0) then
       write (6,30) il,ierr,mm
       write (6,31) lne(idx1(1):idx2(mm))
 30    format('_p2deb_ ',3i6)
 31    format('[',a,']')
      endif
      if(ierr.eq.-1) goto 1
      if(ierr.eq.1)  goto 1000
      if(ierr.eq.2)  goto 1001
c     parse line ...
      call pline(iu,lne,il,fmts,itmp,ftmp,atmp,is1,is2,keys,keylen,nargs
     $     ,ierr,isil,ikw,nkey,idx1,idx2,mm)
      if(ierr.eq.0)then
       if(inice.eq.1) then
        call echcmd(keys(ikw)(1:keylen(ikw)))
       endif
       if(ip2deb.gt.0) then
        write (6,*) 'calling p2call ',ikw
        write (6,*) 'first string arg : ',is1(1),is2(1)
        write (6,*) '(''['',a,'']'')',atmp(is1(1):is2(1))
       endif 
       call p2call(ikw,itmp,ftmp,is1,is2,atmp,iend,MXF,iu)
       if(iend.eq.0)then
        goto 1
       else
        goto 1000
       endif
      elseif(ierr.eq.-1)then
       goto 1
      else
       if(isil.eq.1)then
        write (6,10) ierr
       endif
 10    format ('P2: error ',i4)
       if(iimode.eq.0) then
        call fatal('input error')
       else
        goto 1
       endif
      endif
 1000 continue 
      call ectrl2
      write (6,11)
 11   format ('***')
      return
 1001 continue
      write (6,103) il
 103  format('Unbalanced quote on input line ',i6)
      call fatal('input error') 
      end 

      subroutine gpnln2(iu,il,lne,nl,idx1,idx2,mm,lbuf,errmsg)
      implicit double precision (a-h,o-z)
C    wrapper around gpnlne which
C    eats comments and blank lines
C    and reports errors 
C    returns 'end' at end-of-file
      include 'inc/p2_dim.inc'
      character  lne*(MXLNE)
      character  lbuf*(MXLNE)
      dimension idx1(MXF)
      dimension idx2(MXF)
      character errmsg*(*)
 1    continue 
      call gpnlne(iu,il,lne,MXLNE,idx1,idx2,mm,lbuf,ierr)
      if(ierr.eq.0) then
C blank line
       if(mm.eq.0) goto 1
      elseif(ierr.eq.-1) then
C comment line
       goto 1
      elseif(ierr.eq.1) then
C eof
       goto 2000
      else
C some other error
       goto 1000
      endif
      return
 1000 continue 
      write (6,5) il,errmsg
      call fatal('input error')
 5    format ('gpnln2: error at line ',i4,' called from ',a)  
      return
 2000 continue 
      mm=1
      idx1(1)=1
      idx2(1)=3
      lne(1:3)='end'
      return
      end 

      subroutine gpnlne(iu,il,lne,nl,idx1,idx2,mm,lbuf,ierr)
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2b.inc'
C     get & preprocess next input line
c     field sep. " " and ","
c     comments start with cmtchr
c     anticomments start with antch0 (if iant>0)
c     macros start with "$"
c     "\" continues line
      character  lne*(MXLNE)
      character  lbuf*(MXLNE)
      dimension idx1(MXF)
      dimension idx2(MXF)  
c
      ierr=0
      call rdln(iu,lne,il,ieof)
      if(ieof.eq.1) then
       ierr=1
       return
      endif
c     continuation line ??
 7    continue 
      kll=idxln(' ',lne,MXLNE)
      if(kll.gt.0) then
       if(lne(kll:kll).eq.'\\') then
        kl2=kll
        call rdln(iu,lne(kl2:),il,ieof)
        if(ieof.eq.1) then
         ierr=1
         return
        endif
        goto 7
       endif
      endif

c     split line
c      call ssplt(lne,idx1,idx2,MXLNE,MXF,mm,fsep,ifsep)
c   (  sspltq allows singly quoted strings  )
      call sspltq(lne,idx1,idx2,MXLNE,MXF,mm,fsep,ifsep,iqerr)
      if(iqerr.ne.0) then
       ierr=2
       return
      endif
c *** special lines (comment, blank)
c     blank line?
      if(mm.eq.0)then
       ierr=-1
       return
      endif
c     comment? skip?
ctest>
c$$$      write (6,222) lne(idx1(1):idx1(1)),cmtchr
c$$$ 222  format ('P222 [',a1,']  [',a1,']')
ctest<
      if(lne(idx1(1):idx1(1)).eq.cmtchr) then
       if(lne(idx1(1)+1:idx2(1)).eq.skipch) then
        iskip=0
       endif
       ierr=-1
       return
      endif
ccc *** macro expand
      call macex(lne,lbuf,MXLNE,idx1,idx2,MXF,mm,1,il)

      if(iskip.eq.1) then
ctest>
       if(ip2deb.ne.0) then
        write (6,*) 'SKIP ',iskip
       endif
ctest<
       ierr=-1
       return
      endif
ctest>
      if(ip2deb.ne.0) then
       write (6,100) il,mm,lne(idx1(1):idx2(mm))
      endif
 100  format ('### line ',i2,' nf',i3,'[',a,']')
ctest<
ccc *** 
      return
      end 
      
      subroutine rdln(iu,line,lnum,ieof)
C    read next line from input stream 
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      character line*(*)

 1     continue 
       read(iu,'(a)',err=999,end=1000) line
       if(iant.gt.0) then
        if(line(1:1).ne.antchr) then
         goto 1
        else
         line(1:1)=' '
        endif
       endif
       call incpos(iu,lnum,1)
       ieof=0
       return

      return
 999  continue 
      call fatal('rdln: Error while reading line')
      stop
 1000 continue 
      ieof=1
      return
      end 
      subroutine pline(iu,lne,lnr,fmts,itmp,ftmp,atmp,is1,is2,keys,
     $     keylen,nargs,ierr,isil,ikw,NKEY,idx1,idx2,mm)
      implicit double precision (a-h,o-z)
c     process one input line 
c     " as first char 'verboses'
c     "def" defines one-line macro
c     ierr: 0=success, -1=blank line
c           1=keyword not found, 2=conversion error,
c           3=missing argument(s), 4=extra arguments 
c           5=format error 6=variable not found
c           7=error in variable def.
c     isil=1 silent 
c     format string : i,d,a (integer,double,string*MXLNE)
      include 'inc/p2_dim.inc'     
      include 'inc/p2b.inc'     
      include 'inc/p2defv.inc'
      character lne *(*)
      character fmts(NKEY)*(MXIDL)
      dimension nargs(NKEY)
      dimension itmp(MXF)
      dimension ftmp(MXF)
      character atmp*(MXLNE)
      dimension is1(MXF)
      dimension is2(MXF)
      character keys(NKEY)*(MXIDL)
      dimension keylen(NKEY)
      dimension idx1(MXF),idx2(MXF)
      character astr*(MXLNE)
      character lbuf*(MXLNE)
      dimension ikwtmp(MXKEY)
      character pfxstr*(MXIDL)
      character vnm*(MXIDL)
      nl=MXLNE
c
ccc *** special commands 
      if(lne(idx1(1):idx1(1)).eq.'"') then
C    print whole line
       if(idx2(1).eq.idx1(1)) then
        write (6,201) lne(idx1(2):idx2(mm))
       else
        write (6,201) lne((idx1(1)+1):idx2(mm))
       endif
 201   format(a)
       ierr=-1
       return
      elseif(lne(idx1(1):idx2(1)).eq.'sh') then
C    system call (only on DEC-osf1 and linux)
       if(mm.gt.1) then
        write (6,*) ' '
        write (6,*) ' '
        write (6,202) lne(idx1(2):idx2(mm))
        write (6,*) ' '
        call fatal('no shell-command in this version')
        if(iret.ne.0) then
         call fatal('shell command failed')
        endif
       else
        call warn('empty shell command')
       endif
 202   format('executing shell command : ',a)  
       ierr=-1
       return
      elseif(lne(idx1(1):idx1(1)).eq.'&') then
C 'namelist' -like variable input
C expected format: & nam1 val1 nam2 val2 ...
C or               &prefix val1 nam2 val2 ...
       if((mm.lt.3).or.(mod((mm-1),2).ne.0)) then
        ierr=7
        goto 1000
       else
C prefix 
        if(idx2(1).gt.idx1(1)) then
         pfxstr=lne((idx1(1)+1):idx2(1))
         kplen=idx2(1)-idx1(1)
         ipfx=1
c$$$         write (6,888) kplen, pfxstr(1:kplen)
c$$$ 888     format ('prefix ',i8,' [',a,']')
        else
         ipfx=0
        endif
C loop over name-value pairs
        do i=2,mm,2
         if(ipfx.eq.1) then
          write(vnm,789) pfxstr(1:kplen), lne(idx1(i):idx2(i))
 789      format(a,':',a) 
          klen=kplen+(idx2(i)-idx1(i))+2
CC    check if this namelist should only contain
CC    predefined fields and if the current  field is defined
          call p2nmcheck(pfxstr(1:kplen),lne(idx1(i):idx2(i)),ierr)
          if(ierr.ne.0) then
           ierr=8
           goto 1000
          endif
         else
          write(vnm,790) lne(idx1(i):idx2(i))
 790      format(a)
          klen=idx2(i)-idx1(i)+1
         endif
c$$$         write (6,889) klen,vnm(1:klen)
c$$$ 889     format('VAR NAME ',i8, ' [',a,']')
         call p2setv(vnm(1:klen),lne(idx1(i+1):idx2(i+1)))
        enddo
        ierr=-1
        return
       endif
      elseif(lne(idx1(1):idx2(1)).eq.'def') then
C    macro definition
       if(mm.eq.2) then
        call p2setv(lne(idx1(2):idx2(2)),' ')
        ierr=-1
        return
       elseif(mm.gt.2) then
        call unprtct(lne,3,mm,idx1,idx2,MXF)
        call p2setv(lne(idx1(2):idx2(2)),lne(idx1(3):idx2(mm)))
        ierr=-1
        return
       else
        ierr=7
        goto 1000
       endif
      else
C    catch pragmas
       call p2prag(iu,lne(idx1(1):idx2(mm)),itest)
       if(itest.eq.1) then
        ierr=-1
        return
       endif
C    simple keyword search
       iabbr=0
       ikw=0
       len0=idx2(1)-idx1(1)+1
       do k=1,NKEY
C exact match
        if(len0.eq.keylen(k)) then
         if(keys(k)(1:keylen(k)).eq.lne(idx1(1):idx2(1)))then
          ikw=k
          goto 80
         endif
        elseif(len0.lt.keylen(k)) then
C abbrev
         if(ip2ab.eq.1) then
          if(keys(k)(1:len0).eq.lne(idx1(1):idx2(1)))then
           iabbr=iabbr+1
           ikwtmp(iabbr)=k
          endif
ctest>
c$$$         write (6,500) len0,k,keys(k)(1:len0),lne(idx1(1):idx2(1)),
c$$$     $         keys(k)(1:keylen(k))
c$$$         if(iabbr.gt.0) then
c$$$          write (6,501) iabbr,ikwtmp(iabbr)
c$$$         else
c$$$          write (6,*) 'iabbr =0'
c$$$         endif
c$$$ 500     format('_ABR1_ ',i4,i4,'[',a,']  [',a,']  [',a,'] ')
c$$$ 501     format('_ABR2_ ',i4,i4)
ctest<
         endif
        endif
       enddo
       if((ikw.eq.0).and.(ip2ab.eq.1)) then
        if(iabbr.eq.1) then
         ikw=ikwtmp(1)
         goto 80
        endif
       endif 
       ierr=1
       goto 1000
      endif
 80   continue 
c
      iarg=0
      ini=0
      ifl=0
      iaa=0
c    argument list ...
      if(mm.gt.1)then
       do i=2,mm
c       late comment ?
        if(lne(idx1(i):idx1(i)).eq.cmtchr)then
         goto 99
        endif
        iarg=iarg+1
        if(iarg.gt.nargs(ikw))then
         ierr=4
         goto 1000
        endif
        call p2isvar(lne,idx1,idx2,MXF,i,astr,k3,ie)
        if(ie.ne.0) then
         ierr=6
         goto 1000
        endif
        if((fmts(ikw)(iarg:iarg)).eq.'i') then
         read(astr(1:k3),*,err=999) it
         ini=ini+1
         itmp(ini)=it
        elseif((fmts(ikw)(iarg:iarg)).eq.'d') then
         read(astr(1:k3),*,err=999) dt
         ifl=ifl+1
         ftmp(ifl)=dt
        elseif((fmts(ikw)(iarg:iarg)).eq.'a') then
         iaa=iaa+1
         if(iaa.gt.1)then
          is1(iaa)=is2(iaa-1)+1
          is2(iaa)=is1(iaa)+k3-1
         else
          is1(iaa)=1
          is2(iaa)=k3
         endif
         atmp(is1(iaa):is2(iaa))=astr(1:k3)
        else
         ierr=5
         goto 1000
        endif
       enddo
      endif
 99   continue 
      if(iarg.eq.nargs(ikw)) then
       ierr=0
      else
C    insert default values ??
       if((ip2dfl.ne.0).and.(ideflt(ikw).gt.0).and.(ideflt(ikw).le.(iarg
     $      +1))) then
        call p2ddfl(ikw,iarg,ini,ifl,iaa,nargs(ikw),
     $       fmts(ikw),itmp,ftmp,atmp,is1,is2,
     $       idefpp,idefvv,ddefvv,adefvv,
     $       ierr)
       else
        ierr=3
        goto 1000
       endif
      endif
      if(ierr.eq.0) return
C errors ....
 999  continue 
      ierr=2
 1000 continue 
      if(isil.ne.1)then
       k2=idxln(' ',lne,nl)
c       write (6,8) ierr,lnr
c       write (6,9) lne(1:k2)
c 8     format ('P2: error ',i3,' at input line ',i5,' :') 
c 9     format(a)
       write (6,*) ' '
       write (6,7) lnr
 7     format ('*** ERROR at line ',i5,' :')
       if(ierr.ne.1.and.ikw.gt.0) then
        kl=idxf(' ',keys(ikw),MXIDL)
       endif
       if(ierr.eq.1) then 
        if((ip2ab.eq.1).and.(iabbr.gt.1)) then
         write (6,17) lne(idx1(1):idx2(1))
         do i=1,iabbr
          write (6,18) keys(ikwtmp(i))(1:keylen(ikwtmp(i)))
         enddo
 17      format('   keyword ambiguous: ',a)
 18      format('   ?? ',a)
        else
         write (6,10) lne(idx1(1):idx2(1))
 10      format ('   keyword not found : ',a)
        endif
       elseif(ierr.eq.2) then
        write (6,11) keys(ikw)(1:kl),iarg,fmts(ikw)(iarg:iarg),astr
 11     format ('   conversion error: keyw=',a,' arg# ',i4,' type=',a1
     $       ,' got: ',a)
       elseif(ierr.eq.3) then
        write (6,12) keys(ikw)(1:kl),nargs(ikw),iarg
 12     format ('   missing argument(s): keyw=',a,' need ',i4
     $       ,' got only ',i4)
       elseif(ierr.eq.4) then
        write (6,13) keys(ikw)(1:kl),nargs(ikw)
 13     format ('   extra argument(s): keyw=',a,' need only ',i4)
       elseif(ierr.eq.5) then
        write (6,14) keys(ikw)(1:kl),iarg,fmts(ikw)(iarg:iarg)
 14     format ('   format error: ',a,i4,a)
       elseif(ierr.eq.6) then
        write (6,15) lne(idx1(i):idx2(i))
 15     format ('   variable not found: ',a)
       elseif(ierr.eq.7) then
        write (6,16)
 16     format ('   error in variable assignment')
       elseif(ierr.eq.8) then
        write (6,19) lne(idx1(i):idx2(i)),pfxstr(1:kplen)
 19     format ('   Name ',a,' not defined in list ',a)
       else
        write (6,*) ierr
        stop 'P2 !' 
       endif
       write (6,*) ' '
      endif
      end 

      subroutine ssplt(s,idx1,idx2,n,m,mm,c,nc)
      implicit double precision (a-h,o-z)
c     split string s with sep.chars c
      include 'inc/p2_dim.inc'
      character c*(*)
      character s*(*)
      character tmp*(MXLNE)
      dimension idx1(m),idx2(m)
      call strfl(tmp,' ',MXLNE)
      tmp=s
      m0=0
      ifl=1
      do i=1,n
       ls=0
       do k=1,nc
        if(c(k:k).eq.tmp(i:i))then
         ls=1
        endif
       enddo
       if(ls.ne.0)then
        if( ifl.eq.0)then
         idx2(m0)=i-1
c
c         write (6,*) 'idx2 ',m0,idx2(m0),tmp(i:i),tmp(i-1:i-1)
c
        endif
        ifl=1
       else
        if(ifl.eq.1)then
         m0=m0+1
         idx1(m0)=i
c
c        write (6,*) 'idx1 ',m0,idx1(m0),tmp(i:i),tmp(i-1:i-1)
c
        endif
        ifl=0
       endif
      enddo
      mm=m0
      if(m0.gt.m) then
       write (6,*) 'ssplt ', m0,m
       stop  
      endif
      if(ifl.eq.0)then
       idx2(mm)=n
      endif
      end

      function idxf(c,string,n)
c     first occurence of c in string
      character string*(*)
      character c
      do  i=1,n
       if(c.eq.string(i:i)) then
        idxf=i
        return
       endif
      enddo
      idxf=0
      return
      end
      subroutine sspltq(s,idx1,idx2,n,m,mm,c,nc,ierr)
      implicit double precision (a-h,o-z)
c     split string s with sep.chars c
c     keep singly quoted strings together, regardless of whether
c     they contain one of the sep chars or not, strip the qoutes 
      include 'inc/p2_dim.inc'
      character c*(*)
      character s*(*)
      character tmp*(MXLNE)
      character sq
      data sq/''''/
      dimension idx1(m),idx2(m)
      ierr=0
      call strfl(tmp,' ',MXLNE)
      iquote=0
      tmp=s
      m0=0
      ifl=1
      do i=1,n
       ls=0
       if(tmp(i:i).eq.sq) then
        if( iquote.eq.0 ) then
         iquote=1
        else
         iquote=0
        endif
       endif
       if(iquote.eq.0) then
        do k=1,nc
         if(c(k:k).eq.tmp(i:i))then
          ls=1
         endif
        enddo
       endif
       if(ls.ne.0)then
        if( ifl.eq.0)then
         idx2(m0)=i-1
        endif
        ifl=1
       else
        if(ifl.eq.1)then
         m0=m0+1
         idx1(m0)=i
        endif
        ifl=0
       endif
      enddo
      mm=m0
      if(m0.gt.m) then
       write (6,*) 'ssplt ', m0,m
       stop  
      endif
      if(ifl.eq.0)then
       idx2(mm)=n
      endif
      if(iquote.ne.0) then
       ierr=2
      endif
      do i=1,mm
       if((tmp(idx1(i):idx1(i)).eq.sq).and.(tmp(idx2(i):idx2(i)).eq.sq))
     $      then
        idx1(i)=idx1(i)+1
        idx2(i)=idx2(i)-1
       endif
      enddo
      end

      function idxl(c,string,n)
c     last occurence of c in string
      character string*(*)
      character c
      do  i=n,1,-1
       if(c.eq.string(i:i)) then
        idxl=i
        return
       endif
      enddo
      idxl=0
      return
      end

      function idxfn(c,string,n)
c     first occurence of [^c] in string
      character string*(*)
      character c
      integer i
      integer idxfn
      integer n
      do  i=1,n
       if(c.ne.string(i:i)) then
        idxfn=i
        return
       endif
      enddo
      idxfn=0
      return
      end

      function idxln(c,string,n)
c     last occurence of [^c] in string
      character string*(*)
      character c
      do  i=n,1,-1
       if(c.ne.string(i:i)) then
        idxln=i
        return
       endif
      enddo
      idxln=0
      return
      end
      subroutine unprtct(lne,mm0,mm,idx1,idx2,MXF)
      implicit double precision (a-h,o-z)
c     unprotect: change \$ into $ 
c                change $$ into $
      character lne*(*)
      dimension idx1(MXF)
      dimension idx2(MXF)
      do i=mm0,mm
       if((idx2(i)-idx1(i)).gt.1) then
        if(lne(idx1(i):idx1(i)+1).eq.'\\$') then
         lne(idx1(i):idx1(i))=' '       
        elseif(lne(idx1(i):idx1(i)+1).eq.'$$') then
         lne(idx1(i):idx1(i))=' '
        endif
       endif
      enddo
      end 


      subroutine p2arry(iu,nm,iv,n1)
      implicit double precision (a-h,o-z)
CINPUT array  inp a 0 i=1
CINPUT vector inp a 1 i=1
CINPUT table  inp a 2 i=1
C    simple 1d array input
C    nm basename, n1 first index
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      character nm*(*)
      character vnm*64
      dimension idx1(MXF),idx2(MXF)
      idx=n1
      k0=idxln(' ',nm,64)
 1    continue 
      call gpnln2(iu,il,lne,MXLNE,idx1,idx2,mm,lbuf,'p2arry')
      if((mm.eq.1).and.(lne(idx1(1):idx2(1)).eq.'end')) goto 2
      if(iv.eq.0) then
       mm1=1
       mm2=mm
      elseif(iv.eq.1) then
       mm1=1
       mm2=1
       ifc=1
      elseif(iv.eq.2) then
       mm1=1
       mm2=1
       ifc=2
      endif
      do i=mm1,mm2
       if((iv.eq.0).or.(iv.eq.1)) then
        ix2=idx
       else
        ix2=ifroms(lne(idx1(1):idx2(1)),1)
       endif
       if(ix2.lt.10) then
        write(vnm,10) nm,ix2
       elseif(ix2.lt.100) then
        write(vnm,11) nm,ix2
       elseif(ix2.lt.1000) then
        write(vnm,12) nm,ix2
       elseif(ix2.lt.10000) then
        write(vnm,13) nm,ix2
       elseif(ix2.lt.100000) then
        write(vnm,14) nm,ix2
       else
        call fatal('p2arry: too many elements')
       endif
 10    format (a,'.',i1)
 11    format (a,'.',i2)
 12    format (a,'.',i3)
 13    format (a,'.',i4)
 14    format (a,'.',i5)
       k1=idxln(' ',vnm,64)
       if(iv.eq.0) then
        call p2setv(vnm(1:k1),lne(idx1(i):idx2(i)))
       else
        call p2setv(vnm(1:k1),lne(idx1(ifc):idx2(mm)))
       endif
       if(ip2deb.gt.0) then
        if(iv.eq.0) then
         write (6,100) nm,ix2,vnm(1:k1),lne(idx1(i):idx2(i))
        else
         write (6,100) nm,ix2,vnm(1:k1),lne(idx1(ifc):idx2(mm))
        endif
 100    format ('p2arry : ',a,' index=',i5,' nm=',a,'[',a,']')
       endif
       idx=idx+1
      enddo
      goto 1
 2    continue 
C    dimension stored in name.size
      id=idx-1
      write(vnm,20) nm
 20   format(a,'.size')
      k1=idxln(' ',vnm,64)
      call p2sti(id,vnm(1:k1))
      if(ip2deb.gt.0) then
       write (6,101) id,vnm(1:k1)
 101   format('p2arry : size = ',i6,' stored in ',a)
      endif
      end 
      subroutine p2setv(nm,val)
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      character nm*(*)
      character val*(*)
      character tmp*(MXLNE)
      call strfl(tmp,' ',MXLNE)
      tmp=val
      if(tmp(1:1).eq.'"') then
       jx1=2
       jx2=idxln('"',tmp,MXLNE)
      else
       jx1=1
       jx2=idxln(' ',tmp,MXLNE)
      endif
      do k=1,nvtop
       if(varnm(k).eq.nm)then
        varval(k)=tmp(jx1:jx2)
        ivrln(k)=jx2-jx1+1
        goto 99
       endif
      enddo
      nvtop=nvtop+1
      if(nvtop.gt.MXVAR)then
       stop 'P2: variable stack full'
      else
       varnm(nvtop)=nm
       varval(nvtop)=tmp(jx1:jx2)
       ivrln(nvtop)=jx2-jx1+1
       k=nvtop
      endif
 99   continue 
CTEST>
c$$$      write (6,*) 'len ',jx1,jx2 
c$$$      write (6,10) k,nvtop,varnm(k),varval(k)(1:ivrln(k)),ivrln(k)
c$$$ 10   format ('p2setv ',i4,i4,' [',a,']  [',a,'] ',i4)
CTEST< 
      end 

      subroutine p2gtid(nm,iv,idef,iflag)
c     get value of an optional variable
c     if variable is not set, assign idef instead
c     if flag is set, the variable is set
c     (this means that subsequent reads will get *this* default value)
c
c     (interger version)
      implicit double precision (a-h,o-z)
      common /p2gt01/ iwp2gtd
      character nm*(*)
      ierr = 0
      call p2gti(nm,iv,ierr)
      if(ierr.ne.0) then
       iv = idef
       if(iflag.ne.0) then
        call p2sti(iv,nm)
       endif
       if(iwp2gtd.eq.1) then
        write (6,10) iv,nm
 10     format('p2gtid: value= ',i12,' var=',a)
        call warn('p2gtid: default value used.')
       elseif(iwp2gtd.eq.2) then
        write (6,10) iv,nm
        call fatal('p2gtid: default value used.')
       endif
      endif
      end 
      subroutine p2gtfd(nm,v,def,iflag)
c     get value of an optional variable
c     if variable is not set, assign def instead
      implicit double precision (a-h,o-z)
      common /p2gt01/ iwp2gtd
      character nm*(*)
      ierr = 0
      call p2gtf(nm,v,ierr)
      if(ierr.ne.0) then
       v = def
       if(iflag.ne.0) then
        call  p2stf(v,nm)
       endif
       if(iwp2gtd.eq.1) then
        write (6,10) v,nm
 10     format('p2gtfd: value= ',e22.12,' var=',a)
        call warn('p2gtfd: default value used.')
       elseif(iwp2gtd.eq.2) then
        write (6,10) v,nm
        call fatal('p2gtfd: default value used.')
       endif
      endif
      end 

      subroutine p2gtad(nm,v,def,iflag)
c     get value of an optional variable
c     if variable is not set, assign def instead
      implicit double precision (a-h,o-z)
      common /p2gt01/ iwp2gtd
      character nm*(*)
      character v*(*)
      character def*(*)
      ierr = 0
      call p2gta(nm,v,ierr)
      if(ierr.ne.0) then
       v = def
       if(iflag.ne.0) then
        call  p2sta(v,nm)
       endif
       if(iwp2gtd.eq.1) then
        write (6,10) v,nm
 10     format('p2gtfd: value= ',a,' var=',a)
        call warn('p2gtfd: default value used.')
       elseif(iwp2gtd.eq.2) then
        write (6,10) v,nm
        call fatal('p2gtfd: default value used.')
       endif
      endif
      end 


      subroutine p2gti(nm,iv,ierr)
      implicit double precision (a-h,o-z)
C    get integer variable 
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      character nm*(*)
      character val*(MXLNE)
      character tmp*(MXLNE)
      write(tmp,20) nm
      n=idxln(' ',tmp,MXIDL)
      call p2var(tmp,n,val,k3,ie)
      if(ie.eq.0) then
       iv=ifroms(val(1:k3),1)
      else
       if(ierr.ne.0) then
        write (6,10) tmp(1:n)
        call fatal('p2gti: variable not fond')
       endif
      endif
      if(ierr.eq.0) then
       ierr=ie
      endif
 10   format('p2gti: variable name = ',a)
 20   format(a,' ')
      end 
      subroutine p2gtf(nm,v,ierr)
      implicit double precision (a-h,o-z)
C    get fload variable 
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      character nm*(*)
      character val*(MXLNE)
      character tmp*(MXLNE) 
      write(tmp,20) nm
      n=idxln(' ',tmp,MXIDL)
      call p2var(tmp,n,val,k3,ie)
      if(ie.eq.0) then
       v=dfroms(val(1:k3),1)
      else
       if(ierr.ne.0) then
        write (6,10) tmp(1:n)
        call fatal('p2gtf: variable not fond')
       endif
      endif
      if(ierr.eq.0) then
       ierr=ie
      endif
 10   format('p2gtf: variable name = ',a)
 20   format(a,' ')
      end 


      subroutine p2gta(nm,v,ierr)
      implicit double precision (a-h,o-z)
C    get fload variable 
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      character nm*(*)
      character v*(*)
      character val*(MXLNE)
      character tmp*(MXLNE) 
      write(tmp,20) nm
      n=idxln(' ',tmp,MXIDL)
      call p2var(tmp,n,val,k3,ie)
      if(ie.eq.0) then
       v=val(1:k3)
      else
       if(ierr.ne.0) then
        write (6,10) tmp(1:n)
        call fatal('p2gtf: variable not fond')
       endif
      endif
      if(ierr.eq.0) then
       ierr=ie
      endif
 10   format('p2gtf: variable name = ',a)
 20   format(a,' ')
      end 

      subroutine p2var(nm,n,val,k3,ie)
      implicit double precision (a-h,o-z)
      character nm*(*)
      character val*(*)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      k1=n
      do k=1,nvtop
       k0=idxln(' ',varnm(k),MXIDL)
       if(k0.eq.k1) then
        if(varnm(k)(1:k0).eq.nm(1:k0))then
         val=varval(k)
         k3=ivrln(k)
         ie=0
         goto 99
        endif
       endif
      enddo
      ie=1
 99   continue 
CTEST> 
      if(ip2deb.ne.0) then
       if(ie.eq.0) then
        write (6,10) ie,k,nvtop,k3,varnm(k),varval(k)(1:k3)
       else
        write (6,10) ie,k,nvtop,k3,varnm(k),'<not found>'
       endif
 10    format ('p2var ',4i4,' [',a,']  [',a,']')
      endif
CTEST< 
      end 

      subroutine p2vin(fn,ifmt)
CINPUT printmacros a=stdout 0
CINPUT savemacros  a=stdout 1
C list all defined macros
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      character tmp1*(MXIDL)
      character tmp2*(MXLNE)
      character fn*(*)
      if(ifmt.ne.0) then
       call file(iu,fn,'app',1,0)
      endif
      mxw1=6
      mxw2=7
      do i=1,nvtop
       idx=idxln(' ',varnm(i),MXIDL)
       if(idx.gt.mxw1) mxw1=idx
       if(ivrln(i).gt.mxw2) mxw2=ivrln(i)
      enddo
      if(ifmt.eq.0) then
       write (6,9)
       write (6,10) nvtop,MXVAR
       call strfl(tmp1,' ',mxw1)
       call strfl(tmp2,' ',mxw2)
       tmp1(1:6)=' NAME '
       tmp2(1:7)=' VALUE '
       write (6,11) tmp1(1:mxw1),tmp2(1:mxw2)
       write (6,13)
      endif
      do i=1,nvtop
       call strfl(tmp1,' ',mxw1)
       call strfl(tmp2,' ',mxw2)
       idx=idxln(' ',varnm(i),MXIDL)
       tmp1(1:idx)=varnm(i)(1:idx)
       tmp2(1:ivrln(i))=varval(i)(1:ivrln(i))
       if(ifmt.eq.0) then
        write (6,11) tmp1(1:mxw1),tmp2(1:mxw2)
       else
        if(tmp1(1:1).ne.'*') then
         write(iu,15) tmp1(1:mxw1),tmp2(1:mxw2)
        endif
       endif
      enddo
      if(ifmt.eq.0) then
       write (6,13)
       write (6,12)
       write (6,13)
      endif
c NOTE: COMMENTED -> PROBLEMS ON SGI MACHINE WITH f90
c     if(ifmt.ne.0) then
c      call flush(iu)
c     endif
 15   format('def ',a,' ',a)
 9    format('****  LIST OF MACROS ****')
 10   format(i8,' of ',i8,' used')
 11   format(a,' ',a)
 12   format('*************************')
 13   format(' ')
      end 

      subroutine p2ind(s1,s2,n)
      implicit double precision (a-h,o-z)
C    copies s1 to s2 and expands one-chracter macros
C    (index-macros)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      character s1*(*)
      character s2*(*)
      character tmp*(MXIDL)
      character vn
      ii=0
      i=1
 1    continue 
      if(s1(i:i).ne.'$') then
       ii=ii+1
       s2(ii:ii)=s1(i:i)
       i=i+1
      else
       if(i.eq.n) goto 1000
       if(s1(i+1:i+1).eq.'(') goto 1000
       vn=s1(i+1:i+1)
       call p2var(vn,1,tmp,k3,ie)
       if(ie.ne.0) goto 1200
       do k=1,k3
        if(tmp(k:k).eq.'$') then
         goto 1300
        else
         ii=ii+1
         s2(ii:ii)=tmp(k:k)
        endif
       enddo
       i=i+2
      endif
      if(i.le.n) goto 1
      if(ip2deb.ne.0) then
       write (6,100) s1(1:n),s2(1:ii),n,ii
 100   format ('p2ind s1=[',a,'] s2=[',a,'] n=',i3,' ii=',i3)
      endif
      n=ii
      return
 1000 continue 
      call fatal('p2ind : badly formed index name')
 1200 continue 
      write (6,110) vn
 110  format('[',a,']')
      call fatal('p2ind : macro not found')
 1300 continue 
      call fatal('p2ind : recursive index macro')  
      end 

      subroutine p2stf(f,nam)
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2b.inc'
C    put float to variable list
      character tmpstr*64
      character tmp2*64
      character nmstr*25
      character nam*(*)
      tmp2=' '
      write(tmpstr,10) f
 10   format(e22.12) 
      jtmp=0
      do j=1,64
       if(tmpstr(j:j).ne.' ') then
        jtmp=jtmp+1
        tmp2(jtmp:jtmp)=tmpstr(j:j)
       endif
      enddo
      if(nam(1:1).eq.'*') then
       nans=nans+1
       if(nans.lt.10)then
        write(nmstr,12) nans
       elseif(nans.lt.100)then
        write(nmstr,13) nans
       else
        write(nmstr,14) nans
       endif
      else
       nmstr=nam
      endif
 12   format ('ANS',i1)
 13   format ('ANS',i2)
 14   format ('ANS',i3)  
      if(ip2deb.ne.0) then
       call p2setv(nmstr,tmpstr)
       k1=idxln(' ',nmstr,25)
       write (6,20) nmstr(1:k1),tmpstr
      endif
 20   format ('** ',a,' set to ',a)
      end 

      subroutine p2sti(i,nam)
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2b.inc'
C    put integer to variable list
      character nam*(*)
      character tmpstr*64
      character tmp2*64
      character nmstr*25
      tmp2=' '
      write(tmpstr,10) i
 10   format(i60)
      jtmp=0
      do j=1,64
       if(tmpstr(j:j).ne.' ') then
        jtmp=jtmp+1
        tmp2(jtmp:jtmp)=tmpstr(j:j)
       endif
      enddo
      if(nam(1:1).eq.'*') then
       nans=nans+1
       if(nans.lt.10)then
        write(nmstr,12) nans
       elseif(nans.lt.100)then
        write(nmstr,13) nans
       else
        write(nmstr,14) nans
       endif
      else
       nmstr=nam
      endif
 12   format ('ANS',i1)
 13   format ('ANS',i2)
 14   format ('ANS',i3)
      call p2setv(nmstr,tmp2) 
      if(ip2deb.ne.0) then
       k1=idxln(' ',nmstr,25)
       write (6,20) nmstr(1:k1),tmp2(1:jtmp)
 20   format ('** ',a,' set to ',a)
      endif
      end 

      subroutine p2isvar(lne,idx1,idx2,n,i,astr,k3,ie)
      implicit double precision (a-h,o-z)
      character lne*(*)
      character astr*(*)
      dimension idx1(n),idx2(n)
      if(lne(idx1(i):idx1(i)).eq.'$') then
       ll0=idx2(i)-idx1(i)
       call p2var(lne(idx1(i)+1:idx2(i)),ll0,astr,k3,ie)
      else
       astr=lne(idx1(i):idx2(i))
       k3=idx2(i)-idx1(i)+1
       ie=0
      endif
      end 
      subroutine p2sta(a,nam)
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2b.inc'
C    put integer to variable list
      character a*(64)
      character nam*(64)
      character tmpstr*64
      character tmp2*64
      character nmstr*25
      tmp2=' '
      write(tmpstr,10) a
 10   format(a)
      jtmp=0
      do j=1,64
       if(tmpstr(j:j).ne.' ') then
        jtmp=jtmp+1
        tmp2(jtmp:jtmp)=tmpstr(j:j)
       endif
      enddo
      if(nam(1:1).eq.'*') then
       nans=nans+1
       if(nans.lt.10)then
        write(nmstr,12) nans
       elseif(nans.lt.100)then
        write(nmstr,13) nans
       else
        write(nmstr,14) nans
       endif
      else
       nmstr=nam
      endif
 12   format ('ANS',i1)
 13   format ('ANS',i2)
 14   format ('ANS',i3)
      call p2setv(nmstr,tmp2) 
      if(ip2deb.ne.0) then
       k1=idxln(' ',nmstr,25)
       write (6,20) nmstr(1:k1),tmp2(1:jtmp)
 20   format ('** ',a,' set to ',a)
      endif
      end 

      subroutine macex(lbuf1,lbuf2,lbln,idx1,idx2,n,mm,irec,il)
      implicit double precision (a-h,o-z)
c     macro expand from lbuf1 to lbuf2
c     irec=1 allow recursion
      include 'inc/p2_dim.inc'
      include 'inc/p2b.inc'
      character lbuf1*(*)
      character lbuf2*(*)
      character vn*(MXIDL),vn2*(MXIDL)
      character nc
      dimension idx1(n),idx2(n)
      character quote
      quote='"'
      irr=0
      mxrec=100
 1    continue 
      irr=irr+1
      if(irr.gt.mxrec) goto 1200
c
      ib2=1
      k0=0
c     do we have macros ?
      do i=1,mm
       do j=idx1(i),idx2(i)-1
        if(lbuf1(j:j).eq.'$') then
         k0=1
         goto 2
        endif
       enddo
      enddo
      if(k0.eq.0) return
c
 2    continue 
      ixp=0
      do i=1,mm
       j1=idx1(i)
       j2=idx2(i)
       j=j1
       imac=0
       vn=' '
       ivn=0
       iap=1
 3     continue 
       nc=lbuf1(j:j)
       if((nc.eq.quote).and.(lbuf1(j2:j2).eq.quote)) then
        lbuf2(ib2:)=lbuf1(j1:j2)
        ib2=ib2+j2-j1+1
        goto 8
       endif
       if((imac.eq.0).and.(nc.eq.'$')) then
        if(lbuf1(j+1:j+1).eq.'$') then
         lbuf2(ib2:)='$$'
         ib2=ib2+2
         j=j+1
         iap=0
         goto 99
        elseif(lbuf1(j+1:j+1).eq.'(') then
         imac=1
         j=j+1
         goto 99
        else
         imac=2
         ivn=1
         j=j+1
         vn=lbuf1(j:j)
         goto 99
        endif
       else
        if(imac.eq.1) then
         if(nc.eq.')') then          
          imac=2
          goto 99
         else
          ivn=ivn+1
          vn(ivn:ivn)=nc
         endif
        endif
       endif
 99    continue 
       if(imac.eq.0) then
        if(iap.eq.1) then
         lbuf2(ib2:)=nc
         ib2=ib2+1
        endif
        iap=1
       elseif(imac.eq.2) then
C    index inside vn ?? 
C    (macro names can contain one-letter macros like $(foo.$i.$k.bar)) 
        k01=ivn
        call p2ind(vn,vn2,k01)
C
        call p2var(vn2,k01,lbuf2(ib2:),k3,ie)
        if(ie.ne.0) goto 1000
        ib2=ib2+k3
        ixp=ixp+1
        ivn=0
        vn=' '
        imac=0
       endif
       j=j+1
       if(j.le.j2) goto 3
 8     continue 
       lbuf2(ib2:ib2)=' '
       ib2=ib2+1
       if(imac.ne.0) goto 1300
      enddo
c     copy back & resplit
      lbuf1=lbuf2(1:ib2)
      call ssplt(lbuf1,idx1,idx2,lbln,n,mm,fsep,ifsep)
      if((irec.eq.1).and.(ixp.ne.0)) goto 1
      return
 1000 continue 
      write (6,5) il,vn2(1:k01)
 5    format('> line ',i4,' macro: [',a,']')
      call fatal('macex: macro not found')
 1100 continue 
      call fatal('macex: buffer overflow')
 1200 continue 
      call fatal('macex: too many recursions')
 1300 continue 
      call fatal('macex: badly formed macro')
      end 

C------------------------------------------------------------------------
      subroutine echo(iiu,iou)
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
c
      lnum=0
      call rdln(iiu,lne,lnum,ieof)
c$$$      if(lne(1:1).eq.'&') then
c$$$       iimode=1
c$$$       write (6,*) '(interactive mode)'
c$$$       call noxxx(1)
c$$$       return
c$$$      else
       call noxxx(0)
c$$$      endif
      rewind(iiu)
      call setpos(iiu,0)
      write(iou,10)
 10   format ('-------------------------------------------------')
 1    continue 
      lnum=0
      call rdln(iiu,lne,lnum,ieof)
      if(ieof.eq.0) then
       k1=idxln(' ',lne,MXLNE)
       write(iou,'(a)') lne(1:k1)
       goto 1
      else
       write(iou,10)
       rewind(iiu)
       call setpos(iiu,0)
      endif
      end 

      subroutine p2ini0
c     called by p2init
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc' 
      include 'inc/p2b.inc' 
      character tab
      data tab /'\t'/
C    commands in output
      inice=0
C    abbreviation of commands allowed
      ip2ab=1
C    filed separator characters
      ifsep=3
      fsep=' ,'//tab
      nvtop=0
      nans=0
      iskip=0
C    one line comment 
      cmtchr='#'
C    'anticomment'
      antchr='#'
      iant=0
C debugging off
      ip2deb=0
C    nest level 0
      lnest=0
C    init file own handling 
      call setinp(5)
      call file(iutmp,'<input>','old',1,0)
      end 

      subroutine p2ddfl(ikw,iarg,ini,ifl,iaa,nargs,fmt,
     $     itmp,ftmp,atmp,is1,is2,
     $     idefpp,idefvv,ddefvv,adefvv,
     $     ierr)
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      character fmt*(MXIDL)
c$      
      character atmp*(MXLNE)
      dimension itmp(MXF)
      dimension ftmp(MXF)
      dimension is1(MXF)
      dimension is2(MXF)
c$
      dimension idefpp(MXDEF,MXKEY)
      dimension idefvv(MXDEF)
      dimension ddefvv(MXDEF)
      character adefvv(MXDEF)*(MXLNE)
C
      ierr=0
      do i=iarg+1,nargs
       ix=idefpp(i,ikw)
       if(ix.eq.0) then
        write (6,*) '> ',i,ikw
        call fatal('p2defv: invalid index')
       endif
       if(fmt(i:i).eq.'i') then
        ini=ini+1
        itmp(ini)=idefvv(ix)
       elseif(fmt(i:i).eq.'d') then
        ifl=ifl+1
        ftmp(ifl)=ddefvv(ix)
       elseif(fmt(i:i).eq.'a') then
        iaa=iaa+1
        ll=idxln(' ',adefvv(ix),MXLNE)
        if(iaa.gt.1)then
         is1(iaa)=is2(iaa-1)+1
         is2(iaa)=is1(iaa)+ll-1
        else
         is1(iaa)=1
         is2(iaa)=ll
        endif
        atmp(is1(iaa):is2(iaa))=adefvv(ix)(1:ll)
       else
        err=5
       endif
      enddo
      end 

      subroutine p2prag(iu,str,itest)
C pragmas
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      include 'inc/p2defv.inc'
      common /p2gt01/ iwp2gtd
c
      character str*(*)
      if(str(1:1).ne.'.') then
       itest=0
       return
      else
       itest=1
       if(str(2:8).eq.'VERSION') then
        call pver
C allow commands to be abbreviated
       elseif(str(2:7).eq.'ABBREV') then
        ip2ab=1
       elseif(str(2:9).eq.'NOABBREV') then
        ip2ab=0
C print command and separator line
       elseif(str(2:4).eq.'SEP') then
        inice=1
       elseif(str(2:6).eq.'NOSEP') then
        inice=0
C allow use of default arguments 
       elseif(str(2:11).eq.'NODEFAULTS') then
        ip2dfl=0
       elseif(str(2:9).eq.'DEFAULTS') then
        ip2dfl=1
C no warning if default values of variables are used
       elseif(str(2:14).eq.'VARDEF_SILENT') then
        iwp2gtd=0
C  warning if default values of variables are used
       elseif(str(2:12).eq.'VARDEF_WARN') then
        iwp2gtd=1
C    Error if default values of variables are used
       elseif(str(2:13).eq.'VARDEF_ERROR') then
        iwp2gtd=2
C parser debugging 
       elseif(str(2:6).eq.'DEBUG') then
        ip2deb=1
       elseif(str(2:8).eq.'NODEBUG') then
        ip2deb=0
C change character for one-line-comments 
       elseif(str(2:9).eq.'COMMENTS') then
        cmtchr=str(10:10)
C set character with which lines have to begin 
C in order to be parsed  
       elseif(str(2:5).eq.'ONLY') then
        iant=1
        antchr=str(6:6)
C reset this (all lines which are not comments will be parsed)
       elseif(str(2:4).eq.'ALL') then
        iant=0
C set input line number
       elseif(str(2:5).eq.'LINE') then
        ipos=ifroms(str(6:),1)
        call setpos(iu,ipos)
        if(ip2deb.eq.1) then
         write (6,11) iu,ipos
        endif
       else
        write (6,10) str
       endif
      endif
      if(ip2deb.eq.1) then
       call incpos(iu,ii,0)
       write (6,12) ii,str
      endif
 10   format('p2: unnkown option ignored ',a)
 11   format('p2: line number set to ',i6,' (unit=',i4,')')
 12   format('### line ',i6,' PRAGMA ',a) 
      end 

C## simple control rout. (not included by default) ## 
      subroutine skipto(lbl,icond)
      implicit double precision (a-h,o-z)
CINPUT skipto a i=1
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      character lbl*(*)
      if(lbl.eq.'off') then
       iskip=0
       write (6,*) 'skip switched off'
      else
       if(icond.ne.0) then
        skipch=lbl
        iskip=1
        write (6,5) 'skipping to ',skipch
       endif
      endif
 5    format(a,a)
      end 

      subroutine gotol(nl)
CINPUT gotol i
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      if(nl.gt.il) then
       do i=1,nl-il-1
        call rdln(iutmp,lne,il,ieof)
        if(ieof.eq.1) goto 99
       enddo
      elseif(nl.lt.il) then
       do i=1,il-nl+1
        backspace(iutmp,err=99)
        call incpos(iutmp,il,-1)
       enddo
      endif
      return
 99   continue 
      write (6,*) '> ',il,nl
      call fatal('gotol: line out of reach')
      end 

      subroutine loop(nam,ia,ib,is)
CINPUT loop a i i i=1
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      character nam*(*)
      lctrl=lctrl+1
ctest>
c$$$      write (6,*) 'lctrl ',lctrl
ctest<
      if(lctrl.gt.MXCTRL) then
       write (6,*) '> ',lctrl,MXCTRL
       call fatal('control stack full')
      endif
      kctrl(lctrl)=il+1
      ktctrl(lctrl)=1
      ipctl(1,lctrl)=ia
      ipctl(2,lctrl)=ib
      ipctl(3,lctrl)=is
      apctl(1,lctrl)=nam
      call p2sti(ia,nam)
      if(ip2deb.gt.0) then
       write (6,10) il,lctrl,kctrl(lctrl),ktctrl(lctrl),ia,ib,is,nam
 10    format('LOOP ENRTY : ',7(i6,' '),' ',a)
      endif
      end 
      subroutine ectrl2
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      if(lctrl.ne.0) then
       write (6,*) '> line',il
       call fatal('unclosed block at end of input')
      endif
      end 
      subroutine ectrl
CINPUT end
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      if(lctrl.eq.0) then
       write (6,*) '> line',il
       call fatal('unexpected end')
      else
ctest>
       if(ip2deb.gt.0) then
        write (6,10) il,lctrl,ktctrl(lctrl)
 10     format('END ',3i6)
       endif
ctest<
       goto(1,999) ktctrl(lctrl)
       goto 999
 1     continue 
       call eloop
       goto 99
 99    continue 
       return
 999   continue 
       call fatal('error in ectrl')
      endif
      end 
      subroutine eloop
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      character nam*(MXIDL)
      i=ipctl(1,lctrl)
      ib=ipctl(2,lctrl)
      is=ipctl(3,lctrl)
      nam= apctl(1,lctrl)
      i=i+is
      ipctl(1,lctrl)=i
ctest>
      if(ip2deb.gt.0) then
       write (6,10) il,lctrl,kctrl(lctrl),i,ib,is,nam
 10    format('LOOP END : ',6(i6,' '),a)
      endif
ctest<
      if(i.le.ib) then
       call p2sti(i,nam)
       call gotol(kctrl(lctrl))
      else
       ktctrl(lctrl)=0
       lctrl=lctrl-1
       if(lctrl.lt.0) then
        call fatal('endloop: loop missing')
       endif
      endif
      end 
C     prints list of commands
      subroutine cmdlst(icol)
CINPUT ?  3
CINPUT ?? 0
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      write (6,*) '------ commands ------'
      if(icol.gt.0) then
       nl=nkey/3
       m11=nkey-3*nl
       id=1
       do l=1,nl
        write (6,10) keys(id),keys(id+1),keys(id+2)
        id=id+3
       enddo
       if(m11.eq.1) then
        write (6,11) keys(nkey)
       endif
       if(m11.eq.2) then
        write (6,12) keys(nkey-1),keys(nkey)
       endif
      else
       do i=1,nkey
        write (6,20) keys(i),nargs(i),fmts(i)
       enddo
      endif
 10   format(a16,'  ',a16,'  ',a16)
 11   format(a16)
 12   format(a16,'  ',a16)
 20   format(a32,' ',i3,' ',a32)
      end 
CCC load command file ...
      subroutine ldcmdf(fn)
CINPUT load a
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      character fn*(*)
      character fnsave*256
      logical ex
      fnsave=fn
      ifn=LEN(fn)
      if(ip2deb.ne.0) then
       write (6,*) '========== LOAD (ldcmdf) ========== '
       write (6,*) '(''file ['',a,'']'')',fn
       write (6,*) 'Len = ',ifn
      endif
      inquire(file=fnsave(1:ifn),exist=ex)
      if(ex) then
       call file(iu,fnsave(1:ifn),'old',1,0)
       rewind(iu)
       call setpos(iu,0)
       lnest=lnest+1
       if(lnest.gt.MXNEST) then
        call fatal('too many levels for nested load')
       endif
       write (6,10) fnsave(1:ifn)
C    push flags & set defaults
       iio(lnest)=iimode
       co(lnest)=cmtchr
       ianto(lnest)=iant
       iimode=0
       cmtchr='#'
       iant=0
       lio(lnest)=il
       iuu(lnest)=iutmp
       call p2go(iu,0)
       iimode=iio(lnest)
       iant=ianto(lnest)
       cmtchr=co(lnest)
       il=lio(lnest)
       iutmp=iuu(lnest)
       call setinp(iutmp)
       lnest=lnest-1
       call file(iu,fnsave(1:ifn),'cls',1,0)
      else
       write (6,11) fnsave(1:ifn)
       if(iimode.eq.0) then
        call fatal('ldcmdf: cannot load command file')
       endif
      endif
 10   format('executing [',a,'] ... ')
 11   format('file ',a,' not found.')
      end 
      subroutine p2mode(iant0,cmtch0,antch0)
      implicit double precision (a-h,o-z)
      character cmtch0*(*),antch0*(*)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      if(cmtch0(1:4).eq.'hash') then
       cmtchr='#'
      else
       cmtchr=cmtch0(1:1)
      endif
      antchr=antch0(1:1)
      iant=iant0
      end 
      subroutine p2fop(vn,aop,d1,op,d2)
CINPUT fop a a d a d=0
CINPUT @@  a a d a d=0
C    simple operations on integer variables
      implicit double precision (a-h,o-z)
      character vn*(*),op*(*),aop*(*)
      if(aop.eq.'=') then
       if(op.eq.'+') then
        dtmp=d1+d2
       elseif(op.eq.'-') then
        dtmp=d1-d2
       elseif(op.eq.'*') then
        dtmp=d1*d2
       elseif(op.eq.'/') then
        dtmp=d1/d2
       elseif(op.eq.'**') then
        dtmp=d1**d2
       elseif(op.eq.'exp') then
        dtmp=exp(d1)
       else
        call fatal('p2fop: unknown operation')
       endif
       call p2stf(dtmp,vn)
      else
       call fatal('p2iop: unknown assignment op')
      endif
      end       
      subroutine p2iop(vn,aop,i1,op,i2)
CINPUT iop a a i=0 a=x i=0
CINPUT @   a a i=0 a=x i=0
C    simple operations on integer variables
      implicit double precision (a-h,o-z)
      character vn*(*),op*(*),aop*(*)
      call intro('p2iop')
      ierr=0
      iv=0
      call p2gti(vn,iv,ierr)
      if(aop.eq.'=') then
       if(op.eq.'+') then
        itmp=i1+i2
       elseif(op.eq.'-') then
        itmp=i1-i2
       elseif(op.eq.'*') then
        itmp=i1*i2
       elseif(op.eq.'/') then
        itmp=i1/i2
       elseif(op.eq.'mod') then
        itmp=mod(i1,i2)
       elseif(op.eq.'++') then
        itmp=i1+1
       elseif(op.eq.'--') then
        itmp=i1-1
       elseif(op.eq.'x') then
        itmp=i1
       else
        call fatal('p2iop: unknown operation')
       endif
       iv=itmp
      elseif(aop.eq.'++') then
       if(ierr.ne.0) goto 999
       iv=iv+1
      elseif(aop.eq.'--') then
       if(ierr.ne.0) goto 999
       iv=iv-1
      else
       call fatal('p2iop: unknown assignment op')
      endif
      call p2sti(iv,vn) 
      call outro('p2iop')
      return
 999  continue 
      write (6,10) vn
 10   format('Name= ',a)
      call fatal('p2iop: undefined variable')
      end       

      subroutine strfl(s,c,n)
C    string fill
      implicit double precision (a-h,o-z)
      character s*(*)
      character c
      do i=1,n
       s(i:i)=c
      enddo
      end 
cccc############################################################################

C **********************************************************************
      subroutine p2gtfarr(iu,a,np,n,iflag)
C **********************************************************************
c     read line input into array
c     iflag == 1 : warn   if number of elements does not match
c     iflag == 2 : error  if number of elements does not match
c     iu : input unit
c     a  : array to be filled
c     np : physical dimension of a
c     n  : number of elements to fill
c     
C **********************************************************************
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      dimension idx1(MXF),idx2(MXF)
      dimension a(np)
c
      call gpnln2(iu,il,lne,MXLNE,idx1,idx2,mm,lbuf,'p2gtfarr')
      nr=min(n,mm)
c
      if(iflag.ne.0) then
       if(n.ne.mm) then
        write (6,9)
        write (6,10) n,mm
 9      format('p2gtfarr: number of elements does not match')
 10     format('need ',i4,' have ',i4)
        if(iflag.eq.1) then
         call warn('mismatch in p2gtfarr')
        else if(iflag.eq.2) then
         call fatal('mismatch in p2gtfarr')
        endif
       endif
      endif
c 
      do i=1,nr
       a(i) = dfroms(lne(idx1(i):idx2(i)), 0)
      enddo
      end 
C **********************************************************************
      subroutine p2gtiarr(iu,ia,np,n,iflag)
C **********************************************************************
c     read line input into array ** INTEGER VERSION **
c     iflag == 1 : warn   if number of elements does not match
c     iflag == 2 : error  if number of elements does not match
c     iu : input unit
c     ia : array to be filled
c     np : physical dimension of a
c     n  : number of elements to fill
c     
C **********************************************************************
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      dimension idx1(MXF),idx2(MXF)
      dimension ia(np)
c
      call gpnln2(iu,il,lne,MXLNE,idx1,idx2,mm,lbuf,'p2gtiarr')
      nr=min(n,mm)
c
      if(iflag.ne.0) then
       if(n.ne.mm) then
        write (6,9)
        write (6,10) n,mm
 9      format('p2gtiarr: number of elements does not match')
 10     format('need ',i4,' have ',i4)
        if(iflag.eq.1) then
         call warn('mismatch in p2gtiarr')
        else if(iflag.eq.2) then
         call fatal('mismatch in p2gtiarr')
        endif
       endif
      endif
c 
      do i=1,nr
       ia(i) = ifroms(lne(idx1(i):idx2(i)), 0)
      enddo
      end 

      subroutine set_imode(i)
      implicit double precision (a-h,o-z) 
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2b.inc'
      iimode=i
      if(i.eq.0) then
       call noxxx(0)
      else
       call noxxx(1)
      endif
      end 

      subroutine stripquotes(s)
      implicit double precision (a-h,o-z)
      character s*(*)
      character tmp*256
      character sq
      data sq/''''/
      n=LEN(s)
      if(n.gt.256) return 
      ib=idxln(' ',s,n)
      if(s(1:1).eq.sq) then 
       ia=2
      else
       ia=1
      endif
      if(s(ib:ib).eq.sq) then
       ib=ib-1
      endif
      j=0
      do i=ia,ib
       j=j+1
       tmp(j:j)=s(i:i)
      enddo
      s=tmp(1:j)
      end 

      subroutine append_number(s,n,r,nr,idot)
c     append number n to string s (result in r)
      implicit double precision (a-h,o-z)
      character s*(*)
      character r*(*)
      character fmt*32
      lnum=1
      ntmp=n
 1    continue
      ntmp=ntmp/10
      if(ntmp.ne.0) then
       lnum=lnum+1
       goto 1
      endif

      l=LEN(s)

      if(idot.eq.0) then
       write(fmt,100) lnum
       l2=0
      else
       write(fmt,101) lnum
       l2=1
      endif
 100  format('(a,i',i1,')')
 101  format('(a,''.'',i',i1,')')
      write(r,fmt) s(1:l),n
      nr=l+lnum+l2
      end

      block data p2def_data

      implicit double precision (a-h,o-z)

      character p2prmt*16
      common /p26/p2prmt
      common /p30/iimode,idum008
      common /p2gt01/ iwp2gtd

      data p2prmt /'> '/
      data iimode /0/
      data iwp2gtd /0/
      end
