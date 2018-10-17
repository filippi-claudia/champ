      implicit real*8(a-h,o-z)
c removes duplicate lines from walkalize files created by runs that
c did not complete.  ... Cyrus Umrigar (Sep. 2001)
      character*60 string
      character*13 filename
      character*17 filename2

      read(5,*) nproc
      do 60 idtask=0,nproc-1
        if(idtask.le.9) then
          write(filename,'(''walkalize.'',i1)') idtask
          write(filename2,'(''walkalize_new.'',i1)') idtask
         elseif(idtask.le.99) then
          write(filename,'(''walkalize.'',i2)') idtask
          write(filename2,'(''walkalize_new.'',i2)') idtask
         elseif(idtask.le.999) then
          write(filename,'(''walkalize.'',i3)') idtask
          write(filename2,'(''walkalize_new.'',i3)') idtask
         else
          stop 'idtask > 999'
        endif
        open(unit=11,file=filename,status='old')
        open(unit=12,file=filename2)

        nline=0
        do 10 i=1,2000000000
          read(11,fmt=*,end=20)
   10     nline=nline+1
   20   write(6,'(''nline='',i6)') nline
        rewind 11
        read(11,'(a60)') string
        write(12,*) string

        ipass_old=0
        do 30 iline=2,nline-1
          read(11,*) ipass,ffn,w,e,nwalk
          if(ipass.gt.ipass_old) then
            write(12,'(i8,f9.6,f12.5,f11.6,i5)') ipass,ffn,
     &      w,e,nwalk
            ipass_old=ipass
          endif
  30    continue

        read(11,*,err=40) nstep,nblk,nconf,etrial,tau,taueff
        write(12,'(3i5,f11.5,f7.4,f10.7,
     &  '' nstep,nblk,nconf,etrial,tau,taueff'')')
     &  nstep,nblk,nconf,etrial,tau,taueff
        goto 60
  40    backspace 11
        read(11,*) ipass,ffn,w,e,nwalk
        write(12,'(i8,f9.6,f12.5,f11.6,i5)') ipass,ffn,
     &  w,e,nwalk
  60  continue

      stop
      end
