      implicit real*8(a-h,o-z)
c Written by Cyrus Umrigar
c Calculates the mixed and growth energies (using data in fort.11
c created in a dmc run) for various f projection times.  Instead
c of using a certain number of f's it uses all f's but raising them
c to a power, less than one, so that the contributions from the
c ones further back tend to 1.  This gives somewhat smoother curves
c than the old way.       Cyrus Umrigar (Mar 1992, Sep. 2001)
      parameter (MDATA=204000,MPROC=100)
      dimension f(MDATA,MPROC),w(MDATA,MPROC),e(MDATA,MPROC),
     &fprod(0:MDATA,MPROC)

      character*13 filename
      read(5,*) nproc, nf
      if(nproc.gt.MPROC) stop 'nproc.gt.MPROC'

      do iproc=0,nproc-1
        if(iproc.le.9) then
          write(filename,'(''walkalize.'',i1)') iproc
         elseif(iproc.le.99) then
          write(filename,'(''walkalize.'',i2)') iproc
         elseif(iproc.le.999) then
          write(filename,'(''walkalize.'',i3)') iproc
         else
          stop 'iproc > 999'
        endif
        open(11,file=filename,status='old')

        if(iproc.eq.0) then
          read(11,*) nblkeq
          read(11,*) nstep,nblk,nconf,etrial,tau,taueff
          nskip=2*nstep*nblkeq
          ndata=nstep*nblk
          if(nskip+ndata.gt.MDATA) stop 'MDATA exceeded'
         else
          read(11,*)
        endif
        do i=1,nskip+ndata
          read(11,*) ii,f(i,iproc+1),w(i,iproc+1),e(i,iproc+1)
        enddo
      enddo

      do if=0,nf
        esum=0
        wsum=0
        w2sum=0
        wgsum=0
        wdsum=0
        pow=1-1/dfloat(if+1)
        do iproc=1,nproc
          fprod(0,iproc)=1
          do i=1,nskip+ndata
            fprod(i,iproc)=fprod(i-1,iproc)**pow*f(i,iproc)
          enddo
          do i=nskip+1,nskip+ndata
            esum=esum+fprod(i,iproc)*w(i,iproc)*e(i,iproc)
            wsum=wsum+fprod(i,iproc)*w(i,iproc)
            w2sum=w2sum+(fprod(i,iproc)*w(i,iproc))**2
            wgsum=wgsum+fprod(i-1,iproc)*f(i,iproc)*w(i,iproc)
            wdsum=wdsum+fprod(i-1,iproc)*w(i-1,iproc) 
          enddo
        enddo
        emix=esum/wsum
        eig=wgsum/wdsum
        egro=etrial-dlog(eig)/taueff
        data=ndata*nproc
        werr=sqrt(w2sum/data-(wsum/data)**2)
        write(6,'(f8.5,9f14.7)') if*taueff,emix,egro,werr
      enddo
        stop
        end
