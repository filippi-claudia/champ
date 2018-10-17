      implicit real*8(a-h,o-z)
c Written by Cyrus Umrigar (Mar 1992, Sep. 2001)
c Modified by Claudia Filippi (May 2004) to account for glop. pop.
c
c Calculates the mixed and growth energies (using data in fort.11
c created in a dmc run) for various f projection times.  Instead
c of using a certain number of f's it uses all f's but raising them
c to a power, less than one, so that the contributions from the
c ones further back tend to 1.  This gives somewhat smoother curves
c than the old way.
      parameter (MDATA=104000,MPROC=100)
      dimension f(MDATA,MPROC),w(MDATA,MPROC),e(MDATA,MPROC),
     &fprod(0:MDATA,MPROC)

      character*13 filename
c igp=0  separate populations
c igp=1  global population
      read(5,*) nproc, nf,igp
      if(nproc.gt.MPROC) stop 'nproc.gt.MPROC'

      do 5 iproc=0,nproc-1
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
        do 5 i=1,nskip+ndata
    5     read(11,*) ii,f(i,iproc+1),w(i,iproc+1),e(i,iproc+1)

      do 30 if=0,nf
        pow=1-1/dfloat(if+1)
        do 15 iproc=1,nproc
          fprod(0,iproc)=1
          do 15 i=1,nskip+ndata
   15       fprod(i,iproc)=fprod(i-1,iproc)**pow*f(i,iproc)
        esum=0
        wsum=0
        w2sum=0
        wgsum=0
        wdsum=0
        do 20 i=nskip+1,nskip+ndata
          wtmp=0
          do 20 iproc=1,nproc
            esum=esum+fprod(i,iproc)*w(i,iproc)*e(i,iproc)
            wsum=wsum+fprod(i,iproc)*w(i,iproc)
            wgsum=wgsum+fprod(i-1,iproc)*f(i,iproc)*w(i,iproc)
            wdsum=wdsum+fprod(i-1,iproc)*w(i-1,iproc) 
            wtmp=wtmp+w(i,iproc)
            if(igp.eq.0) then
              w2sum=w2sum+(fprod(i,iproc)*w(i,iproc))**2
             elseif(iproc.eq.nproc) then
              w2sum=w2sum+(fprod(i,iproc)*wtmp*wtmp/nproc)**2
            endif 
   20   continue
        emix=esum/wsum
        eig=wgsum/wdsum
        egro=etrial-dlog(eig)/taueff
        data=ndata*nproc
        data2=data
        if(igp.eq.1) data2=ndata
        werr=sqrt(w2sum/data2-(wsum/data)**2)
   30   write(6,'(f8.5,9f14.7)') if*taueff,emix,egro,werr
        stop
        end
