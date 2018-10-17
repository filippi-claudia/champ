      implicit real*8(a-h,o-z)
c Written by Cyrus Umrigar
c Calculates the mixed and growth energies (using data in fort.11
c created in a dmc run) for various f projection times.  Instead
c of using a certain number of f's it uses all f's but raising them
c to a power, less than one, so that the contributions from the
c ones further back tend to 1.  This gives somewhat smoother curves
c than the old way.       Cyrus Umrigar (Mar 1992, Sep. 2001)
      parameter (MDATA=500000)
      dimension f(MDATA),w(MDATA),e(MDATA),fprod(0:MDATA)

      read(5,*) nf
      read(11,*) nstep,nblk,nblkeq,nconf,etrial,tau,taueff
      nskip=2*nstep*nblkeq
      ndata=nstep*nblk
      if(nskip+ndata.gt.MDATA) stop 'MDATA exceeded'
      do 5 i=1,nskip+ndata
    5   read(11,*) istep,f(i),w(i),e(i)

      do 30 if=0,nf
        esum=0
        wsum=0
        w2sum=0
        wgsum=0
        wdsum=0
        pow=1-1/dfloat(if+1)
        fprod(0)=1
        do 15 i=1,nskip+ndata
   15     fprod(i)=fprod(i-1)**pow*f(i)
        do 20 i=nskip+1,nskip+ndata
          esum=esum+fprod(i)*w(i)*e(i)
          wsum=wsum+fprod(i)*w(i)
          w2sum=w2sum+(fprod(i)*w(i))**2
          wgsum=wgsum+fprod(i-1)*f(i)*w(i)
          wdsum=wdsum+fprod(i-1)*w(i-1) 
   20     continue
        emix=esum/wsum
        eig=wgsum/wdsum
        egro=etrial-dlog(eig)/taueff
        werr=sqrt(w2sum/ndata-(wsum/ndata)**2)
   30   write(6,'(f8.5,9f12.7)') if*taueff,emix,egro,werr
        stop
        end
