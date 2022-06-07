      implicit real*8(a-h,o-z)
c Written by Cyrus Umrigar
c Calculates the mixed and growth energies (using data in fort.11
c created in a dmc run) for various f projection times.
c At present it does it in a rather inefiicient ways.  Instead of
c recalculating the various f products could one could update them by
c multiplying and dividing by the appropriate f's.
      parameter (MDATA=500000)
      dimension f(MDATA),w(MDATA),e(MDATA)

      read(5,*) nf
      read(11,*) nstep,nblk,nblkeq,nconf,etrial,tau,taueff
      nskip=2*nstep*nblkeq
      ndata=nstep*nblk
      if(nskip+ndata.gt.MDATA) stop 'MDATA exceeded'
      do i=1,nskip+ndata
        read(11,*) istep,f(i),w(i),e(i)
      enddo

      do if=0,nf
        esum=0
        wsum=0
        w2sum=0
        wgsum=0
        wdsum=0
        do i=nskip+1,nskip+ndata
          fprod=1.d0
          fprodd=1.d0
          do j=0,if-1
            if(i-j.ge.1) fprod=fprod*f(i-j)
            if(i-j-1.ge.1) fprodd=fprodd*f(i-j-1)
          enddo
          fprodg=fprodd*f(i)
          esum=esum+fprod*w(i)*e(i)
          wsum=wsum+fprod*w(i)
          w2sum=w2sum+(fprod*w(i))**2
          wgsum=wgsum+fprodg*w(i)
          wdsum=wdsum+fprodd*w(i-1) 
        enddo
        emix=esum/wsum
        eig=wgsum/wdsum
        egro=etrial-dlog(eig)/taueff
        werr=sqrt(w2sum/ndata-(wsum/ndata)**2)
        write(6,'(f8.5,9f12.7)') if*taueff,emix,egro,werr
      enddo
        stop
        end
