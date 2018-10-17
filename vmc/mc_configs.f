      subroutine mc_configs_start

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /config/ xold(3,MELEC),xnew(3,MELEC),vold(3,MELEC)
     &,vnew(3,MELEC),psi2o(MSTATES,MFORCE),psi2n(MFORCE),eold(MSTATES,MFORCE),enew(MFORCE)
     &,peo(MSTATES),pen,tjfn,tjfo(MSTATES),psido(MSTATES),psijo
     &,rmino(MELEC),rminn(MELEC),rvmino(3,MELEC),rvminn(3,MELEC)
     &,rminon(MELEC),rminno(MELEC),rvminon(3,MELEC),rvminno(3,MELEC)
     &,nearesto(MELEC),nearestn(MELEC),delttn(MELEC)
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent

      character*20 filename

      dimension nsite(MCENT)

      if(irstar.ne.1) then
c check sites flag if one gets initial configuration from sites routine
        if (isite.eq.1) goto 20
        open(unit=9,err=20,file='mc_configs_start')
        rewind 9
        read(9,*,end=20,err=20) ((xold(k,i),k=1,3),i=1,nelec)
        write(6,'(/,''initial configuration from unit 9'')')
        goto 40
   20   l=0
	ntotal_sites=0
	do 25 i=1,ncent
   25     ntotal_sites=ntotal_sites+int(znuc(iwctype(i))+0.5d0)
        icharge_system=ntotal_sites-nelec
        call p2gtid('startend:icharged_atom',icharged_atom,0,1) 
        do 30 i=1,ncent
          nsite(i)=int(znuc(iwctype(i))+0.5d0)
	  if (icharged_atom.eq.i) then
            nsite(i)=int(znuc(iwctype(i))+0.5d0)-icharge_system
	    if (nsite(i).lt.0) call fatal_error('MC_CONFIG: error in icharged_atom')
	  endif
          l=l+nsite(i)
          if (l.gt.nelec) then
            nsite(i)=nsite(i)-(l-nelec)
            l=nelec
          endif
   30   continue
        if (l.lt.nelec) nsite(1)=nsite(1)+(nelec-l)
	
        call sites(xold,nelec,nsite)
        open(unit=9,file='mc_configs_start')
        rewind 9
        write(6,'(/,''initial configuration from sites'')')
   40   continue

c If we are moving one electron at a time, then we need to initialize
c xnew, since only the first electron gets initialized in metrop
        do 50 i=1,nelec
          do 50 k=1,3
   50       xnew(k,i)=xold(k,i)
      endif

c If nconf_new > 0 then we want to dump configurations for a future
c optimization or dmc calculation. So figure out how often we need to write a
c configuration to produce nconf_new configurations. If nconf_new = 0
c then set up so no configurations are written.
      if (nconf_new.ne.0) then
        filename='mc_configs_new'
        open(unit=7,form='formatted',file=filename)
        rewind 7
      endif

      call pcm_qvol(1)

      return

c-----------------------------------------------------------------------
      entry mc_configs_write

      rewind 9
      write(9,*) ((xold(ic,i),ic=1,3),i=1,nelec)
      close(9)

      call fin_reduce

      return
      end
