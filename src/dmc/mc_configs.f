      subroutine mc_configs
      use const, only: delta, deltai, etrial, fbias, hb, imetro, ipr, nelec, pi
      use config, only: d2o, peo_dmc, psido_dmc, psijo_dmc, vold_dmc, xold_dmc
      use mpiconf, only: idtask, nproc, wid, NPROCX

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'
      include 'mpif.h'

      character*25 fmt
      character*20 filename
      character*12 mode

      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk

      common /contr3/ mode

      dimension irn(4)
      save ngfmc

      if(ipr.gt.-2) then
        if(idtask.le.9) then
          write(filename,'(''walkalize.'',i1)') idtask
         elseif(idtask.le.99) then
          write(filename,'(''walkalize.'',i2)') idtask
         elseif(idtask.le.999) then
          write(filename,'(''walkalize.'',i3)') idtask
         else
          call fatal_error('DMC: idtask > 999')
        endif
      endif

c set the random number seed, setrn already called in read_input
      if(irstar.ne.1) then
        if(nproc.gt.1) then
          do 5 id=1,(3*nelec)*idtask
    5       rnd=rannyu(0)
          call savern(irn)
          do 6 i=1,4
    6       irn(i)=mod(irn(i)+int(rannyu(0)*idtask*9999),9999)
          call setrn(irn)
        endif
      endif

      if (irstar.eq.1) then
        open(unit=10,status='old',form='unformatted',file='restart_dmc')
        rewind 10
        call startr
        close (unit=10)
       else
        open(unit=1,status='old',form='formatted',file='mc_configs')
        rewind 1
        do 330 id=0,idtask
          do 330 i=1,nconf
            read(1,fmt=*,end=340) ((xold_dmc(ic,j,i,1),ic=1,3),j=1,nelec)
  330   continue
        goto 345
  340   call fatal_error('DMC: error reading mc_configs')
  345   close (1)
        if(ipr.gt.-2) then
          open(11,file=filename)
          rewind 11
          write(11,'(i3,'' nblkeq to be added to nblock at file end'')')
     &    nblkeq
        endif
      endif

c If nconf_new > 0, dump configurations for a future optimization or dmc calculation. 
c Figure out frequency of configuration writing to produce nconf_new configurations. 
c If nconf_new = 0, then no configurations are written.
      if (nconf_new.eq.0) then
        ngfmc=2*nstep*nblk
       else
        ngfmc=(nstep*nblk+nconf_new-1)*nconf/nconf_new
        if(idtask.lt.10) then
          write(filename,'(i1)') idtask
         elseif(idtask.lt.100) then
          write(filename,'(i2)') idtask
         elseif(idtask.lt.1000) then
          write(filename,'(i3)') idtask
         else
          write(filename,'(i4)') idtask
        endif
        filename='mc_configs_new'//filename(1:index(filename,' ')-1)
        open(unit=7,form='formatted',file=filename)
        rewind 7
      endif

      return

c-----------------------------------------------------------------------
      entry mc_configs_write(iblk,ipass)

c Write out configuration for optimization/dmc/gfmc here
          if (iblk.gt.2*nblkeq .and. (mod(ipass,ngfmc).eq.1 .or.  ngfmc.eq.1)) then
            if(3*nelec.lt.100) then
              write(fmt,'(a1,i2,a21)')'(',3*nelec,'f14.8,i3,d12.4,f12.5)'
             else
              write(fmt,'(a1,i3,a21)')'(',3*nelec,'f14.8,i3,d12.4,f12.5)'
            endif
            do 352 iwalk=1,nwalk
  352         write(7,fmt) ((xold_dmc(ii,jj,iwalk,1),ii=1,3),jj=1,nelec),
     &        int(sign(1.d0,psido_dmc(iwalk,1))),log(dabs(psido_dmc(iwalk,1)))+psijo_dmc(iwalk,1),eold(iwalk,1)
          endif

      return
      end
