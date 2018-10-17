      subroutine mc_configs
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'

      character*25 fmt
      character*12 mode

      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /config/ xold(3,MELEC,MWALK,MFORCE),vold(3,MELEC,MWALK,MFORCE),
     &psido(MWALK,MFORCE),psijo(MWALK,MFORCE),peo(MWALK,MFORCE),d2o(MWALK,MFORCE)
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk

      common /contr3/ mode

      save ngfmc

      if (irstar.eq.1) then
        open(unit=10,status='old',form='unformatted',file='restart_dmc')
        rewind 10
        call startr
        close (unit=10)
       else
        open(unit=1,status='old',form='formatted',file='mc_configs')
        rewind 1
        do 330 i=1,nconf
          read(1,fmt=*,end=340) ((xold(ic,j,i,1),ic=1,3),j=1,nelec)
  330   continue
        goto 345
  340   call fatal_error('MAIN: error reading mc_configs in maindmc')
  345   close (1)
        if(ipr.gt.-2.and.index(mode,'rmc').eq.0) then
          open(11,file='walkalize')
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
        open(unit=7,form='formatted',file='mc_configs_new')
        rewind 7
      endif

      return

c-----------------------------------------------------------------------
      entry mc_configs_write(iblk,ipass)

c Write out configuration for optimization/dmc/gfmc here
c We would like to:
c Reduce each electron to central simulation cell before writing position.
c Warning: This may result in the sign of the wavefunction being wrong depending
c on the k-pt and the number of cells moved.  Test done in si_cub_.5.0.0 shows this
c gives wrong energies and moves, possibly because the velocity is no longer consistent
c with the move, though I would have thought the velocity would be OK since both the
c wavfn and its gradients would change sign simultaneously.  May be I need to move
c in only even multiples of the sim. lattice.
      if (iblk.gt.2*nblkeq .and. (mod(ipass,ngfmc).eq.1 .or.  ngfmc.eq.1)) then
        if(3*nelec.lt.100) then
          write(fmt,'(a1,i2,a21)')'(',3*nelec,'f14.8,i3,d12.4,f12.5)'
         else
          write(fmt,'(a1,i3,a21)')'(',3*nelec,'f14.8,i3,d12.4,f12.5)'
        endif
        do 352 iwalk=1,nwalk
c         if(iperiodic.ne.0) then
c           do 351 jj=1,nelec
c 351         call reduce_sim_cell(xold(1,jj,iwalk,1),rlatt_sim,rlatt_sim_inv)
c         endif
  352     write(7,fmt) ((xold(ii,jj,iwalk,1),ii=1,3),jj=1,nelec),
     &    int(sign(1.d0,psido(iwalk,1))),log(dabs(psido(iwalk,1)))+psijo(iwalk,1),eold(iwalk,1)
      endif

      return
      end
