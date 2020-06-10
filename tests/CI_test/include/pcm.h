      parameter (MCHS=1,MCHV=1,MSPHERE=30)
      character*80 pcmfile_cavity,pcmfile_chs,pcmfile_chv
      common /pcm_cntrl/ ipcm,ichpol,isurf,icall,ipcmprt
      common /pcm_unit/ pcmfile_cavity,pcmfile_chs,pcmfile_chv
      common /pcm_parms/ xpol(3,MCHV),ch(MCHV),xe(MSPHERE),ye(MSPHERE),ze(MSPHERE),re(MSPHERE),re2(MSPHERE)
      common /pcm_parms/ eps_solv,retk,surk,nesph,nchs,nchv,nch,ncopcm,nvopcm,nscv,iscov,nchs1,nchs2
      common /pcm_ameta/amdlg(MCHS),eta(3,MCHS)
      common /pcm_ah/ahca(MCHS,MCHS),bh(MCHS)

      common /pcm_pot/ penups,penupv,penupol,pepol
      common /pcm_fdc/ qfree,qvol,fdc(MCHS),fs,feps,rcolv,rcol,rcolt
      common /pcm_averages/ spcmsum,spcmcum,spcmcm2 
      common /pcm_averages/ vpcmsum,vpcmcum,vpcmcm2,qopcm_sum,qopcm_cum,qopcm_cm2
      common /pcm_averages/ enfpcm_sum(MCHS),enfpcm_cum(MCHS),enfpcm_cm2(MCHS)
      common /pcm_inda/inda(MCHS)



