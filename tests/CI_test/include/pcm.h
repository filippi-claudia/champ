      parameter (MCHS=1,MCHV=1,MSPHERE=30)
      common /pcm_ah/ahca(MCHS,MCHS),bh(MCHS)

      common /pcm_pot/ penups,penupv,penupol,pepol
      common /pcm_fdc/ qfree,qvol,fdc(MCHS),fs,feps,rcolv,rcol,rcolt
      common /pcm_averages/ spcmsum,spcmcum,spcmcm2 
      common /pcm_averages/ vpcmsum,vpcmcum,vpcmcm2,qopcm_sum,qopcm_cum,qopcm_cm2
      common /pcm_averages/ enfpcm_sum(MCHS),enfpcm_cum(MCHS),enfpcm_cm2(MCHS)
      common /pcm_inda/inda(MCHS)



