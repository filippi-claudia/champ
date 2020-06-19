      parameter (MCHMM=1)

      character*80 mmpolfile_sites, mmpolfile_chmm

      common /mmpol_pot/ penu_dp,penu_q,peqq,peq_dp, u_dd, u_self, pepol_dp,pepol_q
      common /mmpol_fdc/ rcolm,a_cutoff,screen1(MCHMM,MCHMM),screen2(MCHMM,MCHMM)
      common /mmpol_inds/inds_pol(MCHMM)
