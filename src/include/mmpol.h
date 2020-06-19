      parameter (MCHMM=1)

      character*80 mmpolfile_sites, mmpolfile_chmm

      common /mmpol_fdc/ rcolm,a_cutoff,screen1(MCHMM,MCHMM),screen2(MCHMM,MCHMM)
      common /mmpol_inds/inds_pol(MCHMM)
