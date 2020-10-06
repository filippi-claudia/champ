C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'atoms','forces','gradients','3dgrid','blocking_dmc',
     $ 'vmc','electrons','optwf','optgeo','pseudo','qmmm','strech',
     $ 'startend','properties','blocking_vmc','jastrow','general',
     $ 'dmc','mstates','ci','periodic'/
      data vars/'nctype','natom','addghostype','nghostcent','istrech',
     $ 'alfstr','nwprod','itausec','ngradnts','igrdtype','delgrdxyz',
     $ 'delgrdbl','delgrdba','delgrdda','stepx','stepy','stepz','x0',
     $ 'y0','z0','xn','yn','zn','nstep','nblk','nblkeq','nconf_new',
     $ 'nconf','imetro','deltar','deltat','delta','fbias',
     $ 'node_cutoff','enode_cutoff','nelec','nup','ioptwf','idl_flag',
     $ 'ilbfgs_flag','ilbfgs_m','method','nopt_iter','ioptjas',
     $ 'ioptorb','ioptci','multiple_adiag','add_diag',
     $ 'ngrad_jas_blocks','nblk_max','nblk_ci','dl_alg','iorbprt',
     $ 'isample_cmat','istddev','limit_cmat','e_shift','save_blocks',
     $ 'force_blocks','iorbsample','iuse_trafo','iuse_orbeigv',
     $ 'ncore','nextorb','no_active','approx','approx_mix',
     $ 'energy_tol','sr_tau','sr_eps','sr_adiag','micro_iter_sr',
     $ 'dl_mom','lin_eps','lin_adiag','lin_nvec','lin_nvecx',
     $ 'lin_jdav','func_omega','omega','n_omegaf','n_omegat',
     $ 'sr_rescale','iforce_analy','iuse_zmat','alfgeo','izvzb',
     $ 'iroot_geo','nloc','nquad','iqmmm','alfstr','idump','irstar',
     $ 'isite','icharged_atom','sample','print','nstep','nblk',
     $ 'nblkeq','nconf_new','nconf','ianalyt_lap','ijas','isc',
     $ 'nspin1','nspin2','ifock','title','unit','mass','iperiodic',
     $ 'ibasis','nforce','nwftype','seed','ipr','pool','basis',
     $ 'pseudopot','i3dsplorb','i3dlagorb','scalecoef','idmc','tau',
     $ 'etrial','nfprod','ipq','itau_eff','iacc_rej','icross',
     $ 'icuspg','idiv_v','icut_br','icut_e','icasula','node_cutoff',
     $ 'enode_cutoff','ibranch_elec','icircular','idrifdifgfunc',
     $ 'mode_dmc','iguiding','iefficiency','iciprt','norb','npoly',
     $ 'np','cutg','cutg_sim','cutg_big','cutg_sim_big'/
      data iaptr/1,5,9,15,24,29,36,38,84,89,91,92,93,97,99,104,110,
     $ 125,144,146,147/
      data ieptr/4,8,14,23,28,35,37,83,88,90,91,92,96,98,103,109,124,
     $ 143,145,146,153/
      nlist=21
      ierr=0
      do i=1,nlist
       if(lists(i).eq.p) then
        do iv=iaptr(i),ieptr(i)
         if(vars(iv).eq.v) then
          return
         endif
        enddo
        ierr=1
        return
       endif
      enddo
      return
      end
