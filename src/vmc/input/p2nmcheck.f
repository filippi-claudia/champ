C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'jastrow','general','strech','qmmm','optgeo',
     $ '3dgrid','periodic','forces','ci','pseudo','blocking_vmc',
     $ 'optwf','dmc','atoms','blocking_dmc','electrons','startend',
     $ 'properties','gradients','vmc','mstates'/
      data vars/'ianalyt_lap','ijas','isc','nspin1','nspin2','ifock',
     $ 'title','unit','mass','iperiodic','ibasis','nforce','nwftype',
     $ 'seed','ipr','pool','basis','pseudopot','i3dsplorb',
     $ 'i3dlagorb','scalecoef','alfstr','iqmmm','iforce_analy',
     $ 'iuse_zmat','alfgeo','izvzb','iroot_geo','stepx','stepy',
     $ 'stepz','x0','y0','z0','xn','yn','zn','norb','npoly','np',
     $ 'cutg','cutg_sim','cutg_big','cutg_sim_big','istrech','alfstr',
     $ 'nwprod','itausec','iciprt','nloc','nquad','nstep','nblk',
     $ 'nblkeq','nconf_new','nconf','ioptwf','idl_flag','ilbfgs_flag',
     $ 'ilbfgs_m','method','nopt_iter','ioptjas','ioptorb','ioptci',
     $ 'multiple_adiag','add_diag','ngrad_jas_blocks','nblk_max',
     $ 'nblk_ci','dl_alg','iorbprt','isample_cmat','istddev',
     $ 'limit_cmat','e_shift','save_blocks','force_blocks',
     $ 'iorbsample','iuse_trafo','iuse_orbeigv','ncore','nextorb',
     $ 'no_active','approx','approx_mix','energy_tol','sr_tau',
     $ 'sr_eps','sr_adiag','micro_iter_sr','dl_mom','lin_eps',
     $ 'lin_adiag','lin_nvec','lin_nvecx','lin_jdav','func_omega',
     $ 'omega','n_omegaf','n_omegat','sr_rescale','idmc','tau',
     $ 'etrial','nfprod','ipq','itau_eff','iacc_rej','icross',
     $ 'icuspg','idiv_v','icut_br','icut_e','icasula','node_cutoff',
     $ 'enode_cutoff','ibranch_elec','icircular','idrifdifgfunc',
     $ 'mode_dmc','nctype','natom','addghostype','nghostcent','nstep',
     $ 'nblk','nblkeq','nconf_new','nconf','nelec','nup','idump',
     $ 'irstar','isite','icharged_atom','sample','print','ngradnts',
     $ 'igrdtype','delgrdxyz','delgrdbl','delgrdba','delgrdda',
     $ 'imetro','deltar','deltat','delta','fbias','node_cutoff',
     $ 'enode_cutoff','iguiding','iefficiency'/
      data iaptr/1,7,22,23,24,29,38,45,49,50,52,57,103,122,126,131,
     $ 133,137,139,145,152/
      data ieptr/6,21,22,23,28,37,44,48,49,51,56,102,121,125,130,132,
     $ 136,138,144,151,153/
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
