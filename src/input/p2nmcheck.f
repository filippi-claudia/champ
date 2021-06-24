C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'periodic','properties','electrons','jastrow',
     $ 'optwf','atoms','general','blocking_vmc','3dgrid','optgeo',
     $ 'strech','qmmm','forces','gradients','blocking_dmc','ci',
     $ 'startend','dmc','pseudo','vmc','mstates'/
      data vars/'norb','npoly','np','cutg','cutg_sim','cutg_big',
     $ 'cutg_sim_big','sample','print','nelec','nup','ianalyt_lap',
     $ 'ijas','isc','nspin1','nspin2','ifock','ioptwf','idl_flag',
     $ 'ilbfgs_flag','ilbfgs_m','method','nopt_iter','ioptjas',
     $ 'ioptorb','ioptci','multiple_adiag','add_diag',
     $ 'ngrad_jas_blocks','nblk_max','nblk_ci','dl_alg','iorbprt',
     $ 'isample_cmat','istddev','limit_cmat','e_shift','save_blocks',
     $ 'force_blocks','iorbsample','iuse_trafo','iuse_orbeigv',
     $ 'ncore','nextorb','no_active','approx','approx_mix',
     $ 'energy_tol','sr_tau','sr_eps','sr_adiag','micro_iter_sr',
     $ 'dl_mom','lin_eps','lin_adiag','lin_nvec','lin_nvecx',
     $ 'lin_jdav','func_omega','omega','n_omegaf','n_omegat',
     $ 'sr_rescale','nctype','natom','addghostype','nghostcent',
     $ 'title','unit','mass','iperiodic','ibasis','nforce','nwftype',
     $ 'seed','ipr','pool','basis','pseudopot','i3dsplorb',
     $ 'i3dlagorb','scalecoef','nstep','nblk','nblkeq','nconf_new',
     $ 'nconf','stepx','stepy','stepz','x0','y0','z0','xn','yn','zn',
     $ 'iforce_analy','iuse_zmat','alfgeo','izvzb','iroot_geo',
     $ 'alfstr','iqmmm','istrech','alfstr','nwprod','itausec',
     $ 'ngradnts','igrdtype','delgrdxyz','delgrdbl','delgrdba',
     $ 'delgrdda','nstep','nblk','nblkeq','nconf_new','nconf',
     $ 'iciprt','idump','irstar','isite','icharged_atom','idmc','tau',
     $ 'etrial','nfprod','ipq','itau_eff','iacc_rej','icross',
     $ 'icuspg','idiv_v','icut_br','icut_e','icasula','node_cutoff',
     $ 'enode_cutoff','ibranch_elec','icircular','idrifdifgfunc',
     $ 'mode_dmc','nloc','nquad','imetro','deltar','deltat','delta',
     $ 'fbias','node_cutoff','enode_cutoff','iguiding','iefficiency'/
      data iaptr/1,8,10,12,18,64,68,83,88,97,102,103,104,108,114,119,
     $ 120,124,143,145,152/
      data ieptr/7,9,11,17,63,67,82,87,96,101,102,103,107,113,118,119,
     $ 123,142,144,151,153/
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
