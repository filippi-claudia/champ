C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'forces','jastrow','blocking_vmc','electrons',
     $ 'gradients','atoms','strech','periodic','properties','vmc',
     $ 'startend','pseudo','3dgrid','general','qmmm','optwf','optgeo',
     $ 'dmc','ci','blocking_dmc','mstates'/
      data vars/'istrech','alfstr','nwprod','itausec','ianalyt_lap',
     $ 'ijas','isc','nspin1','nspin2','ifock','nstep','nblk','nblkeq',
     $ 'nconf_new','nconf','nelec','nup','ngradnts','igrdtype',
     $ 'delgrdxyz','delgrdbl','delgrdba','delgrdda','nctype','natom',
     $ 'addghostype','nghostcent','alfstr','norb','npoly','np','cutg',
     $ 'cutg_sim','cutg_big','cutg_sim_big','sample','print','imetro',
     $ 'deltar','deltat','delta','fbias','node_cutoff','enode_cutoff',
     $ 'idump','irstar','isite','icharged_atom','nloc','nquad',
     $ 'stepx','stepy','stepz','x0','y0','z0','xn','yn','zn','title',
     $ 'unit','mass','iperiodic','ibasis','nforce','nwftype','seed',
     $ 'ipr','pool','basis','pseudopot','i3dsplorb','i3dlagorb',
     $ 'scalecoef','iqmmm','ioptwf','idl_flag','ilbfgs_flag',
     $ 'ilbfgs_m','method','nopt_iter','ioptjas','ioptorb','ioptci',
     $ 'multiple_adiag','add_diag','ngrad_jas_blocks','nblk_max',
     $ 'nblk_ci','dl_alg','iorbprt','isample_cmat','istddev',
     $ 'limit_cmat','e_shift','save_blocks','force_blocks',
     $ 'iorbsample','iuse_trafo','iuse_orbeigv','ncore','nextorb',
     $ 'no_active','approx','approx_mix','energy_tol','sr_tau',
     $ 'sr_eps','sr_adiag','micro_iter_sr','dl_mom','lin_eps',
     $ 'lin_adiag','lin_nvec','lin_nvecx','lin_jdav','func_omega',
     $ 'omega','n_omegaf','n_omegat','sr_rescale','iforce_analy',
     $ 'iuse_zmat','alfgeo','izvzb','iroot_geo','idmc','tau','etrial',
     $ 'nfprod','ipq','itau_eff','iacc_rej','icross','icuspg',
     $ 'idiv_v','icut_br','icut_e','icasula','node_cutoff',
     $ 'enode_cutoff','ibranch_elec','icircular','idrifdifgfunc',
     $ 'mode_dmc','iciprt','nstep','nblk','nblkeq','nconf_new',
     $ 'nconf','iguiding','iefficiency'/
      data iaptr/1,5,11,16,18,24,28,29,36,38,45,49,51,60,75,76,122,
     $ 127,146,147,152/
      data ieptr/4,10,15,17,23,27,28,35,37,44,48,50,59,74,75,121,126,
     $ 145,146,151,153/
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
