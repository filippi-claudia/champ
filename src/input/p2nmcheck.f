C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'mstates','qmmm','vmc','ci','general','dmc',
     $ 'blocking_dmc','optgeo','properties','blocking_vmc','periodic',
     $ 'electrons','atoms','optwf','pseudo','startend','gradients',
     $ 'jastrow','forces','3dgrid','strech'/
      data vars/'iguiding','iefficiency','iqmmm','imetro','deltar',
     $ 'deltat','delta','fbias','node_cutoff','enode_cutoff','iciprt',
     $ 'title','unit','mass','iperiodic','ibasis','nforce','nwftype',
     $ 'seed','ipr','pool','basis','pseudopot','i3dsplorb',
     $ 'i3dlagorb','scalecoef','idmc','tau','etrial','nfprod','ipq',
     $ 'itau_eff','iacc_rej','icross','icuspg','idiv_v','icut_br',
     $ 'icut_e','icasula','node_cutoff','enode_cutoff','ibranch_elec',
     $ 'icircular','idrifdifgfunc','mode_dmc','nstep','nblk','nblkeq',
     $ 'nconf_new','nconf','iforce_analy','iuse_zmat','alfgeo',
     $ 'izvzb','iroot_geo','sample','print','nstep','nblk','nblkeq',
     $ 'nconf_new','nconf','norb','npoly','np','cutg','cutg_sim',
     $ 'cutg_big','cutg_sim_big','nelec','nup','nctype','natom',
     $ 'addghostype','nghostcent','ioptwf','idl_flag','ilbfgs_flag',
     $ 'ilbfgs_m','method','nopt_iter','ioptjas','ioptorb','ioptci',
     $ 'multiple_adiag','add_diag','ngrad_jas_blocks','nblk_max',
     $ 'nblk_ci','dl_alg','iorbprt','isample_cmat','istddev',
     $ 'limit_cmat','e_shift','save_blocks','force_blocks',
     $ 'iorbsample','iuse_trafo','iuse_orbeigv','ncore','nextorb',
     $ 'no_active','approx','approx_mix','energy_tol','sr_tau',
     $ 'sr_eps','sr_adiag','micro_iter_sr','dl_mom','lin_eps',
     $ 'lin_adiag','lin_nvec','lin_nvecx','lin_jdav','func_omega',
     $ 'omega','n_omegaf','n_omegat','sr_rescale','nloc','nquad',
     $ 'idump','irstar','isite','icharged_atom','ngradnts','igrdtype',
     $ 'delgrdxyz','delgrdbl','delgrdba','delgrdda','ianalyt_lap',
     $ 'ijas','isc','nspin1','nspin2','ifock','istrech','alfstr',
     $ 'nwprod','itausec','stepx','stepy','stepz','x0','y0','z0','xn',
     $ 'yn','zn','alfstr'/
      data iaptr/1,3,4,11,12,27,46,51,56,58,63,70,72,76,122,124,128,
     $ 134,140,144,153/
      data ieptr/2,3,10,11,26,45,50,55,57,62,69,71,75,121,123,127,133,
     $ 139,143,152,153/
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
