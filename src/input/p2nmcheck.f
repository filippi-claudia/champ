C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'pseudo','qmmm','startend','vmc','electrons','dmc',
     $ 'optwf','forces','periodic','blocking_dmc','gradients','atoms',
     $ 'properties','mstates','jastrow','general','ci','3dgrid',
     $ 'blocking_vmc','optgeo','strech'/
      data vars/'nloc','nquad','iqmmm','idump','irstar','isite',
     $ 'icharged_atom','imetro','deltar','deltat','delta','fbias',
     $ 'node_cutoff','enode_cutoff','nelec','nup','idmc','tau',
     $ 'etrial','nfprod','ipq','itau_eff','iacc_rej','icross',
     $ 'icuspg','idiv_v','icut_br','icut_e','icasula','node_cutoff',
     $ 'enode_cutoff','ibranch_elec','icircular','idrifdifgfunc',
     $ 'mode_dmc','ioptwf','idl_flag','ilbfgs_flag','ilbfgs_m',
     $ 'method','nopt_iter','ioptjas','ioptorb','ioptci',
     $ 'multiple_adiag','add_diag','ngrad_jas_blocks','nblk_max',
     $ 'nblk_ci','dl_alg','iorbprt','isample_cmat','istddev',
     $ 'limit_cmat','e_shift','save_blocks','force_blocks',
     $ 'iorbsample','iuse_trafo','iuse_orbeigv','ncore','nextorb',
     $ 'no_active','approx','approx_mix','energy_tol','sr_tau',
     $ 'sr_eps','sr_adiag','micro_iter_sr','dl_mom','lin_eps',
     $ 'lin_adiag','lin_nvec','lin_nvecx','lin_jdav','func_omega',
     $ 'omega','n_omegaf','n_omegat','sr_rescale','istrech','alfstr',
     $ 'nwprod','itausec','norb','npoly','np','cutg','cutg_sim',
     $ 'cutg_big','cutg_sim_big','nstep','nblk','nblkeq','nconf_new',
     $ 'nconf','ngradnts','igrdtype','delgrdxyz','delgrdbl',
     $ 'delgrdba','delgrdda','nctype','natom','addghostype',
     $ 'nghostcent','sample','print','iguiding','iefficiency',
     $ 'ianalyt_lap','ijas','isc','nspin1','nspin2','ifock','title',
     $ 'unit','mass','iperiodic','ibasis','nforce','nwftype','seed',
     $ 'ipr','pool','basis','pseudopot','i3dsplorb','i3dlagorb',
     $ 'scalecoef','iciprt','stepx','stepy','stepz','x0','y0','z0',
     $ 'xn','yn','zn','nstep','nblk','nblkeq','nconf_new','nconf',
     $ 'iforce_analy','iuse_zmat','alfgeo','izvzb','iroot_geo',
     $ 'alfstr'/
      data iaptr/1,3,4,8,15,17,36,82,86,93,98,104,108,110,112,118,133,
     $ 134,143,148,153/
      data ieptr/2,3,7,14,16,35,81,85,92,97,103,107,109,111,117,132,
     $ 133,142,147,152,153/
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
