C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'periodic','qmmm','pseudo','strech','jastrow',
     $ 'blocking_dmc','dmc','atoms','optgeo','3dgrid','general',
     $ 'startend','properties','vmc','mstates','electrons',
     $ 'gradients','blocking_vmc','ci','optwf','forces'/
      data vars/'norb','npoly','np','cutg','cutg_sim','cutg_big',
     $ 'cutg_sim_big','iqmmm','nloc','nquad','alfstr','ianalyt_lap',
     $ 'ijas','isc','nspin1','nspin2','ifock','nstep','nblk','nblkeq',
     $ 'nconf_new','nconf','idmc','tau','etrial','nfprod','ipq',
     $ 'itau_eff','iacc_rej','icross','icuspg','idiv_v','icut_br',
     $ 'icut_e','icasula','node_cutoff','enode_cutoff','ibranch_elec',
     $ 'icircular','idrifdifgfunc','mode_dmc','nctype','natom',
     $ 'addghostype','nghostcent','iforce_analy','iuse_zmat','alfgeo',
     $ 'izvzb','iroot_geo','stepx','stepy','stepz','x0','y0','z0',
     $ 'xn','yn','zn','title','unit','mass','iperiodic','ibasis',
     $ 'nforce','nwftype','seed','ipr','pool','basis','pseudopot',
     $ 'i3dsplorb','i3dlagorb','scalecoef','idump','irstar','isite',
     $ 'icharged_atom','sample','print','imetro','deltar','deltat',
     $ 'delta','fbias','node_cutoff','enode_cutoff','iguiding',
     $ 'iefficiency','nelec','nup','ngradnts','igrdtype','delgrdxyz',
     $ 'delgrdbl','delgrdba','delgrdda','nstep','nblk','nblkeq',
     $ 'nconf_new','nconf','iciprt','ioptwf','idl_flag','ilbfgs_flag',
     $ 'ilbfgs_m','method','nopt_iter','ioptjas','ioptorb','ioptci',
     $ 'multiple_adiag','add_diag','ngrad_jas_blocks','nblk_max',
     $ 'nblk_ci','dl_alg','iorbprt','isample_cmat','istddev',
     $ 'limit_cmat','e_shift','save_blocks','force_blocks',
     $ 'iorbsample','iuse_trafo','iuse_orbeigv','ncore','nextorb',
     $ 'no_active','approx','approx_mix','energy_tol','sr_tau',
     $ 'sr_eps','sr_adiag','micro_iter_sr','dl_mom','lin_eps',
     $ 'lin_adiag','lin_nvec','lin_nvecx','lin_jdav','func_omega',
     $ 'omega','n_omegaf','n_omegat','sr_rescale','istrech','alfstr',
     $ 'nwprod','itausec'/
      data iaptr/1,8,9,11,12,18,23,42,46,51,60,75,79,81,88,90,92,98,
     $ 103,104,150/
      data ieptr/7,8,10,11,17,22,41,45,50,59,74,78,80,87,89,91,97,102,
     $ 103,149,153/
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
