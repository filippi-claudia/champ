C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'mstates','strech','startend','periodic','pseudo',
     $ 'vmc','properties','jastrow','optgeo','qmmm','electrons',
     $ 'gradients','3dgrid','ci','optwf','blocking_vmc',
     $ 'blocking_dmc','atoms','forces','dmc','general'/
      data vars/'iguiding','iefficiency','alfstr','idump','irstar',
     $ 'isite','icharged_atom','norb','npoly','np','cutg','cutg_sim',
     $ 'cutg_big','cutg_sim_big','nloc','nquad','imetro','deltar',
     $ 'deltat','delta','fbias','node_cutoff','enode_cutoff','sample',
     $ 'print','ianalyt_lap','ijas','isc','nspin1','nspin2','ifock',
     $ 'iforce_analy','iuse_zmat','alfgeo','izvzb','iroot_geo',
     $ 'iqmmm','nelec','nup','ngradnts','igrdtype','delgrdxyz',
     $ 'delgrdbl','delgrdba','delgrdda','stepx','stepy','stepz','x0',
     $ 'y0','z0','xn','yn','zn','iciprt','ioptwf','idl_flag',
     $ 'ilbfgs_flag','ilbfgs_m','method','nopt_iter','ioptjas',
     $ 'ioptorb','ioptci','multiple_adiag','add_diag',
     $ 'ngrad_jas_blocks','nblk_max','nblk_ci','dl_alg','iorbprt',
     $ 'isample_cmat','istddev','limit_cmat','e_shift','save_blocks',
     $ 'force_blocks','iorbsample','iuse_trafo','iuse_orbeigv',
     $ 'ncore','nextorb','no_active','approx','approx_mix',
     $ 'energy_tol','sr_tau','sr_eps','sr_adiag','micro_iter_sr',
     $ 'dl_mom','lin_eps','lin_adiag','lin_nvec','lin_nvecx',
     $ 'lin_jdav','func_omega','omega','n_omegaf','n_omegat',
     $ 'sr_rescale','nstep','nblk','nblkeq','nconf_new','nconf',
     $ 'nstep','nblk','nblkeq','nconf_new','nconf','nctype','natom',
     $ 'addghostype','nghostcent','istrech','alfstr','nwprod',
     $ 'itausec','idmc','tau','etrial','nfprod','ipq','itau_eff',
     $ 'iacc_rej','icross','icuspg','idiv_v','icut_br','icut_e',
     $ 'icasula','node_cutoff','enode_cutoff','ibranch_elec',
     $ 'icircular','idrifdifgfunc','mode_dmc','title','unit','mass',
     $ 'iperiodic','ibasis','nforce','nwftype','seed','ipr','pool',
     $ 'basis','pseudopot','i3dsplorb','i3dlagorb','scalecoef'/
      data iaptr/1,3,4,8,15,17,24,26,32,37,38,40,46,55,56,102,107,112,
     $ 116,120,139/
      data ieptr/2,3,7,14,16,23,25,31,36,37,39,45,54,55,101,106,111,
     $ 115,119,138,153/
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
