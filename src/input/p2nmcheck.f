C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'gradients','dmc','mstates','qmmm','startend','vmc',
     $ 'jastrow','forces','optwf','ci','general','periodic',
     $ 'blocking_dmc','properties','blocking_vmc','electrons',
     $ 'strech','atoms','pseudo','3dgrid','optgeo'/
      data vars/'ngradnts','igrdtype','delgrdxyz','delgrdbl',
     $ 'delgrdba','delgrdda','idmc','tau','etrial','nfprod','ipq',
     $ 'itau_eff','iacc_rej','icross','icuspg','idiv_v','icut_br',
     $ 'icut_e','icasula','node_cutoff','enode_cutoff','ibranch_elec',
     $ 'icircular','idrifdifgfunc','mode_dmc','iguiding',
     $ 'iefficiency','iqmmm','idump','irstar','isite','icharged_atom',
     $ 'imetro','deltar','deltat','delta','fbias','node_cutoff',
     $ 'enode_cutoff','ianalyt_lap','ijas','isc','nspin1','nspin2',
     $ 'ifock','istrech','alfstr','nwprod','itausec','ioptwf',
     $ 'idl_flag','ilbfgs_flag','ilbfgs_m','method','nopt_iter',
     $ 'ioptjas','ioptorb','ioptci','multiple_adiag','add_diag',
     $ 'ngrad_jas_blocks','nblk_max','nblk_ci','dl_alg','iorbprt',
     $ 'isample_cmat','istddev','limit_cmat','e_shift','save_blocks',
     $ 'force_blocks','iorbsample','iuse_trafo','iuse_orbeigv',
     $ 'ncore','nextorb','no_active','approx','approx_mix',
     $ 'energy_tol','sr_tau','sr_eps','sr_adiag','micro_iter_sr',
     $ 'dl_mom','lin_eps','lin_adiag','lin_nvec','lin_nvecx',
     $ 'lin_jdav','func_omega','omega','n_omegaf','n_omegat',
     $ 'sr_rescale','iciprt','title','unit','mass','iperiodic',
     $ 'ibasis','nforce','nwftype','seed','ipr','pool','basis',
     $ 'pseudopot','i3dsplorb','i3dlagorb','scalecoef','norb','npoly',
     $ 'np','cutg','cutg_sim','cutg_big','cutg_sim_big','nstep',
     $ 'nblk','nblkeq','nconf_new','nconf','sample','print','nstep',
     $ 'nblk','nblkeq','nconf_new','nconf','nelec','nup','alfstr',
     $ 'nctype','natom','addghostype','nghostcent','nloc','nquad',
     $ 'stepx','stepy','stepz','x0','y0','z0','xn','yn','zn',
     $ 'iforce_analy','iuse_zmat','alfgeo','izvzb','iroot_geo'/
      data iaptr/1,7,26,28,29,33,40,46,50,96,97,112,119,124,126,131,
     $ 133,134,138,140,149/
      data ieptr/6,25,27,28,32,39,45,49,95,96,111,118,123,125,130,132,
     $ 133,137,139,148,153/
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
