C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'electrons','pseudo','qmmm','properties',
     $ 'blocking_vmc','startend','jastrow','forces','strech',
     $ 'gradients','blocking_dmc','optgeo','dmc','mstates','periodic',
     $ 'vmc','3dgrid','optwf','general','ci','atoms'/
      data vars/'nelec','nup','nloc','nquad','iqmmm','sample','print',
     $ 'nstep','nblk','nblkeq','nconf_new','nconf','idump','irstar',
     $ 'isite','icharged_atom','ianalyt_lap','ijas','isc','nspin1',
     $ 'nspin2','ifock','istrech','alfstr','nwprod','itausec',
     $ 'alfstr','ngradnts','igrdtype','delgrdxyz','delgrdbl',
     $ 'delgrdba','delgrdda','nstep','nblk','nblkeq','nconf_new',
     $ 'nconf','iforce_analy','iuse_zmat','alfgeo','izvzb',
     $ 'iroot_geo','idmc','tau','etrial','nfprod','ipq','itau_eff',
     $ 'iacc_rej','icross','icuspg','idiv_v','icut_br','icut_e',
     $ 'icasula','node_cutoff','enode_cutoff','ibranch_elec',
     $ 'icircular','idrifdifgfunc','mode_dmc','iguiding',
     $ 'iefficiency','norb','npoly','np','cutg','cutg_sim','cutg_big',
     $ 'cutg_sim_big','imetro','deltar','deltat','delta','fbias',
     $ 'node_cutoff','enode_cutoff','stepx','stepy','stepz','x0','y0',
     $ 'z0','xn','yn','zn','ioptwf','idl_flag','ilbfgs_flag',
     $ 'ilbfgs_m','method','nopt_iter','ioptjas','ioptorb','ioptci',
     $ 'multiple_adiag','add_diag','ngrad_jas_blocks','nblk_max',
     $ 'nblk_ci','dl_alg','iorbprt','isample_cmat','istddev',
     $ 'limit_cmat','e_shift','save_blocks','force_blocks',
     $ 'iorbsample','iuse_trafo','iuse_orbeigv','ncore','nextorb',
     $ 'no_active','approx','approx_mix','energy_tol','sr_tau',
     $ 'sr_eps','sr_adiag','micro_iter_sr','dl_mom','lin_eps',
     $ 'lin_adiag','lin_nvec','lin_nvecx','lin_jdav','func_omega',
     $ 'omega','n_omegaf','n_omegat','sr_rescale','title','unit',
     $ 'mass','iperiodic','ibasis','nforce','nwftype','seed','ipr',
     $ 'pool','basis','pseudopot','i3dsplorb','i3dlagorb','scalecoef',
     $ 'iciprt','nctype','natom','addghostype','nghostcent'/
      data iaptr/1,3,5,6,8,13,17,23,27,28,34,39,44,63,65,72,79,88,134,
     $ 149,150/
      data ieptr/2,4,5,7,12,16,22,26,27,33,38,43,62,64,71,78,87,133,
     $ 148,149,153/
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
