C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'electrons','periodic','vmc','properties',
     $ 'blocking_vmc','general','ci','mstates','strech','jastrow',
     $ 'qmmm','blocking_dmc','dmc','startend','atoms','optgeo',
     $ 'gradients','3dgrid','forces','optwf','pseudo'/
      data vars/'nelec','nup','norb','npoly','np','cutg','cutg_sim',
     $ 'cutg_big','cutg_sim_big','imetro','deltar','deltat','delta',
     $ 'fbias','node_cutoff','enode_cutoff','sample','print','nstep',
     $ 'nblk','nblkeq','nconf_new','nconf','title','unit','mass',
     $ 'iperiodic','ibasis','nforce','nwftype','seed','ipr','pool',
     $ 'basis','pseudopot','i3dsplorb','i3dlagorb','scalecoef',
     $ 'iciprt','iguiding','iefficiency','alfstr','ianalyt_lap',
     $ 'ijas','isc','nspin1','nspin2','ifock','iqmmm','nstep','nblk',
     $ 'nblkeq','nconf_new','nconf','idmc','tau','etrial','nfprod',
     $ 'ipq','itau_eff','iacc_rej','icross','icuspg','idiv_v',
     $ 'icut_br','icut_e','icasula','node_cutoff','enode_cutoff',
     $ 'ibranch_elec','icircular','idrifdifgfunc','mode_dmc','idump',
     $ 'irstar','isite','icharged_atom','nctype','natom',
     $ 'addghostype','nghostcent','iforce_analy','iuse_zmat','alfgeo',
     $ 'izvzb','iroot_geo','ngradnts','igrdtype','delgrdxyz',
     $ 'delgrdbl','delgrdba','delgrdda','stepx','stepy','stepz','x0',
     $ 'y0','z0','xn','yn','zn','istrech','alfstr','nwprod','itausec',
     $ 'ioptwf','idl_flag','ilbfgs_flag','ilbfgs_m','method',
     $ 'nopt_iter','ioptjas','ioptorb','ioptci','multiple_adiag',
     $ 'add_diag','ngrad_jas_blocks','nblk_max','nblk_ci','dl_alg',
     $ 'iorbprt','isample_cmat','istddev','limit_cmat','e_shift',
     $ 'save_blocks','force_blocks','iorbsample','iuse_trafo',
     $ 'iuse_orbeigv','ncore','nextorb','no_active','approx',
     $ 'approx_mix','energy_tol','sr_tau','sr_eps','sr_adiag',
     $ 'micro_iter_sr','dl_mom','lin_eps','lin_adiag','lin_nvec',
     $ 'lin_nvecx','lin_jdav','func_omega','omega','n_omegaf',
     $ 'n_omegat','sr_rescale','nloc','nquad'/
      data iaptr/1,3,10,17,19,24,39,40,42,43,49,50,55,74,78,82,87,93,
     $ 102,106,152/
      data ieptr/2,9,16,18,23,38,39,41,42,48,49,54,73,77,81,86,92,101,
     $ 105,151,153/
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
