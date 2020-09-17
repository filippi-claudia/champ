C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'optgeo','ci','general','3dgrid','startend','vmc',
     $ 'gradients','electrons','blocking_dmc','blocking_vmc','atoms',
     $ 'periodic','mstates','dmc','forces','jastrow','properties',
     $ 'strech','pseudo','optwf','qmmm'/
      data vars/'iforce_analy','iuse_zmat','alfgeo','izvzb',
     $ 'iroot_geo','iciprt','title','unit','mass','iperiodic',
     $ 'ibasis','nforce','nwftype','seed','ipr','pool','basis',
     $ 'pseudopot','i3dsplorb','i3dlagorb','scalecoef','stepx',
     $ 'stepy','stepz','x0','y0','z0','xn','yn','zn','idump','irstar',
     $ 'isite','icharged_atom','imetro','deltar','deltat','delta',
     $ 'fbias','node_cutoff','enode_cutoff','ngradnts','igrdtype',
     $ 'delgrdxyz','delgrdbl','delgrdba','delgrdda','nelec','nup',
     $ 'nstep','nblk','nblkeq','nconf_new','nconf','nstep','nblk',
     $ 'nblkeq','nconf_new','nconf','nctype','natom','addghostype',
     $ 'nghostcent','norb','npoly','np','cutg','cutg_sim','cutg_big',
     $ 'cutg_sim_big','iguiding','iefficiency','idmc','tau','etrial',
     $ 'nfprod','ipq','itau_eff','iacc_rej','icross','icuspg',
     $ 'idiv_v','icut_br','icut_e','icasula','node_cutoff',
     $ 'enode_cutoff','ibranch_elec','icircular','idrifdifgfunc',
     $ 'mode_dmc','istrech','alfstr','nwprod','itausec','ianalyt_lap',
     $ 'ijas','isc','nspin1','nspin2','ifock','sample','print',
     $ 'alfstr','nloc','nquad','ioptwf','idl_flag','ilbfgs_flag',
     $ 'ilbfgs_m','method','nopt_iter','ioptjas','ioptorb','ioptci',
     $ 'multiple_adiag','add_diag','ngrad_jas_blocks','nblk_max',
     $ 'nblk_ci','dl_alg','iorbprt','isample_cmat','istddev',
     $ 'limit_cmat','e_shift','save_blocks','force_blocks',
     $ 'iorbsample','iuse_trafo','iuse_orbeigv','ncore','nextorb',
     $ 'no_active','approx','approx_mix','energy_tol','sr_tau',
     $ 'sr_eps','sr_adiag','micro_iter_sr','dl_mom','lin_eps',
     $ 'lin_adiag','lin_nvec','lin_nvecx','lin_jdav','func_omega',
     $ 'omega','n_omegaf','n_omegat','sr_rescale','iqmmm'/
      data iaptr/1,6,7,22,31,35,42,48,50,55,60,64,71,73,92,96,102,104,
     $ 105,107,153/
      data ieptr/5,6,21,30,34,41,47,49,54,59,63,70,72,91,95,101,103,
     $ 104,106,152,153/
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
