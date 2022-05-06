C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(154)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'vmc','periodic','forces','optgeo','jastrow','dmc',
     $ 'blocking_vmc','atoms','gradients','optwf','blocking_dmc',
     $ 'general','3dgrid','electrons','qmmm','properties','strech',
     $ 'startend','mstates','pseudo','ci'/
      data vars/'imetro','deltar','deltat','delta','fbias',
     $ 'node_cutoff','enode_cutoff','norb','npoly','np','cutg',
     $ 'cutg_sim','cutg_big','cutg_sim_big','istrech','alfstr',
     $ 'nwprod','itausec','iforce_analy','iuse_zmat','alfgeo','izvzb',
     $ 'iroot_geo','ianalyt_lap','ijas','isc','nspin1','nspin2',
     $ 'ifock','idmc','tau','etrial','nfprod','ipq','itau_eff',
     $ 'iacc_rej','icross','icuspg','idiv_v','icut_br','icut_e',
     $ 'icasula','node_cutoff','enode_cutoff','ibranch_elec',
     $ 'icircular','idrifdifgfunc','mode_dmc','nstep','nblk','nblkeq',
     $ 'nconf_new','nconf','nctype','natom','addghostype',
     $ 'nghostcent','ngradnts','igrdtype','delgrdxyz','delgrdbl',
     $ 'delgrdba','delgrdda','ioptwf','idl_flag','ilbfgs_flag',
     $ 'ilbfgs_m','method','nopt_iter','ioptjas','ioptorb','ioptci',
     $ 'multiple_adiag','add_diag','ngrad_jas_blocks','nblk_max',
     $ 'nblk_ci','dl_alg','iorbprt','isample_cmat','istddev',
     $ 'limit_cmat','e_shift','save_blocks','force_blocks',
     $ 'iorbsample','iuse_trafo','iuse_orbeigv','ncore','nextorb',
     $ 'no_active','approx','approx_mix','energy_tol','sr_tau',
     $ 'sr_eps','sr_adiag','alambda','micro_iter_sr','dl_mom',
     $ 'lin_eps','lin_adiag','lin_nvec','lin_nvecx','lin_jdav',
     $ 'func_omega','omega','n_omegaf','n_omegat','sr_rescale',
     $ 'nstep','nblk','nblkeq','nconf_new','nconf','title','unit',
     $ 'mass','iperiodic','ibasis','nforce','nwftype','seed','ipr',
     $ 'pool','basis','pseudopot','i3dsplorb','i3dlagorb','scalecoef',
     $ 'stepx','stepy','stepz','x0','y0','z0','xn','yn','zn','nelec',
     $ 'nup','iqmmm','sample','print','alfstr','idump','irstar',
     $ 'isite','icharged_atom','iguiding','iefficiency','nloc',
     $ 'nquad','iciprt'/
      data iaptr/1,8,15,19,24,30,49,54,58,64,111,116,131,140,142,143,
     $ 145,146,150,152,154/
      data ieptr/7,14,18,23,29,48,53,57,63,110,115,130,139,141,142,
     $ 144,145,149,151,153,154/
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
