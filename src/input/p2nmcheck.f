C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'mstates','dmc','startend','blocking_dmc','3dgrid',
     $ 'atoms','blocking_vmc','ci','gradients','properties','vmc',
     $ 'strech','jastrow','general','pseudo','electrons','qmmm',
     $ 'optwf','optgeo','forces','periodic'/
      data vars/'iguiding','iefficiency','idmc','tau','etrial',
     $ 'nfprod','ipq','itau_eff','iacc_rej','icross','icuspg',
     $ 'idiv_v','icut_br','icut_e','icasula','node_cutoff',
     $ 'enode_cutoff','ibranch_elec','icircular','idrifdifgfunc',
     $ 'mode_dmc','idump','irstar','isite','icharged_atom','nstep',
     $ 'nblk','nblkeq','nconf_new','nconf','stepx','stepy','stepz',
     $ 'x0','y0','z0','xn','yn','zn','nctype','natom','addghostype',
     $ 'nghostcent','nstep','nblk','nblkeq','nconf_new','nconf',
     $ 'iciprt','ngradnts','igrdtype','delgrdxyz','delgrdbl',
     $ 'delgrdba','delgrdda','sample','print','imetro','deltar',
     $ 'deltat','delta','fbias','node_cutoff','enode_cutoff','alfstr',
     $ 'ianalyt_lap','ijas','isc','nspin1','nspin2','ifock','title',
     $ 'unit','mass','iperiodic','ibasis','nforce','nwftype','seed',
     $ 'ipr','pool','basis','pseudopot','i3dsplorb','i3dlagorb',
     $ 'scalecoef','nloc','nquad','nelec','nup','iqmmm','ioptwf',
     $ 'idl_flag','ilbfgs_flag','ilbfgs_m','method','nopt_iter',
     $ 'ioptjas','ioptorb','ioptci','multiple_adiag','add_diag',
     $ 'ngrad_jas_blocks','nblk_max','nblk_ci','dl_alg','iorbprt',
     $ 'isample_cmat','istddev','limit_cmat','e_shift','save_blocks',
     $ 'force_blocks','iorbsample','iuse_trafo','iuse_orbeigv',
     $ 'ncore','nextorb','no_active','approx','approx_mix',
     $ 'energy_tol','sr_tau','sr_eps','sr_adiag','micro_iter_sr',
     $ 'dl_mom','lin_eps','lin_adiag','lin_nvec','lin_nvecx',
     $ 'lin_jdav','func_omega','omega','n_omegaf','n_omegat',
     $ 'sr_rescale','iforce_analy','iuse_zmat','alfgeo','izvzb',
     $ 'iroot_geo','istrech','alfstr','nwprod','itausec','norb',
     $ 'npoly','np','cutg','cutg_sim','cutg_big','cutg_sim_big'/
      data iaptr/1,3,22,26,31,40,44,49,50,56,58,65,66,72,87,89,91,92,
     $ 138,143,147/
      data ieptr/2,21,25,30,39,43,48,49,55,57,64,65,71,86,88,90,91,
     $ 137,142,146,153/
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
