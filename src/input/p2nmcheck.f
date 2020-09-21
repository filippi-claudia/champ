C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'optwf','blocking_dmc','gradients','mstates',
     $ 'forces','qmmm','3dgrid','electrons','blocking_vmc','periodic',
     $ 'jastrow','properties','vmc','strech','dmc','startend','ci',
     $ 'atoms','general','optgeo','pseudo'/
      data vars/'ioptwf','idl_flag','ilbfgs_flag','ilbfgs_m','method',
     $ 'nopt_iter','ioptjas','ioptorb','ioptci','multiple_adiag',
     $ 'add_diag','ngrad_jas_blocks','nblk_max','nblk_ci','dl_alg',
     $ 'iorbprt','isample_cmat','istddev','limit_cmat','e_shift',
     $ 'save_blocks','force_blocks','iorbsample','iuse_trafo',
     $ 'iuse_orbeigv','ncore','nextorb','no_active','approx',
     $ 'approx_mix','energy_tol','sr_tau','sr_eps','sr_adiag',
     $ 'micro_iter_sr','dl_mom','lin_eps','lin_adiag','lin_nvec',
     $ 'lin_nvecx','lin_jdav','func_omega','omega','n_omegaf',
     $ 'n_omegat','sr_rescale','nstep','nblk','nblkeq','nconf_new',
     $ 'nconf','ngradnts','igrdtype','delgrdxyz','delgrdbl',
     $ 'delgrdba','delgrdda','iguiding','iefficiency','istrech',
     $ 'alfstr','nwprod','itausec','iqmmm','stepx','stepy','stepz',
     $ 'x0','y0','z0','xn','yn','zn','nelec','nup','nstep','nblk',
     $ 'nblkeq','nconf_new','nconf','norb','npoly','np','cutg',
     $ 'cutg_sim','cutg_big','cutg_sim_big','ianalyt_lap','ijas',
     $ 'isc','nspin1','nspin2','ifock','sample','print','imetro',
     $ 'deltar','deltat','delta','fbias','node_cutoff','enode_cutoff',
     $ 'alfstr','idmc','tau','etrial','nfprod','ipq','itau_eff',
     $ 'iacc_rej','icross','icuspg','idiv_v','icut_br','icut_e',
     $ 'icasula','node_cutoff','enode_cutoff','ibranch_elec',
     $ 'icircular','idrifdifgfunc','mode_dmc','idump','irstar',
     $ 'isite','icharged_atom','iciprt','nctype','natom',
     $ 'addghostype','nghostcent','title','unit','mass','iperiodic',
     $ 'ibasis','nforce','nwftype','seed','ipr','pool','basis',
     $ 'pseudopot','i3dsplorb','i3dlagorb','scalecoef','iforce_analy',
     $ 'iuse_zmat','alfgeo','izvzb','iroot_geo','nloc','nquad'/
      data iaptr/1,47,52,58,60,64,65,74,76,81,88,94,96,103,104,123,
     $ 127,128,132,147,152/
      data ieptr/46,51,57,59,63,64,73,75,80,87,93,95,102,103,122,126,
     $ 127,131,146,151,153/
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
