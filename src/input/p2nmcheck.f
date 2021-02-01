C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'dmc','vmc','periodic','properties','forces','atoms',
     $ 'mstates','electrons','optwf','jastrow','pseudo','3dgrid',
     $ 'gradients','qmmm','general','startend','strech',
     $ 'blocking_vmc','blocking_dmc','optgeo','ci'/
      data vars/'idmc','tau','etrial','nfprod','ipq','itau_eff',
     $ 'iacc_rej','icross','icuspg','idiv_v','icut_br','icut_e',
     $ 'icasula','node_cutoff','enode_cutoff','ibranch_elec',
     $ 'icircular','idrifdifgfunc','mode_dmc','imetro','deltar',
     $ 'deltat','delta','fbias','node_cutoff','enode_cutoff','norb',
     $ 'npoly','np','cutg','cutg_sim','cutg_big','cutg_sim_big',
     $ 'sample','print','istrech','alfstr','nwprod','itausec',
     $ 'nctype','natom','addghostype','nghostcent','iguiding',
     $ 'iefficiency','nelec','nup','ioptwf','idl_flag','ilbfgs_flag',
     $ 'ilbfgs_m','method','nopt_iter','ioptjas','ioptorb','ioptci',
     $ 'multiple_adiag','add_diag','ngrad_jas_blocks','nblk_max',
     $ 'nblk_ci','dl_alg','iorbprt','isample_cmat','istddev',
     $ 'limit_cmat','e_shift','save_blocks','force_blocks',
     $ 'iorbsample','iuse_trafo','iuse_orbeigv','ncore','nextorb',
     $ 'no_active','approx','approx_mix','energy_tol','sr_tau',
     $ 'sr_eps','sr_adiag','micro_iter_sr','dl_mom','lin_eps',
     $ 'lin_adiag','lin_nvec','lin_nvecx','lin_jdav','func_omega',
     $ 'omega','n_omegaf','n_omegat','sr_rescale','ianalyt_lap',
     $ 'ijas','isc','nspin1','nspin2','ifock','nloc','nquad','stepx',
     $ 'stepy','stepz','x0','y0','z0','xn','yn','zn','ngradnts',
     $ 'igrdtype','delgrdxyz','delgrdbl','delgrdba','delgrdda',
     $ 'iqmmm','title','unit','mass','iperiodic','ibasis','nforce',
     $ 'nwftype','seed','ipr','pool','basis','pseudopot','i3dsplorb',
     $ 'i3dlagorb','scalecoef','idump','irstar','isite',
     $ 'icharged_atom','alfstr','nstep','nblk','nblkeq','nconf_new',
     $ 'nconf','nstep','nblk','nblkeq','nconf_new','nconf',
     $ 'iforce_analy','iuse_zmat','alfgeo','izvzb','iroot_geo',
     $ 'iciprt'/
      data iaptr/1,20,27,34,36,40,44,46,48,94,100,102,111,117,118,133,
     $ 137,138,143,148,153/
      data ieptr/19,26,33,35,39,43,45,47,93,99,101,110,116,117,132,
     $ 136,137,142,147,152,153/
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
