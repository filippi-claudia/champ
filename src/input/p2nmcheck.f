C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'qmmm','forces','dmc','pseudo','ci','atoms','vmc',
     $ 'gradients','optwf','mstates','periodic','general','optgeo',
     $ 'strech','blocking_dmc','jastrow','startend','electrons',
     $ '3dgrid','properties','blocking_vmc'/
      data vars/'iqmmm','istrech','alfstr','nwprod','itausec','idmc',
     $ 'tau','etrial','nfprod','ipq','itau_eff','iacc_rej','icross',
     $ 'icuspg','idiv_v','icut_br','icut_e','icasula','node_cutoff',
     $ 'enode_cutoff','ibranch_elec','icircular','idrifdifgfunc',
     $ 'mode_dmc','nloc','nquad','iciprt','nctype','natom',
     $ 'addghostype','nghostcent','imetro','deltar','deltat','delta',
     $ 'fbias','node_cutoff','enode_cutoff','ngradnts','igrdtype',
     $ 'delgrdxyz','delgrdbl','delgrdba','delgrdda','ioptwf',
     $ 'idl_flag','ilbfgs_flag','ilbfgs_m','method','nopt_iter',
     $ 'ioptjas','ioptorb','ioptci','multiple_adiag','add_diag',
     $ 'ngrad_jas_blocks','nblk_max','nblk_ci','dl_alg','iorbprt',
     $ 'isample_cmat','istddev','limit_cmat','e_shift','save_blocks',
     $ 'force_blocks','iorbsample','iuse_trafo','iuse_orbeigv',
     $ 'ncore','nextorb','no_active','approx','approx_mix',
     $ 'energy_tol','sr_tau','sr_eps','sr_adiag','micro_iter_sr',
     $ 'dl_mom','lin_eps','lin_adiag','lin_nvec','lin_nvecx',
     $ 'lin_jdav','func_omega','omega','n_omegaf','n_omegat',
     $ 'sr_rescale','iguiding','iefficiency','norb','npoly','np',
     $ 'cutg','cutg_sim','cutg_big','cutg_sim_big','title','unit',
     $ 'mass','iperiodic','ibasis','nforce','nwftype','seed','ipr',
     $ 'pool','basis','pseudopot','i3dsplorb','i3dlagorb','scalecoef',
     $ 'iforce_analy','iuse_zmat','alfgeo','izvzb','iroot_geo',
     $ 'alfstr','nstep','nblk','nblkeq','nconf_new','nconf',
     $ 'ianalyt_lap','ijas','isc','nspin1','nspin2','ifock','idump',
     $ 'irstar','isite','icharged_atom','nelec','nup','stepx','stepy',
     $ 'stepz','x0','y0','z0','xn','yn','zn','sample','print','nstep',
     $ 'nblk','nblkeq','nconf_new','nconf'/
      data iaptr/1,2,6,25,27,28,32,39,45,91,93,100,115,120,121,126,
     $ 132,136,138,147,149/
      data ieptr/1,5,24,26,27,31,38,44,90,92,99,114,119,120,125,131,
     $ 135,137,146,148,153/
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
