C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'optwf','strech','electrons','pseudo','startend',
     $ '3dgrid','atoms','periodic','forces','mstates','vmc','general',
     $ 'gradients','qmmm','jastrow','blocking_dmc','optgeo',
     $ 'blocking_vmc','dmc','properties','ci'/
      data vars/'ioptwf','idl_flag','ilbfgs_flag','ilbfgs_m','method',
     $ 'nopt_iter','ioptjas','ioptorb','ioptci','multiple_adiag',
     $ 'add_diag','ngrad_jas_blocks','nblk_max','nblk_ci','dl_alg',
     $ 'iorbprt','isample_cmat','istddev','limit_cmat','e_shift',
     $ 'save_blocks','force_blocks','iorbsample','iuse_trafo',
     $ 'iuse_orbeigv','ncore','nextorb','no_active','approx',
     $ 'approx_mix','energy_tol','sr_tau','sr_eps','sr_adiag',
     $ 'micro_iter_sr','dl_mom','lin_eps','lin_adiag','lin_nvec',
     $ 'lin_nvecx','lin_jdav','func_omega','omega','n_omegaf',
     $ 'n_omegat','sr_rescale','alfstr','nelec','nup','nloc','nquad',
     $ 'idump','irstar','isite','icharged_atom','stepx','stepy',
     $ 'stepz','x0','y0','z0','xn','yn','zn','nctype','natom',
     $ 'addghostype','nghostcent','norb','npoly','np','cutg',
     $ 'cutg_sim','cutg_big','cutg_sim_big','istrech','alfstr',
     $ 'nwprod','itausec','iguiding','iefficiency','imetro','deltar',
     $ 'deltat','delta','fbias','node_cutoff','enode_cutoff','title',
     $ 'unit','mass','iperiodic','ibasis','nforce','nwftype','seed',
     $ 'ipr','pool','basis','pseudopot','i3dsplorb','i3dlagorb',
     $ 'scalecoef','ngradnts','igrdtype','delgrdxyz','delgrdbl',
     $ 'delgrdba','delgrdda','iqmmm','ianalyt_lap','ijas','isc',
     $ 'nspin1','nspin2','ifock','nstep','nblk','nblkeq','nconf_new',
     $ 'nconf','iforce_analy','iuse_zmat','alfgeo','izvzb',
     $ 'iroot_geo','nstep','nblk','nblkeq','nconf_new','nconf','idmc',
     $ 'tau','etrial','nfprod','ipq','itau_eff','iacc_rej','icross',
     $ 'icuspg','idiv_v','icut_br','icut_e','icasula','node_cutoff',
     $ 'enode_cutoff','ibranch_elec','icircular','idrifdifgfunc',
     $ 'mode_dmc','sample','print','iciprt'/
      data iaptr/1,47,48,50,52,56,65,69,76,80,82,89,104,110,111,117,
     $ 122,127,132,151,153/
      data ieptr/46,47,49,51,55,64,68,75,79,81,88,103,109,110,116,121,
     $ 126,131,150,152,153/
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
