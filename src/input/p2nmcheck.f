C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'optwf','vmc','qmmm','forces','pseudo','gradients',
     $ 'dmc','startend','general','atoms','strech','blocking_vmc',
     $ 'optgeo','jastrow','properties','ci','3dgrid','electrons',
     $ 'mstates','periodic','blocking_dmc'/
      data vars/'ioptwf','idl_flag','ilbfgs_flag','ilbfgs_m','method',
     $ 'nopt_iter','ioptjas','ioptorb','ioptci','multiple_adiag',
     $ 'add_diag','ngrad_jas_blocks','nblk_max','nblk_ci','dl_alg',
     $ 'iorbprt','isample_cmat','istddev','limit_cmat','e_shift',
     $ 'save_blocks','force_blocks','iorbsample','iuse_trafo',
     $ 'iuse_orbeigv','ncore','nextorb','no_active','approx',
     $ 'approx_mix','energy_tol','sr_tau','sr_eps','sr_adiag',
     $ 'micro_iter_sr','dl_mom','lin_eps','lin_adiag','lin_nvec',
     $ 'lin_nvecx','lin_jdav','func_omega','omega','n_omegaf',
     $ 'n_omegat','sr_rescale','imetro','deltar','deltat','delta',
     $ 'fbias','node_cutoff','enode_cutoff','iqmmm','istrech',
     $ 'alfstr','nwprod','itausec','nloc','nquad','ngradnts',
     $ 'igrdtype','delgrdxyz','delgrdbl','delgrdba','delgrdda','idmc',
     $ 'tau','etrial','nfprod','ipq','itau_eff','iacc_rej','icross',
     $ 'icuspg','idiv_v','icut_br','icut_e','icasula','node_cutoff',
     $ 'enode_cutoff','ibranch_elec','icircular','idrifdifgfunc',
     $ 'mode_dmc','idump','irstar','isite','icharged_atom','title',
     $ 'unit','mass','iperiodic','ibasis','nforce','nwftype','seed',
     $ 'ipr','pool','basis','pseudopot','i3dsplorb','i3dlagorb',
     $ 'scalecoef','nctype','natom','addghostype','nghostcent',
     $ 'alfstr','nstep','nblk','nblkeq','nconf_new','nconf',
     $ 'iforce_analy','iuse_zmat','alfgeo','izvzb','iroot_geo',
     $ 'ianalyt_lap','ijas','isc','nspin1','nspin2','ifock','sample',
     $ 'print','iciprt','stepx','stepy','stepz','x0','y0','z0','xn',
     $ 'yn','zn','nelec','nup','iguiding','iefficiency','norb',
     $ 'npoly','np','cutg','cutg_sim','cutg_big','cutg_sim_big',
     $ 'nstep','nblk','nblkeq','nconf_new','nconf'/
      data iaptr/1,47,54,55,59,61,67,86,90,105,109,110,115,120,126,
     $ 128,129,138,140,142,149/
      data ieptr/46,53,54,58,60,66,85,89,104,108,109,114,119,125,127,
     $ 128,137,139,141,148,153/
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
