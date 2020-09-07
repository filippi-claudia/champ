C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(21)*(12)
      character vars(153)*(16)
      dimension iaptr(21),ieptr(21)
      data lists/'jastrow','optwf','blocking_vmc','ci','strech',
     $ 'gradients','forces','pseudo','3dgrid','properties','mstates',
     $ 'periodic','general','qmmm','startend','optgeo','vmc','dmc',
     $ 'blocking_dmc','electrons','atoms'/
      data vars/'ianalyt_lap','ijas','isc','nspin1','nspin2','ifock',
     $ 'ioptwf','idl_flag','ilbfgs_flag','ilbfgs_m','method',
     $ 'nopt_iter','ioptjas','ioptorb','ioptci','multiple_adiag',
     $ 'add_diag','ngrad_jas_blocks','nblk_max','nblk_ci','dl_alg',
     $ 'iorbprt','isample_cmat','istddev','limit_cmat','e_shift',
     $ 'save_blocks','force_blocks','iorbsample','iuse_trafo',
     $ 'iuse_orbeigv','ncore','nextorb','no_active','approx',
     $ 'approx_mix','energy_tol','sr_tau','sr_eps','sr_adiag',
     $ 'micro_iter_sr','dl_mom','lin_eps','lin_adiag','lin_nvec',
     $ 'lin_nvecx','lin_jdav','func_omega','omega','n_omegaf',
     $ 'n_omegat','sr_rescale','nstep','nblk','nblkeq','nconf_new',
     $ 'nconf','iciprt','alfstr','ngradnts','igrdtype','delgrdxyz',
     $ 'delgrdbl','delgrdba','delgrdda','istrech','alfstr','nwprod',
     $ 'itausec','nloc','nquad','stepx','stepy','stepz','x0','y0',
     $ 'z0','xn','yn','zn','sample','print','iguiding','iefficiency',
     $ 'norb','npoly','np','cutg','cutg_sim','cutg_big',
     $ 'cutg_sim_big','title','unit','mass','iperiodic','ibasis',
     $ 'nforce','nwftype','seed','ipr','pool','basis','pseudopot',
     $ 'i3dsplorb','i3dlagorb','scalecoef','iqmmm','idump','irstar',
     $ 'isite','icharged_atom','iforce_analy','iuse_zmat','alfgeo',
     $ 'izvzb','iroot_geo','imetro','deltar','deltat','delta','fbias',
     $ 'node_cutoff','enode_cutoff','idmc','tau','etrial','nfprod',
     $ 'ipq','itau_eff','iacc_rej','icross','icuspg','idiv_v',
     $ 'icut_br','icut_e','icasula','node_cutoff','enode_cutoff',
     $ 'ibranch_elec','icircular','idrifdifgfunc','mode_dmc','nstep',
     $ 'nblk','nblkeq','nconf_new','nconf','nelec','nup','nctype',
     $ 'natom','addghostype','nghostcent'/
      data iaptr/1,7,53,58,59,60,66,70,72,81,83,85,92,107,108,112,117,
     $ 124,143,148,150/
      data ieptr/6,52,57,58,59,65,69,71,80,82,84,91,106,107,111,116,
     $ 123,142,147,149,153/
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
