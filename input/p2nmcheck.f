C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(27)*(12)
      character vars(186)*(16)
      dimension iaptr(27),ieptr(27)
      data lists/'jastrow','general','optmethod','strech','qmmm',
     $ 'optgeo','embedding','3dgrid','periodic','forces','ci',
     $ 'pseudo','blocking_vmc','optwf','dmc','atoms','blocking_dmc',
     $ 'electrons','fit','rmc','startend','properties','jconstraint',
     $ 'gradients','blocking_rmc','vmc','mstates'/
      data vars/'ianalyt_lap','ijas','isc','nspin1','nspin2','ifock',
     $ 'title','unit','mass','iperiodic','ibasis','nforce','nwftype',
     $ 'seed','ipr','pool','basis','pseudopot','i3dsplorb',
     $ 'i3dlagorb','scalecoef','nsig','ncalls','pmarquardt','nstep',
     $ 'ibold','tau','noutput','analytic','cholesky','alfstr','iqmmm',
     $ 'iforce_analy','iuse_zmat','alfgeo','izvzb','iroot_geo',
     $ 'ihtree','ihpol','ixchng','ixpol','stepx','stepy','stepz','x0',
     $ 'y0','z0','xn','yn','zn','norb','npoly','np','cutg','cutg_sim',
     $ 'cutg_big','cutg_sim_big','istrech','alfstr','nwprod',
     $ 'itausec','iciprt','nloc','nquad','nstep','nblk','nblkeq',
     $ 'nconf_new','nconf','ioptwf','method','nopt_iter','ioptjas',
     $ 'ioptorb','ioptci','multiple_adiag','add_diag',
     $ 'ngrad_jas_blocks','nblk_max','nblk_ci','iorbprt',
     $ 'isample_cmat','istddev','limit_cmat','e_shift','save_blocks',
     $ 'force_blocks','iorbsample','iuse_trafo','iuse_orbeigv',
     $ 'ncore','nextorb','no_active','approx','approx_mix',
     $ 'energy_tol','sr_tau','sr_eps','sr_adiag','micro_iter_sr',
     $ 'lin_eps','lin_adiag','lin_nvec','lin_nvecx','lin_jdav','idmc',
     $ 'tau','etrial','nfprod','ipq','itau_eff','iacc_rej','icross',
     $ 'icuspg','idiv_v','icut_br','icut_e','icasula','node_cutoff',
     $ 'enode_cutoff','ibranch_elec','icircular','idrifdifgfunc',
     $ 'nctype','natom','addghostype','nghostcent','nstep','nblk',
     $ 'nblkeq','nconf_new','nconf','nelec','nup','ndata','nparm',
     $ 'eguess','iopt','icusp','icusp2','ipr','irewgt','iaver',
     $ 'istrch','tau','ntau','etrial','if_dmc_rej','ibranch_iel',
     $ 'istep_quad','icasula','if_gsym','idump','irstar','isite',
     $ 'icharged_atom','sample','print','ipos','idcds','idcdr',
     $ 'idcdt','id2cds','id2cdr','id2cdt','idbds','idbdr','idbdt',
     $ 'ngradnts','igrdtype','delgrdxyz','delgrdbl','delgrdba',
     $ 'delgrdda','nstep','nblk','nblkeq','imetro','deltar','deltat',
     $ 'delta','fbias','node_cutoff','enode_cutoff','iguiding',
     $ 'iefficiency'/
      data iaptr/1,7,22,31,32,33,38,42,51,58,62,63,65,70,106,124,128,
     $ 133,135,145,153,157,159,169,175,178,185/
      data ieptr/6,21,30,31,32,37,41,50,57,61,62,64,69,105,123,127,
     $ 132,134,144,152,156,158,168,174,177,184,186/
      nlist=27
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
