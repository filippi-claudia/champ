C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(27)*(12)
      character vars(190)*(16)
      dimension iaptr(27),ieptr(27)
      data lists/'blocking_vmc','gradients','blocking_rmc','startend',
     $ 'mstates','3dgrid','general','optgeo','dmc','ci','jconstraint',
     $ 'electrons','periodic','blocking_dmc','embedding','forces',
     $ 'fit','properties','strech','atoms','jastrow','optwf',
     $ 'optmethod','vmc','pseudo','rmc','qmmm'/
      data vars/'nstep','nblk','nblkeq','nconf_new','nconf',
     $ 'ngradnts','igrdtype','delgrdxyz','delgrdbl','delgrdba',
     $ 'delgrdda','nstep','nblk','nblkeq','idump','irstar','isite',
     $ 'icharged_atom','iguiding','iefficiency','stepx','stepy',
     $ 'stepz','x0','y0','z0','xn','yn','zn','title','unit','mass',
     $ 'iperiodic','ibasis','nforce','nwftype','seed','ipr','pool',
     $ 'basis','pseudopot','i3dsplorb','i3dlagorb','scalecoef',
     $ 'iforce_analy','iuse_zmat','alfgeo','izvzb','idmc','tau',
     $ 'etrial','nfprod','ipq','itau_eff','iacc_rej','icross',
     $ 'icuspg','idiv_v','icut_br','icut_e','icasula','node_cutoff',
     $ 'enode_cutoff','ibranch_elec','icircular','idrifdifgfunc',
     $ 'iciprt','ipos','idcds','idcdr','idcdt','id2cds','id2cdr',
     $ 'id2cdt','idbds','idbdr','idbdt','nelec','nup','norb','npoly',
     $ 'np','cutg','cutg_sim','cutg_big','cutg_sim_big','nstep',
     $ 'nblk','nblkeq','nconf_new','nconf','ihtree','ihpol','ixchng',
     $ 'ixpol','istrech','alfstr','nwprod','itausec','ndata','nparm',
     $ 'eguess','iopt','icusp','icusp2','ipr','irewgt','iaver',
     $ 'istrch','sample','print','alfstr','nctype','natom',
     $ 'addghostype','nghostcent','ianalyt_lap','ijas','isc','nspin1',
     $ 'nspin2','ifock','ioptwf','idl_flag','method','nopt_iter',
     $ 'ioptjas','ioptorb','ioptci','multiple_adiag','add_diag',
     $ 'ngrad_jas_blocks','nblk_max','nblk_ci','dl_alg','iorbprt',
     $ 'isample_cmat','istddev','limit_cmat','e_shift','save_blocks',
     $ 'force_blocks','iorbsample','iuse_trafo','iuse_orbeigv',
     $ 'ncore','nextorb','no_active','approx','approx_mix',
     $ 'energy_tol','sr_tau','sr_eps','sr_adiag','micro_iter_sr',
     $ 'dl_mom','lin_eps','lin_adiag','lin_nvec','lin_nvecx',
     $ 'lin_jdav','func_omega','omega','nsig','ncalls','pmarquardt',
     $ 'nstep','ibold','tau','noutput','analytic','cholesky','imetro',
     $ 'deltar','deltat','delta','fbias','node_cutoff','enode_cutoff',
     $ 'nloc','nquad','tau','ntau','etrial','if_dmc_rej',
     $ 'ibranch_iel','istep_quad','icasula','if_gsym','iqmmm'/
      data iaptr/1,6,12,15,19,21,30,45,49,67,68,78,80,87,92,96,100,
     $ 110,112,113,117,123,164,173,180,182,190/
      data ieptr/5,11,14,18,20,29,44,48,66,67,77,79,86,91,95,99,109,
     $ 111,112,116,122,163,172,179,181,189,190/
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
