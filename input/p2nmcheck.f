C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(27)*(12)
      character vars(190)*(16)
      dimension iaptr(27),ieptr(27)
      data lists/'properties','startend','rmc','optgeo','3dgrid',
     $ 'vmc','optwf','strech','jconstraint','blocking_rmc',
     $ 'blocking_vmc','atoms','pseudo','qmmm','forces','fit',
     $ 'embedding','dmc','gradients','periodic','optmethod','ci',
     $ 'mstates','jastrow','electrons','general','blocking_dmc'/
      data vars/'sample','print','idump','irstar','isite',
     $ 'icharged_atom','tau','ntau','etrial','if_dmc_rej',
     $ 'ibranch_iel','istep_quad','icasula','if_gsym','iforce_analy',
     $ 'iuse_zmat','alfgeo','izvzb','stepx','stepy','stepz','x0','y0',
     $ 'z0','xn','yn','zn','imetro','deltar','deltat','delta','fbias',
     $ 'node_cutoff','enode_cutoff','ioptwf','idl_flag','method',
     $ 'nopt_iter','ioptjas','ioptorb','ioptci','multiple_adiag',
     $ 'add_diag','ngrad_jas_blocks','nblk_max','nblk_ci','dl_alg',
     $ 'iorbprt','isample_cmat','istddev','limit_cmat','e_shift',
     $ 'save_blocks','force_blocks','iorbsample','iuse_trafo',
     $ 'iuse_orbeigv','ncore','nextorb','no_active','approx',
     $ 'approx_mix','energy_tol','sr_tau','sr_eps','sr_adiag',
     $ 'micro_iter_sr','dl_mom','lin_eps','lin_adiag','lin_nvec',
     $ 'lin_nvecx','lin_jdav','func_omega','omega','alfstr','ipos',
     $ 'idcds','idcdr','idcdt','id2cds','id2cdr','id2cdt','idbds',
     $ 'idbdr','idbdt','nstep','nblk','nblkeq','nstep','nblk',
     $ 'nblkeq','nconf_new','nconf','nctype','natom','addghostype',
     $ 'nghostcent','nloc','nquad','iqmmm','istrech','alfstr',
     $ 'nwprod','itausec','ndata','nparm','eguess','iopt','icusp',
     $ 'icusp2','ipr','irewgt','iaver','istrch','ihtree','ihpol',
     $ 'ixchng','ixpol','idmc','tau','etrial','nfprod','ipq',
     $ 'itau_eff','iacc_rej','icross','icuspg','idiv_v','icut_br',
     $ 'icut_e','icasula','node_cutoff','enode_cutoff','ibranch_elec',
     $ 'icircular','idrifdifgfunc','ngradnts','igrdtype','delgrdxyz',
     $ 'delgrdbl','delgrdba','delgrdda','norb','npoly','np','cutg',
     $ 'cutg_sim','cutg_big','cutg_sim_big','nsig','ncalls',
     $ 'pmarquardt','nstep','ibold','tau','noutput','analytic',
     $ 'cholesky','iciprt','iguiding','iefficiency','ianalyt_lap',
     $ 'ijas','isc','nspin1','nspin2','ifock','nelec','nup','title',
     $ 'unit','mass','iperiodic','ibasis','nforce','nwftype','seed',
     $ 'ipr','pool','basis','pseudopot','i3dsplorb','i3dlagorb',
     $ 'scalecoef','nstep','nblk','nblkeq','nconf_new','nconf'/
      data iaptr/1,3,7,15,19,28,35,76,77,87,90,95,99,101,102,106,116,
     $ 120,138,144,151,160,161,163,169,171,186/
      data ieptr/2,6,14,18,27,34,75,76,86,89,94,98,100,101,105,115,
     $ 119,137,143,150,159,160,162,168,170,185,190/
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
