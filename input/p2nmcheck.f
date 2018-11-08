C ---------- predefined namlist check -------
C this file is auto generated, do not edit   
      subroutine p2nmcheck(p,v,ierr)
      character p*(*),v*(*)
      character lists(25)*(12)
      character vars(151)*(16)
      dimension iaptr(25),ieptr(25)
      data lists/'jastrow','general','optmethod','strech','qmmm',
     $ 'embedding','3dgrid','periodic','forces','ci','pseudo',
     $ 'blocking_vmc','optwf','dmc','atoms','blocking_dmc',
     $ 'electrons','fit','rmc','startend','properties','jconstraint',
     $ 'blocking_rmc','vmc','mstates'/
      data vars/'ianalyt_lap','ijas','isc','nspin1','nspin2','ifock',
     $ 'title','unit','mass','iperiodic','ibasis','nforce','nwftype',
     $ 'seed','ipr','pool','basis','pseudopot','i3dsplorb',
     $ 'i3dlagorb','nsig','ncalls','pmarquardt','nstep','ibold','tau',
     $ 'noutput','analytic','cholesky','alfstr','iqmmm','ihtree',
     $ 'ihpol','ixchng','ixpol','stepx','stepy','stepz','x0','y0',
     $ 'z0','xn','yn','zn','norb','npoly','np','cutg','cutg_sim',
     $ 'cutg_big','cutg_sim_big','istrech','alfstr','nwprod',
     $ 'itausec','iciprt','nloc','nquad','nstep','nblk','nblkeq',
     $ 'nconf_new','nconf','ioptwf','method','nopt_iter','ioptjas',
     $ 'ioptorb','ioptci','multiple_adiag','add_diag',
     $ 'ngrad_jas_blocks','nblk_max','nblk_ci','iorbprt',
     $ 'isample_cmat','istddev','limit_cmat','e_shift','save_blocks',
     $ 'force_blocks','iorbsample','iuse_trafo','iuse_orbeigv',
     $ 'ncore','nextorb','no_active','approx','approx_mix','idmc',
     $ 'tau','etrial','nfprod','ipq','itau_eff','iacc_rej','icross',
     $ 'icuspg','idiv_v','icut_br','icut_e','icasula','nctype',
     $ 'natom','addghostype','nghostcent','nstep','nblk','nblkeq',
     $ 'nconf_new','nconf','nelec','nup','ndata','nparm','eguess',
     $ 'iopt','icusp','icusp2','ipr','irewgt','iaver','istrch','tau',
     $ 'ntau','etrial','idump','irstar','isite','sample','print',
     $ 'ipos','idcds','idcdr','idcdt','id2cds','id2cdr','id2cdt',
     $ 'idbds','idbdr','idbdt','nstep','nblk','nblkeq','imetro',
     $ 'deltar','deltat','delta','fbias','iguiding','iefficiency'/
      data iaptr/1,7,21,30,31,32,36,45,52,56,57,59,64,90,103,107,112,
     $ 114,124,127,130,132,142,145,150/
      data ieptr/6,20,29,30,31,35,44,51,55,56,58,63,89,102,106,111,
     $ 113,123,126,129,131,141,144,149,151/
      nlist=25
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
