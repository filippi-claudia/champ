      subroutine optwf_matrix_corsamp
c written by Claudia Filippi

      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use forcepar, only: deltot, istrech, nforce
      use numbas, only: arg, d2rwf, igrid, iwrwf, nr, nrbas, numr, r0, rwf

      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      use optwf_corsam, only: add_diag, add_diag_tmp, energy, energy_err, force, force_err
      use optwf_parms, only: nparmd, nparme, nparmg, nparmj, nparml, nparms
      use wfsec, only: iwf, iwftype, nwftype
      use contrl, only: idump, irstar, isite, n_conf, nblk, nblkeq, nconf_new, nstep
      implicit real*8(a-h,o-z)








      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optjas.h'
      include 'optci.h'
      include 'optorb.h'
      include 'numbas.h'

      parameter(MPARMALL=MPARMJ+MXCIREDUCED+MXREDUCED,MPARMALL2=MPARMALL*(MPARMALL+1)/2)
      parameter(MWORK=50*MPARMALL)





      common /gradhess_all/ grad(MPARMALL),h(MPARMALL,MPARMALL),s(MPARMALL,MPARMALL)



      dimension grad_sav(MPARMALL),h_sav(MPARMALL,MPARMALL),s_sav(MPARMALL2)
      dimension work(MWORK),work2(MPARMALL,MPARMALL)

      call p2gtid('optwf:ioptjas',ioptjas,0,1)
      call p2gtid('optwf:ioptorb',ioptorb,0,1)
      call p2gtid('optwf:ioptci',ioptci,0,1)

c No dump/restart if optimizing wave function
      irstar=0
      idump=0
c Set up basis functions for test run
      do 1 iwft=2,3
   1    iwftype(iwft)=iwft
      if(numr.gt.0) then
        do 2 iwft=2,3
   2      call read_bas_num(iwft)
       else
        do 3 iwft=2,3
   3      call copy_zex(iwft)
      endif
      call set_displace_zero(3)

c Number of iterations
      call p2gtid('optwf:nopt_iter',nopt_iter,6,1)
      write(6,'(/,''Number of iterations'',i2)') nopt_iter
      if(ioptci.eq.1.and.ioptjas.eq.0.and.ioptorb.eq.0) then
        nopt_iter=1
        write(6,'(''Reset number of iterations to 1'')')
      endif
c Max number of blocks
      call p2gtid('optwf:nblk_max',nblk_max,nblk,1)
      write(6,'(/,''Maximum number of blocks'',i4)') nblk_max
c Compute multiple adiag
      call p2gtid('optwf:multiple_adiag',multiple_adiag,0,1)
      write(6,'(/,''Perform test run with multiple adiag'',i2)') multiple_adiag
c Tolerance on energy
      call p2gtfd('optwf:energy_tol',energy_tol,1.d-3,1)
      write(6,'(/,''Energy tolerance'',d12.2)') energy_tol
  
      call p2gtfd('optwf:add_diag',add_diag(1),1.d-6,1)
      if(ioptjas.eq.0.and.ioptorb.eq.0) add_diag(1)=-1

c Set dparm_norm_min
      call p2gtfd('optwf:dparm_norm_min',dparm_norm_min,1.0d0,1)
      write(6,'(''Starting dparm_norm_min'',g12.4)') dparm_norm_min

      ioptjas_sav=ioptjas
      ioptorb_sav=ioptorb
      ioptci_sav=ioptci
      call save_nparms

      increase_nblk=0
      energy_plus_err_best=1.d99

c CI step for state average of multiple states (optimal CI for input Jastrow and LCAO)
      if(ioptci.ne.0.and.nstates.gt.1.and.(ioptorb+ioptjas.gt.0)) then
        write(6,'(/,''Perform CI run for SA calculation'')') 
        ioptjas=0
        ioptorb=0
        add_diag_sav=add_diag(1)
        add_diag(1)=-1

        call set_nparms

        nblk_sav=nblk
        call p2gtid('optwf:nblk_ci',nblk_ci,nblk,1)
        nblk=nblk_ci
        call qmc
        nblk=nblk_sav

        call combine_derivatives

        call save_wf

        call setup_optimization(nparm,MPARMALL,MWORK,lwork,h,h_sav,s,s_sav,work,work2,add_diag(1),iter)

        write(6,'(/,''Compute CI parameters'',/)')
        call compute_dparm(nparm,MPARMALL,lwork_ci_save,grad,h,h_sav,s,s_sav,work,work2,
     &                     add_diag(1),energy(1),energy_err(1))

        call compute_parameters(grad,iflag,1)

        call write_wf(1,0)

        add_diag(1)=add_diag_sav
      endif

c Iterate optimization
      do 900 iter=1,nopt_iter

      
      write(6,'(/,''Optimization iteration'',i2)') iter
      iadd_diag_loop1=0

 100  ioptjas=ioptjas_sav
      ioptorb=ioptorb_sav
      ioptci=ioptci_sav

      if(ioptci.ne.0.and.nstates.gt.1.and.(ioptorb+ioptjas.gt.0)) ioptci=0
      call set_nparms

      nforce=1
      nwftype=1
c Generate gradient, hessian
 200  call qmc
      call combine_derivatives

      if(iter.ge.2) then
       denergy=energy(1)-energy_sav
       denergy_err=sqrt(energy_err(1)**2+energy_err_sav**2)

c For multiple states, this should never happen as you have checked the energy in the CI step and
c the CI step is unlikely to go wrong (unless the CI run is too short)
       if(denergy.gt.3*denergy_err) then
         iadd_diag_loop1=iadd_diag_loop1+1
         if(iadd_diag_loop1.gt.5) call fatal_error('OPTWF: energy went up a lot and iadd_diag_loop1 > 5')

         add_diag(1)=200*add_diag(1)
         write(6,'(/,''Iteration '',i4,'' sampling run to generate new parms '')') iter
         write(6,'(''old energy'',2f12.5)') energy_sav,energy_err_sav
         write(6,'(''new energy'',2f12.5)') energy(1),energy_err(1)
         write(6,'(/,''Energy is worse, increase adiag to'',1pd11.4)') add_diag(1)
         call restore_wf(1)
         call compute_dparm(nparm,MPARMALL,lwork_all_save,grad,h,h_sav,s,s_sav,work,work2,
     &                     add_diag(1),energy_sav,energy_err_sav)
         call compute_parameters(grad,iflag,1)
c In case starting config is very bad, reset configuration by calling sites
         isite=1
         call reset_configs_start
         if(ioptci_sav.ne.0.and.nstates.gt.1.and.(ioptorb.ne.0.or.ioptjas.ne.0)) then
c This case should never happen
           call fatal_error('OPTWF: Multiple state - Energy already checked: CI run too short?')
          else
           goto 200
         endif
       endif
      endif

c Save current energy and sigma
      energy_sav=energy(1)
      energy_err_sav=energy_err(1)

      write(6,'(/,''Current energy = '',f12.7,'' +- '',f11.7)') energy_sav,energy_err_sav
      energy_plus_err=energy(1)+2*energy_err(1)
      if(energy_plus_err.lt.energy_plus_err_best) then
        write(6,'(/,''Current best energy + 2*error = '',f11.4)') energy_plus_err
        energy_plus_err_best=energy_plus_err
        call save_wf_best(ioptjas_sav,ioptorb_sav,ioptci_sav)
      endif

      call save_wf

      call setup_optimization(nparm,MPARMALL,MWORK,lwork,h,h_sav,s,s_sav,work,work2,add_diag(1),iter)
      if(iter.eq.1) lwork_all_save=lwork

c Compute corrections to parameters
    6 write(6,'(/,''Compute parameters 1'',/)')
      call compute_dparm(nparm,MPARMALL,lwork_all_save,grad,h,h_sav,s,s_sav,work,work2,
     &                     add_diag(1),energy_sav,energy_err_sav)
 
      call test_solution_parm(nparm,grad,dparm_norm,dparm_norm_min,add_diag(1),iflag)
      if(iflag.ne.0) then
       write(6,'(''Warning: add_diag_1 has dparm_norm>1'')')
       add_diag(1)=10*add_diag(1)
       write(6,'(''adiag_1 increased to '',g12.5)') add_diag(1)
       go to 6
      endif

c     write(6,'(/,''change in parameters 1'')')
c     write(6,'(''-x='',9f15.9)') (-grad(i),i=1,nparm)

c Compute new parameters
      call compute_parameters(grad,iflag,1)
      if(iflag.ne.0) then
        write(6,'(''Warning: add_diag_1 has problems with a2 and/or b2'')')
        call restore_wf(1)
        add_diag(1)=10*add_diag(1)
        write(6,'(''adiag_1 increased to '',g12.5)') add_diag(1)
        go to 6
      endif

      if(multiple_adiag.ne.0) then

       nforce=3
       nwftype=3
       call setup_wf

       do 300 iadiag=2,3

c add_diag=add_diag*10
        add_diag(iadiag)=10**(iadiag-1)*add_diag(1)

        call restore_wf(iadiag)
        write(6,'(/,''Compute parameters '',i1,/)') iadiag
   10   call compute_dparm(nparm,MPARMALL,lwork_all_save,grad,h,h_sav,s,s_sav,work,work2,
     &                     add_diag(iadiag),energy_sav,energy_err_sav)

        call test_solution_parm(nparm,grad,dparm_norm,dparm_norm_min,add_diag(iadiag),iflag)
        if(iflag.ne.0) then
          write(6,'(''Warning: adiag_'',i1,'' has dparm_norm>1'')') iadiag
          add_diag(1)=2*10**(iadiag-1)*add_diag(1)
          write(6,'(''adiag_1 increased to '',g12.5)') add_diag(1)
          go to 6
        endif
c       write(6,'(/,''change in parameters '',i1)') iadiag
c       write(6,'(''-x='',9f15.9)') (-grad(i),i=1,nparm)
        call compute_parameters(grad,iflag,iadiag)
        if(iflag.ne.0) call fatal_error('OPTWF: adiag_1 or 2 still has problems')
 300   enddo

       write(6,'(/,''adiag1,adiag2,adiag3'',1p3g15.8,/)') (add_diag(i),i=1,3)
       write(6,'(/,''Correlated sampling test run for adiag'',/)')
 
       ioptjas=0
       ioptorb=0
       ioptci=0

c Test run for adiag_1,2,3 with correlated sampling
       nblk_sav=nblk
       nblk=max(2,nblk/2)
c      nblk=max(2,nblk/10)

       call qmc

       ioptjas=ioptjas_sav
       ioptorb=ioptorb_sav
       ioptci=ioptci_sav
       if(ioptci.ne.0.and.nstates.gt.1.and.(ioptorb.ne.0.or.ioptjas.ne.0)) ioptci=0

       nblk=nblk_sav

c Check if something is very wrong in correlated sampling run
       denergy_min=1.d+99
       denergy_max=0
       do k=1,3
         if(energy_err(k).gt.denergy_max) then
           denergy_max=energy_err(k)
           k_demax=k
         endif
         if(energy_err(k).lt.denergy_min) then
           denergy_min=energy_err(k)
           k_demin=k
         endif
       enddo
       if(denergy_max/denergy_min.gt.10) then
         write(6,'(/,''Problem with correlated sampling run'')')
         write(6,'(''e,demin,e,demax'',2(f12.5,'' +- '',f12.5))') 
     &   energy(k_demin),denergy_min,energy(k_demax),denergy_max

         if(k_demax.eq.1) then
           add_diag(1)=add_diag(1)*20
          else
           add_diag(1)=add_diag(1)*200
         endif
         write(6,'(''adiag_1 increased to '',g12.5)') add_diag(1)
         write(6,'(''generate again parameters for correlated sampling'')')

         call restore_wf(1)
c In case starting config is very bad, reset configuration by calling sites
         isite=1
         call reset_configs_start

         energy(1)=energy_sav
         go to 6
       endif

       de_worse1=energy(1)-energy_sav
       de_worse2=energy(2)-energy_sav
       de_worse3=energy(3)-energy_sav
       de_worse_err1=sqrt(energy_err(1)**2+energy_err_sav**2)
       de_worse_err2=sqrt(energy_err(2)**2+energy_err_sav**2)
       de_worse_err3=sqrt(energy_err(3)**2+energy_err_sav**2)
       write(6,'(/,''adiag, correlated energies, forces'')')
       do k=1,3
         if(k.eq.1) then
           write(6,'(g15.5,f12.5,'' +- '',f12.5)') add_diag(k),energy(k),energy_err(k)
          else
           write(6,'(g15.5,f12.5,'' +- '',f12.5,e12.5,'' +- '',e12.5)') add_diag(k),energy(k),energy_err(k),force(k),force_err(k)
         endif
       enddo
       write(6,'(''old energy'',2f12.5)') energy_sav,energy_err_sav
       write(6,'(''energy_adiag1-energy_old'',3f12.5)') de_worse1,de_worse_err1,3*de_worse_err1
       write(6,'(''energy_adiag2-energy_old'',3f12.5)') de_worse2,de_worse_err2,3*de_worse_err2
       write(6,'(''energy_adiag3-energy_old'',3f12.5)') de_worse3,de_worse_err3,3*de_worse_err3

c      if(de_worse2.gt.3*de_worse_err2.or.de_worse3.gt.3*de_worse_err3) then
       if(de_worse3.gt.3*de_worse_err3) then
c        write(6,'(/,''energy_adiag2_3 is much worse than old energy'')')
         write(6,'(/,''energy_adiag3 is much worse than old energy'')')
 
         add_diag(1)=add_diag(1)*200
         write(6,'(''adiag_1 increased to '',g12.5)') add_diag(1)
         write(6,'(''generate again parameters for correlated sampling'')')

         call restore_wf(1)
c In case starting config is very bad, reset configuration by calling sites
         isite=1
         call reset_configs_start
 
         energy(1)=energy_sav
         go to 6
       endif

       write(6,'(/,''find optimal adiag'')')
c Find optimal a_diag
       call quad_min

       call restore_wf(1)

   7   call compute_dparm(nparm,MPARMALL,lwork_all_save,grad,h,h_sav,s,s_sav,work,work2,
     &                     add_diag(1),energy_sav,energy_err_sav)

       call test_solution_parm(nparm,grad,dparm_norm,dparm_norm_min,add_diag(1),iflag)
       if(iflag.ne.0) then
        write(6,'(''Warning: adiag_optimal has dparm_norm>1'')')
        add_diag(1)=200*add_diag(1)
        write(6,'(''adiag_1 increased to '',g12.5)') add_diag(1)
        call restore_wf(1)
        go to 7
       endif

c      write(6,'(/,''Optimal change in parameters'')')
c      write(6,'(''-x='',9f15.9)') (-grad(i),i=1,nparm)

c Compute new parameters
       write(6,'(/,''Compute parameters for optimal adiag'')')
       call compute_parameters(grad,iflag,1)
       if(iflag.ne.0) then
         write(6,'(''Warning: adiag_optimal has problems'')')
         add_diag(1)=200*add_diag(1)
         write(6,'(''adiag_1 increased to '',g12.5)') add_diag(1)
         call restore_wf(1)
         go to 7
       endif

       call write_wf(1,iter)

       call check_length_run(iter,increase_nblk,nblk,nblk_max,denergy,denergy_err,energy_err_sav,energy_tol)

       add_diag(1)=0.1d0*add_diag(1)
      else
       call write_wf(1,iter)
c endif for multiple_adiag
      endif

      ioptjas=ioptjas_sav
      ioptorb=ioptorb_sav
      ioptci=ioptci_sav

c CI step for state average of multiple states
      if(ioptci.ne.0.and.nstates.gt.1.and.(ioptorb+ioptjas.gt.0)) then
  800   ioptjas=0
        ioptorb=0
        ioptci=ioptci_sav

        nforce=1
        nwftype=1

        add_diag_sav=add_diag(1)
        add_diag(1)=-1

        call set_nparms

        nblk_sav=nblk
        call p2gtid('optwf:nblk_ci',nblk_ci,nblk,1)
        nblk=nblk_ci
        call qmc
        nblk=nblk_sav
        
        call combine_derivatives

        if(iter.ge.1) then
          denergy=energy(1)-energy_sav
          denergy_err=sqrt(energy_err(1)**2+energy_err_sav**2)

          if(denergy.gt.3*denergy_err) then
            iadd_diag_loop1=iadd_diag_loop1+1
            if(iadd_diag_loop1.gt.5) call fatal_error('OPTWF: energy went up a lot and iadd_diag_loop1 > 5')

c add_diag used to generate current Jastrow+orbitals was divided by 10 before the end of the loop 900
c           add_diag(1)=20*add_diag_sav
            add_diag(1)=200*add_diag_sav
            write(6,'(/,''Iteration '',i4,'' sampling run to generate new CI coefs'')') iter
            write(6,'(''old energy'',2f12.5)') energy_sav,energy_err_sav
            write(6,'(''new energy'',2f12.5)') energy(1),energy_err(1)
            write(6,'(/,''Energy is worse, increase adiag to '',1pd11.4)') add_diag(1)

c Jastrow and orbital parameters give worse energy 
            ioptjas=ioptjas_sav
            ioptorb=ioptorb_sav
            ioptci=0

            call set_nparms

            call restore_wf(1)
            call compute_dparm(nparm,MPARMALL,lwork_all_save,grad,h,h_sav,s,s_sav,work,work2,
     &                     add_diag(1),energy_sav,energy_err_sav)
            call compute_parameters(grad,iflag,1)
c In case starting config is very bad, reset configuration by calling sites
            isite=1
            call reset_configs_start

            call write_wf(1,iter)
            goto 800
          endif
        endif

        call setup_optimization(nparm,MPARMALL,MWORK,lwork,h,h_sav,s,s_sav,work,work2,add_diag(1),iter)
        if(iter.eq.1) lwork_ci_save=lwork

        write(6,'(/,''Compute CI parameters'',/)')
        call compute_dparm(nparm,MPARMALL,lwork_ci_save,grad,h,h_sav,s,s_sav,work,work2,
     &                     add_diag(1),energy(1),energy_err(1))

        call compute_parameters(grad,iflag,1)

c save CI coefficients
        call save_wf
        call write_wf(1,iter)

c save orb and jastrow 
        ioptjas=ioptjas_sav
        ioptorb=ioptorb_sav
        ioptci=0
        call save_wf

        add_diag(1)=add_diag_sav

        call set_nparms
c endif CI step for multiple states
      endif

c end of optimization loop
 900  continue

 950  nforce=1

      ioptjas=0
      ioptorb=0
      ioptci=0

      call p2gtid('optwf:ilastvmc',ilastvmc,1,1)
      if(ilastvmc.eq.0) go to 970

      call qmc
      write(6,'(/,''Current energy = '',f12.7,'' +- '',f11.7)') energy(1),energy_err(1)
      energy_plus_err=energy(1)+2*energy_err(1)
      if(energy_plus_err.lt.energy_plus_err_best) then
        write(6,'(/,''Current best energy + 2*error = '',f11.4)') energy_plus_err
        energy_plus_err_best=energy_plus_err
        call save_wf_best(ioptjas_sav,ioptorb_sav,ioptci_sav)
      endif

 970  ioptjas=ioptjas_sav
      ioptorb=ioptorb_sav
      ioptci=ioptci_sav

      call write_wf_best

      return
      end
c-----------------------------------------------------------------------
      subroutine check_length_run(iter,increase_nblk,nblk,nblk_max,denergy,denergy_err,energy_err_sav,energy_tol)

      implicit real*8(a-h,o-z)

c Increase nblk if near convergence to value needed to get desired statistical error
      increase_nblk=increase_nblk+1

c Check if energies for 3 different values of a_diag are less than the tolerance
c     energy_min= 1.d99
c     energy_max=-1.d99
c     do 5 k=1,3
c       energy_min=min(energy_min,energy(k))
c   5   energy_max=max(energy_max,energy(k))
c     e_diff=energy_max-energy_min
c     write(6,'(''iter,e_diff='',i4,d12.4)') iter,e_diff
c
c     dforce2=force_err(2)-energy_tol
c     dforce3=force_err(2)-energy_tol
c     if(e_diff.lt.energy_tol.and.dforce2.lt.0.and.dforce3.lt.0) then
c       nblk_new=nblk*max(1.d0,(energy_err_sav/energy_tol)**2)
c       nblk_new=min(nblk_new,nblk_max)
c       if(nblk_new.gt.nblk) then
c         increase_nblk=0
c         nblk=nblk_new
c         write(6,'(''nblk reset to'',i8,9d12.4)') nblk,energy_err(1),energy_tol
c       endif
c       write(6,'(''energy differences for different add_diag converged to'',d12.4)') energy_tol
c       goto 950
c     endif

c Increase if subsequent energies are within errorbar
      if(iter.gt.2.and.dabs(denergy).lt.3*denergy_err) then
        nblk_new=nblk*max(1.d0,(energy_err_sav/energy_tol)**2)
        nblk_new=min(nblk_new,nblk_max)
        if(nblk_new.gt.nblk) then
          increase_nblk=0
          nblk=nblk_new
          write(6,'(''nblk reset to'',i8,9d12.4)') nblk,dabs(denergy),energy_tol
        endif
      endif
c Always increase nblk by a factor of 2 every other iteration
      if(increase_nblk.eq.2.and.nblk.lt.nblk_max) then
        increase_nblk=0
        nblk=min(2*nblk,nblk_max)
        write(6,'(''nblk reset to'',i8,9d12.4)') nblk
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine quad_min

      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      use optwf_corsam, only: add_diag, add_diag_tmp, energy, energy_err, force, force_err
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'

      parameter(MFUNC=3)


      dimension add_diag_log(MFUNC),a(MFUNC,MFUNC),b(MFUNC)

      npts=3
      nfunc=3

      do 5 k=1,npts
    5   add_diag_log(k)=dlog10(add_diag(k))

      do 30 i=1,nfunc
        b(i)=0
        do 10 k=1,npts
   10     b(i)=b(i)+energy(k)*add_diag_log(k)**(i-1)
        do 30 j=1,i
          a(i,j)=0
          do 20 k=1,npts
   20       a(i,j)=a(i,j)+add_diag_log(k)**(i+j-2)
   30     a(j,i)=a(i,j)

c Do cholesky decomposition
      call chlsky(a,nfunc,MFUNC,ierr)
      if(ierr.ne.0) stop 'ierr ne 0 in chlsky'

c Symmetrize decomposed matrix (needs to be done before calling uxb
c or need to modify uxb)
      do 40 i=1,nfunc
        do 40 j=i+1,nfunc
   40     a(i,j)=a(j,i)

c Solve linear equations
      call lxb(a,nfunc,MFUNC,b)
      call uxb(a,nfunc,MFUNC,b)

      write(6,'(''polinomial coeffcients b1+b2*adiag+b3*adiag^2'',f12.5,1p2e12.4)') (b(i),i=1,nfunc)
      energy_min= 1.d99
      energy_max=-1.d99
      rms=0
      do 50 k=1,npts
        ee=b(1)+b(2)*add_diag_log(k)+b(3)*add_diag_log(k)**2
        write(6,'(''fit log(adiag),e_fit,e '',3f12.5)') add_diag_log(k),ee,energy(k)
        if(energy(k).lt.energy_min) then
          k_min=k
          energy_min=energy(k)
        endif
        energy_max=max(energy_max,energy(k))
   50   rms=rms+(ee-energy(k))**2
      rms=dsqrt(rms/npts)
      write(6,'(''rms error in fit of energy to get optimal add_diag is'',d12.4)') rms

      energy_var=energy_max-energy_min
      if(b(3).gt.0.and.abs(force(2)).gt.3*force_err(2).and.abs(force(3)).gt.3*force_err(3)) then
        iwadd_diag=0
        add_diag_log_min=-0.5d0*b(2)/b(3)
        add_diag_log_min=min(max(add_diag_log_min,add_diag_log(1)-1),add_diag_log(1)+3)
        write(6,'(/,''computed optimal adiag '',g12.4)') 10**add_diag_log_min
        eopt=b(1)+b(2)*add_diag_log_min+b(3)*add_diag_log_min**2
        write(6,'(/,''computed optimal energy'',f12.5)') eopt
       elseif(energy(1).lt.energy(2)+force_err(2).and.energy_var.lt.1.d-3*abs(energy_max)) then
        iwadd_diag=1
        add_diag_log_min=add_diag_log(1)
        if(energy(1).lt.energy(2).and.add_diag_log_min.ge.0.d0) add_diag_log_min=add_diag_log_min-1.d0
       elseif(energy(2).lt.energy(3)+force_err(3).and.energy_var.lt.1.d-2*abs(energy_max).and.k_min.eq.3) then
        iwadd_diag=2
        add_diag_log_min=add_diag_log(2)
       else
        iwadd_diag=k_min
        write(6,'(/,''b3 < 0 or error on one force too large'')')
        add_diag_log_min=add_diag_log(k_min)
        if(k_min.eq.1.and.energy(1).lt.energy(2).and.add_diag_log_min.ge.0.d0) add_diag_log_min=add_diag_log_min-1.d0
      endif
      add_diag_log_min=max(add_diag_log_min,-6*1.d0)
      add_diag_min=10**add_diag_log_min
      write(6,'(/,''optimal adiag '',i2,g12.4,/)') iwadd_diag,add_diag_min

      add_diag(1)=add_diag_min

      return
      end
c-----------------------------------------------------------------------
      subroutine combine_derivatives

      use gradhess_ci, only: grad_ci, h_ci, s_ci
      use gradhess_jas, only: grad_jas, h_jas, s_jas
      use gradhess_mix_jas_ci, only: h_mix_jas_ci, s_mix_jas_ci
      use gradhess_mix_jas_orb, only: h_mix_jas_orb, s_mix_jas_orb
      use gradhess_mix_orb_ci, only: h_mix_ci_orb, s_mix_ci_orb
      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      use optwf_parms, only: nparmd, nparme, nparmg, nparmj, nparml, nparms
      implicit real*8(a-h,o-z)







      include 'vmc.h'
      include 'force.h'
      include 'optjas.h'
      include 'optci.h'
      include 'optorb.h'
        
      parameter(MPARMALL=MPARMJ+MXCIREDUCED+MXREDUCED)

c     common /gradhess_orb/ grad_orb(MXORBOP),h_orb(MXMATDIM),s_orb(MXMATDIM)


      common /gradhess_all/ grad(MPARMALL),h(MPARMALL,MPARMALL),s(MPARMALL,MPARMALL)


c Note: we do not vary the first (i0) CI coefficient unless full CI

      if(method.eq.'linear') then
        
       is=1
       i0=0
       ishift=1
       if(ioptjas.eq.0) go to 115 

c Jastrow Hamiltonian
       do 110 j=1,nparmj+is
         do 110 i=1,nparmj+is
           h(i,j)=h_jas(i,j)
  110      s(i,j)=s_jas(i,j)
        
c      do 111 i=1,nparmj+1
c 111    write(6,'(''h1= '',1000d12.5)') (h(i,j),j=1,nparmj+1)
c      do 112 i=1,nparmj+1
c 112    write(6,'(''h1= '',1000d12.5)') (s(i,j),j=1,nparmj+1)

       ishift=nparmj+is

  115  continue

       if(ioptci.eq.0) go to 135

       if(ioptjas.eq.0.and.ioptorb.eq.0) then
        is=0
        i0=0
        ishift=0
       else
        i0=1
        h(1,1)=h_ci(1,1)
        s(1,1)=s_ci(1,1)
        do 125 i=1,nciterm-i0
          h(ishift+i,1)=h_ci(i+i0+is,1)
          h(1,ishift+i)=h_ci(1,i+i0+is)
          s(ishift+i,1)=s_ci(i+i0+is,1)
  125     s(1,ishift+i)=s_ci(1,i+i0+is)
       endif

c CI Hamiltonian
       do 120 j=1,nciterm-i0
         do 120 i=1,nciterm-i0
           h(ishift+i,ishift+j)=h_ci(i+i0+is,j+i0+is)
  120      s(ishift+i,ishift+j)=s_ci(i+i0+is,j+i0+is)

c      write(6,'(''h2 shift ='',i4)') ishift
c      do 121 i=1,nciterm-i0
c 121    write(6,'(''h2= '',1000f12.5)') (h_ci(i+i0+is,j+i0+is),j=1,nciterm-i0)

c Jastrow-CI Hamiltonian
       do 130 j=1,nciterm-i0
         do 130 i=1,nparmj
           h(i+1,j+ishift)=h_mix_jas_ci(i,j+i0)
           h(j+ishift,i+1)=h_mix_jas_ci(i+nparmj,j+i0)
           s(i+1,j+ishift)=s_mix_jas_ci(i,j+i0)
  130      s(j+ishift,i+1)=s_mix_jas_ci(i,j+i0)

c      do 131 i=1,nparmj
c        write(6,'(''h3= '',1000f12.5)') (h_mix_jas_ci(i,j+i0),j=1,nciterm-i0)
c 131    write(6,'(''h3= '',1000f12.5)') (h_mix_jas_ci(i+nparmj,j+i0),j=1,nciterm-i0)

       ishift=ishift+nciterm-i0

  135  continue

       if(ioptorb.eq.0) go to 175

c      h(1,1)=h_orb(1)
c      s(1,1)=s_orb(1)
c      do 140 i=1,nreduced
c        ik=i*(nreduced+1)
c        h(ishift+i,1)=h_orb(i+1)
c        h(1,ishift+i)=h_orb(ik+1)
c        s(ishift+i,1)=s_orb(i+1)
c 140    s(1,ishift+i)=s_orb(ik+1)

c ORB Hamiltonian
c     do 150 j=1,nreduced
c       jk=j*(nreduced+1)
c        do 150 i=1,nreduced
c          h(ishift+i,ishift+j)=h_orb(i+1+jk)
c 150      s(ishift+i,ishift+j)=s_orb(i+1+jk)

c Jastrow-ORB Hamiltonian
       do 160 j=1,nreduced
         do 160 i=1,nparmj
           h(i+1,j+ishift)=h_mix_jas_orb(i,j)
           h(j+ishift,i+1)=h_mix_jas_orb(i+nparmj,j)
           s(i+1,j+ishift)=s_mix_jas_orb(i,j)
  160      s(j+ishift,i+1)=s_mix_jas_orb(i,j)

c ORB-CI Hamiltonian
       do 170 j=1,nreduced
         do 170 i=1,nciterm-i0
           h(i+nparmj+1,j+ishift)=h_mix_ci_orb(i+i0,j)
           h(j+ishift,i+nparmj+1)=h_mix_ci_orb(i+nciterm+i0,j)
           s(i+nparmj+1,j+ishift)=s_mix_ci_orb(i+i0,j)
  170      s(j+ishift,i+nparmj+1)=s_mix_ci_orb(i+i0,j)

  175  nparm=nparmj+nciterm+nreduced-i0

       write(6,'(/,''number of parms: total, Jastrow, CI, orbitals= '',4i5)') 
     & nparm,nparmj,nciterm,nreduced

c      do 180 i=1,nparm+1
c 180    write(6,'(''h= '',1000d12.5)') (h(i,j),j=1,nparm+1)
c      do 185 i=1,nparm+1
c 185    write(6,'(''s= '',1000d12.5)') (s(i,j),j=1,nparm+1)

      endif

      return
      end
c-----------------------------------------------------------------------
      block data optprt_count

      implicit real*8(a-h,o-z)

      common /icount_ci/ icount_ci
      common /icount_orb/ icount_orb
      common /icount_prop/ icount_prop
      data icount_ci /1/
      data icount_orb /1/
      data icount_prop /1/

      end
