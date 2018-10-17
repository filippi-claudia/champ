      subroutine optwf
c written by Claudia Filippi

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

      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar

      common /numbas/ arg(MCTYPE),r0(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS,MCTYPE)

      common /wfsec/ iwftype(MFORCE),iwf,nwftype
      common /forcepar/ deltot(MFORCE),nforce,istrech

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_corsam/ add_diag(MFORCE),energy(MFORCE),energy_err(MFORCE),force(MFORCE),force_err(MFORCE)

      common /gradhess_all/ grad(MPARMALL),h(MPARMALL,MPARMALL),s(MPARMALL,MPARMALL)

      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj

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

        call setup_optimization(nparm,MPARMALL,MWORK,lwork,grad,grad_sav,h,h_sav,s,s_sav,work,work2,add_diag(1),iter)

        write(6,'(/,''Compute CI parameters'',/)')
        call compute_dparm(nparm,MPARMALL,lwork_ci_save,grad,grad_sav,h,h_sav,s,s_sav,work,work2,
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
         call compute_dparm(nparm,MPARMALL,lwork_all_save,grad,grad_sav,h,h_sav,s,s_sav,work,work2,
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

      call setup_optimization(nparm,MPARMALL,MWORK,lwork,grad,grad_sav,h,h_sav,s,s_sav,work,work2,add_diag(1),iter)
      if(iter.eq.1) lwork_all_save=lwork

c Compute corrections to parameters
    6 write(6,'(/,''Compute parameters 1'',/)')
      call compute_dparm(nparm,MPARMALL,lwork_all_save,grad,grad_sav,h,h_sav,s,s_sav,work,work2,
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
   10   call compute_dparm(nparm,MPARMALL,lwork_all_save,grad,grad_sav,h,h_sav,s,s_sav,work,work2,
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

   7   call compute_dparm(nparm,MPARMALL,lwork_all_save,grad,grad_sav,h,h_sav,s,s_sav,work,work2,
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
            call compute_dparm(nparm,MPARMALL,lwork_all_save,grad,grad_sav,h,h_sav,s,s_sav,work,work2,
     &                     add_diag(1),energy_sav,energy_err_sav)
            call compute_parameters(grad,iflag,1)
c In case starting config is very bad, reset configuration by calling sites
            isite=1
            call reset_configs_start

            call write_wf(1,iter)
            goto 800
          endif
        endif

        call setup_optimization(nparm,MPARMALL,MWORK,lwork,grad,grad_sav,h,h_sav,s,s_sav,work,work2,add_diag(1),iter)
        if(iter.eq.1) lwork_ci_save=lwork

        write(6,'(/,''Compute CI parameters'',/)')
        call compute_dparm(nparm,MPARMALL,lwork_ci_save,grad,grad_sav,h,h_sav,s,s_sav,work,work2,
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
      subroutine setup_optimization(nparm,mparmx,MWORK,lwork,grad,grad_sav,h,h_sav,s,s_sav,work,eig_vec,add_diag,iter)

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optjas.h'
      include 'optci.h'
      include 'optorb.h'
      include 'optorb_cblk.h'

      parameter(MPARMALL=MPARMJ+MXCIREDUCED+MXREDUCED,eps=1.d-12)

      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj
      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm_save
      common /optwf_corsam/ add_diag_tmp(MFORCE),energy(MFORCE),energy_err(MFORCE),force(MFORCE),force_err(MFORCE)

      common /linear_norm/ oav(MXCITERM)

      dimension h(mparmx,*),grad(*),s(mparmx,*)
      dimension h_sav(mparmx,*),grad_sav(*),s_sav(*)

      dimension eig(MPARMALL),eigi(MPARMALL),seig_inv(MPARMALL)
      dimension eig_vec(MPARMALL,*),hmod(MPARMALL,MPARMALL)
      dimension isort(MPARMALL)

      dimension work(*)

      save eig_min

      write(6,'(/,''Setup starting adiag'',/)') 
      eig_min=0
      if(method.eq.'hessian') then

c Symmetrize the hessian, save the hessian and the gradient
      do 1 i=1,nparm
        h_sav(i,i)=h(i,i)
        grad_sav(i)=grad(i)
        do 1 j=1,i-1
          h_sav(i,j)=h(i,j)
          h_sav(j,i)=h(i,j)
    1     h(j,i)=h(i,j)

      write(6,'(''grad='',9g12.4)') (grad(i),i=1,nparm)

      grad_norm=0.d0
      do 2 i=1,nparm
    2   grad_norm=grad_norm+grad(i)*grad(i)
      grad_norm=dsqrt(grad_norm)
      write(6,'(''grad_norm='',g10.5,/)') grad_norm

c Make hessian positive definite: add min(0,-lowest eigenvalue) to diagonal
      call posdef_hessian(nparm,mparmx,h,h_sav,work,MWORK,eig_min)

      elseif(method.eq.'linear') then

      nparmd=max(nciterm-1,0)

c If we optimize jas, orb, or both, we solve secular super-CI equation on basis of psi_0,dpsi_0 for a total dimension of nparm+1
      is=1
      if(ioptjas.eq.0.and.ioptorb.eq.0) then
c If CI only, we solve secular equation on basis of J*CSF
        is=0
       else
c Save the overlap and hamiltonian
        idx=0
        do 3 i=1,nparm+is
          do 3 j=1,i
            idx=idx+1
            s_sav(idx)=s(i,j)
            h_sav(i,j)=h(i,j)
   3        h_sav(j,i)=h(j,i)
      endif

c Symmetrize the overlap
      do 4 i=1,nparm+is
        do 4 j=1,i
   4      s(j,i)=s(i,j)

c     do i=1,nparm+is
c       write(6,*) 'h =',(h(i,j),j=1,nparm+is)
c     enddo
c     do i=1,nparm+is
c       write(6,*) 's =',(s(i,j),j=1,nparm+is)
c     enddo

      if(add_diag.gt.0.and.iter.eq.1) then

      write(6,'(''Determine energy gap in super-CI hamiltonian at step 1'',/)')
      call regularize_geneig(nparm+is,mparmx,h,s,work,seig_inv,hmod)
      call solve_geneig(nparm+is,mparmx,hmod,s,seig_inv,work,eig,eigi,eig_vec)

      imag=0
      do 5 j=1,nparm+is
       if (eigi(j).ne.0.d0) imag=1
    5 continue
      if(imag.eq.1) write(6,'(''Warning: imaginary eigenvalues'')')

c Sort the eigenvalues
      call sort(nparm+is,eig,isort)

      ireal=0
      do 6 ii=1,nparm+is
        i=isort(ii)
        if(eigi(i).eq.0) then 
         ireal=ireal+1
         if(ireal.le.5)
     &   write(6,'(''eigenvalue '',i4,'' = '',f15.5,'' + i* '',f15.5)') i,eig(i),eigi(i)
        endif
    6 continue
      write(6,*)

c use semiorthogonal basis and compute minimum norm in direction orthogonal to psi0
      dmult=1.d0
      scale=.1d0
    9 de_range=scale*dabs(energy(1))
      emin=energy(1)-de_range
      emax=energy(1)+dmult*3*energy_err(1)

      ireal=0
      anorm_orth_min=1.d+99
      do 20 ii=1,nparm+is
        i=isort(ii)
        if(eigi(i).ne.0.or.eig(i).lt.emin) go to 20
        if(eig(i).gt.emax) go to 25
        ireal=ireal+1
        bot=1.d0
        do 10 j=1,nparmd
   10     bot=bot-eig_vec(j+nparmj+is,i)*oav(j+1)/eig_vec(1,i)
        idx=0
        anorm_orth=0.d0
        do 15 j=1,nparm+is
          do 15 k=1,j
            idx=idx+1
c 21/8 - ERROR?
c           if(j.ne.1.and.k.ne.1) then
            if(j.ne.1.or.k.ne.1) then
              anorm_orth=anorm_orth+eig_vec(j,i)*eig_vec(k,i)*s_sav(idx)
              if(j.ne.k) anorm_orth=anorm_orth+eig_vec(k,i)*eig_vec(j,i)*s_sav(idx)
            endif
   15   continue
        anorm_orth=sqrt(dabs(anorm_orth))/abs(bot*eig_vec(1,i))
        if(anorm_orth.lt.anorm_orth_min) then
          anorm_orth_min=anorm_orth
          i_ovr=i
          isort_ovr=ii
        endif
        if(ireal.le.5) then
          write(6,'(i4,'' eigenvalue '',i4,'' = '',f15.5,'' + i* '',f15.5)') ii,i,eig(i),eigi(i)
          write(6,'(4x,'' ortho dir  '',i4,'' = '',f15.5)') i,anorm_orth
        endif
   20 continue
   25 write(6,'(i4,'' real roots in energy range '',f15.5,'','',f15.5)') ireal,emin,emax
      if(ireal.eq.0) then
        dmult=dmult*2
        scale=scale*2
        goto 9
      endif
      write(6,'(''maximum overlap = '',i5,f15.5,'' + i * '',2f15.5)') i_ovr,eig(i_ovr),eigi(i_ovr),anorm_orth_min

      i0=i_ovr
      if(isort_ovr.lt.nparm+is) then
        i1=isort(isort_ovr+1)
        eig_min=eig(i1)-eig(i0)
      else
        eig_min=0
      endif
 
      write(6,'(''energy gap = '',f15.5)') eig_min
      endif

      endif

      call p2gtid('optwf:multiple_adiag',multiple_adiag,0,1)
c Set add_diag
      if(add_diag.gt.0.or.multiple_adiag.eq.1) then
        add_diag=max(add_diag,1.d-6)
        add_diag=max(add_diag,1.d-2*eig_min)
c       add_diag=max(add_diag,1.d-2*abs(eig_min))
       else
c       add_diag=0
        add_diag=dabs(add_diag)
      endif
      if(ioptjas.eq.0.and.ioptorb.eq.0) add_diag=0
      write(6,'(/,''starting adiag'',g12.4,/)') add_diag
       

      return
      end
c-----------------------------------------------------------------------
      subroutine posdef_hessian(nparm,mparmx,h,h_sav,work,MWORK,eig_min)

      implicit real*8(a-h,o-z)

      include 'optjas.h'

c TEMPORARY: complaints on novamaris with eigv(mparmx)
      dimension h(mparmx,*),h_sav(mparmx,*)
      dimension work(*),eigv(MPARMJ)

c Calculate eigenvalues
      call dsyev('V','U',nparm,h,mparmx,eigv,work,-1,info)
      if(work(1).gt.MWORK) then
        write(6,'(''work(1), MWORK'',2i7)') int(work(1)),MWORK
        call fatal_error('POSDEF_HESS: work(1).gt.MWORK')
      endif
      if(info.ne.0) call fatal_error('POSDEF_HESS: info from dggev != 0')
      lwork=int(work(1))
      call dsyev('V','U',nparm,h,mparmx,eigv,work,lwork,info)
      write(6,'(''eigs='',1p9g8.1)') (eigv(i),i=1,nparm)

      eig_min=1.d99
      do 1 i=1,nparm
    1   eig_min=min(eig_min,eigv(i))
      write(6,'(''eig_min='',(1p20g8.1))') eig_min

c Add min(0,-lowest eigenvalue) to diagonal of hessian
c Another possibility is first diagonalizing hess by a unitary transform
c taking the inverse of the diagonalized matrix to construct the
c inverse of A, but resetting the diagonal elements of the diagonalized A
c to a minimum positive value if they are less than this value.
      do 5 i=1,nparm
        h_sav(i,i)=h_sav(i,i)+max(-eig_min,0.d0)
    5 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_dparm(nparm,mparmx,lwork,grad,grad_sav,h,h_sav,s,s_sav,work,eig_vec,
     &                     add_diag,energy_sav,energy_err_sav)

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optjas.h'
      include 'optci.h'
      include 'optorb.h'
      include 'optorb_cblk.h'

      parameter(MPARMALL=MPARMJ+MXCIREDUCED+MXREDUCED,eps=1.d-12)

      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj
      common /optwf_wjas/ iwjasa(83,MCTYP3X),iwjasb(83,3),iwjasc(83,MCTYPE),iwjasf(15,MCTYPE)
      common /optwf_nparmj/ nparma(MCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE)

      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparmsav

      common /linear_norm/ oav(MXCITERM)

      dimension grad(*),grad_sav(*)
      dimension h(mparmx,*),h_sav(mparmx,*)
      dimension s(mparmx,*),s_sav(*)
      dimension work(*)

      dimension eig(MPARMALL),eigi(MPARMALL),seig_inv(MPARMALL)
      dimension eig_vec(MPARMALL,*),hmod(MPARMALL,MPARMALL)
      dimension s_norm(MXCIMATDIM)
      dimension isort(MPARMALL)
      dimension cdelta(MPARMALL)
      dimension overlap(MXCITERM)
    
      if(method.eq.'hessian') then

      do 8 i=1,nparm
        grad(i)=grad_sav(i)
        h(i,i)=h_sav(i,i)+add_diag
        do 8 j=1,i-1
          h(i,j)=h_sav(i,j)
          h(j,i)=h_sav(j,i)
    8 continue

      call p2gtid('optwf:cholesky',ichlsky,1,1)

      if(ichlsky.gt.0) then
c Do cholesky decomposition
       call chlsky(h,nparm,mparmx,ierr)
       if(ierr.ne.0) goto 99

c Symmetrize decomposed matrix (needs to be done before calling uxb
c or need to modify uxb)
       do 20 i=1,nparm
         do 20 j=i+1,nparm
   20      h(i,j)=h(j,i)

c Solve linear equations
       call lxb(h,nparm,mparmx,grad)
       call uxb(h,nparm,mparmx,grad)

       return
   99  stop 'ierr.ne.0 in chlsky'
      else
       call p2gtfd('optwf:svd_threshold',ftol,1.0d-4,1)
       call svdchk(h,mparmx,nparm,grad,cdelta,nkept,ftol,1)
       do 100 i=1,nparm
  100    grad(i)=cdelta(i)
      endif

      elseif(method.eq.'linear') then

      nparmd=max(nciterm-1,0)

      is=1
      if(ioptjas.eq.0.and.ioptorb.eq.0) then 
	is=0
        idx=0
        do 105 i=1,nparm+is
          do 105 j=1,i
            idx=idx+1
  105       s_norm(idx)=s(i,j)
       else
        idx=0
        do 107 i=1,nparm+is
          do 107 j=1,i
            idx=idx+1
            s(i,j)=s_sav(idx)
  107       s(j,i)=s_sav(idx)
        do 108 i=1,nparm+is
          do 108 j=1,nparm+is
            h(i,j)=h_sav(i,j)
  108   continue
        do 109 i=2,nparm+is
          h(i,i)=h(i,i)+add_diag
  109   continue
      endif

c     do i=1,nparm+is
c     write(6,'(''h '',1000e25.15)') (h(i,j),j=1,nparm+is)
c     enddo
c     do i=1,nparm+is
c     write(6,'(''s '',1000e25.15)') (s(i,j),j=1,nparm+is)
c     enddo

      call regularize_geneig(nparm+is,mparmx,h,s,work,seig_inv,hmod)
      call solve_geneig(nparm+is,mparmx,hmod,s,seig_inv,work,eig,eigi,eig_vec)

      imag=0
      do 110 j=1,nparm+is
       if (eigi(j).ne.0.d0) imag=1
  110 continue

      if(imag.eq.1) write(6,'(''Warning: imaginary eigenvalues'')')

c Sort the eigenvalues
      call sort(nparm+is,eig,isort)

      i_min=isort(1)
      write(6,'(''generalized eigenvalue problem H c= E S c'')')
      write(6,'(''lowest eigenval = '',i5,f15.5,'' + i * '',f15.5)') i_min,eig(i_min),eigi(i_min)

c use semiorthogonal basis and compute minimum norm in direction orthogonal to psi0
      dmult=1.d0
      scale=.1d0

      no_real_found=0
  116 de_range=scale*dabs(energy_sav)
      emin=energy_sav-de_range
      emax=energy_sav+dmult*3*energy_err_sav

      ireal=0
      anorm_orth_min=1.d+99
      do 119 ii=1,nparm+is
        i=isort(ii)
        if(eigi(i).ne.0.or.eig(i).lt.emin) go to 119
        if(eig(i).gt.emax) go to 120
        ireal=ireal+1

        if(ioptjas.ne.0.or.ioptorb.ne.0) then

          bot=1.d0
          do 117 j=1,nparmd
  117        bot=bot-eig_vec(nparmj+is+j,i)*oav(j+1)/eig_vec(1,i)
          idx=0
          anorm_orth=0.d0
          do 118 j=1,nparm+is
            do 118 k=1,j
              idx=idx+1
c 21/8 - ERROR?
c             if(j.ne.1.and.k.ne.1) then
              if(j.ne.1.or.k.ne.1) then
                anorm_orth=anorm_orth+eig_vec(j,i)*eig_vec(k,i)*s_sav(idx)
                if(j.ne.k) anorm_orth=anorm_orth+eig_vec(k,i)*eig_vec(j,i)*s_sav(idx)
              endif
  118     continue
          anorm_orth=sqrt(dabs(anorm_orth))/abs(bot*eig_vec(1,i))
          if(anorm_orth.lt.anorm_orth_min) then
            anorm_orth_min=anorm_orth
            i_ovr=i
          endif
          if(ireal.le.5) then
            write(6,'(i4,'' eigenvalue '',i4,'' = '',f15.5,'' + i* '',f15.5)') ii,i,eig(i),eigi(i)
            write(6,'(4x,'' ortho dir  '',i4,'' = '',f15.5)') i,anorm_orth
          endif
        else
          if(ireal.le.5) write(6,'(i4,'' eigenvalue '',i4,'' = '',f15.5,'' + i* '',f15.5)') ii,i,eig(i),eigi(i)
        endif
  119 continue
  120 write(6,'(i4,'' real roots in energy range '',f15.5,'','',f15.5)') ireal,emin,emax
      if(ireal.eq.0) then
        dmult=dmult*2
        scale=scale*2

        no_real_found=no_real_found+1
        if(no_real_found.gt.nparmd) 
     &      call fatal_error('OPTWF: Cannot find real eigenvalues')
        goto 116
      endif

      if(ioptjas.eq.1.or.ioptorb.eq.1) then
        write(6,'(''maximum overlap = '',i5,f15.5,'' + i * '',2f15.5)') i_ovr,eig(i_ovr),eigi(i_ovr),anorm_orth_min

        if(i_ovr.ne.i_min) write(6,'(''Warning: max overlap not for min eigenvalue'')')

        i_good=i_ovr
        do 130 i=2,nparm+is
  130    cdelta(i-1)=eig_vec(i,i_good)/eig_vec(1,i_good)

        write(6,'(''dp  ='',1000f10.5)') (cdelta(i),i=1,nparm)

        bot=1
        do 145 i=1,nparmd
  145     bot=bot-cdelta(nparmj+i)*oav(i+1)

c minus sign because variation is subtracted when computing new parameters
        do 150 i=1,nparm
  150     cdelta(i)=-cdelta(i)/bot

        write(6,'(''dpn ='',1000f10.5)') (-cdelta(i),i=1,nparm)

        do 160 i=1,nparm
  160     grad(i)=cdelta(i)

       else
        i0=0
        do 180 jj=1,nstates

          write(6,'(''State '',i4)') jj
          dnorm_jj=0
          idx=0

          if(ncsf.gt.0) then
            do 162 i=1,nparm
              do 162 k=1,i
                idx=idx+1
                dmul=1.d0
                if(i.ne.k) dmul=2.d0
  162            dnorm_jj=dnorm_jj+dmul*ccsf(i,jj,1)*ccsf(k,jj,1)*s_norm(idx)
          else
            do 163 i=1,nparm
              do 163 k=1,i
                idx=idx+1
                dmul=1.d0
                if(i.ne.k) dmul=2.d0
  163            dnorm_jj=dnorm_jj+dmul*cdet(i,jj,1)*cdet(k,jj,1)*s_norm(idx)
          endif
          dnorm_jj=1.d0/dsqrt(dnorm_jj)

          write(6,'(''state '',i4,'' norm '',1p1e12.5)') jj,dnorm_jj

          do 172 j=i0+1,nparm

            overlap(j)=0.d0

            jsort=isort(j)

            if(eig(jsort).ne.0.and.eigi(jsort).eq.0.d0) then
              idx=0
              if(ncsf.gt.0) then
                do 166 i=1,nparm
                  do 166 k=1,i
                    idx=idx+1
                    if(i.ne.k) then
                      overlap(j)=overlap(j)+(eig_vec(i,jsort)*ccsf(k,jj,1)+eig_vec(k,jsort)*ccsf(i,jj,1))*s_norm(idx)
                     else
                      overlap(j)=overlap(j)+eig_vec(i,jsort)*ccsf(k,jj,1)*s_norm(idx)
                    endif
  166           continue
               else
                do 167 i=1,nparm
                  do 167 k=1,i
                    idx=idx+1
                    if(i.ne.k) then
                      overlap(j)=overlap(j)+(eig_vec(i,jsort)*cdet(k,jj,1)+eig_vec(k,jsort)*cdet(i,jj,1))*s_norm(idx)
                     else
                      overlap(j)=overlap(j)+eig_vec(i,jsort)*cdet(k,jj,1)*s_norm(idx)
                    endif
  167           continue
              endif

              dnorm=0.d0
              idx=0
              do 170 i=1,nparm
                do 170 k=1,i
                  idx=idx+1
                  dmul=1.d0
                  if(i.ne.k) dmul=2.d0
  170             dnorm=dnorm+dmul*eig_vec(i,jsort)*eig_vec(k,jsort)*s_norm(idx)
              dnorm=1.d0/dsqrt(dnorm)

              overlap(j)=dabs(overlap(j))*dnorm*dnorm_jj

              write(6,'(i4,'' eigenvalue '',i4,'' = '',f15.5,'' + i* '',f15.5)') j,jsort,eig(jsort),eigi(jsort)
              write(6,'('' overlap state,eigenstate '',2i4,'' = '',f15.5)') jj,j,overlap(j)
  
            endif
  172     continue

          target_overlap=-1.d0
          do 175 j=i0+1,nparm
            if(overlap(j).gt.target_overlap) then
              i0=j
              target_overlap=overlap(j)
            endif
  175     continue

          do 178 i=1,nparm
  178       grad(i+nparm*(jj-1))=eig_vec(i,isort(i0))*dnorm
          write(6,'(''state '',i4,'' norm'',1p1e12.5,'' overlap '',1p1e12.5)') jj,dnorm,overlap(i0)
          write(6,'(''pn  ='',1000f10.5)') (grad(i+nparm*(jj-1)),i=1,nparm)

          if(nstates.gt.1.and.jj.ne.nstates.and.eig(isort(i0+1)).eq.0.d0) 
     &      call fatal_error('OPTWF: Overlap with state 1 for highest eigenvalue >0')
  180   continue
       endif

      elseif(method.eq.'perturbative') then

       idx=0
       do 208 i=1,nparm
         grad(i)=grad_sav(i)
         do 208 j=1,i
           idx=idx+1
           s(i,j)=s_sav(idx)
           s(j,i)=s_sav(idx)
  208  continue

      call p2gtfd('optwf:svd_threshold',ftol,1.0d-4,1)
CVARDOC Threshold for keeping eigenvalues in singular value decomposition
      call svdchk(s,mparmx,nparm,grad,cdelta,nkept,ftol,1)
      write(6,'(''eigenvalues kept: '',i4,'' of '',i4)') nkept,nparm

c Scale orbital variations
      call p2gtid('optwf:mscaling',mscaling,0,1)
      call p2gtfd('optwf:scale',scale,1.d0,1)
      do 310 i=1,nparm
        cdelta(i)=cdelta(i)*scale
        if(mscaling.gt.0) then
          cdelta(i)=0.d0
        endif
  310 continue

      dwf=0
c     do 320 i=1,nparm
c       do 320 j=1,nparm
c 320     dwf=dwf+cdelta(i)*cdelta(j)*(s_sav(i,j)-s_sav(i,1)*s_sav(j,1))
      idx=0
      do 320 i=1,nparm
        do 319 j=1,i-1
          idx=idx+1
  319     dwf=dwf+2*cdelta(i)*cdelta(j)*s_sav(idx)
          idx=idx+1
  320     dwf=dwf+cdelta(i)*cdelta(i)*s_sav(idx)
      idx=1
      cdelta_s=0
      do 321 i=1,nparm
        idx=idx+i-1
  321   delta_s=cdelta_s+cdelta(i)*s_sav(idx)
      dwf=dwf-cdelta_s*cdelta_s

      dwf=sqrt(dwf)
      write(6,'(''wave function variation  ='',f10.5)') dwf
      call p2gtfd('optwf:wfscale',dwf_max,dwf,1)
      if(dwf_max.lt.dwf) then
        s1=sqrt(dwf_max/dwf)
        do 16 i=1,nparm
          cdelta(i)=cdelta(i)*s1
 16     enddo
        write(6,'(''wave function variation reset ='',f10.5)') dwf_max
      endif

      do 360 i=1,nparm
  360   grad(i)=cdelta(i)

      endif

      end
c-----------------------------------------------------------------------
      subroutine regularize_geneig(n,mparmx,h,s,work,seig_valinv,hmod)
      
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optjas.h'
      include 'optci.h'
      include 'optorb.h'
      include 'numbas.h'
      parameter(MPARMALL=MPARMJ+MXCIREDUCED+MXREDUCED,eps=1.d-12)
      parameter(MWORK=50*MPARMALL)
      parameter(eps_eigval=1.d-14)
      
      dimension h(mparmx,*),s(mparmx,*)
      dimension seig_vals(MPARMALL),seig_valinv(*)
      dimension hmod(MPARMALL,*),work(*)

      call cpu_time(t0)

c call dsyev to determine lworks
      call dsyev('V','U',n,s,mparmx,seig_vals,work,-1,isdinfo)
      lworks=work(1)
c     write(6,*) 'S diag opt lwork=',lworks

c diagonalize s=S -> S_diag=U^T S U -> in output, s contains the unitary matrix U
      call dsyev('V','U',n,s,mparmx,seig_vals,work,lworks,isdinfo)
c     call cpu_time(t)
c     t_sdiag=t 
c     write(6,*) 'elapsed time for diagonalization:',t_sdiag-t0

      icut=1
      ineg=1
      write(6,'(''overlap matrix: Maximum eigenvalue, threshold: '',1p2e12.4)') seig_vals(n),eps_eigval*seig_vals(n)
      do 10 i=1,n
        if((ineg.eq.1).and.(seig_vals(i).gt.0.d0)) then
          ineg=0
          write(6,'(''first positive eigenvalue'',t41,i6)') i
        endif
        if((icut.eq.1).and.(seig_vals(i).gt.eps_eigval*seig_vals(n))) then
          icut=0
          write(6,'(''first eigenvalue larger than threshold'',t41,i6)') i
        endif
  10  continue
      do 20 i=1,n
        if(seig_vals(i)/seig_vals(n).gt.eps_eigval) then
          seig_valinv(i)=1.0d0/dsqrt(seig_vals(i))
         else 
          seig_valinv(i)=0.0d0
        endif
  20  continue

c I = A^T S A where A_ij= U_ij/sqrt(s_eigval(j))
c s in input contains U -> s is overwritten and contains A 
      do 30 i=1,n
        do 30 j=1,n
          s(i,j)=s(i,j)*seig_valinv(j)
  30  continue

c Compute work_mat=A^T*H*A
      do 35 i=1,n
        do 35 l=1,n
  35      hmod(l,i)=0.d0
       
      do 50 l=1,n
        do 40 m=1,n
          work(m)=0.d0
          do 40 k=1,n
  40        work(m)=work(m)+s(k,l)*h(k,m)
        do 50 i=1,n
          do 50 m=1,n
  50        hmod(l,i)=hmod(l,i)+work(m)*s(m,i)

      call cpu_time(t)
c     t_hmod=t
c     write(6,*) 'elapsed time to build Hmod:',t_hmod-t_sdiag

      return
      end

c-----------------------------------------------------------------------
      subroutine solve_geneig(n,mparmx,hmod,s,seig_valinv,work,eig,eigi,eig_vec)
      
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optjas.h'
      include 'optci.h'
      include 'optorb.h'
      include 'numbas.h'
      parameter(MPARMALL=MPARMJ+MXCIREDUCED+MXREDUCED,eps=1.d-12)
      parameter(MWORK=50*MPARMALL)
      
      dimension seig_valinv(*)
      dimension hmod(mparmx,*),s(mparmx,*)

      dimension eig_vec(MPARMALL,*),eig_vecl(MPARMALL,1)
      dimension work(*)

c s_fordiag: a copy of S for diagonalization. 
c hmod: the modified Hamiltonian matrix (in the end, S^-1*U*H*U^T)
c s: overlap matrix, h: hamiltonian, eigenvec: eigenvectors, 

      call cpu_time(t0)

c MISSING
c hmod+adiag/s_diag

c Determine ilwork
      call dgeev('N','V',n,hmod,MPARMALL,eig,eigi,eig_vecl,
     &        MPARMALL,eig_vec,MPARMALL,work,-1,isdinfo)
      ilwork=work(1)
c     write(6,*) 'isdinfo, optimal lwork=',isdinfo,ilwork

c Diagonalize
      call dgeev('N','V',n,hmod,MPARMALL,eig,eigi,eig_vecl,
     &        MPARMALL,eig_vec,MPARMALL,work,ilwork,isdinfo)
c     write(6,*) 'isdinfo=',isdinfo
c     call cpu_time(t)
c     t_hmdiag=t
c     write(6,*) 'elapsed time to diagonalize Hmod:',t-t0

      do 20 k=1,n
        do 10 i=1,n
          work(i)=0.d0
          do 10 j=1,n
  10        work(i)=work(i)+s(i,j)*eig_vec(j,k)
        do 20 i=1,n
  20      eig_vec(i,k)=work(i)
      
c     call cpu_time(t)
c     t_eigvec=t
c     write(6,*) 'elapsed time to get eigenvectors:',t_eigvec-t_hmdiag

      return
      end
c-----------------------------------------------------------------------
      subroutine test_solution_parm(nparm,grad,
     &              dparm_norm,dparm_norm_min,add_diag,iflag)

      implicit real*8(a-h,o-z)

      include 'vmc.h'

      dimension grad(*)

      iflag=0
      if(add_diag.le.0.d0) return

c Calculate rms change in parameters
      dparm_norm=0
      do 30 i=1,nparm
  30    dparm_norm=dparm_norm+grad(i)**2
      dparm_norm=sqrt(dparm_norm/nparm)

      write(6,'(''dparm_norm,adiag ='',3g12.5)') 
     &dparm_norm,add_diag

      if(dparm_norm.gt.dparm_norm_min) iflag=1
      
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
      subroutine setup_wf
      implicit real*8(a-h,o-z)
  
      do 10 k=2,3
        call copy_jastrow(k)
        call copy_lcao(k)
  10    call copy_ci(k)

      return
      end
c-----------------------------------------------------------------------
      subroutine save_wf
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      if(ioptjas.ne.0) call save_jastrow
      if(ioptorb.ne.0) call save_lcao 
      if(ioptci.ne.0) call save_ci

      return
      end
c-----------------------------------------------------------------------
      subroutine restore_wf(iadiag)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      if(ioptjas.ne.0) call restore_jastrow(iadiag)
      if(ioptorb.ne.0) call restore_lcao(iadiag)
      if(ioptci.ne.0) call restore_ci(iadiag)

      return
      end
c-----------------------------------------------------------------------
      subroutine save_wf_best(ioptjas,ioptorb,ioptci)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      if(ioptjas.ne.0) call save_jastrow_best
      if(ioptorb.ne.0) call save_lcao_best
      if(ioptci.ne.0) call save_ci_best

      return
      end
c-----------------------------------------------------------------------
      subroutine save_jastrow
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(83,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc
      common /bparm/ nspin2b,nocuspb

      dimension a4_save(MORDJ1,MCTYPE,MWF),b_save(MORDJ1,2,MWF),
     &c_save(83,MCTYPE,MWF)

      save a4_save,b_save,c_save
      save mparmja,mparmjb,mparmjc

c Save parameters corresponding to run generating hessian

      mparmja=2+max(0,norda-1)
      mparmjb=2+max(0,nordb-1)
      mparmjc=nterms4(nordc)

      do 50 ict=1,nctype
        do 50 i=1,mparmja
   50     a4_save(i,ict,1)=a4(i,ict,1)
      do 60 i=1,mparmjb
   60   b_save(i,1,1)=b(i,1,1)
      do 70 ict=1,nctype
        do 70 i=1,mparmjc
   70     c_save(i,ict,1)=c(i,ict,1)

      return

      entry restore_jastrow(iadiag)

c Restore parameters corresponding to run generating hessian
      do 80 ict=1,nctype
        do 80 i=1,mparmja
   80     a4(i,ict,iadiag)=a4_save(i,ict,1)
      do 90 i=1,mparmjb
   90   b(i,1,iadiag)=b_save(i,1,1)
      do 100 ict=1,nctype
        do 100 i=1,mparmjc
  100     c(i,ict,iadiag)=c_save(i,ict,1)

      return
      end

c-----------------------------------------------------------------------
      subroutine save_lcao
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb

      dimension coef_save(MBASIS,MORB,MWF)

      save coef_save

      do 10 i=1,norb
       do 10 j=1,nbasis
   10   coef_save(j,i,1)=coef(j,i,1)

      return

      entry restore_lcao(iadiag)

      do 20 i=1,norb
       do 20 j=1,nbasis
   20   coef(j,i,iadiag)=coef_save(j,i,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine save_ci
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb

      dimension cdet_save(MDET,MSTATES),ccsf_save(MDET,MSTATES)
      save cdet_save,ccsf_save

      do 10 j=1,nstates
        do 10 i=1,ndet
   10     cdet_save(i,j)=cdet(i,j,1)

      do 20 j=1,nstates
       do 20 icsf=1,ncsf
   20   ccsf_save(icsf,j)=ccsf(icsf,j,1)

      return

      entry restore_ci(iadiag)

      do 30 j=1,nstates
        do 30 i=1,ndet
   30     cdet(i,j,iadiag)=cdet_save(i,j)

      do 40 j=1,nstates
       do 40 icsf=1,ncsf
   40   ccsf(icsf,j,iadiag)=ccsf_save(icsf,j)

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_jastrow(iadiag)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(83,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc
      common /bparm/ nspin2b,nocuspb

      mparmja=2+max(0,norda-1)
      mparmjb=2+max(0,nordb-1)
      mparmjc=nterms4(nordc)

      scalek(iadiag)=scalek(1)
      do 80 ict=1,nctype
        do 80 i=1,mparmja
   80     a4(i,ict,iadiag)=a4(i,ict,1)
      do 90 i=1,mparmjb
   90   b(i,1,iadiag)=b(i,1,1)
      do 100 ict=1,nctype
        do 100 i=1,mparmjc
  100     c(i,ict,iadiag)=c(i,ict,1)

      return
      end

c-----------------------------------------------------------------------
      subroutine copy_lcao(iadiag)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb

      do 20 i=1,norb+nadorb
       do 20 j=1,nbasis
   20   coef(j,i,iadiag)=coef(j,i,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_ci(iadiag)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      do 30 j=1,nstates
        do 30 i=1,ndet
   30     cdet(i,j,iadiag)=cdet(i,j,1)

      do 40 j=1,nstates
       do 40 icsf=1,ncsf
   40   ccsf(icsf,j,iadiag)=ccsf(icsf,j,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_zex(iadiag)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'basis.h'

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb

      do 20 i=1,nbasis
   20   zex(i,iadiag)=zex(i,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine save_jastrow_best
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(83,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc
      common /bparm/ nspin2b,nocuspb

      dimension a4_best(MORDJ1,MCTYPE,MWF),b_best(MORDJ1,2,MWF),
     &c_best(83,MCTYPE,MWF)

      save a4_best,b_best,c_best
      save mparmja,mparmjb,mparmjc

c Save parameters corresponding to run generating hessian

      mparmja=2+max(0,norda-1)
      mparmjb=2+max(0,nordb-1)
      mparmjc=nterms4(nordc)

      do 50 ict=1,nctype
        do 50 i=1,mparmja
   50     a4_best(i,ict,1)=a4(i,ict,1)
      do 60 i=1,mparmjb
   60   b_best(i,1,1)=b(i,1,1)
      do 70 ict=1,nctype
        do 70 i=1,mparmjc
   70     c_best(i,ict,1)=c(i,ict,1)

      return

      entry restore_jastrow_best

c Restore parameters corresponding to run generating hessian
      do 80 ict=1,nctype
        do 80 i=1,mparmja
   80     a4(i,ict,1)=a4_best(i,ict,1)
      do 90 i=1,mparmjb
   90   b(i,1,1)=b_best(i,1,1)
      do 100 ict=1,nctype
        do 100 i=1,mparmjc
  100     c(i,ict,1)=c_best(i,ict,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine save_lcao_best
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      dimension coef_best(MBASIS,MORB,MWF)

      save coef_best

      do 10 i=1,norb
       do 10 j=1,nbasis
   10   coef_best(j,i,1)=coef(j,i,1)

      return

      entry restore_lcao_best

      if(ioptorb.eq.0) return

      do 20 i=1,norb
       do 20 j=1,nbasis
   20   coef(j,i,1)=coef_best(j,i,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine save_ci_best
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      dimension cdet_best(MDET,MSTATES),ccsf_best(MDET,MSTATES)
      save cdet_best,ccsf_best

      do 10 j=1,nstates
        do 10 i=1,ndet
   10     cdet_best(i,j)=cdet(i,j,1)

      do 20 j=1,nstates
       do 20 icsf=1,ncsf
   20   ccsf_best(icsf,j)=ccsf(icsf,j,1)

      return

      entry restore_ci_best

      if(ioptci.eq.0) return

      do 30 j=1,nstates
        do 30 i=1,ndet
   30     cdet(i,j,1)=cdet_best(i,j)

      do 40 j=1,nstates
       do 40 icsf=1,ncsf
   40   ccsf(icsf,j,1)=ccsf_best(icsf,j)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_parameters(grad,iflag,iadiag)
      implicit real*8(a-h,o-z)

      dimension grad(*)

      iflag=0
      call compute_jastrow(grad,iflag,iadiag)

      if(iflag.ne.0) return

      call compute_lcao(grad,iadiag)

      call compute_ci(grad,iadiag)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_jastrow(grad,iflag,iadiag)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(83,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc
      common /bparm/ nspin2b,nocuspb

      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj
      common /optwf_wjas/ iwjasa(83,MCTYP3X),iwjasb(83,3),iwjasc(83,MCTYPE),iwjasf(15,MCTYPE)
      common /optwf_nparmj/ nparma(MCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE)

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      dimension grad(*)

      if(ioptjas.eq.0) return

c Set up cusp conditions
      call cuspinit4(0)

c Add change to old parameters
      iparm=0
      do 50 ict=1,nctype
        do 50 i=1,nparma(ict)
          iparm=iparm+1
   50     a4(iwjasa(i,ict),ict,iadiag)=a4(iwjasa(i,ict),ict,iadiag)-grad(iparm)
      do 60 i=1,nparmb(1)
        iparm=iparm+1
   60   b(iwjasb(i,1),1,iadiag)=b(iwjasb(i,1),1,iadiag)-grad(iparm)
      do 70 ict=1,nctype
        do 70 i=1,nparmc(ict)
          iparm=iparm+1
   70     c(iwjasc(i,ict),ict,iadiag)=c(iwjasc(i,ict),ict,iadiag)-grad(iparm)
      call cuspexact4(0,iadiag)

c Check parameters a2 and b2 > -scalek
      call check_parms_jas(iflag)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_lcao(grad,iadiag)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'
      include 'optorb_cblk.h'

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj

      dimension acoef(MBASIS,MORB),grad(*)

      if(ioptorb.eq.0) return

      do 10 i=1,norb
       do 10 j=1,nbasis
 10     acoef(j,i)=coef(j,i,iadiag)

c Update the orbitals
      do 30 i=1,norbterm
       io=ideriv(1,i)
       jo=ideriv(2,i)
       do 30 j=1,nbasis
 30     acoef(j,io)=acoef(j,io)-grad(i+nparmj+nparmd)*coef(j,jo,iadiag)

      do 50 i=1,norb
       do 50 j=1,nbasis
 50     coef(j,i,iadiag)=acoef(j,i)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_ci(grad,iadiag)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj

      dimension grad(*)

      if(ioptci.eq.0) return

c Update the ci coef
      if((method.eq.'linear'.or.method.eq.'lin_d').and.ioptjas+ioptorb.eq.0) then
        do 31 k=1,nstates

          if(ncsf.eq.0) then
            do 10 idet=1,ndet
              cdet(idet,k,1)=grad(idet+ndet*(k-1))
 10         continue
           else
            do 15 j=1,ndet
 15           cdet(j,k,1)=0
            do 30 icsf=1,ncsf
              do 20 j=iadet(icsf),ibdet(icsf)
                jx=icxdet(j)
                cdet(jx,k,1)=cdet(jx,k,1)+grad(icsf+ncsf*(k-1))*cxdet(j)
 20           continue
              ccsf(icsf,k,1)=grad(icsf+ncsf*(k-1))
 30         continue
          endif

 31     continue
       else
         if(ncsf.eq.0) then
           do 35 idet=2,ndet
             cdet(idet,1,iadiag)=cdet(idet,1,iadiag)-grad(idet-1+nparmj)
 35        continue
          else
           do 50 icsf=2,ncsf
             do 40 j=iadet(icsf),ibdet(icsf)
               jx=icxdet(j)
               cdet(jx,1,iadiag)=cdet(jx,1,iadiag)-grad(icsf-1+nparmj)*cxdet(j)
 40          continue
             ccsf(icsf,1,iadiag)=ccsf(icsf,1,iadiag)-grad(icsf-1+nparmj)
 50        continue
         endif
      endif

c     do 90 j=1,nstates
c90     write(6,'(''csf ='',1000f20.15)') (ccsf(i,j,iadiag),i=1,ncsf)

      return
      end
c-----------------------------------------------------------------------
      subroutine check_parms_jas(iflag)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(83,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc
      common /bparm/ nspin2b,nocuspb

      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj
      common /optwf_wjas/ iwjasa(83,MCTYP3X),iwjasb(83,3),iwjasc(83,MCTYPE),iwjasf(15,MCTYPE)
      common /optwf_nparmj/ nparma(MCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE)

      iflag=0
      iflaga=0
      iflagb=0

      scalem=-scalek(1)
      do 50 ict=1,nctype
        do 50 i=1,nparma(ict)
   50     if(iwjasa(i,ict).eq.2.and.a4(2,ict,1).le.scalem) iflaga=1
      if(iflaga.eq.1) then
        do 55 ict=1,nctype
   55     write(6,'(''a2 < -scalek'',f10.5)') a4(2,ict,1)
      endif
      do 60 i=1,nparmb(1)
   60   if(iwjasb(i,1).eq.2.and.b(2,1,1).le.scalem) iflagb=1
      if(iflagb.eq.1) write(6,'(''b2 < -scalek'',f10.5)') b(2,1,1)
      
      if(iflaga.eq.1.or.iflagb.eq.1) iflag=1

      return
      end
c-----------------------------------------------------------------------
      subroutine quad_min

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      parameter(MFUNC=3)

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_corsam/ add_diag(MFORCE),energy(MFORCE),energy_err(MFORCE),force(MFORCE),force_err(MFORCE)

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
      subroutine write_wf(iwf_fit,iter)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      character*40 filetype,wf,itn

      if(iter.lt.0) then
        filetype='_optimal.'//wf(1:index(wf,' ')-1)
       else
        write(wf,'(i1)') iwf_fit
        if(iter.lt.10) then
          write(itn,'(i1)') iter
         elseif(iter.lt.100) then
          write(itn,'(i2)') iter
         elseif(iter.lt.1000) then
          write(itn,'(i3)') iter
        endif
        filetype='_optimal.'//wf(1:index(wf,' ')-1)//'.iter'//itn(1:index(itn,' ')-1)
      endif

      call write_jastrow(iwf_fit,filetype)
      call write_lcao(iwf_fit,filetype)
      call write_ci(iwf_fit,filetype)

      return
      end
c-----------------------------------------------------------------------
      subroutine write_wf_best
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      call restore_jastrow_best
      call restore_lcao_best
      call restore_ci_best

      call write_wf(1,-1)

      return
      end
c-----------------------------------------------------------------------
      subroutine write_jastrow(iwf_fit,filetype)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      character*50 fmt
      character*40 filename,filetype

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(83,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc
      common /bparm/ nspin2b,nocuspb

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj
      common /optwf_wjas/ iwjasa(83,MCTYP3X),iwjasb(83,3),iwjasc(83,MCTYPE),iwjasf(15,MCTYPE)
      common /optwf_nparmj/ nparma(MCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE)

      if(ioptjas.eq.0) return

      filename='jastrow'//filetype(1:index(filetype,' ')-1)

      open(2,file=filename,status='unknown')

      call p2gtid('jastrow:ianalyt_lap',ianalyt_lap,1,1)
      call p2gti('jastrow:ijas',ijas,1)
      call p2gti('jastrow:isc',isc,1)
      call p2gtid('jastrow:nspin1',nspin1,1,1)
      call p2gtid('jastrow:nspin2',nspin2,1,1)
      call p2gtid('jastrow:ifock',ifock,0,1)
      write(2,'(''&jastrow ianalyt_lap'',i2,'' ijas'',i2,'' isc'',i2,
     &'' nspin1'',i2,'' nspin2'',i2,'' ifock'',i2)') ianalyt_lap,ijas,isc,nspin1,nspin2,ifock
      write(2,*)
      write(2,'(''jastrow_parameter'',i4)') iwf_fit
      write(2,'(3i3,a28)') norda,nordb,nordc,' norda,nordb,nordc'
c tmp
      a21=0
      write(2,'(2f13.8,a15)') scalek(1),a21,' scalek,a21'
      mparmja=2+max(0,norda-1)
      mparmjb=2+max(0,nordb-1)
      mparmjc=nterms4(nordc)
      if(mparmja.gt.0) then
        write(fmt,'(''(''i2,''f13.8,a28)'')') mparmja
       else
        write(fmt,'(''(a28)'')')
      endif
      do 80 ict=1,nctype
   80   write(2,fmt) (a4(i,ict,1),i=1,mparmja),' (a(iparmj),iparmj=1,nparma)'

      if(mparmjb.gt.0) then
        write(fmt,'(''(''i2,''f13.8,a28)'')') mparmjb
       else
        write(fmt,'(''(a28)'')')
      endif
      write(2,fmt) (b(i,1,1),i=1,mparmjb),' (b(iparmj),iparmj=1,nparmb)'

      if(mparmjc.gt.0) then
        write(fmt,'(''(''i2,''f13.8,a28)'')') mparmjc
       else
        write(fmt,'(''(a28)'')')
      endif
      do 90 ict=1,nctype
   90   write(2,fmt) (c(i,ict,1),i=1,mparmjc),' (c(iparmj),iparmj=1,nparmc)'
      write(2,'(''end'')')
      close(2)

      return
      end
c-----------------------------------------------------------------------
      subroutine write_lcao(iwf_fit,filetype)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'numbas.h'

      common /numbas/ arg(MCTYPE),r0(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS,MCTYPE)

      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      character*40 filename,filetype

      dimension anorm(MBASIS)

      if(ioptorb.eq.0) return

      filename='orbitals'//filetype(1:index(filetype,' ')-1)
      open(2,file=filename,status='unknown')
      write(2,'(''lcao '',3i4)') norb+nadorb,nbasis,iwf_fit

      call p2gtfd('general:scalecoef',scalecoef,1.0d0,1)
      if(numr.gt.0) then
        do 20 i=1,norb+nadorb
   20     write(2,'(1000e20.8)') (coef(j,i,1)/scalecoef,j=1,nbasis)
      else
        call basis_norm(1,anorm,1)
        do 40 i=1,norb+nadorb
   40     write(2,'(1000e20.8)') (coef(j,i,1)/(anorm(j)*scalecoef),j=1,nbasis)
      endif

      write(2,'(''end'')')
      close(2)

      return
      end
c-----------------------------------------------------------------------
      subroutine write_ci(iwf_fit,filetype)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      character*40 filename,filetype

      common /dorb/ iworbd(MELEC,MDET)

      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj

      if(ioptci.eq.0) return

      filename='det'//filetype(1:index(filetype,' ')-1)
      open(2,file=filename,status='unknown')

      call p2gti('electrons:nelec',nelec,1)
      call p2gti('electrons:nup',nup,1)
      write(2,'(''&electrons nelec '',i4,'' nup '',i4)') nelec,nup
      do 1 istate=1,nstates
      write(2,'(''# State '',i4)') istate
      write(2,'(''determinants'',i10,i4)') ndet,iwf_fit
      write(2,'(100f15.8)') (cdet(i,istate,1),i=1,ndet)
      do 1 k=1,ndet
   1   write(2,'(100i4)') (iworbd(i,k),i=1,nelec)
 
      write(2,'(''end'')')

      if(ncsf.ne.0) then
        write(2,'(''csf '',i10,i4)') ncsf,nstates
        do i=1,nstates
          write(2,'(100f15.8)') (ccsf(j,i,1),j=1,ncsf)
        enddo
        write(2,'(''end'')')
c
        nmap=0
        do 5 i=1,ncsf
   5      nmap=nmap+ibdet(i)-iadet(i)+1
        write(2,'(''csfmap'')') 
        write(2,'(3i10)') ncsf,ndet,nmap
        nptr=0
        do 10 i=1,ncsf
         nterm=ibdet(i)-iadet(i)+1
         write(2,'(i10)') ibdet(i)-iadet(i)+1
         do 12 j=1,nterm
           nptr=nptr+1
           write(2,'(i10,f10.6)') icxdet(nptr),cxdet(nptr)
  12     enddo
  10    enddo
        write(2,'(''end'')')
      endif

      close(2)

      return
      end
c-----------------------------------------------------------------------
      subroutine combine_derivatives

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'optjas.h'
      include 'optci.h'
      include 'optorb.h'
        
      parameter(MPARMALL=MPARMJ+MXCIREDUCED+MXREDUCED)

      common /gradhess_jas/ grad_jas(MPARMJ),h_jas(MPARMJ,MPARMJ),s_jas(MPARMJ,MPARMJ)
      common /gradhess_ci/  grad_ci(MXCITERM),h_ci(MXCITERM,MXCIREDUCED),s_ci(MXCITERM,MXCIREDUCED)
c     common /gradhess_orb/ grad_orb(MXORBOP),h_orb(MXMATDIM),s_orb(MXMATDIM)

      common /gradhess_mix_jas_ci/  h_mix_jas_ci(2*MPARMJ,MXCITERM),s_mix_jas_ci(MPARMJ,MXCITERM)
      common /gradhess_mix_jas_orb/ h_mix_jas_orb(2*MPARMJ,MXREDUCED),s_mix_jas_orb(MPARMJ,MXREDUCED)
      common /gradhess_mix_orb_ci/  h_mix_ci_orb(2*MXCITERM,MXREDUCED),s_mix_ci_orb(MXCITERM,MXREDUCED)

      common /gradhess_all/ grad(MPARMALL),h(MPARMALL,MPARMALL),s(MPARMALL,MPARMALL)

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj

c Note: we do not vary the first (i0) CI coefficient unless full CI

      if(method.eq.'hessian') then

        nparmd=max(nciterm-1,0)

c Jastrow Hessian
        do 10 j=1,nparmj
          grad(j)=grad_jas(j)
          do 10 i=1,nparmj
  10        h(i,j)=h_jas(i,j)
        
c CI Hessian
        do 20 j=1,nparmd
          grad(j+nparmj)=grad_ci(j+1)
          do 20 i=1,nparmd
  20        h(i+nparmj,j+nparmj)=h_ci(i+1,j+1)

c Jastrow-CI Hessian
        do 30 j=1,nparmd
          do 30 i=1,nparmj
            h(i,j+nparmj)=h_mix_jas_ci(i,j+1)
  30        h(j+nparmj,i)=h_mix_jas_ci(i,j+1)

c ORB Hessian
c       do 40 j=1,nreduced
c         jk=(j-1)*nreduced
c         grad(j+nparmj+nparmd)=grad_orb(j)
c         do 40 i=1,nreduced
c 40        h(i+nparmj+nparmd,j+nparmj+nparmd)=h_orb(i+jk)

c Jastrow-ORB Hessian
        do 50 j=1,nreduced
          do 50 i=1,nparmj
            h(i,j+nparmj+nparmd)=h_mix_jas_orb(i,j)
  50        h(j+nparmj+nparmd,i)=h_mix_jas_orb(i,j)

c CI-ORB Hessian
        do 60 j=1,nreduced
          do 60 i=1,nparmd
            h(i+nparmj,j+nparmj+nparmd)=h_mix_ci_orb(i+1,j)
  60        h(j+nparmj+nparmd,i+nparmj)=h_mix_ci_orb(i+1,j)

        nparm=nparmj+nparmd+nreduced

c       do 70 i=1,nparm
c 70      write(6,'(''h = '',1000f10.5)') (h(i,j),j=1,nparm)

      elseif(method.eq.'linear') then
        
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

      elseif(method.eq.'perturbative') then

c ORB gradient and overlap
c       do 220 j=1,nreduced
c         grad(j)=grad_orb(j)
c         jk=(j-1)*nreduced
c         do 220 i=1,nreduced
c 220       s(i,j)=s_orb(i+jk)

       nparm=nreduced
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine save_nparms

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'optci.h'
      include 'optorb.h'

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj

      save nparmj_sav,norbterm_sav,nciterm_sav,nparmd_sav,nreduced_sav

      nparmj_sav=nparmj
      norbterm_sav=norbterm
      nreduced_sav=nreduced
      nciterm_sav=nciterm
      nparmd=max(nciterm-1,0)
      nparmd_sav=nparmd

      write(6,'(''Saved max number of parameters, nparmj,norb,nciterm,nciterm-1: '',5i5)') nparmj,norbterm,nciterm,nparmd
      return

      entry set_nparms

      nparmj=nparmj_sav
      nparmd=nparmd_sav
      norbterm=norbterm_sav
      nreduced=nreduced_sav
      nciterm=nciterm_sav

      if(ioptjas.eq.0) nparmj=0
      if(ioptorb.eq.0) then
        norbterm=0
        nreduced=0
      endif
      if(ioptci.eq.0) then
        nciterm=0
        nparmd=0
      endif

      write(6,'(''Max number of parameters set to nparmj,norb,nciterm,nciterm-1: '',5i5)') nparmj,norbterm,nciterm,nparmd
      call set_nparms_tot

      return
      end
c-----------------------------------------------------------------------
      subroutine set_nparms_tot

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'optjas.h'
      include 'optci.h'
      include 'optorb.h'
        
      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj

c Note: we do not vary the first (i0) CI coefficient unless a run where we only optimize the CI coefs

      if(method.eq.'hessian'.or.method.eq.'sr_n') then

        nparmd=max(nciterm-1,0)
        nparm=nparmj+nparmd+norbterm

      elseif(method.eq.'linear'.or.method.eq.'lin_d') then
        
       i0=0
       if(ioptci.ne.0) i0=1
       if(ioptjas.eq.0.and.ioptorb.eq.0) i0=0

       nparmd=max(nciterm-1,0)
       nparm=nparmj+norbterm+nciterm-i0

      elseif(method.eq.'perturbative') then

       nparm=norbterm
      endif

      write(6,'(/,''number of parms: total, Jastrow, CI, orbitals= '',4i5)') 
     & nparm,nparmj,nciterm,norbterm

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
