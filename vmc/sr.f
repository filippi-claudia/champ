      subroutine sr_optwf

      implicit real*8 (a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'sr.h'

      character*20 method_sav

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contrl/ nstep,nblk,nblkeq,nconf_old,nconf_new,isite,idump,irstar
      common /optwf_corsam/ add_diag(MFORCE),energy(MFORCE),energy_err(MFORCE),force(MFORCE),force_err(MFORCE)
      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_func/ omega,ifunc_omega

      common /force_analy/ iforce_analy,iuse_zmat,alfgeo
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      dimension grad(MPARM*MSTATES)

      save method_sav

      if(method.ne.'sr_n'.and.method.ne.'lin_d'.and.method.ne.'mix_n')return

      if(method.eq.'mix_n'.and.(nstates.eq.1.or.ioptci.eq.0)) call fatal_error('SR_OPTWF: no need mix_n')

      method_sav=method

      call set_nparms_tot

      if(nparm.gt.MPARM)call fatal_error('SR_OPTWF: nparmtot gt MPARM')

      call p2gtid('optwf:nopt_iter',nopt_iter,6,1)
      call p2gtid('optwf:nblk_max',nblk_max,nblk,1)
      call p2gtfd('optwf:energy_tol',energy_tol,1.d-3,1)

      call p2gtfd('optwf:dparm_norm_min',dparm_norm_min,1.0d0,1)
      write(6,'(''Starting dparm_norm_min'',g12.4)') dparm_norm_min

      call p2gtfd('optwf:sr_tau',sr_tau,0.02,1)
      call p2gtfd('optwf:sr_adiag',sr_adiag,0.01,1)
      call p2gtfd('optwf:sr_eps',sr_eps,0.001,1)

      call p2gtid('optwf:lin_nvec',nvec,5,1)
      call p2gtid('optwf:lin_nvecx',nvecx,MVEC,1)
      call p2gtfd('optwf:lin_adiag',alin_adiag,0.01,1)
      call p2gtfd('optwf:lin_eps',alin_eps,0.001,1)

      call p2gtid('optwf:func_omega',ifunc_omega,0,1)
      if(ifunc_omega.gt.0) then
       call p2gtfd('optwf:omega',omega,0.d0,1)
       write(6,'(/,''LIN_D omega: '',f10.5)') omega
      endif

      call p2gtid('optwf:micro_iter_sr',micro_iter_sr,1,1)

      call p2gtid('optgeo:izvzb',izvzb,0,1)

      if(method.eq.'mix_n'.and.iforce_analy.gt.0) call p2gtid('optgeo:iroot_geo',iroot_geo,0,0)

      if(method.eq.'sr_n'.or.method.eq.'mix_n') then

        write(6,'(/,''SR adiag: '',f10.5)') sr_adiag
        write(6,'(''SR tau:   '',f10.5)') sr_tau
        write(6,'(''SR eps:   '',f10.5)') sr_eps
      endif

      if(method.eq.'lin_d'.or.method.eq.'mix_n') then

        if(nvecx.gt.MVEC) call fatal_error('SR_OPTWF: nvecx > MVEC')
        write(6,'(/,''LIN_D adiag: '',f10.5)') alin_adiag
        write(6,'(''LIN_D ethr:  '',f10.5)') alin_eps
        write(6,'(''LIN_D nvec:  '',i4)') nvec
        write(6,'(''LIN_D nvecx: '',i4)') nvecx

        if(nstates.gt.1.and.nvec.lt.nstates) call fatal_error('SR_OPTWF: nvec < nstates')
      endif

      inc_nblk=0
      
      sr_adiag_sav=sr_adiag
      alin_adiag_sav=alin_adiag

      nstates_sav=nstates
      iforce_analy_sav=iforce_analy

      ioptjas_sav=ioptjas
      ioptorb_sav=ioptorb
      ioptci_sav=ioptci
      call save_nparms

      call write_geometry(0)

c do iteration
      do iter=1,nopt_iter
        write(6,'(/,''Optimization iteration'',i5,'' of'',i5)')iter,nopt_iter

        if(ipr.gt.1) write(88,'(/,''Optimization iteration'',i5,'' of'',i5)')iter,nopt_iter
 
        iforce_analy=0

c do micro_iteration
        do miter=1,micro_iter_sr

          if(micro_iter_sr.gt.1) write(6,'(/,''Micro iteration'',i5,'' of'',i5)')miter,micro_iter_sr

          if(method_sav.eq.'mix_n') then
            method='sr_n'
            ioptci=0
            ioptorb=ioptorb_sav
            ioptjas=ioptjas_sav

            if(miter.eq.micro_iter_sr) then
              ioptci=ioptci_sav
              ioptorb=0
              ioptjas=0
              method='lin_d'
            endif
            call set_nparms

           else
            if(miter.eq.micro_iter_sr) iforce_analy=iforce_analy_sav
          endif

          call qmc

          write(6,'(/,''Completed sampling'')')

   6      continue

          if(method.eq.'sr_n') then
            call sr(nparm,grad,sr_adiag,sr_eps,i)
            call dscal(nparm,-sr_tau,grad,1)
           else
            call lin_d(nparm,nvec,nvecx,grad,alin_adiag,alin_eps)
            if(nstates.eq.1) call dscal(nparm,-1.d0,grad,1)
          endif

          if(method.eq.'sr_n'.or.(method.eq.'lin_d'.and.ioptorb+ioptjas.gt.0)) then
            if(method.eq.'sr_n') then
              adiag=sr_adiag
             else
              adiag=alin_adiag
            endif
            call test_solution_parm(nparm,grad,dparm_norm,dparm_norm_min,adiag,iflag)
            write(6,'(''Norm of parm variation '',g12.5)') dparm_norm
            if(iflag.ne.0) then
              write(6,'(''Warning: dparm_norm>1'')')
              adiag=10*adiag
              write(6,'(''adiag increased to '',f10.5)') adiag

              sr_adiag=adiag
              alin_adiag=adiag
              go to 6
             else
              sr_adiag=sr_adiag_sav
              alin_adiag=alin_adiag_sav
            endif
          endif
 
          call compute_parameters(grad,iflag,1)
          call write_wf(1,iter)

          call save_wf

c enters only if method_sav.not.mix_n
          if(iforce_analy.gt.0) then

            if(izvzb.gt.0) call forces_zvzb(nparm)

            call compute_positions
            call write_geometry(iter)
          endif
        enddo
c enddo micro_iteration
 
        if(method_sav.eq.'mix_n'.and.iforce_analy_sav.eq.1) then

          call select_ci_root(iroot_geo)

          iforce_analy=iforce_analy_sav
          iguiding_sav=iguiding

          nstates=1
          iguiding=0

          ioptjas=0
          ioptorb=0
          ioptci=0

          call set_nparms
          call qmc
          
          call compute_positions
          call write_geometry(iter)

          nstates=nstates_sav
          iguiding=iguiding_sav

          ioptjas=ioptjas_sav
          ioptorb=ioptorb_sav
          ioptci=ioptci_sav

          call restore_wf(1)
        endif
          
        if(iter.ge.2) then
          denergy=energy(1)-energy_sav
          denergy_err=sqrt(energy_err(1)**2+energy_err_sav**2)
c         call check_length_run_sr(iter,inc_nblk,nblk,nblk_max,denergy,denergy_err,energy_err_sav,energy_tol)
          nblk=nblk*1.2
          nblk=min(nblk,nblk_max)
c         if(-denergy.gt.3*denergy_err) alfgeo=alfgeo/1.2
        endif
        write(6,'(''nblk = '',i6)') nblk
        write(6,'(''alfgeo = '',f10.4)') alfgeo

        energy_sav=energy(1)
        energy_err_sav=energy_err(1)
      enddo
c enddo iteration

      write(6,'(/,''Check last iteration'')')

      ioptjas=0
      ioptorb=0
      ioptci=0
      iforce_analy=0

      call set_nparms

      call qmc

      call write_wf(1,-1)
      call write_geometry(-1)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sr(nparm,deltap,sr_adiag,sr_eps,i)
c solve S*deltap=h_sr (call in optwf)

      implicit real*8(a-h,o-z)
      include 'mstates.h'
      include 'sr.h'

      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES)
     &,s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,nconf

      dimension deltap(nparm)

      call sr_hs(nparm,sr_adiag)

      imax=nparm          ! max n. iterations conjugate gradients
      imod=50             ! inv. freq. of calc. r=b-Ax vs. r=r-\alpha q (see pcg)
      do i=1,nparm   
       deltap(i)=0.d0     ! initial guess of solution
      enddo
      call pcg(nparm,h_sr,deltap,i,imax,imod,sr_eps)
      write(6,*) 'CG iter ',i

      return              ! deltap
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sr_store(l,wt,psid,energy)
c store elocal and derivatives of psi for each configuration (call in vmc)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optjas.h'
      include 'optorb.h'
      include 'optorb_cblk.h'
      include 'optci.h'
      include 'optci_cblk.h'
      include 'sr.h'

      common /force_analy/ iforce_analy,iuse_zmat,alfgeo

      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent

      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj
      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      common /derivjas/ gvalue(MPARMJ),g(3,MELEC,MPARMJ)
     &,d2g(MPARMJ),go(MELEC,MELEC,MPARMJ)
      common /deloc_dj/ denergy(MPARMJ,MSTATES)

      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES)
     &,s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,nconf

      dimension tmp_ho(MPARMJ),wt(*),psid(*),energy(*)

      call p2gtid('optgeo:izvzb',izvzb,0,1)

      if(iforce_analy.gt.0.and.izvzb.eq.1) call force_store(l)

      if((method.ne.'sr_n'.and.method.ne.'lin_d').or.ioptjas+ioptorb+ioptci.eq.0)return

      i0=1
      if(method.eq.'lin_d'.and.ioptjas+ioptorb.eq.0) i0=0

      if(l.gt.MCONF) call fatal_error('SR_STORE: l gt MCONF')

      call dcopy(nparmj,gvalue,1,sr_o(1,l),1)

      ntmp=max(nciterm-i0,0)
      call dcopy(ntmp,ci_o(1+i0),1,sr_o(nparmj+1,l),1)

      ijasci=nparmj+ntmp
      if(ijasci+nstates*norbterm+nstates.gt.MPARM) call fatal_error('SR_STORE: iparm gt MPARM')

      do istate=1,nstates
        ii=ijasci+(istate-1)*norbterm
        call dcopy(norbterm,orb_o(1,istate),1,sr_o(ii+1,l),1)
        elocal(l,istate)=energy(istate)
        wtg(l,istate)=wt(istate)
      enddo

      ii=ijasci+nstates*norbterm
      do istate=1,nstates
        sr_o(ii+istate,l)=psid(istate)
      enddo
      
      nconf=l

      if(method.eq.'sr_n'.and.izvzb.eq.0) return

c TO FIX: we are assuming optjas.ne.0 or optorb.ne.0 -> Otherwise, standard secular problem

      do 10 j=1,nparmj
  10    tmp_ho(j)=denergy(j,1)+gvalue(j)*energy(1)

      call dcopy(nparmj,tmp_ho,1,sr_ho(1,l),1)

      call dcopy(ntmp,ci_e(1+i0),1,sr_ho(nparmj+1,l),1)

      call dcopy(norbterm,orb_ho(1,1),1,sr_ho(nparmj+ntmp+1,l),1)
      
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine force_store(l)

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'sr.h'

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent

      common /da_energy_now/ da_energy(3,MCENT),da_psi(3,MCENT)

      common /force_mat_n/ force_o(6*MCENT,MCONF)

      ii=0
      do 10 i=1,ncent
        do 10 k=1,3
          ii=ii+1
  10      force_o(ii,l)=da_psi(k,i)

      do 30 i=1,ncent
        do 30 k=1,3
          ii=ii+1
  30      force_o(ii,l)=da_energy(k,i)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine check_length_run_sr(iter,increase_nblk,nblk,nblk_max,denergy,denergy_err,energy_err_sav,energy_tol)

      implicit real*8(a-h,o-z)

c Increase nblk if near convergence to value needed to get desired statistical error
      increase_nblk=increase_nblk+1

c Increase if subsequent energies are within errorbar
      if(dabs(denergy).lt.3*denergy_err.and.energy_tol.gt.0.d0) then
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
        nbkl=1.2*nblk
        nblk=min(nblk,nblk_max)
        write(6,'(''nblk reset to'',i8,9d12.4)') nblk
      endif

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write_geometry(iter)

      implicit real*8(a-h,o-z)
      include 'vmc.h'

      character*40 filename,itn

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent

      if(iter.lt.0) then
        filename='geo_optimal_final'
       else
        if(iter.lt.10) then
          write(itn,'(i1)') iter
         elseif(iter.lt.100) then
          write(itn,'(i2)') iter
         elseif(iter.lt.1000) then
          write(itn,'(i3)') iter
        endif
        filename='geo_optimal.iter'//itn(1:index(itn,' ')-1)
      endif

      open(2,file=filename,status='unknown')

      write(2,'(''# geometry iter '',i5)') iter
      write(2,'(''&atoms nctype '',i5,'' natom '',i5)') nctype,ncent
      write(2,'(''geometry'')')
      do 20 i=1,ncent
  20    write(2,'(3f14.6,i5)') (cent(k,i),k=1,3),iwctype(i)
      write(2,'(''end'')')

      close(2)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      subroutine compute_positions
!      implicit real*8(a-h,o-z)
!
!      include 'vmc.h'
!      include 'force.h'
!
!      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent,iwctype(MCENT),nctype,ncent
!
!      common /force_analy/ iforce_analy,iuse_zmat,alfgeo
!      common /force_fin/ da_energy_ave(3,MCENT),da_energy_err(3)
!      common /zmatrix/ czcart(3,MCENT),czint(3,MCENT),
!     &                 czcart_ref(3,3),izcmat(3,MCENT),
!     &                 izmatrix
!      common /grdnthes/ hessian_zmat(3,MCENT)
!
!      dimension cent_ref(3,MCENT)
!
!      if(iforce_analy.eq.0)return
!
!      call compute_position_bcast
!
!      if(iuse_zmat.eq.0) then
!
!        do ic=1,ncent
!          do k=1,3
!            cent(k,ic)=cent(k,ic)-alfgeo*da_energy_ave(k,ic)
!          enddo
!          write(6,*)'CENT ',(cent(k,ic),k=1,3)
!        enddo
!
!      else
!
!        do 10 ic=1,3
!          do 10 k=1,3
!   10       cent_ref(k,ic)=cent(k,ic)
!
!        call cart2zmat(ncent,cent,izcmat,czint)
!
!        do ic=1,ncent
!          do k=1,3
!            czint(k,ic)=czint(k,ic)-alfgeo*da_energy_ave(k,ic)/hessian_zmat(k,ic)
!          enddo
!          write(6,*)'CENT ',(czint(k,ic),k=1,3)
!        enddo
!
!        call zmat2cart_rc(ncent,izcmat,czint,cent,cent_ref)
!
!      endif
!
!      return
!      end


      subroutine compute_positions
        use coords_int
        implicit real*8(a-h,o-z)
      
        include 'vmc.h'
        include 'force.h'
        
        common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent,iwctype(MCENT),nctype,ncent
        
        common /force_analy/ iforce_analy,iuse_zmat,alfgeo
        common /force_fin/ da_energy_ave(3,MCENT),da_energy_err(3)
        common /zmatrix/ czcart(3,MCENT),czint(3,MCENT),
     &                   czcart_ref(3,3),izcmat(3,MCENT),
     &                   izmatrix
        
        if (iforce_analy.eq.0) return !TODO why is this here?
        
        call compute_position_bcast
        
        if(iuse_zmat.eq.0) then
      
          call init (MCENT)
          call transform_gradients (da_energy_ave)
          call compute_step_int (alfgeo)
          call do_step (czint, cent, izcmat)
        endif
        
        return
      end


