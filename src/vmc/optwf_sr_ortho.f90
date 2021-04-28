!------------------------------------------------------------------------------
!        Optimization routine using stochastic reconfiguration
!------------------------------------------------------------------------------
!> @author
!> Claudia Filippi, Modified by: Ramón L. Panadés-Barrueta
!
! DESCRIPTION:
!> Optimize multiple states using the penalty method introduced by Pathak et al.
!> J. Chem. Phys. 154, 034101 (2021)
!
! URL           : https://github.com/filippi-claudia/champ
!---------------------------------------------------------------------------

module optwf_sr_ortho_mod

  use precision_kinds, only: dp
  use optwf_contrl, only: ioptci, ioptjas, ioptorb
  use force_analy, only: iforce_analy
  use contrl, only: nblk_max
  use optwf_contrl, only: energy_tol, nopt_iter, micro_iter_sr, dparm_norm_min
  use optwf_contrl, only: sr_tau , sr_adiag, sr_eps 

  integer :: ioptjas_sav, ioptorb_sav, ioptci_sav, iforce_analy_sav
  integer :: i_sr_rescale, izvzb
  real(dp), dimension(:, :), allocatable :: deltap
  real(dp) :: sr_adiag_sav

  private
  public :: optwf_sr_ortho, atimes_n_ortho
  save

contains

  subroutine optwf_sr_ortho

    use precision_kinds, only: dp
    use sr_mod, only: MPARM
    use csfs, only: nstates
    use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
    use mstates_mod, only: MSTATES
    use optwf_corsam, only: energy, energy_err, force
    use optwf_func, only: ifunc_omega, omega0, n_omegaf, n_omegat, omega_hes
    use contrl, only: nblk
    use force_analy, only: alfgeo
    use optwf_contrl, only: nparm
    use method_opt, only: method
    use optwf_sr_mod, only: forces_zvzb

    implicit real*8(a-h, o-z)

    allocate(deltap(MPARM,MSTATES))

    call set_nparms_tot

    if(nparm.gt.MPARM) call fatal_error('SR_OPTWF: nparmtot gt MPARM')

    write(6,'(''Starting dparm_norm_min'',g12.4)') dparm_norm_min
    write(6,'(/,''SR adiag: '',f10.5)') sr_adiag
    write(6,'(''SR tau:   '',f10.5)') sr_tau
    write(6,'(''SR eps:   '',f10.5)') sr_eps

    call save_params()
    call save_nparms
    call write_geometry(0)

    do iter=1,nopt_iter
       write(6,'(/,''Optimization iteration'',i5,'' of'',i5)')iter,nopt_iter
       iforce_analy=0

       do miter=1,micro_iter_sr
          if(micro_iter_sr.gt.1) then
             write(6,'(/,''Micro iteration'',i5,'' of'',i5)')miter,micro_iter_sr
          end if

          if(miter.eq.micro_iter_sr) iforce_analy=iforce_analy_sav

          call qmc
          write(6,'(/,''Completed sampling'')')

6         continue

          call sr_ortho(nparm,deltap,sr_adiag,sr_eps,i)

          adiag=sr_adiag
	  iflagin=0
          do istate=1,nstates
             call dscal(nparm,-sr_tau,deltap(:,istate),1)
             call test_solution_parm(nparm,deltap(:,istate),&
		     dparm_norm,dparm_norm_min,adiag,iflag)
             write(6,'(''Norm of parm variation '',d12.5)') dparm_norm
             if(iflag.ne.0) iflagin=1
          end do

          if(iflagin.ne.0) then
             write(6,'(''Warning: dparm_norm>1'')')
             adiag=10*adiag
             write(6,'(''adiag increased to '',f10.5)') adiag
             sr_adiag=adiag
             go to 6
          else
             sr_adiag=sr_adiag_sav
          endif

	  call compute_norm_lin(nparm,deltap)
          call compute_parameters(deltap,iflag,1)
          call write_wf(1,iter)
          call save_wf

          if(iforce_analy.gt.0) then
             if(izvzb.gt.0) call forces_zvzb(nparm)
             call compute_positions
             call write_geometry(iter)
          endif
       enddo

       if(iter.ge.2) then
          denergy=energy(1)-energy_sav
          denergy_err=sqrt(energy_err(1)**2+energy_err_sav**2)
          nblk=nblk*1.2
          nblk=min(nblk,nblk_max)
       endif

       write(6,'(''nblk = '',i6)') nblk
       write(6,'(''alfgeo = '',f10.4)') alfgeo

       energy_sav=energy(1)
       energy_err_sav=energy_err(1)
       sigma_sav=sigma
    end do

    write(6,'(/,''Check last iteration'')')

    ioptjas=0
    ioptorb=0
    ioptci=0
    iforce_analy=0

    call set_nparms
    call qmc
    call write_wf(1,-1)
    call write_geometry(-1)

    deallocate(deltap)

  end subroutine optwf_sr_ortho

  subroutine sr_ortho(nparm,deltap,sr_adiag,sr_eps,i)
    use const, only: ipr
    use csfs, only: nstates
    use sr_mat_n, only: h_sr, istat_curr
    use mpiconf, only: idtask
    use optwf_sr_mod, only: sr_rescale_deltap

    implicit real*8(a-h,o-z)
    integer, intent(in) :: nparm
    real(dp), dimension(:,:), intent(out) :: deltap
    real(dp), intent(in) :: sr_adiag
    real(dp), intent(in) :: sr_eps
    integer, intent(inout) :: i

    call compute_gradient_sr_ortho(nparm,sr_adiag)

    imax=nparm                ! max n. iterations conjugate gradients
    imod=50                   ! inv. freq. of calc. r=b-Ax vs. r=r-\alpha q (see pcg)
    deltap=0.0d0              ! initial guess of solution

    if (idtask.eq.0.and.ipr.gt.4) then
       print *, "Gradient state 1"
       print *, h_sr(1:nparm,1)
       print *, "Gradient state 2"
       print *, h_sr(1:nparm,2)
    end if

    do j=1,nstates
       istat_curr = j
       write(6,*) 'Orthogonal optimization state ', istat_curr
       call pcg(nparm,h_sr(1:nparm,istat_curr),deltap(1:nparm,istat_curr),i,imax,imod,sr_eps)
       call sr_rescale_deltap(nparm,deltap(1:nparm,istat_curr))
    enddo

  end subroutine sr_ortho

  subroutine compute_norm_lin(nparm,deltap)
    use mpi
    use mpiconf, only: idtask
    use sr_mod, only: MOBS
    use mstates_mod, only: MSTATES
    use sr_mod, only: MPARM
    use csfs, only: nstates
    use sr_mat_n, only: elocal, h_sr, jefj, jfj, jhfj, nconf_n, obs, s_diag, s_ii_inv, sr_ho
    use sr_mat_n, only: sr_o, wtg, obs_tot
    use estcum, only: iblk
    use contrl, only: nstep
    use optorb_cblock, only: norbterm
    use config, only: anormo

    implicit real*8 (a-h,o-z)

    integer, intent(in) :: nparm
    real(dp), dimension(MPARM,MSTATES), intent(in) :: deltap

    real(dp), dimension(MOBS,MSTATES) :: obs_norm
    real(dp), dimension(MOBS,MSTATES) :: obs_norm_tot

    jwtg=1
    jwfj=1
    jsqfj=2
    n_obs=2

    obs_norm=0.0d0
    passes=dfloat(iblk*nstep)

    do istate=1,nstates
       do iconf=1,nconf_n
          tmp=0.0d0
          do i=1,nparm
             tmp=tmp+deltap(i,istate)*sr_o(i,iconf,istate)
          enddo
          obs_norm(jwfj,istate)=obs_norm(jwfj,istate)+tmp*wtg(iconf,istate)
          obs_norm(jsqfj,istate)=obs_norm(jsqfj,istate)+tmp*tmp*wtg(iconf,istate)
       enddo
       obs_norm(jwfj,istate)=2.0d0*obs_norm(jwfj,istate)
       call MPI_REDUCE(obs_norm(1,istate),obs_norm_tot(1,istate),&
       	       n_obs,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ier)
    enddo

    if(idtask.eq.0) then
       do istate=1,nstates
          anormo(istate)=obs_tot(jwtg,istate)+obs_norm_tot(jwfj,istate)&
		  +obs_norm_tot(jsqfj,istate)
       enddo
       anormo=anormo/passes
       print *, "ANORMO", anormo
    endif
    call MPI_BCAST(anormo(1),nstates,MPI_REAL8,0,MPI_COMM_WORLD,i)

  end subroutine

  subroutine compute_gradient_sr_ortho(nparm,sr_adiag)
    use mpi
    use sr_mod, only: MOBS
    use csfs, only: nstates
    use mstates_mod, only: MSTATES
    use mpiconf, only: idtask
    use optwf_func, only: ifunc_omega, omega
    use sr_index, only: jelo, jelo2, jelohfj
    use sr_mat_n, only: elocal, h_sr, jefj, jfj, jhfj, nconf_n, obs, s_diag, s_ii_inv, sr_ho
    use sr_mat_n, only: sr_o, wtg, obs_tot
    use optorb_cblock, only: norbterm
    use method_opt, only: method

    implicit real*8 (a-h,o-z)

    integer, intent(in) :: nparm
    real(dp), dimension(nconf_n) :: aux

    jwtg=1
    n_obs=nstates
    jelo=n_obs+1
    n_obs=n_obs+1
    jfj=n_obs+1
    n_obs=n_obs+nparm
    jefj=n_obs+1
    n_obs=n_obs+nparm
    jfifj=n_obs+1
    n_obs=n_obs+nparm

    if(nstates.gt.1) then
       jfjsi=n_obs+1
       n_obs=n_obs+nparm
    endif

    if(n_obs.gt.MOBS) call fatal_error('SR_HS LIN: n_obs_ortho BS)')

    do iconf=1,nconf_n
       aux(iconf) = 0.0d0
    enddo

    do istate=1,nstates
       do i=1,n_obs
          obs(i,istate) = 0.0d0
       enddo
       do iconf=1,nconf_n
          obs(jwtg,istate)=obs(jwtg,istate)+wtg(iconf,istate)

          do jstate=1,istate-1
             obs(jwtg+jstate,istate)=obs(jwtg+jstate,istate)&
		     +sr_o(nparm+2,iconf,jstate)*sr_o(nparm+2,iconf,istate)
          enddo

          obs(jelo,istate)=obs(jelo,istate)&
		  +elocal(iconf,istate)*wtg(iconf,istate)

          do i=1,nparm
             obs(jfj+i-1,istate)=obs(jfj+i-1,istate)&              
		     +sr_o(i,iconf,istate)*wtg(iconf,istate)
             obs(jefj+i-1,istate)=obs(jefj+i-1,istate)&
       		     +elocal(iconf,istate)*sr_o(i,iconf,istate)*wtg(iconf,istate)
             obs(jfifj+i-1,istate)=obs(jfifj+i-1,istate)&
       		     +sr_o(i,iconf,istate)*sr_o(i,iconf,istate)*wtg(iconf,istate)
          enddo
       enddo

       n_obs_reduce=n_obs
       if(nstates.gt.1) n_obs_reduce=n_obs-nparm

       call MPI_REDUCE(obs(1,istate),obs_tot(1,istate),&
       	       n_obs_reduce,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ier)

       call MPI_BCAST(obs_tot(1,istate),nstates,MPI_REAL8,0,MPI_COMM_WORLD,i)

       if(istate.gt.1) then
          do iconf=1,nconf_n
             aux(iconf)=aux(iconf)+sr_o(nparm+2,iconf,istate-1)&
       		     *obs_tot(jwtg+istate-1,istate)/obs_tot(jwtg,istate-1)
             do i=1,nparm
                obs(jfjsi+i-1,istate)=obs(jfjsi+i-1,istate)&
       			+sr_o(i,iconf,istate)*sr_o(nparm+2,iconf,istate)*aux(iconf)
             enddo
          enddo
          call MPI_REDUCE(obs(jfjsi,istate),obs_tot(jfjsi,istate),&
       		  nparm,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ier)
       endif
    enddo

    if(idtask.eq.0) then
       alambda=1.0d0
       aux1=0.0d0
       do istate=1,nstates
          do i=2,n_obs
             obs_tot(i,istate)=obs_tot(i,istate)/obs_tot(jwtg,istate)
          enddo

          do k=1,nparm
             aux2=obs_tot(jfifj+k-1,istate)&
       		     -obs_tot(jfj+k-1,istate)*obs_tot(jfj+k-1,istate)
             s_diag(k,istate)=aux2*sr_adiag
             s_ii_inv(k,istate)=aux2+s_diag(k,istate)
             h_sr(k,istate)=-2.0d0*(obs_tot(jefj+k-1,istate)&
       		     -obs_tot(jfj+k-1,istate)*obs_tot(jelo,istate))
          enddo

          if(istate.gt.1) then
             aux1=aux1+obs_tot(jwtg+istate-1,istate)*obs_tot(jwtg+istate-1,istate)&
       	             *obs_tot(jwtg,istate)/obs_tot(jwtg,istate-1)

             do k=1,nparm
                h_sr(k,istate)=h_sr(k,istate)&
       	        	-alambda*2.d0*(obs_tot(jfjsi+k-1,istate)-aux1*obs_tot(jfj+k-1,istate))
             enddo
          endif

          smax=0.0d0
          do k=1,nparm
             if(s_ii_inv(k,istate).gt.smax) smax=s_ii_inv(k,istate)
          enddo
          write(6,'(''max S diagonal element '',t41,d8.2)') smax

          kk=0
          do k=1,nparm
             if(s_ii_inv(k,istate)/smax.gt.eps_eigval) then
                kk=kk+1
                s_ii_inv(k,istate)=1.0d0/s_ii_inv(k,istate)
             else
                s_ii_inv(k,istate)=0.0d0
             endif
          enddo
          write(6,'(''nparm, non-zero S diag'',t41,2i5)') nparm,kk
       enddo
    endif

  end subroutine compute_gradient_sr_ortho

  subroutine atimes_n_ortho(n,z,r)
    use mpi
    use sr_mod, only: MOBS,MCONF,MPARM
    use csfs, only: nstates
    use mstates_mod, only: MSTATES
    use mpiconf, only: idtask
    use optwf_func, only: ifunc_omega, omega
    use sr_index, only: jelo, jelo2, jelohfj
    use sr_mat_n, only: elocal, h_sr, jefj, jfj, jhfj, nconf_n, obs, s_diag, s_ii_inv, sr_ho
    use sr_mat_n, only: sr_o, wtg, obs_tot, istat_curr
    use optorb_cblock, only: norbterm
    use method_opt, only: method

    implicit real*8 (a-h,o-z)

    dimension z(*),r(*),aux(0:MCONF),rloc(MPARM)

    call MPI_BCAST(z,n,MPI_REAL8,0,MPI_COMM_WORLD,i)

    do iconf=1,nconf_n
       aux(iconf)=ddot(n,z,1,sr_o(1,iconf,istat_curr),1)*wtg(iconf,istat_curr)
    enddo

    do i=1,n
       rloc(i)=ddot(nconf_n,aux(1),1,sr_o(i,1,istat_curr),MPARM)
    enddo

    call MPI_REDUCE(rloc,r,n,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,i)

    if(idtask.eq.0)then
       aux0=ddot(n,z,1,obs_tot(jfj,istat_curr),1)
       do i=1,n
          r(i)=r(i)/obs_tot(1,istat_curr)&
       		  -obs_tot(jfj+i-1,istat_curr)*aux0+s_diag(i,istat_curr)*z(i)
       enddo
    endif

  end subroutine atimes_n_ortho

  subroutine save_params()
      sr_adiag_sav = sr_adiag
      iforce_analy_sav = iforce_analy
      ioptjas_sav = ioptjas
      ioptorb_sav = ioptorb
      ioptci_sav = ioptci
  end subroutine save_params

end module optwf_sr_ortho_mod

