ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sr_hs(nparm,sr_adiag)
c <elo>, <o_i>, <elo o_i>, <o_i o_i>; s_diag, s_ii_inv, h_sr

      use force, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use csfs, only: nstates
      use mstates_mod, only: MSTATES
      use mpiconf, only: idtask
      use optwf_func, only: ifunc_omega, omega
      use sa_weights, only: weights
      use sr_index, only: jelo, jelo2, jelohfj
      use sr_mat_n, only: elocal, h_sr, jefj, jfj, jhfj, nconf_n, obs, s_diag, s_ii_inv, sr_ho
      use sr_mat_n, only: sr_o, wtg, obs_tot
      use optorb_cblock, only: norbterm

      use method_opt, only: method

      implicit real*8(a-h,o-z)


      include 'mpif.h'
      include 'sr.h'
      include 'optorb.h'





      dimension obs_wtg(MSTATES),obs_wtg_tot(MSTATES)

      call p2gtid('optgeo:izvzb',izvzb,0,1)
      call p2gtid('optwf:sr_rescale',i_sr_rescale,0,1)

      nstates_eff=nstates
      if(method.eq.'lin_d') nstates_eff=1

      jwtg=1
      jelo=2
      n_obs=2
      jfj=n_obs+1
      n_obs=n_obs+nparm
      jefj=n_obs+1
      n_obs=n_obs+nparm
      jfifj=n_obs+1
      n_obs=n_obs+nparm


      jhfj=n_obs+1
      n_obs=n_obs+nparm
      jfhfj=n_obs+1
      n_obs=n_obs+nparm

c for omega functional
      jelo2=n_obs+1
      n_obs=n_obs+1
      jelohfj=n_obs+1
      n_obs=n_obs+nparm

      if(n_obs.gt.MOBS) call fatal_error('SR_HS LIN: n_obs > MOBS)')

      do k=1,nparm
        h_sr(k)=0.d0
        s_ii_inv(k)=0.d0
      enddo

      nparm_jasci=max(nparm-norbterm,0)

      do istate=1,nstates
        obs(jwtg,istate)=0.d0
        do iconf=1,nconf_n
          obs(jwtg,istate)=obs(jwtg,istate)+wtg(iconf,istate)
        enddo
        obs_wtg(istate)=obs(jwtg,istate)
      enddo

      call MPI_REDUCE(obs_wtg,obs_wtg_tot,nstates,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ier)
      do istate=1,nstates
        obs_tot(jwtg,istate)=obs_wtg_tot(istate)
      enddo

      do istate=1,nstates_eff
        do i=2,n_obs
         obs(i,istate)=0.d0
        enddo

        ish=(istate-1)*norbterm
        do iconf=1,nconf_n
          obs(jelo,istate)=obs(jelo,istate)+elocal(iconf,istate)*wtg(iconf,istate)
          do i=1,nparm_jasci
            obs(jfj +i-1,istate)=obs(jfj +i-1,istate)+sr_o(i,iconf)*wtg(iconf,istate)
            obs(jefj+i-1,istate)=obs(jefj+i-1,istate)+elocal(iconf,istate)*sr_o(i,iconf)*wtg(iconf,istate)
            obs(jfifj+i-1,istate)=obs(jfifj+i-1,istate)+sr_o(i,iconf)*sr_o(i,iconf)*wtg(iconf,istate)
          enddo
          do i=nparm_jasci+1,nparm
            obs(jfj +i-1,istate)=obs(jfj +i-1,istate)+sr_o(ish+i,iconf)*wtg(iconf,istate)
            obs(jefj+i-1,istate)=obs(jefj+i-1,istate)+elocal(iconf,istate)*sr_o(ish+i,iconf)*wtg(iconf,istate)
            obs(jfifj+i-1,istate)=obs(jfifj+i-1,istate)+sr_o(ish+i,iconf)*sr_o(ish+i,iconf)*wtg(iconf,istate)
          enddo
        enddo

        call MPI_REDUCE(obs(1,istate),obs_tot(1,istate),n_obs,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ier)
      enddo

      if(idtask.eq.0) then
        do istate=1,nstates_eff
          wts=weights(istate)
          if(method.eq.'lin_d') wts=1.d0

          do i=2,n_obs
            obs_tot(i,istate)=obs_tot(i,istate)/obs_tot(1,istate)
          enddo

          do k=1,nparm
            aux=obs_tot(jfifj+k-1,istate)-obs_tot(jfj+k-1,istate)*obs_tot(jfj+k-1,istate)
            s_diag(k,istate)=aux*sr_adiag
            s_ii_inv(k)=s_ii_inv(k)+wts*(aux+s_diag(k,istate))
            h_sr(k)=h_sr(k)-2*wts*(obs_tot(jefj+k-1,istate)-obs_tot(jfj +k-1,istate)*obs_tot(jelo,istate))
          enddo
        enddo

        smax=0.d0
        do k=1,nparm
          if(s_ii_inv(k).gt.smax) smax=s_ii_inv(k)
        enddo 
        write(6,'(''max S diagonal element '',t41,d8.2)') smax

        kk=0
        do k=1,nparm
          if(s_ii_inv(k)/smax.gt.eps_eigval) then
            kk=kk+1
            s_ii_inv(k)=1.d0/s_ii_inv(k)
           else
            s_ii_inv(k)=0.d0
          endif
        enddo
        write(6,'(''nparm, non-zero S diag'',t41,2i5)') nparm,kk

      endif

      if(method.eq.'sr_n'.and.i_sr_rescale.eq.0.and.izvzb.eq.0.and.ifunc_omega.eq.0) return

      if(method.ne.'sr_n') then
        s_diag(1,1)=sr_adiag !!!

        do k=1,nparm
         h_sr(k)=-0.5d0*h_sr(k)
        enddo
      elseif(ifunc_omega.ne.0) then
        s_diag(1,1)=sr_adiag !!!
      endif

      if(n_obs.gt.MOBS) call fatal_error('SR_HS LIN: n_obs > MOBS)')

      do i=jhfj,n_obs
       obs(i,1)=0.d0
      enddo
      do iconf=1,nconf_n
       obs(jelo2,1)=obs(jelo2,1)+elocal(iconf,1)*elocal(iconf,1)*wtg(iconf,1)
       do i=1,nparm
         obs(jhfj+i-1,1)=obs(jhfj+i-1,1)+sr_ho(i,iconf)*wtg(iconf,1)
         obs(jfhfj+i-1,1)=obs(jfhfj+i-1,1)+sr_o(i,iconf)*sr_ho(i,iconf)*wtg(iconf,1)
         obs(jelohfj+i-1,1)=obs(jelohfj+i-1,1)+elocal(iconf,1)*sr_ho(i,iconf)*wtg(iconf,1)
       enddo
      enddo

      nreduce=n_obs-jhfj+1
      call MPI_REDUCE(obs(jhfj,1),obs_tot(jhfj,1),nreduce,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,j)

      if(idtask.eq.0)then
        do i=jhfj,n_obs
         obs_tot(i,1)=obs_tot(i,1)/obs_tot(1,1)
        enddo

        if(ifunc_omega.eq.1) then
c variance
          var=obs_tot(jelo2,1)-obs_tot(jelo,1)**2
          do k=1,nparm
           h_sr(k)=-2*(obs_tot(jelohfj+k-1,1)-(obs_tot(jhfj+k-1,1)-obs_tot(jefj+k-1,1))*obs_tot(jelo,1)
     &             -obs_tot(jfj+k-1,1)*obs_tot(jelo2,1) 
     &             -2*obs_tot(jelo,1)*(obs_tot(jefj+k-1,1)-obs_tot(jfj+k-1,1)*obs_tot(jelo,1)))
          enddo
        elseif(ifunc_omega.eq.2) then
c variance with fixed average energy (omega)
          var=omega*omega+obs_tot(jelo2,1)-2*omega*obs_tot(jelo,1)
          dum1=-2
          do k=1,nparm
           h_sr(k)=dum1*(omega*omega*obs_tot(jfj+k-1,1)+obs_tot(jelohfj+k-1,1)-omega*(obs_tot(jhfj+k-1,1)+obs_tot(jefj+k-1,1))
     &     -var*obs_tot(jfj+k-1,1))
c adding a term which intergrates to zero
c    &     -(obs_tot(jelo,1)-omega)*(obs_tot(jhfj+k-1,1)-obs_tot(jefj+k-1,1)))
          enddo

        elseif(ifunc_omega.eq.3.and.method.eq.'sr_n') then
c Neuscamman's functional
          den=omega*omega+obs_tot(jelo2,1)-2*omega*obs_tot(jelo,1)
          dum1=-2/den
          dum2=(omega-obs_tot(jelo,1))/den
          do k=1,nparm
           h_sr(k)=dum1*(omega*obs_tot(jfj+k-1,1)-obs_tot(jefj+k-1,1)
     &     -dum2*(omega*omega*obs_tot(jfj+k-1,1)+obs_tot(jelohfj+k-1,1)-omega*(obs_tot(jhfj+k-1,1)+obs_tot(jefj+k-1,1))))
          enddo
        endif

      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine pcg(n,b,x,i,imax,imod,eps)
c one-shot preconditioned conjugate gradients; convergence thr is residual.lt.initial_residual*eps**2 (after J.R.Shewchuck)

      use mpiconf, only: idtask
      implicit real*8(a-h,o-z)

      include 'mpif.h'
      integer m_parm_opt
      parameter(m_parm_opt=59000)
      integer n,imax,imod,i,j
      real*8 b(*),x(*),eps
      real*8 r(m_parm_opt),d(m_parm_opt),q(m_parm_opt),s(m_parm_opt)
      real*8 delta_0,delta_new,delta_old,alpha,beta,ddot

      if(n.gt.m_parm_opt) stop 'nparm > m_parm_opt'

      call atimes_n(n,x,r)         ! r=Ax neuscamman
      if(idtask.eq.0)then
       call daxpy(n,-1.d0,b,1,r,1)       ! r=r-b
       call dscal(n,-1.d0,r,1)           ! r=b-r
       call asolve(n,r,d)                ! d=M^{-1}r preconditioner
       delta_new=ddot(n,d,1,r,1)           ! \delta_new=r^T d
       print*,'delta0 = ',delta_new
      endif
      call MPI_BCAST(delta_new,1,MPI_REAL8,0,MPI_COMM_WORLD,j)
      delta_0=delta_new*eps**2            ! convergence thr
      do i=0,imax-1
c      write(*,*)i,idtask,'ECCO ',delta_0,delta_new 
       if(delta_new.lt.delta_0)then
        if(idtask.eq.0)print*,'CG iter ',i
c     write(*,*)'ECCO pcg esce ',idtask
        call MPI_BCAST(x,n,MPI_REAL8,0,MPI_COMM_WORLD,j)
        return
       endif
       call atimes_n(n,d,q)        ! q=Ad neuscamman
       if(idtask.eq.0)then
        alpha=delta_new/ddot(n,d,1,q,1)  ! \alpha=\delta_new/(d^T q)
        call daxpy(n,alpha,d,1,x,1)      ! x=x+\alpha d
       endif
       if(mod(i,imod).eq.0)then
        call atimes_n(n,x,r)       ! r=Ax neuscamman
        if(idtask.eq.0)then
         call daxpy(n,-1.d0,b,1,r,1)     ! r=r-b
         call dscal(n,-1.d0,r,1)         ! r=b-r
        endif
       else
        if(idtask.eq.0)call daxpy(n,-alpha,q,1,r,1)    ! r=r-\alpha q
       endif
       if(idtask.eq.0)then
        call asolve(n,r,s)               ! s=M^{-1}r preconditioner
        delta_old=delta_new              ! \delta_old=\delta_new
        delta_new=ddot(n,r,1,s,1)        ! \delta_new=r^T s
        print*,'delta_new ',delta_new
        beta=delta_new/delta_old         ! \beta=\delta_new/\delta_old
        call dscal(n,beta,d,1)           ! d=\beta d
        call daxpy(n,1.d0,s,1,d,1)       ! d=s+d
       endif
       call MPI_BCAST(delta_new,1,MPI_REAL8,0,MPI_COMM_WORLD,j)
      enddo

      if(idtask.eq.0)print*,'CG iter ',i
      call MPI_BCAST(x,n,MPI_REAL8,0,MPI_COMM_WORLD,j)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine asolve(n,b,x)
c x(i)=b(i)/s(i,i) (preconditioning with diag(S))

      use sr_mat_n, only: s_ii_inv
      implicit real*8(a-h,o-z)

      include 'sr.h'


      dimension x(*),b(*)

      do i=1,n
       x(i)=b(i)*s_ii_inv(i)
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine atimes_n(n,z,r)
c r=a*z, i cicli doppi su n e nconf_n sono parallelizzati

      use force, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use csfs, only: nstates

      use optwf_func, only: ifunc_omega, omega, omega_hes
      use sa_weights, only: weights
      use sr_index, only: jelo, jelo2, jelohfj
      use sr_mat_n, only: jefj, jfj, jhfj, nconf_n, s_diag, sr_ho
      use sr_mat_n, only: sr_o, wtg, obs_tot
      use optorb_cblock, only: norbterm
      implicit real*8(a-h,o-z)

      include 'mpif.h'
      include 'optorb.h'
      include 'sr.h'




      dimension z(*),r(*),aux(0:MCONF),aux1(0:MCONF),rloc(MPARM),r_s(MPARM),oz_jasci(MCONF)
      dimension tmp(MPARM),tmp2(MPARM)

      call MPI_BCAST(z,n,MPI_REAL8,0,MPI_COMM_WORLD,i)

      nparm_jasci=max(n-norbterm,0)

      do i=1,n
        r(i)=0.d0
      enddo

      if(ifunc_omega.eq.0) then 

      do iconf=1,nconf_n
        oz_jasci(iconf)=ddot(nparm_jasci,z,1,sr_o(1,iconf),1)
      enddo

      do istate=1,nstates
        wts=weights(istate)

        i0=nparm_jasci+(istate-1)*norbterm+1
        do iconf=1,nconf_n
          oz_orb=ddot(norbterm,z(nparm_jasci+1),1,sr_o(i0,iconf),1)
          aux(iconf)=(oz_jasci(iconf)+oz_orb)*wtg(iconf,istate)
        enddo

        do i=1,nparm_jasci
          rloc(i)=ddot(nconf_n,aux(1),1,sr_o(i,1),MPARM)
        enddo
        do i=nparm_jasci+1,n
          i0=i+(istate-1)*norbterm
          rloc(i)=ddot(nconf_n,aux(1),1,sr_o(i0,1),MPARM)
        enddo
        call MPI_REDUCE(rloc,r_s,n,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,i)
        
        if(idtask.eq.0)then
         aux0=ddot(n,z,1,obs_tot(jfj,istate),1)
         do i=1,n
          r(i)=r(i)+wts*(r_s(i)/obs_tot(1,istate)-obs_tot(jfj+i-1,istate)*aux0+s_diag(i,istate)*z(i))
         enddo
        endif

      enddo

      else
c ifunc_omega.gt.0

      if(ifunc_omega.eq.1.or.ifunc_omega.eq.2) omega_hes=omega

      do iconf=1,nconf_n
        hoz=ddot(n,z,1,sr_ho(1,iconf),1)
        oz=ddot(n,z,1,sr_o(1,iconf),1)
        aux(iconf)=(hoz-omega_hes*oz)*wtg(iconf,1)
      enddo
      do i=1,n
        rloc(i)=ddot(nconf_n,aux(1),1,sr_ho(i,1),MPARM)
        rloc(i)=rloc(i)-omega_hes*ddot(nconf_n,aux(1),1,sr_o(i,1),MPARM)
      enddo
      call MPI_REDUCE(rloc,r,n,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,i)

      if(idtask.eq.0) then

        do i=1,n
          r(i)=r(i)/obs_tot(1,1)+s_diag(1,1)*z(i)
        enddo

        var=omega_hes*omega_hes+obs_tot(jelo2,1)-2*omega_hes*obs_tot(jelo,1)

        do k=1,n
          tmp(k)=obs_tot(jelohfj+k-1,1)-omega_hes*obs_tot(jefj+k-1,1)-omega_hes*(obs_tot(jhfj+k-1,1)-omega_hes*obs_tot(jfj+k-1,1))
        enddo

        aux0=ddot(n,z,1,tmp(1),1)
        aux2=ddot(n,z,1,obs_tot(jfj,1),1)
        do i=1,n
          r(i)=r(i)-tmp(i)*aux2-obs_tot(jfj+i-1,1)*aux0+var*obs_tot(jfj+i-1,1)*aux2
        enddo

        do k=1,n
          tmp(k)=obs_tot(jhfj+k-1,1)-obs_tot(jelo,1)*obs_tot(jfj+k-1,1)
        enddo
        aux3=ddot(n,z,1,tmp(1),1)
        do i=1,n
          r(i)=r(i)-tmp(i)*aux3
        enddo

        do k=1,n
          tmp2(k)=obs_tot(jefj+k-1,1)-obs_tot(jelo,1)*obs_tot(jfj+k-1,1)
        enddo
        aux4=ddot(n,z,1,tmp2(1),1)
        do i=1,n
          r(i)=r(i)-tmp(i)*aux4-tmp2(i)*aux3
        enddo

      endif
c endif idtask.eq.0

      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sr_rescale_deltap(nparm,deltap)

      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use mpiconf, only: idtask
      use sr_mat_n, only: jefj, jfj, jhfj
      use sr_mat_n, only: obs_tot
    
      implicit real*8(a-h,o-z)



      include 'mpif.h'
      include 'sr.h'


      dimension deltap(*)

      call p2gtid('optwf:sr_rescale',i_sr_rescale,0,1)
      if(i_sr_rescale.eq.0) return

      jwtg=1
      jelo=2
      n_obs=2
      jfj=n_obs+1
      n_obs=n_obs+nparm
      jefj=n_obs+1
      n_obs=n_obs+nparm
      jfifj=n_obs+1
      n_obs=n_obs+nparm

      jhfj=n_obs+1
      n_obs=n_obs+nparm
      jfhfj=n_obs+1
      n_obs=n_obs+nparm

      jelo2=n_obs+1
      n_obs=n_obs+1
      jelohfj=n_obs+1
      n_obs=n_obs+nparm

      if(idtask.eq.0) then
        do i=1,nparm
          write(6,*) 'CIAO',obs_tot(jfhfj+i-1,1)/obs_tot(jfifj+i-1,1),obs_tot(jelo,1),
     &    obs_tot(jfhfj+i-1,1)/obs_tot(jfifj+i-1,1)-obs_tot(jelo,1)
          deltap(i)=deltap(i)/(obs_tot(jfhfj+i-1,1)/obs_tot(jfifj+i-1,1)-obs_tot(jelo,1))
        enddo
      endif

      call MPI_BCAST(deltap,nparm,MPI_REAL8,0,MPI_COMM_WORLD,j)

      return 
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine compute_position_bcast

      use force, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use atom, only: ncent
      use force_fin, only: da_energy_ave
      use force_analy, only: iforce_analy

      implicit real*8(a-h,o-z)

      include 'mpif.h'


      if(iforce_analy.eq.0)return

      call MPI_BCAST(da_energy_ave,3*ncent,MPI_REAL8,0,MPI_COMM_WORLD,i)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine forces_zvzb(nparm)

      use force, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use atom, only: ncent

      use force_fin, only: da_energy_ave
      use force_mat_n, only: force_o
      use mpiconf, only: idtask
      use sr_mat_n, only: elocal, jefj, jfj, jhfj, nconf_n, obs, sr_ho
      use sr_mat_n, only: sr_o, wtg

      implicit real*8(a-h,o-z)


      include 'mpif.h'
      include 'sr.h'

      parameter (MTEST=1500)
      dimension cloc(MTEST,MTEST),c(MTEST,MTEST),oloc(MPARM),o(MPARM),p(MPARM),tmp(MPARM)
      dimension ipvt(MTEST),work(MTEST)

      if(nparm.gt.MTEST) stop 'mparm>MTEST'

      jwtg=1
      jelo=2
      n_obs=2
      jfj=n_obs+1
      n_obs=n_obs+nparm
      jefj=n_obs+1
      n_obs=n_obs+nparm
      jfifj=n_obs+1
      n_obs=n_obs+nparm

      jhfj=n_obs+1
      n_obs=n_obs+nparm
      jfhfj=n_obs+1
      n_obs=n_obs+nparm

      do 10 i=1,nparm
        do 10 j=i,nparm
  10      cloc(i,j)=0.d0

      do l=1,nconf_n
        do i=1,nparm
          tmp(i)=(sr_ho(i,l)-elocal(l,1)*sr_o(i,l))*sqrt(wtg(l,1))
        enddo
        do k=1,nparm
          do j=k,nparm
            cloc(k,j)=cloc(k,j)+tmp(k)*tmp(j)
          enddo
        enddo
      enddo

      call MPI_REDUCE(cloc,c,MTEST*nparm,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,i)

      if(idtask.eq.0) then

        wtoti=1.d0/obs(1,1)
        do 20 i=1,nparm
          dum=(obs(jhfj+i-1,1)-obs(jefj+i-1,1))
          c(i,i)=c(i,i)*wtoti-dum*dum
          do 20 j=i+1,nparm
            c(i,j)=c(i,j)*wtoti-dum*(obs(jhfj+j-1,1)-obs(jefj+j-1,1))
  20        c(j,i)=c(i,j)

        call dgetrf(nparm,nparm,c,MTEST,ipvt,info)
        if(info.gt.0) then
          write(6,'(''MATINV: u(k,k)=0 with k= '',i5)') info
          call fatal_error('MATINV: info ne 0 in dgetrf')
        endif
        call dgetri(nparm,c,MTEST,ipvt,work,MTEST,info)

c ZVZB
c       do 30 iparm=1,nparm
c 30      tmp(iparm)=obs(jhfj+iparm-1,1)+obs(jefj+iparm-1,1)-2*obs(2,1)*obs(jfj+iparm-1,1)

c ZV
        do 30 iparm=1,nparm
  30      tmp(iparm)=obs(jhfj+iparm-1,1)-obs(jefj+iparm-1,1)
     
      endif

      energy_tot=obs(2,1)

      call MPI_BCAST(energy_tot,1,MPI_REAL8,0,MPI_COMM_WORLD,j)

      ia=0
      ish=3*ncent
      do icent=1,ncent
        write(6,'(''FORCE before'',i4,3e15.7)') icent,(da_energy_ave(k,icent),k=1,3)
        do k=1,3
          ia=ia+1

c         test=0.d0
c         do l=1,nconf_n
c           test=test+(force_o(ia+ish,l)-2*obs(2,1)*force_o(ia,l))*wtg(l,1)*wtoti
c         enddo
c         write(6,*) 'TEST ',test
         
          do i=1,nparm
            oloc(i)=0.d0
            do l=1,nconf_n
              oloc(i)=oloc(i)+(sr_ho(i,l)-elocal(l,1)*sr_o(i,l))*
     &                        (force_o(ia+ish,l)-2*energy_tot*force_o(ia,l))*wtg(l,1)
            enddo
          enddo

          call MPI_REDUCE(oloc,o,nparm,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,i)

          if(idtask.eq.0) then
            do i=1,nparm
              o(i)=o(i)*wtoti-(obs(jhfj+i-1,1)-obs(jefj+i-1,1))*da_energy_ave(k,icent)
            enddo
            do iparm=1,nparm
              p(iparm)=0.d0
              do jparm=1,nparm
                p(iparm)=p(iparm)+c(iparm,jparm)*o(jparm)
              enddo
              p(iparm)=-0.5*p(iparm)
            enddo

            force_tmp=da_energy_ave(k,icent)
            do iparm=1,nparm
              force_tmp=force_tmp+p(iparm)*tmp(iparm)
            enddo
            da_energy_ave(k,icent)=force_tmp

          endif
        enddo
        write(6,'(''FORCE after '',i4,3e15.7)') icent,(da_energy_ave(k,icent),k=1,3)
      enddo
          
      return
      end

