ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sr_hs(nparm,sr_adiag)
c <elo>, <o_i>, <elo o_i>, <o_i o_i>; s_diag, s_ii_inv, h_sr

      use sr_mod, only: MOBS
      use csfs, only: nstates
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





      dimension obs_wtg(nstates),obs_wtg_tot(nstates)

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



      dimension x(*),b(*)

      do i=1,n
       x(i)=b(i)*s_ii_inv(i)
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine atimes_n(n,z,r)
c r=a*z, i cicli doppi su n e nconf_n sono parallelizzati

      use sr_mod, only: MPARM, MCONF
      use csfs, only: nstates

      use optwf_func, only: ifunc_omega, omega, omega_hes
      use sa_weights, only: weights
      use sr_index, only: jelo, jelo2, jelohfj
      use sr_mat_n, only: jefj, jfj, jhfj, nconf_n, s_diag, sr_ho
      use sr_mat_n, only: sr_o, wtg, obs_tot
      use optorb_cblock, only: norbterm

      ! as not in master ... 
      use mpiconf, only: idtask

      implicit real*8(a-h,o-z)

      include 'mpif.h'




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
