module sr_more
      interface !LAPACK interface
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        DOUBLE PRECISION FUNCTION ddot(N,DX,INCX,DY,INCY)
!*  -- Reference BLAS level1 routine --
          INTEGER incx,incy,n
          DOUBLE PRECISION dx(*),dy(*)
        END FUNCTION
        SUBROUTINE dgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!*  -- Reference BLAS level2 routine --
          DOUBLE PRECISION ALPHA,BETA
          INTEGER INCX,INCY,LDA,M,N
          CHARACTER TRANS
          DOUBLE PRECISION A(LDA,*),X(*),Y(*)
        END SUBROUTINE
        SUBROUTINE daxpy(N,DA,DX,INCX,DY,INCY)
! -- Reference BLAS level1 routine --
          DOUBLE PRECISION DA
          INTEGER INCX,INCY,N
          DOUBLE PRECISION DX(*),DY(*)
        END SUBROUTINE
        SUBROUTINE dscal(N,DA,DX,INCX)
! -- Reference BLAS level1 routine --
          DOUBLE PRECISION DA
          INTEGER INCX,N
          DOUBLE PRECISION DX(*)
        END SUBROUTINE
        end interface
contains
      subroutine pcg(n,b,x,i,imax,imod,eps)
! one-shot preconditioned conjugate gradients; convergence thr is residual.lt.initial_residual*eps**2 (after J.R.Shewchuck)

      use contrl_file, only: ounit
      use mpi
      use mpiconf, only: idtask
      use sr_mat_n, only: ortho
      implicit none

      integer, parameter :: m_parm_opt = 59000
      integer :: n, imax, imod, i, j
      real(8) b(*),x(*),eps
      real(8), save :: r(m_parm_opt), d(m_parm_opt), q(m_parm_opt), s(m_parm_opt)
      real(8) delta_0,delta_new,delta_old,alpha,beta

      if(n.gt.m_parm_opt) stop 'nparm > m_parm_opt'

      call atimes_n(n,x,r)         ! r=Ax neuscamman
      if(idtask.eq.0)then
       call daxpy(n,-1.d0,b,1,r,1)       ! r=r-b
       call dscal(n,-1.d0,r,1)           ! r=b-r
       call asolve(n,r,d)                ! d=M^{-1}r preconditioner
       delta_new=ddot(n,d,1,r,1)           ! \delta_new=r^T d
       write(ounit,'(a12,f24.16)') 'delta0 = ',delta_new
      endif
      call MPI_BCAST(delta_new,1,MPI_REAL8,0,MPI_COMM_WORLD,j)
      delta_0=delta_new*eps**2            ! convergence thr
      do i=0,imax-1
!      write(*,*)i,idtask,'ECCO ',delta_0,delta_new
       if(delta_new.lt.delta_0)then
        if(idtask.eq.0) write(ounit,*) 'CG iter ',i
!     write(*,*)'ECCO pcg esce ',idtask
        call MPI_BCAST(x,n,MPI_REAL8,0,MPI_COMM_WORLD,j)
        return
       endif
       call atimes_n(n,d,q)         ! r=Ax neuscamman
       if(idtask.eq.0)then
        alpha=delta_new/ddot(n,d,1,q,1)  ! \alpha=\delta_new/(d^T q)
        call daxpy(n,alpha,d,1,x,1)      ! x=x+\alpha d
       endif
       if(mod(i,imod).eq.0)then
        call atimes_n(n,x,r)         ! r=Ax neuscamman
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
        write(ounit,'(a12,f24.16)') 'delta_new ',delta_new
        beta=delta_new/delta_old         ! \beta=\delta_new/\delta_old
        call dscal(n,beta,d,1)           ! d=\beta d
        call daxpy(n,1.d0,s,1,d,1)       ! d=s+d
       endif
       call MPI_BCAST(delta_new,1,MPI_REAL8,0,MPI_COMM_WORLD,j)
      enddo

      if(idtask.eq.0) write(ounit,*) 'CG iter ',i
      call MPI_BCAST(x,n,MPI_REAL8,0,MPI_COMM_WORLD,j)
      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine asolve(n,b,x)
! x(i)=b(i)/s(i,i) (preconditioning with diag(S))

      use sr_mat_n, only: s_ii_inv, sr_state
      use precision_kinds, only: dp
      use contrl_file, only: ounit

      implicit none

      integer :: i, n
      real(dp), dimension(*) :: x
      real(dp), dimension(*) :: b

      do i=1,n
       x(i)=b(i)*s_ii_inv(i, sr_state)
      enddo

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine atimes_n(n,z,r)
! r=a*z, i cicli doppi su n e nconf_n sono parallelizzati

      use sr_mod, only: mparm, mconf
      use csfs, only: nstates
      use optwf_func, only: ifunc_omega, omega, omega_hes
      use sa_weights, only: weights
      use sr_index, only: jelo, jelo2, jelohfj
      use sr_mat_n, only: jefj, jfj, jhfj, nconf_n, s_diag, sr_ho
      use sr_mat_n, only: sr_o, wtg, obs_tot, ortho, sr_state
      use optorb_cblock, only: norbterm
      use mpiconf, only: idtask
      use mpi
      use mpiconf, only: idtask
      use optorb_cblock, only: norbterm
      use optwf_func, only: ifunc_omega,omega,omega_hes
      use precision_kinds, only: dp
      use contrl_file, only: ounit
      use sr_mat_n, only: elocal, h_sr 

      implicit none

      integer :: i, i0, i1, iconf, istate
      integer :: k, n, nmparm_jasci, nparm_jasci
      real(dp) :: aux0, aux2, aux3, aux4
      real(dp) :: hoz, oz, oz_orb, var
      real(dp) :: wts
      real(dp), dimension(*) :: z
      real(dp), dimension(*) :: r
      real(dp), dimension(0:mconf) :: aux
      real(dp), dimension(0:mconf) :: aux1
      real(dp), dimension(mparm) :: rloc
      real(dp), dimension(mparm) :: r_s
      real(dp), dimension(mconf) :: oz_jasci
      real(dp), dimension(mparm) :: tmp
      real(dp), dimension(mparm) :: tmp2


      call MPI_BCAST(z,n,MPI_REAL8,0,MPI_COMM_WORLD,i)

      nparm_jasci=max(n-norbterm,0)

      if (ortho.eq.0) then

        do i=1,n
          r(i)=0.d0
        enddo

        if(ifunc_omega.eq.0) then

          do iconf=1,nconf_n
            oz_jasci(iconf)=ddot(nparm_jasci,z,1,sr_o(1,iconf,1),1)
          enddo

          do istate=1,nstates
            wts=weights(istate)

            i0=nparm_jasci+(istate-1)*norbterm+1
            do iconf=1,nconf_n
              oz_orb=ddot(norbterm,z(nparm_jasci+1),1,sr_o(i0,iconf,1),1)
              aux(iconf)=(oz_jasci(iconf)+oz_orb)*wtg(iconf,istate)
            enddo

!       Following three lines commented and replaced by dgemv after profiling
!        do i=1,nparm_jasci
!          rloc(i)=ddot(nconf_n,aux(1),1,sr_o(i,1),mparm)
!        enddo

            call dgemv('N', nparm_jasci, nconf_n, 1.0d0, sr_o(1,1,1), mparm, aux(1), 1, 0.0d0, rloc(1), 1)


!       Following code commented and replaced by dgemv after profiling
!        do i=nparm_jasci+1,n
!          i0=i+(istate-1)*norbterm
!          rloc(i)=ddot(nconf_n,aux(1),1,sr_o(i0,1),mparm)
!        enddo

            i0 = nparm_jasci + 1 +(istate-1)*norbterm
            i1 = nparm_jasci + 1
            call dgemv('N', n - nparm_jasci, nconf_n, 1.0d0, sr_o(i0,1,1), mparm, aux(1), 1, 0.0d0, rloc(i1), 1)

            call MPI_REDUCE(rloc,r_s,n,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,i)

            if(idtask.eq.0)then
              aux0=ddot(n,z,1,obs_tot(jfj,istate),1)
              do i=1,n
                r(i)=r(i)+wts*(r_s(i)/obs_tot(1,istate)-obs_tot(jfj+i-1,istate)*aux0+s_diag(i,istate)*z(i))
              enddo
            endif

          enddo

        else
! ifunc_omega.gt.0

          if(ifunc_omega.eq.1.or.ifunc_omega.eq.2) omega_hes=omega

          do iconf=1,nconf_n
            hoz=ddot(n,z,1,sr_ho(1,iconf),1)
            oz=ddot(n,z,1,sr_o(1,iconf,1),1)
            aux(iconf)=(hoz-omega_hes*oz)*wtg(iconf,1)
          enddo
          do i=1,n
            rloc(i)=ddot(nconf_n,aux(1),1,sr_ho(i,1),mparm)
            rloc(i)=rloc(i)-omega_hes*ddot(nconf_n,aux(1),1,sr_o(i,1,1),mparm)
          enddo
          call MPI_REDUCE(rloc,r,n,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,i)

          if(idtask.eq.0) then

            do i=1,n
              r(i)=r(i)/obs_tot(1,1)+s_diag(1,1)*z(i)
            enddo

            var=omega_hes*omega_hes+obs_tot(jelo2,1)-2*omega_hes*obs_tot(jelo,1)

            do k=1,n
              tmp(k)=obs_tot(jelohfj+k-1,1)-omega_hes*obs_tot(jefj+k-1,1) &
                     -omega_hes*(obs_tot(jhfj+k-1,1)-omega_hes*obs_tot(jfj+k-1,1))
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
! endif idtask.eq.0

        endif
! endif omega.eq.1.or.2

      else ! (ortho.ne.0)

        do iconf=1,nconf_n
          aux(iconf)=ddot(n,z,1,sr_o(1,iconf,sr_state),1)*wtg(iconf,sr_state)
        enddo

        do i=1,n
          rloc(i)=ddot(nconf_n,aux(1),1,sr_o(i,1,sr_state),mparm)
        enddo

        call MPI_REDUCE(rloc,r,n,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,i)

        if(idtask.eq.0)then
          aux0=ddot(n,z,1,obs_tot(jfj,sr_state),1)
          do i=1,n
            r(i)=r(i)/obs_tot(1,sr_state) &
               -obs_tot(jfj+i-1,sr_state)*aux0+s_diag(i,sr_state)*z(i)
          enddo
        endif

      endif !end ortho if

      return
      end

end module
