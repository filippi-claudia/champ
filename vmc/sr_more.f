cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sr_hs(nparm,sr_adiag)
c <elo>, <o_i>, <elo o_i>, <o_i o_i>; s_diag, s_ii_inv, h_sr

      implicit real*8 (a-h,o-z)
      include 'sr.h'
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'
      parameter(eps_eigval=1.d-14)

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm_sav
      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES)
     &,s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates
      common /sa_weights/ weights(MSTATES),iweight(MSTATES),nweight

      common /optwf_func/ omega,ifunc_omega

      call p2gtid('optgeo:izvzb',izvzb,0,1)

      write(6,*) 'NPARM',nparm

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

      print *, 'n_obs = ', n_obs

      if(n_obs.gt.MOBS) call fatal_error('SR_HS LIN: n_obs > MOBS)')

      do k=1,nparm
        h_sr(k)=0.d0
        s_ii_inv(k)=0.d0
      enddo

      nparm_jasci=max(nparm-norbterm,0)

      do istate=1,nstates
        obs(jwtg,istate)=0.d0
      enddo

      do istate=1,nstates
        do iconf=1,nconf
          obs(jwtg,istate)=obs(jwtg,istate)+wtg(iconf,istate)
        enddo
      enddo

      do istate=1,nstates_eff
        wts=weights(istate)
        if(method.eq.'lin_d') wts=1.d0

        do i=2,n_obs
         obs(i,istate)=0.d0
        enddo

        ish=(istate-1)*norbterm
        do iconf=1,nconf
c         obs(jwtg,istate)=obs(jwtg,istate)+wtg(iconf,istate)
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

        do i=2,n_obs
          obs(i,istate)=obs(i,istate)/obs(1,istate)
        enddo
        do k=1,nparm
         aux=obs(jfifj+k-1,istate)-obs(jfj+k-1,istate)*obs(jfj+k-1,istate)
         s_diag(k,istate)=aux*sr_adiag
         s_ii_inv(k)=s_ii_inv(k)+wts*(aux+s_diag(k,istate))
         h_sr(k)=h_sr(k)-2*wts*(obs(jefj+k-1,istate)-obs(jfj +k-1,istate)*obs(jelo,istate))
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
        
      if(method.eq.'sr_n'.and.izvzb.eq.0.and.ifunc_omega.eq.0) return

      if(method.ne.'sr_n') then
        s_diag(1,1)=sr_adiag !!!
        do k=1,nparm
         h_sr(k)=-0.5d0*h_sr(k)
        enddo
      endif

      do i=jhfj,n_obs
       obs(i,1)=0.d0
      enddo
      do iconf=1,nconf
       obs(jelo2,1)=obs(jelo2,1)+elocal(iconf,1)*elocal(iconf,1)*wtg(iconf,1)
       do i=1,nparm
         obs(jhfj+i-1,1)=obs(jhfj+i-1,1)+sr_ho(i,iconf)*wtg(iconf,1)
         obs(jfhfj+i-1,1)=obs(jfhfj+i-1,1)+sr_o(i,iconf)*sr_ho(i,iconf)*wtg(iconf,1)
         obs(jelohfj+i-1,1)=obs(jelohfj+i-1,1)+elocal(iconf,1)*sr_ho(i,iconf)*wtg(iconf,1)
       enddo
      enddo
      do i=jhfj,n_obs
       obs(i,1)=obs(i,1)/obs(1,1)
      enddo

      do k=1,nparm
        obs(jfhfj+i-1,1)=obs(jfhfj+i-1,1) !!! +sr_adiag
      enddo

      if(ifunc_omega.ne.0) then
        den=omega*omega+obs(jelo2,1)-2*omega*obs(jelo,1)
        dum1=-2/den
        dum2=(omega-obs(jelo,1))/den
        do k=1,nparm
         h_sr(k)=dum1*(omega*obs(jfj+k-1,1)-obs(jefj+k-1,1)
     &   -dum2*(omega*omega*obs(jfj+k-1,1)+obs(jelohfj+k-1,1)-omega*(obs(jhfj+k-1,1)+obs(jefj+k-1,1))))
        enddo
      endif

c     do k=1,nparm
c      s_kk=obs(jfifj+k-1)-obs(jfj+k-1)*obs(jfj+k-1)
c      h_kk=obs(jfhfj+k-1)-obs(jfj+k-1)*(obs(jhfj+k-1)*(1-obs(jelo))+obs(jefj+k-1))
c      s_ii_inv(k)=1.d0/(h_kk-ene*s_kk)
c     enddo

c     write(6,*) 'end of SR_HS'
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine forces_zvzb(nparm)

      implicit real*8 (a-h,o-z)

      include 'sr.h'
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent

      common /force_fin/ da_energy_ave(3,MCENT),da_energy_err(3)

      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES)
     &,s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf
 
      common /force_mat_n/ force_o(6*MCENT,MCONF)

      parameter (MTEST=1500)
      dimension c(MTEST,MTEST),o(MPARM),p(MPARM),tmp(MPARM)
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
  10      c(i,j)=0.d0

      wtoti=1.d0/obs(1,1)
      do l=1,nconf
        do i=1,nparm
          tmp(i)=(sr_ho(i,l)-elocal(l,1)*sr_o(i,l))*sqrt(wtg(l,1)*wtoti)
        enddo
        do k=1,nparm
          do j=k,nparm
            c(k,j)=c(k,j)+tmp(k)*tmp(j)
          enddo
        enddo
      enddo

      do 20 i=1,nparm
        dum=(obs(jhfj+i-1,1)-obs(jefj+i-1,1))
        c(i,i)=c(i,i)-dum*dum
        do 20 j=i+1,nparm
          c(i,j)=c(i,j)-dum*(obs(jhfj+j-1,1)-obs(jefj+j-1,1))
  20      c(j,i)=c(i,j)

c     do k=1,nparm
c       write(6,*) 'C ',(c(k,j),j=1,nparm)
c     enddo

      call dgetrf(nparm,nparm,c,MTEST,ipvt,info)
      if(info.gt.0) then
        write(6,'(''MATINV: u(k,k)=0 with k= '',i5)') info
        call fatal_error('MATINV: info ne 0 in dgetrf')
      endif
      call dgetri(nparm,c,MTEST,ipvt,work,MTEST,info)

c ZVZB
c     do 30 iparm=1,nparm
c 30    tmp(iparm)=obs(jhfj+iparm-1,1)+obs(jefj+iparm-1,1)-2*obs(2,1)*obs(jfj+iparm-1,1)

c ZV
      do 30 iparm=1,nparm
  30    tmp(iparm)=obs(jhfj+iparm-1,1)-obs(jefj+iparm-1,1)

      ia=0
      ish=3*ncent
      do icent=1,ncent
        write(6,'(''FORCE before'',i4,3e15.7)') icent,(da_energy_ave(k,icent),k=1,3)
        do k=1,3
          ia=ia+1

c         test=0.d0
c         do l=1,nconf
c           test=test+(force_o(ia+ish,l)-2*obs(2,1)*force_o(ia,l))*wtg(l,1)*wtoti
c         enddo
c         write(6,*) 'TEST ',test
         
          do i=1,nparm
            o(i)=0.d0
            do l=1,nconf
              o(i)=o(i)+(sr_ho(i,l)-elocal(l,1)*sr_o(i,l))*
     &                  (force_o(ia+ish,l)-2*obs(2,1)*force_o(ia,l))*wtg(l,1)*wtoti
            enddo
            o(i)=o(i)-(obs(jhfj+i-1,1)-obs(jefj+i-1,1))*da_energy_ave(k,icent)
            
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

        enddo
        write(6,'(''FORCE after '',i4,3e15.7)') icent,(da_energy_ave(k,icent),k=1,3)
      enddo
          
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine pcg(n,b,x,i,imax,imod,eps)
c one-shot preconditioned conjugate gradients; convergence thr is residual.lt.initial_residual*eps**2 (after J.R.Shewchuck)

      implicit none
      integer m_parm_opt
      parameter(m_parm_opt=59000)
      integer n,imax,imod,i,j
      real*8 b(*),x(*),eps
      real*8 r(m_parm_opt),d(m_parm_opt),q(m_parm_opt),s(m_parm_opt)
      real*8 delta_0,delta_new,delta_old,alpha,beta,ddot

      if(n.gt.m_parm_opt) stop 'nparm > m_parm_opt'

      call atimes_n(n,x,r)         ! r=Ax neuscamman
      call daxpy(n,-1.d0,b,1,r,1)       ! r=r-b
      call dscal(n,-1.d0,r,1)           ! r=b-r
      call asolve(n,r,d)                ! d=M^{-1}r preconditioner
      delta_new=ddot(n,d,1,r,1)           ! \delta_new=r^T d
      print*,'delta0 = ',delta_new
      delta_0=delta_new*eps**2            ! convergence thr
      do i=0,imax-1
       if(delta_new.lt.delta_0)return
       call atimes_n(n,d,q)        ! q=Ad neuscamman
       alpha=delta_new/ddot(n,d,1,q,1)  ! \alpha=\delta_new/(d^T q)
       call daxpy(n,alpha,d,1,x,1)      ! x=x+\alpha d
       if(mod(i,imod).eq.0)then
        call atimes_n(n,x,r)       ! r=Ax neuscamman
        call daxpy(n,-1.d0,b,1,r,1)     ! r=r-b
        call dscal(n,-1.d0,r,1)         ! r=b-r
       else
        call daxpy(n,-alpha,q,1,r,1)    ! r=r-\alpha q
       endif
       call asolve(n,r,s)               ! s=M^{-1}r preconditioner
       delta_old=delta_new              ! \delta_old=\delta_new
       delta_new=ddot(n,r,1,s,1)        ! \delta_new=r^T s
       print*,'delta_new ',delta_new
       beta=delta_new/delta_old         ! \beta=\delta_new/\delta_old
       call dscal(n,beta,d,1)           ! d=\beta d
       call daxpy(n,1.d0,s,1,d,1)       ! d=s+d
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine asolve(n,b,x)
c x(i)=b(i)/s(i,i) (preconditioning with diag(S))

      implicit real*8 (a-h,o-z)
      include 'sr.h'
      include 'mstates.h'
      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES)
     &,s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf

      dimension x(*),b(*)
      do i=1,n
       x(i)=b(i)*s_ii_inv(i)
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine atimes_n(n,z,r)
c r=a*z, i cicli doppi su n e nconf sono parallelizzati

      implicit real*8 (a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'
      include 'sr.h'

      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES)
     &,s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates
      common /sa_weights/ weights(MSTATES),iweight(MSTATES),nweight

      dimension z(*),r(*),aux(0:MCONF),oz_jasci(MCONF)

      nparm_jasci=max(n-norbterm,0)

      do i=1,n
        r(i)=0.d0
      enddo

      do iconf=1,nconf
        oz_jasci(iconf)=ddot(nparm_jasci,z,1,sr_o(1,iconf),1)
      enddo

      do istate=1,nstates
        wts=weights(istate)

        i0=nparm_jasci+(istate-1)*norbterm+1
        do iconf=1,nconf
          oz_orb=ddot(norbterm,z(nparm_jasci+1),1,sr_o(i0,iconf),1)
          aux(iconf)=(oz_jasci(iconf)+oz_orb)*wtg(iconf,istate)
        enddo
        aux0=ddot(n,z,1,obs(jfj,istate),1)
        do i=1,nparm_jasci
          r(i)=r(i)+wts*(ddot(nconf,aux(1),1,sr_o(i,1),MPARM)/obs(1,istate)
     &                  -obs(jfj+i-1,istate)*aux0+s_diag(i,istate)*z(i))
        enddo
        do i=nparm_jasci+1,n
          i0=i+(istate-1)*norbterm
          r(i)=r(i)+wts*(ddot(nconf,aux(1),1,sr_o(i0,1),MPARM)/obs(1,istate)
     &                  -obs(jfj+i-1,istate)*aux0+s_diag(i,istate)*z(i))
        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine compute_position_bcast
      implicit real*8 (a-h,o-z)

      return
      end
