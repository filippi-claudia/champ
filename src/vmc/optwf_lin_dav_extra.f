      module optwf_lin_dav_extra
      use sr_more, only: ddot
      contains
      subroutine h_psi_energymin(ndim,nvec,psi,hpsi )

      use contrl_file, only: errunit,ounit
      use mpi
      use mpiconf, only: idtask
      use optwf_control, only: ioptjas,ioptorb,nparm,nvecx
      use precision_kinds, only: dp
      use sr_mat_n, only: h_sr,jefj,jfj,jhfj,nconf_n,obs_tot,s_diag
      use sr_mat_n, only: sr_ho,sr_o,wtg
      use sr_mod,  only: mconf,mparm
      ! these were not called in the master
      ! but they seem to be needed
      ! use sr_index, only: jelo, jelo2, jelohfj

      implicit none

      integer :: i, i0, iconf, ier, ivec
      integer :: jelo, jfifj, jwtg, n_obs
      integer :: ndim, nvec
      real(dp) :: aux0, aux1, aux2
      real(dp), dimension(mparm,*) :: psi
      real(dp), dimension(mparm,*) :: hpsi
      real(dp), dimension(mconf) :: aux
      real(dp), dimension(mparm,nvecx) :: hpsiloc

      i0=1
      if(ioptorb.eq.0.and.ioptjas.eq.0) i0=0
      nparm=ndim-i0

c     write(ounit,*) 'HPSI_LIN',ndim,nvec

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


c loop vec
      do ivec=1,nvec

      call MPI_BCAST(psi(1,ivec),ndim,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      do iconf=1,nconf_n
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_ho(1,iconf),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       hpsiloc(i+i0,ivec)=0.5d0*ddot(nconf_n,aux(1),1,sr_o(i,1,1),mparm)
      enddo

      do iconf=1,nconf_n
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_o(1,iconf,1),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       hpsiloc(i+i0,ivec)=hpsiloc(i+i0,ivec)+0.5d0*ddot(nconf_n,aux(1),1,sr_ho(i,1),mparm)
      enddo

      call MPI_REDUCE(hpsiloc(1+i0,ivec),hpsi(1+i0,ivec),nparm,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ier)

      if(idtask.eq.0) then

        do i=1,nparm
         hpsi(i+i0,ivec)=hpsi(i+i0,ivec)/obs_tot(1,1)
        enddo

        if(i0.eq.1)then

          aux0=ddot(nparm,psi(i0+1,ivec),1,obs_tot(jfj,1),1)
          aux1=ddot(nparm,psi(i0+1,ivec),1,obs_tot(jefj,1),1)
          aux2=ddot(nparm,psi(i0+1,ivec),1,obs_tot(jhfj,1),1)
          do i=1,nparm
           hpsi(i+i0,ivec)=hpsi(i+i0,ivec)+aux0*obs_tot(jelo,1)*obs_tot(jfj+i-1,1)
     &                            -0.5d0*(aux0*(obs_tot(jefj+i-1,1)+obs_tot(jhfj+i-1,1))
     &                                   +obs_tot(jfj+i-1,1)*(aux1+aux2))
          enddo

          hpsi(1,ivec)=obs_tot(jelo,1)*psi(1,ivec)+ddot(nparm,h_sr,1,psi(2,ivec),1)
          do i=1,nparm
           hpsi(i+i0,ivec)=hpsi(i+i0,ivec)+h_sr(i,1)*psi(1,ivec)+s_diag(1,1)*psi(i+i0,ivec) !!!
          enddo

        endif

      endif

c end loop vec
      enddo

c     do i=1,nparm+1
c       write(ounit,'(''HPSI_LIN'',100e12.3)')(hpsi(i,ivec),ivec=1,nvec)
c     enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine s_psi_energymin(ndim,nvec,psi,spsi )

      use contrl_file, only: errunit,ounit
      use mpi
      use mpiconf, only: idtask
      use optwf_control, only: ioptjas,ioptorb,nparm,nvecx
      use precision_kinds, only: dp
      use sr_mat_n, only: jefj,jfj,jhfj,nconf_n,obs_tot,sr_o,wtg
      use sr_mod,  only: mconf,mparm
      ! these were not called in the master
      ! but they seem to be needed
      ! use sr_index, only: jelo, jelo2, jelohfj

      implicit none

      integer :: i, i0, iconf, ier, ivec
      integer :: jelo, jfifj, jwtg, n_obs
      integer :: ndim, nvec
      real(dp) :: aux0
      real(dp), dimension(mparm,*) :: psi
      real(dp), dimension(mparm,*) :: spsi
      real(dp), dimension(mparm,nvecx) :: spsiloc
      real(dp), dimension(mconf) :: aux

      i0=1
      if(ioptorb.eq.0.and.ioptjas.eq.0) i0=0
      nparm=ndim-i0

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

c loop vec
      do ivec=1,nvec

      call MPI_BCAST(psi(1,ivec),ndim,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      do iconf=1,nconf_n
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_o(1,iconf,1),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       spsiloc(i+i0,ivec)=ddot(nconf_n,aux(1),1,sr_o(i,1,1),mparm)
      enddo

      call MPI_REDUCE(spsiloc(1+i0,ivec),spsi(1+i0,ivec),nparm,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ier)

      if(idtask.eq.0) then

        do i=1,nparm
         spsi(i+i0,ivec)=spsi(i+i0,ivec)/obs_tot(1,1)
        enddo

        if(i0.eq.1) then
          aux0=ddot(nparm,psi(1+i0,ivec),1,obs_tot(jfj,1),1)
          do i=1,nparm
            spsi(i+i0,ivec)=spsi(i+i0,ivec)-aux0*obs_tot(jfj+i-1,1)
          enddo

          spsi(1,ivec)=psi(1,ivec)
        endif
      endif

c end loop vec
      enddo

c     do i=1,nparm+1
c       write(ounit,'(''SPSI_LIN'',100e12.3)')(spsi(i,ivec),ivec=1,nvec)
c     enddo
c     STOP

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine h_psi_omegamin(ndim,nvec,psi,hpsi )

      use contrl_file, only: errunit,ounit
      use mpi
      use mpiconf, only: idtask
      use optwf_control, only: ioptjas,ioptorb,nparm,nvecx
      use optwf_func, only: omega
      use precision_kinds, only: dp
      use sr_mat_n, only: h_sr,jefj,jfj,jhfj,nconf_n,obs_tot,s_diag
      use sr_mat_n, only: sr_ho,sr_o,wtg
      use sr_mod,  only: mconf,mparm
      ! these were not called in the master
      ! but they seem to be needed
      ! use sr_index, only: jelo, jelo2, jelohfj

      implicit none

      integer :: i, i0, iconf, ier, ivec
      integer :: jelo, jfifj, jwtg, n_obs
      integer :: ndim, nvec
      real(dp) :: aux0, aux1, aux2
      real(dp), dimension(mparm,*) :: psi
      real(dp), dimension(mparm,*) :: hpsi
      real(dp), dimension(mparm,nvecx) :: hpsiloc
      real(dp), dimension(mconf) :: aux

      i0=1
      if(ioptorb.eq.0.and.ioptjas.eq.0) i0=0
      nparm=ndim-i0

c     do i=1,nparm+i0
c       write(ounit,'(''PSI NEW'',100e12.3)')(psi(i,ivec),ivec=1,nvec)
c     enddo

c     write(ounit,*) 'HPSI_LIN',ndim,nvec

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

c loop vec
      do ivec=1,nvec

      call MPI_BCAST(psi(1,ivec),ndim,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      do iconf=1,nconf_n
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_ho(1,iconf),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       hpsiloc(i+i0,ivec)=-0.5d0*ddot(nconf_n,aux(1),1,sr_o(i,1,1),mparm)
      enddo

c     do i=1,nparm+1
c       write(ounit,*) 'H1 ',hpsi(i,ivec)
c     enddo

      do iconf=1,nconf_n
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_o(1,iconf,1),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       hpsiloc(i+i0,ivec)=hpsiloc(i+i0,ivec)-0.5d0*ddot(nconf_n,aux(1),1,sr_ho(i,1),mparm)
     &                                      +omega*ddot(nconf_n,aux(1),1,sr_o(i,1,1),mparm)
      enddo

      call MPI_REDUCE(hpsiloc(1+i0,ivec),hpsi(1+i0,ivec),nparm,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ier)

c     do i=1,nparm+1
c       write(ounit,*) 'H2 ',hpsi(i,ivec)
c     enddo

      if(idtask.eq.0) then

        do i=1,nparm
         hpsi(i+i0,ivec)=hpsi(i+i0,ivec)/obs_tot(1,1)
        enddo
        if(i0.eq.1) then

          aux0=ddot(nparm,psi(i0+1,ivec),1,obs_tot(jfj,1),1)
          aux1=ddot(nparm,psi(i0+1,ivec),1,obs_tot(jefj,1),1)
          aux2=ddot(nparm,psi(i0+1,ivec),1,obs_tot(jhfj,1),1)
          do i=1,nparm
           hpsi(i+i0,ivec)=hpsi(i+i0,ivec)-aux0*obs_tot(jelo,1)*obs_tot(jfj+i-1,1)
     &                              +0.5d0*(aux0*(obs_tot(jefj+i-1,1)+obs_tot(jhfj+i-1,1))
     &                                     +obs_tot(jfj+i-1,1)*(aux1+aux2))
          enddo

          hpsi(1,ivec)=-(obs_tot(jelo,1)*psi(1,ivec)+ddot(nparm,h_sr,1,psi(2,ivec),1))
          do i=1,nparm
           hpsi(i+i0,ivec)=hpsi(i+i0,ivec)-h_sr(i,1)*psi(1,ivec)+s_diag(1,1)*psi(i+i0,ivec) !!!
          enddo

          do i=1,nparm
           hpsi(i+i0,ivec)=hpsi(i+i0,ivec)-omega*aux0*obs_tot(jfj+i-1,1)
          enddo
          hpsi(1,ivec)=hpsi(1,ivec)+omega*psi(1,ivec)

        endif
      endif

c end loop vec
      enddo

c     do i=1,nparm+i0
c       write(ounit,'(''HPSI_LIN'',100e12.3)')(hpsi(i,ivec),ivec=1,nvec)
c     enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine s_psi_omegamin(ndim,nvec,psi,spsi )

      use contrl_file, only: errunit,ounit
      use mpi
      use mpiconf, only: idtask
      use optwf_control, only: ioptjas,ioptorb,nparm,nvecx
      use optwf_func, only: omega
      use precision_kinds, only: dp
      use sr_mat_n, only: h_sr,jefj,jfj,jhfj,nconf_n,obs_tot,sr_ho,sr_o
      use sr_mat_n, only: wtg
      use sr_mod,  only: mconf,mparm
      ! these were not called in the master
      ! but they seem to be needed
      ! use sr_index, only: jelo, jelo2, jelohfj

      implicit none

      integer :: i, i0, iconf, ier, ivec
      integer :: jelo, jelo2, jelohfj, jfhfj
      integer :: jfifj, jwtg, n_obs, ndim
      integer :: nvec
      real(dp) :: aux0, aux1, aux2, aux3
      real(dp), dimension(mparm,*) :: psi
      real(dp), dimension(mparm,*) :: spsi
      real(dp), dimension(mparm,nvecx) :: spsiloc
      real(dp), dimension(mconf) :: aux
      real(dp), dimension(mparm) :: h_sr_sym

      i0=1
      if(ioptorb.eq.0.and.ioptjas.eq.0) i0=0
      nparm=ndim-i0

      jwtg=1
      jelo=2
      n_obs=2
      jfj=n_obs+1
      n_obs=n_obs+nparm
      jefj=n_obs+1
      n_obs=n_obs+nparm
      jfifj=n_obs+1
      n_obs=n_obs+nparm

c for lin_d
      jhfj=n_obs+1
      n_obs=n_obs+nparm
      jfhfj=n_obs+1
      n_obs=n_obs+nparm

c for omega functional
      jelo2=n_obs+1
      n_obs=n_obs+1
      jelohfj=n_obs+1
      n_obs=n_obs+nparm

      if(idtask.eq.0) then
        do i=1,nparm
          h_sr_sym(i)=0.5d0*(h_sr(i,1)+obs_tot(jhfj+i-1,1)-obs_tot(jfj+i-1,1)*obs_tot(jelo,1))
        enddo
      endif

c loop vec
      do ivec=1,nvec

      call MPI_BCAST(psi(1,ivec),ndim,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      do iconf=1,nconf_n
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_o(1,iconf,1),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       spsiloc(i+i0,ivec)=omega*omega*ddot(nconf_n,aux(1),1,sr_o(i,1,1),mparm)
     &                         -omega*ddot(nconf_n,aux(1),1,sr_ho(i,1),mparm)
      enddo

      do iconf=1,nconf_n
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_ho(1,iconf),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       spsiloc(i+i0,ivec)=spsiloc(i+i0,ivec)-omega*ddot(nconf_n,aux(1),1,sr_o(i,1,1),mparm)
     &                                            +ddot(nconf_n,aux(1),1,sr_ho(i,1),mparm)
      enddo

      call MPI_REDUCE(spsiloc(1+i0,ivec),spsi(1+i0,ivec),nparm,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ier)

      if(idtask.eq.0) then

        do i=1,nparm
         spsi(i+i0,ivec)=spsi(i+i0,ivec)/obs_tot(1,1)
        enddo

        if(i0.eq.1) then
          aux0=ddot(nparm,psi(1+i0,ivec),1,obs_tot(jfj,1),1)
          do i=1,nparm
           spsi(i+i0,ivec)=spsi(i+i0,ivec)-omega*omega*aux0*obs_tot(jfj+i-1,1)
          enddo

          spsi(1,ivec)=omega*omega*psi(1,ivec)

          aux1=ddot(nparm,psi(i0+1,ivec),1,obs_tot(jefj,1),1)
          aux2=ddot(nparm,psi(i0+1,ivec),1,obs_tot(jhfj,1),1)
          aux3=ddot(nparm,psi(i0+1,ivec),1,obs_tot(jelohfj,1),1)
          do i=1,nparm
           spsi(i+i0,ivec)=spsi(i+i0,ivec)-2*omega*(aux0*obs_tot(jelo,1)*obs_tot(jfj+i-1,1)
     &                                    -0.5d0*(aux0*(obs_tot(jefj+i-1,1)+obs_tot(jhfj+i-1,1))
     &                                    +obs_tot(jfj+i-1,1)*(aux1+aux2)))
           spsi(i+i0,ivec)=spsi(i+i0,ivec)-(aux3*obs_tot(jfj+i-1,1)+aux0*obs_tot(jelohfj+i-1,1))
     &                                    +obs_tot(jelo2,1)*aux0*obs_tot(jfj+i-1,1)
          enddo

          spsi(1,ivec)=spsi(1,ivec)-2*omega*(obs_tot(jelo,1)*psi(1,ivec)+ddot(nparm,h_sr_sym,1,psi(2,ivec),1))
     &                             +obs_tot(jelo2,1)*psi(1,ivec)+(aux3-obs_tot(jelo2,1)*aux0)
          do i=1,nparm
           spsi(i+i0,ivec)=spsi(i+i0,ivec)-2*omega*h_sr_sym(i)*psi(1,ivec)
     &                                    +(obs_tot(jelohfj+i-1,1)-obs_tot(jelo2,1)*obs_tot(jfj+i-1,1))*psi(1,ivec)
          enddo
        endif

      endif
c end loop vec
      enddo

c     do i=1,nparm+i0
c       write(ounit,'(''SPSI_LIN'',100e12.3)')(spsi(i,ivec),ivec=1,nvec)
c     enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine h_psi_varmin(ndim,nvec,psi,hpsi )
      use contrl_file, only: errunit,ounit
      use mpi
      use mpiconf, only: idtask
      use optwf_control, only: ioptjas,ioptorb,nparm,nvecx
      use optwf_func, only: ifunc_omega
      use precision_kinds, only: dp
      use sr_mat_n, only: elocal,h_sr,jefj,jfj,jhfj,nconf_n,obs_tot
      use sr_mat_n, only: s_diag,sr_ho,sr_o,wtg
      use sr_mod,  only: mconf,mparm
      ! these were not called in the master
      ! but they seem to be needed
      ! use sr_index, only: jelo, jelo2, jelohfj

      implicit none

      integer :: i, i0, iconf, ier, ivec
      integer :: jelo, jelo2, jelohfj, jfhfj
      integer :: jfifj, jwtg, k, n_obs
      integer :: ndim, nvec
      real(dp) :: auxx0, auxx2, auxx3, hoz
      real(dp) :: oz, var
      real(dp), dimension(mparm,*) :: psi
      real(dp), dimension(mparm,*) :: hpsi
      real(dp), dimension(mparm,nvecx) :: hpsiloc
      real(dp), dimension(mconf) :: aux0
      real(dp), dimension(mconf) :: aux1
      real(dp), dimension(mconf) :: aux2
      real(dp), dimension(mparm) :: grad_ene

      i0=1
      if(ioptorb.eq.0.and.ioptjas.eq.0) i0=0
      nparm=ndim-i0

c     do i=1,nparm+i0
c       write(ounit,'(''PSI NEW'',100e12.3)')(psi(i,ivec),ivec=1,nvec)
c     enddo

c     write(ounit,*) 'HPSI_LIN',ndim,nvec

      jwtg=1
      jelo=2
      n_obs=2
      jfj=n_obs+1
      n_obs=n_obs+nparm
      jefj=n_obs+1
      n_obs=n_obs+nparm
      jfifj=n_obs+1
      n_obs=n_obs+nparm

c for lin_d
      jhfj=n_obs+1
      n_obs=n_obs+nparm
      jfhfj=n_obs+1
      n_obs=n_obs+nparm

c for omega functional
      jelo2=n_obs+1
      n_obs=n_obs+1
      jelohfj=n_obs+1
      n_obs=n_obs+nparm

      if(idtask.eq.0) then

        var=obs_tot(jelo2,1)-obs_tot(jelo,1)*obs_tot(jelo,1)
        do k=1,nparm
          grad_ene(k)=2*(obs_tot(jefj+k-1,1)-obs_tot(jfj+k-1,1)*obs_tot(jelo,1))
        enddo

      endif
      call MPI_BCAST(var,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)

c loop vec
      do ivec=1,nvec

      do i=1,ndim
        hpsi(i,ivec)=0
      enddo

      call MPI_BCAST(psi(1,ivec),ndim,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      do iconf=1,nconf_n
        hoz=ddot(nparm,psi(i0+1,ivec),1,sr_ho(1,iconf),1)
        oz =ddot(nparm,psi(i0+1,ivec),1,sr_o(1,iconf,1),1)
        aux0(iconf)=(hoz-oz*elocal(iconf,1))*wtg(iconf,1)
        aux1(iconf)=(oz*elocal(iconf,1)**2-hoz*elocal(iconf,1))*wtg(iconf,1)
        aux2(iconf)=oz*wtg(iconf,1)
      enddo
      do i=1,nparm
        hpsiloc(i+i0,ivec)=ddot(nconf_n,aux0(1),1,sr_ho(i,1),mparm)
     &                    +ddot(nconf_n,aux2(1),1,sr_o(i,1,1),mparm)*var
     &                    +ddot(nconf_n,aux1(1),1,sr_o(i,1,1),mparm)
      enddo
      call MPI_REDUCE(hpsiloc(1+i0,ivec),hpsi(1+i0,ivec),nparm,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,i)

      if(idtask.eq.0) then

        do i=1,nparm
         hpsi(i+i0,ivec)=hpsi(i+i0,ivec)/obs_tot(1,1)
        enddo

        if(ifunc_omega.eq.1) then
          auxx0=ddot(nparm,psi(1+i0,ivec),1,grad_ene(1),1)
          auxx2=ddot(nparm,psi(1+i0,ivec),1,obs_tot(jhfj,1),1)
          auxx3=ddot(nparm,psi(1+i0,ivec),1,obs_tot(jefj,1),1)
          do i=1,nparm
           hpsi(i+i0,ivec)=hpsi(i+i0,ivec)
     &              +grad_ene(i)*auxx0-(obs_tot(jhfj+i-1,1)-obs_tot(jefj+i-1,1))*auxx0
     &              -grad_ene(i)*(auxx2-auxx3)
          enddo
        endif

        auxx0=ddot(nparm,psi(1+i0,ivec),1,obs_tot(jfj,1),1)*var
        do i=1,nparm
          hpsi(i+i0,ivec)=hpsi(i+i0,ivec)-auxx0*obs_tot(jfj+i-1,1)
        enddo

        hpsi(1,ivec)=(var*psi(1,ivec)-0.5*ddot(nparm,h_sr(1,1),1,psi(2,ivec),1))
        do i=1,nparm
          hpsi(i+i0,ivec)=hpsi(i+i0,ivec)-0.5*h_sr(i,1)*psi(1,ivec)+s_diag(1,1)*psi(i+i0,ivec) !!!
        enddo

      endif

c end loop vec
      enddo

      do i=1,nparm+i0
        write(ounit,'(''HPSI_LIN'',100e12.3)')(hpsi(i,ivec),ivec=1,nvec)
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine g_psi_lin_d( ndim, nvec, nb1, psi, ew )

      use mpi
      use optwf_control, only: ioptjas,ioptorb
      use precision_kinds, only: dp
      use sr_mat_n, only: jefj,jfj,jhfj,obs_tot,s_diag
      use sr_mod,  only: mparm

      ! these were not called in the master
      ! but they seem to be needed
      ! use sr_index, only: jelo

      ! this was not in master but is clearly needed

      implicit none

      integer :: i, i0, ivec, jelo, jfhfj
      integer :: jfifj, jwtg, k, n_obs
      integer :: nb1, ndim, nparm, nvec
      real(dp), dimension(mparm,*) :: psi
      real(dp), dimension(*) :: ew
      real(dp), dimension(mparm) :: s
      real(dp), dimension(mparm) :: h

      i0=1
      if(ioptorb.eq.0.and.ioptjas.eq.0) i0=0
      nparm=ndim-i0

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

      if(i0.eq.1) then

      s(1)=1
      h(1)=obs_tot(jelo,1)

      do k=1,nparm
       s(k+i0)=obs_tot(jfifj+k-1,1)-obs_tot(jfj+k-1,1)*obs_tot(jfj+k-1,1)
       h(k+i0)=obs_tot(jfhfj+k-1,1)-obs_tot(jfj+k-1,1)*(obs_tot(jhfj+k-1,1)+obs_tot(jefj+k-1,1)-obs_tot(jelo,1)*obs_tot(jfj+k-1,1))
      enddo

      do ivec=1,nvec
        do i=1,ndim
c         if(i.ne.ivec+nb1-1) psi(i,ivec)=psi(i,ivec)/(h(i)+s_diag(1,1)-ew(ivec)*s(i))
          psi(i,ivec)=psi(i,ivec)/(h(i)+s_diag(1,1)-ew(ivec)*s(i))
        enddo
      enddo

      else

      do k=1,nparm
       s(k)=obs_tot(jfifj+k-1,1)
       h(k)=obs_tot(jfhfj+k-1,1)
      enddo

      do ivec=1,nvec
        do i=1,ndim
          psi(i,ivec)=psi(i,ivec)/(h(i)+s_diag(1,1)-ew(ivec)*s(i))
        enddo
      enddo

      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine compute_overlap_psi(ndim,nvec,psi,overlap_psi,anorm)

      use csfs,    only: nstates
      use mpi
      use mpiconf, only: idtask,nproc
      use mstates_mod, only: MSTATES
      use optwf_control, only: ioptjas,ioptorb,nparm,nvecx
      use precision_kinds, only: dp
      use sr_mat_n, only: nconf_n,obs_tot,sr_o,wtg
      use sr_mod,  only: mparm

      implicit none

      integer :: i0, iconf, ier, istate, ivec
      integer :: ndim, nvec
      real(dp) :: den, dum, ratio
      real(dp), dimension(mparm,*) :: psi
      real(dp), dimension(nvecx,*) :: overlap_psi
      real(dp), dimension(*) :: anorm
      real(dp), dimension(nvecx,MSTATES) :: overlap_psiloc
      real(dp), dimension(nvecx) :: anorm_loc

      i0=1
      if(ioptjas+ioptorb.eq.0) i0=0
      nparm=ndim-i0
      if (nproc > 1) then
        do ivec=1,nvec
          call MPI_BCAST(psi(1,ivec),ndim,MPI_REAL8,0,MPI_COMM_WORLD,ier)
        enddo
      endif

      ratio=1.d0
      do ivec=1,nvec
        anorm_loc(ivec)=0.d0
        do istate=1,nstates
          overlap_psiloc(ivec,istate)=0.d0
        enddo
        do iconf=1,nconf_n
          dum=ddot(nparm,psi(i0+1,ivec),1,sr_o(1,iconf,1),1)
          anorm_loc(ivec)=anorm_loc(ivec)+dum*dum*wtg(iconf,1)
          overlap_psiloc(ivec,1)=overlap_psiloc(ivec,1)+dum*wtg(iconf,1)

          do istate=2,nstates
            if(i0.eq.0) ratio=sr_o(nparm+1,iconf,1)/sr_o(nparm+istate,iconf,1)
            !write(ounit,'(A,i4,f15.10)') "state, sr_o(nparm+istate", istate, sr_o(nparm+istate,iconf,1)
            overlap_psiloc(ivec,istate)=overlap_psiloc(ivec,istate)+dum*ratio*wtg(iconf,istate)
          enddo
        enddo
      enddo

      call MPI_REDUCE(overlap_psiloc,overlap_psi,nvecx*nstates,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ier)
      call MPI_REDUCE(anorm_loc,anorm,nvec,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ier)

      if(idtask.eq.0) then

        do ivec=1,nvec

          if(i0.eq.1) then
            anorm(ivec)=anorm(ivec)/obs_tot(1,1)
            overlap_psi(ivec,1)=overlap_psi(ivec,1)/obs_tot(1,1)

            den=dabs(psi(i0,ivec))
            overlap_psi(ivec,1)=dsqrt(dabs(anorm(ivec)-overlap_psi(ivec,1)*overlap_psi(ivec,1)))/den
           else
            anorm(ivec)=dsqrt(anorm(ivec)/obs_tot(1,1))
            do istate=1,nstates
              den=dsqrt(obs_tot(1,istate)*obs_tot(1,1))*anorm(ivec)
              overlap_psi(ivec,istate)=dabs(overlap_psi(ivec,istate))/den
            enddo
          endif
        enddo

      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine h_psi_lin_d(ndim,nvec,psi,hpsi)

      use optwf_func, only: ifunc_omega
      use precision_kinds, only: dp
      use sr_mod,  only: mparm

      implicit none

      integer :: ndim, nvec
      real(dp), dimension(mparm,*) :: psi
      real(dp), dimension(mparm,*) :: hpsi

      if(ifunc_omega.eq.0) then
        call h_psi_energymin(ndim,nvec,psi,hpsi)
       elseif(ifunc_omega.le.2) then
        call h_psi_varmin(ndim,nvec,psi,hpsi)
       else
        call h_psi_omegamin(ndim,nvec,psi,hpsi)
      endif

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine s_psi_lin_d(ndim,nvec,psi,spsi)

      use optwf_func, only: ifunc_omega
      use precision_kinds, only: dp
      use sr_mod,  only: mparm

      implicit none

      integer :: ndim, nvec
      real(dp), dimension(mparm,*) :: psi
      real(dp), dimension(mparm,*) :: spsi

      if(ifunc_omega.le.2) then
        call s_psi_energymin(ndim,nvec,psi,spsi )
       else
        call s_psi_omegamin(ndim,nvec,psi,spsi )
      endif

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine select_ci_root(iroot)

      use csfs,    only: ccsf,ncsf
      use slater,  only: cdet,ndet

      implicit none

      integer :: i, iadiag, icsf, iroot


      do i=1,ndet
        cdet(i,1,1)=cdet(i,iroot,1)
      enddo

      do icsf=1,ncsf
        ccsf(icsf,1,iadiag)=ccsf(icsf,iroot,1)
      enddo

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end module
