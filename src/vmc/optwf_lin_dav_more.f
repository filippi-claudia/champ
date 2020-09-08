      subroutine lin_d(nparm,nvec,nvecx,deltap,deltap_more,adiag,ethr)

      use force, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use const, only: ipr
      use mstates_mod, only: MSTATES
      use csfs, only: nstates

      use mpiconf, only: idtask
      use optwf_contrl, only: ioptjas, ioptorb
      use optwf_corsam, only: energy, force
      use optwf_parms, only: nparmd, nparmj
      use sr_mat_n, only: jfj
      use sr_mat_n, only: obs_tot
    
      implicit real*8(a-h,o-z)

      include 'mpif.h'
      include 'sr.h'




      dimension e(MVEC),evc(MPARM,MVEC),itype(MVEC),overlap_psi(MVEC,MSTATES),index_overlap(MVEC),anorm(MVEC)
      dimension deltap(*),deltap_more(MPARM*MSTATES,5)

      call p2gtid('optwf:lin_jdav',lin_jdav,0,1)

      write(6,*) 'LIN_D NPARM',nparm

      call sr_hs(nparm,adiag)

      i0=1
      if(ioptorb+ioptjas.eq.0) i0=0
      nparm_p1=nparm+i0

      do 10 ivec=1,nvec
        do 10 iparm=1,nparm_p1
 10       evc(iparm,ivec)=0.d0

      do 20 ivec=1,nvec
        itype(ivec)=1
 20     evc(ivec,ivec)=1.d0

      if(lin_jdav.eq.0) then

       write(6,*) "USING OLD REGTERG"

        call regterg( nparm_p1, MPARM, nvec, nvecx, evc, ethr,
     &                e, itype, notcnv, idav_iter, ipr, idtask )

        write(6,'(''LIN_D: no. iterations'',i4)') idav_iter
        write(6,'(''LIN_D: no. not converged roots '',i4)') notcnv

       elseif(lin_jdav.eq.1) then
       write(6,*) "USING DAVIDSON WRAP: FREE VERSION"
        call davidson_wrap( nparm_p1, MPARM, nvec, nvecx, MVEC, evc, 
     &       ethr, e, itype, notcnv, idav_iter, ipr)
       else
         call fatal_error('LIND: lin_jdav < 2')
      endif

      call my_second(2,'david ')
      call compute_overlap_psi(nparm_p1,nvec,evc,overlap_psi,anorm)
c idtask.eq.0
      if(idtask.eq.0)  then

        do istate=1,nstates
          do ivec=1,nvec
            write(6,'(''LIN_D: state, vec, energy'',2i4,2f12.5)') istate,ivec,e(ivec),overlap_psi(ivec,istate)
          enddo
        enddo

        if(i0.eq.1) then
          ortho_min=1.d+99
          do ivec=1,nvec
            if(overlap_psi(ivec,1).lt.ortho_min) then
              ortho_min=overlap_psi(ivec,1)
              i_ortho_min=ivec
            endif
          enddo
          write(6,'(''LIN_D: max overlap ivec'',i4)') i_ortho_min

          do i=1,nparm
            deltap(i)=evc(i+1,i_ortho_min)/evc(1,i_ortho_min)
          enddo

          bot=1
          do i=1,nparmd
            bot=bot-deltap(nparmj+i)*obs_tot(jfj+nparmj +i-1,1)
          enddo

          do i=1,nparm
            deltap(i)=deltap(i)/bot
          enddo

         else                   
c else means if i optimize jastrow and or orbitals

          do istate=1,nstates
            call sort(nvec,overlap_psi(1,istate),index_overlap)
            i_overlap_max=index_overlap(nvec)
            write(6,'(''LIN_D: state, max overlap ivec'',2i4)') istate,i_overlap_max

            do i=1,nparm
              deltap(i+nparm*(istate-1))=evc(i,i_overlap_max)/anorm(i_overlap_max)
            enddo

c save 5 additional vectors with large overlap
            do ivec=1,5
              idx_ivec=index_overlap(nvec-ivec)
              do i=1,nparm
                deltap_more(i+nparm*(istate-1),ivec)=evc(i,idx_ivec)/anorm(idx_ivec)
              enddo
            enddo
          enddo

        endif

c endif idtask.eq.0
      endif

c     do istate=1,nstates
c       call MPI_BCAST(deltap(1+nparm*(istate-1)),nparm,MPI_REAL8,0,MPI_COMM_WORLD,ier)
c     enddo

      call MPI_BCAST(deltap,nparm*nstates,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      if(i0.eq.0) then
        do ivec=1,5
          call MPI_BCAST(deltap_more(1,ivec),nparm*nstates,MPI_REAL8,0,MPI_COMM_WORLD,ier)
        enddo
      endif


      return              ! deltap
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine amul(n,q,r)

      use jd_scratch, only: qr, rr
      implicit real*8(a-h,o-z)


      include 'mpif.h'
      include 'sr.h'


      complex*16 q(n),r(n)
      do i=1,n
        qr(i)=real(q(i))
      enddo

      call h_psi_lin_d(n,1,qr,rr)

      call MPI_BCAST(rr,n,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      do i=1,n
        q(i)=cmplx(qr(i),0.d0)
        r(i)=cmplx(rr(i),0.d0)
      enddo

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bmul(n,q,r)

      use jd_scratch, only: qr, rr
      implicit real*8(a-h,o-z)


      include 'mpif.h'
      include 'sr.h'


      complex*16 q(n),r(n)
      do i=1,n
        qr(i)=real(q(i))
      enddo

      call s_psi_lin_d(n,1,qr,rr)

      call MPI_BCAST(rr,n,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      do i=1,n
        q(i)=cmplx(qr(i),0.d0)
        r(i)=cmplx(rr(i),0.d0)
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine precon(n,q)
      implicit real*8 (a-h,o-z)

      complex*16 q(n)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine h_psi_energymin(ndim,nvec,psi,hpsi )
      use mpiconf, only: idtask
      use optwf_contrl, only: ioptjas, ioptorb, nparm
      use sr_mat_n, only: h_sr, jefj, jfj, jhfj, nconf_n, s_diag, sr_ho
      use sr_mat_n, only: sr_o, wtg, obs_tot
      implicit real*8(a-h,o-z)




      include 'mpif.h'
      include 'sr.h'



      dimension psi(MPARM,*),hpsi(MPARM,*),aux(MCONF),hpsiloc(MPARM,MVEC)


      i0=1
      if(ioptorb.eq.0.and.ioptjas.eq.0) i0=0
      nparm=ndim-i0

c     write(6,*) 'HPSI_LIN',ndim,nvec

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
       hpsiloc(i+i0,ivec)=0.5d0*ddot(nconf_n,aux(1),1,sr_o(i,1),MPARM)
      enddo

      do iconf=1,nconf_n
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_o(1,iconf),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       hpsiloc(i+i0,ivec)=hpsiloc(i+i0,ivec)+0.5d0*ddot(nconf_n,aux(1),1,sr_ho(i,1),MPARM)
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
           hpsi(i+i0,ivec)=hpsi(i+i0,ivec)+h_sr(i)*psi(1,ivec)+s_diag(1,1)*psi(i+i0,ivec) !!!
          enddo

        endif

      endif

c end loop vec
      enddo

c     do i=1,nparm+1
c       write(6,'(''HPSI_LIN'',100e12.3)')(hpsi(i,ivec),ivec=1,nvec)
c     enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine s_psi_energymin(ndim,nvec,psi,spsi )
      use mpiconf, only: idtask
      use optwf_contrl, only: ioptjas, ioptorb, nparm
      use sr_mat_n, only: jefj, jfj, jhfj, nconf_n
      use sr_mat_n, only: sr_o, wtg, obs_tot
      implicit real*8(a-h,o-z)




      include 'mpif.h'
      include 'sr.h'



      dimension psi(MPARM,*),spsi(MPARM,*),spsiloc(MPARM,MVEC),aux(MCONF)

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
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_o(1,iconf),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       spsiloc(i+i0,ivec)=ddot(nconf_n,aux(1),1,sr_o(i,1),MPARM)
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
c       write(6,'(''SPSI_LIN'',100e12.3)')(spsi(i,ivec),ivec=1,nvec)
c     enddo
c     STOP

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine h_psi_omegamin(ndim,nvec,psi,hpsi )
      use mpiconf, only: idtask
      use optwf_contrl, only: ioptjas, ioptorb, nparm
      use optwf_func, only: omega
      use sr_mat_n, only: h_sr, jefj, jfj, jhfj, nconf_n, s_diag, sr_ho
      use sr_mat_n, only: sr_o, wtg, obs_tot
      implicit real*8(a-h,o-z)





      include 'mpif.h'
      include 'sr.h'
 




      dimension psi(MPARM,*),hpsi(MPARM,*),hpsiloc(MPARM,MVEC),aux(MCONF)

      i0=1
      if(ioptorb.eq.0.and.ioptjas.eq.0) i0=0
      nparm=ndim-i0

c     do i=1,nparm+i0
c       write(6,'(''PSI NEW'',100e12.3)')(psi(i,ivec),ivec=1,nvec)
c     enddo

c     write(6,*) 'HPSI_LIN',ndim,nvec

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
       hpsiloc(i+i0,ivec)=-0.5d0*ddot(nconf_n,aux(1),1,sr_o(i,1),MPARM)
      enddo

c     do i=1,nparm+1
c       write(6,*) 'H1 ',hpsi(i,ivec)
c     enddo

      do iconf=1,nconf_n
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_o(1,iconf),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       hpsiloc(i+i0,ivec)=hpsiloc(i+i0,ivec)-0.5d0*ddot(nconf_n,aux(1),1,sr_ho(i,1),MPARM)
     &                                      +omega*ddot(nconf_n,aux(1),1,sr_o(i,1),MPARM)
      enddo

      call MPI_REDUCE(hpsiloc(1+i0,ivec),hpsi(1+i0,ivec),nparm,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ier)

c     do i=1,nparm+1
c       write(6,*) 'H2 ',hpsi(i,ivec)
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
           hpsi(i+i0,ivec)=hpsi(i+i0,ivec)-h_sr(i)*psi(1,ivec)+s_diag(1,1)*psi(i+i0,ivec) !!!
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
c       write(6,'(''HPSI_LIN'',100e12.3)')(hpsi(i,ivec),ivec=1,nvec)
c     enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine s_psi_omegamin(ndim,nvec,psi,spsi )
      use mpiconf, only: idtask
      use optwf_contrl, only: ioptjas, ioptorb, nparm
      use optwf_func, only: omega
      use sr_mat_n, only: h_sr, jefj, jfj, jhfj, nconf_n, sr_ho
      use sr_mat_n, only: sr_o, wtg, obs_tot
      implicit real*8(a-h,o-z)





      include 'mpif.h'
      include 'sr.h'





      dimension psi(MPARM,*),spsi(MPARM,*),spsiloc(MPARM,MVEC),aux(MCONF),h_sr_sym(MPARM)

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
          h_sr_sym(i)=0.5d0*(h_sr(i)+obs_tot(jhfj+i-1,1)-obs_tot(jfj+i-1,1)*obs_tot(jelo,1))
        enddo
      endif

c loop vec
      do ivec=1,nvec

      call MPI_BCAST(psi(1,ivec),ndim,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      do iconf=1,nconf_n
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_o(1,iconf),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       spsiloc(i+i0,ivec)=omega*omega*ddot(nconf_n,aux(1),1,sr_o(i,1),MPARM)
     &                         -omega*ddot(nconf_n,aux(1),1,sr_ho(i,1),MPARM)
      enddo

      do iconf=1,nconf_n
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_ho(1,iconf),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       spsiloc(i+i0,ivec)=spsiloc(i+i0,ivec)-omega*ddot(nconf_n,aux(1),1,sr_o(i,1),MPARM)
     &                                            +ddot(nconf_n,aux(1),1,sr_ho(i,1),MPARM)
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
c       write(6,'(''SPSI_LIN'',100e12.3)')(spsi(i,ivec),ivec=1,nvec)
c     enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine h_psi_varmin(ndim,nvec,psi,hpsi )
      use mpiconf, only: idtask
      use optwf_contrl, only: ioptjas, ioptorb, nparm
      use optwf_func, only: ifunc_omega, omega
      use sr_mat_n, only: elocal, h_sr, jefj, jfj, jhfj, nconf_n, s_diag, sr_ho
      use sr_mat_n, only: sr_o, wtg, obs_tot
      implicit real*8(a-h,o-z)





      include 'mpif.h'
      include 'sr.h'

 


      dimension psi(MPARM,*),hpsi(MPARM,*),hpsiloc(MPARM,MVEC),aux0(MCONF),aux1(MCONF),aux2(MCONF)
      dimension grad_ene(MPARM)

      i0=1
      if(ioptorb.eq.0.and.ioptjas.eq.0) i0=0
      nparm=ndim-i0

c     do i=1,nparm+i0
c       write(6,'(''PSI NEW'',100e12.3)')(psi(i,ivec),ivec=1,nvec)
c     enddo

c     write(6,*) 'HPSI_LIN',ndim,nvec

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
        oz =ddot(nparm,psi(i0+1,ivec),1,sr_o(1,iconf),1)
        aux0(iconf)=(hoz-oz*elocal(iconf,1))*wtg(iconf,1)
        aux1(iconf)=(oz*elocal(iconf,1)**2-hoz*elocal(iconf,1))*wtg(iconf,1)
        aux2(iconf)=oz*wtg(iconf,1)
      enddo
      do i=1,nparm
        hpsiloc(i+i0,ivec)=ddot(nconf_n,aux0(1),1,sr_ho(i,1),MPARM)
     &                    +ddot(nconf_n,aux2(1),1,sr_o(i,1),MPARM)*var
     &                    +ddot(nconf_n,aux1(1),1,sr_o(i,1),MPARM)
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

        hpsi(1,ivec)=(var*psi(1,ivec)-0.5*ddot(nparm,h_sr(1),1,psi(2,ivec),1))
        do i=1,nparm
          hpsi(i+i0,ivec)=hpsi(i+i0,ivec)-0.5*h_sr(i)*psi(1,ivec)+s_diag(1,1)*psi(i+i0,ivec) !!!
        enddo

      endif

c end loop vec
      enddo

      do i=1,nparm+i0
        write(6,'(''HPSI_LIN'',100e12.3)')(hpsi(i,ivec),ivec=1,nvec)
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine g_psi_lin_d( ndim, nvec, nb1, psi, ew )

      use sr_mat_n, only: jefj, jfj, jhfj, s_diag
      use sr_mat_n, only: obs_tot
      implicit real*8(a-h,o-z)


      include 'mpif.h'
      include 'sr.h'


      dimension psi(MPARM,*),ew(*)
      dimension s(MPARM),h(MPARM)

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
      use force, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use csfs, only: nstates
      use mstates_mod, only: MSTATES

      use mpiconf, only: idtask, nproc
      use optwf_contrl, only: ioptjas, ioptorb, nparm
      use sr_mat_n, only: nconf_n
      use sr_mat_n, only: sr_o, wtg, obs_tot
      implicit real*8(a-h,o-z)





      include 'mpif.h'
      include 'sr.h'


      dimension psi(MPARM,*),overlap_psi(MVEC,*),anorm(*),overlap_psiloc(MVEC,MSTATES),anorm_loc(MVEC)

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
          dum=ddot(nparm,psi(i0+1,ivec),1,sr_o(1,iconf),1)
          anorm_loc(ivec)=anorm_loc(ivec)+dum*dum*wtg(iconf,1)
          overlap_psiloc(ivec,1)=overlap_psiloc(ivec,1)+dum*wtg(iconf,1)

          do istate=2,nstates
            if(i0.eq.0) ratio=sr_o(nparm+1,iconf)/sr_o(nparm+istate,iconf)
            overlap_psiloc(ivec,istate)=overlap_psiloc(ivec,istate)+dum*ratio*wtg(iconf,istate)
          enddo
        enddo
      enddo

      call MPI_REDUCE(overlap_psiloc,overlap_psi,MVEC*nstates,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ier)
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
