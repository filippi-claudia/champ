      subroutine lin_d(nparm,nvec,nvecx,deltap,adiag,ethr)

      implicit real*8(a-h,o-z)

      include 'sr.h'
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm_sav
      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj
      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES)
     &,s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      common /optwf_corsam/ add_diag(MFORCE),energy(MFORCE),energy_err(MFORCE),force(MFORCE),force_err(MFORCE)

      dimension e(MVEC),evc(MPARM,MVEC),itype(MVEC),overlap_psi(MVEC,MSTATES),anorm(MVEC)
      dimension deltap(*)

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

        call regterg( nparm_p1, MPARM, nvec, nvecx, evc, ethr,
     &                e, itype, notcnv, idav_iter , ipr)

      else

c e0,nvec_e0
        e0=energy(1)
        nvec_e0=nvec/2
        write(6,'(''LIN_D : target energy'',f14.6)') e0

        call jdqz_driver( nparm_p1, nvec_e0, nvec, nvecx, evc, ethr,
     &                    e, e0, itype, notcnv, idav_iter , ipr )

      endif

c TMP
c     STOP

      call compute_overlap_psi(nparm_p1,nvec,evc,overlap_psi,anorm)

      write(6,'(''LIN_D: no. iterations'',i4)') idav_iter
      write(6,'(''LIN_D: no. not converged roots '',i4)') notcnv
      
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
          bot=bot-deltap(nparmj+i)*obs(jfj+nparmj +i-1,1)
        enddo

        do i=1,nparm
          deltap(i)=deltap(i)/bot
        enddo
       else
       
        do istate=1,nstates
          overlap_max=0.d0
          do ivec=1,nvec
            if(overlap_psi(ivec,istate).gt.overlap_max) then
              overlap_max=overlap_psi(ivec,istate)
              i_overlap_max=ivec
            endif
          enddo
          write(6,'(''LIN_D: state, max overlap ivec'',2i4)') istate,i_overlap_max

          do i=1,nparm
            deltap(i+nparm*(istate-1))=evc(i,i_overlap_max)/anorm(i_overlap_max)
          enddo
        enddo
      endif
 

      return              ! deltap
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine amul(n,q,r)

      implicit real*8 (a-h,o-z)

      include 'sr.h'

      common /jd_scratch/ qr(MPARM),rr(MPARM)

      complex*16 q(n),r(n)
      do i=1,n
        qr(i)=real(q(i))
      enddo

      call h_psi_lin_d(n,1,qr,rr)

      do i=1,n
        r(i)=cmplx(rr(i),0.d0)
      enddo

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bmul(n,q,r)

      implicit real*8 (a-h,o-z)

      include 'sr.h'

      common /jd_scratch/ qr(MPARM),rr(MPARM)

      complex*16 q(n),r(n)
      do i=1,n
        qr(i)=real(q(i))
      enddo

      call s_psi_lin_d(n,1,qr,rr)

      do i=1,n
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
      implicit real*8 (a-h,o-z)

      include 'sr.h'
      include 'mstates.h'
 
      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm_sav
      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES)
     &,s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf


      dimension psi(MPARM,*),hpsi(MPARM,*),aux(MCONF)

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

      do iconf=1,nconf
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_ho(1,iconf),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       hpsi(i+i0,ivec)=0.5d0*ddot(nconf,aux(1),1,sr_o(i,1),MPARM)/obs(1,1)
      enddo

c     do i=1,nparm+1
c       write(6,*) 'H1 ',hpsi(i,ivec)
c     enddo

      do iconf=1,nconf
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_o(1,iconf),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       hpsi(i+i0,ivec)=hpsi(i+i0,ivec)+0.5d0*ddot(nconf,aux(1),1,sr_ho(i,1),MPARM)/obs(1,1)
      enddo

c     do i=1,nparm+1
c       write(6,*) 'H2 ',hpsi(i,ivec)
c     enddo

c TEST
      if(i0.eq.1) then

        aux0=ddot(nparm,psi(i0+1,ivec),1,obs(jfj,1),1)
        aux1=ddot(nparm,psi(i0+1,ivec),1,obs(jefj,1),1)
        aux2=ddot(nparm,psi(i0+1,ivec),1,obs(jhfj,1),1)
        do i=1,nparm
         hpsi(i+i0,ivec)=hpsi(i+i0,ivec)+aux0*obs(jelo,1)*obs(jfj+i-1,1)
     &                            -0.5d0*(aux0*(obs(jefj+i-1,1)+obs(jhfj+i-1,1))
     &                                   +obs(jfj+i-1,1)*(aux1+aux2))
        enddo

        hpsi(1,ivec)=obs(jelo,1)*psi(1,ivec)+ddot(nparm,h_sr,1,psi(2,ivec),1)
        do i=1,nparm
         hpsi(i+i0,ivec)=hpsi(i+i0,ivec)+h_sr(i)*psi(1,ivec)+s_diag(1,1)*psi(i+i0,ivec) !!!
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

      subroutine s_psi_energymin(ndim,nvec,psi,spsi )
      implicit real*8 (a-h,o-z)

      include 'sr.h'
      include 'mstates.h'

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm_sav
      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES)
     &,s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf


      dimension psi(MPARM,*),spsi(MPARM,*),aux(MCONF)

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

      do iconf=1,nconf
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_o(1,iconf),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       spsi(i+i0,ivec)=ddot(nconf,aux(1),1,sr_o(i,1),MPARM)/obs(1,1)
      enddo

      if(i0.eq.1) then
        aux0=ddot(nparm,psi(1+i0,ivec),1,obs(jfj,1),1)
        do i=1,nparm
         spsi(i+i0,ivec)=spsi(i+i0,ivec)-aux0*obs(jfj+i-1,1)
        enddo

        spsi(1,ivec)=psi(1,ivec)

      endif
c end loop vec
      enddo

      do i=1,nparm+i0
        write(6,'(''SPSI_LIN'',100e12.3)')(spsi(i,ivec),ivec=1,nvec)
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine h_psi_omegamin(ndim,nvec,psi,hpsi )
      implicit real*8 (a-h,o-z)

      include 'sr.h'
      include 'mstates.h'
 
      common /optwf_func/ omega,ifunc_omega

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm_sav
      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES)
     &,s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf

      dimension psi(MPARM,*),hpsi(MPARM,*),aux(MCONF)

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

      do iconf=1,nconf
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_ho(1,iconf),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       hpsi(i+i0,ivec)=-0.5d0*ddot(nconf,aux(1),1,sr_o(i,1),MPARM)/obs(1,1)
      enddo

c     do i=1,nparm+1
c       write(6,*) 'H1 ',hpsi(i,ivec)
c     enddo

      do iconf=1,nconf
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_o(1,iconf),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       hpsi(i+i0,ivec)=hpsi(i+i0,ivec)-0.5d0*ddot(nconf,aux(1),1,sr_ho(i,1),MPARM)/obs(1,1)
     &                                +omega*ddot(nconf,aux(1),1,sr_o(i,1),MPARM)/obs(1,1)
      enddo

c     do i=1,nparm+1
c       write(6,*) 'H2 ',hpsi(i,ivec)
c     enddo

c LATER
      if(i0.eq.1) then

        aux0=ddot(nparm,psi(i0+1,ivec),1,obs(jfj,1),1)
        aux1=ddot(nparm,psi(i0+1,ivec),1,obs(jefj,1),1)
        aux2=ddot(nparm,psi(i0+1,ivec),1,obs(jhfj,1),1)
        do i=1,nparm
         hpsi(i+i0,ivec)=hpsi(i+i0,ivec)-aux0*obs(jelo,1)*obs(jfj+i-1,1)
     &                            +0.5d0*(aux0*(obs(jefj+i-1,1)+obs(jhfj+i-1,1))
     &                                   +obs(jfj+i-1,1)*(aux1+aux2))
        enddo

        hpsi(1,ivec)=-(obs(jelo,1)*psi(1,ivec)+ddot(nparm,h_sr,1,psi(2,ivec),1))
        do i=1,nparm
         hpsi(i+i0,ivec)=hpsi(i+i0,ivec)-h_sr(i)*psi(1,ivec)+s_diag(1,1)*psi(i+i0,ivec) !!!
        enddo

        do i=1,nparm
         hpsi(i+i0,ivec)=hpsi(i+i0,ivec)-omega*aux0*obs(jfj+i-1,1)
        enddo
        hpsi(1,ivec)=hpsi(1,ivec)+omega*psi(1,ivec)

      endif

c end loop vec
      enddo

      do i=1,nparm+i0
        write(6,'(''HPSI_LIN'',100e12.3)')(hpsi(i,ivec),ivec=1,nvec)
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine s_psi_omegamin(ndim,nvec,psi,spsi )
      implicit real*8 (a-h,o-z)

      include 'sr.h'
      include 'mstates.h'

      common /optwf_func/ omega,ifunc_omega

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm_sav
      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES)
     &,s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf

      dimension psi(MPARM,*),spsi(MPARM,*),aux(MCONF)

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

c loop vec
      do ivec=1,nvec

      do iconf=1,nconf
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_o(1,iconf),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       spsi(i+i0,ivec)=omega*omega*ddot(nconf,aux(1),1,sr_o(i,1),MPARM)/obs(1,1)
     &                      -omega*ddot(nconf,aux(1),1,sr_ho(i,1),MPARM)/obs(1,1)
      enddo

      do iconf=1,nconf
       aux(iconf)=ddot(nparm,psi(i0+1,ivec),1,sr_ho(1,iconf),1)*wtg(iconf,1)
      enddo
      do i=1,nparm
       spsi(i+i0,ivec)=spsi(i+i0,ivec)-omega*ddot(nconf,aux(1),1,sr_o(i,1),MPARM)/obs(1,1)
     &                                      +ddot(nconf,aux(1),1,sr_ho(i,1),MPARM)/obs(1,1)
      enddo

c LATER
      if(i0.eq.1) then
        aux0=ddot(nparm,psi(1+i0,ivec),1,obs(jfj,1),1)
        do i=1,nparm
         spsi(i+i0,ivec)=spsi(i+i0,ivec)-omega*omega*aux0*obs(jfj+i-1,1)
        enddo

        spsi(1,ivec)=omega*omega*psi(1,ivec)

        aux1=ddot(nparm,psi(i0+1,ivec),1,obs(jefj,1),1)
        aux2=ddot(nparm,psi(i0+1,ivec),1,obs(jhfj,1),1)
        aux3=ddot(nparm,psi(i0+1,ivec),1,obs(jelohfj,1),1)
        do i=1,nparm
         spsi(i+i0,ivec)=spsi(i+i0,ivec)-2*omega*(aux0*obs(jelo,1)*obs(jfj+i-1,1)
     &                                  -0.5d0*(aux0*(obs(jefj+i-1,1)+obs(jhfj+i-1,1))
     &                                  +obs(jfj+i-1,1)*(aux1+aux2)))
         spsi(i+i0,ivec)=spsi(i+i0,ivec)-(aux3*obs(jfj+i-1,1)+aux0*obs(jelohfj+i-1,1))
     &                                  +obs(jelo2,1)*aux0*obs(jfj+i-1,1)
        enddo

        spsi(1,ivec)=spsi(1,ivec)-2*omega*(obs(jelo,1)*psi(1,ivec)+ddot(nparm,h_sr,1,psi(2,ivec),1))
     &                           +obs(jelo2,1)*psi(1,ivec)+(aux3-obs(jelo2,1)*aux0)
        do i=1,nparm
         spsi(i+i0,ivec)=spsi(i+i0,ivec)-2*omega*h_sr(i)*psi(1,ivec)
     &                                  +(obs(jelohfj+i-1,1)-obs(jelo2,1)*obs(jfj+i-1,1))*psi(1,ivec)
        enddo

      endif
c end loop vec
      enddo

      do i=1,nparm+i0
        write(6,'(''SPSI_LIN'',100e12.3)')(spsi(i,ivec),ivec=1,nvec)
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine h_psi_varmin(ndim,nvec,psi,hpsi )
      implicit real*8 (a-h,o-z)

      include 'sr.h'
      include 'mstates.h'

      common /optwf_func/ omega,ifunc_omega

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm_sav
      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES)
     &,s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf

      dimension psi(MPARM,*),hpsi(MPARM,*)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine g_psi_lin_d( ndim, nvec, nb1, psi, ew )

      implicit real*8 (a-h,o-z)

      include 'sr.h'
      include 'mstates.h'

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm_sav
      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES)
     &,s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf


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

c     do i=1,nparm+1
c       write(6,'(''PSI NEW BEFORE G'',100e12.3)')(psi(i,ivec),ivec=1,nvec)
c     enddo

      if(i0.eq.1) then

      s(1)=1
      h(1)=obs(jelo,1)

      do k=1,nparm
       s(k+i0)=obs(jfifj+k-1,1)-obs(jfj+k-1,1)*obs(jfj+k-1,1)
       h(k+i0)=obs(jfhfj+k-1,1)-obs(jfj+k-1,1)*(obs(jhfj+k-1,1)+obs(jefj+k-1,1)-obs(jelo,1)*obs(jfj+k-1,1))
      enddo

      do ivec=1,nvec
        do i=1,ndim
          if(i.ne.ivec+nb1-1) psi(i,ivec)=psi(i,ivec)/(h(i)+s_diag(1,1)-ew(ivec)*s(i)) !!!
        enddo
      enddo

      else

      do k=1,nparm
       s(k)=obs(jfifj+k-1,1)
       h(k)=obs(jfhfj+k-1,1)
      enddo

      do ivec=1,nvec
        do i=1,ndim
          if(i.ne.ivec+nb1-1) psi(i,ivec)=psi(i,ivec)/(h(i)+s_diag(1,1)-ew(ivec)*s(i)) !!!
        enddo
      enddo

      endif

c     do i=1,nparm+1
c       write(6,'(''PSI NEW AFTER G'',100e12.3)')(psi(i,ivec),ivec=1,nvec)
c     enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine compute_overlap_psi(ndim,nvec,psi,overlap_psi,anorm)
      implicit real*8 (a-h,o-z)

      include 'sr.h'
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm_sav
      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES)
     &,s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      dimension psi(MPARM,*),overlap_psi(MVEC,*),anorm(*)

      i0=1
      if(ioptjas+ioptorb.eq.0) i0=0
      nparm=ndim-i0

      ratio=1.d0
      do ivec=1,nvec

        anorm(ivec)=0.d0
        do istate=1,nstates
          overlap_psi(ivec,istate)=0.d0
        enddo

        do iconf=1,nconf
          dum=ddot(nparm,psi(i0+1,ivec),1,sr_o(1,iconf),1)

          anorm(ivec)=anorm(ivec)+dum*dum*wtg(iconf,1)
          overlap_psi(ivec,1)=overlap_psi(ivec,1)+dum*wtg(iconf,1)

          do istate=2,nstates
            if(i0.eq.0) ratio=sr_o(nparm+1,iconf)/sr_o(nparm+istate,iconf)
            overlap_psi(ivec,istate)=overlap_psi(ivec,istate)+dum*ratio*wtg(iconf,istate)
          enddo
        enddo

        if(i0.eq.1) then
          anorm(ivec)=anorm(ivec)/obs(1,1)
          overlap_psi(ivec,1)=dabs(overlap_psi(ivec,1))/obs(1,1)

          den=dabs(psi(i0,ivec))
          overlap_psi(ivec,1)=dsqrt(dabs(anorm(ivec)-overlap_psi(ivec,1)*overlap_psi(ivec,1)))/den

         else
          anorm(ivec)=dsqrt(anorm(ivec)/obs(1,1))
          do istate=1,nstates
            den=dsqrt(obs(1,istate)*obs(1,1))*anorm(ivec)
            overlap_psi(ivec,istate)=dabs(overlap_psi(ivec,istate))/den
          enddo
        endif

      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
