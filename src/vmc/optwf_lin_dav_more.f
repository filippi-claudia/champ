      module optwf_lin_dav_more
      use error,   only: fatal_error
      use sr_more, only: ddot
      contains
      subroutine lin_d(nparm,nvec,nvecx,deltap,deltap_more,index_more,adiag,ethr)

      use contrl_file, only: errunit,ounit
      use control, only: ipr
      use csfs,    only: nstates
      use davidson_wrap_mod, only: davidson_wrap
      use mpi
      use mpiconf, only: idtask
      use mpitimer, only: time,time_check1,time_check2,time_start
      use mstates_mod, only: MSTATES
      use optwf_control, only: ioptjas,ioptorb,lin_jdav
      use optwf_corsam, only: energy
      use optwf_lib, only: sort
      use optwf_lin_dav_extra, only: compute_overlap_psi
      use optwf_parms, only: nparmd,nparmj
      use optwf_sr_mod, only: sr_hs
      use precision_kinds, only: dp
      use regterg_mod, only: regterg
      use sr_mat_n, only: jfj,obs_tot
      use sr_mod,  only: mparm
      implicit none

      integer :: i, i0, i_ortho_min, i_overlap_max, idav_iter
      integer :: idx_ivec, ier, iparm, istate
      integer :: ivec, notcnv, nparm, nparm_p1
      integer :: nvec, nvecx
      integer, dimension(nvecx) :: itype
      integer, dimension(nvecx) :: index_overlap
      integer, dimension(5,MSTATES) :: index_more
      real(dp) :: adiag, bot, ethr, ortho_min
      real(dp), dimension(nvecx) :: e
      real(dp), dimension(mparm,nvecx) :: evc
      real(dp), dimension(nvecx,MSTATES) :: overlap_psi
      real(dp), dimension(nvecx) :: anorm
      real(dp), dimension(*) :: deltap
      real(dp), dimension(mparm*MSTATES,5) :: deltap_more

      ! include 'mpif.h'

      time_check1 = time()

      write(ounit,*) 'LIN_D NPARM',nparm

      call sr_hs(nparm,adiag)

      i0=1
      if(ioptorb+ioptjas.eq.0) i0=0
      nparm_p1=nparm+i0

      do ivec=1,nvec
        do iparm=1,nparm_p1
          evc(iparm,ivec)=0.d0
        enddo
      enddo

      do ivec=1,nvec
        itype(ivec)=1
        evc(ivec,ivec)=1.d0
      enddo

      ! regterg
      if(lin_jdav.eq.0) then
       write(ounit,*) "USING OLD REGTERG"

        call regterg( nparm_p1, mparm, nvec, nvecx, evc, ethr,
     &                e, itype, notcnv, idav_iter, ipr, idtask )

       ! Davidson DPR
       elseif(lin_jdav.eq.1) then
        write(ounit,*) "USING DPR DAVIDSON"
        call davidson_wrap( nparm_p1, mparm, nvec, nvecx, nvecx, evc,
     &       ethr, e, itype, notcnv, idav_iter, ipr, "DPR")

       ! Davidson JOCC
       elseif(lin_jdav.eq.2) then
        write(ounit,*) "USING GJD DAVIDSON"
         call davidson_wrap( nparm_p1, mparm, nvec, nvecx, nvecx, evc,
     &       ethr, e, itype, notcnv, idav_iter, ipr, "GJD")

       else
         call fatal_error('LIND: lin_jdav must be 0, 1 or 2')

      endif

      write(ounit,'(''LIN_D: no. iterations'',i4)') idav_iter
      write(ounit,'(''LIN_D: no. not converged roots '',i4)') notcnv

      time_check2 = time()
      write(ounit, '(a,t40, f12.3, f12.3)') "END OF David, REAL TIME IS", time_check2 - time_start, time_check2 - time_check1
      time_check1 = time_check2


      call compute_overlap_psi(nparm_p1,nvec,evc,overlap_psi,anorm)
c idtask.eq.0
      if(idtask.eq.0)  then

        do istate=1,nstates
          do ivec=1,nvec
            write(ounit,'(''LIN_D: state, vec, energy'',2i4,2f12.5)') istate,ivec,e(ivec),overlap_psi(ivec,istate)
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
          write(ounit,'(''LIN_D: max overlap ivec'',i4)') i_ortho_min

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
c elseif I do not optimize jastrow and or orbitals

          do istate=1,nstates
            call sort(nvec,overlap_psi(1,istate),index_overlap)
            i_overlap_max=index_overlap(nvec)
            write(ounit,'(''LIN_D: state, max overlap ivec'',2i4)') istate,i_overlap_max
            do i=1,nparm
              deltap(i+nparm*(istate-1))=evc(i,i_overlap_max)/anorm(i_overlap_max)
            enddo
c Save 5 additional vectors with large overlap
            do ivec=1,5
              if (nvec-ivec > 0) then
                idx_ivec=index_overlap(nvec-ivec)
                index_more(ivec,istate)=idx_ivec
                do i=1,nparm
                  deltap_more(i+nparm*(istate-1),ivec)=evc(i,idx_ivec)/anorm(idx_ivec)
                enddo
              endif
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
          ! NR why a BCAST in a loop ?
          call MPI_BCAST(deltap_more(1,ivec),nparm*nstates,MPI_REAL8,0,MPI_COMM_WORLD,ier)
        enddo

        call MPI_BCAST(index_more,5*nstates,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      endif
      return              ! deltap
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine amul(n,q,r)

      use jd_scratch, only: qr,rr
      use mpi
      use optwf_lin_dav_extra, only: h_psi_lin_d

      implicit none

      integer :: i, ier, n
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

      use jd_scratch, only: qr,rr
      use mpi
      use optwf_lin_dav_extra, only: s_psi_lin_d

      implicit none

      integer :: i, ier, n
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

      implicit none

      integer :: n
      complex*16 q(n)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end module
