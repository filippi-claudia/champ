      module optwf_lin_matrix
      use error, only: fatal_error
      use optwf_lib, only: sort
      interface !LAPACK interface
      SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
     $                  LDVR, WORK, LWORK, INFO )
*  -- LAPACK driver routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
      DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WI( * ), WORK( * ), WR( * )
      END SUBROUTINE
      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
*  -- LAPACK driver routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
      END SUBROUTINE
      end interface
      contains
      subroutine setup_optimization(nparm,mparmx,MWORK,lwork,h,h_sav,s,s_sav,work,eig_vec,add_diag,iter)

      use linear_norm, only: ci_oav
      use optwf_control, only: ioptjas, ioptorb, multiple_adiag
      use optwf_corsam, only: energy, energy_err
      use optwf_parms, only: nparmd, nparmj
      use gradhess_all, only: nparmall
      use ci000, only: nciterm
      use optwf_control, only: method
      use precision_kinds, only: dp
      use contrl_file,    only: ounit

      implicit none

      integer :: i, i0, i1, i_ovr, idx
      integer :: ii, imag, ireal, is
      integer :: isort_ovr, iter, j, k
      integer :: lwork, mparmx, nparm
      integer, dimension(nparmall) :: isort
      integer  :: MWORK
      real(dp) :: add_diag, anorm_orth, anorm_orth_min, bot
      real(dp) :: de_range, dmult, eig_min
      real(dp) :: emax, emin, scale
      real(dp), dimension(mparmx,*) :: h
      real(dp), dimension(mparmx,*) :: s
      real(dp), dimension(mparmx,*) :: h_sav
      real(dp), dimension(*) :: s_sav
      real(dp), dimension(nparmall) :: eig
      real(dp), dimension(nparmall) :: eigi
      real(dp), dimension(nparmall) :: seig_inv
      real(dp), dimension(nparmall,*) :: eig_vec
      real(dp), dimension(nparmall,nparmall) :: hmod
      real(dp), dimension(*) :: work
      real(dp), parameter :: eps = 1.d-12



      save eig_min

      write(ounit,'(/,''Setup starting adiag'',/)')
      eig_min=0

      if(method.eq.'linear') then

      nparmd=max(nciterm-1,0)

c If we optimize jas, orb, or both, we solve secular super-CI equation on basis of psi_0,dpsi_0 for a total dimension of nparm+1
      is=1
      if(ioptjas.eq.0.and.ioptorb.eq.0) then
c If CI only, we solve secular equation on basis of J*CSF
        is=0
       else
c Save the overlap and hamiltonian
        idx=0
        do i=1,nparm+is
          do j=1,i
            idx=idx+1
            s_sav(idx)=s(i,j)
            h_sav(i,j)=h(i,j)
            h_sav(j,i)=h(j,i)
          enddo
        enddo
      endif

c Symmetrize the overlap
      do i=1,nparm+is
        do j=1,i
          s(j,i)=s(i,j)
        enddo
      enddo

      ! do i=1,nparm+is
      !   write(ounit,*) 'h =',(h(i,j),j=1,nparm+is)
      ! enddo
      ! do i=1,nparm+is
      !   write(ounit,*) 's =',(s(i,j),j=1,nparm+is)
      ! enddo


      if(add_diag.gt.0.and.iter.eq.1) then

      write(ounit,'(''Determine energy gap in super-CI hamiltonian at step 1'',/)')
      call regularize_geneig(nparm+is,mparmx,h,s,work,seig_inv,hmod)
      call solve_geneig(nparm+is,mparmx,hmod,s,seig_inv,work,eig,eigi,eig_vec)

      imag=0
      do j=1,nparm+is
       if (eigi(j).ne.0.d0) imag=1
      enddo
      if(imag.eq.1) write(ounit,'(''Warning: imaginary eigenvalues'')')

c Sort the eigenvalues
      call sort(nparm+is,eig,isort)

      ireal=0
      do ii=1,nparm+is
        i=isort(ii)
        if(eigi(i).eq.0) then
         ireal=ireal+1
         if(ireal.le.5)
     &   write(ounit,'(''eigenvalue '',i4,'' = '',f15.5,'' + i* '',f15.5)') i,eig(i),eigi(i)
        endif
      enddo
      write(ounit,*)

c use semiorthogonal basis and compute minimum norm in direction orthogonal to psi0
      dmult=1.d0
      scale=.1d0
    9 de_range=scale*dabs(energy(1))
      emin=energy(1)-de_range
      emax=energy(1)+dmult*3*energy_err(1)

      ireal=0
      anorm_orth_min=1.d+99
      do ii=1,nparm+is
        i=isort(ii)
        if(eigi(i).ne.0.or.eig(i).lt.emin) go to 20
        if(eig(i).gt.emax) go to 25
        ireal=ireal+1
        bot=1.d0
        do j=1,nparmd
          bot=bot-eig_vec(j+nparmj+is,i)*ci_oav(j+1)/eig_vec(1,i)
        enddo
        idx=0
        anorm_orth=0.d0
        do j=1,nparm+is
          do k=1,j
            idx=idx+1
c 21/8 - ERROR?
c           if(j.ne.1.and.k.ne.1) then
            if(j.ne.1.or.k.ne.1) then
              anorm_orth=anorm_orth+eig_vec(j,i)*eig_vec(k,i)*s_sav(idx)
              if(j.ne.k) anorm_orth=anorm_orth+eig_vec(k,i)*eig_vec(j,i)*s_sav(idx)
            endif
          enddo
        enddo
        anorm_orth=sqrt(dabs(anorm_orth))/abs(bot*eig_vec(1,i))
        if(anorm_orth.lt.anorm_orth_min) then
          anorm_orth_min=anorm_orth
          i_ovr=i
          isort_ovr=ii
        endif
        if(ireal.le.5) then
          write(ounit,'(i4,'' eigenvalue '',i4,'' = '',f15.5,'' + i* '',f15.5)') ii,i,eig(i),eigi(i)
          write(ounit,'(4x,'' ortho dir  '',i4,'' = '',f15.5)') i,anorm_orth
        endif
   20 continue
      enddo
   25 write(ounit,'(i4,'' real roots in energy range '',f15.5,'','',f15.5)') ireal,emin,emax
      if(ireal.eq.0) then
        dmult=dmult*2
        scale=scale*2
        goto 9
      endif
      write(ounit,'(''maximum overlap = '',i5,f15.5,'' + i * '',2f15.5)') i_ovr,eig(i_ovr),eigi(i_ovr),anorm_orth_min

      i0=i_ovr
      if(isort_ovr.lt.nparm+is) then
        i1=isort(isort_ovr+1)
        eig_min=eig(i1)-eig(i0)
      else
        eig_min=0
      endif

      write(ounit,'(''energy gap = '',f15.5)') eig_min
      endif

      endif

c Set add_diag
      if(add_diag.gt.0.or.multiple_adiag.eq.1) then
        add_diag=max(add_diag,1.d-6)
        add_diag=max(add_diag,1.d-2*eig_min)
c       add_diag=max(add_diag,1.d-2*abs(eig_min))
       else
c       add_diag=0
        add_diag=dabs(add_diag)
      endif
      if(ioptjas.eq.0.and.ioptorb.eq.0) add_diag=0
      write(ounit,'(/,''starting adiag'',g12.4,/)') add_diag


      return
      end
c-----------------------------------------------------------------------
      subroutine regularize_geneig(n,mparmx,h,s,work,seig_valinv,hmod)

      use gradhess_all, only: nparmall
      use precision_kinds, only: dp
      use contrl_file,    only: ounit
      implicit none

      integer :: i, icut, ineg, isdinfo, j
      integer :: k, l, lworks, m
      integer :: mparmx, n
      real(dp) :: t, t0
      real(dp), dimension(mparmx,*) :: h
      real(dp), dimension(mparmx,*) :: s
      real(dp), dimension(nparmall) :: seig_vals
      real(dp), dimension(*) :: seig_valinv
      real(dp), dimension(nparmall,*) :: hmod
      real(dp), dimension(*) :: work
      real(dp), parameter :: eps = 1.d-12
      real(dp), parameter :: eps_eigval = 1.d-14

      ! parameter(MWORK=50*nparmall)

      call cpu_time(t0)
c call dsyev to determine lworks
      call dsyev('V','U',n,s,mparmx,seig_vals,work,-1,isdinfo)
      lworks=work(1)

c     write(ounit,*) 'S diag opt lwork=',lworks

c diagonalize s=S -> S_diag=U^T S U -> in output, s contains the unitary matrix U
      call dsyev('V','U',n,s,mparmx,seig_vals,work,lworks,isdinfo)
      if (isdinfo.gt.0) call fatal_error('Eigenvalue issues in regularize_geneig')

c     call cpu_time(t)
c     t_sdiag=t
c     write(ounit,*) 'elapsed time for diagonalization:',t_sdiag-t0

      icut=1
      ineg=1
      write(ounit,'(''overlap matrix: Maximum eigenvalue, threshold: '',1p2e12.4)') seig_vals(n),eps_eigval*seig_vals(n)
      do i=1,n
        if((ineg.eq.1).and.(seig_vals(i).gt.0.d0)) then
          ineg=0
          write(ounit,'(''first positive eigenvalue'',t41,i6)') i
        endif
        if((icut.eq.1).and.(seig_vals(i).gt.eps_eigval*seig_vals(n))) then
          icut=0
          write(ounit,'(''first eigenvalue larger than threshold'',t41,i6)') i
        endif
      enddo
      do i=1,n
        if(seig_vals(i)/seig_vals(n).gt.eps_eigval) then
          seig_valinv(i)=1.0d0/dsqrt(seig_vals(i))
         else
          seig_valinv(i)=0.0d0
        endif
  20  continue
      enddo

c I = A^T S A where A_ij= U_ij/sqrt(s_eigval(j))
c s in input contains U -> s is overwritten and contains A
      do i=1,n
        do j=1,n
          s(i,j)=s(i,j)*seig_valinv(j)
        enddo
      enddo

c Compute work_mat=A^T*H*A
      do i=1,n
        do l=1,n
          hmod(l,i)=0.d0
        enddo
      enddo

      do l=1,n
        do m=1,n
          work(m)=0.d0
          do k=1,n
            work(m)=work(m)+s(k,l)*h(k,m)
          enddo
        enddo
        do i=1,n
          do m=1,n
            hmod(l,i)=hmod(l,i)+work(m)*s(m,i)
          enddo
        enddo
      enddo

      call cpu_time(t)
c     t_hmod=t
c     write(ounit,*) 'elapsed time to build Hmod:',t_hmod-t_sdiag

      return
      end

c-----------------------------------------------------------------------
      subroutine solve_geneig(n,mparmx,hmod,s,seig_valinv,work,eig,eigi,eig_vec)

      use gradhess_all, only: nparmall
      use precision_kinds, only: dp

      implicit none

      integer :: i, ilwork, isdinfo, j, k
      integer :: mparmx, n
      real(dp) :: t0
      real(dp), dimension(*) :: seig_valinv
      real(dp), dimension(mparmx,*) :: hmod
      real(dp), dimension(mparmx,*) :: s
      real(dp), dimension(nparmall,*) :: eig_vec
      real(dp), dimension(nparmall,nparmall) :: eig_vecl
      real(dp), dimension(*) :: work
      real(dp), dimension(nparmall) :: eig
      real(dp), dimension(nparmall) :: eigi
      real(dp), parameter :: eps = 1.d-12
      ! parameter(MWORK=50*nparmall)
      ! dimension eig_vecl(nparmall,1)

c s_fordiag: a copy of S for diagonalization.
c hmod: the modified Hamiltonian matrix (in the end, S^-1*U*H*U^T)
c s: overlap matrix, h: hamiltonian, eigenvec: eigenvectors,

      call cpu_time(t0)


c MISSING
c hmod+adiag/s_diag

c Determine ilwork
      call dgeev('N','V',n,hmod,nparmall,eig,eigi,eig_vecl,
     &        nparmall,eig_vec,nparmall,work,-1,isdinfo)
      ilwork=work(1)
c     write(ounit,*) 'isdinfo, optimal lwork=',isdinfo,ilwork

c Diagonalize
      call dgeev('N','V',n,hmod,nparmall,eig,eigi,eig_vecl,
     &        nparmall,eig_vec,nparmall,work,ilwork,isdinfo)
c     write(ounit,*) 'isdinfo=',isdinfo
c     call cpu_time(t)
c     t_hmdiag=t
c     write(ounit,*) 'elapsed time to diagonalize Hmod:',t-t0

      do k=1,n
        do i=1,n
          work(i)=0.d0
          do j=1,n
            work(i)=work(i)+s(i,j)*eig_vec(j,k)
          enddo
        enddo
        do i=1,n
  20      eig_vec(i,k)=work(i)
        enddo
      enddo

c     call cpu_time(t)
c     t_eigvec=t
c     write(ounit,*) 'elapsed time to get eigenvectors:',t_eigvec-t_hmdiag

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_dparm(nparm,mparmx,lwork,dparm,h,h_sav,s,s_sav,work,eig_vec,
     &                     add_diag,energy_sav,energy_err_sav)

      use optci, only: mxciterm, ncimatdim
      use csfs, only: ccsf, ncsf, nstates
      use slater, only: cdet
      use linear_norm, only: ci_oav
      use optwf_control, only: ioptjas, ioptorb
      use optwf_parms, only: nparmd, nparmj
      use gradhess_all, only: nparmall
      use ci000, only: nciterm
      use optwf_control, only: method
      use precision_kinds, only: dp
      use contrl_file,    only: ounit

      implicit none

      integer :: i, i0, i_good, i_min, i_ovr
      integer :: idx, ii, imag, ireal
      integer :: is, j, jj, jsort
      integer :: k, lwork, mparmx, no_real_found
      integer :: nparm
      integer, dimension(nparmall) :: isort
      real(dp) :: add_diag, anorm_orth, anorm_orth_min, bot
      real(dp) :: de_range, dmul, dmult, dnorm
      real(dp) :: dnorm_jj, emax, emin, energy_err_sav
      real(dp) :: energy_sav, scale, target_overlap
      real(dp), dimension(*) :: dparm
      real(dp), dimension(mparmx,*) :: h
      real(dp), dimension(mparmx,*) :: h_sav
      real(dp), dimension(mparmx,*) :: s
      real(dp), dimension(*) :: s_sav
      real(dp), dimension(*) :: work
      real(dp), dimension(nparmall) :: eig
      real(dp), dimension(nparmall) :: eigi
      real(dp), dimension(nparmall) :: seig_inv
      real(dp), dimension(nparmall,*) :: eig_vec
      real(dp), dimension(nparmall,nparmall) :: hmod
      real(dp), dimension(ncimatdim) :: s_norm
      real(dp), dimension(nparmall) :: cdelta
      real(dp), dimension(mxciterm) :: overlap
      real(dp), parameter :: eps = 1.d-12

      if(method.eq.'linear') then

      nparmd=max(nciterm-1,0)

      is=1
      if(ioptjas.eq.0.and.ioptorb.eq.0) then
	is=0
        idx=0
        do i=1,nparm+is
          do j=1,i
            idx=idx+1
            s_norm(idx)=s(i,j)
          enddo
        enddo
       else
        idx=0
        do i=1,nparm+is
          do j=1,i
            idx=idx+1
            s(i,j)=s_sav(idx)
            s(j,i)=s_sav(idx)
          enddo
        enddo
        do i=1,nparm+is
          do j=1,nparm+is
            h(i,j)=h_sav(i,j)
          enddo
        enddo
        do i=2,nparm+is
          h(i,i)=h(i,i)+add_diag
        enddo
      endif

c     do i=1,nparm+is
c     write(ounit,'(''h '',1000e25.15)') (h(i,j),j=1,nparm+is)
c     enddo
c     do i=1,nparm+is
c     write(ounit,'(''s '',1000e25.15)') (s(i,j),j=1,nparm+is)
c     enddo

      call regularize_geneig(nparm+is,mparmx,h,s,work,seig_inv,hmod)
      call solve_geneig(nparm+is,mparmx,hmod,s,seig_inv,work,eig,eigi,eig_vec)

      imag=0
      do j=1,nparm+is
       if (eigi(j).ne.0.d0) imag=1
      enddo

      if(imag.eq.1) write(ounit,'(''Warning: imaginary eigenvalues'')')

c Sort the eigenvalues
      call sort(nparm+is,eig,isort)

      i_min=isort(1)
      write(ounit,'(''generalized eigenvalue problem H c= E S c'')')
      write(ounit,'(''lowest eigenval = '',i5,f15.5,'' + i * '',f15.5)') i_min,eig(i_min),eigi(i_min)

c use semiorthogonal basis and compute minimum norm in direction orthogonal to psi0
      dmult=1.d0
      scale=.1d0

      no_real_found=0
  116 de_range=scale*dabs(energy_sav)
      emin=energy_sav-de_range
      emax=energy_sav+dmult*3*energy_err_sav

      ireal=0
      anorm_orth_min=1.d+99
      do ii=1,nparm+is
        i=isort(ii)
        if(eigi(i).ne.0.or.eig(i).lt.emin) go to 119
        if(eig(i).gt.emax) go to 120
        ireal=ireal+1

        if(ioptjas.ne.0.or.ioptorb.ne.0) then

          bot=1.d0
          do j=1,nparmd
             bot=bot-eig_vec(nparmj+is+j,i)*ci_oav(j+1)/eig_vec(1,i)
          enddo
          idx=0
          anorm_orth=0.d0
          do j=1,nparm+is
            do k=1,j
              idx=idx+1
c 21/8 - ERROR?
c             if(j.ne.1.and.k.ne.1) then
              if(j.ne.1.or.k.ne.1) then
                anorm_orth=anorm_orth+eig_vec(j,i)*eig_vec(k,i)*s_sav(idx)
                if(j.ne.k) anorm_orth=anorm_orth+eig_vec(k,i)*eig_vec(j,i)*s_sav(idx)
              endif
            enddo
          enddo
          anorm_orth=sqrt(dabs(anorm_orth))/abs(bot*eig_vec(1,i))
          if(anorm_orth.lt.anorm_orth_min) then
            anorm_orth_min=anorm_orth
            i_ovr=i
          endif
          if(ireal.le.5) then
            write(ounit,'(i4,'' eigenvalue '',i4,'' = '',f15.5,'' + i* '',f15.5)') ii,i,eig(i),eigi(i)
            write(ounit,'(4x,'' ortho dir  '',i4,'' = '',f15.5)') i,anorm_orth
          endif
        else
          if(ireal.le.5) write(ounit,'(i4,'' eigenvalue '',i4,'' = '',f15.5,'' + i* '',f15.5)') ii,i,eig(i),eigi(i)
        endif
  119 continue
      enddo
  120 write(ounit,'(i4,'' real roots in energy range '',f15.5,'','',f15.5)') ireal,emin,emax
      if(ireal.eq.0) then
        dmult=dmult*2
        scale=scale*2

        no_real_found=no_real_found+1
        if(no_real_found.gt.nparmd)
     &      call fatal_error('OPTWF: Cannot find real eigenvalues')
        goto 116
      endif

      if(ioptjas.eq.1.or.ioptorb.eq.1) then
        write(ounit,'(''maximum overlap = '',i5,f15.5,'' + i * '',2f15.5)') i_ovr,eig(i_ovr),eigi(i_ovr),anorm_orth_min

        if(i_ovr.ne.i_min) write(ounit,'(''Warning: max overlap not for min eigenvalue'')')

        i_good=i_ovr
        do i=2,nparm+is
         cdelta(i-1)=eig_vec(i,i_good)/eig_vec(1,i_good)
        enddo

        write(ounit,'(''dp  ='',1000f10.5)') (cdelta(i),i=1,nparm)

        bot=1
        do i=1,nparmd
          bot=bot-cdelta(nparmj+i)*ci_oav(i+1)
        enddo

c minus sign because variation is subtracted when computing new parameters
        do i=1,nparm
          cdelta(i)=-cdelta(i)/bot
        enddo

        write(ounit,'(''dpn ='',1000f10.5)') (-cdelta(i),i=1,nparm)

        do i=1,nparm
          dparm(i)=cdelta(i)
        enddo

       else
        i0=0
        do jj=1,nstates

          write(ounit,'(''State '',i4)') jj
          dnorm_jj=0
          idx=0

          if(ncsf.gt.0) then
            do i=1,nparm
              do k=1,i
                idx=idx+1
                dmul=1.d0
                if(i.ne.k) dmul=2.d0
                 dnorm_jj=dnorm_jj+dmul*ccsf(i,jj,1)*ccsf(k,jj,1)*s_norm(idx)
              enddo
            enddo
          else
            do i=1,nparm
              do k=1,i
                idx=idx+1
                dmul=1.d0
                if(i.ne.k) dmul=2.d0
                 dnorm_jj=dnorm_jj+dmul*cdet(i,jj,1)*cdet(k,jj,1)*s_norm(idx)
              enddo
            enddo
          endif
          dnorm_jj=1.d0/dsqrt(dnorm_jj)

          write(ounit,'(''state '',i4,'' norm '',1p1e12.5)') jj,dnorm_jj

          do j=i0+1,nparm

            overlap(j)=0.d0

            jsort=isort(j)

            if(eig(jsort).ne.0.and.eigi(jsort).eq.0.d0) then
              idx=0
              if(ncsf.gt.0) then
                do i=1,nparm
                  do k=1,i
                    idx=idx+1
                    if(i.ne.k) then
                      overlap(j)=overlap(j)+(eig_vec(i,jsort)*ccsf(k,jj,1)+eig_vec(k,jsort)*ccsf(i,jj,1))*s_norm(idx)
                     else
                      overlap(j)=overlap(j)+eig_vec(i,jsort)*ccsf(k,jj,1)*s_norm(idx)
                    endif
                  enddo
                enddo
               else
                do i=1,nparm
                  do k=1,i
                    idx=idx+1
                    if(i.ne.k) then
                      overlap(j)=overlap(j)+(eig_vec(i,jsort)*cdet(k,jj,1)+eig_vec(k,jsort)*cdet(i,jj,1))*s_norm(idx)
                     else
                      overlap(j)=overlap(j)+eig_vec(i,jsort)*cdet(k,jj,1)*s_norm(idx)
                    endif
                  enddo
                enddo
              endif

              dnorm=0.d0
              idx=0
              do i=1,nparm
                do k=1,i
                  idx=idx+1
                  dmul=1.d0
                  if(i.ne.k) dmul=2.d0
                  dnorm=dnorm+dmul*eig_vec(i,jsort)*eig_vec(k,jsort)*s_norm(idx)
                enddo
              enddo
              dnorm=1.d0/dsqrt(dnorm)

              overlap(j)=dabs(overlap(j))*dnorm*dnorm_jj

              write(ounit,'(i4,'' eigenvalue '',i4,'' = '',f15.5,'' + i* '',f15.5)') j,jsort,eig(jsort),eigi(jsort)
              write(ounit,'('' overlap state,eigenstate '',2i4,'' = '',f15.5)') jj,j,overlap(j)

            endif
          enddo

          target_overlap=-1.d0
          do j=i0+1,nparm
            if(overlap(j).gt.target_overlap) then
              i0=j
              target_overlap=overlap(j)
            endif
          enddo

          do i=1,nparm
            dparm(i+nparm*(jj-1))=eig_vec(i,isort(i0))*dnorm
          enddo
          write(ounit,'(''state '',i4,'' norm'',1p1e12.5,'' overlap '',1p1e12.5)') jj,dnorm,overlap(i0)
          write(ounit,'(''pn  ='',1000f10.5)') (dparm(i+nparm*(jj-1)),i=1,nparm)

          if(nstates.gt.1.and.jj.ne.nstates.and.eig(isort(i0+1)).eq.0.d0)
     &      call fatal_error('OPTWF: Overlap with state 1 for highest eigenvalue >0')
        enddo
       endif

      endif

      end
c-----------------------------------------------------------------------
      end module
