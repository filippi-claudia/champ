      subroutine setup_optimization(nparm,mparmx,MWORK,lwork,h,h_sav,s,s_sav,work,eig_vec,add_diag,iter)

      use force, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use linear_norm, only: oav
      use optwf_contrl, only: ioptjas, ioptorb
      use optwf_corsam, only: energy, energy_err, force
      use optwf_parms, only: nparmd, nparmj
      use gradhess_all, only: MPARMALL

      use ci000, only: nciterm

      use method_opt, only: method

      implicit real*8(a-h,o-z)



      include 'optjas.h'
      include 'optci.h'
      include 'optorb.h'

      parameter(eps=1.d-12)

      dimension h(mparmx,*),s(mparmx,*)
      dimension h_sav(mparmx,*),s_sav(*)
      dimension eig(MPARMALL),eigi(MPARMALL),seig_inv(MPARMALL)
      dimension eig_vec(MPARMALL,*),hmod(MPARMALL,MPARMALL)
      dimension isort(MPARMALL)
      dimension work(*)

      save eig_min

      write(6,'(/,''Setup starting adiag'',/)') 
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
        do 3 i=1,nparm+is
          do 3 j=1,i
            idx=idx+1
            s_sav(idx)=s(i,j)
            h_sav(i,j)=h(i,j)
   3        h_sav(j,i)=h(j,i)
      endif

c Symmetrize the overlap
      do 4 i=1,nparm+is
        do 4 j=1,i
   4      s(j,i)=s(i,j)

c     do i=1,nparm+is
c       write(6,*) 'h =',(h(i,j),j=1,nparm+is)
c     enddo
c     do i=1,nparm+is
c       write(6,*) 's =',(s(i,j),j=1,nparm+is)
c     enddo

      if(add_diag.gt.0.and.iter.eq.1) then

      write(6,'(''Determine energy gap in super-CI hamiltonian at step 1'',/)')
      call regularize_geneig(nparm+is,mparmx,h,s,work,seig_inv,hmod)
      call solve_geneig(nparm+is,mparmx,hmod,s,seig_inv,work,eig,eigi,eig_vec)

      imag=0
      do 5 j=1,nparm+is
       if (eigi(j).ne.0.d0) imag=1
    5 continue
      if(imag.eq.1) write(6,'(''Warning: imaginary eigenvalues'')')

c Sort the eigenvalues
      call sort(nparm+is,eig,isort)

      ireal=0
      do 6 ii=1,nparm+is
        i=isort(ii)
        if(eigi(i).eq.0) then 
         ireal=ireal+1
         if(ireal.le.5)
     &   write(6,'(''eigenvalue '',i4,'' = '',f15.5,'' + i* '',f15.5)') i,eig(i),eigi(i)
        endif
    6 continue
      write(6,*)

c use semiorthogonal basis and compute minimum norm in direction orthogonal to psi0
      dmult=1.d0
      scale=.1d0
    9 de_range=scale*dabs(energy(1))
      emin=energy(1)-de_range
      emax=energy(1)+dmult*3*energy_err(1)

      ireal=0
      anorm_orth_min=1.d+99
      do 20 ii=1,nparm+is
        i=isort(ii)
        if(eigi(i).ne.0.or.eig(i).lt.emin) go to 20
        if(eig(i).gt.emax) go to 25
        ireal=ireal+1
        bot=1.d0
        do 10 j=1,nparmd
   10     bot=bot-eig_vec(j+nparmj+is,i)*oav(j+1)/eig_vec(1,i)
        idx=0
        anorm_orth=0.d0
        do 15 j=1,nparm+is
          do 15 k=1,j
            idx=idx+1
c 21/8 - ERROR?
c           if(j.ne.1.and.k.ne.1) then
            if(j.ne.1.or.k.ne.1) then
              anorm_orth=anorm_orth+eig_vec(j,i)*eig_vec(k,i)*s_sav(idx)
              if(j.ne.k) anorm_orth=anorm_orth+eig_vec(k,i)*eig_vec(j,i)*s_sav(idx)
            endif
   15   continue
        anorm_orth=sqrt(dabs(anorm_orth))/abs(bot*eig_vec(1,i))
        if(anorm_orth.lt.anorm_orth_min) then
          anorm_orth_min=anorm_orth
          i_ovr=i
          isort_ovr=ii
        endif
        if(ireal.le.5) then
          write(6,'(i4,'' eigenvalue '',i4,'' = '',f15.5,'' + i* '',f15.5)') ii,i,eig(i),eigi(i)
          write(6,'(4x,'' ortho dir  '',i4,'' = '',f15.5)') i,anorm_orth
        endif
   20 continue
   25 write(6,'(i4,'' real roots in energy range '',f15.5,'','',f15.5)') ireal,emin,emax
      if(ireal.eq.0) then
        dmult=dmult*2
        scale=scale*2
        goto 9
      endif
      write(6,'(''maximum overlap = '',i5,f15.5,'' + i * '',2f15.5)') i_ovr,eig(i_ovr),eigi(i_ovr),anorm_orth_min

      i0=i_ovr
      if(isort_ovr.lt.nparm+is) then
        i1=isort(isort_ovr+1)
        eig_min=eig(i1)-eig(i0)
      else
        eig_min=0
      endif
 
      write(6,'(''energy gap = '',f15.5)') eig_min
      endif

      endif

      call p2gtid('optwf:multiple_adiag',multiple_adiag,0,1)
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
      write(6,'(/,''starting adiag'',g12.4,/)') add_diag
       

      return
      end
c-----------------------------------------------------------------------
      subroutine regularize_geneig(n,mparmx,h,s,work,seig_valinv,hmod)

      use force, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use gradhess_all, only: MPARMALL
      
      implicit real*8(a-h,o-z)

      include 'optjas.h'
      include 'optci.h'
      include 'optorb.h'
      include 'numbas.h'

      parameter(eps=1.d-12)
      parameter(MWORK=50*MPARMALL)
      parameter(eps_eigval=1.d-14)
      
      dimension h(mparmx,*),s(mparmx,*)
      dimension seig_vals(MPARMALL),seig_valinv(*)
      dimension hmod(MPARMALL,*),work(*)

      call cpu_time(t0)

c call dsyev to determine lworks
      call dsyev('V','U',n,s,mparmx,seig_vals,work,-1,isdinfo)
      lworks=work(1)
c     write(6,*) 'S diag opt lwork=',lworks

c diagonalize s=S -> S_diag=U^T S U -> in output, s contains the unitary matrix U
      call dsyev('V','U',n,s,mparmx,seig_vals,work,lworks,isdinfo)
c     call cpu_time(t)
c     t_sdiag=t 
c     write(6,*) 'elapsed time for diagonalization:',t_sdiag-t0

      icut=1
      ineg=1
      write(6,'(''overlap matrix: Maximum eigenvalue, threshold: '',1p2e12.4)') seig_vals(n),eps_eigval*seig_vals(n)
      do 10 i=1,n
        if((ineg.eq.1).and.(seig_vals(i).gt.0.d0)) then
          ineg=0
          write(6,'(''first positive eigenvalue'',t41,i6)') i
        endif
        if((icut.eq.1).and.(seig_vals(i).gt.eps_eigval*seig_vals(n))) then
          icut=0
          write(6,'(''first eigenvalue larger than threshold'',t41,i6)') i
        endif
  10  continue
      do 20 i=1,n
        if(seig_vals(i)/seig_vals(n).gt.eps_eigval) then
          seig_valinv(i)=1.0d0/dsqrt(seig_vals(i))
         else 
          seig_valinv(i)=0.0d0
        endif
  20  continue

c I = A^T S A where A_ij= U_ij/sqrt(s_eigval(j))
c s in input contains U -> s is overwritten and contains A 
      do 30 i=1,n
        do 30 j=1,n
          s(i,j)=s(i,j)*seig_valinv(j)
  30  continue

c Compute work_mat=A^T*H*A
      do 35 i=1,n
        do 35 l=1,n
  35      hmod(l,i)=0.d0
       
      do 50 l=1,n
        do 40 m=1,n
          work(m)=0.d0
          do 40 k=1,n
  40        work(m)=work(m)+s(k,l)*h(k,m)
        do 50 i=1,n
          do 50 m=1,n
  50        hmod(l,i)=hmod(l,i)+work(m)*s(m,i)

      call cpu_time(t)
c     t_hmod=t
c     write(6,*) 'elapsed time to build Hmod:',t_hmod-t_sdiag

      return
      end

c-----------------------------------------------------------------------
      subroutine solve_geneig(n,mparmx,hmod,s,seig_valinv,work,eig,eigi,eig_vec)
      
      use force, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use gradhess_all, only: MPARMALL

      implicit real*8(a-h,o-z)

      include 'optjas.h'
      include 'optci.h'
      include 'optorb.h'
      include 'numbas.h'

      parameter(eps=1.d-12)
      parameter(MWORK=50*MPARMALL)
      
      dimension seig_valinv(*)
      dimension hmod(mparmx,*),s(mparmx,*)

      dimension eig_vec(MPARMALL,*),eig_vecl(MPARMALL,1)
      dimension work(*)

c s_fordiag: a copy of S for diagonalization. 
c hmod: the modified Hamiltonian matrix (in the end, S^-1*U*H*U^T)
c s: overlap matrix, h: hamiltonian, eigenvec: eigenvectors, 

      call cpu_time(t0)

c MISSING
c hmod+adiag/s_diag

c Determine ilwork
      call dgeev('N','V',n,hmod,MPARMALL,eig,eigi,eig_vecl,
     &        MPARMALL,eig_vec,MPARMALL,work,-1,isdinfo)
      ilwork=work(1)
c     write(6,*) 'isdinfo, optimal lwork=',isdinfo,ilwork

c Diagonalize
      call dgeev('N','V',n,hmod,MPARMALL,eig,eigi,eig_vecl,
     &        MPARMALL,eig_vec,MPARMALL,work,ilwork,isdinfo)
c     write(6,*) 'isdinfo=',isdinfo
c     call cpu_time(t)
c     t_hmdiag=t
c     write(6,*) 'elapsed time to diagonalize Hmod:',t-t0

      do 20 k=1,n
        do 10 i=1,n
          work(i)=0.d0
          do 10 j=1,n
  10        work(i)=work(i)+s(i,j)*eig_vec(j,k)
        do 20 i=1,n
  20      eig_vec(i,k)=work(i)
      
c     call cpu_time(t)
c     t_eigvec=t
c     write(6,*) 'elapsed time to get eigenvectors:',t_eigvec-t_hmdiag

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_dparm(nparm,mparmx,lwork,dparm,h,h_sav,s,s_sav,work,eig_vec,
     &                     add_diag,energy_sav,energy_err_sav)

      use force, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use csfs, only: ccsf, ncsf, nstates
      use dets, only: cdet
      use linear_norm, only: oav
      use optwf_contrl, only: ioptjas, ioptorb
      use optwf_parms, only: nparmd, nparmj
      use gradhess_all, only: MPARMALL

      use ci000, only: nciterm

      use method_opt, only: method

      implicit real*8(a-h,o-z)



      include 'optjas.h'
      include 'optci.h'
      include 'optorb.h'

      parameter(eps=1.d-12)

      dimension dparm(*)
      dimension h(mparmx,*),h_sav(mparmx,*)
      dimension s(mparmx,*),s_sav(*)
      dimension work(*)

      dimension eig(MPARMALL),eigi(MPARMALL),seig_inv(MPARMALL)
      dimension eig_vec(MPARMALL,*),hmod(MPARMALL,MPARMALL)
      dimension s_norm(MXCIMATDIM)
      dimension isort(MPARMALL)
      dimension cdelta(MPARMALL)
      dimension overlap(MXCITERM)
    
      if(method.eq.'linear') then

      nparmd=max(nciterm-1,0)

      is=1
      if(ioptjas.eq.0.and.ioptorb.eq.0) then 
	is=0
        idx=0
        do 105 i=1,nparm+is
          do 105 j=1,i
            idx=idx+1
  105       s_norm(idx)=s(i,j)
       else
        idx=0
        do 107 i=1,nparm+is
          do 107 j=1,i
            idx=idx+1
            s(i,j)=s_sav(idx)
  107       s(j,i)=s_sav(idx)
        do 108 i=1,nparm+is
          do 108 j=1,nparm+is
            h(i,j)=h_sav(i,j)
  108   continue
        do 109 i=2,nparm+is
          h(i,i)=h(i,i)+add_diag
  109   continue
      endif

c     do i=1,nparm+is
c     write(6,'(''h '',1000e25.15)') (h(i,j),j=1,nparm+is)
c     enddo
c     do i=1,nparm+is
c     write(6,'(''s '',1000e25.15)') (s(i,j),j=1,nparm+is)
c     enddo

      call regularize_geneig(nparm+is,mparmx,h,s,work,seig_inv,hmod)
      call solve_geneig(nparm+is,mparmx,hmod,s,seig_inv,work,eig,eigi,eig_vec)

      imag=0
      do 110 j=1,nparm+is
       if (eigi(j).ne.0.d0) imag=1
  110 continue

      if(imag.eq.1) write(6,'(''Warning: imaginary eigenvalues'')')

c Sort the eigenvalues
      call sort(nparm+is,eig,isort)

      i_min=isort(1)
      write(6,'(''generalized eigenvalue problem H c= E S c'')')
      write(6,'(''lowest eigenval = '',i5,f15.5,'' + i * '',f15.5)') i_min,eig(i_min),eigi(i_min)

c use semiorthogonal basis and compute minimum norm in direction orthogonal to psi0
      dmult=1.d0
      scale=.1d0

      no_real_found=0
  116 de_range=scale*dabs(energy_sav)
      emin=energy_sav-de_range
      emax=energy_sav+dmult*3*energy_err_sav

      ireal=0
      anorm_orth_min=1.d+99
      do 119 ii=1,nparm+is
        i=isort(ii)
        if(eigi(i).ne.0.or.eig(i).lt.emin) go to 119
        if(eig(i).gt.emax) go to 120
        ireal=ireal+1

        if(ioptjas.ne.0.or.ioptorb.ne.0) then

          bot=1.d0
          do 117 j=1,nparmd
  117        bot=bot-eig_vec(nparmj+is+j,i)*oav(j+1)/eig_vec(1,i)
          idx=0
          anorm_orth=0.d0
          do 118 j=1,nparm+is
            do 118 k=1,j
              idx=idx+1
c 21/8 - ERROR?
c             if(j.ne.1.and.k.ne.1) then
              if(j.ne.1.or.k.ne.1) then
                anorm_orth=anorm_orth+eig_vec(j,i)*eig_vec(k,i)*s_sav(idx)
                if(j.ne.k) anorm_orth=anorm_orth+eig_vec(k,i)*eig_vec(j,i)*s_sav(idx)
              endif
  118     continue
          anorm_orth=sqrt(dabs(anorm_orth))/abs(bot*eig_vec(1,i))
          if(anorm_orth.lt.anorm_orth_min) then
            anorm_orth_min=anorm_orth
            i_ovr=i
          endif
          if(ireal.le.5) then
            write(6,'(i4,'' eigenvalue '',i4,'' = '',f15.5,'' + i* '',f15.5)') ii,i,eig(i),eigi(i)
            write(6,'(4x,'' ortho dir  '',i4,'' = '',f15.5)') i,anorm_orth
          endif
        else
          if(ireal.le.5) write(6,'(i4,'' eigenvalue '',i4,'' = '',f15.5,'' + i* '',f15.5)') ii,i,eig(i),eigi(i)
        endif
  119 continue
  120 write(6,'(i4,'' real roots in energy range '',f15.5,'','',f15.5)') ireal,emin,emax
      if(ireal.eq.0) then
        dmult=dmult*2
        scale=scale*2

        no_real_found=no_real_found+1
        if(no_real_found.gt.nparmd) 
     &      call fatal_error('OPTWF: Cannot find real eigenvalues')
        goto 116
      endif

      if(ioptjas.eq.1.or.ioptorb.eq.1) then
        write(6,'(''maximum overlap = '',i5,f15.5,'' + i * '',2f15.5)') i_ovr,eig(i_ovr),eigi(i_ovr),anorm_orth_min

        if(i_ovr.ne.i_min) write(6,'(''Warning: max overlap not for min eigenvalue'')')

        i_good=i_ovr
        do 130 i=2,nparm+is
  130    cdelta(i-1)=eig_vec(i,i_good)/eig_vec(1,i_good)

        write(6,'(''dp  ='',1000f10.5)') (cdelta(i),i=1,nparm)

        bot=1
        do 145 i=1,nparmd
  145     bot=bot-cdelta(nparmj+i)*oav(i+1)

c minus sign because variation is subtracted when computing new parameters
        do 150 i=1,nparm
  150     cdelta(i)=-cdelta(i)/bot

        write(6,'(''dpn ='',1000f10.5)') (-cdelta(i),i=1,nparm)

        do 160 i=1,nparm
  160     dparm(i)=cdelta(i)

       else
        i0=0
        do 180 jj=1,nstates

          write(6,'(''State '',i4)') jj
          dnorm_jj=0
          idx=0

          if(ncsf.gt.0) then
            do 162 i=1,nparm
              do 162 k=1,i
                idx=idx+1
                dmul=1.d0
                if(i.ne.k) dmul=2.d0
  162            dnorm_jj=dnorm_jj+dmul*ccsf(i,jj,1)*ccsf(k,jj,1)*s_norm(idx)
          else
            do 163 i=1,nparm
              do 163 k=1,i
                idx=idx+1
                dmul=1.d0
                if(i.ne.k) dmul=2.d0
  163            dnorm_jj=dnorm_jj+dmul*cdet(i,jj,1)*cdet(k,jj,1)*s_norm(idx)
          endif
          dnorm_jj=1.d0/dsqrt(dnorm_jj)

          write(6,'(''state '',i4,'' norm '',1p1e12.5)') jj,dnorm_jj

          do 172 j=i0+1,nparm

            overlap(j)=0.d0

            jsort=isort(j)

            if(eig(jsort).ne.0.and.eigi(jsort).eq.0.d0) then
              idx=0
              if(ncsf.gt.0) then
                do 166 i=1,nparm
                  do 166 k=1,i
                    idx=idx+1
                    if(i.ne.k) then
                      overlap(j)=overlap(j)+(eig_vec(i,jsort)*ccsf(k,jj,1)+eig_vec(k,jsort)*ccsf(i,jj,1))*s_norm(idx)
                     else
                      overlap(j)=overlap(j)+eig_vec(i,jsort)*ccsf(k,jj,1)*s_norm(idx)
                    endif
  166           continue
               else
                do 167 i=1,nparm
                  do 167 k=1,i
                    idx=idx+1
                    if(i.ne.k) then
                      overlap(j)=overlap(j)+(eig_vec(i,jsort)*cdet(k,jj,1)+eig_vec(k,jsort)*cdet(i,jj,1))*s_norm(idx)
                     else
                      overlap(j)=overlap(j)+eig_vec(i,jsort)*cdet(k,jj,1)*s_norm(idx)
                    endif
  167           continue
              endif

              dnorm=0.d0
              idx=0
              do 170 i=1,nparm
                do 170 k=1,i
                  idx=idx+1
                  dmul=1.d0
                  if(i.ne.k) dmul=2.d0
  170             dnorm=dnorm+dmul*eig_vec(i,jsort)*eig_vec(k,jsort)*s_norm(idx)
              dnorm=1.d0/dsqrt(dnorm)

              overlap(j)=dabs(overlap(j))*dnorm*dnorm_jj

              write(6,'(i4,'' eigenvalue '',i4,'' = '',f15.5,'' + i* '',f15.5)') j,jsort,eig(jsort),eigi(jsort)
              write(6,'('' overlap state,eigenstate '',2i4,'' = '',f15.5)') jj,j,overlap(j)
  
            endif
  172     continue

          target_overlap=-1.d0
          do 175 j=i0+1,nparm
            if(overlap(j).gt.target_overlap) then
              i0=j
              target_overlap=overlap(j)
            endif
  175     continue

          do 178 i=1,nparm
  178       dparm(i+nparm*(jj-1))=eig_vec(i,isort(i0))*dnorm
          write(6,'(''state '',i4,'' norm'',1p1e12.5,'' overlap '',1p1e12.5)') jj,dnorm,overlap(i0)
          write(6,'(''pn  ='',1000f10.5)') (dparm(i+nparm*(jj-1)),i=1,nparm)

          if(nstates.gt.1.and.jj.ne.nstates.and.eig(isort(i0+1)).eq.0.d0) 
     &      call fatal_error('OPTWF: Overlap with state 1 for highest eigenvalue >0')
  180   continue
       endif

      endif

      end
c-----------------------------------------------------------------------
