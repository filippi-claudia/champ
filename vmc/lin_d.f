ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine h_psi_lin_d(ndim,nvec,psi,hpsi )
      implicit real*8 (a-h,o-z)

      include 'sr.h'
      include 'mstates.h'
 
      common /optwf_func/ omega,ifunc_omega

      dimension psi(MPARM,*),hpsi(MPARM,*)

      if(ifunc_omega.eq.0) then
        call h_psi_energymin(ndim,nvec,psi,hpsi )
       else
        call h_psi_omegamin(ndim,nvec,psi,hpsi )
      endif
      
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine s_psi_lin_d(ndim,nvec,psi,spsi )
      implicit real*8 (a-h,o-z)

      include 'sr.h'
      include 'mstates.h'
 
      common /optwf_func/ omega,ifunc_omega

      dimension psi(MPARM,*),spsi(MPARM,*)

      if(ifunc_omega.eq.0) then
        call s_psi_energymin(ndim,nvec,psi,spsi )
       else
        call s_psi_omegamin(ndim,nvec,psi,spsi )
      endif
      
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine select_ci_root(iroot)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      do 30 i=1,ndet
   30   cdet(i,1,1)=cdet(i,iroot,1)

      do 40 icsf=1,ncsf
   40   ccsf(icsf,1,iadiag)=ccsf(icsf,iroot,1)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine jdqz_driver( n, kmax, jmin, jmax, evc, eps,
     &                        e, e0, itype, notcnv, idav_iter , ipr )
      implicit real*8 (a-h,o-z)
      include 'sr.h'
c     parameter(lwork=10+6*MVEC+5*MVEC+3*MVEC)
      parameter(lwork=MPARM*100)
      dimension e(MVEC),evc(MPARM,MVEC),itype(MVEC)
      complex*16 alpha(MVEC),beta(MVEC),eivec(MPARM*MVEC),zwork(lwork),tmp(MPARM)
      complex*16 target,residu
      logical wanted

      wanted=.true.
      target=cmplx(e0,0.d0)

      method=1
      mxmv=100
      maxstep=100
      alock=eps
      iorder=0
      itestspace=3

      call JDQZ(alpha,beta,eivec,wanted,n,target,eps
     &         ,kmax,jmax,jmin,method,jmax,0,mxmv,maxstep,alock,iorder
     &         ,itestspace,zwork,lwork)

      write(6,'(''converged roots : '',i4)') kmax
      do j=1,kmax
        ish=n*(j-1)
        write(6,'(''norm : '',e11.4)') dznrm2( n, eivec(1+ish), 1 )
        call amul(n,eivec(1+ish),residu)
        call zscal(n,beta(j),residu,1)
        call bmul(n,eivec(1+ish),tmp)
        call zaxpy(n,-alpha(j),tmp,1,residu,1)
        write(6,'(''lambda('',i2,''): ('',1p,e11.4,'','',e11.4,'' )'')') j,alpha(j)/beta(j)
        write(6,'(a30,d13.6)') '||beta Ax - alpha Bx||:', dznrm2( n, residu, 1 )
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
