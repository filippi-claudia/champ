      module multideterminante_mod
      contains
      subroutine multideterminante(iel)

      use vmc_mod, only: norb_tot
      use vmc_mod, only: MEXCIT
      use csfs, only: nstates
      use dets, only: ndet
      use system, only: ndn, nup
      use multidet, only: irepcol_det, ireporb_det, ivirt, iwundet, kref, numrep_det, k_det, ndetiab, ndet_req
      use multidet, only: k_det2, ndetiab2, k_aux, ndetsingle
      use slatn, only: slmin
      use ycompactn, only: ymatn
      use coefs, only: norb
      use multimatn, only: aan, wfmatn
      use multislatern, only: detn, orbn
      use system, only: nelec
      use orbval, only: orb
      use multislater, only: detiab
      use precision_kinds, only: dp
      use contrl_file,    only: ounit
      use matinv_mod, only: matinv
      use multideterminant_mod, only: compute_ymat
      implicit none

      integer :: i, iab, iel, index_det, iorb
      integer :: irep, ish, istate, jj
      integer :: jorb, jrep, k, ndim, ndim2, kun, kw, kk
      integer :: nel
      real(dp) :: det, dum1
      real(dp), dimension(nelec, norb_tot, 3) :: gmat
      real(dp), dimension(MEXCIT**2, 3) :: gmatn
      real(dp), dimension(norb_tot, 3) :: b
      real(dp), dimension(3) :: ddx_mdet
      real(dp), dimension(norb_tot) :: orb_sav
      real(dp), dimension(ndet_req) :: ddetn




      if(ndet.eq.1) return

      iab=1
      nel=nup
      ish=0
      if(iel.gt.nup) then
        iab=2
        nel=ndn
        ish=nup
      endif

c temporarely copy orbn to orb
      do iorb=1,norb
        orb_sav(iorb)=orb(iel,iorb)
        orb(iel,iorb)=orbn(iorb)
      enddo

      do jrep=ivirt(iab),norb
        do irep=1,nel

          dum1=0.d0
          do i=1,nel
           dum1=dum1+slmin(irep+(i-1)*nel)*orb(i+ish,jrep)
          enddo
          aan(irep+nelec*(jrep-1))=dum1
        enddo
      enddo

c compute wave function
c     loop inequivalent determinants
c     loop over single exitations
      do k=1,ndetsingle(iab)
         
         jorb=ireporb_det(1,k,iab)
         iorb=irepcol_det(1,k,iab)
         wfmatn(k,1)=aan(iorb+nelec*(jorb-1))
         ddetn(k)=wfmatn(k,1)
         wfmatn(k,1)=1.0d0/wfmatn(k,1)
      enddo

c     loop over single exitations      
c         jorb=ireporb_det(1,1:ndetsingle(iab),iab)
c         iorb=irepcol_det(1,1Kndetsingle(iab),iab)
      
c      wfmatn(1:ndetsingle(iab),1)=aan(irepcol_det(1,1:ndetsingle(iab),iab)+nelec*(ireporb_det(1,1:ndetsingle(iab),iab)-1))
c      ddetn(1:ndetsingle(iab))=wfmatn(1:ndetsingle(iab),1)
c      wfmatn(1:ndetsingle(iab),1)=1.0d0/wfmatn(1:ndetsingle(iab),1)

      
c     loop over multiple exitations      
      do k=ndetsingle(iab)+1,ndetiab(iab)
     
         ndim=numrep_det(k,iab)
         ndim2=ndim*ndim
     
         jj=0
         do jrep=1,ndim
            jorb=ireporb_det(jrep,k,iab)
            do irep=1,ndim
               iorb=irepcol_det(irep,k,iab)
               jj=jj+1
               wfmatn(k,jj)=aan(iorb+nelec*(jorb-1))
            enddo
         enddo
         call matinv(wfmatn(k,1:ndim2),ndim,ddetn(k))
      enddo
      
      

c     Unrolling determinats different to kref
      detn=detn(kref)
      do kk=1,ndetiab2(iab)
         k=k_det2(kk,iab)
         kw=k_aux(kk,iab)  
         detn(k)=detn(k)*ddetn(kw)
c     print *, "k ",k,"detn(k) ",detn(k)
      enddo
c      k_det2(1:ndetiab2(iab),iab)
c      k_aux(1:ndetiab2(iab),iab)  
c      detn(k_det2(1:ndetiab2(iab),iab))=detn(k_det2(1:ndetiab2(iab),iab))*ddetn(k_aux(1:ndetiab2(iab),iab))
      
c     do istate=1,nstates
c        if(iab.eq.1) call compute_ymat(iab,detn,detiab(1,2),wfmatn,ymatn(1,1,istate),istate)
c     if(iab.eq.2) call compute_ymat(iab,detiab(1,1),detn,wfmatn,ymatn(1,1,istate),istate)
c     enddo
      
      if (iab.eq.1) then
         do istate=1,nstates
            call compute_ymat(iab,detn,detiab(1,2),wfmatn,ymatn(1,1,istate),istate)
         enddo
      else
         do istate=1,nstates
            call compute_ymat(iab,detiab(1,1),detn,wfmatn,ymatn(1,1,istate),istate)
         enddo
      endif
      
      do iorb=1,norb
        orb(iel,iorb)=orb_sav(iorb)
      enddo
      
      return
      end

c-----------------------------------------------------------------------
      subroutine multideterminante_grad(iel,b,norbs,detratio,slmi,aa,ymat,velocity)

      use precision_kinds, only: dp
      use vmc_mod, only: norb_tot
      use vmc_mod, only: nmat_dim
      use vmc_mod, only: MEXCIT
      use dets, only: ndet
      use system, only: ndn, nup
      use multidet, only: iactv, ivirt, kref
      use coefs, only: norb
      use dorb_m, only: iworbd
      use system, only: nelec

      implicit none

      integer :: iab, iel, iorb, irep, ish, norbs
      integer :: j, jel, jrep, k
      integer :: kk, nel
      real(dp) :: detratio, dum
      real(dp), dimension(nelec*norb_tot) :: aa
      real(dp), dimension(norb_tot, nelec) :: ymat
      real(dp), dimension(norbs, 3) :: b
      real(dp), dimension(nelec, norb_tot, 3) :: gmat
      real(dp), dimension(3) :: velocity
      real(dp), dimension(nmat_dim) :: slmi




      do k=1,3
        velocity(k)=0.d0
      enddo
      if(ndet.eq.1) return

      if(iel.le.nup) then
        iab=1
        nel=nup
        ish=0
       else
        iab=2
        nel=ndn
        ish=nup
      endif

      jel=iel-ish


      do kk=1,3

        do jrep=ivirt(iab),norb
          dum=0
          do j=1,nel
             dum=dum+b(iworbd(j+ish,kref),kk)*aa(j+nelec*(jrep-1))
          enddo
          dum=b(jrep,kk)-dum
          do irep=iactv(iab),nel
            gmat(irep,jrep,kk)=dum*slmi(irep+(jel-1)*nel)
          enddo
        enddo

      enddo

c     if(iab.eq.2) write(ounit,*) 'gmat ',(((gmat(irep,jrep,kk),irep=iactv(iab),nel),jrep=ivirt(iab),norb),kk=1,3)

      do kk=1,3
        dum=0
        do jrep=ivirt(iab),norb
          do irep=iactv(iab),nel
            dum=dum+ymat(jrep,irep)*gmat(irep,jrep,kk)
          enddo
        enddo
        velocity(kk)=dum*detratio
      enddo

c     if(iab.eq.2) write(ounit,*) 'ymat ',((ymat(jrep,irep),irep=iactv(iab),nel),jrep=ivirt(iab),norb)
      return
      end
c-----------------------------------------------------------------------
      end module
