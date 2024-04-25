module multideterminante_mod
contains
      subroutine multideterminante(iel)

      use contrl_file,    only: ounit
      use csfs, only: nstates
      use matinv_mod, only: matinv
      use multidet, only: irepcol_det, ireporb_det, ivirt, numrep_det, k_det, ndetiab, ndet_req
      use multidet, only: k_aux2, k_det2, ndetiab2, ndetsingle, ndetdouble
      use multideterminant_mod, only: compute_ymat
      use multimatn, only: aan, wfmatn
      use multislater, only: detiab
      use multislatern, only: detn, orbn
      use orbval, only: orb
      use precision_kinds, only: dp
      use slater, only: iwundet, kref, ndet, norb
      use slatn, only: slmin
      use system, only: ndn, nup, nelec
      use vmc_mod, only: norb_tot, MEXCIT, stoo
      use ycompactn, only: ymatn
      implicit none

      integer :: i, iab, iel, index_det, iorb
      integer :: irep, ish, istate, jj
      integer :: jorb, jrep, k, ndim, ndim2, kun, kw, kk, kcum
      integer :: nel, o
      real(dp) :: det, dum1, deti, auxdet
      real(dp), dimension(nelec, norb_tot, 3) :: gmat
      real(dp), dimension(MEXCIT**2, 3) :: gmatn
      real(dp), dimension(norb_tot, 3) :: b
      real(dp), dimension(3) :: ddx_mdet
      real(dp), dimension(norb_tot) :: orb_sav
      real(dp), dimension(ndet_req) :: ddetn

      if(ndet.eq.1) return

      do istate=1,nstates
        o=stoo(istate)
        iab=1
        nel=nup
        ish=0
        if(iel.gt.nup) then
          iab=2
          nel=ndn
          ish=nup
        endif

! temporarely copy orbn to orb
        do iorb=1,norb
          orb_sav(iorb)=orb(iel,iorb,o)
          orb(iel,iorb,o)=orbn(iorb,o)
        enddo

        do jrep=ivirt(iab),norb
          do irep=1,nel

            dum1=0.d0
            do i=1,nel
             dum1=dum1+slmin(irep+(i-1)*nel,o)*orb(i+ish,jrep,o)
            enddo
            aan(irep,jrep,o)=dum1
          enddo
        enddo

! compute wave function
! loop inequivalent determinants
!
! loop over single exitations
        if(ndetsingle(iab).ge.1) then
          do k=1,ndetsingle(iab)
            jorb=ireporb_det(1,k,iab)
            iorb=irepcol_det(1,k,iab)
            wfmatn(k,1,o)=aan(iorb,jorb,o)
            ddetn(k)=wfmatn(k,1,o)
            wfmatn(k,1,o)=1.0d0/wfmatn(k,1,o)
          enddo
        endif

        kcum=ndetsingle(iab)+ndetdouble(iab)

! loop over double exitations
        if(ndetdouble(iab).ge.1)then
          do k=ndetsingle(iab)+1,kcum
            jorb=ireporb_det(1,k,iab)
            iorb=irepcol_det(1,k,iab)
            wfmatn(k,1,o)=aan(iorb,jorb,o)
            iorb=irepcol_det(2,k,iab)
            wfmatn(k,2,o)=aan(iorb,jorb,o)
            jorb=ireporb_det(2,k,iab)
            iorb=irepcol_det(1,k,iab)
            wfmatn(k,3,o)=aan(iorb,jorb,o)
            iorb=irepcol_det(2,k,iab)
            wfmatn(k,4,o)=aan(iorb,jorb,o)
            
            ddetn(k)=wfmatn(k,1,o)*wfmatn(k,4,o)-wfmatn(k,3,o)*wfmatn(k,2,o)
            deti=1.d0/ddetn(k)
            auxdet=wfmatn(k,1,o)
            wfmatn(k,1,o)=wfmatn(k,4,o)*deti
            wfmatn(k,2,o)=-wfmatn(k,2,o)*deti
            wfmatn(k,3,o)=-wfmatn(k,3,o)*deti
            wfmatn(k,4,o)=auxdet*deti
                  
          enddo
        endif

! loop over multiple exitations      
        if(kcum.lt.ndetiab(iab))then
          do k=kcum+1,ndetiab(iab)

            ndim=numrep_det(k,iab)
            ndim2=ndim*ndim

            jj=0
            do jrep=1,ndim
              jorb=ireporb_det(jrep,k,iab)
              do irep=1,ndim
                iorb=irepcol_det(irep,k,iab)
                jj=jj+1
                wfmatn(k,jj,o)=aan(iorb,jorb,o)
              enddo
            enddo
            call matinv(wfmatn(k,1:ndim2,o),ndim,ddetn(k))
          enddo
        endif 

! unrolling determinats different from kref
        detn(:,o)=detn(kref,o)
        do kk=1,ndetiab2(iab)
          k=k_det2(kk,iab)
          kw=k_aux2(kk,iab)  
          detn(k,o)=detn(k,o)*ddetn(kw)
        enddo

        if (iab.eq.1) then
          call compute_ymat(iab,detn(1,o),detiab(1,2,o),wfmatn(1,1,o),ymatn(1,1,istate),istate)
        else
          call compute_ymat(iab,detiab(1,1,o),detn(1,o),wfmatn(1,1,o),ymatn(1,1,istate),istate)
        endif
        do iorb=1,norb
          orb(iel,iorb,o)=orb_sav(iorb)
        enddo

      enddo

      return
      end

!-----------------------------------------------------------------------
      subroutine multideterminante_grad(iel,b,norbs,detratio,slmi,aa,ymat,velocity)

      use dorb_m,  only: iworbd
      use multidet, only: iactv,ivirt
      use precision_kinds, only: dp
      use vmc_mod, only: norb_tot, nmat_dim, MEXCIT
      use slater, only: ndet
      use system, only: ndn, nup, nelec
      use multidet, only: iactv, ivirt
      use slater, only: kref
      use slater, only: norb
      use dorb_m, only: iworbd
      use contrl_file, only: ounit

      implicit none

      integer :: iab, iel, iorb, irep, ish, norbs
      integer :: j, jel, jrep, k
      integer :: kk, nel
      real(dp) :: detratio, dum
      real(dp), dimension(nelec,norb_tot) :: aa
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
          dum=0.0d0
          do j=1,nel
             dum=dum+b(iworbd(j+ish,kref),kk)*aa(j,jrep)
          enddo
          dum=b(jrep,kk)-dum
          do irep=iactv(iab),nel
            gmat(irep,jrep,kk)=dum*slmi(irep+(jel-1)*nel)
          enddo
        enddo

      enddo

!     if(iab.eq.2) write(ounit,*) 'gmat ',(((gmat(irep,jrep,kk),irep=iactv(iab),nel),jrep=ivirt(iab),norb),kk=1,3)

      do kk=1,3
        dum=0
        do jrep=ivirt(iab),norb
          do irep=iactv(iab),nel
            dum=dum+ymat(jrep,irep)*gmat(irep,jrep,kk)
          enddo
        enddo
        velocity(kk)=dum*detratio
      enddo

!     if(iab.eq.2) write(ounit,*) 'ymat ',((ymat(jrep,irep),irep=iactv(iab),nel),jrep=ivirt(iab),norb)
      return
      end
!-----------------------------------------------------------------------
end module
