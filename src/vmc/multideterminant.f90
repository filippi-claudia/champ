module multideterminant_mod
contains
      subroutine multideterminant_hpsi(vj,ekin_det,vpsp_det,eloc_det)

      use Bloc,    only: b, bkin,tildem,tildemkin,xmat, xmatkin
      use bxmatrices, only: bxmatrix
      use constants, only: hb
      use csfs, only: nstates
      use contrldmc, only: icut_e
      use denergy_det_m, only: denergy_det, allocate_denergy_det_m
      use fragments, only: eloc_i, elocfrag, ifragelec, nfrag
      use m_force_analytic, only: iforce_analy
      use matinv_mod, only: matinv
      use multidet, only: irepcol_det, ireporb_det, numrep_det, ndetiab, k_det, ndet_req
      use multidet, only: k_aux2, k_det2, ndetiab2, ndetsingle, ndetdouble
      use multimat, only: aa, wfmat
      use multislater, only: detiab
      use optwf_control, only: ioptorb
      use orbval, only: nadorb, orb
      use ycompact, only: dymat, ymat
      use zcompact, only: aaz, dzmat, emz, zmat
      use slater, only: norb
      use optwf_control, only: method
      use vmc_mod, only: nwftypeorb, nwftypejas, stoo, stoj, stobjx, nbjx
      use precision_kinds, only: dp
      use slater, only: ndet, iwundet, kref, d2dx2, ddx, slmi
      use system, only: nelec, ndn, nup
      use contrl_file, only: ounit

      implicit none

      integer :: i, iab, iel, index_det, iorb, kun, kw
      integer :: irep, ish, istate, jorb, j, o, x
      integer :: jrep, k, ndim, nel, ndim2, kk, kcum
      real(dp) :: det, dum1, dum2, dum3, dum4, dum5, deti, auxdet, tmpe
      real(dp), dimension(ndet, 2, nbjx) :: eloc_det
      real(dp), dimension(3, nelec, nwftypejas) :: vj
      real(dp), dimension(2, nbjx) :: vpsp_det
      real(dp), dimension(2, nbjx) :: ekin_det
      real(dp), dimension(nelec**2, 2) :: btemp
      real(dp), dimension(ndet_req,2) :: ddetiab
      real(dp), dimension(ndet_req,2) :: ddenergy_det

      do istate=1,nstates 
        j=stoj(istate) 
        o=stoo(istate)
        x=stobjx(istate)
        nel=nup
        ish=0
        do iab=1,2
          if(iab.eq.2) then
            nel=ndn
            ish=nup
          endif
          ekin_det(iab,x)=0.d0
          do i=1,nel
            tmpe=-hb*(d2dx2(i+ish,o)+2.d0*(vj(1,i+ish,j)*ddx(1,i+ish,o)+vj(2,i+ish,j)*ddx(2,i+ish,o)+vj(3,i+ish,j)*ddx(3,i+ish,o)))
            if (icut_e.lt.0) then
              eloc_i(i+ish)=eloc_i(i+ish)+tmpe
            endif
            if (nfrag.gt.1) then
              elocfrag(ifragelec(i+ish))=elocfrag(ifragelec(i+ish))+tmpe
            endif
            ekin_det(iab,x)=ekin_det(iab,x)+tmpe
          enddo
          eloc_det(kref,iab,x)=ekin_det(iab,x) + vpsp_det(iab,x)
        enddo

        if(ndet.ne.1.or.iforce_analy.ne.0.or.ioptorb.ne.0) then
          call bxmatrix(kref,xmat(1,1,x),xmat(1,2,x),b(1,1,x),x)
          call bxmatrix(kref,xmatkin(1,1,x),xmatkin(1,2,x),bkin(1,1,x),x)
        endif

        if(ndet.eq.1.and.ioptorb.eq.0) return

        nel=nup
        iel=0
        do iab=1,2
          if(iab.eq.2) then
            nel=ndn
            iel=nup
          endif

!         ish=-nel
!         do 110 i=1,nel
!           ish=ish+nel
!           do 110 j=1,nel
! 110         btemp(j+ish,iab)=b(iworbd(j+iel,kref),i+iel)

!         do jrep=ivirt(iab),norb+nadorb
          do jrep=1,norb+nadorb

            do irep=1,nel
              dum1=0.d0
              dum2=0.d0
              dum3=0.d0
              dum4=0.d0 
              dum5=0.d0
              do i=1,nel
                dum1=dum1+slmi(irep+(i-1)*nel,iab,o)*orb(i+iel,jrep,o)
                dum2=dum2+slmi(irep+(i-1)*nel,iab,o)*b(jrep,i+iel,x)
                dum3=dum3+xmat(i+(irep-1)*nel,iab,x)*orb(i+iel,jrep,o)

                dum4=dum4+slmi(irep+(i-1)*nel,iab,o)*bkin(jrep,i+iel,x)
                dum5=dum5+xmatkin(i+(irep-1)*nel,iab,x)*orb(i+iel,jrep,o)
              enddo
              aa(irep,jrep,iab,o)=dum1
              tildem(irep,jrep,iab,x)=dum2-dum3
              tildemkin(irep,jrep,iab,x)=dum4-dum5
            enddo

!           do irep=1,nel
!             dum1=0.d0
!             do i=1,nel
!               dum4=0.d0
!               do kk=1,nel
!                 dum4=dum4+btemp(kk+nel*(i-1),iab)*aa(kk,jrep,iab)
!               enddo
!               dum1=dum1+slmi(irep+(i-1)*nel,iab)*(b(jrep,i+iel)-dum4)
!             enddo
!             tildem(irep,jrep,iab)=dum1
!           enddo

          enddo

        enddo

!       if(kref.ne.1) then
!         do irep=1,13
!           write(ounit,'(''SLM  '',15f7.2)') (slmi(irep+(i-1)*ndn,2),i=1,13)
!         enddo
!         do irep=1,13
!           write(ounit,'(''AA-2 '',15f7.2)') (aa(irep,jrep,2),jrep=1,15)
!         enddo
!       endif

        call allocate_denergy_det_m()
        denergy_det=0.0d0
        ddenergy_det=0.0d0

        if(ndet.eq.1) return
        do iab=1,2
         
! loop inequivalent determinants
!
! determinants with single exitations
          if(ndetsingle(iab).ge.1)then
            do k=1,ndetsingle(iab)
              iorb=irepcol_det(1,k,iab)
              jorb=ireporb_det(1,k,iab)
              ddetiab(k,iab)=aa(iorb,jorb,iab,o)
              wfmat(k,1,iab,o)=1.0d0/ddetiab(k,iab)
              ddenergy_det(k,iab)=wfmat(k,1,iab,o)*tildem(iorb,jorb,iab,x)
            enddo
          endif

          kcum=ndetsingle(iab)+ndetdouble(iab)
         
! determinants double exitations
          if(ndetdouble(iab).ge.1)then
            do k=ndetsingle(iab)+1,kcum
           
!               ndim=numrep_det(k,iab)
!               do irep=1,ndim
!                  iorb=irepcol_det(irep,k,iab)
!                  do jrep=1,ndim
!                     jorb=ireporb_det(jrep,k,iab)
!                     wfmat(k,irep+(jrep-1)*ndim,iab)=aa(iorb,jorb,iab)
!                  enddo
!               enddo
!               ndim2=ndim*ndim
!               call matinv(wfmat(k,1:ndim2,iab),ndim,det)
!               ddetiab(k,iab)=det

              iorb=irepcol_det(1,k,iab)
              jorb=ireporb_det(1,k,iab)
              wfmat(k,1,iab,o)=aa(iorb,jorb,iab,o)
              jorb=ireporb_det(2,k,iab)
              wfmat(k,3,iab,o)=aa(iorb,jorb,iab,o)
              iorb=irepcol_det(2,k,iab)
              jorb=ireporb_det(1,k,iab)
              wfmat(k,2,iab,o)=aa(iorb,jorb,iab,o)
              jorb=ireporb_det(2,k,iab)
              wfmat(k,4,iab,o)=aa(iorb,jorb,iab,o)

!             call matinv(wfmat(k,1:4,iab),2,det)             
!             ddetiab(k,iab)=det
              ddetiab(k,iab)=wfmat(k,1,iab,o)*wfmat(k,4,iab,o) &
                   -wfmat(k,3,iab,o)*wfmat(k,2,iab,o)
              deti=1.d0/ddetiab(k,iab)
              auxdet=wfmat(k,1,iab,o)
              wfmat(k,1,iab,o)=wfmat(k,4,iab,o)*deti
              wfmat(k,2,iab,o)=-wfmat(k,2,iab,o)*deti
              wfmat(k,3,iab,o)=-wfmat(k,3,iab,o)*deti
              wfmat(k,4,iab,o)=auxdet*deti

!               do irep=1,ndim
!                  iorb=irepcol_det(irep,k,iab)
!                  do jrep=1,ndim
!                     jorb=ireporb_det(jrep,k,iab)
!                     ddenergy_det(k,iab)=ddenergy_det(k,iab)+wfmat(k,jrep+(irep-1)*ndim,iab)*tildem(iorb,jorb,iab)
!                  enddo
!               enddo

              iorb=irepcol_det(1,k,iab)
              jorb=ireporb_det(1,k,iab)
              ddenergy_det(k,iab)=ddenergy_det(k,iab)+wfmat(k,1,iab,o)*tildem(iorb,jorb,iab,x)
              jorb=ireporb_det(2,k,iab)
              ddenergy_det(k,iab)=ddenergy_det(k,iab)+wfmat(k,2,iab,o)*tildem(iorb,jorb,iab,x)
              iorb=irepcol_det(2,k,iab)
              jorb=ireporb_det(1,k,iab)
              ddenergy_det(k,iab)=ddenergy_det(k,iab)+wfmat(k,3,iab,o)*tildem(iorb,jorb,iab,x)
              jorb=ireporb_det(2,k,iab)
              ddenergy_det(k,iab)=ddenergy_det(k,iab)+wfmat(k,4,iab,o)*tildem(iorb,jorb,iab,x)
               
            enddo
          endif

          if(kcum.lt.ndetiab(iab))then
! determinants multiple exitations
            do k=kcum+1,ndetiab(iab)
           
              ndim=numrep_det(k,iab)
              do irep=1,ndim
                iorb=irepcol_det(irep,k,iab)
                do jrep=1,ndim
                  jorb=ireporb_det(jrep,k,iab)
                  wfmat(k,irep+(jrep-1)*ndim,iab,o)=aa(iorb,jorb,iab,o)
                  !write(ounit,*) k,iab,iorb,jorb,aa(iorb,jorb,iab,o)
                enddo
              enddo
           
              ndim2=ndim*ndim
              call matinv(wfmat(k,1:ndim2,iab,o),ndim,det)
              ddetiab(k,iab)=det
              !write(ounit,*) k,iab,det
           
              do irep=1,ndim
                iorb=irepcol_det(irep,k,iab)
                do jrep=1,ndim
                  jorb=ireporb_det(jrep,k,iab)
                  ddenergy_det(k,iab)=ddenergy_det(k,iab)+wfmat(k,jrep+(irep-1)*ndim,iab,o)*tildem(iorb,jorb,iab,x)
                enddo
              enddo
            enddo
          endif
        
! unrolling determinants different to kref
          detiab(:,iab,o)=detiab(kref,iab,o)
          eloc_det(:,iab,x)=eloc_det(kref,iab,x)
          do kk=1,ndetiab2(iab)
            k=k_det2(kk,iab)
            kw=k_aux2(kk,iab)
            detiab(k,iab,o)=detiab(k,iab,o)*ddetiab(kw,iab)
            denergy_det(k,iab,x)=ddenergy_det(kw,iab)
            eloc_det(k,iab,x)=eloc_det(k,iab,x)+denergy_det(k,iab,x)
          enddo

        enddo ! end iab loop
         
! compute Ymat for future use

        call compute_ymat(1,detiab(1,1,o),detiab(1,2,o),wfmat(1,1,1,o),ymat(1,1,1,istate),istate)
!        if(iforce_analy.gt.0.or.(ioptorb.gt.0.and.(method(1:3) == 'lin'))) call compute_dymat(1,dymat(1,1,1,istate))
        if(iforce_analy.gt.0.or.ioptorb.gt.0) call compute_dymat(1,dymat(1,1,1,istate),istate)

        if(ndn.gt.0) then
          call compute_ymat(2,detiab(1,1,o),detiab(1,2,o),wfmat(1,1,2,o),ymat(1,1,2,istate),istate)
!          if(iforce_analy.gt.0.or.(ioptorb.gt.0.and.(method(1:3) == 'lin'))) call compute_dymat(2,dymat(1,1,2,istate))
          if(iforce_analy.gt.0.or.ioptorb.gt.0) call compute_dymat(2,dymat(1,1,2,istate),istate)
        endif

!        if(iforce_analy.gt.0.or.(ioptorb.gt.0.and.(method(1:3) == 'lin'))) call compute_zmat(ymat(1,1,1,istate),dymat(1,1,1,istate)
        if(iforce_analy.gt.0.or.ioptorb.gt.0) call compute_zmat(ymat(1,1,1,istate),dymat(1,1,1,istate) &
                   ,zmat(1,1,1,istate),dzmat(1,1,1,istate),emz(1,1,1,istate),aaz(1,1,1,istate),istate)

      enddo ! end of istate loop from start of routine

      return
      end

!-----------------------------------------------------------------------
      subroutine compute_ymat(iab,detu,detd,wfmat,ymat,istate)

      use denergy_det_m, only: denergy_det
      use multidet, only: irepcol_det, ireporb_det, numrep_det, k_det, ndetiab
      use multidet, only: k_aux2, k_det2, ndetiab2, ndetsingle, ndetdouble
      use multiple_geo, only: iwf
      use precision_kinds, only: dp
      use slater, only: ndet, cdet, iwundet, kref, cdet_equiv, dcdet_equiv
      use slater, only: norb
      use system, only: nelec
      use vmc_mod, only: MEXCIT, norb_tot, stoo, nwftypeorb, stobjx
      use contrl_file, only: ounit


      implicit none

      integer :: i, iab, iorb, irep, istate
      integer :: j, jorb, jrep, k, kun, kw
      integer :: kk, ndim, ndim2, kcum, iwf_save
      real(dp) :: detall, detrefi
      real(dp), dimension(ndet) :: detu
      real(dp), dimension(ndet) :: detd
      real(dp), dimension(ndet, MEXCIT**2) :: wfmat
      real(dp), dimension(norb_tot*nelec) :: ymat
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: half = 0.5d0


      real(dp), dimension(ndetiab2(iab)) :: detallv
      real(dp), dimension(ndetiab2(iab)) :: sumde

      detrefi=1.d0/(detu(kref)*detd(kref))

      ymat=0
      cdet_equiv=0.0d0
      dcdet_equiv=0.0d0

      iwf_save=iwf
      if(nwftypeorb.gt.1) iwf=1

! unroling determinants different to kref
      do kk=1,ndetiab2(iab)
         k=k_det2(kk,iab)
         kw=k_aux2(kk,iab)
         detall=detrefi*detu(k)*detd(k)*cdet(k,istate,iwf)
         cdet_equiv(kw)=cdet_equiv(kw)+detall
         dcdet_equiv(kw)=dcdet_equiv(kw)+detall*(denergy_det(k,1,stobjx(istate))+denergy_det(k,2,stobjx(istate))) 
      enddo

! loop over single exitations
      if(ndetsingle(iab).ge.1)then
         do kk=1,ndetsingle(iab)
            iorb=irepcol_det(1,kk,iab)
            jorb=ireporb_det(1,kk,iab)
            ymat(jorb+norb_tot*(iorb-1))=ymat(jorb+norb_tot*(iorb-1))+cdet_equiv(kk)*wfmat(kk,1)
         enddo
      endif

      kcum=ndetsingle(iab)+ndetdouble(iab)
      if(ndetdouble(iab).ge.1)then

         do kk=ndetsingle(iab)+1,kcum !ndetiab(iab)
         
!            ndim=numrep_det(kk,iab)
!            do irep=1,ndim
!               iorb=irepcol_det(irep,kk,iab)
!               do jrep=1,ndim
!                  jorb=ireporb_det(jrep,kk,iab)
!                  ymat(jorb+norb_tot*(iorb-1))=ymat(jorb+norb_tot*(iorb-1))+cdet_equiv(kk)*wfmat(kk,jrep+(irep-1)*ndim)
!               enddo  
!            enddo

            iorb=irepcol_det(1,kk,iab)
            jorb=ireporb_det(1,kk,iab)
            ymat(jorb+norb_tot*(iorb-1))=ymat(jorb+norb_tot*(iorb-1))+cdet_equiv(kk)*wfmat(kk,1)
            jorb=ireporb_det(2,kk,iab)
            ymat(jorb+norb_tot*(iorb-1))=ymat(jorb+norb_tot*(iorb-1))+cdet_equiv(kk)*wfmat(kk,2)
            
            iorb=irepcol_det(2,kk,iab)
            jorb=ireporb_det(1,kk,iab)
            ymat(jorb+norb_tot*(iorb-1))=ymat(jorb+norb_tot*(iorb-1))+cdet_equiv(kk)*wfmat(kk,3)
            jorb=ireporb_det(2,kk,iab)
            ymat(jorb+norb_tot*(iorb-1))=ymat(jorb+norb_tot*(iorb-1))+cdet_equiv(kk)*wfmat(kk,4)
         enddo

      endif


      if(kcum.lt.ndetiab(iab))then

         do kk=kcum+1,ndetiab(iab)
            ndim=numrep_det(kk,iab)
            do irep=1,ndim
               iorb=irepcol_det(irep,kk,iab)
               do jrep=1,ndim
                  jorb=ireporb_det(jrep,kk,iab)
                  ymat(jorb+norb_tot*(iorb-1))=ymat(jorb+norb_tot*(iorb-1))+cdet_equiv(kk)*wfmat(kk,jrep+(irep-1)*ndim)
               enddo   
            enddo
         enddo
      endif


      return
      end
!-----------------------------------------------------------------------
      subroutine compute_dymat(iab,dymat,istate)

      use vmc_mod, only: norb_tot, MEXCIT, stoo, stobjx
      use system, only: nelec
      use multidet, only: irepcol_det, ireporb_det, numrep_det, ndetiab, ndetsingle, ndetdouble
      use slater, only: iwundet, kref, norb, ndet, cdet_equiv, dcdet_equiv
      use Bloc, only: tildem
      use multimat, only: wfmat
      use precision_kinds, only: dp
      use slater,  only: cdet_equiv,dcdet_equiv,iwundet,kref,ndet,norb
      use system,  only: nelec
      use vmc_mod, only: MEXCIT,norb_tot


      implicit none

      integer :: i, iab, iorb, irep, j, istate, o, x
      integer :: jj, jorb, jrep, kk, kcum
      integer :: ll, lorb, lrep, ndim

      real(dp), dimension(norb_tot, nelec) :: dymat
      real(dp), dimension(MEXCIT*MEXCIT) :: dmat1
      real(dp), dimension(MEXCIT*MEXCIT) :: dmat2

      dymat=0.0d0
      o=stoo(istate)
      x=stobjx(istate)
! loop over single exitations      
      if(ndetsingle(iab).ge.1) then
         do kk=1,ndetsingle(iab)
            iorb=ireporb_det(1,kk,iab)
            jorb=irepcol_det(1,kk,iab)
            dmat1(1)=wfmat(kk,1,iab,o)*tildem(jorb,iorb,iab,x)
            dmat2(1)=dmat1(1)*wfmat(kk,1,iab,o)
            dymat(iorb,jorb)=dymat(iorb,jorb)+wfmat(kk,1,iab,o)*dcdet_equiv(kk)-cdet_equiv(kk)*dmat2(1)
         enddo
      endif

      kcum=ndetsingle(iab)+ndetdouble(iab)

! double excitations
      if(ndetdouble(iab).ge.1) then
         do kk=ndetsingle(iab)+1,kcum

            dmat1(1:4)=0.d0
            
            iorb=ireporb_det(1,kk,iab)
            lorb=irepcol_det(1,kk,iab)
            dmat1(1)=dmat1(1)+wfmat(kk,1,iab,o)*tildem(lorb,iorb,iab,x)
            dmat1(2)=dmat1(2)+wfmat(kk,2,iab,o)*tildem(lorb,iorb,iab,x)
            lorb=irepcol_det(2,kk,iab)
            dmat1(1)=dmat1(1)+wfmat(kk,3,iab,o)*tildem(lorb,iorb,iab,x)
            dmat1(2)=dmat1(2)+wfmat(kk,4,iab,o)*tildem(lorb,iorb,iab,x)

            iorb=ireporb_det(2,kk,iab)
            lorb=irepcol_det(1,kk,iab)
            dmat1(3)=dmat1(3)+wfmat(kk,1,iab,o)*tildem(lorb,iorb,iab,x)
            dmat1(4)=dmat1(4)+wfmat(kk,2,iab,o)*tildem(lorb,iorb,iab,x)
            lorb=irepcol_det(2,kk,iab)
            dmat1(3)=dmat1(3)+wfmat(kk,3,iab,o)*tildem(lorb,iorb,iab,x)
            dmat1(4)=dmat1(4)+wfmat(kk,4,iab,o)*tildem(lorb,iorb,iab,x)

            dmat2(1:4)=0.d0
            dmat2(1)=dmat1(1)*wfmat(kk,1,iab,o)+dmat1(3)*wfmat(kk,2,iab,o)
            dmat2(2)=dmat1(2)*wfmat(kk,1,iab,o)+dmat1(4)*wfmat(kk,2,iab,o)
            dmat2(3)=dmat1(1)*wfmat(kk,3,iab,o)+dmat1(3)*wfmat(kk,4,iab,o)
            dmat2(4)=dmat1(2)*wfmat(kk,3,iab,o)+dmat1(4)*wfmat(kk,4,iab,o)
           
            iorb=irepcol_det(1,kk,iab)
            jorb=ireporb_det(1,kk,iab)
            dymat(jorb,iorb)=dymat(jorb,iorb)+wfmat(kk,1,iab,o)*dcdet_equiv(kk)-cdet_equiv(kk)*dmat2(1)
            jorb=ireporb_det(2,kk,iab)
            dymat(jorb,iorb)=dymat(jorb,iorb)+wfmat(kk,2,iab,o)*dcdet_equiv(kk)-cdet_equiv(kk)*dmat2(2)

            iorb=irepcol_det(2,kk,iab)
            jorb=ireporb_det(1,kk,iab)
            dymat(jorb,iorb)=dymat(jorb,iorb)+wfmat(kk,3,iab,o)*dcdet_equiv(kk)-cdet_equiv(kk)*dmat2(3)
            jorb=ireporb_det(2,kk,iab)
            dymat(jorb,iorb)=dymat(jorb,iorb)+wfmat(kk,4,iab,o)*dcdet_equiv(kk)-cdet_equiv(kk)*dmat2(4)
            
         enddo
      endif

! multiple excitations
      if(kcum.lt.ndetiab(iab)) then
        do kk=kcum+1,ndetiab(iab)
         
          ndim=numrep_det(kk,iab)
         
          do irep=1,ndim
            iorb=ireporb_det(irep,kk,iab)
            do jrep=1,ndim
              jj=jrep+(irep-1)*ndim
              dmat1(jj)=0.d0
              do lrep=1,ndim
                lorb=irepcol_det(lrep,kk,iab)
                dmat1(jj)=dmat1(jj)+wfmat(kk,jrep+(lrep-1)*ndim,iab,o)*tildem(lorb,iorb,iab,x)
              enddo
            enddo
          enddo
         
          do irep=1,ndim
            do jrep=1,ndim
              jj=jrep+(irep-1)*ndim
              dmat2(jj)=0.d0
              do lrep=1,ndim
                ll=jrep+(lrep-1)*ndim
                dmat2(jj)=dmat2(jj)+dmat1(ll)*wfmat(kk,lrep+(irep-1)*ndim,iab,o)
              enddo
            enddo
          enddo
         
          do irep=1,ndim
            iorb=irepcol_det(irep,kk,iab)
            do jrep=1,ndim
              jorb=ireporb_det(jrep,kk,iab)                 
              jj=jrep+(irep-1)*ndim
              dymat(jorb,iorb)=dymat(jorb,iorb)+wfmat(kk,jj,iab,o)*dcdet_equiv(kk)-cdet_equiv(kk)*dmat2(jj)
           enddo
          enddo
         
        enddo
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine compute_zmat(ymat,dymat,zmat,dzmat,emz,aaz,istate)

      use vmc_mod, only: norb_tot, stoo, stobjx
      use multidet, only: iactv, ivirt
      use slater, only: norb
      use Bloc, only: tildem, xmat
      use multimat, only: aa
      use slater, only: slmi
      use system, only: nelec, ndn, nup

      use precision_kinds, only: dp
      use slater,  only: norb,slmi
      use system,  only: ndn,nelec,nup
      use vmc_mod, only: norb_tot

      implicit none

      integer :: iab, irep, ish, jrep, krep
      integer :: nel, istate, x, o

      real(dp), dimension(norb_tot, nelec, 2) :: ymat
      real(dp), dimension(norb_tot, nelec, 2) :: dymat
      real(dp), dimension(norb_tot, nelec, 2) :: zmat
      real(dp), dimension(norb_tot, nelec, 2) :: dzmat
      real(dp), dimension(nelec, nelec, 2) :: emz
      real(dp), dimension(nelec, nelec, 2) :: aaz

      o=stoo(istate)
      x=stobjx(istate)

      do iab=1,2
        if(iab.eq.2.and.ndn.eq.0) goto 100

        if(iab.eq.1) then
          ish=0
          nel=nup
         else
          ish=nup
          nel=ndn
        endif

        do irep=1,nel
          do jrep=ivirt(iab),norb
            zmat(jrep,irep,iab)=0
            dzmat(jrep,irep,iab)=0
            do krep=iactv(iab),nel
              zmat(jrep,irep,iab)=zmat(jrep,irep,iab)+ymat(jrep,krep,iab)*slmi(krep+(irep-1)*nel,iab,o)
              dzmat(jrep,irep,iab)=dzmat(jrep,irep,iab)+dymat(jrep,krep,iab)*slmi(krep+(irep-1)*nel,iab,o) &
              -ymat(jrep,krep,iab)*xmat(irep+(krep-1)*nel,iab,x)
            enddo
          enddo
        enddo

        do irep=1,nel
          do jrep=1,nel
            emz(jrep,irep,iab)=0
            aaz(jrep,irep,iab)=0
            do krep=ivirt(iab),norb
              emz(jrep,irep,iab)=emz(jrep,irep,iab)+tildem(jrep,krep,iab,x)*zmat(krep,irep,iab) &
                                 +aa(jrep,krep,iab,o)*dzmat(krep,irep,iab)
              aaz(jrep,irep,iab)=aaz(jrep,irep,iab)+aa(jrep,krep,iab,o)*zmat(krep,irep,iab)
            enddo
          enddo
        enddo

      100 continue
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine update_ymat(iel)

      use system, only: nelec, ndn, nup
      use csfs, only: nstates
      use ycompact, only: ymat
      use multimat, only: wfmat
      use vmc_mod, only: stoo

      use multislater, only: detiab
      use system,  only: ndn,nelec,nup
      use ycompact, only: ymat

      implicit none

      integer :: iab, iel, istate

      if((iel.ne.nup.and.iel.ne.nelec).or.ndn.eq.0) return

      if(iel.eq.nup) then
        iab=2
       elseif(iel.eq.nelec) then
        iab=1
      endif

      do istate=1,nstates
         call compute_ymat(iab,detiab(1,1,stoo(istate)),detiab(1,2,stoo(istate)), &
                  wfmat(:,:,iab,stoo(istate)),ymat(1,1,iab,istate),istate)
      enddo

      return
      end

!-----------------------------------------------------------------------
      function idiff0(j,i,iab)
! used when setting up the determinant list and numrep, irepcol, ireporb available for all dets

      use multidet, only: irepcol_det, ireporb_det, numrep_det
      use contrl_file, only: ounit

      implicit none

      integer :: i, iab, j, k
      integer :: idiff0

      idiff0=1

      if(numrep_det(i,iab).ne.numrep_det(j,iab)) return
      do k=1,numrep_det(i,iab)
        if(irepcol_det(k,j,iab).ne.irepcol_det(k,i,iab)) return
        if(ireporb_det(k,j,iab).ne.ireporb_det(k,i,iab)) return
      enddo

      idiff0=0

      return
      end

!-----------------------------------------------------------------------
      function idiff(j,i,iab)
! used when iwundet available

      use slater, only: iwundet
      use contrl_file, only: ounit

      implicit none

      integer :: i, iab, j
      integer :: idiff

      idiff=1
      if(iwundet(i,iab).ne.iwundet(j,iab)) return

      idiff=0

      return
      end

!-----------------------------------------------------------------------
end module
