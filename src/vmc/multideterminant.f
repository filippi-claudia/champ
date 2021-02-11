      subroutine multideterminant_hpsi(vj,vpsp_det,eloc_det)

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use const, only: hb, nelec
      use csfs, only: nstates
      use dets, only: ndet
      use elec, only: ndn, nup
      use multidet, only: irepcol_det, ireporb_det, ivirt, iwundet, kref, numrep_det
      use optwf_contrl, only: ioptorb
      use ycompact, only: dymat, ymat
      use zcompact, only: aaz, dzmat, emz, zmat
      use coefs, only: norb
      use Bloc, only: b, tildem, xmat
      use denergy_det_m, only: denergy_det
      use dorb_m, only: iworbd
      use multimat, only: aa, wfmat
      use force_analy, only: iforce_analy
      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb
      use slater, only: d2dx2, ddx, fp, fpp, slmi
      use array_resize_utils, only: resize_matrix, resize_tensor
      use multislater, only: detiab

      implicit real*8(a-h,o-z)



      parameter (one=1.d0,half=0.5d0)

c note that the dimension of the slater matrices is assumed
c to be given by MMAT_DIM = (MELEC/2)**2, that is there are
c as many ups as downs. If this is not true then be careful if
c nelec is close to MELEC. The Slater matrices must be
c dimensioned at least max(nup**2,ndn**2)


      dimension eloc_det(MDET,2)
      dimension vj(3,nelec),vpsp_det(*)

      dimension btemp(nelec**2,2)

      ! call resize_matrix(b, norb+nadorb, 1)
      ! call resize_matrix(orb, norb+nadorb, 2)
      ! call resize_tensor(tildem, norb+nadorb, 2)
      ! call resize_tensor(aa, norb+nadorb, 2)

      nel=nup
      ish=0
      do iab=1,2
        if(iab.eq.2) then
          nel=ndn
          ish=nup
        endif
        eloc_det(kref,iab)=vpsp_det(iab)
        do i=1,nel
          eloc_det(kref,iab)=eloc_det(kref,iab)
     &    -hb*(d2dx2(i+ish)+2.d0*(vj(1,i+ish)*ddx(1,i+ish)+vj(2,i+ish)*ddx(2,i+ish)+vj(3,i+ish)*ddx(3,i+ish)))
        enddo
      enddo

c     write(6,*) 'eloc_ref',eloc_det(kref,1),eloc_det(kref,2)

      if(ndet.ne.1.or.iforce_analy.ne.0.or.ioptorb.ne.0) call bxmatrix(kref,xmat(1,1),xmat(1,2),b)
     
      if(ndet.eq.1.and.ioptorb.eq.0) return

      nel=nup
      iel=0
      do iab=1,2
        if(iab.eq.2) then
          nel=ndn
          iel=nup
        endif

c       ish=-nel
c       do 110 i=1,nel
c         ish=ish+nel
c         do 110 j=1,nel
c 110       btemp(j+ish,iab)=b(iworbd(j+iel,kref),i+iel)

c       do jrep=ivirt(iab),norb+nadorb
        do jrep=1,norb+nadorb

          do irep=1,nel
            dum1=0.d0
            dum2=0.d0
            dum3=0.d0
            do i=1,nel
              dum1=dum1+slmi(irep+(i-1)*nel,iab)*orb(i+iel,jrep)
              dum2=dum2+slmi(irep+(i-1)*nel,iab)*b(jrep,i+iel)
              dum3=dum3+xmat(i+(irep-1)*nel,iab)*orb(i+iel,jrep)
            enddo
            aa(irep,jrep,iab)=dum1
            tildem(irep,jrep,iab)=dum2-dum3
          enddo

c         do irep=1,nel
c           dum1=0.d0
c           do i=1,nel
c             dum4=0.d0
c             do kk=1,nel
c               dum4=dum4+btemp(kk+nel*(i-1),iab)*aa(kk,jrep,iab)
c             enddo
c             dum1=dum1+slmi(irep+(i-1)*nel,iab)*(b(jrep,i+iel)-dum4)
c           enddo
c           tildem(irep,jrep,iab)=dum1
c         enddo

        enddo

      enddo

c     if(kref.ne.1) then
c       do irep=1,13
c         write(6,'(''SLM  '',15f7.2)') (slmi(irep+(i-1)*ndn,2),i=1,13)
c       enddo
c       do irep=1,13
c         write(6,'(''AA-2 '',15f7.2)') (aa(irep,jrep,2),jrep=1,15)
c       enddo
c     endif

      denergy_det(kref,1)=0
      denergy_det(kref,2)=0

      if(ndet.eq.1) return

      do 200 k=1,ndet

        if(k.eq.kref) then
c         write(6,*) 'energy_det',eloc_det(k,1),eloc_det(k,2)
          goto 200
        endif

        do iab=1,2

          if(iwundet(k,iab).eq.k) then

            iel=0
            nel=nup
            if(iab.eq.2) then
              iel=nup
              nel=ndn
            endif
            ndim=numrep_det(k,iab)

            do irep=1,ndim
              iorb=irepcol_det(irep,k,iab)
              do jrep=1,ndim
                jorb=ireporb_det(jrep,k,iab)

                wfmat(irep+(jrep-1)*ndim,k,iab)=aa(iorb,jorb,iab)
              enddo
            enddo

c           if(kref.ne.1.and.k.eq.405) then
c           write(6,'('' AA det'',3i6)') k, iab,ndim
c           do irep=1,ndim
c             write(6,'(''AA'',10d12.4)') (wfmat(irep+(jrep-1)*ndim,k,iab),jrep=1,ndim)
c           enddo
c           endif

c           write(6,*) 'B HELLO',k,ndim,(irepcol_det(irep,k,iab),irep=1,ndim),(ireporb_det(jrep,k,iab),jrep=1,ndim)
            call matinv(wfmat(1,k,iab),ndim,det)
c           write(6,*) 'A HELLO',k,det
            detiab(k,iab)=det

c           if(kref.ne.1.and.k.eq.405) then
c           write(6,'('' AA det'',i6,d12.4)') k, det
c           do irep=1,ndim
c             write(6,'(''AA'',10d12.4)') (wfmat(irep+(jrep-1)*ndim,k,iab),jrep=1,ndim)
c           enddo
c           endif

            denergy_det(k,iab)=0
            do irep=1,ndim
              iorb=irepcol_det(irep,k,iab)
              do jrep=1,ndim
                jorb=ireporb_det(jrep,k,iab)
                denergy_det(k,iab)=denergy_det(k,iab)+wfmat(jrep+(irep-1)*ndim,k,iab)*tildem(iorb,jorb,iab)
              enddo
            enddo

          else
            index_det=iwundet(k,iab)

            denergy_det(k,iab)=denergy_det(index_det,iab)
            detiab(k,iab)=detiab(index_det,iab)

          endif

        enddo

c       write(6,*) 'denergy_det',denergy_det(k,1),denergy_det(k,2)

        eloc_det(k,1)=denergy_det(k,1)+eloc_det(kref,1)
        eloc_det(k,2)=denergy_det(k,2)+eloc_det(kref,2)

c       write(6,*) 'energy_det',eloc_det(k,1),eloc_det(k,2)
 200  continue

      do 400 k=1,ndet
        if(k.eq.kref) goto 400
        do 300 iab=1,2
          if(iwundet(k,iab).ne.kref) then
            detiab(k,iab)=detiab(k,iab)*detiab(kref,iab)
          endif
 300  continue
 400  continue


c compute Ymat for future use

      do 800 istate=1,nstates

        call compute_ymat(1,detiab(1,1),detiab(1,2),wfmat(1,1,1),ymat(1,1,1,istate),istate)
        if(iforce_analy.gt.0.or.ioptorb.gt.0) call compute_dymat(1,dymat(1,1,1,istate))

        if(ndn.gt.0) then
          call compute_ymat(2,detiab(1,1),detiab(1,2),wfmat(1,1,2),ymat(1,1,2,istate),istate)
          if(iforce_analy.gt.0.or.ioptorb.gt.0) call compute_dymat(2,dymat(1,1,2,istate))
        endif

        if(iforce_analy.gt.0.or.ioptorb.gt.0) call compute_zmat(ymat(1,1,1,istate),dymat(1,1,1,istate)
     &    ,zmat(1,1,1,istate),dzmat(1,1,1,istate),emz(1,1,1,istate),aaz(1,1,1,istate))
 800  continue

      return
      end

c-----------------------------------------------------------------------
      subroutine compute_ymat(iab,detu,detd,wfmat,ymat,istate)

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use const, only: nelec
      use dets, only: cdet, ndet
      use dets_equiv, only: cdet_equiv, dcdet_equiv
      use multidet, only: irepcol_det, ireporb_det, iwundet, kref, numrep_det
      use wfsec, only: iwf
      use coefs, only: norb
      use denergy_det_m, only: denergy_det

      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb
      use slater, only: d2dx2, ddx, fp, fpp, slmi

      implicit real*8(a-h,o-z)




      parameter (one=1.d0,half=0.5d0)


      dimension detu(MDET),detd(MDET),wfmat(MEXCIT**2,MDET),ymat(MORB,nelec)

      detrefi=1.d0/(detu(kref)*detd(kref))

      do 10 i=1,nelec
        do 10 j=1,norb
 10       ymat(j,i)=0

      do 300 k=1,ndet
        if(k.eq.kref) goto 300

        cdet_equiv(k)=0
        dcdet_equiv(k)=0

        if(iwundet(k,iab).eq.kref) goto 300

        kk=k
        if(iwundet(k,iab).ne.k) kk=iwundet(k,iab)

        if(iab.eq.1) then
          detall=detrefi*detu(kk)*detd(k)
         else
          detall=detrefi*detd(kk)*detu(k)
        endif

        cdet_equiv(kk)=cdet_equiv(kk)+cdet(k,istate,iwf)*detall
        dcdet_equiv(kk)=dcdet_equiv(kk)+cdet(k,istate,iwf)*detall*(denergy_det(k,1)+denergy_det(k,2))

 300  continue

      do 400 kk=1,ndet

        if(kk.eq.kref.or.iwundet(kk,iab).ne.kk) goto 400

        ndim=numrep_det(kk,iab)

        do irep=1,ndim
          iorb=irepcol_det(irep,kk,iab)
          do jrep=1,ndim
            jorb=ireporb_det(jrep,kk,iab)

            ymat(jorb,iorb)=ymat(jorb,iorb)+cdet_equiv(kk)*wfmat(jrep+(irep-1)*ndim,kk)
          enddo
        enddo

 400  continue

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_dymat(iab,dymat)

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use const, only: nelec
      use dets, only: ndet
      use dets_equiv, only: cdet_equiv, dcdet_equiv
      use multidet, only: irepcol_det, ireporb_det, iwundet, kref, numrep_det
      use coefs, only: norb
      use Bloc, only: tildem

      use multimat, only: wfmat

      implicit real*8(a-h,o-z)


      dimension dymat(MORB,nelec),dmat1(MEXCIT*MEXCIT),dmat2(MEXCIT*MEXCIT)

      do 10 i=1,nelec
        do 10 j=1,norb
 10       dymat(j,i)=0

      do 400 kk=1,ndet

        if(kk.eq.kref.or.iwundet(kk,iab).ne.kk) goto 400

        ndim=numrep_det(kk,iab)
        do irep=1,ndim
          iorb=ireporb_det(irep,kk,iab)
          do jrep=1,ndim
             jj=jrep+(irep-1)*ndim
             dmat1(jj)=0.d0
             do lrep=1,ndim
               lorb=irepcol_det(lrep,kk,iab)
               dmat1(jj)=dmat1(jj)+wfmat(jrep+(lrep-1)*ndim,kk,iab)*tildem(lorb,iorb,iab)
             enddo
          enddo
        enddo
        do irep=1,ndim
           do jrep=1,ndim
              jj=jrep+(irep-1)*ndim
              dmat2(jj)=0.d0
              do lrep=1,ndim
                 ll=jrep+(lrep-1)*ndim
                 dmat2(jj)=dmat2(jj)+dmat1(ll)*wfmat(lrep+(irep-1)*ndim,kk,iab)
              enddo
           enddo
        enddo

        do irep=1,ndim
          iorb=irepcol_det(irep,kk,iab)
          do jrep=1,ndim
            jorb=ireporb_det(jrep,kk,iab)

            jj=jrep+(irep-1)*ndim
            dymat(jorb,iorb)=dymat(jorb,iorb)+wfmat(jj,kk,iab)*dcdet_equiv(kk)-cdet_equiv(kk)*dmat2(jj)
          enddo
        enddo

 400  continue

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_zmat(ymat,dymat,zmat,dzmat,emz,aaz)

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use elec, only: ndn, nup
      use multidet, only: iactv, ivirt
      use coefs, only: norb
      use Bloc, only: tildem, xmat
      use multimat, only: aa

      use slater, only: d2dx2, ddx, fp, fpp, slmi

      implicit real*8(a-h,o-z)




      dimension ymat(MORB,nelec,2),dymat(MORB,nelec,2)
      dimension zmat(MORB,nelec,2),dzmat(MORB,nelec,2),emz(nelec,nelec,2),aaz(nelec,nelec,2)

      do 100 iab=1,2
        if(iab.eq.2.and.ndn.eq.0) goto 100

        if(iab.eq.1) then
          ish=0
          nel=nup
         else
          ish=nup
          nel=ndn
        endif

        do irep=1,nel
c         do jrep=ivirt(iab),norb+nadorb
          do jrep=ivirt(iab),norb
            zmat(jrep,irep,iab)=0
            dzmat(jrep,irep,iab)=0
            do krep=iactv(iab),nel
              zmat(jrep,irep,iab)=zmat(jrep,irep,iab)+ymat(jrep,krep,iab)*slmi(krep+(irep-1)*nel,iab)
              dzmat(jrep,irep,iab)=dzmat(jrep,irep,iab)+dymat(jrep,krep,iab)*slmi(krep+(irep-1)*nel,iab)
     &                                                 -ymat(jrep,krep,iab)*xmat(irep+(krep-1)*nel,iab)
            enddo
          enddo
        enddo

        do irep=1,nel
          do jrep=1,nel
            emz(jrep,irep,iab)=0
            aaz(jrep,irep,iab)=0
c           do krep=ivirt(iab),norb+nadorb
            do krep=ivirt(iab),norb
              emz(jrep,irep,iab)=emz(jrep,irep,iab)+tildem(jrep,krep,iab)*zmat(krep,irep,iab)
     &                           +aa(jrep,krep,iab)*dzmat(krep,irep,iab)
              aaz(jrep,irep,iab)=aaz(jrep,irep,iab)+aa(jrep,krep,iab)*zmat(krep,irep,iab)
            enddo
          enddo
        enddo

  100 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine update_ymat(iel)

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use const, only: nelec
      use csfs, only: nstates
      use elec, only: ndn, nup
      use ycompact, only: ymat
      use multimat, only: wfmat

      use multislater, only: detiab
      implicit real*8(a-h,o-z)




      if((iel.ne.nup.and.iel.ne.nelec).or.ndn.eq.0) return

      if(iel.eq.nup) then
        iab=2
       elseif(iel.eq.nelec) then
        iab=1
      endif

      do 100 istate=1,nstates
 100    call compute_ymat(iab,detiab(1,1),detiab(1,2),wfmat(1,1,iab),ymat(1,1,iab,istate),istate)

c     write(6,*) 'DU',(detiab(k,1),k=1,56)
c     write(6,*) 'DD',(detiab(k,2),k=1,56)
c     write(6,*) 'WF',((wfmat(i,k,iab),i=1,9),k=1,56)
c     do j=1,13
c     if(iab.eq.2) write(6,*) j,'YMAT 1',(ymat(i,j,iab,1),i=1,96)
c     enddo
c     do j=1,13
c     if(iab.eq.2) write(6,*) j,'YMAT 2',(ymat(i,j,iab,2),i=1,96)
c     enddo

      return
      end

c-----------------------------------------------------------------------
 
c-----------------------------------------------------------------------
      function idiff(j,i,iab)
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use multidet, only: irepcol_det, ireporb_det, numrep_det

      implicit real*8(a-h,o-z)

      idiff=1
      if(numrep_det(i,iab).ne.numrep_det(j,iab))return
      do k=1,numrep_det(i,iab)
        if(irepcol_det(k,j,iab).ne.irepcol_det(k,i,iab))return
        if(ireporb_det(k,j,iab).ne.ireporb_det(k,i,iab))return
      enddo
      idiff=0
      return
      end
           
c-----------------------------------------------------------------------
