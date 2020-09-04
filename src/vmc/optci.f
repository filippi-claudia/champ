      subroutine optci_deloc(eloc_det,e_other,psid,energy)

      use csfs, only: cxdet, iadet, ibdet, icxdet, ncsf
      use optwf_contrl, only: ioptci
      use ci000, only: iciprt, nciprim, nciterm
      use ci001_blk, only: ci_o, ci_oe
      use ci003_blk, only: ci_e, ci_e_old
      use ci004_blk, only: ci_de, ci_de_old

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'pseudo.h'
      include 'optci.h'

      common /multislater/ detiab(MDET,2)

      dimension ciprim(MDET),cieprim(MDET)
      dimension eloc_det(MDET,2)
      
      if(ioptci.eq.0) return 
      
      psidi=1.d0/psid

      do 1 k=1,nciprim
        ciprim(k)=detiab(k,1)*detiab(k,2)*psidi
        cieprim(k)=(eloc_det(k,1)+eloc_det(k,2)+e_other)*ciprim(k)
   1  continue

c Update <Oi>,<Ei>,<dEi>,<Oi*Ej> 
c Correlation matrix <Oi*Oj> is computed in ci_sum

      if(ncsf.eq.0) then
        do 10 i=1,nciprim
         ci_o(i)=ciprim(i)
         ci_e(i)=cieprim(i)
  10     ci_de(i)=cieprim(i)-ciprim(i)*energy
       else
        do 30 icsf=1,ncsf
          ci_o_csf=0
          ci_e_csf=0
          do 20 ix=iadet(icsf),ibdet(icsf)
            idet=icxdet(ix)
            ci_o_csf=ci_o_csf+ciprim(idet)*cxdet(ix)
  20        ci_e_csf=ci_e_csf+cieprim(idet)*cxdet(ix)
        ci_o(icsf)=ci_o_csf
        ci_e(icsf)=ci_e_csf
  30    ci_de(icsf)=ci_e_csf-ci_o_csf*energy

      endif

      if(method.eq.'sr_n'.or.method.eq.'lin_d') return

      if(ncsf.eq.0) then
       do 50 i=1,nciprim
         do 50 j=1,nciprim
  50       ci_oe(i,j) = ciprim(i)*cieprim(j)
      else
        do 60 icsf=1,ncsf
          do 60 jcsf=1,ncsf
  60        ci_oe(icsf,jcsf)=ci_o(icsf)*ci_e(jcsf)
      endif

      end
c-----------------------------------------------------------------------
      subroutine optci_init(iflg)

      use optwf_contrl, only: ioptci
      use ci000, only: iciprt, nciprim, nciterm
      use ci005_blk, only: ci_o_cum, ci_o_sum
      use ci006_blk, only: ci_de_cum, ci_de_sum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum, ci_oe_sum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum, ci_oo_sum
      use ci010_blk, only: ci_ooe_cum, ci_ooe_sum

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'optci.h'

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do 10 i=1,nciterm
       ci_o_sum(i) =0.d0
       ci_de_sum(i) =0.d0
       do 10 j=1,nciterm
  10    ci_oe_sum(j,i) =0.d0

      idx=0
      do 20 i=1,nciterm
       do 20 j=1,i
        idx=idx+1
        ci_oo_sum(idx)=0.d0
  20    ci_ooe_sum(idx)=0.d0

C$ iflg = 0: init *cum, *cm2 as well
      if(iflg.gt.0) return

      guid_weight=0.d0
      guid_weight_sq=0.d0

      do 30 i=1,nciterm
       ci_o_cum(i)=0.d0
       ci_de_cum(i)=0.d0
       do 30 j=1,nciterm
        ci_oe_cum(j,i)=0.d0
  30    ci_oe_cm2(j,i)=0.d0

      idx=0
      do 40 i=1,nciterm
       do 40 j=1,i
        idx=idx+1
        ci_oo_cum(idx)=0.d0
        ci_oo_cm2(idx)=0.d0
  40    ci_ooe_cum(idx)=0.d0

      end

c-----------------------------------------------------------------------
      subroutine optci_save

      use optwf_contrl, only: ioptci
      use ci000, only: iciprt, nciprim, nciterm
      use ci001_blk, only: ci_o, ci_oe
      use ci002_blk, only: ci_o_old, ci_oe_old
      use ci003_blk, only: ci_e, ci_e_old
      use ci004_blk, only: ci_de, ci_de_old

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'optci.h'

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do 10 i=1,nciterm
       ci_o_old(i)=ci_o(i)
       ci_e_old(i)=ci_e(i)
       ci_de_old(i)=ci_de(i)
       do 10 j=1,nciterm
  10    ci_oe_old(j,i)=ci_oe(j,i)

      end
c-----------------------------------------------------------------------
      subroutine optci_restore

      use optwf_contrl, only: ioptci
      use ci000, only: iciprt, nciprim, nciterm
      use ci001_blk, only: ci_o, ci_oe
      use ci002_blk, only: ci_o_old, ci_oe_old
      use ci003_blk, only: ci_e, ci_e_old
      use ci004_blk, only: ci_de, ci_de_old

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'optci.h'

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do 10 i=1,nciterm
       ci_o(i)=ci_o_old(i)
       ci_e(i)=ci_e_old(i)
       ci_de(i)=ci_de_old(i)
       do 10 j=1,nciterm
  10    ci_oe(j,i)=ci_oe_old(j,i)

      end
c-----------------------------------------------------------------------
      subroutine optci_sum(p,q,enew,eold)

      use optwf_contrl, only: ioptci
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      use ci000, only: iciprt, nciprim, nciterm
      use ci001_blk, only: ci_o, ci_oe
      use ci002_blk, only: ci_o_old, ci_oe_old
      use ci004_blk, only: ci_de, ci_de_old
      use ci005_blk, only: ci_o_cum, ci_o_sum
      use ci006_blk, only: ci_de_cum, ci_de_sum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum, ci_oe_sum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum, ci_oo_sum
      use ci010_blk, only: ci_ooe_cum, ci_ooe_sum

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'optci.h'

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do 10 j=1,nciterm
       ci_o_sum(j) =ci_o_sum(j)+p*ci_o(j)+q*ci_o_old(j)
       ci_de_sum(j) =ci_de_sum(j)+p*ci_de(j)+q*ci_de_old(j)
       do 10 i=1,nciterm
  10    ci_oe_sum(i,j) =ci_oe_sum(i,j)+p*ci_oe(i,j)+q*ci_oe_old(i,j)

      idx=0
      do 20 i=1,nciterm
       do 20 j=1,i
        idx=idx+1
        ci_oo_new=ci_o(i)*ci_o(j)
        ci_oo_old=ci_o_old(i)*ci_o_old(j)
        ci_oo_sum(idx)=ci_oo_sum(idx)  +p*ci_oo_new+q*ci_oo_old
  20    ci_ooe_sum(idx)=ci_ooe_sum(idx)+p*ci_oo_new*enew+q*ci_oo_old*eold

      end
c-----------------------------------------------------------------------
      subroutine optci_cum(wsum)

      use optwf_contrl, only: ioptci
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      use ci000, only: iciprt, nciprim, nciterm

      use ci005_blk, only: ci_o_cum, ci_o_sum
      use ci006_blk, only: ci_de_cum, ci_de_sum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum, ci_oe_sum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum, ci_oo_sum
      use ci010_blk, only: ci_ooe_cum, ci_ooe_sum

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'optci.h'

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do 10 i=1,nciterm
       ci_o_cum(i)=ci_o_cum(i)+ci_o_sum(i)
       ci_de_cum(i)=ci_de_cum(i)+ci_de_sum(i)
       do 10 j=1,nciterm
        ci_oe_now=ci_oe_sum(i,j)/wsum
        ci_oe_cum(i,j)=ci_oe_cum(i,j)+ci_oe_sum(i,j)
  10    ci_oe_cm2(i,j)=ci_oe_cm2(i,j)+ci_oe_sum(i,j)*ci_oe_now

      idx=0
      do 20 i=1,nciterm
       do 20 j=1,i
        idx=idx+1
        ci_oo_cum(idx)=ci_oo_cum(idx)+ci_oo_sum(idx)
        ci_oo_cm2(idx)=ci_oo_cm2(idx)+ci_oo_sum(idx)*ci_oo_sum(idx)/wsum
  20    ci_ooe_cum(idx)=ci_ooe_cum(idx)+ci_ooe_sum(idx)

      end
c-----------------------------------------------------------------------
      subroutine optci_dump(iu)

      use optwf_contrl, only: ioptci
      use ci000, only: iciprt, nciprim, nciterm
      use ci005_blk, only: ci_o_cum, ci_o_sum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum, ci_oe_sum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum, ci_oo_sum
      use ci010_blk, only: ci_ooe_cum, ci_ooe_sum

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'optci.h'


      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      matdim=nciterm*(nciterm+1)/2
      write(iu) nciprim,nciterm
      write(iu) (ci_o_cum(i),i=1,nciterm)
      write(iu) ((ci_oe_cum(i,j),ci_oe_cm2(i,j),i=1,nciterm),j=1,nciterm)
      write(iu) (ci_oo_cum(i),ci_oo_cm2(i),i=1,matdim)
      write(iu) (ci_ooe_cum(i),i=1,matdim)

      end
c-----------------------------------------------------------------------
      subroutine optci_rstrt(iu)
      use optwf_contrl, only: ioptci
      use ci000, only: iciprt, nciprim, nciterm

      use ci005_blk, only: ci_o_cum, ci_o_sum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum, ci_oe_sum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum, ci_oo_sum
      use ci010_blk, only: ci_ooe_cum, ci_ooe_sum

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'optci.h'

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      read(iu) mciprim,mciterm
      if(mciprim.ne.nciprim) then
       write (6,*) 'wrong number of primitive ci terms!'
       write (6,*) 'old ',mciprim,' new ',nciprim
       call fatal_error('CI: Restart, inconsistent CI information')
      endif
      if(mciterm.ne.nciterm) then
       write (6,*) 'wrong number of ci terms!'
       write (6,*) 'old ',mciterm,' new ',nciterm
       call fatal_error('CI: Restart, inconsistent CI information')
      endif

      read(iu) (ci_o_cum(i),i=1,nciterm)
      read(iu) ((ci_oe_cum(i,j),ci_oe_cm2(i,j),i=1,nciterm),j=1,nciterm)
      matdim=nciterm*(nciterm+1)/2 
      read(iu) (ci_oo_cum(i),ci_oo_cm2(i),i=1,matdim)
      read(iu) (ci_ooe_cum(i),i=1,matdim)

      end
c-----------------------------------------------------------------------
      subroutine optci_avrg(wcum,iblk,oav,deav,oeav,oeerr,ooav,ooerr,ooeav)

      use optwf_contrl, only: ioptci
      use ci000, only: iciprt, nciprim, nciterm
      use ci005_blk, only: ci_o_cum, ci_o_sum
      use ci006_blk, only: ci_de_cum, ci_de_sum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum, ci_oe_sum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum, ci_oo_sum
      use ci010_blk, only: ci_ooe_cum, ci_ooe_sum

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'optci.h'

      dimension oav(MXCITERM),deav(MXCITERM)
      dimension oeav(MXCITERM,MXCIREDUCED),oeerr(MXCITERM,MXCIREDUCED)
      dimension ooav(MXCIMATDIM),ooerr(MXCIMATDIM)
      dimension ooeav(MXCIMATDIM)

      err(x,x2)=dsqrt(abs(x2/wcum-(x/wcum)**2)/iblk)

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do 20 i=1,nciterm
       oav(i)=ci_o_cum(i)/wcum
       deav(i)=ci_de_cum(i)/wcum
       do 10 j=1,nciterm
  10    oeav(i,j)=ci_oe_cum(i,j)/wcum

       do 20 j=1,nciterm
  20    oeerr(i,j)=err(ci_oe_cum(i,j),ci_oe_cm2(i,j))

      idx=0
      do 30 i=1,nciterm
       do 30 j=1,i
        idx=idx+1
        ooav(idx)=ci_oo_cum(idx)/wcum
        ooeav(idx)=ci_ooe_cum(idx)/wcum
  30    ooerr(idx)=err(ci_oo_cum(idx),ci_oo_cm2(idx))

      end
c-----------------------------------------------------------------------
      subroutine optci_fin(iblk,passes,etot)

      use csfs, only: ccsf, ncsf
      use dets, only: cdet
      use gradhess_ci, only: grad_ci, h_ci, s_ci
      use linear_norm, only: oav, ci_oav
      use optwf_contrl, only: ioptci, ioptjas, ioptorb
      use ci000, only: iciprt, nciprim, nciterm

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optci.h'

      dimension deav(MXCITERM)
      dimension oeav(MXCITERM,MXCIREDUCED),oeerr(MXCITERM,MXCIREDUCED)
      dimension ooav(MXCIMATDIM),ooerr(MXCIMATDIM)
      dimension ooeav(MXCIMATDIM)
      dimension oelocav(MXCITERM),eav(MXCITERM)

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      iciprt_sav=iciprt
      iciprt=-1
      call optci_avrg(passes,iblk,oav,deav,oeav,oeerr,ooav,ooerr,ooeav)
      iciprt=iciprt_sav
     
      if(ncsf.eq.0) then
        do 20 i=1,nciterm
          ci_oav(i)=oav(i)
          oelocav(i)=0
          eav(i)=0
          do 20 j=1,nciterm
            oelocav(i)=oelocav(i)+oeav(i,j)*cdet(j,1,1)
  20        eav(i)=eav(i)+oeav(j,i)*cdet(j,1,1)
       else
        do 25 i=1,ncsf
          ci_oav(i)=oav(i)
          oelocav(i)=0
          eav(i)=0
          do 25 j=1,ncsf
            oelocav(i)=oelocav(i)+oeav(i,j)*ccsf(j,1,1)
  25        eav(i)=eav(i)+oeav(j,i)*ccsf(j,1,1)
      endif

c Compute the gradient (also for first coefficient)
      do 30 i=1,nciterm
  30    grad_ci(i)=2*(oelocav(i)-etot*oav(i))

      if(method.eq.'hessian') then

c Compute the hessian (also for first coefficient)
        idx=0
        do 40 i=1,nciterm
          do 40 j=1,i
            idx=idx+1
            h_ci(i,j)=2*(ooeav(idx)-etot*ooav(idx)-oav(i)*grad_ci(j)-oav(j)*grad_ci(i))
     &               +oeav(i,j)+oeav(j,i)-2*ooeav(idx)
     &               +oav(i)*oelocav(j)+oav(j)*oelocav(i)-oav(i)*eav(j)-oav(j)*eav(i)
  40        h_ci(j,i)=h_ci(i,j)

        write(6,'(''opening grad_hess.ci.dat'')')
        open(43,file='grad_hess.ci.dat',status='unknown')

        write(43,*) nciterm
        write(43,*) (grad_ci(i),i=1,nciterm)
        write(43,*) ((h_ci(i,j),j=1,i),i=1,nciterm)

       else

        if(ioptjas.eq.0.and.ioptorb.eq.0) then
          is=0
          do 50 i=1,nciterm
            do 50 j=1,nciterm
  50          h_ci(i+is,j+is)=oeav(i,j)

          idx=0
          do 60 i=1,nciterm
            do 60 j=1,i
              idx=idx+1
              s_ci(i+is,j+is)=ooav(idx)
  60          s_ci(j+is,i+is)=ooav(idx)

         else
          is=1

c h_0,0, h_0,ci, h_ci,0, s_0,ci, s_ci,0
          h_ci(1,1)=etot
          s_ci(1,1)=1
          do 70 j=1,nciterm
            h_ci(1,j+is)=eav(j)-etot*oav(j)
            h_ci(j+is,1)=oelocav(j)-etot*oav(j)
            s_ci(1,j+is)=0
            s_ci(j+is,1)=0
            do 70 i=1,nciterm
  70          h_ci(i+is,j+is)=oeav(i,j)+etot*oav(i)*oav(j)-oav(i)*eav(j)-oav(j)*oelocav(i)

        idx=0
        do 80 i=1,nciterm
          do 80 j=1,i
            idx=idx+1
            s_ci(i+is,j+is)=ooav(idx)-oav(i)*oav(j)
  80        s_ci(j+is,i+is)=s_ci(i+is,j+is)

        endif

      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine optci_prt(w,iblk,iu)
      use linear_norm, only: oav
      use optwf_contrl, only: ioptci
      use ci000, only: iciprt, nciprim, nciterm

      implicit real*8(a-h,o-z)

c compute averages and print then out
      include 'vmc.h'
      include 'force.h'
      include 'optci.h'

      common /icount_ci/ icount_ci

      dimension deav(MXCITERM)
      dimension oeav(MXCITERM,MXCIREDUCED),oeerr(MXCITERM,MXCIREDUCED)
      dimension ooav(MXCIMATDIM),ooerr(MXCIMATDIM)
      dimension ooeav(MXCIMATDIM)
      dimension itemp_print(5), temp_print(5)

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      if(iciprt.eq.0) return

c iciprt 0 no printout
c         1 each iteration full printout
c        >1 after iciprt iterations reduced printout
c        -1 force printout

      if(iciprt.gt.0.and.icount_ci.ne.iciprt) then
        icount_ci=icount_ci+1
        return
      endif

      icount_ci=1

      call optci_avrg(w,iblk,oav,deav,oeav,oeerr,ooav,ooerr,ooeav)

c     print the Ok
      write (45,*) 
      write (45,'(''CI Operators'')')
      write (45,'(''------------'')')
      write (45,*) 

c     print the OkOl
      write (45,'(''Overlap'')')
      write (45,*) 
      idx=0
      do k=1,nciterm
       write (45,*)
       write (45,*) '< Ok O',k,' >'
       i=1
       do while (i.lt.min(k+1,nciterm))
        jmax=5
        if (i+jmax.gt.k) then
         jmax=mod(k-1,5)+1
        endif
        do j=1,jmax
          idx=idx+1
          temp_print(j) = ooav(idx)
          itemp_print(j) = int(1000000*ooerr(idx))
        enddo
        write (45,22) i, i+jmax-1, (temp_print(j), itemp_print(j), j=1,jmax)
        i=i+jmax
       enddo
      enddo

c     print the OkEL
      write (45,*) 
      write (45,'(''Hamiltonian'')')
      do k=1,nciterm
       write (45,*)
       write (45,*) '< Ok E',k,' >'
       i=1
       do while (i.lt.nciterm)
        jmax=5
        if (i+jmax.gt.nciterm) then
         jmax=nciterm-i+1
        endif
        do j=1,jmax
          temp_print(j) = oeav(i+j-1,k)
          itemp_print(j) = int(1000000*oeerr(i+j-1,k))
        enddo
        write (45,22) i, i+jmax-1, (temp_print(j), itemp_print(j), j=1,jmax)
        i=i+jmax
       enddo
      enddo

      write(iu,'(''Operators correspond to primitive set'')')

 22   format('k= ',i4,'->',i4,2X,5(f11.5,'(',i6.0,')',1X))

      end

c-----------------------------------------------------------------------
      subroutine optci_define

      use csfs, only: ncsf
      use dets, only: ndet
      use optwf_contrl, only: ioptjas, ioptorb
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      use ci000, only: iciprt, nciprim, nciterm

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'optci.h'
      include 'mstates.h'

      nciprim = ndet

      ici_def = 1
      if(ncsf.eq.0) then
        nciterm=nciprim
       else
        is=1
        if((method.eq.'linear'.or.method.eq.'lin_d').and.ioptjas.eq.0.and.ioptorb.eq.0) is=0
        nciterm=ncsf-is
      endif

      end
c-----------------------------------------------------------------------
