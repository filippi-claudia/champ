      module optci_mod
      use error, only: fatal_error
      contains
      subroutine optci_deloc(eloc_det,e_other,psid,energy)

      use slater, only: ndet
      use csfs, only: cxdet, iadet, ibdet, icxdet, ncsf, nstates
      use optwf_control, only: ioptci
      use ci000, only: nciprim
      use ci001_blk, only: ci_o, ci_oe
      use ci003_blk, only: ci_e
      use ci004_blk, only: ci_de
      use vmc_mod, only: nwftypeorb

      use optwf_control, only: method

      use multislater, only: detiab
      use precision_kinds, only: dp
      use contrl_file,    only: ounit

      implicit none

      integer :: i, icsf, idet, ix, j
      integer :: jcsf, k, istate
      real(dp) :: ci_e_csf, ci_o_csf, e_other
      real(dp) :: psidi
      real(dp), dimension(ndet,nwftypeorb) :: ciprim
      real(dp), dimension(ndet,nwftypeorb) :: cieprim
      real(dp), dimension(ndet, 2) :: eloc_det
      real(dp), dimension(nstates) :: psid
      real(dp), dimension(nstates) :: energy






      if(ioptci.eq.0) return

      do istate=1,nstates !STU check state mapping
        psidi=1.d0/psid(istate)

        do k=1,nciprim
          ciprim(k,istate)=detiab(k,1,istate)*detiab(k,2,istate)*psidi
          cieprim(k,istate)=(eloc_det(k,1)+eloc_det(k,2)+e_other)*ciprim(k,istate)
        enddo

c Update <Oi>,<Ei>,<dEi>,<Oi*Ej>
c Correlation matrix <Oi*Oj> is computed in ci_sum

        if(ncsf.eq.0) then
          do i=1,nciprim
            ci_o(i,istate)=ciprim(i,istate)
            ci_e(i,istate)=cieprim(i,istate)
            ci_de(i,istate)=cieprim(i,istate)-ciprim(i,istate)*energy(istate)
          enddo
        else
          do icsf=1,ncsf
            ci_o_csf=0
            ci_e_csf=0
            do ix=iadet(icsf),ibdet(icsf)
              idet=icxdet(ix)
              ci_o_csf=ci_o_csf+ciprim(idet,istate)*cxdet(ix)
              ci_e_csf=ci_e_csf+cieprim(idet,istate)*cxdet(ix)
            enddo
            ci_o(icsf,istate)=ci_o_csf
            ci_e(icsf,istate)=ci_e_csf
            ci_de(icsf,istate)=ci_e_csf-ci_o_csf*energy(istate)
          enddo

        endif
      enddo

      if(method.eq.'sr_n'.or.method.eq.'lin_d') return !STU this line requires 2 local variables and 2 global to be given an extra index.

      do istate=1,nstates !STU also here
        if(ncsf.eq.0) then
          do i=1,nciprim
            do j=1,nciprim
              ci_oe(i,j,istate) = ciprim(i,istate)*cieprim(j,istate)
            enddo
          enddo
        else
          do icsf=1,ncsf
            do jcsf=1,ncsf
              ci_oe(icsf,jcsf,istate)=ci_o(icsf,istate)*ci_e(jcsf,istate)
            enddo
          enddo
        endif
      enddo

      end
c-----------------------------------------------------------------------
      subroutine optci_init(iflg)

      use optwf_control, only: ioptci
      use ci000, only: nciterm
      use ci005_blk, only: ci_o_cum, ci_o_sum
      use ci006_blk, only: ci_de_cum, ci_de_sum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum, ci_oe_sum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum, ci_oo_sum
      use ci010_blk, only: ci_ooe_cum, ci_ooe_sum
      use vmc_mod, only: nwftypeorb

      use optwf_control, only: method

      use precision_kinds, only: dp
      implicit none

      integer :: i, idx, iflg, j, k
      real(dp) :: guid_weight, guid_weight_sq



      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do k=1,nwftypeorb
        do i=1,nciterm
          ci_o_sum(i,k) =0.d0
          ci_de_sum(i,k) =0.d0
          do j=1,nciterm
            ci_oe_sum(j,i,k) =0.d0
          enddo
        enddo

        idx=0
        do i=1,nciterm
          do j=1,i
            idx=idx+1
            ci_oo_sum(idx,k)=0.d0
            ci_ooe_sum(idx,k)=0.d0
          enddo
        enddo
      enddo

C$ iflg = 0: init *cum, *cm2 as well
      if(iflg.gt.0) return
      
      do k=1,nwftypeorb
        guid_weight=0.d0
        guid_weight_sq=0.d0
        do i=1,nciterm
          ci_o_cum(i,k)=0.d0
          ci_de_cum(i,k)=0.d0
          do j=1,nciterm
            ci_oe_cum(j,i,k)=0.d0
            ci_oe_cm2(j,i,k)=0.d0
          enddo
        enddo

        idx=0
        do i=1,nciterm
          do j=1,i
            idx=idx+1
            ci_oo_cum(idx,k)=0.d0
            ci_oo_cm2(idx,k)=0.d0
            ci_ooe_cum(idx,k)=0.d0
          enddo
        enddo
      enddo

      end

c-----------------------------------------------------------------------
      subroutine optci_save

      use optwf_control, only: ioptci
      use ci000, only: nciterm
      use ci001_blk, only: ci_o, ci_oe
      use ci002_blk, only: ci_o_old, ci_oe_old
      use ci003_blk, only: ci_e, ci_e_old
      use ci004_blk, only: ci_de, ci_de_old
      use vmc_mod, only: nwftypeorb

      use optwf_control, only: method

      implicit none

      integer :: i, j, k




      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do k=1, nwftypeorb
        do i=1,nciterm
          ci_o_old(i,k)=ci_o(i,k)
          ci_e_old(i,k)=ci_e(i,k)
          ci_de_old(i,k)=ci_de(i,k)
          do j=1,nciterm
            ci_oe_old(j,i,k)=ci_oe(j,i,k)
          enddo
        enddo
      enddo

      end
c-----------------------------------------------------------------------
      subroutine optci_restore

      use optwf_control, only: ioptci
      use ci000, only: nciterm
      use ci001_blk, only: ci_o, ci_oe
      use ci002_blk, only: ci_o_old, ci_oe_old
      use ci003_blk, only: ci_e, ci_e_old
      use ci004_blk, only: ci_de, ci_de_old
      use vmc_mod, only: nwftypeorb

      use optwf_control, only: method

      implicit none

      integer :: i, j, k




      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do k=1,nwftypeorb
        do i=1,nciterm
          ci_o(i,k)=ci_o_old(i,k)
          ci_e(i,k)=ci_e_old(i,k)
          ci_de(i,k)=ci_de_old(i,k)
          do j=1,nciterm
            ci_oe(j,i,k)=ci_oe_old(j,i,k)
          enddo
        enddo
      enddo

      end
c-----------------------------------------------------------------------
      subroutine optci_sum(p,q,enew,eold)

      use optwf_control, only: ioptci

      use ci000, only: nciterm
      use ci001_blk, only: ci_o, ci_oe
      use ci002_blk, only: ci_o_old, ci_oe_old
      use ci004_blk, only: ci_de, ci_de_old
      use ci005_blk, only: ci_o_sum
      use ci006_blk, only: ci_de_sum
      use ci008_blk, only: ci_oe_sum
      use ci009_blk, only: ci_oo_sum
      use ci010_blk, only: ci_ooe_sum
      use vmc_mod, only: nwftypeorb

      use optwf_control, only: method

      use precision_kinds, only: dp
      implicit none

      integer :: i, idx, j, k
      real(dp) :: ci_oo_new, ci_oo_old, enew, eold, p
      real(dp) :: q



      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do k=1,nwftypeorb !STU orb or max? check other subroutines in here
        do j=1,nciterm
          ci_o_sum(j,k) =ci_o_sum(j,k)+p*ci_o(j,k)+q*ci_o_old(j,k)
          ci_de_sum(j,k) =ci_de_sum(j,k)+p*ci_de(j,k)+q*ci_de_old(j,k)
          do i=1,nciterm
            ci_oe_sum(i,j,k)=ci_oe_sum(i,j,k)+p*ci_oe(i,j,k)+q*ci_oe_old(i,j,k)
          enddo
        enddo

        idx=0
        do i=1,nciterm
          do j=1,i
            idx=idx+1
            ci_oo_new=ci_o(i,k)*ci_o(j,k)
            ci_oo_old=ci_o_old(i,k)*ci_o_old(j,k)
            ci_oo_sum(idx,k)=ci_oo_sum(idx,k)  +p*ci_oo_new+q*ci_oo_old
            ci_ooe_sum(idx,k)=ci_ooe_sum(idx,k)+p*ci_oo_new*enew+q*ci_oo_old*eold
          enddo
        enddo
      enddo

      end
c-----------------------------------------------------------------------
      subroutine optci_cum(wsum)

      use optwf_control, only: ioptci

      use ci000, only: nciterm

      use ci005_blk, only: ci_o_cum, ci_o_sum
      use ci006_blk, only: ci_de_cum, ci_de_sum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum, ci_oe_sum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum, ci_oo_sum
      use ci010_blk, only: ci_ooe_cum, ci_ooe_sum

      use optwf_control, only: method
      use vmc_mod, only: nwftypeorb

      use precision_kinds, only: dp
      implicit none

      integer :: i, idx, j, k
      real(dp) :: ci_oe_now, wsum



      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do k=1,nwftypeorb
        do i=1,nciterm
          ci_o_cum(i,k)=ci_o_cum(i,k)+ci_o_sum(i,k)
          ci_de_cum(i,k)=ci_de_cum(i,k)+ci_de_sum(i,k)
          do j=1,nciterm
            ci_oe_now=ci_oe_sum(i,j,k)/wsum
            ci_oe_cum(i,j,k)=ci_oe_cum(i,j,k)+ci_oe_sum(i,j,k)
            ci_oe_cm2(i,j,k)=ci_oe_cm2(i,j,k)+ci_oe_sum(i,j,k)*ci_oe_now
          enddo
        enddo

        idx=0
        do i=1,nciterm
          do j=1,i
            idx=idx+1
            ci_oo_cum(idx,k)=ci_oo_cum(idx,k)+ci_oo_sum(idx,k)
            ci_oo_cm2(idx,k)=ci_oo_cm2(idx,k)+ci_oo_sum(idx,k)*ci_oo_sum(idx,k)/wsum
            ci_ooe_cum(idx,k)=ci_ooe_cum(idx,k)+ci_ooe_sum(idx,k)
          enddo
        enddo
      enddo

      end
c-----------------------------------------------------------------------
      subroutine optci_dump(iu)

      use optwf_control, only: ioptci
      use ci000, only: nciprim, nciterm
      use ci005_blk, only: ci_o_cum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum
      use ci010_blk, only: ci_ooe_cum

      use optwf_control, only: method

      implicit none

      integer :: i, iu, j, matdim, k




      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return
      
      k=1 !STU decode how to print multiple state info later

      matdim=nciterm*(nciterm+1)/2
      write(iu) nciprim,nciterm
      write(iu) (ci_o_cum(i,k),i=1,nciterm)
      write(iu) ((ci_oe_cum(i,j,k),ci_oe_cm2(i,j,k),i=1,nciterm),j=1,nciterm)
      write(iu) (ci_oo_cum(i,k),ci_oo_cm2(i,k),i=1,matdim)
      write(iu) (ci_ooe_cum(i,k),i=1,matdim)

      end
c-----------------------------------------------------------------------
      subroutine optci_rstrt(iu)
      use optwf_control, only: ioptci
      use ci000, only: nciprim, nciterm

      use ci005_blk, only: ci_o_cum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum
      use ci010_blk, only: ci_ooe_cum

      use optwf_control, only: method
      use contrl_file,    only: ounit
      implicit none

      integer :: i, iu, j, matdim, mciprim, k
      integer :: mciterm

      


      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      k=1 !STU decide how to setup restart file with multiple state

      read(iu) mciprim,mciterm
      if(mciprim.ne.nciprim) then
       write(ounit,*) 'wrong number of primitive ci terms!'
       write(ounit,*) 'old ',mciprim,' new ',nciprim
       call fatal_error('CI: Restart, inconsistent CI information')
      endif
      if(mciterm.ne.nciterm) then
       write(ounit,*) 'wrong number of ci terms!'
       write(ounit,*) 'old ',mciterm,' new ',nciterm
       call fatal_error('CI: Restart, inconsistent CI information')
      endif

      read(iu) (ci_o_cum(i,k),i=1,nciterm)
      read(iu) ((ci_oe_cum(i,j,k),ci_oe_cm2(i,j,k),i=1,nciterm),j=1,nciterm)
      matdim=nciterm*(nciterm+1)/2
      read(iu) (ci_oo_cum(i,k),ci_oo_cm2(i,k),i=1,matdim)
      read(iu) (ci_ooe_cum(i,k),i=1,matdim)

      end
c-----------------------------------------------------------------------
      subroutine optci_avrg(wcum,iblk,oav,deav,oeav,oeerr,ooav,ooerr,ooeav,k)

      use optci, only: mxciterm, mxcireduced, ncimatdim
      use optwf_control, only: ioptci
      use ci000, only: nciterm
      use ci005_blk, only: ci_o_cum
      use ci006_blk, only: ci_de_cum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum
      use ci010_blk, only: ci_ooe_cum

      use optwf_control, only: method

      use precision_kinds, only: dp
      implicit none

      integer :: i, iblk, idx, j, k
      real(dp) :: err, wcum, x, x2
      real(dp), dimension(mxciterm) :: oav
      real(dp), dimension(mxciterm) :: deav
      real(dp), dimension(mxciterm, mxcireduced) :: oeav
      real(dp), dimension(mxciterm, mxcireduced) :: oeerr
      real(dp), dimension(ncimatdim) :: ooav
      real(dp), dimension(ncimatdim) :: ooerr
      real(dp), dimension(ncimatdim) :: ooeav




      err(x,x2)=dsqrt(abs(x2/wcum-(x/wcum)**2)/iblk)

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do i=1,nciterm
       oav(i)=ci_o_cum(i,k)/wcum
       deav(i)=ci_de_cum(i,k)/wcum
       do j=1,nciterm
        oeav(i,j)=ci_oe_cum(i,j,k)/wcum
       enddo

       do j=1,nciterm
        oeerr(i,j)=err(ci_oe_cum(i,j,k),ci_oe_cm2(i,j,k))
       enddo
      enddo

      idx=0
      do i=1,nciterm
       do j=1,i
        idx=idx+1
        ooav(idx)=ci_oo_cum(idx,k)/wcum
        ooeav(idx)=ci_ooe_cum(idx,k)/wcum
        ooerr(idx)=err(ci_oo_cum(idx,k),ci_oo_cm2(idx,k))
       enddo
      enddo

      end
c-----------------------------------------------------------------------
      subroutine optci_fin(iblk,passes,etot)

      use optci, only: mxciterm, mxcireduced, ncimatdim
      use csfs, only: ccsf, ncsf
      use slater, only: cdet
      use gradhess_ci, only: grad_ci, h_ci, s_ci
      use linear_norm, only: ci_oav
      use optwf_control, only: ioptci, ioptjas, ioptorb
      use ci000, only: iciprt, nciterm

      use optwf_control, only: method
      use contrl_file,    only: ounit
      use precision_kinds, only: dp
      implicit none

      integer :: i, iblk, iciprt_sav, idx, is
      integer :: j, k
      real(dp) :: etot, passes
      real(dp), dimension(mxciterm) :: deav
      real(dp), dimension(mxciterm, mxcireduced) :: oeav
      real(dp), dimension(mxciterm, mxcireduced) :: oeerr
      real(dp), dimension(ncimatdim) :: ooav
      real(dp), dimension(ncimatdim) :: ooerr
      real(dp), dimension(ncimatdim) :: ooeav
      real(dp), dimension(mxciterm) :: oelocav
      real(dp), dimension(mxciterm) :: eav
      real(dp), dimension(mxciterm) :: oav




      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      k=1 !STU not sure what type of calculation calls this

      iciprt_sav=iciprt
      iciprt=-1
      call optci_avrg(passes,iblk,oav,deav,oeav,oeerr,ooav,ooerr,ooeav,k)
      iciprt=iciprt_sav

      if(ncsf.eq.0) then
        do i=1,nciterm
          ci_oav(i)=oav(i)
          oelocav(i)=0
          eav(i)=0
          do j=1,nciterm
            oelocav(i)=oelocav(i)+oeav(i,j)*cdet(j,1,1)
            eav(i)=eav(i)+oeav(j,i)*cdet(j,1,1)
          enddo
        enddo
       else
        do i=1,ncsf
          ci_oav(i)=oav(i)
          oelocav(i)=0
          eav(i)=0
          do j=1,ncsf
            oelocav(i)=oelocav(i)+oeav(i,j)*ccsf(j,1,1)
            eav(i)=eav(i)+oeav(j,i)*ccsf(j,1,1)
          enddo
        enddo
      endif

c Compute the gradient (also for first coefficient)
      do i=1,nciterm
        grad_ci(i)=2*(oelocav(i)-etot*oav(i))
      enddo

      if(method.eq.'hessian') then

c Compute the hessian (also for first coefficient)
        idx=0
        do i=1,nciterm
          do j=1,i
            idx=idx+1
            h_ci(i,j)=2*(ooeav(idx)-etot*ooav(idx)-oav(i)*grad_ci(j)-oav(j)*grad_ci(i))
     &               +oeav(i,j)+oeav(j,i)-2*ooeav(idx)
     &               +oav(i)*oelocav(j)+oav(j)*oelocav(i)-oav(i)*eav(j)-oav(j)*eav(i)
            h_ci(j,i)=h_ci(i,j)
          enddo
        enddo

        write(ounit,'(''opening grad_hess.ci.dat'')')
        open(43,file='grad_hess.ci.dat',status='unknown')

        write(43,*) nciterm
        write(43,*) (grad_ci(i),i=1,nciterm)
        write(43,*) ((h_ci(i,j),j=1,i),i=1,nciterm)

       else

        if(ioptjas.eq.0.and.ioptorb.eq.0) then
          is=0
          do i=1,nciterm
            do j=1,nciterm
              h_ci(i+is,j+is)=oeav(i,j)
            enddo
          enddo

          idx=0
          do i=1,nciterm
            do j=1,i
              idx=idx+1
              s_ci(i+is,j+is)=ooav(idx)
              s_ci(j+is,i+is)=ooav(idx)
            enddo
          enddo

         else
          is=1

c h_0,0, h_0,ci, h_ci,0, s_0,ci, s_ci,0
          h_ci(1,1)=etot
          s_ci(1,1)=1
          do j=1,nciterm
            h_ci(1,j+is)=eav(j)-etot*oav(j)
            h_ci(j+is,1)=oelocav(j)-etot*oav(j)
            s_ci(1,j+is)=0
            s_ci(j+is,1)=0
            do i=1,nciterm
              h_ci(i+is,j+is)=oeav(i,j)+etot*oav(i)*oav(j)-oav(i)*eav(j)-oav(j)*oelocav(i)
            enddo
          enddo

        idx=0
        do i=1,nciterm
          do j=1,i
            idx=idx+1
            s_ci(i+is,j+is)=ooav(idx)-oav(i)*oav(j)
            s_ci(j+is,i+is)=s_ci(i+is,j+is)
          enddo
        enddo

        endif

      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine optci_prt(w,iblk,iu)
      use optci, only: mxciterm, mxcireduced, ncimatdim
      use optwf_control, only: ioptci
      use ci000, only: iciprt, nciterm
      use m_icount, only: icount_ci
      use optwf_control, only: method

      use precision_kinds, only: dp
      implicit none

      integer :: i, iblk, idx, iu, j
      integer :: jmax, k, istate
      integer, dimension(5) :: itemp_print
      real(dp) :: w
      real(dp), dimension(mxciterm) :: deav
      real(dp), dimension(mxciterm) :: oav
      real(dp), dimension(mxciterm, mxcireduced) :: oeav
      real(dp), dimension(mxciterm, mxcireduced) :: oeerr
      real(dp), dimension(ncimatdim) :: ooav
      real(dp), dimension(ncimatdim) :: ooerr
      real(dp), dimension(ncimatdim) :: ooeav
      real(dp), dimension(5) :: temp_print


c compute averages and print then out



      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      if(iciprt.eq.0) return

      istate=1 !STU when is this called? will change

c iciprt 0 no printout
c         1 each iteration full printout
c        >1 after iciprt iterations reduced printout
c        -1 force printout

      if(iciprt.gt.0.and.icount_ci.ne.iciprt) then
        icount_ci=icount_ci+1
        return
      endif

      icount_ci=1

      call optci_avrg(w,iblk,oav,deav,oeav,oeerr,ooav,ooerr,ooeav,istate)

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
      use slater, only: ndet
      use optwf_control, only: ioptjas, ioptorb
      use inputflags, only: ici_def

      use ci000, only: nciprim, nciterm

      use optwf_control, only: method

      implicit none

      integer :: is




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
      end module
