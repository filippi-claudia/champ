      subroutine optci_deloc(eloc_det,e_other,psid,energy)
      use vmc_mod, only: MDET
      use mstates_mod, only: MSTATES
      use csfs, only: cxdet, iadet, ibdet, icxdet, ncsf
      use optwf_contrl, only: ioptci
      use ci000, only: nciprim
      use ci001_blk, only: ci_o, ci_oe
      use ci003_blk, only: ci_e
      use ci004_blk, only: ci_de
      use method_opt, only: method
      use multislater, only: detiab
      use csfs, only: nstates

      implicit real*8(a-h,o-z)

      dimension ciprim(MDET,MSTATES),cieprim(MDET,MSTATES)
      dimension eloc_det(MDET,2),psid(MSTATES),energy(MSTATES)
      
      if(ioptci.eq.0) return 

      do istate=1,nstates
         psidi=1.0d0/psid(istate)
         do k=1,nciprim
            ciprim(k,istate)=detiab(k,1,istate)*detiab(k,2,istate)*psidi
            cieprim(k,istate)=(eloc_det(k,1)+eloc_det(k,2)+e_other)*ciprim(k,istate)
         enddo

c     Update <Oi>,<Ei>,<dEi>,<Oi*Ej> 
c     Correlation matrix <Oi*Oj> is computed in ci_sum

         if(ncsf.eq.0) then
            do i=1,nciprim
               ci_o(i,istate)=ciprim(i,istate)
               ci_e(i,istate)=cieprim(i,istate)
               ci_de(i,istate)=cieprim(i,istate)-ciprim(i,istate)*energy(istate)
            enddo
         else
            do icsf=1,ncsf
               ci_o_csf=0.0d0
               ci_e_csf=0.0d0
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

      if(method.eq.'sr_n'.or.method.eq.'lin_d') return

      do istate=1,nstates
         if(ncsf.eq.0) then
            do i=1,nciprim
               do j=1,nciprim
                  ci_oe(i,j,istate)=ciprim(i,istate)*cieprim(j,istate)
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

      end subroutine

c-----------------------------------------------------------------------

      subroutine optci_init(iflg)
      use optwf_contrl, only: ioptci
      use ci000, only: nciterm
      use ci005_blk, only: ci_o_cum, ci_o_sum
      use ci006_blk, only: ci_de_cum, ci_de_sum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum, ci_oe_sum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum, ci_oo_sum
      use ci010_blk, only: ci_ooe_cum, ci_ooe_sum
      use method_opt, only: method
      use csfs, only: nstates

      implicit real*8(a-h,o-z)

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do istate=1,nstates
         do i=1,nciterm
            ci_o_sum(i,istate)=0.0d0
            ci_de_sum(i,istate)=0.0d0
            do j=1,nciterm
               ci_oe_sum(j,i,istate)=0.0d0
            enddo
         enddo

         idx=0
         do i=1,nciterm
            do j=1,i
               idx=idx+1
               ci_oo_sum(idx,istate)=0.0d0
               ci_ooe_sum(idx,istate)=0.0d0
            enddo
         enddo
      enddo

C     $ iflg = 0: init *cum, *cm2 as well
      if(iflg.gt.0) return

      do istate=1,nstates
         guid_weight=0.0d0
         guid_weight_sq=0.0d0

         do i=1,nciterm
            ci_o_cum(i,istate)=0.0d0
            ci_de_cum(i,istate)=0.0d0
            do j=1,nciterm
               ci_oe_cum(j,i,istate)=0.0d0
               ci_oe_cm2(j,i,istate)=0.0d0
            enddo
         enddo

         idx=0
         do i=1,nciterm
            do j=1,i
               idx=idx+1
               ci_oo_cum(idx,istate)=0.0d0
               ci_oo_cm2(idx,istate)=0.0d0
               ci_ooe_cum(idx,istate)=0.0d0
            enddo
         enddo
      enddo

      end subroutine

c-----------------------------------------------------------------------

      subroutine optci_save
      use optwf_contrl, only: ioptci
      use ci000, only: nciterm
      use ci001_blk, only: ci_o, ci_oe
      use ci002_blk, only: ci_o_old, ci_oe_old
      use ci003_blk, only: ci_e, ci_e_old
      use ci004_blk, only: ci_de, ci_de_old
      use method_opt, only: method
      use csfs, only: nstates

      implicit real*8(a-h,o-z)

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do istate=1,nstates
         do i=1,nciterm
            ci_o_old(i,istate)=ci_o(i,istate)
            ci_e_old(i,istate)=ci_e(i,istate)
            ci_de_old(i,istate)=ci_de(i,istate)
            do j=1,nciterm
               ci_oe_old(j,i,istate)=ci_oe(j,i,istate)
            enddo
         enddo
      enddo

      end subroutine

c-----------------------------------------------------------------------

      subroutine optci_restore
      use optwf_contrl, only: ioptci
      use ci000, only: nciterm
      use ci001_blk, only: ci_o, ci_oe
      use ci002_blk, only: ci_o_old, ci_oe_old
      use ci003_blk, only: ci_e, ci_e_old
      use ci004_blk, only: ci_de, ci_de_old
      use method_opt, only: method
      use csfs, only: nstates

      implicit real*8(a-h,o-z)

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do istate=1,nstates
         do i=1,nciterm
            ci_o(i,istate)=ci_o_old(i,istate)
            ci_e(i,istate)=ci_e_old(i,istate)
            ci_de(i,istate)=ci_de_old(i,istate)
            do j=1,nciterm
               ci_oe(j,i,istate)=ci_oe_old(j,i,istate)
            enddo
         enddo
      enddo

      end subroutine

c-----------------------------------------------------------------------

      subroutine optci_sum(p,q,enew,eold)
      use optwf_contrl, only: ioptci
      use ci000, only: nciterm
      use ci001_blk, only: ci_o, ci_oe
      use ci002_blk, only: ci_o_old, ci_oe_old
      use ci004_blk, only: ci_de, ci_de_old
      use ci005_blk, only: ci_o_sum
      use ci006_blk, only: ci_de_sum
      use ci008_blk, only: ci_oe_sum
      use ci009_blk, only: ci_oo_sum
      use ci010_blk, only: ci_ooe_sum
      use method_opt, only: method
      use csfs, only: nstates

      implicit real*8(a-h,o-z)

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do istate=1,nstates
         do j=1,nciterm
            ci_o_sum(j,istate)=ci_o_sum(j,istate)+p*ci_o(j,istate)+q*ci_o_old(j,istate)
            ci_de_sum(j,istate)=ci_de_sum(j,istate)+p*ci_de(j,istate)+q*ci_de_old(j,istate)
            do i=1,nciterm
               ci_oe_sum(i,j,istate)=ci_oe_sum(i,j,istate)
     &              +p*ci_oe(i,j,istate)+q*ci_oe_old(i,j,istate)
            enddo
         enddo

         idx=0
         do i=1,nciterm
            do j=1,i
               idx=idx+1
               ci_oo_new=ci_o(i,istate)*ci_o(j,istate)
               ci_oo_old=ci_o_old(i,istate)*ci_o_old(j,istate)
               ci_oo_sum(idx,istate)=ci_oo_sum(idx,istate)+p*ci_oo_new+q*ci_oo_old
               ci_ooe_sum(idx,istate)=ci_ooe_sum(idx,istate)
     &              +p*ci_oo_new*enew+q*ci_oo_old*eold
            enddo
         enddo
      enddo

      end subroutine

c-----------------------------------------------------------------------

      subroutine optci_cum(wsum)
      use optwf_contrl, only: ioptci
      use ci000, only: nciterm
      use ci005_blk, only: ci_o_cum, ci_o_sum
      use ci006_blk, only: ci_de_cum, ci_de_sum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum, ci_oe_sum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum, ci_oo_sum
      use ci010_blk, only: ci_ooe_cum, ci_ooe_sum
      use method_opt, only: method
      use csfs, only: nstates

      implicit real*8(a-h,o-z)

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do istate=1,nstates
         do i=1,nciterm
            ci_o_cum(i,istate)=ci_o_cum(i,istate)+ci_o_sum(i,istate)
            ci_de_cum(i,istate)=ci_de_cum(i,istate)+ci_de_sum(i,istate)
            do j=1,nciterm
               ci_oe_now=ci_oe_sum(i,j,istate)/wsum
               ci_oe_cum(i,j,istate)=ci_oe_cum(i,j,istate)+ci_oe_sum(i,j,istate)
               ci_oe_cm2(i,j,istate)=ci_oe_cm2(i,j,istate)+ci_oe_sum(i,j,istate)*ci_oe_now
            enddo
         enddo

         idx=0
         do i=1,nciterm
            do j=1,i
               idx=idx+1
               ci_oo_cum(idx,istate)=ci_oo_cum(idx,istate)+ci_oo_sum(idx,istate)
               ci_oo_cm2(idx,istate)=ci_oo_cm2(idx,istate)
     &              +ci_oo_sum(idx,istate)*ci_oo_sum(idx,istate)/wsum
               ci_ooe_cum(idx,istate)=ci_ooe_cum(idx,istate)+ci_ooe_sum(idx,istate)
            enddo
         enddo
      enddo

      end subroutine

c-----------------------------------------------------------------------

      subroutine optci_dump(iu)
      use optwf_contrl, only: ioptci
      use ci000, only: nciprim, nciterm
      use ci005_blk, only: ci_o_cum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum
      use ci010_blk, only: ci_ooe_cum
      use method_opt, only: method

      implicit real*8(a-h,o-z)

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      istate=1

      matdim=nciterm*(nciterm+1)/2
      write(iu) nciprim,nciterm
      write(iu) (ci_o_cum(i,istate),i=1,nciterm)
      write(iu) ((ci_oe_cum(i,j,istate),ci_oe_cm2(i,j,istate),i=1,nciterm),j=1,nciterm)
      write(iu) (ci_oo_cum(i,istate),ci_oo_cm2(i,istate),i=1,matdim)
      write(iu) (ci_ooe_cum(i,istate),i=1,matdim)

      end subroutine

c-----------------------------------------------------------------------

      subroutine optci_rstrt(iu)
      use optwf_contrl, only: ioptci
      use ci000, only: nciprim, nciterm
      use ci005_blk, only: ci_o_cum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum
      use ci010_blk, only: ci_ooe_cum
      use method_opt, only: method

      implicit real*8(a-h,o-z)

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      istate=1

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

      read(iu) (ci_o_cum(i,istate),i=1,nciterm)
      read(iu) ((ci_oe_cum(i,j,istate),ci_oe_cm2(i,j,istate),i=1,nciterm),j=1,nciterm)
      matdim=nciterm*(nciterm+1)/2 
      read(iu) (ci_oo_cum(i,istate),ci_oo_cm2(i,istate),i=1,matdim)
      read(iu) (ci_ooe_cum(i,istate),i=1,matdim)

      end subroutine

c-----------------------------------------------------------------------

      subroutine optci_avrg(wcum,iblk,oav,deav,oeav,oeerr,ooav,ooerr,ooeav,istate)
      use optci, only: MXCITERM, MXCIREDUCED, MXCIMATDIM
      use optwf_contrl, only: ioptci
      use ci000, only: nciterm
      use ci005_blk, only: ci_o_cum
      use ci006_blk, only: ci_de_cum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum
      use ci010_blk, only: ci_ooe_cum
      use method_opt, only: method
      use csfs, only: nstates

      implicit real*8(a-h,o-z)

      dimension oav(MXCITERM),deav(MXCITERM)
      dimension oeav(MXCITERM,MXCIREDUCED),oeerr(MXCITERM,MXCIREDUCED)
      dimension ooav(MXCIMATDIM),ooerr(MXCIMATDIM)
      dimension ooeav(MXCIMATDIM)

      err(x,x2)=dsqrt(abs(x2/wcum-(x/wcum)**2)/iblk)

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do i=1,nciterm
         oav(i)=ci_o_cum(i,istate)/wcum
         deav(i)=ci_de_cum(i,istate)/wcum
         do j=1,nciterm
            oeav(i,j)=ci_oe_cum(i,j,istate)/wcum
         enddo

         do j=1,nciterm
            oeerr(i,j)=err(ci_oe_cum(i,j,istate),ci_oe_cm2(i,j,istate))
         enddo
      enddo

      idx=0
      do i=1,nciterm
         do j=1,i
            idx=idx+1
            ooav(idx)=ci_oo_cum(idx,istate)/wcum
            ooeav(idx)=ci_ooe_cum(idx,istate)/wcum
            ooerr(idx)=err(ci_oo_cum(idx,istate),ci_oo_cm2(idx,istate))
         enddo
      enddo

      end subroutine

c-----------------------------------------------------------------------

      subroutine optci_fin(iblk,passes,etot)
      use optci, only: MXCITERM, MXCIREDUCED, MXCIMATDIM
      use csfs, only: ccsf, ncsf
      use dets, only: cdet
      use gradhess_ci, only: grad_ci, h_ci, s_ci
      use linear_norm, only: oav, ci_oav
      use optwf_contrl, only: ioptci, ioptjas, ioptorb
      use ci000, only: iciprt, nciterm
      use method_opt, only: method

      implicit real*8(a-h,o-z)

      dimension deav(MXCITERM)
      dimension oeav(MXCITERM,MXCIREDUCED)
      dimension oeerr(MXCITERM,MXCIREDUCED)
      dimension ooav(MXCIMATDIM),ooerr(MXCIMATDIM)
      dimension ooeav(MXCIMATDIM)
      dimension oelocav(MXCITERM),eav(MXCITERM)

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

c     RLPB: Warning! At the moment this is done for only one state
      istate=1

      iciprt_sav=iciprt
      iciprt=-1
      call optci_avrg(passes,iblk,oav,deav,oeav,oeerr,ooav,ooerr,ooeav,istate)
      iciprt=iciprt_sav
      
      if(ncsf.eq.0) then
         do i=1,nciterm
            ci_oav(i)=oav(i)
            oelocav(i)=0.0d0
            eav(i)=0.0d0
            do j=1,nciterm
               oelocav(i)=oelocav(i)+oeav(i,j)*cdet(j,1,1)
               eav(i)=eav(i)+oeav(j,i)*cdet(j,1,1)
            enddo
         enddo
      else
         do i=1,ncsf
            ci_oav(i)=oav(i)
            oelocav(i)=0.0d0
            eav(i)=0.0d0
            do j=1,ncsf
               oelocav(i)=oelocav(i)+oeav(i,j)*ccsf(j,1,1)
               eav(i)=eav(i)+oeav(j,i)*ccsf(j,1,1)
            enddo
         enddo
      endif

c     Compute the gradient (also for first coefficient)
      do i=1,nciterm
         grad_ci(i)=2*(oelocav(i)-etot*oav(i))
      enddo

      if(method.eq.'hessian') then
c     Compute the hessian (also for first coefficient)
         idx=0
         do i=1,nciterm
            do j=1,i
               idx=idx+1
               h_ci(i,j)=2*(ooeav(idx)-etot*ooav(idx)
     &              -oav(i)*grad_ci(j)-oav(j)*grad_ci(i))
     &              +oeav(i,j)+oeav(j,i)-2*ooeav(idx)
     &              +oav(i)*oelocav(j)+oav(j)*oelocav(i)
     &              -oav(i)*eav(j)-oav(j)*eav(i)
               h_ci(j,i)=h_ci(i,j)
            enddo
         enddo

         write(6,'(''opening grad_hess.ci.dat'')')
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
c     h_0,0, h_0,ci, h_ci,0, s_0,ci, s_ci,0
            is=1
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

      end subroutine

c-----------------------------------------------------------------------

      subroutine optci_prt(w,iblk,iu)
c     compute averages and print them out
      use optci, only: MXCITERM, MXCIREDUCED, MXCIMATDIM
      use linear_norm, only: oav
      use optwf_contrl, only: ioptci
      use ci000, only: iciprt, nciterm
      use m_icount, only: icount_ci
      use method_opt, only: method

      implicit real*8(a-h,o-z)

      dimension deav(MXCITERM)
      dimension oeav(MXCITERM,MXCIREDUCED),oeerr(MXCITERM,MXCIREDUCED)
      dimension ooav(MXCIMATDIM),ooerr(MXCIMATDIM)
      dimension ooeav(MXCIMATDIM)
      dimension itemp_print(5), temp_print(5)

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      if(iciprt.eq.0) return

c     RLPB: Warning! At the moment this is done for only one state
      istate=1

c     iciprt 0 no printout
c     1 each iteration full printout
c     >1 after iciprt iterations reduced printout
c     -1 force printout

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

      end subroutine

c-----------------------------------------------------------------------

      subroutine optci_define
      use csfs, only: ncsf
      use dets, only: ndet
      use optwf_contrl, only: ioptjas, ioptorb
      use inputflags, only: ici_def
      use ci000, only: nciprim, nciterm
      use method_opt, only: method

      implicit real*8(a-h,o-z)

      nciprim = ndet

      ici_def = 1
      if(ncsf.eq.0) then
         nciterm=nciprim
      else
         is=1
         if((method.eq.'linear'.or.method.eq.'lin_d')
     &        .and.ioptjas.eq.0.and.ioptorb.eq.0) is=0
         nciterm=ncsf-is
      endif

      end subroutine
