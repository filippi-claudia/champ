      subroutine fetch_parameters(p)
c     this was not in master but I think it's needed
c     this is so confusng ... 
      use optwf_contrl, only: nparm
      
      implicit real*8(a-h,o-z)

      dimension p(*)

      n=0
      ip=1
      call fetch_jastrow(p(1),n)
      nparm=n

      call fetch_lcao(p(nparm+1),n)
      nparm=nparm+n

      call fetch_ci(p(nparm+1),n)
      nparm=nparm+n

      end subroutine

c-----------------------------------------------------------------------

      subroutine fetch_jastrow(p,n)
      use atom, only: nctype
      use jaspar3, only: a, b, c
      use jaspar4, only: a4
      use optwf_contrl, only: ioptjas
      use optwf_nparmj, only: nparma, nparmb, nparmc
      use optwf_wjas, only: iwjasa, iwjasb, iwjasc

      implicit real*8(a-h,o-z)

      dimension p(*)

      if(ioptjas.eq.0) return

      iparm=0
      do ict=1,nctype
         do i=1,nparma(ict)
            iparm=iparm+1
            p(iparm)=a4(iwjasa(i,ict),ict,1)
         enddo
      enddo

      do i=1,nparmb(1)
         iparm=iparm+1
         p(iparm)=b(iwjasb(i,1),1,1)
      enddo

      do ict=1,nctype
         do i=1,nparmc(ict)
            iparm=iparm+1
            p(iparm)=c(iwjasc(i,ict),ict,1)
         enddo
      enddo

      n=iparm

      end subroutine

c-----------------------------------------------------------------------

      subroutine fetch_lcao(p,n)
      use optwf_contrl, only: ioptorb
      use optorb_cblock, only: norbterm

      implicit real*8(a-h,o-z)

      dimension p(*)

      if(ioptorb.eq.0) return

      do i=1,norbterm
         p(i)=0.d0
      enddo
      n=norbterm

      end subroutine

c-----------------------------------------------------------------------

      subroutine fetch_ci(p,n)
      use csfs, only: ccsf, ncsf, nstates
      use dets, only: cdet, ndet
      use optwf_contrl, only: ioptci

      implicit real*8(a-h,o-z)

      dimension p(*)

      if(ioptci.eq.0) return

      if(ncsf.eq.0) then
         do idet=2,ndet
            p(idet-1)=cdet(idet,1,1)
         enddo
         n=ndet-1
      else
         do icsf=2,ncsf
            p(icsf-1)=ccsf(icsf,1,1)
         enddo
         n=ncsf-1
      endif

      end subroutine
