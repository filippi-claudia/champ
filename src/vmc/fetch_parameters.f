      module fetch_parameters_mod
      contains
c-----------------------------------------------------------------------
      subroutine fetch_parameters(p)

      ! this was not in master but I think it's needed
      ! this is so confusng ...
      use optwf_control, only: nparm

      use precision_kinds, only: dp
      implicit none

      integer :: ip, n

      real(dp), dimension(*) :: p



      n=0

      ip=1
      call fetch_jastrow(p(1),n)
      nparm=n

      call fetch_lcao(p(nparm+1),n)
      nparm=nparm+n

      call fetch_ci(p(nparm+1),n)
      nparm=nparm+n

      return
      end
c-----------------------------------------------------------------------
      subroutine fetch_jastrow(p,n)
      use system, only: nctype
      use jaspar3, only: b, c

      use jaspar4, only: a4
      use optwf_control, only: ioptjas
      use optwf_nparmj, only: nparma, nparmb, nparmc
      use optwf_wjas, only: iwjasa, iwjasb, iwjasc
      use precision_kinds, only: dp
      implicit none

      integer :: i, ict, iparm, n

      real(dp), dimension(*) :: p



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

      return
      end
c-----------------------------------------------------------------------
      subroutine fetch_lcao(p,n)

      use optwf_control, only: ioptorb
      use optorb_cblock, only: norbterm

      use precision_kinds, only: dp
      implicit none

      integer :: i, n

      real(dp), dimension(*) :: p





      if(ioptorb.eq.0) return

      do i=1,norbterm
       p(i)=0.d0
      enddo
      n=norbterm

      return
      end
c-----------------------------------------------------------------------
      subroutine fetch_ci(p,n)
      use csfs, only: ccsf, ncsf
      use contrl_file, only: ounit
      use dets, only: cdet, ndet
      use optwf_control, only: ioptci
      use precision_kinds, only: dp
      implicit none

      integer :: i, iadiag, icsf, idet, j
      integer :: n
      real(dp) :: c90
      real(dp), dimension(*) :: p








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

c     do 90 j=1,nstates
c90     write(ounit,'(''csf ='',1000f20.15)') (ccsf(i,j,iadiag),i=1,ncsf)

      return
      end
c-----------------------------------------------------------------------
      end module
