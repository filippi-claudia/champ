module fetch_parameters_mod
contains
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
      subroutine fetch_jastrow(p,n)
      use jastrow, only: a4,b,c
      use optwf_control, only: ioptjas
      use optwf_nparmj, only: nparma,nparmb,nparmc
      use optwf_wjas, only: iwjasa,iwjasb,iwjasc
      use precision_kinds, only: dp
      use system,  only: nctype
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
!-----------------------------------------------------------------------
      subroutine fetch_lcao(p,n)

      use optorb_cblock, only: norbterm
      use optwf_control, only: ioptorb
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
!-----------------------------------------------------------------------
      subroutine fetch_ci(p,n)
      use contrl_file, only: ounit
      use csfs, only: ccsf,maxcsf,ncsf
      use optwf_control, only: ioptci
      use precision_kinds, only: dp
      use slater,  only: cdet,ndet
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
        do icsf=1,maxcsf(1)-1
          p(icsf)=ccsf(icsf,1,1)
        enddo
        do icsf=maxcsf(1)+1,ncsf
          p(icsf-1)=ccsf(icsf,1,1)
        enddo
        n=ncsf-1
      endif

!     do 90 j=1,nstates
!90     write(ounit,'(''csf ='',1000f20.15)') (ccsf(i,j,iadiag),i=1,ncsf)

      return
      end
!-----------------------------------------------------------------------
end module
