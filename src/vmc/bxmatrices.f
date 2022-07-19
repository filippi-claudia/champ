      module bxmatrices
      contains
      subroutine bxmatrix(kref,xmatu,xmatd,b)

      use vmc_mod, only: norb_tot
      use elec, only: ndn, nup
      use dorb_m, only: iworbd
      use slater, only: slmi

      use precision_kinds, only: dp
      use multiply_slmi_mderiv, only: multiply_slmi_mderiv_simple
      use system, only: nelec

      implicit none

      integer :: i, iab, iel, ish, j
      integer :: kref, nel

      real(dp), dimension(norb_tot, nelec) :: b
      real(dp), dimension(nelec**2, 2) :: btemp
      real(dp), dimension(nelec**2) :: xmatu
      real(dp), dimension(nelec**2) :: xmatd
      real(dp), dimension(nelec) :: work





      do iab=1,2
        if(iab.eq.1) then
          iel=0
          nel=nup
         else
          iel=nup
          nel=ndn
        endif
        ish=-nel
        do i=1,nel
          ish=ish+nel
          do j=1,nel
            btemp(j+ish,iab)=b(iworbd(j+iel,kref),i+iel)
          enddo
        enddo
      enddo

      call multiply_slmi_mderiv_simple(nup,btemp(1,1),work,slmi(1,1),xmatu)
      call multiply_slmi_mderiv_simple(ndn,btemp(1,2),work,slmi(1,2),xmatd)

      return
      end
c-----------------------------------------------------------------------
      end module
