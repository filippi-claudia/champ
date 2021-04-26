      subroutine multideterminant_tmove(psid,iel_move)

      use vmc_mod, only: MELEC, MORB
      use const, only: nelec
      use atom, only: ncent
      use qua, only: nquad
      use b_tmove, only: b_t, iskip
      use casula, only: icasula, t_vpsp
      use slater, only: slmi
      use elec, only: ndn, nup
      use dorb_m, only: iworbd
      use coefs, only: norb
      use ycompact, only: ymat
      use multislater, only: detd, detu
      use multidet, only: iactv, ivirt, kref
      use multimat, only: aa

      use precision_kinds, only: dp
      implicit none

      integer :: i1, i2, iab, ic, iel
      integer :: iel_move, iq, irep, ish
      integer :: j, jel, jrep, nel
      real(dp) :: detratio, dum, psid
      real(dp), dimension(MELEC, MORB) :: gmat
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: half = 0.5d0



      if(icasula.gt.0)then
        i1=iel_move
        i2=iel_move
       else
        i1=1
        i2=nelec
      endif

      do iel=i1,i2

      do ic=1,ncent
        
      if(iskip(iel,ic).eq.0) then

      if(iel.le.nup) then
        iab=1
        nel=nup
        ish=0
       else
        iab=2
        nel=ndn
        ish=nup
      endif

      detratio=detu(kref)*detd(kref)/psid

      jel=iel-ish

      do iq=1,nquad

        do jrep=ivirt(iab),norb
          dum=0
          do j=1,nel
            dum=dum+b_t(iworbd(j+ish,kref),iq,ic,iel)*aa(j,jrep,iab)
          enddo
          dum=b_t(jrep,iq,ic,iel)-dum

          do irep=iactv(iab),nel
            gmat(irep,jrep)=dum*slmi(irep+(jel-1)*nel,iab)
          enddo
        enddo

c     t_vpsp(ic,iq,iel)=t_vpsp_ref

      dum=0
      do jrep=ivirt(iab),norb
        do irep=iactv(iab),nel
          dum=dum+ymat(jrep,irep,iab,1)*gmat(irep,jrep)
        enddo
      enddo
      t_vpsp(ic,iq,iel)=t_vpsp(ic,iq,iel)+dum*detratio

      enddo

      endif

      enddo
      enddo

      return
      end
