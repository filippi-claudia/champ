      subroutine psie(iel,coord,psid,psij,ipass,iflag)
c Written by Claudia Filippi by modifying hpsi

      use csfs, only: nstates
      use estpsi, only: apsi, aref
      use multidet, only: kref
      use wfsec, only: iwf, iwftype
      use contr2, only: ianalyt_lap
      use velocity_jastrow, only: vjn

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'pseudo.h'
      include 'force.h'
      include 'optjas.h'
      include 'optci.h'

c Calculates wave function

      common /multislatern/ detn(MDET)
     &,orbn(MORB),dorbn(3,MORB),ddorbn(MORB)
      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)

      dimension coord(3,*),psid(*)

      iwf=iwftype(1)

      call distances(iel,coord)

      if(ianalyt_lap.eq.1) then
        call jastrowe(iel,coord,vjn,d2j,psij,iflag)
       else
        call fatal_error('HPSIE: numerical one-electron move not implemented')
      endif

c compute all determinants 
      call determinante(iel,x,rvec_en,r_en,iflag)

      if(detn(kref).eq.0.d0) then
        do 1 istate=1,nstates
   1      psid(istate)=0.d0
        return
      endif

      call multideterminante(iel)

c combine determinantal quantities to obtain trial wave function
      do 10 istate=1,nstates
   10   call determinante_psit(iel,psid(istate),istate)

      if(ipass.gt.2) then

        check_apsi_min=1.d+99
        do 20 istate=1,nstates
          apsi_now=apsi(istate)/(ipass-1)
          check_apsi=abs(psid(istate))/apsi_now

   20    check_apsi_min=min(check_apsi,check_apsi_min)

        aref_now=aref/(ipass-1)
        check_dref=abs(detn(kref))/aref_now

      endif

      return
      end
