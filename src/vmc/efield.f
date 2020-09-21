      subroutine efield_extpot_ene(coord,nelec,efield_pot)

      use efield_mod, only: MCHARGES
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use efield_blk, only: ascreen, bscreen, qcharge, xcharge, ycharge, zcharge

      use efield, only: iscreen, ncharges

      implicit real*8(a-h,o-z)



      dimension coord(3,*)

      atmp=0
      efield_pot=0
      do 20 i=1,ncharges
         do 20 j=1,nelec
           rij2=(xcharge(i)-coord(1,j))**2+(ycharge(i)-coord(2,j))**2+(zcharge(i)-coord(3,j))**2
           rij=sqrt(rij2)
           if(iscreen.gt.0) atmp=ascreen(i)*dexp(-bscreen(i)*rij**2)
  20       efield_pot=efield_pot-qcharge(i)*(1.d0-atmp)/rij

      return
      end
c-----------------------------------------------------------------------
      subroutine efield_compute_extint

      use efield_mod, only: MCHARGES
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use atom, only: znuc, cent, iwctype, ncent
      use efield_blk, only: ascreen, bscreen, qcharge, xcharge, ycharge, zcharge

      use efield, only: iscreen, ncharges

      implicit real*8(a-h,o-z)




      efield_nucext=0
      atmp=0
      do 20 i=1,ncharges
c         write(6,*) 'HELLO',i,ascreen(i),bscreen(i)
         do 20 j=1,ncent
           rij2=(xcharge(i)-cent(1,j))**2+(ycharge(i)-cent(2,j))**2+(zcharge(i)-cent(3,j))**2
           rij=sqrt(rij2)
           if(iscreen.gt.0) atmp=ascreen(i)*dexp(-bscreen(i)*rij**2)
  20       efield_nucext=efield_nucext+znuc(iwctype(j))*qcharge(i)*(1.d0-atmp)/rij

      efield_extext=0
      do 30 i=1,ncharges
         do 30 j=i+1,ncharges
           rij2=(xcharge(i)-xcharge(j))**2+(ycharge(i)-ycharge(j))**2+(zcharge(i)-zcharge(j))**2
           rij=sqrt(rij2)
  30       efield_extext=efield_extext+qcharge(i)*qcharge(j)/rij

      write(6,'(''Number of charges: '',i5)') ncharges
      write(6,'(''Number of centers: '',i5)') ncent
     
c      write(6,*) 'znuc=',znuc(1:8)

      write(6,'(''Nuclear-external interaction: '',f15.9)') efield_nucext
      write(6,'(''External-external interaction:'',f15.9)') efield_extext
      return
      end
c-----------------------------------------------------------------------
