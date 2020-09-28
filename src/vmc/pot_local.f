      subroutine pot_local(pe)
      use vmc_mod, only: MELEC, MCENT
      use vmc_mod, only: MMAT_DIM2
      use atom, only: znuc, pecent, iwctype, ncent
      use ghostatom, only: nghostcent
      use const, only: nelec, ipr
      use contrl_per, only: iperiodic
      use distance_mod, only: rshift, r_en, rvec_en
      use pseudo, only: nloc

      implicit real*8(a-h,o-z)



      ! common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)
      ! common /distance/ rshift(3,MELEC,MCENT), rvec_ee(3,MMAT_DIM2), r_ee(MMAT_DIM2)
      ! common /distance/ rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
      common /distance/ rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)
c  pe from nucleus-nucleus repulsion
      pe=pecent
      pe_ee=0.d0
      pe_en=0.d0
      if(iperiodic.eq.0) then
        do 40 i=1,nelec
          do 40 ic=1,ncent+nghostcent
   40       if(nloc.eq.0.and.ic.le.ncent) pe=pe-znuc(iwctype(ic))/r_en(i,ic)
        ij=0
        do 50 i=2,nelec
          do 50 j=1,i-1
            ij=ij+1
   50       pe=pe+1/r_ee(ij)
       else 
        call pot_en_ewald(x,pe_en)
        call pot_ee_ewald(x,pe_ee)
        pe=pe+pe_en+pe_ee
      endif 
      if(ipr.ge.3) write(6,'(''pe,pe_en(loc),pe_ee'',9f9.5)') pe,pe_en,pe_ee
        
      return
      end 
