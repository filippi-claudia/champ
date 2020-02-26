      subroutine distances(iel,x)
c Written by Cyrus Umrigar
c calculate interparticle distances
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use ghostatom, only: newghostype, nghostcent
      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use distances_sav, only: r_ee_sav, r_en_sav, rshift_sav, rvec_ee_sav, rvec_en_sav
      implicit real*8(a-h,o-z)




      include 'vmc.h'
      include 'force.h'
      include 'pseudo.h'

      common /contrl_per/ iperiodic,ibasis
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,lpot(MCTYPE),nloc
      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)

      dimension x(3,*)

      if(iel.eq.0) then
        i1=1
        i2=nelec
       else
        i1=iel
        i2=iel

        do 10 ic=1,ncent+nghostcent
          r_en_sav(ic)=r_en(iel,ic)
          do 10 m=1,3
            rshift_sav(m,ic)=rshift(m,iel,ic)
   10       rvec_en_sav(m,ic)=rvec_en(m,iel,ic)
        ij=0
        do 20 jj=1,nelec
          if(jj.eq.iel) goto 20
          if(jj.lt.iel) then
            i=iel
            j=jj
           else
            i=jj
            j=iel
          endif
          ij=((i-1)*(i-2))/2+j
          
          r_ee_sav(jj)=r_ee(ij)
          do 15 m=1,3
   15       rvec_ee_sav(m,jj)=rvec_ee(m,ij)
   20 continue
      endif

c Calculate e-N inter-particle distances
      do 27 i=i1,i2
        do 27 ic=1,ncent+nghostcent
          do 25 m=1,3
   25       rvec_en(m,i,ic)=x(m,i)-cent(m,ic)
          if(iperiodic.eq.0) then
            r_en(i,ic)=0
            do 26 m=1,3
   26         r_en(i,ic)=r_en(i,ic)+rvec_en(m,i,ic)**2
            r_en(i,ic)=dsqrt(r_en(i,ic))
           else
            call find_image4(rshift(1,i,ic),rvec_en(1,i,ic),r_en(i,ic))
          endif
   27 continue

c Calculate e-e inter-particle distances
      do 30 i=1,nelec
        if(iel.eq.0) i2=i-1
        do 30 j=i1,i2
          if(i.eq.j) goto 30
          if(i.lt.j) then
            ii=j
            jj=i
           else
            ii=i
            jj=j
          endif
          ij=((ii-1)*(ii-2))/2+jj
          do 28 m=1,3
   28       rvec_ee(m,ij)=x(m,ii)-x(m,jj)
          if(iperiodic.eq.0) then
            r_ee(ij)=0
            do 29 m=1,3
   29         r_ee(ij)=r_ee(ij)+rvec_ee(m,ij)**2
            r_ee(ij)=dsqrt(r_ee(ij))
           else
            call find_image3(rvec_ee(1,ij),r_ee(ij))
          endif
   30   continue

c     write(6,*) 'in distances'
c     write(6,'(''r_en(i,j)'',9f9.5)') ((r_en(i,j),i=1,nelec),j=1,2)
c     write(6,'(''r_ee(ij)'',9f9.5)') (r_ee(ij),ij=1,nelec*(nelec-1)/2)
c     write(6,'(''rvec_ee(k,ij)'',9f12.4)') ((rvec_ee(k,ij),k=1,3),ij=1,nelec*(nelec-1)/2)

      return
      end
c-----------------------------------------------------------------------
      subroutine distancese_restore(iel)
c Written by Cyrus Umrigar
c restore interparticle distances (called if move rejected)

      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use ghostatom, only: newghostype, nghostcent
      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use distances_sav, only: r_ee_sav, r_en_sav, rshift_sav, rvec_ee_sav, rvec_en_sav
      implicit real*8(a-h,o-z)




      include 'vmc.h'

      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)

c Calculate e-N inter-particle distances
      do 25 ic=1,ncent+nghostcent
        r_en(iel,ic)=r_en_sav(ic)
        do 25 m=1,3
          rshift(m,iel,ic)=rshift_sav(m,ic)
   25     rvec_en(m,iel,ic)=rvec_en_sav(m,ic)

c Calculate e-e inter-particle distances
      do 29 jj=1,nelec

        if(jj.eq.iel) goto 29
        if(jj.lt.iel) then
          i=iel
          j=jj
         else
          i=jj
          j=iel
        endif
        ij=((i-1)*(i-2))/2+j

        r_ee(ij)=r_ee_sav(jj)
        do 28 m=1,3
   28     rvec_ee(m,ij)=rvec_ee_sav(m,jj)
   29 continue

      return
      end
c-----------------------------------------------------------------------
