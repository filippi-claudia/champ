      subroutine distances(iel,x)
c Written by Cyrus Umrigar
c calculate interparticle distances
      use atom, only: cent, ncent
      use ghostatom, only: nghostcent
      use const, only: nelec
      use distances_sav, only: r_ee_sav, r_en_sav, rshift_sav, rvec_ee_sav, rvec_en_sav
      use contrl_per, only: iperiodic
      use distance_mod, only: rshift, r_en, rvec_en, r_ee, rvec_ee
      use contrl_file,    only: ounit
      use precision_kinds, only: dp
      implicit none

      integer :: i, i1, i2, ic, iel
      integer :: ii, ij, j, jj
      integer :: m

      real(dp), dimension(3, *) :: x


      if(iel.eq.0) then
        i1=1
        i2=nelec
       else
        i1=iel
        i2=iel

        do ic=1,ncent+nghostcent
          r_en_sav(ic)=r_en(iel,ic)
          do m=1,3
            rshift_sav(m,ic)=rshift(m,iel,ic)
            rvec_en_sav(m,ic)=rvec_en(m,iel,ic)
          enddo
        enddo
        ij=0
        do jj=1,nelec
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
          do m=1,3
            rvec_ee_sav(m,jj)=rvec_ee(m,ij)
          enddo
   20 continue
        enddo
      endif

c Calculate e-N inter-particle distances
      do i=i1,i2
        do ic=1,ncent+nghostcent
          do m=1,3
            rvec_en(m,i,ic)=x(m,i)-cent(m,ic)
          enddo
          if(iperiodic.eq.0) then
            r_en(i,ic)=0
            do m=1,3
              r_en(i,ic)=r_en(i,ic)+rvec_en(m,i,ic)**2
            enddo
            r_en(i,ic)=dsqrt(r_en(i,ic))
           else
            call find_image4(rshift(1,i,ic),rvec_en(1,i,ic),r_en(i,ic))
          endif
        enddo
      enddo

c Calculate e-e inter-particle distances
      do i=1,nelec
        if(iel.eq.0) i2=i-1
        do j=i1,i2
          if(i.eq.j) goto 30
          if(i.lt.j) then
            ii=j
            jj=i
           else
            ii=i
            jj=j
          endif
          ij=((ii-1)*(ii-2))/2+jj
          do m=1,3
            rvec_ee(m,ij)=x(m,ii)-x(m,jj)
          enddo
          if(iperiodic.eq.0) then
            r_ee(ij)=0
            do m=1,3
              r_ee(ij)=r_ee(ij)+rvec_ee(m,ij)**2
            enddo
            r_ee(ij)=dsqrt(r_ee(ij))
           else
            call find_image3(rvec_ee(1,ij),r_ee(ij))
          endif
   30   continue
        enddo
      enddo

c     write(ounit,*) 'in distances'
c     write(ounit,'(''r_en(i,j)'',9f9.5)') ((r_en(i,j),i=1,nelec),j=1,2)
c     write(ounit,'(''r_ee(ij)'',9f9.5)') (r_ee(ij),ij=1,nelec*(nelec-1)/2)
c     write(ounit,'(''rvec_ee(k,ij)'',9f12.4)') ((rvec_ee(k,ij),k=1,3),ij=1,nelec*(nelec-1)/2)

      return
      end
c-----------------------------------------------------------------------
      subroutine distancese_restore(iel)
c Written by Cyrus Umrigar
c restore interparticle distances (called if move rejected)

      use atom, only: ncent
      use ghostatom, only: nghostcent
      use const, only: nelec
      use distance_mod, only: rshift, r_en, rvec_en
      use distances_sav, only: r_ee_sav, r_en_sav, rshift_sav, rvec_ee_sav, rvec_en_sav
      use distance_mod, only: rshift, r_en, rvec_en, r_ee, rvec_ee
      implicit none

      integer :: i, ic, iel, ij, j
      integer :: jj, m


c Calculate e-N inter-particle distances
      do ic=1,ncent+nghostcent
        r_en(iel,ic)=r_en_sav(ic)
        do m=1,3
          rshift(m,iel,ic)=rshift_sav(m,ic)
          rvec_en(m,iel,ic)=rvec_en_sav(m,ic)
        enddo
      enddo

c Calculate e-e inter-particle distances
      do jj=1,nelec

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
        do m=1,3
          rvec_ee(m,ij)=rvec_ee_sav(m,jj)
        enddo
   29 continue
      enddo

      return
      end
c-----------------------------------------------------------------------
