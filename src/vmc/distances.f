      module distances_mod
      contains
      subroutine distances(iel,x)
c Written by Cyrus Umrigar
c calculate interparticle distances
      use contrl_file,    only: ounit
      use contrl_per, only: iperiodic
      use distance_mod, only: r_en, rvec_en, r_ee, rvec_ee
      use distances_sav, only: r_ee_sav, r_en_sav, rvec_ee_sav, rvec_en_sav
      use precision_kinds, only: dp
      use find_pimage, only: find_image3, find_image_pbc
      use system, only: cent, ncent, nghostcent, nelec
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
               rvec_en_sav(m,ic)=rvec_en(m,iel,ic)
            enddo
         enddo
         
         ij=0
         
         do jj=1,iel-1

            ij=((iel-1)*(iel-2))/2+jj
            
            r_ee_sav(jj)=r_ee(ij)
            do m=1,3
               rvec_ee_sav(m,jj)=rvec_ee(m,ij)
            enddo
            
         enddo

         do jj=iel+1,nelec

            ij=((jj-1)*(jj-2))/2+iel

            r_ee_sav(jj)=r_ee(ij)
            do m=1,3
               rvec_ee_sav(m,jj)=rvec_ee(m,ij)
            enddo

         enddo
         
      endif



      if(iperiodic.eq.0) then


c     Calculate e-N inter-particle distances
         do i=i1,i2
            do ic=1,ncent+nghostcent
            
               do m=1,3
                  rvec_en(m,i,ic)=x(m,i)-cent(m,ic)
               enddo
         
               r_en(i,ic)=0
               do m=1,3
                  r_en(i,ic)=r_en(i,ic)+rvec_en(m,i,ic)**2
               enddo
               r_en(i,ic)=dsqrt(r_en(i,ic))

            enddo
         enddo




c     Calculate e-e inter-particle distances      
         if(iel.eq.0) then

            do i=2,nelec
               do j=1,i-1
                  
                  ij=((i-1)*(i-2))/2+j
               
                  do m=1,3
                     rvec_ee(m,ij)=x(m,i)-x(m,j)
                  enddo
                  
                  r_ee(ij)=0
                  do m=1,3
                     r_ee(ij)=r_ee(ij)+rvec_ee(m,ij)**2
                  enddo
                  r_ee(ij)=dsqrt(r_ee(ij))

               enddo
            enddo
         
            
         else

c     iel!=0         
         
            do i=1,iel-1
            
               ij=((iel-1)*(iel-2))/2+i
               
               do m=1,3
                  rvec_ee(m,ij)=x(m,iel)-x(m,i)
               enddo
               
               r_ee(ij)=0
               do m=1,3
                  r_ee(ij)=r_ee(ij)+rvec_ee(m,ij)**2
               enddo
               r_ee(ij)=dsqrt(r_ee(ij))
            
         
            enddo
         

            do i=iel+1,nelec
                        
               ij=((i-1)*(i-2))/2+iel
            
               do m=1,3
                  rvec_ee(m,ij)=x(m,i)-x(m,iel)
               enddo
               
               
               r_ee(ij)=0
               do m=1,3
                  r_ee(ij)=r_ee(ij)+rvec_ee(m,ij)**2
               enddo
               r_ee(ij)=dsqrt(r_ee(ij))
               
         enddo
         
         

         
      endif



                  
         
      else


!     periodic systems

         
c     Calculate e-N inter-particle distances
         do i=i1,i2
            do ic=1,ncent+nghostcent
               
               do m=1,3
                  rvec_en(m,i,ic)=x(m,i)-cent(m,ic)
               enddo
               
               call find_image_pbc(rvec_en(1,i,ic),r_en(i,ic))
               
            enddo
         enddo
         



c     Calculate e-e inter-particle distances      
         if(iel.eq.0) then

            do i=2,nelec
               do j=1,i-1

                  ij=((i-1)*(i-2))/2+j
               
                  do m=1,3
                     rvec_ee(m,ij)=x(m,i)-x(m,j)
                  enddo
                  
                  call find_image_pbc(rvec_ee(1,ij),r_ee(ij))
                 
               enddo
            enddo
         
            
         else

c     iel!=0         
            
            do i=1,iel-1
               
               ij=((iel-1)*(iel-2))/2+i
               
               do m=1,3
                  rvec_ee(m,ij)=x(m,iel)-x(m,i)
               enddo
               
               call find_image_pbc(rvec_ee(1,ij),r_ee(ij))
               
            enddo
         

            do i=iel+1,nelec
                        
               ij=((i-1)*(i-2))/2+iel
            
               do m=1,3
                  rvec_ee(m,ij)=x(m,i)-x(m,iel)
               enddo
            
               call find_image_pbc(rvec_ee(1,ij),r_ee(ij))
                        
            enddo

         

         
         endif
         !!iel enif


         
         
      endif
      !! periodic endif




      

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

      use system, only: ncent, nghostcent, nelec
      use distances_sav, only: r_ee_sav, r_en_sav, rvec_ee_sav, rvec_en_sav
      use distance_mod, only: r_en, rvec_en, r_ee, rvec_ee
      implicit none

      integer :: i, ic, iel, ij, j
      integer :: jj, m


c Calculate e-N inter-particle distances
      do ic=1,ncent+nghostcent
        r_en(iel,ic)=r_en_sav(ic)
        do m=1,3
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
      end module
