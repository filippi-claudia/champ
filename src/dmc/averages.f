      module averages
      contains
      subroutine deriv(wtg,eold,pwt,ajac,psid,psij,idrifdifgfunc,iw,mwalk)

      use c_averages, only: prop,wprop
      use c_averages_index, only: jderiv
      use multiple_geo, only: nforce
      use precision_kinds, only: dp

      implicit none

      integer :: idrifdifgfunc, ifr, iw, mwalk
      real(dp) :: ro, wtg
      real(dp), dimension(mwalk, *) :: eold
      real(dp), dimension(mwalk, *) :: pwt
      real(dp), dimension(mwalk, *) :: ajac
      real(dp), dimension(mwalk, *) :: psij
      real(dp), dimension(mwalk, *) :: psid


      do ifr=1,nforce
        if(idrifdifgfunc.gt.0) then
          ro=0.d0
         else
          ro=2*(log(abs(psid(iw,ifr)))+psij(iw,ifr))
          if(ifr.gt.1) ro=ro+log(ajac(iw,ifr))
        endif
        prop(jderiv(1,ifr))=eold(iw,ifr)
        prop(jderiv(2,ifr))=eold(iw,1)*(pwt(iw,ifr)+ro)
        prop(jderiv(3,ifr))=pwt(iw,ifr)+ro
        wprop(jderiv(1,ifr))=wtg
        wprop(jderiv(2,ifr))=wtg
        wprop(jderiv(3,ifr))=wtg
      enddo
      return
      end

      subroutine init_averages_index

      use c_averages_index, only: jderiv,jeloc
      use multiple_geo, only: nforce
      use prp000,  only: nprop

      implicit none

      integer :: ifr, j, k


      nprop=0
c elocal
      j=nprop+1
      jeloc=j
      nprop=j
c deriv
      j=nprop
      do ifr=1,nforce
       do k=1,3
        j=j+1
        jderiv(k,ifr)=j
       enddo
      enddo
      nprop=j
      return
      end

      subroutine average(ido)

      use c_averages, only: cum_av,cum_av2,cum_w,mprop,prop,wprop
      use precision_kinds, only: dp
      use prp000,  only: nprop

      implicit none

      integer :: i, ido

      real(dp), dimension(mprop) :: sum_av
      real(dp), dimension(mprop) :: sum_w

      if(ido.eq.0)then
       do i=1,nprop
        cum_av(i)=0
        cum_av2(i)=0
        cum_w(i)=0
        sum_av(i)=0
        sum_w(i)=0
       enddo
      elseif(ido.eq.1)then
       do i=1,nprop
        sum_av(i)=sum_av(i)+prop(i)*wprop(i)
        sum_w(i)=sum_w(i)+wprop(i)
       enddo
      else
       do i=1,nprop
        cum_av(i)=cum_av(i)+sum_av(i)
        cum_av2(i)=cum_av2(i)+sum_av(i)**2*sum_w(i)
        cum_w(i)=cum_w(i)+sum_w(i)
        sum_av(i)=0
        sum_w(i)=0
       enddo
      endif
      return
      end

      subroutine average_write

      use c_averages, only: cum_av,cum_w
      use c_averages_index, only: jderiv
      use multiple_geo, only: nforce
      use precision_kinds, only: dp

      implicit none

      integer :: ifr
      real(dp) :: derivtotave, egave

      egave=cum_av(jderiv(1,1))/cum_w(jderiv(1,1))
      do ifr=2,nforce
       derivtotave=-(cum_av(jderiv(1,ifr))-cum_av(jderiv(1,1))
     &              +cum_av(jderiv(2,ifr))-cum_av(jderiv(2,1))
     &              -egave*(cum_av(jderiv(3,ifr))-cum_av(jderiv(3,1))))
       derivtotave=derivtotave/cum_w(jderiv(1,1))
c      write(ounit,*)'test deriv ',ifr,derivtotave
      enddo
      return
      end
      end module
