c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##                                                         ##
c     #############################################################
c
c
c
c
      subroutine fit_force(istep,dt)

      use md_fit
      use md_var
      use atom, only: ncent

      implicit none
      integer j,k,m
      integer ifit,degr,ndeg
      integer ierror
      integer istep,iopt
      real*8 dt, eps, yp
      real*8 tauH,tauS
      real*8,dimension(1000):: Acoef
      real*8, allocatable :: acc_fit(:,:)
      real*8, dimension(50) :: xfit, afit, weights, Rfit
      
      allocate(acc_fit(3,ncent))

      write(6,*) "istep", istep, dt

      if (istep .le. (nfit-1)) then
      do k = 1, ncent
        do j = 1, 3
           a_save(j,k,istep)    = acc_new(j,k)
           a_save(:,:,istep+1:) = 0.d0
           acc_fit(j,k)         = acc_new(j,k)
           aY(j,k)              = 0.d0 
           x_save(j,k,istep)    = pos(j,k)
        enddo
      enddo

      else

      do k = 1, ncent
        do j = 1, 3
           if (istep .eq. nfit) then
             a_save_up(j,k,:) = a_save(j,k,:)
             x_save_up(j,k,:) = x_save(j,k,:)
           else
             do m = 1, nfit-1
               a_save_up(j,k,m) = a_save(j,k,m+1)
               x_save_up(j,k,m) = x_save(j,k,m+1)
             enddo
           endif
           a_save_up(j,k,nfit) = acc_new(j,k) 
           a_save(j,k,:)       = a_save_up(j,k,:)
           x_save(j,k,nfit)    = pos(j,k)
           x_save(j,k,:)       = x_save_up(j,k,:)
        enddo
      enddo

      degr = 3
      eps = 0.
c     deg is the degree of the polynomial 
      tauH = 1.00d0
      tauS = 1.50d0

      if(istep==1) then
       write (6,*) "tauH", tauH
       write (6,*) "tauS", tauS
      endif

      do k = 1, ncent
        do j = 1, 3
         acc_fit(j,k) = 0.0 
         do ifit = 1, nfit
           if(k .ge. 6) then
              weights(ifit) = exp(-1.0*(nfit - ifit )/tauH)
           else
              weights(ifit) = exp(-1.0*(nfit - ifit)/tauS)
           endif 
c           xfit(ifit) = x_save(j,k,ifit)
           xfit(ifit) = ifit * dt
           afit(ifit) = a_save(j,k,ifit) 
         enddo

        call dpolft((nfit), xfit, afit, weights, degr,
     &              ndeg , eps, Rfit, ierror, Acoef)

        call dp1vlu(ndeg, 0.0,xfit(nfit),acc_fit(j,k),yp,Acoef)

         aY(j,k)      = acc_new(j,k) - acc_fit(j,k) 
         acc_new(j,k) = acc_fit(j,k)

        enddo
      enddo
      write(6,*) "ACC UPDATED WITH FIT"

      endif

      return
      end
