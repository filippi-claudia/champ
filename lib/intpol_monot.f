      function intpol_monot(x,y,xi,delta_c3_max) ! author: Cyrus Umrigar
c **Warning: This routine is prone to numerical instabilities
c **Warning: This routine works for all cases I have tried it on but it
c is not guanranteed to work.  Since the Newton root finding method has
c good local convergence but poor global convergence, I have augmented
c it in a non-standard way.

c Interpolates a monotonic function through 3 points using the form
c y=c1+c2*exp(c3*x)
c x and y are arrays of length 3 with the original data.
c xi is the point at which the interpolated value is needed.

      implicit real*8(a-h,o-z)
      real*8 intpol_monot
      dimension x(*),y(*)

c Make initial guess for exponent
      c3=(2/(x(1)-x(3)))*
     &log((y(1)-y(2))*(x(2)-x(3))/((y(2)-y(3))*(x(1)-x(2))))

c Iterate with Newton to determine exponent
      iter_newton=0
      delta_c3_mx=delta_c3_max
   10 iter_newton=iter_newton+1
      fun=(exp(c3*x(1))-exp(c3*x(3)))/(exp(c3*x(2))-exp(c3*x(3)))
     &-((y(1)-y(3))/(y(2)-y(3)))
      dfun_dc3=((exp(c3*x(2))-exp(c3*x(3)))*
     &(x(1)*exp(c3*x(1))-x(3)*exp(c3*x(3)))-
     &(exp(c3*x(1))-exp(c3*x(3)))*
     &(x(2)*exp(c3*x(2))-x(3)*exp(c3*x(3))))/
     &(exp(c3*x(2))-exp(c3*x(3)))**2
      delta_c3=-fun/dfun_dc3
      delta_c3_newton=delta_c3
      if(iter_newton.ge.2 .and. delta_c3/delta_c3_old.gt.0.25d0) then
c       delta_c3=delta_c3*delta_c3_old/(delta_c3_old-delta_c3)
        delta_c3=delta_c3*(1+delta_c3/delta_c3_old)
      endif
      if(abs(delta_c3).gt.delta_c3_mx) then
        delta_c3=max(-delta_c3_mx,min(delta_c3_mx,-fun/dfun_dc3))
        delta_c3_mx=delta_c3_mx/2
      endif
      c3=c3+delta_c3

      if(abs(delta_c3).gt.1.d-6) then
        if(iter_newton.lt.50) then
          delta_c3_old=delta_c3_newton
          goto 10
         else
          stop 'intpol_monot did not converge in 50 steps'
        endif
      endif

c Determine other 2 coeffs.
      c2=(y(2)-y(1))/(exp(c3*x(2))-exp(c3*x(1)))
      c1=y(1)-c2*exp(c3*x(1))

c Determine interpolated value
      intpol_monot=c1+c2*exp(c3*xi)

      return
      end
