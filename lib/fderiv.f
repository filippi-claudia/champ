      subroutine fderiv(rho,gradrho,del2rho,fder,eps,h,n,ipr)
      implicit real*8(a-h,o-z)
c Calculates the functional derivative for a spherically symmetric rho
c using 5 pt formulae on uniform radial mesh.
c The integrand of the functional may depend on rho, gradrho and del2rho
c The 5 pt formulae can be found in Abramowitz and Stegun pg 914
c                                             Cyrus Feb 1992

c f       = The functional is the integral of this function. Fderiv
c         = calculates the functional derivative of the integral.
c rho     = the charge density
c gradrho = grad(rho)
c del2rho = laplacian(rho)
c fder    = functional derivative
c eps     = initial radial coordinate
c h       = step size in radial coordinate
c n       = number of points
c ipr     = print switch
c The radial mesh is uniform i.e.  r(i) = eps+(i-1)*h

c Presently uses auxilliary arrays of dim MRAD.
c We could use one work array to calculate various c quantities and
c keep adding them to fder.
c I don't do that, because I want to be able to print out the individual
c quantities at the end of the program.

      parameter(
     1 c00=  1.d0, c01= -8.d0, c02=  0.d0, c03= 8.d0, c04=-1.d0,
     1 c10= -3.d0, c11=-10.d0, c12= 18.d0, c13=-6.d0, c14= 1.d0,
     1 c20=-25.d0, c21= 48.d0, c22=-36.d0, c23=16.d0, c24=-3.d0)
      parameter (one=1.d0,small=1.d-3)
      parameter (MRAD=2001)
      dimension rho(*),gradrho(*),del2rho(*),fder(*)
     1,df_dgradrho(MRAD) ,div_df_dgradrho(MRAD)
     1,df_ddel2rho(MRAD) ,del2_df_ddel2rho(MRAD)
      external f

      if(n.gt.MRAD) stop 'n must be smaller than MRAD in fderiv'

      hinv12=1.d0/(h*12.d0)

      if(ipr.ge.1) write(6,*) '  r      rho        gradrho     del2rho
     1  df_drho   df_dgradrho df_ddel2rho -div_df_dgradrho del2_df_ddel2
     1rho fder'

      do 10 i=1,n
        eps1=small*abs(gradrho(i))
c       eps1=small*max(1.d-3,abs(gradrho(i)))
        eps2=eps1*2
        eps1inv12=1.d0/(12.d0*eps1)
        df_dgradrho(i)=eps1inv12*
     1  (c00*f(rho(i),abs(gradrho(i)-eps2),del2rho(i))
     1  +c01*f(rho(i),abs(gradrho(i)-eps1),del2rho(i))
     1  +c02*f(rho(i),abs(gradrho(i)     ),del2rho(i))
     1  +c03*f(rho(i),abs(gradrho(i)+eps1),del2rho(i))
     1  +c04*f(rho(i),abs(gradrho(i)+eps2),del2rho(i)))

        eps1=small*abs(del2rho(i))
c       eps1=small*max(1.d-3,abs(del2rho(i)))
        eps2=eps1*2
        eps1inv12=1.d0/(12.d0*eps1)
        df_ddel2rho(i)=eps1inv12*
     1  (c00*f(rho(i),abs(gradrho(i)),del2rho(i)-eps2)
     1  +c01*f(rho(i),abs(gradrho(i)),del2rho(i)-eps1)
     1  +c02*f(rho(i),abs(gradrho(i)),del2rho(i)     )
     1  +c03*f(rho(i),abs(gradrho(i)),del2rho(i)+eps1)
     1  +c04*f(rho(i),abs(gradrho(i)),del2rho(i)+eps2))
   10   continue

      call divergence(df_dgradrho,div_df_dgradrho,eps,h,n)
      call laplac(df_ddel2rho,del2_df_ddel2rho,eps,h,n)


      do 20 i=1,n
        r=eps+(i-1)*h

        eps1=small*rho(i)
        eps2=eps1*2
        eps1inv12=1.d0/(12.d0*eps1)
        df_drho=eps1inv12*
     1  (c00*f(rho(i)-eps2,abs(gradrho(i)),del2rho(i))
     1  +c01*f(rho(i)-eps1,abs(gradrho(i)),del2rho(i))
     1  +c02*f(rho(i)     ,abs(gradrho(i)),del2rho(i))
     1  +c03*f(rho(i)+eps1,abs(gradrho(i)),del2rho(i))
     1  +c04*f(rho(i)+eps2,abs(gradrho(i)),del2rho(i)))

        fder(i)=df_drho-div_df_dgradrho(i)+del2_df_ddel2rho(i)
        if(ipr.ge.1) write(6,'(f7.3,9d12.5)')
     1  r,rho(i),gradrho(i),del2rho(i),
     1  df_drho,df_dgradrho(i),df_ddel2rho(i),
     1  -div_df_dgradrho(i),del2_df_ddel2rho(i),fder(i)
   20   continue
      return
      end
