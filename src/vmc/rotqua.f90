module rotqua_mod
contains
        subroutine gesqua(nq,xq,yq,zq,wq)
! Written by Lubos Mitas

      use precision_kinds, only: dp
      implicit none

      integer :: i20, i24, i86, iabs, ic
      integer :: icase, icasem, icosa, ii
      integer :: iocta, ip, iq, itetr
      integer :: j, jj, jmax, k
      integer :: k0 = 0, kk, kk1, kk2
      integer :: kk3, kk4, l, lp
      integer :: lq, m, mindif, ncase
      integer :: newmin, npoint, nq, ntypa
      integer :: ntype, ntypm
      real(dp) :: cf, cp, crk2, crk21
      real(dp) :: csth1, csth2, cstha, den
      real(dp) :: fi, fi0, hold, p
      real(dp) :: pi, pp, q, qq
      real(dp) :: r, rk2, rk21, rl1
      real(dp) :: rl2, s, s5, sf
      real(dp) :: si, sj, sk, sm1
      real(dp) :: sm2, snth1, snth2, sntha
      real(dp) :: sp, srk2, srk21, u
      real(dp) :: v, w, wq, ww
      real(dp) :: x, xq, y, yq
      real(dp) :: z, zq

        parameter (ncase=10,npoint=86,ntype=4,ntypm=6)

        dimension wq(*),xq(*),yq(*),zq(*), &
              iocta(ntype),icosa(ntype),ww(ncase,ntypm),iq(ncase), &
        itetr(ntype),i20(ntype),i24(ntype),i86(ntypm)

        data iq/4,6,12,18,26,32,50,24,20,86/
        data itetr/4,0,0,0/
        data iocta/6,12,8,24/
        data icosa/2,10,20,0/
        data i24/24,0,0,0/
        data i20/20,0,0,0/
        data i86/6,0,8,24,24,24/
        data pi/3.14159265359d0/

        icasem=1
        mindif=iabs(nq-iq(1))
        do ic=2,ncase
          newmin=iabs(nq-iq(ic))
          if (mindif.gt.newmin) then
            mindif=newmin
            icasem=ic
          endif
        end do
        icase=icasem
        nq=iq(icase)

      do ic=1,ncase
       do k=1,ntypm
         ww(ic,k)=0.d0
        end do
      end do
        do j=1,npoint
          xq(j)=0.d0
          yq(j)=0.d0
          zq(j)=0.d0
          wq(j)=0.d0
        end do

        ww(1,1)=1.d0/4.d0
        ww(2,1)=1.d0/6.d0
        ww(3,1)=1.d0/12.d0
        ww(4,1)=1.d0/30.d0
        ww(5,1)=1.d0/21.d0
        ww(6,1)=5.d0/168.d0
        ww(7,1)=4.d0/315.d0
        ww(3,2)=1.d0/12.d0
        ww(4,2)=1.d0/15.d0
        ww(5,2)=4.d0/105.d0
        ww(6,2)=5.d0/168.d0
        ww(7,2)=64.d0/2835.d0
        ww(5,3)=27.d0/840.d0
        ww(6,3)=27.d0/840.d0
        ww(7,3)=27.d0/1280.d0
        ww(7,4)=14641.d0/725760.d0
        ww(8,1)=1.d0/24.d0
        ww(9,1)=1.d0/20.d0
        ww(10,1)=.0115440115441d0
        ww(10,3)=.0119439090859d0
        ww(10,4)=.0111105557106d0
        ww(10,5)=.0118765012945d0
        ww(10,6)=.0118123037469d0

         jj=0
        ntypa=ntype
        if(nq.eq.86) ntypa=ntypm
        do k=1,ntypa
          jmax=iocta(k)
          if (nq.eq.4) jmax=itetr(k)
          if (nq.eq.12 .or. nq.eq.32) jmax=icosa(k)
          if (nq.eq.20)        jmax=i20(k)
          if (nq.eq.24)        jmax=i24(k)
          if (nq.eq.86)        jmax=i86(k)
          do j=1,jmax
            jj=jj+1
            if (jj.le.nq)        wq(jj)=ww(icase,k)
          end do
        end do
        p=1.d0/dsqrt(2.d0)
        q=1.d0/dsqrt(3.d0)
        r=1.d0/dsqrt(11.d0)
        s=3.d0/dsqrt(11.d0)

        select case(icase)
          case(1)
            goto 10
          case(2)
            goto 20
          case(3)
            goto 130
          case(4)
            goto 20
          case(5)
            goto 20
          case(6)
            goto 130
          case(7)
            goto 20
          case(8)
            goto 230
          case(9)
            goto 130
          case(10)
            goto 530
          case default
        end select

!       tetrahedron symmetry quadrature

      10     continue
        do j=1,4
          xq(j)=q
          if (j.gt.2) xq(j)=-q
          yq(j)=q
          if (j.eq.2 .or. j.eq.4) yq(j)=-q
          zq(j)=q
          if (j.eq.2 .or. j.eq.3) zq(j)=-q
        end do
        return

!       octahedron symmetry quadrature

      20      continue
        ip=0
        lp=0
        si=-1.d0
        do ii=1,2
            si=-si
            xq(ii)=si
            yq(2+ii)=si
            zq(4+ii)=si
            if (nq.gt.6) then
              sj=-1.d0
              do jj=1,2
                sj=-sj
                ip=ip+1
                xq(6+ip)=si*p
                yq(6+ip)=sj*p
                yq(10+ip)=si*p
                zq(10+ip)=sj*p
                zq(14+ip)=si*p
                xq(14+ip)=sj*p
                if (nq.gt.18) then
                  sk=-1.d0
                  do kk=1,2
                    sk=-sk
                    lp=lp+1
                    xq(18+lp)=si*q
                    yq(18+lp)=sj*q
                    zq(18+lp)=sk*q
                    if (nq.gt.26) then
                      xq(26+lp)=si*r
                      yq(26+lp)=sj*r
                      zq(26+lp)=sk*s
                      xq(34+lp)=si*r
                      yq(34+lp)=sj*s
                      zq(34+lp)=sk*r
                      xq(42+lp)=si*s
                      yq(42+lp)=sj*r
                      zq(42+lp)=sk*r
                    endif
                  end do
                endif
              end do
            endif
        end do
        return

!       icosahedron symmetry quadrature

      130       zq(1)=1.d0
        zq(2)=-1.d0
        fi0=pi/5.d0
        s5=dsqrt(5.d0)
        cstha=1.d0/s5
        sntha=2.d0/s5
        den=dsqrt(15.d0+6.d0*s5)
        csth1=(2.d0+s5)/den
        snth1=dsqrt(1.d0-csth1**2)
        csth2=1.d0/den
        snth2=dsqrt(1.d0-csth2**2)

        do j=1,5
          rk2=2*dble(j-1)
          crk2=dcos(rk2*fi0)
          srk2=dsin(rk2*fi0)
          if (nq.ne.20)        then
            xq(j+2)=sntha*crk2
            yq(j+2)=sntha*srk2
            zq(j+2)=cstha
          endif
          rk21=rk2+1.d0
          crk21=dcos(rk21*fi0)
          srk21=dsin(rk21*fi0)
          if (nq.ne.20)        then
            xq(j+7)=sntha*crk21
            yq(j+7)=sntha*srk21
            zq(j+7)=-cstha
          endif
          if (nq.gt.12)        then
            if (nq.eq.20) k0=0
            if (nq.eq.32) k0=12
            kk1=j+k0
            kk2=j+k0+5
            kk3=j+k0+10
            kk4=j+k0+15
            xq(kk1)=snth1*crk21
            yq(kk1)=snth1*srk21
            zq(kk1)=csth1
            xq(kk2)=snth2*crk21
            yq(kk2)=snth2*srk21
            zq(kk2)=csth2
            xq(kk3)=snth1*crk2
            yq(kk3)=snth1*srk2
            zq(kk3)=-csth1
            xq(kk4)=snth2*crk2
            yq(kk4)=snth2*srk2
            zq(kk4)=-csth2
          endif
        end do
        return

      230       continue
        lq=0
        do l=1,3
          if (l.eq.1) then
            x=.8662468181078d0
            y=.4225186537611d0
              z=.2666354015167d0
          endif
          if (l.eq.2) then
            hold=y
            y=x
            x=z
            z=hold
          endif
          if (l.eq.3) then
            hold=z
            z=y
            y=x
            x=hold
          endif
          cp=dcos(-pi/2.d0)
          sp=dsin(-pi/2.d0)
          u=cp*x-sp*z
          v=y
          w=sp*x+cp*z
          do m=1,4
            fi=pi*dble(m-1)/2
            cf=dcos(fi)
            sf=dsin(fi)
            lq=lq+1
            xq(lq)=x*cf-y*sf
            yq(lq)=x*sf+y*cf
            zq(lq)=z
            lq=lq+1
            xq(lq)=u*cf-v*sf
            yq(lq)=u*sf+v*cf
            zq(lq)=w
          end do
        end do
        return

      530       continue
        ip=0
        lp=0
        pp=.927330657151d0
        qq=.374243039090d0
        rl1=.369602846454d0
        sm1=.852518311701d0
        rl2=.694354006603d0
        sm2=.189063552885d0
        si=-1.d0
        do ii=1,2
          si=-si
          xq(ii)=si
          yq(2+ii)=si
          zq(4+ii)=si
          sj=-1.d0
          do jj=1,2
            sj=-sj
            ip=ip+1
            xq(62+ip)=si*pp
            yq(62+ip)=sj*qq
            xq(66+ip)=si*pp
            zq(66+ip)=sj*qq
            yq(70+ip)=si*pp
            zq(70+ip)=sj*qq
            xq(74+ip)=si*qq
            yq(74+ip)=sj*pp
            xq(78+ip)=si*qq
            zq(78+ip)=sj*pp
            yq(82+ip)=si*qq
            zq(82+ip)=sj*pp
            sk=-1.d0
            do kk=1,2
              sk=-sk
              lp=lp+1
              xq(6+lp)=si*q
              yq(6+lp)=sj*q
              zq(6+lp)=sk*q
              xq(14+lp)=si*rl1
              yq(14+lp)=sj*rl1
              zq(14+lp)=sk*sm1
              xq(22+lp)=si*rl1
              yq(22+lp)=sj*sm1
              zq(22+lp)=sk*rl1
              xq(30+lp)=si*sm1
              yq(30+lp)=sj*rl1
              zq(30+lp)=sk*rl1
              xq(38+lp)=si*rl2
              yq(38+lp)=sj*rl2
              zq(38+lp)=sk*sm2
              xq(46+lp)=si*rl2
              yq(46+lp)=sj*sm2
              zq(46+lp)=sk*rl2
              xq(54+lp)=si*sm2
              yq(54+lp)=sj*rl2
              zq(54+lp)=sk*rl2
            end do
          end do
        end do 
        return
        end
!-----------------------------------------------------------------------

        subroutine rotqua
! Written by Lubos Mitas
      use qua, only: nquad, wq, xq, xq0, yq, yq0, zq, zq0
      use precision_kinds, only: dp
      use random_mod, only: random_dp
      implicit none

      integer :: iq
      real(dp) :: cfi, sfi
      real(dp) :: sthet, theta, u1, u2
      real(dp) :: usum, uu, x1, x2
      real(dp) :: x3, xsum, xsum2, y1
      real(dp) :: y2, y3, yy1, yy2
      real(dp) :: yy3, z1, z2, z3
      real(dp) :: zz1, zz2, zz3


      2     x1=1.d0-2.d0*random_dp()
        x2=1.d0-2.d0*random_dp()
        xsum=x1*x1+x2*x2
        if (xsum.ge.1.d0) goto 2
        xsum2=2.d0*dsqrt(dabs(1.d0-xsum))
        x1=x1*xsum2
        x2=x2*xsum2
        x3=1.d0-2.d0*xsum

        theta=dacos(x3)
        sthet=dsin(theta)
        if (sthet.lt.1.d-05) then
          cfi=1.d0
          sfi=0.d0
          x1=0.d0
          x2=0.d0
          x3=1.d0
          sthet=0.d0
         else
          cfi=x1/sthet
          sfi=x2/sthet
        endif

        yy1=x3*cfi
        yy2=x3*sfi
        yy3=-sthet
        zz1=-sfi
        zz2=cfi
        zz3=0.d0

      3   u1=random_dp()*2.d0-1.d0
        u2=random_dp()*2.d0-1.d0
        usum=u1*u1+u2*u2
        if(usum.ge.1.d0) goto 3
        uu=dsqrt(usum)
        u1=u1/uu
        u2=u2/uu
!       yu1=yy1*u1
!       yu2=yy2*u1
!       yu3=yy3*u1
!       zu1=zz1*u2
!       zu2=zz2*u2
!       zu3=zz3*u2
        y1=yy1*u1+zz1*u2
        y2=yy2*u1+zz2*u2
        y3=yy3*u1+zz3*u2
        z1=yy1*u2-zz1*u1
        z2=yy2*u2-zz2*u1
        z3=yy3*u2-zz3*u1
        do iq=1,nquad
          xq(iq)=xq0(iq)*x1+yq0(iq)*y1+zq0(iq)*z1
          yq(iq)=xq0(iq)*x2+yq0(iq)*y2+zq0(iq)*z2
          zq(iq)=xq0(iq)*x3+yq0(iq)*y3+zq0(iq)*z3
        end do
        return
        end
end module 
