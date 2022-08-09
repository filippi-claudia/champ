      module trexio_basis_fns_mod

      contains
      subroutine trexio_basis_fns(ie1,ie2,rvec_en,r_en,ider)
c calculate the values of the basis functions and their derivatives
c ider = 0 -> value
c ider = 1 -> value, gradient
c ider = 2 -> value, gradient, laplacian
c ider = 3 -> value, gradient, laplacian, forces

      use numbas_mod, only: MRWF
      use system, only: iwctype, ncent, ncent_tot, nghostcent, nelec
      use numbas, only: iwrwf, nrbas!, rmax
      use numbas1, only: iwlbas, nbastyp
      use phifun, only: phin, dphin, d2phin, d2phin_all, d3phin, n0_nbasis
      use multiple_geo, only: iwf
      use m_force_analytic, only: iforce_analy
      use splfit_mod, only: splfit
      use slm_mod, only: slm
#if defined(TREXIO_FOUND)
      use m_trexio_basis,     only: slm_per_l, index_slm, num_rad_per_cent
      use m_trexio_basis,     only: num_ao_per_cent, ao_radial_index
      use m_trexio_basis,     only: basis_num_shell, basis_shell_ang_mom
#endif

      use precision_kinds, only: dp
      implicit none

      integer :: it, ic, ider, irb
      integer :: iwlbas0
      integer :: j, k, nrbasit, nbastypit, i
      integer :: ie1, ie2, l, ilm, num_slms, maxlval


      real(dp) :: r, r2, ri, ri2
      real(dp), dimension(3, nelec, ncent_tot) :: rvec_en
      real(dp), dimension(nelec, ncent_tot) :: r_en
      real(dp), dimension(4, MRWF) :: wfv
      real(dp), dimension(3) :: xc
      real(dp), parameter :: one = 1.d0

      ! Temporary arrays for basis function values and derivatives
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: ddy_lap(:)
      real(dp), allocatable :: dy(:,:)
      real(dp), allocatable :: ddy(:,:,:)
      real(dp), allocatable :: dlapy(:,:)

      integer                       :: upper_range, lower_range
      integer                       :: upper_rad_range, lower_rad_range

#ifndef VECTORIZATION
      do j=ie1,ie2
        n0_nbasis(j)=0
      enddo
#endif

#if defined(TREXIO_FOUND)

      lower_range = 1
      lower_rad_range = 1


c loop through centers
      do ic=1,ncent+nghostcent
        upper_range = lower_range + num_ao_per_cent(ic) -1
        upper_rad_range = lower_rad_range + num_rad_per_cent(ic) -1

        it=iwctype(ic)
        nrbasit   = num_rad_per_cent(ic)
        nbastypit = num_ao_per_cent(ic)

        ! Get the maximum angular momentum for this center
        maxlval = maxval(basis_shell_ang_mom(lower_rad_range:upper_rad_range))
        ! num_slms will give number of slms needed to evaluate per atom
        num_slms = sum(slm_per_l(1:maxlval+1))

        allocate (y(num_slms))
        allocate (ddy_lap(num_slms))
        allocate (dy(3,num_slms))
        allocate (ddy(3,3,num_slms))
        allocate (dlapy(3,num_slms))

c     numerical atomic orbitals
        do k=ie1,ie2

c get distance to center

          xc(1)=rvec_en(1,k,ic)
          xc(2)=rvec_en(2,k,ic)
          xc(3)=rvec_en(3,k,ic)

          r=r_en(k,ic)
          r2=r*r
          ri=one/r
          ri2=ri*ri

          do irb=1,nrbasit
          ! only evaluate for r <= rmax
          ! if (r <= rmax(irb,it)) then
            call splfit(r,irb,it,iwf,wfv(1,irb),ider)
          ! else
          !   wfv(1:4,irb)=0.d0
          ! endif
          enddo

          ! Get the Slm evaluated and store them arrays
          do i=1, num_slms
            call slm(i,xc,r2,y(i),dy(1,i),ddy(1,1,i),ddy_lap(i),dlapy(1,i),ider)
          enddo


          l = 1
          iwlbas0=0
          ! Run a loop over all the AOs in this center
          do ilm=lower_range, upper_range
            iwlbas0=index_slm(ilm)

!     compute sml and combine to generate molecular orbitals
            irb = ao_radial_index(ilm)
            call trexio_phi_combine(iwlbas0,xc,ri,ri2,wfv(1,irb),
     &            y(iwlbas0),
     &            dy(:,iwlbas0),
     &            ddy(:,:,iwlbas0),
     &            ddy_lap(iwlbas0),
     &            dlapy(:,iwlbas0),
     &            phin(ilm,k),
     &            dphin(ilm,k,:),
     &            d2phin(ilm,k),
     &            d2phin_all(1,1,ilm,k),
     &            d3phin(1,ilm,k),
     &            ider)

#ifndef VECTORIZATION
            ! localization
            call n0_inc(l,k,ic)
#endif
            l = l + 1
          enddo
        enddo ! loop over electrons
        lower_range = upper_range + 1
        lower_rad_range = upper_rad_range + 1

!     deallocate temporary arrays
        deallocate (y)
        deallocate (ddy_lap)
        deallocate (dy)
        deallocate (ddy)
        deallocate (dlapy)

      enddo ! loop over all atoms
#endif
      return
      end
c-------------------------------------------------------------------
      subroutine trexio_phi_combine(l,xc,ri,ri2,wfv,y,dy,ddy,ddy_lap,dlapy,phi,dphi,d2phi,d2phi_all,d3phi,ider)
      use precision_kinds, only: dp
      implicit none

      integer :: iforce_analy, ider, ii, jj, l
      real(dp) :: d2phi, ddy_lap, dum, dum1, phi
      real(dp) :: prod, ri, ri2, ri3
      real(dp) :: y
      real(dp), dimension(3) :: xc
      real(dp), dimension(3) :: xcri
      real(dp), dimension(4) :: wfv
      real(dp), dimension(3) :: dy
      real(dp), dimension(3, 3) :: ddy
      real(dp), dimension(3) :: dphi
      real(dp), dimension(3, 3) :: d2phi_all
      real(dp), dimension(3) :: d3phi
      real(dp), dimension(3) :: dlapy
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: three = 3.d0
      real(dp), parameter :: four = 4.d0
      real(dp), parameter :: five = 5.d0
      real(dp), parameter :: six = 6.d0

c     phi is computed for all ider values
      phi=y*wfv(1)

      if(ider.eq.1) then

        xcri(1)=xc(1)*ri
        xcri(2)=xc(2)*ri
        xcri(3)=xc(3)*ri

        do jj=1,3
          dphi(jj)=y*xcri(jj)*wfv(2)+dy(jj)*wfv(1)
        enddo

      elseif(ider.ge.2) then


        xcri(1)=xc(1)*ri
        xcri(2)=xc(2)*ri
        xcri(3)=xc(3)*ri

        d2phi=y*wfv(3)+y*two*ri*wfv(2)+ddy_lap*wfv(1)
        dum=0.d0
        do jj=1,3
          dphi(jj)=y*xcri(jj)*wfv(2)+dy(jj)*wfv(1)
          dum=dum+dy(jj)*xcri(jj)
        enddo
        d2phi=d2phi+two*dum*wfv(2)


        if(ider.eq.3) then

          do jj=1,3
            dum1=0
            do ii=1,3
              dum1=dum1+ddy(jj,ii)*xcri(ii)
            enddo
            d3phi(jj)=wfv(4)*y*xcri(jj)
     &               +wfv(3)*(dy(jj)+two*xcri(jj)*(y*ri+dum))
     &               +wfv(2)*(xcri(jj)*(ddy_lap-two*ri*(dum+y*ri))+two*(dum1+two*dy(jj)*ri))
     &               +wfv(1)*dlapy(jj)
          enddo

          do jj=1,3
            do ii=jj,3
              prod=xcri(jj)*xcri(ii)
              d2phi_all(ii,jj)=ddy(ii,jj)*wfv(1)+wfv(2)*(dy(ii)*xcri(jj)+dy(jj)*xcri(ii)-y*ri*prod)+wfv(3)*y*prod
              d2phi_all(jj,ii)=d2phi_all(ii,jj)
            enddo
            d2phi_all(jj,jj)=d2phi_all(jj,jj)+y*ri*wfv(2)
          enddo

        endif

      endif

      return
      end
c-------------------------------------------------------------------
      subroutine n0_inc(l,k,ic)

      use phifun, only: phin, dphin, n0_ibasis, n0_ic, n0_nbasis
      implicit none

      integer :: ic, k, l


      if(abs(phin(l,k))+abs(dphin(l,k,1))+abs(dphin(l,k,2))+abs(dphin(l,k,3)).gt.1.d-20)then

       n0_nbasis(k)=n0_nbasis(k)+1
       n0_ibasis(n0_nbasis(k),k)=l
       n0_ic(n0_nbasis(k),k)=ic
      endif

      return
      end
      end module
