      subroutine orbitals_pw(x,orb,dorb,ddorb)
c Written by Cyrus Umrigar
c Calculate pw orbitals, gradient and laplacian.
c isortg could be used to map g-vectors from iv to ig and
c isortk could be used to map k-vectors.
c At present it is assumed that both g- and k-vectors are in the correct order.

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use ewald_mod, only: IVOL_RATIO
      use ewald_mod, only: NGVECX
      use const, only: nelec, ipr
      use periodic, only: glatt
      use periodic, only: glatt_sim, gnorm, gvec, igmult, igvec
      use periodic, only: ireal_imag, k_inv, kvec, nband, ng1d, ng1d_sim, ngnorm_orb
      use periodic, only: ngvec_orb, nkvec
      use periodic, only: rknorm, rkvec, rkvec_shift
      use pworbital, only: c_im, c_ip, c_rm, c_rp, isortg, isortk, ngorb
      use contrl_file,    only: ounit
      use coefs, only: norb
      use precision_kinds, only: dp
      implicit none

      integer :: iband, iel, ig, ik, ikvec
      integer :: iorb, iv, jorb, k
      real(dp) :: cos_g, cos_ip, cos_k, cos_rp, dcos_g
      real(dp) :: dcos_k, ddcos_g, ddcos_ip, ddcos_k
      real(dp) :: ddcos_rp, ddsin_g, ddsin_im, ddsin_k
      real(dp) :: ddsin_rm, dsin_g, dsin_k, sin_g
      real(dp) :: sin_im, sin_k, sin_rm
      real(dp), dimension(3,nelec) :: x
      real(dp), dimension(nelec,*) :: orb
      real(dp), dimension(3,nelec,*) :: dorb
      real(dp), dimension(nelec,*) :: ddorb
      real(dp), dimension(3) :: dcos_rp
      real(dp), dimension(3) :: dsin_rm
      real(dp), dimension(3) :: dcos_ip
      real(dp), dimension(3) :: dsin_im








     &,cos_g(NGVECX),sin_g(NGVECX),dcos_g(3,NGVECX),dsin_g(3,NGVECX)
     &,ddcos_g(NGVECX),ddsin_g(NGVECX)
     &,cos_k(IVOL_RATIO),sin_k(IVOL_RATIO),dcos_k(3,IVOL_RATIO),dsin_k(3,IVOL_RATIO)
     &,ddcos_k(IVOL_RATIO),ddsin_k(IVOL_RATIO)

      do iorb=1,norb
        do iel=1,nelec
          orb(iel,iorb)=0
          dorb(1,iel,iorb)=0
          dorb(2,iel,iorb)=0
          dorb(3,iel,iorb)=0
          ddorb(iel,iorb)=0
        enddo
      enddo

      do iel=1,nelec

c compute cos(g.r), sin(g.r) and derivatives
c     call cossin_psi_g(glatt,gnorm,igmult,ngnorm_orb,gvec,igvec,ngvec_orb,x,nelec,ng1d,cos_g,sin_g
      call cossin_psi_g(glatt,gnorm,igmult,ngnorm_orb,gvec,igvec,ngvec_orb,x(1,iel),iel,ng1d,cos_g,sin_g
     &,dcos_g,dsin_g,ddcos_g,ddsin_g,rkvec_shift)

c     write(ounit,'(''cos_g,sin_g,dcos_g,dsin_g,ddcos_g,ddsin_g='',30f9.4)')
c    &cos_g(1,1),sin_g(1,1),(dcos_g(k,1,1),k=1,3),(dsin_g(k,1,1),k=1,3),ddcos_g(1,1),ddsin_g(1,1)
c     write(ounit,'(''cos_g,sin_g,dcos_g,dsin_g,ddcos_g,ddsin_g='',30f9.4)')
c    &cos_g(1,2),sin_g(1,2),(dcos_g(k,1,2),k=1,3),(dsin_g(k,1,2),k=1,3),ddcos_g(1,2),ddsin_g(1,2)

c compute cos(k.r), sin(k.r) and derivatives
c     call cossin_psi_k(glatt_sim,rknorm,rkvec,kvec,nkvec,x,nelec,ng1d_sim,cos_k,sin_k
      call cossin_psi_k(glatt_sim,rknorm,rkvec,kvec,nkvec,x(1,iel),iel,ng1d_sim,cos_k,sin_k
     &,dcos_k,dsin_k,ddcos_k,ddsin_k,rkvec_shift)

      if(ipr.ge.4) then
        write(ounit,'(''rkvec'',9f9.5)') ((rkvec(k,ikvec),k=1,3),ikvec=1,nkvec)
        write(ounit,'(''cos_k,sin_k,dcos_k,dsin_k,ddcos_k,ddsin_k='',10f9.4)')
     &  (cos_k(ik),sin_k(ik),(dcos_k(k,ik),k=1,3),(dsin_k(k,ik),k=1,3),ddcos_k(ik),ddsin_k(ik),ik=1,nkvec)
      endif

        iorb=0
        jorb=0
        do ikvec=1,nkvec
          do iband=1,nband(ikvec)
            jorb=jorb+1

            cos_rp=0
            sin_rm=0
            cos_ip=0
            sin_im=0
            do k=1,3
              dcos_rp(k)=0
              dsin_rm(k)=0
              dcos_ip(k)=0
              dsin_im(k)=0
            enddo
            ddcos_rp=0
            ddsin_rm=0
            ddcos_ip=0
            ddsin_im=0
c           do 80 iv=2,ngorb(ikvec)
c             ig=isortg(iv,ikvec)
            do iv=2,ngvec_orb
              ig=iv
              cos_rp=cos_rp+cos_g(ig)*c_rp(iv,jorb)
              sin_rm=sin_rm+sin_g(ig)*c_rm(iv,jorb)
              cos_ip=cos_ip+cos_g(ig)*c_ip(iv,jorb)
              sin_im=sin_im+sin_g(ig)*c_im(iv,jorb)
c             write(ounit,'(''iel,ig,jorb,cos_rp,cos_g(ig),c_rp(iv,jorb)'',3i5,9d12.4)')
c    & iel,ig,jorb,cos_rp,cos_g(ig),c_rp(iv,jorb)
              do k=1,3
                dcos_rp(k)=dcos_rp(k)+dcos_g(k,ig)*c_rp(iv,jorb)
                dsin_rm(k)=dsin_rm(k)+dsin_g(k,ig)*c_rm(iv,jorb)
                dcos_ip(k)=dcos_ip(k)+dcos_g(k,ig)*c_ip(iv,jorb)
                dsin_im(k)=dsin_im(k)+dsin_g(k,ig)*c_im(iv,jorb)
              enddo
              ddcos_rp=ddcos_rp+ddcos_g(ig)*c_rp(iv,jorb)
              ddsin_rm=ddsin_rm+ddsin_g(ig)*c_rm(iv,jorb)
              ddcos_ip=ddcos_ip+ddcos_g(ig)*c_ip(iv,jorb)
              ddsin_im=ddsin_im+ddsin_g(ig)*c_im(iv,jorb)
            enddo

c           write(ounit,'(''dcos_k(k,ikvec),dsin_k(k,ikvec),dcos_rp(k),dsin_rm(k),dsin_im(k),dcos_ip(k)'',30f9.5)')
c    &(dcos_k(k,ikvec),k=1,3),(dsin_k(k,ikvec),k=1,3),(dcos_rp(k),k=1,3),(dsin_rm(k),k=1,3),(dsin_im(k),k=1,3),(dcos_ip(k),k=1,3
c    &)

            if(k_inv(ikvec).eq.2. .or. ireal_imag(iorb+1).eq.1) then

            iorb=iorb+1

            if(ipr.ge.4) write(ounit,'(''1iorb,ireal_imag(iorb)'',9i5)') iorb,ireal_imag(iorb)

            orb(iel,iorb)=cos_k(ikvec)*(c_rp(1,jorb)+cos_rp-sin_im)
     &                 -sin_k(ikvec)*(c_ip(1,jorb)+sin_rm+cos_ip)
            if(ipr.ge.4) write(ounit,'(''1orb(iel,iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb)
     &,cos_rp,cos_ip,sin_rm,sin_im='',2i5,20d12.4)')
     & iel,iorb,orb(iel,iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb),cos_rp,cos_ip,sin_rm,sin_im

            do k=1,3
              dorb(k,iel,iorb)=dcos_k(k,ikvec)*(c_rp(1,jorb)+cos_rp-sin_im)
     &                        -dsin_k(k,ikvec)*(c_ip(1,jorb)+sin_rm+cos_ip)
     &                        +cos_k(ikvec)*(dcos_rp(k)-dsin_im(k))
     &                        -sin_k(ikvec)*(dsin_rm(k)+dcos_ip(k))
            enddo

            ddorb(iel,iorb)=ddcos_k(ikvec)*(c_rp(1,jorb)+cos_rp-sin_im)
     &                     -ddsin_k(ikvec)*(c_ip(1,jorb)+sin_rm+cos_ip)
     &                     +cos_k(ikvec)*(ddcos_rp-ddsin_im)
     &                     -sin_k(ikvec)*(ddsin_rm+ddcos_ip)
c      write(ounit,'(''ddcos_k(ikvec),ddsin_k(ikvec),cos_k(ikvec),ddcos_rp,sin_k(ikvec),ddsin_rm='',9f9.4)')
c    &ddcos_k(ikvec),ddsin_k(ikvec),cos_k(ikvec),ddcos_rp,sin_k(ikvec),ddsin_rm

c           write(ounit,'(''orb'',2i5,9d12.4)') iel,iorb,orb(iel,iorb),(dorb(k,iel,iorb),k=1,3)
            do k=1,3
              ddorb(iel,iorb)=ddorb(iel,iorb)
     &                       +2*(dcos_k(k,ikvec)*(dcos_rp(k)-dsin_im(k))
     &                          -dsin_k(k,ikvec)*(dsin_rm(k)+dcos_ip(k)))
            enddo

            if(k_inv(ikvec).eq.1) goto 130
            endif

            if(iorb.lt.norb) then
c           if(k_inv(ikvec).eq.2. .or. ireal_imag(iorb+1).eq.2) then

            iorb=iorb+1


            if(ipr.ge.4) write(ounit,'(''2iorb,ireal_imag(iorb)'',9i5)') iorb,ireal_imag(iorb)

            orb(iel,iorb)=cos_k(ikvec)*(c_ip(1,jorb)+sin_rm+cos_ip)
     &                   +sin_k(ikvec)*(c_rp(1,jorb)+cos_rp-sin_im)
            if(ipr.ge.4) write(ounit,'(''2orb(iel,iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb)
     &,cos_rp,cos_ip,sin_rm,sin_im='',2i5,20d12.4)')
     & iel,iorb,orb(iel,iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb),cos_rp,cos_ip,sin_rm,sin_im

            do k=1,3
              dorb(k,iel,iorb)=dcos_k(k,ikvec)*(c_ip(1,jorb)+sin_rm+cos_ip)
     &                        +dsin_k(k,ikvec)*(c_rp(1,jorb)+cos_rp-sin_im)
     &                        +cos_k(ikvec)*(dsin_rm(k)+dcos_ip(k))
     &                        +sin_k(ikvec)*(dcos_rp(k)-dsin_im(k))
            enddo

            ddorb(iel,iorb)=ddcos_k(ikvec)*(c_ip(1,jorb)+sin_rm+cos_ip)
     &                     +ddsin_k(ikvec)*(c_rp(1,jorb)+cos_rp-sin_im)
     &                     +cos_k(ikvec)*(ddsin_rm+ddcos_ip)
     &                     +sin_k(ikvec)*(ddcos_rp-ddsin_im)
c           write(ounit,'(''orb2'',2i5,9d12.4)') iel,iorb,orb(iel,iorb),(dorb(k,iel,iorb),k=1,3)
            do k=1,3
              ddorb(iel,iorb)=ddorb(iel,iorb)
     &                       +2*(dcos_k(k,ikvec)*(dsin_rm(k)+dcos_ip(k))
     &                          +dsin_k(k,ikvec)*(dcos_rp(k)-dsin_im(k)))
            enddo
c           endif
            endif

  130       continue
          enddo
        enddo
      enddo

      if(ipr.ge.4) write(ounit,'(i4,'' electrons placed in'',i4,'' orbitals'')') nelec,iorb

      return
      end
c-----------------------------------------------------------------------

      subroutine orbitals_pw_grade(iel,x,orb,dorb,ddorb)
c Calculate pw orbitals, gradient and laplacian for electron iel.
c isortg could be used to map g-vectors from iv to ig and
c isortk could be used to map k-vectors.
c At present it is assumed that both g- and k-vectors are in the correct order.

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use ewald_mod, only: IVOL_RATIO
      use ewald_mod, only: NGVECX
      use const, only: nelec, ipr
      use periodic, only: glatt
      use periodic, only: glatt_sim, gnorm, gvec, igmult, igvec
      use periodic, only: ireal_imag, k_inv, kvec, nband, ng1d, ng1d_sim, ngnorm_orb
      use periodic, only: ngvec_orb, nkvec
      use periodic, only: rknorm, rkvec, rkvec_shift
      use pworbital, only: c_im, c_ip, c_rm, c_rp, isortg, isortk, ngorb
      use contrl_file,    only: ounit
      use coefs, only: norb
      use precision_kinds, only: dp
      use vmc_mod, only: norb_tot
      implicit none

      integer :: iband, iel, ig, ikvec, iorb
      integer :: iv, jorb, k
      real(dp) :: cos_ip, cos_rp, ddcos_ip, ddcos_rp, ddsin_im
      real(dp) :: ddsin_rm, sin_im, sin_rm
      real(dp), dimension(3) :: x
      real(dp), dimension(*) :: orb
      real(dp), dimension(norb_tot,3) :: dorb
      real(dp), dimension(*) :: ddorb
      real(dp), dimension(3) :: dcos_rp
      real(dp), dimension(3) :: dsin_rm
      real(dp), dimension(3) :: dcos_ip
      real(dp), dimension(3) :: dsin_im
      real(dp), dimension(NGVECX) :: cos_g
      real(dp), dimension(NGVECX) :: sin_g
      real(dp), dimension(3,NGVECX) :: dcos_g
      real(dp), dimension(3,NGVECX) :: dsin_g
      real(dp), dimension(NGVECX) :: ddcos_g
      real(dp), dimension(NGVECX) :: ddsin_g
      real(dp), dimension(IVOL_RATIO) :: cos_k
      real(dp), dimension(IVOL_RATIO) :: sin_k
      real(dp), dimension(3,IVOL_RATIO) :: dcos_k
      real(dp), dimension(3,IVOL_RATIO) :: dsin_k
      real(dp), dimension(IVOL_RATIO) :: ddcos_k
      real(dp), dimension(IVOL_RATIO) :: ddsin_k


      do k=1,3
          do iorb=1,norb
            dorb(iorb,k)=0
          enddo
      enddo

      do iorb=1,norb
         orb(iorb)=0
         ddorb(iorb)=0
      enddo

c     do 130 iel=1,nelec

c compute cos(g.r), sin(g.r) and derivatives
c     call cossin_psi_g(glatt,gnorm,igmult,ngnorm_orb,gvec,igvec,ngvec_orb,x,nelec,ng1d,cos_g,sin_g
      call cossin_psi_g(glatt,gnorm,igmult,ngnorm_orb,gvec,igvec,ngvec_orb,x,iel,ng1d,cos_g,sin_g
     &,dcos_g,dsin_g,ddcos_g,ddsin_g,rkvec_shift)

c     write(ounit,'(''cos_g,sin_g,dcos_g,dsin_g,ddcos_g,ddsin_g='',30f9.4)')
c    &cos_g(1,1),sin_g(1,1),(dcos_g(k,1,1),k=1,3),(dsin_g(k,1,1),k=1,3),ddcos_g(1,1),ddsin_g(1,1)
c     write(ounit,'(''cos_g,sin_g,dcos_g,dsin_g,ddcos_g,ddsin_g='',30f9.4)')
c    &cos_g(1,2),sin_g(1,2),(dcos_g(k,1,2),k=1,3),(dsin_g(k,1,2),k=1,3),ddcos_g(1,2),ddsin_g(1,2)

c compute cos(k.r), sin(k.r) and derivatives
c     call cossin_psi_k(glatt_sim,rknorm,rkvec,kvec,nkvec,x,nelec,ng1d_sim,cos_k,sin_k
      call cossin_psi_k(glatt_sim,rknorm,rkvec,kvec,nkvec,x,iel,ng1d_sim,cos_k,sin_k
     &,dcos_k,dsin_k,ddcos_k,ddsin_k,rkvec_shift)

      if(ipr.ge.4) then
        write(ounit,'(''cos_k,sin_k,dcos_k,dsin_k,ddcos_k,ddsin_k='',30f9.4)')
     &  cos_k(1),sin_k(1),(dcos_k(k,1),k=1,3),(dsin_k(k,1),k=1,3),ddcos_k(1),ddsin_k(1)
      endif

c     do 130 iel=1,nelec

        iorb=0
        jorb=0
        do ikvec=1,nkvec
          do iband=1,nband(ikvec)
            jorb=jorb+1

            cos_rp=0
            sin_rm=0
            cos_ip=0
            sin_im=0
            do k=1,3
              dcos_rp(k)=0
              dsin_rm(k)=0
              dcos_ip(k)=0
              dsin_im(k)=0
            enddo
            ddcos_rp=0
            ddsin_rm=0
            ddcos_ip=0
            ddsin_im=0
c           do 80 iv=2,ngorb(ikvec)
c             ig=isortg(iv,ikvec)
            do iv=2,ngvec_orb
              ig=iv
              cos_rp=cos_rp+cos_g(ig)*c_rp(iv,jorb)
              sin_rm=sin_rm+sin_g(ig)*c_rm(iv,jorb)
              cos_ip=cos_ip+cos_g(ig)*c_ip(iv,jorb)
              sin_im=sin_im+sin_g(ig)*c_im(iv,jorb)
c             write(ounit,'(''iel,ig,jorb,cos_rp,cos_g(ig),c_rp(iv,jorb)'',3i5,9d12.4)')
c    & iel,ig,jorb,cos_rp,cos_g(ig),c_rp(iv,jorb)
              do k=1,3
                dcos_rp(k)=dcos_rp(k)+dcos_g(k,ig)*c_rp(iv,jorb)
                dsin_rm(k)=dsin_rm(k)+dsin_g(k,ig)*c_rm(iv,jorb)
                dcos_ip(k)=dcos_ip(k)+dcos_g(k,ig)*c_ip(iv,jorb)
                dsin_im(k)=dsin_im(k)+dsin_g(k,ig)*c_im(iv,jorb)
              enddo
              ddcos_rp=ddcos_rp+ddcos_g(ig)*c_rp(iv,jorb)
              ddsin_rm=ddsin_rm+ddsin_g(ig)*c_rm(iv,jorb)
              ddcos_ip=ddcos_ip+ddcos_g(ig)*c_ip(iv,jorb)
              ddsin_im=ddsin_im+ddsin_g(ig)*c_im(iv,jorb)
            enddo

c           write(ounit,'(''dcos_k(k,ikvec),dsin_k(k,ikvec),dcos_rp(k),dsin_rm(k),dsin_im(k),dcos_ip(k)'',30f9.5)')
c    &(dcos_k(k,ikvec),k=1,3),(dsin_k(k,ikvec),k=1,3),(dcos_rp(k),k=1,3),(dsin_rm(k),k=1,3),(dsin_im(k),k=1,3),(dcos_ip(k),k=1,3
c    &)

            if(k_inv(ikvec).eq.2. .or. ireal_imag(iorb+1).eq.1) then

            iorb=iorb+1

            orb(iorb)=cos_k(ikvec)*(c_rp(1,jorb)+cos_rp-sin_im)
     &               -sin_k(ikvec)*(c_ip(1,jorb)+sin_rm+cos_ip)
            if(ipr.ge.4) write(ounit,'(''1orb(iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb)
     &,cos_rp,cos_ip,sin_rm,sin_im='',2i5,20d12.4)')
     & iel,iorb,orb(iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb),cos_rp,cos_ip,sin_rm,sin_im

            do k=1,3
              dorb(iorb,k)=dcos_k(k,ikvec)*(c_rp(1,jorb)+cos_rp-sin_im)
     &                    -dsin_k(k,ikvec)*(c_ip(1,jorb)+sin_rm+cos_ip)
     &                    +cos_k(ikvec)*(dcos_rp(k)-dsin_im(k))
     &                    -sin_k(ikvec)*(dsin_rm(k)+dcos_ip(k))
            enddo

            ddorb(iorb)=ddcos_k(ikvec)*(c_rp(1,jorb)+cos_rp-sin_im)
     &                 -ddsin_k(ikvec)*(c_ip(1,jorb)+sin_rm+cos_ip)
     &                 +cos_k(ikvec)*(ddcos_rp-ddsin_im)
     &                 -sin_k(ikvec)*(ddsin_rm+ddcos_ip)
c      write(ounit,'(''ddcos_k(ikvec),ddsin_k(ikvec),cos_k(ikvec),ddcos_rp,sin_k(ikvec),ddsin_rm='',9f9.4)')
c    &ddcos_k(ikvec),ddsin_k(ikvec),cos_k(ikvec),ddcos_rp,sin_k(ikvec),ddsin_rm

c           write(ounit,'(''orb'',2i5,9d12.4)') iel,iorb,orb(iorb),(dorb(iorb,k),k=1,3)
            do k=1,3
              ddorb(iorb)=ddorb(iorb)
     &                   +2*(dcos_k(k,ikvec)*(dcos_rp(k)-dsin_im(k))
     &                      -dsin_k(k,ikvec)*(dsin_rm(k)+dcos_ip(k)))
            enddo

            if(k_inv(ikvec).eq.1) goto 130

            endif

            if(iorb.lt.norb) then
            if(k_inv(ikvec).eq.2. .or. ireal_imag(iorb+1).eq.2) then
c           if(k_inv(ikvec).eq.1) goto 130

            iorb=iorb+1

            orb(iorb)=cos_k(ikvec)*(c_ip(1,jorb)+sin_rm+cos_ip)
     &               +sin_k(ikvec)*(c_rp(1,jorb)+cos_rp-sin_im)

            do k=1,3
              dorb(iorb,k)=dcos_k(k,ikvec)*(c_ip(1,jorb)+sin_rm+cos_ip)
     &                    +dsin_k(k,ikvec)*(c_rp(1,jorb)+cos_rp-sin_im)
     &                    +cos_k(ikvec)*(dsin_rm(k)+dcos_ip(k))
     &                    +sin_k(ikvec)*(dcos_rp(k)-dsin_im(k))
            enddo

            ddorb(iorb)=ddcos_k(ikvec)*(c_ip(1,jorb)+sin_rm+cos_ip)
     &                 +ddsin_k(ikvec)*(c_rp(1,jorb)+cos_rp-sin_im)
     &                 +cos_k(ikvec)*(ddsin_rm+ddcos_ip)
     &                 +sin_k(ikvec)*(ddcos_rp-ddsin_im)
c           write(ounit,'(''orb2'',2i5,9d12.4)') iel,iorb,orb(iorb),(dorb(iorb,k),k=1,3)
            do k=1,3
              ddorb(iorb)=ddorb(iorb)
     &                   +2*(dcos_k(k,ikvec)*(dsin_rm(k)+dcos_ip(k))
     &                      +dsin_k(k,ikvec)*(dcos_rp(k)-dsin_im(k)))
            enddo

            endif
            endif

  130       continue
          enddo
        enddo

      if(ipr.ge.4) write(ounit,'(i4,'' electrons placed in'',i4,'' orbitals'')') nelec,iorb

      return
      end
