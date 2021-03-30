      subroutine orbitals_pwe(iel,x,orb)
c Written by Cyrus Umrigar
c Calculate pw orbitals.
c isortg could be used to map g-vectors from iv to ig and
c isortk could be used to map k-vectors.
c At present it is assumed that both g- and k-vectors are in the correct order.

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use ewald_mod, only: IVOL_RATIO
      use ewald_mod, only: NGVECX
      use vmc_mod, only: MELEC
      use const, only: nelec, ipr
      use periodic, only: glatt
      use periodic, only: glatt_sim, gnorm, gvec, igmult, igvec
      use periodic, only: ireal_imag, k_inv, kvec, nband, ng1d, ng1d_sim, ngnorm_orb
      use periodic, only: ngvec_orb, nkvec
      use periodic, only: rknorm, rkvec, rkvec_shift
      use pworbital, only: c_im, c_ip, c_rm, c_rp, isortg, isortk, ngorb
      use coefs, only: norb

      implicit real*8(a-h,o-z)



      dimension x(3),orb(*)
      dimension cos_g(NGVECX),sin_g(NGVECX),dcos_g(3,NGVECX),dsin_g(3,NGVECX)
     &,ddcos_g(NGVECX),ddsin_g(NGVECX)
     &,cos_k(IVOL_RATIO),sin_k(IVOL_RATIO),dcos_k(3,IVOL_RATIO),dsin_k(3,IVOL_RATIO)
     &,ddcos_k(IVOL_RATIO),ddsin_k(IVOL_RATIO)

      do 5 iorb=1,norb
c       do 5 iel=1,nelec
    5     orb(iorb)=0

c     do 130 iel=1,nelec

c compute cos(g.r), sin(g.r) and derivatives
c     call cossin_psi_g(glatt,gnorm,igmult,ngnorm_orb,gvec,igvec,ngvec_orb,x,nelec,ng1d,cos_g,sin_g
c    &,dcos_g,dsin_g,ddcos_g,ddsin_g,rkvec_shift,0)
      call cossin_psi_g(glatt,gnorm,igmult,ngnorm_orb,gvec,igvec,ngvec_orb,x,iel,ng1d,cos_g,sin_g
     &,dcos_g,dsin_g,ddcos_g,ddsin_g,rkvec_shift)

c     write(6,'(''cos_g,sin_g,dcos_g,dsin_g,ddcos_g,ddsin_g='',30f9.4)')
c    &cos_g(1,1),sin_g(1,1),(dcos_g(k,1,1),k=1,3),(dsin_g(k,1,1),k=1,3),ddcos_g(1,1),ddsin_g(1,1)
c     write(6,'(''cos_g,sin_g,dcos_g,dsin_g,ddcos_g,ddsin_g='',30f9.4)')
c    &cos_g(1,2),sin_g(1,2),(dcos_g(k,1,2),k=1,3),(dsin_g(k,1,2),k=1,3),ddcos_g(1,2),ddsin_g(1,2)

c compute cos(k.r), sin(k.r) and derivatives
c     call cossin_psi_k(glatt_sim,rknorm,rkvec,kvec,nkvec,x,nelec,ng1d_sim,cos_k,sin_k
      call cossin_psi_k(glatt_sim,rknorm,rkvec,kvec,nkvec,x,iel,ng1d_sim,cos_k,sin_k
     &,dcos_k,dsin_k,ddcos_k,ddsin_k,rkvec_shift)

c     write(6,'(''cos_k,sin_k,dcos_k,dsin_k,ddcos_k,ddsin_k='',30f9.4)')
c    &cos_k(1,1),sin_k(1,1),(dcos_k(k,1,1),k=1,3),(dsin_k(k,1,1),k=1,3),ddcos_k(1,1),ddsin_k(1,1)


        iorb=0
        jorb=0
        do 130 ikvec=1,nkvec
          do 130 iband=1,nband(ikvec)
            jorb=jorb+1

            cos_rp=0
            sin_rm=0
            cos_ip=0
            sin_im=0
c           do 80 iv=2,ngorb(ikvec)
c             ig=isortg(iv,ikvec)
            do 80 iv=2,ngvec_orb
              ig=iv
              cos_rp=cos_rp+cos_g(ig)*c_rp(iv,jorb)
              sin_rm=sin_rm+sin_g(ig)*c_rm(iv,jorb)
              cos_ip=cos_ip+cos_g(ig)*c_ip(iv,jorb)
              sin_im=sin_im+sin_g(ig)*c_im(iv,jorb)
c             write(6,'(''iel,ig,jorb,cos_rp,cos_g(ig),c_rp(iv,jorb)'',3i5,9d12.4)')
c    & iel,ig,jorb,cos_rp,cos_g(ig),c_rp(iv,jorb)
   80         continue

c           write(6,'(''dcos_k(k,iel,ikvec),dsin_k(k,iel,ikvec),dcos_rp(k),dsin_rm(k),dsin_im(k),dcos_ip(k)'',30f9.5)')
c    &(dcos_k(k,iel,ikvec),k=1,3),(dsin_k(k,iel,ikvec),k=1,3),(dcos_rp(k),k=1,3),(dsin_rm(k),k=1,3),(dsin_im(k),k=1,3),(dcos_ip(k),k=1,3
c    &)

            if(k_inv(ikvec).eq.2. .or. ireal_imag(iorb+1).eq.1) then

            iorb=iorb+1

            orb(iorb)=cos_k(ikvec)*(c_rp(1,jorb)+cos_rp-sin_im)
     &               -sin_k(ikvec)*(c_ip(1,jorb)+sin_rm+cos_ip)
            if(ipr.ge.4) then
            write(6,'(''1x='',3f9.5)') (x(k),k=1,3)
            write(6,'(''21orb(iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb)
     &,cos_rp,cos_ip,sin_rm,sin_im='',2i5,20d12.4)')
     & iel,iorb,orb(iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb),cos_rp,cos_ip,sin_rm,sin_im
            endif

c      write(6,'(''ddcos_k(ikvec),ddsin_k(ikvec),cos_k(ikvec),ddcos_rp,sin_k(ikvec),ddsin_rm='',9f9.4)')
c    &ddcos_k(ikvec),ddsin_k(ikvec),cos_k(ikvec),ddcos_rp,sin_k(ikvec),ddsin_rm

c           write(6,'(''orb'',2i5,9d12.4)') iel,iorb,orb(iorb),(dorb(k,iel,iorb),k=1,3)

            if(k_inv(ikvec).eq.1) goto 130

            endif

            if(iorb.lt.norb) then
c           if(k_inv(ikvec).eq.2. .or. ireal_imag(iorb+1).eq.2) then

            iorb=iorb+1

            orb(iorb)=cos_k(ikvec)*(c_ip(1,jorb)+sin_rm+cos_ip)
     &               +sin_k(ikvec)*(c_rp(1,jorb)+cos_rp-sin_im)
            if(ipr.ge.4) write(6,'(''22orb(iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb)
     &,cos_rp,cos_ip,sin_rm,sin_im='',2i5,20d12.4)')
     & iel,iorb,orb(iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb),cos_rp,cos_ip,sin_rm,sin_im

c           endif
            endif

  130       continue

      if(ipr.ge.4) write(6,'(i4,'' electrons placed in'',i4,'' orbitals'')') nelec,iorb

      return
      end
