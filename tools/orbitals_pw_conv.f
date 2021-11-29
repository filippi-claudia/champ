      program orbitals_pw_conv
c Written by Cyrus Umrigar
c Convert dft pw orbitals into ones suitable for input to qmc

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'ewald.h'

      common /pworbital/c_rp(NGVECX,MORB),c_rm(NGVECX,MORB),c_ip(NGVECX,MORB)
     &,c_im(NGVECX,MORB),ngorb(MORB),isortg(NGVECX,MORB),isortk(IVOL_RATIO),icmplx
      common /periodic/ rlatt(3,3),glatt(3,3),rlatt_sim(3,3),glatt_sim(3,3)
     &,rlatt_inv(3,3),rlatt_sim_inv(3,3),glatt_inv(3,3)
     &,cutr,cutr_sim,cutg,cutg_sim,cutg_big,cutg_sim_big
     &,igvec(3,NGVEC_BIGX),gvec(3,NGVEC_BIGX),gnorm(NGNORM_BIGX),igmult(NGNORM_BIGX)
     &,igvec_sim(3,NGVEC_SIM_BIGX),gvec_sim(3,NGVEC_SIM_BIGX),gnorm_sim(NGNORM_SIM_BIGX),igmult_sim(NGNORM_SIM_BIGX)
     &,rkvec_shift(3),kvec(3,IVOL_RATIO),rkvec(3,IVOL_RATIO),rknorm(IVOL_RATIO)
     &,k_inv(IVOL_RATIO),nband(IVOL_RATIO),ireal_imag(MORB)
     &,znuc_sum,znuc2_sum,vcell,vcell_sim
     &,ngnorm,ngvec,ngnorm_sim,ngvec_sim,ngnorm_orb,ngvec_orb,nkvec
     &,ngnorm_big,ngvec_big,ngnorm_sim_big,ngvec_sim_big
     &,ng1d(3),ng1d_sim(3),npoly,ncoef,np,isrange

      dimension igvec_dft(3,NGVEC_BIGX),iwgvec(NGVEC_BIGX),c_real(NGVEC_BIGX),c_imag(NGVEC_BIGX)

c Warning: For the moment we assume that orbitals_pw_tm contains only the bands we want to keep.
c for each k-pt.

      open(1,file='gvectors_qmc')
      read(1,*) ngvec
      if(ngvec.gt.NGVECX) call fatal_error ('ngvec>NGVECX')
      do i=1,ngvec
        read(1,*) (igvec(k,i),k=1,3)
      enddo
c  10   write(6,*) (igvec(k,i),k=1,3)

c     open(2,file='gvectors')
      open(30,file='orbitals_pw_tm')
      read(30,*) ngvec_dft
      do i=1,ngvec_dft
        read(30,*,end=22) (igvec_dft(k,i),k=1,3)
      enddo
c  20   read(30,*,end=22) junk,junk,(igvec_dft(k,i),k=1,3)
c  20   write(6,*) (igvec_dft(k,i),k=1,3)

      read(30,*) nkvec
      open(3,file='orbitals_pw')
      write(3,'(2i4,'' nkvec,ngvec'')') nkvec,ngvec
      if(nkvec.gt.IVOL_RATIO) call fatal_error ('nkvec>IVOL_RATIO')

   22 do i=1,nkvec
        read(30,*) ikpt,nband(i),ng,(rkvec(k,i),k=1,3)
        write(3,'(2i4,3f9.5'' ikpt, nband, rkvec(in recip. lat. units)'')') i,nband(i),(rkvec(k,i),k=1,3)
        if(ng.gt.NGVEC_BIGX) call fatal_error ('ng>NGVEC_BIGX')
        read(30,*) (iwgvec(j),j=1,ng)
        do iband=1,nband(i)
          read(30,*) iband_tmp,eig
          read(30,*) (c_real(ig),c_imag(ig),ig=1,ng)
          do igv=1,ngvec
            c_rp(igv,i)=0
            c_rm(igv,i)=0
            c_ip(igv,i)=0
            c_im(igv,i)=0
          enddo
          ifound=0
          do igv=1,ngvec
            if(igvec(1,igv).eq.0 .and. igvec(2,igv).eq.0 .and.  igvec(3,igv).eq.0) then
              isign_min=1
             else
              isign_min=-1
            endif
            do ig=1,ng
              do isign=1,isign_min,-2
                norm=0
c               write(6,*) (igvec_dft(k,iwgvec(ig)),k=1,3),(igvec(k,igv),k=1,3)
                do k=1,3
                  norm=norm+(igvec_dft(k,iwgvec(ig))-isign*igvec(k,igv))**2
                enddo
                if(norm.eq.0) then
                  ifound=ifound+1
                  c_rp(igv,i)=c_rp(igv,i)+c_real(ig)
                  c_rm(igv,i)=c_rm(igv,i)+c_real(ig)*isign
                  c_ip(igv,i)=c_ip(igv,i)+c_imag(ig)
                  c_im(igv,i)=c_im(igv,i)+c_imag(ig)*isign
                endif
              enddo
            enddo
          enddo
         if(ifound.ne.ng) then
           write(3,'(''ng,ifound='',2i5)') ng,ifound
           if(ifound.lt.ng) call fatal_error ('probably need more g-vectors in gvectors_qmc file')
           if(ifound.gt.ng) call fatal_error ('probably need more g-vectors in gvectors file')
         endif
         write(3,'(2i5,f10.6,'' ikpt, iband, eig (Ry)'')') i,iband,eig
         write(3,'(1p4d22.14)') (c_rp(igv,i),c_rm(igv,i),c_ip(igv,i),c_im(igv,i),igv=1,ngvec)
        enddo
      enddo

      call fatal_error('')
      end
