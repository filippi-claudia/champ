c$$$      subroutine pw_setup_input
c$$$      use periodic, only: cutg, cutg_big, cutg_sim, cutg_sim_big, cutr, cutr_sim, glatt,
c$$$    &glatt_inv, glatt_sim, gnorm, gnorm_sim, gvec, gvec_sim, igmult, igmult_sim, igvec, igvec_sim,
c$$$    &ireal_imag, isrange, k_inv, kvec, nband, ncoef_per, ng1d, ng1d_sim, ngnorm, ngnorm_big, ngnorm_orb,
c$$$    &ngnorm_sim, ngnorm_sim_big, ngvec, ngvec_big, ngvec_orb, ngvec_sim, ngvec_sim_big, nkvec,
c$$$    &np, npoly, rknorm, rkvec, rkvec_shift, rlatt, rlatt_inv, rlatt_sim, rlatt_sim_inv, vcell,
c$$$    &vcell_sim, znuc2_sum, znuc_sum
c$$$     implicit real*8(a-h,o-z)

c$$$
c$$$      include 'vmc.h'
c$$$      include 'ewald.h'
c$$$      include 'inputflags.h'
c$$$
c$$$c  npoly,np,cutg,cutg_sim,cutg_big,cutg_sim_big
c$$$c  alattice   lattice constant to multiply rlatt
c$$$c  rlatt      lattice vectors of primitive cell
c$$$c  rlatt_sim  lattice vectors of simulation cell
c$$$c  rkvec_shift k-shift for generating k-vector lattice
c$$$
c$$$c$$$     &,rlatt_inv(3,3),rlatt_sim_inv(3,3),glatt_inv(3,3)
c$$$     &,cutr,cutr_sim,cutg,cutg_sim,cutg_big,cutg_sim_big
c$$$     &,igvec(3,NGVEC_BIGX),gvec(3,NGVEC_BIGX),gnorm(NGNORM_BIGX),igmult(NGNORM_BIGX)
c$$$     &,igvec_sim(3,NGVEC_SIM_BIGX),gvec_sim(3,NGVEC_SIM_BIGX),gnorm_sim(NGNORM_SIM_BIGX),igmult_sim(NGNORM_SIM_BIGX)
c$$$     &,rkvec_shift(3),kvec(3,IVOL_RATIO),rkvec(3,IVOL_RATIO),rknorm(IVOL_RATIO)
c$$$     &,k_inv(IVOL_RATIO),nband(IVOL_RATIO),ireal_imag(norb_tot)
c$$$     &,znuc_sum,znuc2_sum,vcell,vcell_sim
c$$$     &,ngnorm,ngvec,ngnorm_sim,ngvec_sim,ngnorm_orb,ngvec_orb,nkvec
c$$$     &,ngnorm_big,ngvec_big,ngnorm_sim_big,ngvec_sim_big
c$$$     &,ng1d(3),ng1d_sim(3),npoly,ncoef_per,np,isrange
c$$$
c$$$c note that if iperiodic=0, norb is fetched in read_lcao
c$$$      call p2gti('periodic:norb',norb,1)
c$$$
c$$$c npoly is the polynomial order for short-range part
c$$$      call p2gti('periodic:npoly',npoly,1)
c$$$      call p2gti('periodic:np',np,1)
c$$$      ncoef_per=npoly+1
c$$$      if(ncoef_per.gt.NCOEFX) call fatal_error('INPUT: ncoef_per gt NCOEFX')
c$$$
c$$$      call p2gtf('periodic:cutg',cutg,1)
c$$$      call p2gtf('periodic:cutg_sim',cutg_sim,1)
c$$$      call p2gtf('periodic:cutg_big',cutg_big,1)
c$$$      call p2gtf('periodic:cutg_sim_big',cutg_sim_big,1)
c$$$      write(ounit,'(/,''Npoly,np,cutg,cutg_sim,cutg_big,cutg_sim_big'',2i4,9f8.2)')
c$$$     & npoly,np,cutg,cutg_sim,cutg_big,cutg_sim_big
c$$$
c$$$c Lattice vectors fetched in read_lattice
c$$$      write(ounit,'(/,''Lattice basis vectors'',3(/,3f10.6))')
c$$$     & ((rlatt(k,j),k=1,3),j=1,3)
c$$$      write(ounit,'(/,''Simulation lattice basis vectors'',3(/,3f10.6))')
c$$$     & ((rlatt_sim(k,j),k=1,3),j=1,3)
c$$$
c$$$      write(ounit,'(/,''center positions in primitive lattice vector units '',
c$$$     &  ''and in cartesian coordinates'')')
c$$$c Convert center positions from primitive lattice vector units to cartesian coordinates
c$$$      do 5 ic=1,ncent
c$$$        do 3 k=1,3
c$$$    3     cent_tmp(k)=cent(k,ic)
c$$$        do 4 k=1,3
c$$$          cent(k,ic)=0
c$$$          do 4 i=1,3
c$$$    4       cent(k,ic)=cent(k,ic)+cent_tmp(i)*rlatt(k,i)
c$$$    5   write(ounit,'(''center'',i4,1x,''('',3f9.5,'')'',1x,''('',3f9.5,'')'')')
c$$$     &  ic,(cent_tmp(k),k=1,3),(cent(k,ic),k=1,3)
c$$$
c$$$c cutjas already fetched by read_jastrow_parameter, set_ewald would overwrite cutjas
c$$$      cutjas_tmp=cutjas
c$$$      call set_ewald
c$$$
c$$$      return
c$$$      end
c----------------------------------------------------------------------
      module pw_read
      contains
      subroutine do_read_lattice(iu)

      use periodic, only: rkvec_shift, rlatt, rlatt_sim
      use inputflags, only: ilattice
      use contrl_file, only: ounit
      use precision_kinds, only: dp
      implicit none

      integer :: i, itmp, iu, k
      real(dp) :: alattice



      call incpos(iu,itmp,1)
      read(iu,*) alattice
      do i=1,3
        call incpos(iu,itmp,1)
        read(iu,*) (rlatt(k,i),k=1,3)
        do k=1,3
          rlatt(k,i)=rlatt(k,i)*alattice
        enddo
      enddo

c Read the dimensions of the simulation 'cube'
        do i=1,3
          call incpos(iu,itmp,1)
          read(iu,*) (rlatt_sim(k,i),k=1,3)
          do k=1,3
            rlatt_sim(k,i)=rlatt_sim(k,i)*alattice
          enddo
        enddo

c Read k-shift for generating k-vector lattice
      call incpos(iu,itmp,1)
      read(iu,*) (rkvec_shift(k),k=1,3)

      ilattice=1
      call p2chkend(iu, 'lattice')

      end
c----------------------------------------------------------------------
      subroutine read_orb_pw
c Written by Cyrus Umrigar
c Reads in pw basis orbitals that have already been converted to be real.
c Presently not used.

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use ewald_mod, only: IVOL_RATIO
      use periodic, only: nband
      use periodic, only: ngvec, nkvec
      use periodic, only: rkvec
      use pworbital, only: c_im, c_ip, c_rm, c_rp
      use contrl_file, only: ounit
      use coefs, only: norb
      use precision_kinds, only: dp
      implicit none

      integer :: i, i3, ikvec, jorb, k
      real(dp), dimension(3) :: rkvec_tmp







      open(3,file='orbitals_pw')
      read(3,*) nkvec,ngvec
      if(nkvec.gt.IVOL_RATIO) call fatal_error ('nkvec>IVOL_RATIO in read_orb_pw')

      jorb=0
      do i=1,nkvec
        read(3,*) ikvec,nband(i),(rkvec_tmp(k),k=1,3)
        do k=1,3
          if(abs(rkvec_tmp(k)-rkvec(k,i)).gt.1.d-5) call fatal_error ('rkvec_tmp != rkvec in read_orb_pw')
        enddo
        do iband=1,nband(i)
          jorb=jorb+1
          read(3,*) ikvec,ib,eig
          read(3,*) (c_rp(igv,jorb),c_rm(igv,jorb),c_ip(igv,jorb),c_im(igv,jorb),igv=1,ngvec)
        enddo
      enddo

      write(ounit,'(i3,'' orbitals read in from file orbitals_pw'')') jorb
      if(jorb.lt.norb) then
        write(ounit,'(''jorb,norb='',2i5)') jorb,norb
        call fatal_error ('jorb < norb in read_orb_pw')
      endif
      close(3)

      return
      end
c-----------------------------------------------------------------------

      subroutine read_orb_pw_tm
c Written by Cyrus Umrigar
c Reads in pw basis orbitals and convert them to real ones suitable for qmc.
c Warning: At present NGVECX is used to dimension not only quantities that are
c used for the Ewald sum, but also quantities used for the wavefunction pw coefficients.
c There is no reason for the latter to be tied to the former.

c I write out orbitals_pw so that in future runs I can just read that file in
c rather than orbitals_pw_tm, but at the moment I am not using that feature.
c Also, I first write out a temporary fort.3 and then delete it just because
c it is only after one has processed all the k-pts that one knows how big ngvec_orb is.
c However, that causes problems when running with mpi, so comment out that part.

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use ewald_mod, only: IVOL_RATIO
      use ewald_mod, only: NGVECX
      use ewald_mod, only: NGVEC_BIGX
      use vmc_mod, only: norb_tot
      use periodic, only: glatt
      use periodic, only: igmult, igvec
      use periodic, only: ireal_imag, k_inv, kvec, nband, ngnorm, ngnorm_orb
      use periodic, only: ngvec, ngvec_orb, nkvec
      use periodic, only: rknorm, rkvec, rkvec_shift
      use pworbital, only: c_im, c_ip, c_rm, c_rp, icmplx

      use tempor_test, only: c_imag, c_real, igvec_dft, iwgvec, ngg, ngvec_dft, rkvec_tmp, rkvec_tmp2
      use contrl_file, only: ounit
      use coefs, only: norb
      use precision_kinds, only: dp
      use system, only: nelec
      implicit none

      integer :: i, i1, i3, ib, iband
      integer :: iband_tmp, ifound, ig, ig_min
      integer :: igv, ikv, ikvec, in
      integer :: isign, isign_min, j, jorb
      integer :: k, nelec_sav, ng, nkvec_tmp
      integer :: norm, nsum
      real(dp) :: eig, sum
      real(dp) :: sum_abs, units
      real(dp), dimension(3) :: r
      real(dp), dimension(nelec,norb_tot) :: orb
      real(dp), dimension(norb_tot,nelec,3) :: dorb
      real(dp), dimension(norb_tot,nelec) :: ddorb








c Warning: Temporary
c     dimension igvec_dft(3,NGVEC_BIGX),iwgvec(NGVEC_BIGX),c_real(NGVEC_BIGX),c_imag(NGVEC_BIGX)
c    &,rkvec_tmp(3),rkvec_tmp2(3)

c Warning: Temporary
c Warning: why do I print out zeros if I dimension to norb_tot rather than IVOL_RATIO?
c It should be norb_tot


c Warning: For the moment we assume that orbitals_pw_tm contains only the bands we want to keep.
c for each k-pt.

c     open(1,file='gvectors_qmc')
c     read(1,*) ngvec
c     if(ngvec.gt.NGVECX) call fatal_error ('ngvec>NGVECX')
c     do 10 i=1,ngvec
c  10   read(1,*) (igvec(k,i),k=1,3)

c     open(2,file='gvectors')
      open(30,file='orbitals_pw_tm')
      read(30,*) icmplx,ngvec_dft
      do i=1,ngvec_dft
        read(30,*) (igvec_dft(k,i),k=1,3)
      enddo
c  20   read(30,*) junk,junk,(igvec_dft(k,i),k=1,3)
c  20   write(ounit,*) (igvec_dft(k,i),k=1,3)

      read(30,*) nkvec_tmp
      if(nkvec_tmp.ne.nkvec) then
        write(ounit,'(''nkvec_tmp,nkvec='',9i5)') nkvec_tmp,nkvec
        call fatal_error ('nkvec_tmp != nkvec in read_orb_pw_tm')
      endif

      write(3,'(2i6,'' nkvec,ngvec'')') nkvec,ngvec
c     if(nkvec.gt.IVOL_RATIO) call fatal_error ('nkvec>IVOL_RATIO')

      jorb=0
      ngvec_orb=0
      do ikv=1,nkvec
        read(30,*) ikvec,nband(ikv),ng,(rkvec_tmp(k),k=1,3)
c Warning temp
        ngg(ikv)=ng

c       write(3,'(2i4,3f9.5'' ikvec, nband, rkvec(in cartesian units)'')') ikv,nband(ikv),(rkvec(k,ikv),k=1,3)
        write(3,'(2i4,3f9.5'' ikvec, nband, rkvec(in recip. lat. units)'')') ikv,nband(ikv),(rkvec_tmp(k),k=1,3)
        do k=1,3
        rkvec_tmp2(k)=rkvec_tmp(k)
        enddo
        do k=1,3
          rkvec_tmp(k)=0
          do i=1,3
            rkvec_tmp(k)=rkvec_tmp(k)+rkvec_tmp2(i)*glatt(k,i)
          enddo
        enddo
        write(ounit,'(/,''rkvec_tmp in recip. lat. vec. units'',9f9.4)') (rkvec_tmp2(k),k=1,3)
        write(ounit,'(''rkvec_tmp in cartesian coodinates'',9f9.4)') (rkvec_tmp(k),k=1,3)
        do k=1,3
c Warning: tmp commented out
c         if(abs(rkvec_tmp(k)-rkvec(k,ikv)).gt.1.d-5) then
c           write(ounit,'(''kvec='',9f9.6)') (kvec(kk,ikv)+rkvec_shift(kk),kk=1,3)
c           write(ounit,'(''rkvec_tmp,rkvec='',9f9.6)') (rkvec_tmp(kk),kk=1,3),(rkvec(kk,ikv),kk=1,3)
c           call fatal_error ('rkvec_tmp != rkvec in read_orb_pw')
c         endif
        enddo
        if(ng.gt.NGVEC_BIGX) call fatal_error ('ng>NGVEC_BIGX')
        read(30,*) (iwgvec(j),j=1,ng)
        do iband=1,nband(ikv)
          jorb=jorb+1
          read(30,*) iband_tmp,eig
          if(icmplx.ne.0) then
            read(30,*) (c_real(ig),c_imag(ig),ig=1,ng)
           else
            read(30,*) (c_real(ig),ig=1,ng)
          endif
          if(icmplx.eq.0) then
            do ig=1,ng
              c_imag(ig)=0
            enddo
          endif
c If there is only one linearly indep. state formed from psi_k and psi_-k then
c determine if that state is real or imaginary.
          if(k_inv(ikv).eq.1) then
            if(rknorm(ikv).eq.0.d0) then
              ig_min=2
             else
              ig_min=1
            endif
            sum=0
            sum_abs=0
            do ig=ig_min,ng
              sum=sum+c_real(ig)
              sum_abs=sum_abs+abs(c_real(ig))
            enddo
            if(abs(sum/sum_abs).gt.1.d-6) then
              ireal_imag(jorb)=1
             else
              ireal_imag(jorb)=2
            endif
            write(ounit,'(''ikv,iband,ireal_imag,sum,sum_abs='',3i4,9d12.4)')
     &      ikv,iband,ireal_imag(jorb),sum,sum_abs
           else
            ireal_imag(jorb)=0
          endif
          do igv=1,ngvec
            c_rp(igv,jorb)=0
            c_rm(igv,jorb)=0
            c_ip(igv,jorb)=0
            c_im(igv,jorb)=0
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
c               write(ounit,*) (igvec_dft(k,iwgvec(ig)),k=1,3),(igvec(k,igv),k=1,3)
                do k=1,3
                  norm=norm+(igvec_dft(k,iwgvec(ig))-isign*igvec(k,igv))**2
                enddo
                if(norm.eq.0) then
                  ifound=ifound+1
                  ngvec_orb=max(ngvec_orb,igv)
                  c_rp(igv,jorb)=c_rp(igv,jorb)+c_real(ig)
                  c_rm(igv,jorb)=c_rm(igv,jorb)+c_real(ig)*isign
                  c_ip(igv,jorb)=c_ip(igv,jorb)+c_imag(ig)
                  c_im(igv,jorb)=c_im(igv,jorb)+c_imag(ig)*isign
                endif
              enddo
            enddo
          enddo
          if(ifound.ne.ng) then
            write(3,'(''ng,ifound='',2i5)') ng,ifound
            write(ounit,'(''ng,ifound='',2i5)') ng,ifound
            if(ifound.lt.ng) call fatal_error ('probably need to generate more g-vectors')
            call fatal_error ('ifound != ng in read__orb_pw_tm')
          endif
          write(3,'(2i5,f10.6,'' ikvec, iband, eig (Ha)'')') ikv,iband,eig
          if(k_inv(ikv).eq.1) then
            write(3,'(1p2d22.14)') (c_rp(igv,jorb),c_rm(igv,jorb),igv=1,ngvec)
           else
            write(3,'(1p4d22.14)') (c_rp(igv,jorb),c_rm(igv,jorb),c_ip(igv,jorb),c_im(igv,jorb),igv=1,ngvec)
          endif

c Test to see if orbitals and derivs. calculated correctly
          r(1)=.1d0
          r(2)=.2d0
          r(3)=.3d0
          call orbitals_pw2(ikv,iband,jorb,r,orb,dorb,ddorb)
          if(k_inv(ikv).eq.1) then
            if(abs(orb(1,jorb)/orb(1,jorb+1)).gt.1.d0) then
              ireal_imag(jorb)=1
              write(ounit,'(''ireal_imag,orb2_sim='',i1,5d15.8)') ireal_imag(jorb),orb(1,jorb),ddorb(jorb,1),(dorb(jorb,1,k),k=1,3)
             else
              ireal_imag(jorb)=2
              write(ounit,'(''ireal_imag,orb2_sim='',i1,5d15.8)')
     &        ireal_imag(jorb),orb(1,jorb+1),ddorb(jorb+1,1),(dorb(jorb+1,1,k),k=1,3)
            endif
           else
            ireal_imag(jorb)=0
            write(ounit,'(''ireal_imag,orb2_sim='',i1,5d15.8)') ireal_imag(jorb),orb(1,jorb),ddorb(jorb,1),(dorb(jorb,1,k),k=1,3)
            write(ounit,'(''ireal_imag,orb2_sim='',i1,5d15.8)')
     &      ireal_imag(jorb),orb(1,jorb+1),ddorb(jorb+1,1),(dorb(jorb+1,1,k),k=1,3)
          endif

        enddo
      enddo

      write(ounit,'(i3,'' orbitals read in from file orbitals_pw'')') jorb
      if(jorb.lt.norb) then
        write(ounit,'(''jorb,norb='',2i5)') jorb,norb
        call fatal_error ('jorb < norb in read_orb_pw_tm')
      endif

      close(30)

      nsum=0
      do i=1,ngnorm
        nsum=nsum+igmult(i)
        if(nsum.ge.ngvec_orb) then
          ngnorm_orb=i
          goto 70
        endif
      enddo
   70 write(ounit,'(''ngnorm_orb,ngvec_orb='',9i6)') ngnorm_orb,ngvec_orb

c Test to see if orbitals and derivs. calculated correctly by comparing to above
      nelec_sav=nelec
      nelec=1
      call orbitals_pw(r,orb,dorb,ddorb)
      write(ounit,'(''orb2_com='',5d15.8)') (orb(1,ib),ddorb(ib,1),(dorb(ib,1,k),k=1,3),ib=1,4)
      nelec=nelec_sav

c Write file for subsequent read-in, though at present I am not using it.
c Comment it out since it causes problems with mpi.
c     open(4,file='orbitals_pw')
c     write(4,'(i1,3i5,'' icmplx,nkvec,ngnorm_orb,ngvec_orb'')') icmplx,nkvec,ngnorm_orb,ngvec_orb
c     rewind 3
c     read(3,*) nkvec_tmp,ngvec_tmp
c     do 80 ikv=1,nkvec
c       read(3,'(2i4,3f9.5)') ikv_tmp,nband(ikv),(rkvec_tmp(k),k=1,3)
c       write(4,'(2i4,3f9.5'' ikvec, nband, rkvec(in recip. lat. units)'')') ikv,nband(ikv),(rkvec_tmp(k),k=1,3)
c       do 80 iband=1,nband(ikv)
c         read(3,'(2i5,f10.6)') ikv_tmp,iband_tmp,eig
c         write(4,'(2i5,f10.6,'' ikvec, iband, eig (Ha)'')') ikv,iband,eig
c         if(k_inv(ikv).eq.1) then
c           read(3,'(1p2d22.14)') (c_rp(igv,jorb),c_rm(igv,jorb),igv=1,ngvec)
c           write(4,'(1p2d22.14)') (c_rp(igv,jorb),c_rm(igv,jorb),igv=1,ngvec_orb)
c          else
c           read(3,'(1p4d22.14)') (c_rp(igv,jorb),c_rm(igv,jorb),c_ip(igv,jorb),c_im(igv,jorb),igv=1,ngvec)
c           write(4,'(1p4d22.14)') (c_rp(igv,jorb),c_rm(igv,jorb),c_ip(igv,jorb),c_im(igv,jorb),igv=1,ngvec_orb)
c         endif
c  80     continue

c     close(3,status='delete')
      close(4)
      return
      end
c-----------------------------------------------------------------------
      subroutine orbitals_pw2(ikv,iband,jorb,x,orb,dorb,ddorb)
c Written by Cyrus Umrigar
c Calculate pw orbitals.
c isortg is used to map g-vectors from iv to ig.
c At present it is assumed that k-vectors are in the correct order, but
c if not one could use isortk to map iorb.
c This is the straightforward evaluation for checking purposes only.

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use ewald_mod, only: IVOL_RATIO
      use ewald_mod, only: NGVECX
      use ewald_mod, only: NGVEC_BIGX
      use periodic, only: glatt
      use periodic, only: gvec
      use periodic, only: ireal_imag, k_inv, nband
      use periodic, only: nkvec
      use periodic, only: rkvec
      use tempor_test, only: c_imag, c_real, igvec_dft, iwgvec, ngg, ngvec_dft
      use contrl_file, only: ounit
      use coefs, only: norb
      use precision_kinds, only: dp
      use vmc_mod, only: norb_tot
      use control, only: ipr
      use system, only: nelec
      implicit none

      integer :: i, iband, ig, ig2, ikv
      integer :: iorb, jorb, k, nelecc
      real(dp) :: cos, dot, sin
      real(dp), dimension(3) :: x
      real(dp), dimension(nelec,*) :: orb
      real(dp), dimension(norb_tot,nelec,3) :: dorb
      real(dp), dimension(norb_tot,nelec) :: ddorb
      real(dp), dimension(3,NGVEC_BIGX) :: gvec_dft
      real(dp), dimension(NGVEC_BIGX) :: gnorm_dft








c     dimension dcos_rp(3),dsin_rm(3),dcos_ip(3),dsin_im(3)
c    &,cos_g(nelec,NGVECX),sin_g(nelec,NGVECX),dcos_g(3,nelec,NGVECX),dsin_g(3,nelec,NGVECX)
c    &,ddcos_g(nelec,NGVECX),ddsin_g(nelec,NGVECX)
c    &,cos_k(nelec,IVOL_RATIO),sin_k(nelec,IVOL_RATIO),dcos_k(3,nelec,IVOL_RATIO),dsin_k(3,nelec,IVOL_RATIO)
c    &,ddcos_k(nelec,IVOL_RATIO),ddsin_k(nelec,IVOL_RATIO)

      write(ounit,'(''nelec,norb,nkvec in orbitals_pw'',9i5)') nelec,norb,nkvec

      do ig=1,ngvec_dft
        do k=1,3
          gvec_dft(k,ig)=0
          do i=1,3
            gvec_dft(k,ig)=gvec_dft(k,ig)+igvec_dft(i,ig)*glatt(k,i)
          enddo
        enddo
c       write(ounit,'(/,''igvec_dft in recip. lat. vec. units'',9i4)') (igvec_dft(i,ig),i=1,3)
c  26   write(ounit,'(''gvec_dft in cartesian coodinates'',9f9.4)') (gvec_dft(k,ig),k=1,3)
      enddo

c Warning: c_real and c_imag are only stored for particular k-vector and band
c Warning: for the moment just use this for testing 1 electron
      nelecc=1
      do i=1,nelecc
c     iorb=0
c     do 80 ikv=1,nkvec
c       do 80 iband=1,nband(ikv)

        iorb=jorb
c       if(ireal_imag(iorb).eq.0 .or. ireal_imag(iorb).eq.1) then

c       iorb=iorb+1
        orb(i,iorb)=0
        ddorb(iorb,i)=0
        do k=1,3
          dorb(iorb,i,k)=0
        enddo
        do ig=1,ngg(ikv)
          ig2=iwgvec(ig)
          dot=0
          gnorm_dft(ig2)=0
          do k=1,3
            gnorm_dft(ig2)=gnorm_dft(ig2)+(rkvec(k,ikv)+gvec_dft(k,ig2))**2
            dot=dot+(rkvec(k,ikv)+gvec_dft(k,ig2))*x(k)
          enddo
          orb(i,iorb)=orb(i,iorb)+c_real(ig)*cos(dot)-c_imag(ig)*sin(dot)
          if(ipr.ge.4 .and. ig.le.22) write(ounit,'(''rkvec+gvec'',2i4,7f9.4,f18.12)')
     &    ig,ig2,(rkvec(k,ikv)+gvec_dft(k,ig2),k=1,3),c_real(ig),dot,cos(dot),sin(dot),orb(i,iorb)
          do k=1,3
            dorb(iorb,i,k)=dorb(iorb,i,k)+(rkvec(k,ikv)+gvec_dft(k,ig2))*(-c_real(ig)*sin(dot)-c_imag(ig)*cos(dot))
          enddo
          ddorb(iorb,i)=ddorb(iorb,i)-gnorm_dft(ig2)*(c_real(ig)*cos(dot)-c_imag(ig)*sin(dot))
        enddo

        write(ounit,'(''ikv,iband,nband(ikv),k_inv(ikv)'',9i5)') ikv,iband,nband(ikv),k_inv(ikv)
        write(ounit,'(''real orb='',i5,9d12.4)') iorb,orb(i,iorb)

c       if(k_inv(ikv).eq.1) goto 80
        iorb=iorb+1
c       endif

c       if(iorb.lt.norb) then
c       if(ireal_imag(iorb).eq.0 .or. ireal_imag(iorb).eq.2) then

        orb(i,iorb)=0
        ddorb(iorb,i)=0
        do k=1,3
          dorb(iorb,i,k)=0
        enddo
        do ig=1,ngg(ikv)
          ig2=iwgvec(ig)
          dot=0
          gnorm_dft(ig2)=0
          do k=1,3
            gnorm_dft(ig2)=gnorm_dft(ig2)+(rkvec(k,ikv)+gvec_dft(k,ig2))**2
            dot=dot+(rkvec(k,ikv)+gvec_dft(k,ig2))*x(k)
          enddo
          orb(i,iorb)=orb(i,iorb)+c_real(ig)*sin(dot)+c_imag(ig)*cos(dot)
c         if(ipr.ge.4 .and. ig.le.22) write(ounit,'(''rkvec+gvec'',2i4,7f9.4,f18.12)')
c    &    ig,ig2,(rkvec(k,ikv)+gvec_dft(k,ig2),k=1,3),c_real(ig),dot,cos(dot),sin(dot),orb(i,iorb)
          if(ipr.ge.4 .and. ig.le.22) write(ounit,'(''rkvec,gvec'',8f9.4)') (rkvec(k,ikv),gvec_dft(k,ig2),k=1,3)
          if(ipr.ge.4 .and. ig.le.22) write(ounit,'(''rkvec+gvec'',2i4,8f9.4,f18.12)')
     &    ig,ig2,(rkvec(k,ikv)+gvec_dft(k,ig2),k=1,3),c_real(ig),x(1),dot,cos(dot),sin(dot),orb(i,iorb)
          do k=1,3
            dorb(iorb,i,k)=dorb(iorb,i,k)+(rkvec(k,ikv)+gvec_dft(k,ig2))*(c_real(ig)*cos(dot)-c_imag(ig)*sin(dot))
          enddo
   70     ddorb(iorb,i)=ddorb(iorb,i)-gnorm_dft(ig2)*(c_real(ig)*sin(dot)+c_imag(ig)*cos(dot))
        enddo

        write(ounit,'(''ikv,iband,nband(ikv),k_inv(ikv)'',9i5)') ikv,iband,nband(ikv),k_inv(ikv)
        write(ounit,'(''imag orb='',i5,9d12.4)') iorb,orb(i,iorb)

c       endif
c       endif

   80   continue
      enddo

      return
      end
      end module
