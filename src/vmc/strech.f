      subroutine strech(x,xstrech,ajacob,ifr,istrech_el)
c Written by Cyrus Umrigar and Claudia Filippi
c Modified by A. Amovilli for forces in PCM
c Uses the coordinate transform described in:
c 1) Two Aspects of Quantum Monte Carlo: Determination of Accurate Wavefunctions and
c    Determination of Potential Energy Surfaces of Molecules, C.J. Umrigar,
c    Int. J. Quant. Chem. Symp., 23, 217 (1989).
c 2) Correlated sampling in quantum Monte Carlo: A route to forces,
c    Claudia Filippi and C. J. Umrigar, Phys. Rev. B., 61, R16291, (2000).

c stretch space so that electrons close to a nucleus move almost
c rigidly with that nucleus
      use precision_kinds, only: dp
      use pcm, only: MCHS, MCHV
      use force_mod, only: MFORCE, MFORCE_WT_PRD
      use forcepar, only: istrech, alfstr
      use vmc_mod, only: MELEC, MCENT
      use atom, only: znuc, cent, pecent, iwctype, ncent, ncent_tot
      use const, only: nelec
      use force_dmc, only: itausec, nwprod
      use forcepar, only: deltot, istrech, nforce
      use forcestr, only: delc
      use pcm_force, only: sch_s
      use wfsec, only: iwftype
      use contr3, only: mode
      use pcm_cntrl, only: ipcm
      use pcm_parms, only: ch, nch, nchs
      use pcm_parms, only: nesph
      use pcm_parms, only: xpol
      use pcm_ameta, only: eta
      use pcm_pot, only: penups, penupv
      use pcm_inda, only: inda
      use optwf_contrl, only: ioptwf

      implicit real*8(a-h,o-z)

      parameter (zero=0.d0,one=1.d0)

      real(dp), ALLOCATABLE, save :: centsav(:,:)
      real(dp), ALLOCATABLE, save :: pecentn(:)
      real(dp), ALLOCATABLE, save :: xpolsav(:,:)

      dimension x(3,nelec),xstrech(3,nelec)
      dimension wt(ncent_tot),dvol(3,3),dwt(3,ncent_tot),dwtsm(3)
      dimension cent_str(3,ncent_tot)
      dimension q_strech(MCHS),efsol(MCHS),wt_pcm(ncent_tot)

      if(.not.allocated(centsav)) allocate(centsav(3, ncent_tot))
      if(.not.allocated(pecentn)) allocate(pecentn(MFORCE))
      if(.not.allocated(xpolsav)) allocate(xpolsav(3,MCHV))


c set center and n-n potential for secondary geometries
      pecent=pecentn(ifr)
      do 1 icent=1,ncent
        do 1 k=1,3
    1     cent(k,icent)=centsav(k,icent)+delc(k,icent,ifr)

      ajacob=one

c PCM
      if(ipcm.eq.3.and.ioptwf.eq.0) then

c positions of surface charges on spheres rigidly displaced
          do 3 j=1,nchs
            is=inda(j)
            do 3 k=1,3
   3         xpol(k,j)=xpolsav(k,j)+delc(k,is,ifr)
c positions of volume charges space warped
          do 6 j=nchs+1,nch
            wtsm=zero
            do 5 icent=1,ncent
              dist2=zero
              do 4 k=1,3
    4           dist2=dist2+(xpolsav(k,j)-centsav(k,icent))**2
                dist=dsqrt(dist2)
                if(istrech.eq.1) wt_pcm(icent)=dexp(-alfstr*dist)
                if(istrech.eq.2) wt_pcm(icent)=one/dist**alfstr
                if(istrech.eq.3) wt_pcm(icent)=dexp(alfstr/dist)
                wtsm=wtsm+wt_pcm(icent)
    5       continue
            wtsmi=one/wtsm
            do 6 icent=1,ncent
              wt_pcm(icent)=wt_pcm(icent)*wtsmi
              do 6 k=1,3
    6           xpol(k,j)=xpolsav(k,j)+wt_pcm(icent)*delc(k,icent,ifr)

          do 7 i=1,nchs
    7        ch(i)=sch_s(i,ifr)

c endif PCM
      endif

      if(istrech_el.eq.0) return

      do 8 i=1,nelec
        do 8 k=1,3
    8     xstrech(k,i)=x(k,i)

      if(istrech.eq.0) return

      do 50 i=1,nelec

        wtsm=zero
c initialize volume change matrix
        do 10 k=1,3
          dwtsm(k)=zero
          do 10 j=1,3
            dvol(j,k)=zero
            if(j.eq.k) dvol(j,k)=one
   10   continue

        do 20 icent=1,ncent
          dist2=zero
          do 15 k=1,3
   15       dist2=dist2+(x(k,i)-centsav(k,icent))**2
            dist=dsqrt(dist2)
            if(istrech.eq.1) wt(icent)=dexp(-alfstr*dist)
            if(istrech.eq.2) wt(icent)=one/dist**alfstr
            if(istrech.eq.3) wt(icent)=dexp(alfstr/dist)
            wtsm=wtsm+wt(icent)
            do 20 k=1,3
              if(istrech.eq.1) dwt(k,icent)=-alfstr*dexp(-alfstr*dist)*
     &        (x(k,i)-centsav(k,icent))/dist
              if(istrech.eq.2) dwt(k,icent)=-alfstr*
     &        (x(k,i)-centsav(k,icent))/dist**(alfstr+2)
              if(istrech.eq.3) dwt(k,icent)=-alfstr*dexp(alfstr/dist)*
     &        (x(k,i)-centsav(k,icent))/(dist2*dist)
              dwtsm(k)=dwtsm(k)+dwt(k,icent)
   20     continue
          wtsmi=one/wtsm

          do 40 icent=1,ncent
            do 40 k=1,3
              dwt(k,icent)=(wtsm*dwt(k,icent)-wt(icent)*dwtsm(k))/
     &        wtsm**2
              dvol(1,k)=dvol(1,k)+dwt(k,icent)*delc(1,icent,ifr)
              dvol(2,k)=dvol(2,k)+dwt(k,icent)*delc(2,icent,ifr)
   40         dvol(3,k)=dvol(3,k)+dwt(k,icent)*delc(3,icent,ifr)
          call matinv(dvol,3,det)
          ajacob=ajacob*det
          do 50 icent=1,ncent
            wt(icent)=wt(icent)*wtsmi
            do 50 k=1,3
              xstrech(k,i)=xstrech(k,i)+wt(icent)*delc(k,icent,ifr)
c end loop over electrons
   50 continue

      return

c Set up n-n potential energy (and PCM related quantities) at displaced positions 
      entry setup_force

      write(6,'(''istrech,alfstr ='',i4,2f10.5)') istrech,alfstr

      do 60 i=1,nforce
        write(6,*) '--------------'
        do 60 ic=1,ncent
   60     write(6,'(''center '',i2,'' conf '',i2,'' displace '',
     &    3f15.7)') ic,i,(delc(k,ic,i),k=1,3)
      write(6,'(''iwftypes'',20i2)') (iwftype(i),i=1,nforce)

      if(index(mode,'dmc').ne.0) then
        if(nwprod.gt.MFORCE_WT_PRD) call fatal_error('STRETCH: nwprod gt MFORCE_WT_PRD')
        write(6,'(''nwprod,itausec='',2i4)') nwprod,itausec
      endif

      do 65 icent=1,ncent
        do 65 k=1,3
   65     centsav(k,icent)=cent(k,icent)

c' PCM
      if(ipcm.eq.3) then
        do 66 j=1,nch
          do 66 k=1,3
   66       xpolsav(k,j)=xpol(k,j)

c Interatomic forces (and not wave function optimization) 
        if(ioptwf.eq.0) then

          open(54,file='field',status='old',form='formatted')
          rewind(54)
          do 67 i=1,nchs
            do 67 j=1,nesph
   67        read(54,*)
          do 68 i=1,nchs
   68       read(54,*) xi,yi,zi,efsol(i)
          close(54)
          do 71 k=1,nchs
            enk=0.0d0
            do 69 l=1,ncent
              xx=xpolsav(1,k)-cent(1,l)
              yy=xpolsav(2,k)-cent(2,l)
              zz=xpolsav(3,k)-cent(3,l)
              rr2=xx**2+yy**2+zz**2
              rr3=rr2**1.5d0
              cc1=xx*eta(1,k)
              cc2=yy*eta(2,k)
              cc3=zz*eta(3,k)
              cc=cc1+cc2+cc3
   69         enk=enk+znuc(iwctype(l))*cc/rr3
            env=0.0d0
c           do 70 l=nchs+1,nch
c             xx=xpolsav(1,k)-xpolsav(1,l)
c             yy=xpolsav(2,k)-xpolsav(2,l)
c             zz=xpolsav(3,k)-xpolsav(3,l)
c             rr2=xx**2+yy**2+zz**2
c             rr3=rr2**1.5d0
c             cc1=xx*eta(1,k)
c             cc2=yy*eta(2,k)
c             cc3=zz*eta(3,k)
c             cc=cc1+cc2+cc3
c  70         env=env+ch(l)*cc/rr3
   71       efsol(k)=efsol(k)+enk+env
c endif interatomic forces are being computed
        endif
c endif PCM
      endif

c loop over geometries (if wf optimization, geometries for different adiag are equal)
      do 200 ifl=1,nforce

        do 80 i=1,ncent
            do 81 k=1,3
   81         cent_str(k,i)=cent(k,i)+delc(k,i,ifl)
   80     call pot_nn(cent_str,znuc,iwctype,ncent,pecentn(ifl))

c PCM
c for wave function optimization or energy calculation, positions/charges unchanged
         if(ipcm.eq.3.and.ioptwf.eq.0) then

            do 82 j=1,nchs
              is=inda(j)
              do 82 k=1,3
   82          xpol(k,j)=xpolsav(k,j)+delc(k,is,ifl)

            call sigma_R(efsol,q_strech)

            do 83 i=1,nchs
   83         sch_s(i,ifl)=2*q_strech(i)

            delta_qs=0.d0
            do 84 i=1,nchs
   84         delta_qs=delta_qs+dabs(q_strech(i)-ch(i))
c check deviation of surface charges from charges of primary geometry
            write (6,'(''Geometry'',i4,'' : Deviation of surface charges from primary charges'',1p1d14.5)') ifl,delta_qs
c printout charges if deviation is  big
            if (delta_qs.gt.1.d-3) then
              do 85 i=1,nchs
   85           write(6,'(''Warning: Large deviation in surface charges'',2d8.4)') q_strech(i),ch(i)
            endif

            do 89 j=nchs+1,nch
              wtsm=zero
              do 87 icent=1,ncent
                dist2=zero
                do 86 k=1,3
   86             dist2=dist2+(xpolsav(k,j)-centsav(k,icent))**2
                  dist=dsqrt(dist2)
                  if(istrech.eq.1) wt_pcm(icent)=dexp(-alfstr*dist)
                  if(istrech.eq.2) wt_pcm(icent)=one/dist**alfstr
                  if(istrech.eq.3) wt_pcm(icent)=dexp(alfstr/dist)
                  wtsm=wtsm+wt_pcm(icent)
   87         continue
              wtsmi=one/wtsm
              do 89 icent=1,ncent
                wt_pcm(icent)=wt_pcm(icent)*wtsmi
                do 89 k=1,3
   89             xpol(k,j)=xpolsav(k,j)+wt_pcm(icent)*delc(k,icent,ifl)

              penups_fc=0.d0
              penupv_fc=0.d0
                do 96 i=1,ncent
                  do 92 j=1,nchs
                    js=inda(j)
                    rnp2=0.d0
                    do 90 k=1,3
   90                 rnp2=rnp2+(xpol(k,j)-cent_str(k,i))**2.0d0
                    rnp=dsqrt(rnp2)
                    penups_fc=penups_fc+0.5d0*znuc(iwctype(i))*sch_s(j,ifl)/rnp
   92           continue

                do 96 j=nchs+1,nch
                  rnp2=0.d0
                  do 94 k=1,3
   94               rnp2=rnp2+(xpol(k,j)-cent_str(k,i))**2.0d0
                  rnp=dsqrt(rnp2)
                  penupv_fc=penupv_fc+0.5d0*znuc(iwctype(i))*ch(j)/rnp
   96         continue
              delta_gpol_fc=penups_fc-penups+penupv_fc-penupv
              write(6,'(''nuclear delta_gpol contribution to force'',i4,'' ='',1p1d14.5)') ifl,delta_gpol_fc
              pecentn(ifl)=pecentn(ifl)+delta_gpol_fc
c endif PCM
          endif
          
c end loop forces
  200 continue

      write(6,'(''n-n potential energies '',10f10.5)') (pecentn(ifl),ifl=1,nforce)

      do 300 ifl=1,nforce
        deltot(ifl)=zero
        rsq=zero
        do 350 jc=1,ncent
          do 350 k=1,3
            rcm=zero
            do 150 ic=1,ncent
              rcm=rcm+delc(k,ic,ifl)
  150         rsq=rsq+
     &        (cent(k,ic)+delc(k,ic,ifl)-cent(k,jc)-delc(k,jc,ifl))**2
            rcm=rcm/ncent
  350       deltot(ifl)=deltot(ifl)+(delc(k,jc,ifl)-rcm)**2
        if(ifl.eq.1) rsq1=rsq
c Warning: TEMPORARY: multiplication by ncent right for diatomics
c        deltot(ifl)=sign(dsqrt(deltot(ifl)*ncent),rsq-rsq1)
        deltot(ifl)=1.d0
  300   if(deltot(ifl).eq.0) deltot(ifl)=1.d0


      write(6,'(''deltot '',10f10.5)') (deltot(ifl),ifl=1,nforce)

      return
      end
