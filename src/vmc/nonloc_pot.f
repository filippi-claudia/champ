      subroutine nonloc_pot(x,rshift,rvec_en,r_en,pe,vpsp_det,dvpsp_dj,t_vpsp,i_vpsp,ifr)
c Written by Claudia Filippi; modified by Cyrus Umrigar
c Calculates the local and nonlocal components of the pseudopotential
c Calculates non-local potential derivatives
c pe_en(loc) is computed in distances and pe_en(nonloc) here in nonloc_pot if nloc !=0 and iperiodic!=0.
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'pseudo.h'

      common /contrl_per/ iperiodic,ibasis
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,lpot(MCTYPE),nloc
      common /dets/ cdet(MDET,MSTATES,MWF),ndet

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      dimension x(3,*),rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
     &,vpsp_det(*),dvpsp_dj(*),t_vpsp(MCENT,MPS_QUAD,*)

      if(i_vpsp.gt.0)then
        i1=i_vpsp
        i2=i_vpsp
       else
        i1=1
        i2=nelec
      endif
      do 20 i=i1,i2
        if(nloc.eq.1) then
          call getvps(r_en,i)
         elseif(nloc.eq.2.or.nloc.eq.3) then
          call getvps_tm(r_en,i)
         elseif(nloc.eq.4) then
          call getvps_gauss(rvec_en,r_en,i)
         elseif(nloc.eq.5) then
          call getvps_champ(r_en,i)
        endif
   20 continue

c local component (highest angular momentum)
      if(iperiodic.eq.0) then
        do 30 ic=1,ncent
          do 30 i=i1,i2
   30       pe=pe+vps(i,ic,lpot(iwctype(ic)))
      endif

c non-local component (division by the Jastrow already in nonloc)
      call nonloc(x,rshift,rvec_en,r_en,vpsp_det,dvpsp_dj,t_vpsp,i_vpsp)

      return
      end
