      subroutine readps
c Written by Claudia Filippi
      use pseudo_mod, only: MPS_GRID
      use vmc_mod, only: nrad
      use atom, only: znuc, iwctype, nctype

      use pseudo_fahy, only: drad, dradl, nlrad, npotl, potl, ptnlc, rcmax

      use pseudo, only: lpot

      use qua, only: nquad, wq, xq, xq0, yq, yq0, zq, zq0

      implicit real*8(a-h,o-z)


      character*20 filename,atomtyp

c nquad = number of quadrature points
c nlang = number of non-local potentials
c rcmax = cutoff radius for non-local potential
c npotl = number of mesh point for local potential
c dradl = spacing of uniform mesh for local potential

      do 20 ic=1,nctype

      if(ic.lt.10) then
        write(atomtyp,'(i1)') ic
       elseif(ic.lt.100) then
        write(atomtyp,'(i2)') ic
       else
        call fatal_error('READPS: problem atomtyp')
      endif

      filename='ps.data.'//atomtyp(1:index(atomtyp,' ')-1)
      open(3,file=filename,status='old',form='formatted')

      read(3,*) nquad
      write(6,'(''quadrature points'',i4)') nquad

      read(3,*) nlang,rcmax(ic)
      lpot(ic)=nlang+1

c local potential
      read(3,*)
      read(3,*) npotl(ic),nzion,dradl(ic)
      if(npotl(ic).gt.MPS_GRID) call fatal_error('READPS: npotl gt MPS_GRID')
      if(nzion.ne.int(znuc(iwctype(ic)))) call fatal_error('READPS: nzion ne znuc')

      read(3,*) (potl(i,ic),i=1,npotl(ic))

      do 5 i=1,npotl(ic)
  5     write(33,*) (i-1)*dradl(ic),potl(i,ic)

c non-local potential
      read(3,*)
      read(3,*) nlrad(ic),drad(ic)
      if(nlrad(ic).gt.MPS_GRID) call fatal_error('READPS: nrad gt MPS_GRID')

      if(drad(ic)*(nlrad(ic)-1).le.rcmax(ic)) then
        write(6,'(''non-local table max radius = '',
     &  f10.5,'' too small for cut-off = '',f10.5)')
     &  drad(ic)*(nlrad(ic)-1),rcmax(ic)
        call fatal_error('READPS')
      endif

      do 10 l=1,nlang
        read(3,*)
        read(3,*) (ptnlc(i,ic,l),i=1,nlrad(ic))
      do 10 i=1,nlrad(ic)
 10     write(34,*) (i-1)*drad(ic),ptnlc(i,ic,l)

      close(3)
 20   continue


      call gesqua (nquad,xq0,yq0,zq0,wq)
c     call gesqua (nquad,xq,yq,zq,wq)

      write(6,'(''quadrature points'')')
      do 30 i=1,nquad
 30     write(6,'(''xyz,w'',4f10.5)') xq0(i),yq0(i),zq0(i),wq(i)

      return
      end

c-----------------------------------------------------------------------
      subroutine getvps(rad,iel)
c Written by Claudia Filippi
      use vmc_mod, only: MELEC, MCENT
      use atom, only: znuc, iwctype, ncent, ncent_tot
      use const, only: nelec
      use pseudo_fahy, only: drad, dradl, npotl, potl, ptnlc, rcmax

      use pseudo, only: lpot, vps

      implicit real*8(a-h,o-z)





      dimension rad(nelec,ncent_tot)

      do 10 ic=1,ncent
        ict=iwctype(ic)
        r=rad(iel,ic)
c local potential
        if(r.lt.(npotl(ict)-1)*dradl(ict)) then
          ri=r/dradl(ict)
          ir=int(ri)
          ri=ri-dfloat(ir)
          ir=ir+1
          vps(iel,ic,lpot(ict))=potl(ir+1,ict)*ri+(1.d0-ri)*potl(ir,ict)
         else
          vps(iel,ic,lpot(ict))=-znuc(ict)/r
        endif
c non-local pseudopotential
        do 10 l=1,lpot(ict)-1
          if(r.lt.rcmax(ict)) then
            ri=r/drad(ict)
            ir=int(ri)
            ri=ri-dfloat(ir)
            ir=ir+1
            vps(iel,ic,l)=ptnlc(ir+1,ict,l)*ri+(1.d0-ri)*ptnlc(ir,ict,l)
           else
            vps(iel,ic,l)=0.0d0
          endif
   10 continue

      return
      end
