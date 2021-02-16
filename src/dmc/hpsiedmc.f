      subroutine psiedmc(iel,iw,coord,psid,psij,iflag)
c Written by Claudia Filippi

      use const, only: delta, deltai, etrial, fbias, hb, imetro, ipr, nelec, pi
      use forcest, only: fgcm2, fgcum
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'
      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /force_dmc/ itausec,nwprod

      common /config/ xold(3,MELEC,MWALK,MFORCE),vold(3,MELEC,MWALK,MFORCE),
     &psido(MWALK,MFORCE),psijo(MWALK,MFORCE),peo(MWALK,MFORCE),d2o(MWALK,MFORCE)

      dimension coord(3),x(3,MELEC)

      do 10 ic=1,3
      do 10 i=1,iel-1
  10    x(ic,i)=xold(ic,i,iw,1)

      do 20 ic=1,3
  20    x(ic,iel)=coord(ic)

      do 30 ic=1,3
      do 30 i=iel+1,nelec
  30    x(ic,i)=xold(ic,i,iw,1)

      idum=1
      call psie(iel,x,psid,psij,idum,iflag)

      return
      end
