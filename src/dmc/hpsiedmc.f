      subroutine psiedmc(iel,iw,coord,psid,psij,iflag)
c Written by Claudia Filippi

      use const, only: delta, deltai, etrial, fbias, hb, imetro, ipr, nelec, pi
      use forcest, only: fgcm2, fgcum
      use forcepar, only: deltot, istrech, nforce
      use config, only: d2o, peo_dmc, psido_dmc, psijo_dmc, vold_dmc, xold_dmc

      implicit real*8(a-h,o-z)




      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'
      common /force_dmc/ itausec,nwprod


      dimension coord(3),x(3,MELEC)

      do 10 ic=1,3
      do 10 i=1,iel-1
  10    x(ic,i)=xold_dmc(ic,i,iw,1)

      do 20 ic=1,3
  20    x(ic,iel)=coord(ic)

      do 30 ic=1,3
      do 30 i=iel+1,nelec
  30    x(ic,i)=xold_dmc(ic,i,iw,1)

      idum=1
      call psie(iel,x,psid,psij,idum,iflag)

      return
      end
