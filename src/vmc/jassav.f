      subroutine jassav(iel,iflag,istate)
c     Written by Claudia Filippi
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use const, only: nelec
      use jaso, only: d2ijo, d2o, fijo, fjo, fso, fsumo
      use velocity_jastrow, only: vj, vjn
      use jasn, only: d2ijn, d2n, fijn, fjn, fsn, fsumn

      implicit real*8(a-h,o-z)

      fsumo(istate)=fsumn(istate)
      do i=1,nelec
         fjo(1,i,istate)=fjn(1,i,istate)
         fjo(2,i,istate)=fjn(2,i,istate)
         fjo(3,i,istate)=fjn(3,i,istate)
      enddo

      do j=1,iel
         fso(iel,j,istate)=fsn(iel,j,istate)
      enddo

      do j=iel+1,nelec
         fso(j,iel,istate)=fsn(j,iel,istate)
      enddo

      do j=1,nelec
         fijo(1,iel,j,istate)=fijn(1,iel,j,istate)
         fijo(2,iel,j,istate)=fijn(2,iel,j,istate)
         fijo(3,iel,j,istate)=fijn(3,iel,j,istate)
         fijo(1,j,iel,istate)=fijn(1,j,iel,istate)
         fijo(2,j,iel,istate)=fijn(2,j,iel,istate)
         fijo(3,j,iel,istate)=fijn(3,j,iel,istate)
      enddo

      if(iflag.gt.0) then
         d2o=d2n
         do j=1,iel
            d2ijo(iel,j,istate)=d2ijn(iel,j,istate)
         enddo

         do j=iel+1,nelec
            d2ijo(j,iel,istate)=d2ijn(j,iel,istate)
         enddo
      endif

      do i=1,nelec
         vj(1,i,istate)=vjn(1,i,istate)
         vj(2,i,istate)=vjn(2,i,istate)
         vj(3,i,istate)=vjn(3,i,istate)
      enddo

      end subroutine
