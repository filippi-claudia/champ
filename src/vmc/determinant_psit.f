      subroutine determinant_psit(determ,istate)

      use force, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use dets, only: cdet, ndet
      use wfsec, only: iwf

      implicit real*8(a-h,o-z)


      common /multislater/ detiab(MDET,2)

      determ=0
      do 110 k=1,ndet
  110   determ=determ+detiab(k,1)*detiab(k,2)*cdet(k,istate,iwf)


      return
      end
