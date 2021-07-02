      subroutine determinant_psit(determ,istate)

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use dets, only: cdet, ndet
      use wfsec, only: iwf

      use multislater, only: detiab
      implicit real*8(a-h,o-z)


      determ=0.0d0
      do 110 k=1,ndet
  110   determ=determ+detiab(k,1)*detiab(k,2)*cdet(k,istate,iwf)


      return
      end
