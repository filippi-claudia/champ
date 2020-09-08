c-----------------------------------------------------------------------
      subroutine determinante_psit(iel,determ,istate)

      use force, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use dets, only: cdet, ndet
      use elec, only: nup
      use wfsec, only: iwf
      implicit real*8(a-h,o-z)


      common /multislater/ detiab(MDET,2)
      common /multislatern/ detn(MDET)
     &,orbn(MORB),dorbn(3,MORB),ddorbn(MORB)


      determ=0

      do 110 k=1,ndet
        if(iel.le.nup) then
          det=detn(k)*detiab(k,2)
         else
          det=detiab(k,1)*detn(k)
        endif

        determ=determ+det*cdet(k,istate,iwf)
  110 continue

      return
      end
