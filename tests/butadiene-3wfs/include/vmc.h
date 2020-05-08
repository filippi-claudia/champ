c MELEC  >= number of electrons
c MORB   >= number of orbitals
c MBASIS >= number of basis functions
c MDET   >= number of determinants
c MCENT  >= number of centers
c MCTYPE >= number of center types
c MCTYP3X=max(3,MCTYPE)

c Slater matrices are dimensioned (MELEC/2)**2 assuming 
c equal numbers of up and down spins. MELEC has to be 
c correspondingly larger if spin polarized calculations 
c are attempted.

      integer MELEC,MORB,MBASIS,MDET,MCENT,MCTYPE,MCTYP3X,
     &NSPLIN,nrad,MORDJ,MORDJ1,MMAT_DIM,MMAT_DIM2

      real*8 radmax,delri

      character*20 method

      parameter(MELEC=50,MORB=500,MBASIS=500,MDET=40000,MCENT=10,MCTYPE=5,
     &MCTYP3X=5,NSPLIN=1001,MORDJ=7,radmax=10.d0,nrad=3001,
     &MMAT_DIM=(MELEC*MELEC)/4,MMAT_DIM2=(MELEC*(MELEC-1))/2,
     &MORDJ1=MORDJ+1,delri=(nrad-1)/radmax)
! PLT@eScienceCenter(2020) Moved the parameter here:
! "For Jastrow4 NEQSX=2*(MORDJ-1) is sufficient.
!  For Jastrow3 NEQSX=2*MORDJ should be sufficient.
!  I am setting NEQSX=6*MORDJ simply because that is how it was for
!  Jastrow3 for reasons I do not understand."
!     parameter(NEQSX=2*(MORDJ-1),MTERMS=55)
      parameter(NEQSX=6*MORDJ,MTERMS=55)

      common /method_opt/ method
