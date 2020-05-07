! MELEC  >= number of electrons
! MORB   >= number of orbitals
! MBASIS >= number of basis functions
! MDET   >= number of determinants
! MCENT  >= number of centers
! MCTYPE >= number of center types
! MCTYP3X=max(3,MCTYPE)

! Slater matrices are dimensioned (MELEC/2)**2 assuming 
! equal numbers of up and down spins. MELEC has to be 
! correspondingly larger if spin polarized calculations 
! are attempted.

      integer MELEC,MORB,MBASIS,MDET,MCENT,MCTYPE,MCTYP3X
      integer NSPLIN,nrad,MORDJ,MORDJ1,MMAT_DIM,MMAT_DIM2

      real*8 radmax,delri

      character*20 method

      parameter(MELEC=50,MORB=550,MBASIS=550,MDET=15000,MCENT=20,MCTYPE=3)
      parameter(MCTYP3X=5,NSPLIN=1001,MORDJ=7,radmax=10.d0,nrad=3001)
      parameter(MMAT_DIM=(MELEC*MELEC)/4,MMAT_DIM2=(MELEC*(MELEC-1))/2)
      parameter(MORDJ1=MORDJ+1,delri=(nrad-1)/radmax)

      common /method_opt/ method
