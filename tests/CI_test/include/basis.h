!   /basis/
!        ncent  = number of centers
!        zex    = screening constants for each basis function
!        cent   = positions of each center
!        pecent = potential energy of the centers
!        znuc   = charge of the nuclei (centers)
!        n1s    = number of 1s functions at each center
!        n2s    = number of 2s functions at each center
!        n2p    = number of 2p functions of each type at each center
!        n3s    = number of 3s functions at each center
!        n3p    = number of 3p functions of each type at each center
!        n3dzr  = number of z**2-r**2 d functions at each center
!        n3dx2  = number of x**2-y**2 d functions at each center
!        n3dxy  = number of xy d functions at each center
!        n3dxz  = number of xz d functions at each center
!        n3dyz  = number of yz d functions at each center
!        n4s    = number of 4s functions at each center
!        n4p    = number of 4p functions of each type at each center
!
      common /basis/ zex(MBASIS,MWF), betaq, n1s(MCTYPE)
      common /basis/ n2s(MCTYPE),n2p(3,MCTYPE),n3s(MCTYPE),n3p(3,MCTYPE)
      common /basis/ n3dzr(MCTYPE),n3dx2(MCTYPE),n3dxy(MCTYPE),n3dxz(MCTYPE),n3dyz(MCTYPE)
      common /basis/ n4s(MCTYPE),n4p(3,MCTYPE)
      common /basis/ n4fxxx(MCTYPE),n4fyyy(MCTYPE),n4fzzz(MCTYPE),n4fxxy(MCTYPE),n4fxxz(MCTYPE)
      common /basis/ n4fyyx(MCTYPE),n4fyyz(MCTYPE),n4fzzx(MCTYPE),n4fzzy(MCTYPE),n4fxyz(MCTYPE)
      common /basis/ nsa(MCTYPE),npa(3,MCTYPE)
      common /basis/ ndzra(MCTYPE),ndx2a(MCTYPE),ndxya(MCTYPE),ndxza(MCTYPE),ndyza(MCTYPE)
