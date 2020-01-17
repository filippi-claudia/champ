c   /basis/
c        ncent  = number of centers
c        zex    = screening constants for each basis function
c        cent   = positions of each center
c        pecent = potential energy of the centers
c        znuc   = charge of the nuclei (centers)
c        n1s    = number of 1s functions at each center
c        n2s    = number of 2s functions at each center
c        n2p    = number of 2p functions of each type at each center
c        n3s    = number of 3s functions at each center
c        n3p    = number of 3p functions of each type at each center
c        n3dzr  = number of z**2-r**2 d functions at each center
c        n3dx2  = number of x**2-y**2 d functions at each center
c        n3dxy  = number of xy d functions at each center
c        n3dxz  = number of xz d functions at each center
c        n3dyz  = number of yz d functions at each center
c        n4s    = number of 4s functions at each center
c        n4p    = number of 4p functions of each type at each center
c
      common /basis/ zex(MBASIS,MWF),betaq
     &,n1s(MCTYPE)
     &,n2s(MCTYPE),n2p(3,MCTYPE),n3s(MCTYPE),n3p(3,MCTYPE)
     &,n3dzr(MCTYPE),n3dx2(MCTYPE),n3dxy(MCTYPE),n3dxz(MCTYPE),n3dyz(MCTYPE)
     &,n4s(MCTYPE),n4p(3,MCTYPE)
     &,n4fxxx(MCTYPE),n4fyyy(MCTYPE),n4fzzz(MCTYPE),n4fxxy(MCTYPE),n4fxxz(MCTYPE)
     &,n4fyyx(MCTYPE),n4fyyz(MCTYPE),n4fzzx(MCTYPE),n4fzzy(MCTYPE),n4fxyz(MCTYPE)
     &,nsa(MCTYPE),npa(3,MCTYPE)
     &,ndzra(MCTYPE),ndx2a(MCTYPE),ndxya(MCTYPE),ndxza(MCTYPE),ndyza(MCTYPE)
