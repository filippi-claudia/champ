!     common blocks for p2 (part2)
! Revision: 1.1 $
! various flags
      common /p32/ ip2deb,inice,ip2ab,idum2
      character fsep*8
      common /sep00/ fsep,ifsep,id011
      character cmtchr,antchr,idum3*6
      common /p27/ cmtchr,antchr,idum3    
      character skipch*(MXIDL)
      common /skpp1/ skipch
      common /skpp2/ iskip,idum
      common /p2sto/ nans,id017
