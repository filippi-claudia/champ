!     tables for hartree and exchange potentials 
!     max interpolation order
      parameter (MXTABI=6)
!     max number of  hartree, exchangetables 
      parameter (MXHTAB=1,MXXTAB=1,MXTAB=MXHTAB+MXXTAB)
!     max number of gridpoints in hartree table
      parameter (MXHGRID=1)
!     max number of gridpoints in exchange tables
      parameter (MXXGRID=1)
!     flags whether tables are used and interpolation order
      common /tab000/ ihtree,ixchng,ihpol,ixpol
!     number of tables,  table sizes
      common /tab001/ nhtab,nxtab,ndimh(3,MXHTAB),ndimx(3,MXXTAB)
!     table bounds and spacing
      common /tab002/ rminh(3,MXHTAB), rmaxh(3,MXHTAB),drh(3,MXHTAB)
      common /tab003/ rminx(3,MXXTAB), rmaxx(3,MXXTAB),drx(3,MXXTAB)
!     tables
      common /tab004/ thtree(MXHGRID*MXHGRID*MXHGRID,MXHTAB)
      common /tab005/ txchng(MXXGRID*MXXGRID*MXXGRID,MXXTAB)
!     evaluation count
      common /tab006/ nevalu(2,MXTAB)
!     used part of tables 
      common /tab007/ tmin(3,MXTAB),tmax(3,MXTAB)
