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
!     number of tables,  table sizes
!     table bounds and spacing
!     tables
!     evaluation count
!     used part of tables 
