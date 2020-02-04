c     common blocks for p2 (part1)
C$Revision: 1.1 $
      character sfile*64
      character fmts(MXKEY)*(MXIDL)
      dimension nargs(MXKEY)
      character keys(MXKEY)*(MXIDL)
      dimension keylen(MXKEY)
c
      character lne*(MXLNE)
      character lbuf*(MXLNE)
      dimension itmp(MXF)
      dimension ftmp(MXF)
      character atmp*(MXLNE)
      dimension is1(MXF)
      dimension is2(MXF)
c
      common /p21/ fmts,keys,nargs,nkey,idum006
      common /p211/  keylen
      common /p22/ lne,atmp
      common /p23/ itmp,is1,is2
      common /p24/ ftmp
      common /p25/ il,iutmp
c
      character p2prmt*16
      common /p26/p2prmt
      common /p28/ iant,idum4
c
      common /p30/iimode,idum008
      common /p31/sfile
c
      character varnm(MXVAR)*(MXIDL)
      character varval(MXVAR)*(MXLNE)
      dimension ivrln(MXVAR)
      common /vars1/ varnm,varval
      common /vars2/ nvtop,ivrln,id020
c
c control
C    types 
C    1: i-loop
C 
      common /p2ct1/  lctrl,ktctrl(MXCTRL),kctrl(MXCTRL),
     $     ipctl(MXCP,MXCTRL),id014
      character apctl(MXCP,MXCTRL)*(MXIDL)
      common /p2ct2/  apctl
      common /p2ct3/  fpctl(MXCP,MXCTRL)
C nesting
      parameter(MXNEST=8)
      character co(MXNEST)
      common /p2nst1/co
      common /p2nst2/lnest,iio(MXNEST),ianto(MXNEST),
     $     lio(MXNEST),iuu(MXNEST)
