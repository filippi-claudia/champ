c     open files
      parameter(MXFNL=4096)
      parameter(MXFLS=1024)
      dimension ifu(MXFLS)
      dimension istat(MXFLS)
      dimension itype(MXFLS)
      character ifn(MXFLS)*4096
      common /fle1/ nfle,iul,ifu,istat,itype
      common /fle2/ ifn
      common /fle2b/ inptu
C range of unit numbers for 'sr file'	 
      parameter(MINUNT=50,MXUNT=90)
C    current positions
      common /fle3/ icpos(MXUNT)
C    other flags
      common /fle4/ irewnd
C    internal buffers	
      parameter(IBUF1=7,IBUF2=8,NBUFFER=2,MBUFLNE=256,KBUFLNE=256)
      common /ibuf00/  ibufline(NBUFFER),kbufline(NBUFFER)
      character cbuffer(MBUFLNE,NBUFFER)*(KBUFLNE)
      common /ibuf01/ cbuffer
      common /ibuf02/ ibuflen(MBUFLNE,NBUFFER)
