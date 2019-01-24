      program main

      implicit real*8(a-h,o-z)
      character*12 mode

      common /contr3/ mode

      mode='dmc         '

      open(45,file='output.log',status='unknown')

      call read_input

      call p2gtid('optwf:ioptwf',ioptwf,0,1)

      if(ioptwf.gt.0) then
       call optwf_matrix_corsamp
      else
       call dmc
      endif

      close(5)
      close(6)
      close(45)

      end
