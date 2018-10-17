      program main

      implicit real*8(a-h,o-z)
      character*12 mode
      character method*20

      common /contr3/ mode

      mode='vmc         '

      open(45,file='output.log',status='unknown')

      call read_input

      call p2gtid('optwf:ioptwf',ioptwf,0,1)
      call p2gtad('optwf:method',method,'linear',1)
      call p2gtid('optwf:idl_flag',idl_flag,0,1)

      if(ioptwf.gt.0) then
        if(idl_flag.gt.0) then
          call dl_optwf
        elseif(method.eq.'sr_n'.or.method.eq.'lin_d'.or.method.eq.'mix_n') then
          call sr_optwf
        else
         call optwf
       endif
      else
       call vmc
      endif

      close(5)
      close(6)
      close(45)

      end
