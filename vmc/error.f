      subroutine fatal_error(msg)

      implicit double precision (a-h,o-z)
      character msg*(*)

      write(6,'(''Fatal error: '',a)') msg
      stop

      return
      end
