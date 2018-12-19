      subroutine fatal_error(msg)

      character msg*(*)

      write(6,'(''Fatal error: '',a)') msg
      stop

      return
      end
