cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     placeholders for wave function optimization related subroutines
c     (c) Friedemann Schautz, 2002
c     Modified by C. Filippi
c     dummy versions for use in DMC code
c     where optwf makes no sense 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine do_read_multiple_cistates(iu,ns)
      implicit real*8(a-h,o-z)
      character line*800

      do j=1,ns
       read(iu,*) line
       call incpos(iu,itmp,1)
      enddo
      call p2chkend(iu, 'multiple_cistates')

      end
