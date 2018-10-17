      implicit real*8(a-h,o-z)
      character*10 file

c real*8 test
      call clear_locate
      call locate(y)
      y=2.d0
      z=y**2
      file='main'
      line=6
      call monitor(file,line)      
      y=3.d0
      line=9
      call monitor(file,line)      

c integer*4 test
c     call clear_locate
c     call locate_int(i)
      i=5
      j=i**2
c     file='main'
c     line=6
c     call monitor_int(file,line)
      i=7
c     line=9
c     call monitor_int(file,line)

      stop
      end
c-----------------------------------------------------------------------
c Peter's print-anything-from-anywhere routines.
c First call clear_locate to set the number of locations monitored.
c (This is redundant if only one set of locations are monitored since
c blockdata does the same thing.)
c Then call locate for each of the variables to be monitored.
c Then call monitor each time you want to monitor the variables.
c locate_int and monitor_int are the integer versions.

c **Warning: For the real*8 variables it assumes that the locations returned
c are either integers or half-integers when divided by 8.  This will in
c general not be true. In Linux g77 I think one can fix that by including
c -malign-double in FFLAGS.

c The idea is to find the location of any variable as an array element of the array x
c in monitor by using its address as supplied by the routine loc.  The %
c is tricky it has to do with call-by-name (which is FORTRANESE) vs
c call-by-value (C-ESE), depending on which loc function you happen to
c get.  You'll have to experiment with that.

c If you use a compiler that respects the cpp conventions you can use
c __FILE__ and __LINE__ in your code and have the (pre-)compiler
c substitute the file and line info in your code of the place where you
c use this.  That's how I use the arguments of monitor.

c 
c Peter.
c-----------------------------------------------------------------------
      block data locations_block_data
      parameter (MX=2)
      common/locations_block/locations,location(MX)
      data locations/0/
      end
c-----------------------------------------------------------------------
      subroutine clear_locate
      parameter (MX=2)
      common/locations_block/locations,location(MX)
      locations=0
      return
      end
c-----------------------------------------------------------------------
      subroutine locate(x)
      implicit real*8(a-h,o-z)
      parameter (MX=2)
      common/locations_block/locations,location(MX)
      locations=locations+1
      if(locations .gt. MX) stop 'locate: locations .gt. MX'
      location(locations)=%loc(x)
      write(6,*) location(locations),(location(locations)/8)*8
      return
      end
c-----------------------------------------------------------------------
      subroutine locate_int(i)
      implicit real*8(a-h,o-z)
      parameter (MX=2)
      common/locations_block/locations,location(MX)
      locations=locations+1
      if(locations .gt. MX) stop 'locate: locations .gt. MX'
      location(locations)=%loc(i)
      return
      end
c-----------------------------------------------------------------------
      subroutine monitor(file,line)
      implicit real*8(a-h,o-z)
      parameter (MX=2)
      common/locations_block/locations,location(MX)
      parameter (ireal8=8)
      character*(*) file
      dimension x(1)
      do i=1,locations
        print *,'monitor:',file,line,x((location(i)-%loc(x))/ireal8+1)
      enddo
      write(6,*) location(locations),(location(locations)/8)*8
      write(6,*) %loc(x),(%loc(x)/8)*8
      return
      end
c-----------------------------------------------------------------------
      subroutine monitor_int(file,line)
      implicit real*8(a-h,o-z)
      parameter (MX=2)
      common/locations_block/locations,location(MX)
      parameter (int4=4)
      character*(*) file
      dimension i(1)
      do j=1,locations
        print *,'monitor:',file,line,i((location(j)-%loc(i))/int4+1)
      enddo
      return
      end
