      subroutine setrn(iseed)

      use rnyucm, only : switch_rng
      implicit none

      integer, intent(in) :: iseed(4)

      if(switch_rng.eq.0) then
            call setrn_rannyu_reference(iseed)
      else
            call setrn_rannyu_std(iseed)
      endif
      endsubroutine
c-------------------------------------------

      subroutine setrn_rannyu_reference(iseed)
c NYU linear congruential random number generator.
c Uses 48 bits, rather than the usual 32 bits.

      use rnyucm, only: ll
      implicit none

      integer :: i
      integer, intent(in) :: iseed(4)

      do i=1,4
         ll(i)=iseed(i)
      enddo
      ll(4)=2*(ll(4)/2)+1
      return
      end
c----------------------------------------------
      subroutine setrn_rannyu_std(iseed)
      use precision_kinds, only: dp
      use contrl_file,    only: ounit
      implicit none

      integer, intent(in) :: iseed(4)
      integer size
      integer, allocatable :: new_seed(:)
      integer i

      ! determine the seed size
      call random_seed(size = size)

      ! allocate the new seed
      allocate(new_seed(size), source=0)

      ! input the seed
      do i=1,min(size,4)
            new_seed(i) = iseed(i)
      end do

      ! check the new seed
      write(ounit, *) 'seed size', size
      do i=1,size
            write(ounit, *) 'seed     ', i, new_seed(i)
      end do

      ! set the new seed
      call random_seed(put=new_seed)

      ! deallocate
      deallocate(new_seed)
      return
      end

c----------------------------------------------------------

      function rannyu(idum)
      use precision_kinds, only : dp
      use rnyucm, only : switch_rng

      implicit none
      integer, intent(in) :: idum
      real(dp) :: rannyu

      if (switch_rng.eq.0) then
            call rannyu_reference(rannyu, idum)
      else
            call random_number(rannyu)
      endif

      return
      end


c-----------------------------------------------------------------------
      subroutine rannyu_reference(val, idum)
      use rnyucm, only: ll, mm
      use precision_kinds, only: dp
      implicit none

      real(dp), intent(inout) :: val
      integer, intent(in) :: idum
      integer :: i1, i2, i3, i4
      integer, parameter :: itwo12 = 4096

      real(dp), parameter :: two12i = 2.44140625d-4


c     On bat it is more efficient to not precompute 2**12 and 2**-12
c     Therefore use original verion of code instead of this one.
      i1=ll(1)*mm(4)+ll(2)*mm(3)+ll(3)*mm(2)+ll(4)*mm(1)
      i2=ll(2)*mm(4)+ll(3)*mm(3)+ll(4)*mm(2)
      i3=ll(3)*mm(4)+ll(4)*mm(3)
      i4=ll(4)*mm(4)
      ll(4)=mod(i4,itwo12)
      i3=i3+i4/itwo12
      ll(3)=mod(i3,itwo12)
      i2=i2+i3/itwo12
      ll(2)=mod(i2,itwo12)
      ll(1)=mod(i1+i2/itwo12,itwo12)
      val=two12i *(dfloat(ll(1))+
     &       two12i*(dfloat(ll(2))+
     &       two12i*(dfloat(ll(3))+
     &       two12i*(dfloat(ll(4))))))
      return
      end
c-----------------------------------------------------------------------

      subroutine savern(iseed)
      use rnyucm, only: ll
      implicit none

      integer :: i


      integer iseed(4)
      do i=1,4
         iseed(i)=ll(i)
      enddo
      return
      end
