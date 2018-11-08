      subroutine my_second (n,title)
c Prints out cpu and wall-clock time,
c both since beginning of run and since last call.
      implicit real*8 (a-h,o-z)

      character*6 title

      save icall,itim1,itim2,etim1,etim2
      data icall/0/

      call cpu_time(etim)
      call system_clock(itim,irate)
      if(icall.eq.0) then
        icall=1
        itim1=itim
        itim2=itim
        etim1=etim
        etim2=etim
      endif
      itimtot=(itim-itim1)/dfloat(irate)
      itimlast=(itim-itim2)/dfloat(irate)
      itim2=itim
      etimtot=etim-etim1
      etimlast=etim-etim2
      etim2=etim
      if(n.eq.1) write (6,'(''BEGINNING OF '',a6,'' CP, REAL TIME IS '',
     &2f10.2,2i6)') title,etimtot,etimlast,itimtot,itimlast
      if(n.eq.2) write (6,'(''END       OF '',a6,'' CP, REAL TIME IS '',
     &2f10.2,2i6)') title,etimtot,etimlast,itimtot,itimlast
      return
      end
