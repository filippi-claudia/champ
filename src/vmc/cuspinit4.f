      module cuspinit4_mod
      contains
      subroutine cuspinit4(iprin)
c Written by Cyrus Umrigar
      use vmc_mod, only: mterms
      use jaspar4, only: nordc
      use cuspmat4, only: d, iwc4, nterms
      use contrl_file,    only: ounit
      implicit none

      integer :: i, iord, iprin, k, l
      integer :: l_hi, m, n





      if(nordc.eq.0) return

      do n=1,2*(nordc-1)
        iwc4(n)=0
        do i=1,mterms
        d(n,i)=0
        enddo
      enddo

      i=0
      do n=2,nordc
        do k=n-1,0,-1
          if(k.eq.0) then
            l_hi=n-k-2
           else
            l_hi=n-k
          endif
          do l=l_hi,0,-1
            m=(n-k-l)/2
            if(2*m.eq.n-k-l) then
              i=i+1
              if(i.gt.mterms) stop 'nterms>mterms in cuspinit4'
              if(k.eq.1.and.iwc4(n-1).eq.0) iwc4(n-1)=i
              if(k.eq.0.and.iwc4(n+nordc-2).eq.0) iwc4(n+nordc-2)=i
              if(iprin.gt.1) write(ounit,'(9i4)') i,n,k,l,m
              iord=l+2*m
              d(iord,i)=d(iord,i)-2*k
              if(l+m.gt.0) then
                iord=nordc-1+k+m
                d(iord,i)=d(iord,i)-(l+m)
              endif
              if(m.gt.0) then
                iord=nordc-1+k+l+m
                d(iord,i)=d(iord,i)-m
              endif
            endif
          enddo
        enddo
      enddo
      nterms=i
      if(iprin.gt.0) then
        write(ounit,'(''# of e-e-n terms, nterms='',i5)') nterms
        if(iprin.gt.1) then
          write(ounit,'(''d matrix:'')')
          write(ounit,'(55i2)') (i,i=1,nterms)
          do n=1,2*(nordc-1)
            write(ounit,'(55i2)') (nint(d(n,i)),i=1,nterms)
          enddo
        endif
        write(ounit,'(''coefs. fixed by cusp conditions are'')')
        write(ounit,'(55i3)') (iwc4(i),i=1,2*(nordc-1))
      endif

      call checkdepend4(iprin)

      return
      end
c-----------------------------------------------------------------------
      subroutine checkdepend4(iprin)
      use system, only: nctype
      use jaspar4, only: nordc
      use optwf_nparmj, only: nparmc
      use optwf_wjas, only: iwjasc

      use vardep, only: cdep, iwdepend, nvdepend
      use contrl_file,    only: ounit
      use cuspmat4, only: d, iwc4, nterms
      use precision_kinds, only: dp
      implicit none

      integer :: i, iprin, it, j, l
      integer :: neqs
      real(dp) :: factor






      neqs=2*(nordc-1)

      do i=1,neqs
        do it=1,nctype
          nvdepend(i,it)=0
          do l=1,nparmc(it)
            cdep(i,l,it)=0
            iwdepend(i,l,it)=0
          enddo
        enddo
      enddo

c Figure out dependence of all dependent variables except that from
c the 2nd order e-n cusp cond.
      do i=1,neqs
        if(i.eq.nordc) goto 40
        do j=1,nterms
          if(dabs(d(i,j)).gt.1.d-10) then
            do it=1,nctype
              do l=1,nparmc(it)
                if(j.eq.iwjasc(l,it)) then
                  if(j.eq.iwc4(i)) stop 'do not vary dep. variable'
                  nvdepend(i,it)=nvdepend(i,it)+1
                  iwdepend(i,nvdepend(i,it),it)=l
                  cdep(i,nvdepend(i,it),it)=-d(i,j)/d(i,iwc4(i))
                endif
              enddo
            enddo
          endif
        enddo
   40 continue
      enddo

c Now do the one from the 2nd order e-n cusp cond.
c The dep. variable from the 2nd-order e-n cc depends directly on all the
c other dependent variables and only on the other dependent variables.
c Figure out what dependence on the independent variables is implied by that.

c Since it depends directly on all other dependent variables, it depends
c indirectly on all the independent variables that are being varied.
      do it=1,nctype
      do l=1,nparmc(it)
        nvdepend(nordc,it)=nvdepend(nordc,it)+1
        iwdepend(nordc,nvdepend(nordc,it),it)=l
      enddo
      enddo

      do i=1,neqs
        if(i.eq.nordc) goto 70
        factor=-d(nordc,iwc4(i))/d(nordc,iwc4(nordc))
        do j=1,nterms
          do it=1,nctype
            do l=1,nparmc(it)
              if(j.eq.iwjasc(l,it)) then
                cdep(nordc,l,it)= cdep(nordc,l,it)-
     &          factor*(d(i,j)/d(i,iwc4(i)))
              endif
            enddo
          enddo
        enddo
   70 continue
      enddo

      if(iprin.gt.1) then
      write(ounit,'(''i  it nvdepend iwdepend or cdep'')')
      do i=1,neqs
        do it=1,nctype
          write(ounit,'(i2,2i4,2x,50i4)') i,it,nvdepend(i,it),
     &   (iwdepend(i,j,it),j=1,nvdepend(i,it))
          write(ounit,'(i2,2i4,2x,50f4.0)') i,it,nvdepend(i,it),
     &   (cdep(i,j,it),j=1,nvdepend(i,it))
        enddo
      enddo
      endif

      return
      end

      end module
