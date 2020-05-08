      subroutine cuspinit4(iprin)
c Written by Cyrus Umrigar
      use jaspar4, only: a4, norda, nordb, nordc
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      common /cuspmat4/ d(NEQSX,MTERMS),iwc4(NEQSX),nterms

      if(nordc.eq.0) return

      do 10 n=1,2*(nordc-1)
        iwc4(n)=0
        do 10 i=1,MTERMS
   10   d(n,i)=0

      i=0
      do 20 n=2,nordc
        do 20 k=n-1,0,-1
          if(k.eq.0) then
            l_hi=n-k-2
           else
            l_hi=n-k
          endif
          do 20 l=l_hi,0,-1
            m=(n-k-l)/2
            if(2*m.eq.n-k-l) then
              i=i+1
              if(i.gt.MTERMS) stop 'nterms>MTERMS in cuspinit4'
              if(k.eq.1.and.iwc4(n-1).eq.0) iwc4(n-1)=i
              if(k.eq.0.and.iwc4(n+nordc-2).eq.0) iwc4(n+nordc-2)=i
              if(iprin.gt.1) write(6,'(9i4)') i,n,k,l,m
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
   20 continue
      nterms=i
      if(iprin.gt.0) then
        write(6,'(''# of e-e-n terms, nterms='',i5)') nterms
        if(iprin.gt.1) then
          write(6,'(''d matrix:'')')
          write(6,'(55i2)') (i,i=1,nterms)
          do 30 n=1,2*(nordc-1)
   30       write(6,'(55i2)') (nint(d(n,i)),i=1,nterms)
        endif
        write(6,'(''coefs. fixed by cusp conditions are'')')
        write(6,'(55i3)') (iwc4(i),i=1,2*(nordc-1))
      endif

      call checkdepend4(iprin)

      return
      end
c-----------------------------------------------------------------------
      subroutine checkdepend4(iprin)
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use jaspar4, only: a4, norda, nordb, nordc
      use optwf_nparmj, only: nparma, nparmb, nparmc, nparmf
      use optwf_parms, only: nparmd, nparme, nparmg, nparmj, nparml, nparms
      use optwf_wjas, only: iwjasa, iwjasb, iwjasc, iwjasf

      use vardep, only: cdep, iwdepend, nvdepend

      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'

      common /cuspmat4/ d(NEQSX,MTERMS),iwc4(NEQSX),nterms


      neqs=2*(nordc-1)

      do 5 i=1,neqs
        do 5 it=1,nctype
          nvdepend(i,it)=0
          do 5 l=1,nparmc(it)
            cdep(i,l,it)=0
  5         iwdepend(i,l,it)=0

c Figure out dependence of all dependent variables except that from
c the 2nd order e-n cusp cond.
      do 40 i=1,neqs
        if(i.eq.nordc) goto 40
        do 20 j=1,nterms
          if(dabs(d(i,j)).gt.1.d-10) then
            do 10 it=1,nctype
              do 10 l=1,nparmc(it)
                if(j.eq.iwjasc(l,it)) then
                  if(j.eq.iwc4(i)) stop 'do not vary dep. variable'
                  nvdepend(i,it)=nvdepend(i,it)+1
                  iwdepend(i,nvdepend(i,it),it)=l
                  cdep(i,nvdepend(i,it),it)=-d(i,j)/d(i,iwc4(i))
                endif
   10       continue
          endif
   20   continue
   40 continue

c Now do the one from the 2nd order e-n cusp cond.
c The dep. variable from the 2nd-order e-n cc depends directly on all the
c other dependent variables and only on the other dependent variables.
c Figure out what dependence on the independent variables is implied by that.

c Since it depends directly on all other dependent variables, it depends
c indirectly on all the independent variables that are being varied.
      do 50 it=1,nctype
      do 50 l=1,nparmc(it)
        nvdepend(nordc,it)=nvdepend(nordc,it)+1
   50   iwdepend(nordc,nvdepend(nordc,it),it)=l

      do 70 i=1,neqs
        if(i.eq.nordc) goto 70
        factor=-d(nordc,iwc4(i))/d(nordc,iwc4(nordc))
        do 60 j=1,nterms
          do 60 it=1,nctype
            do 60 l=1,nparmc(it)
              if(j.eq.iwjasc(l,it)) then
                cdep(nordc,l,it)= cdep(nordc,l,it)-
     &          factor*(d(i,j)/d(i,iwc4(i)))
              endif
   60   continue
   70 continue

      if(iprin.gt.1) then
      write(6,'(''i  it nvdepend iwdepend or cdep'')')
      do 80 i=1,neqs
        do 80 it=1,nctype
          write(6,'(i2,2i4,2x,50i4)') i,it,nvdepend(i,it),
     &   (iwdepend(i,j,it),j=1,nvdepend(i,it))
   80     write(6,'(i2,2i4,2x,50f4.0)') i,it,nvdepend(i,it),
     &   (cdep(i,j,it),j=1,nvdepend(i,it))
      endif

      return
      end

