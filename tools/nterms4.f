      function nterms4(nord)
      implicit real*8(a-h,o-z)

      i=0
      do 20 n=2,nord
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
            endif
   20 continue
      nterms4=i
      write(6,'(''nterms4='',i5)') nterms4
      return
      end
