      subroutine jdqz (alpha, beta, eivec, wanted,
     $     n, shift, eps, kmax, jmax, jmin,
     $     method, m, l,
     $     maxnmv, maxstep, lock, order, testspace, work, lwork)
c
c     Programmer: Diederik R. Fokkema
c     Modified         : M. van Gijzen
c     Modified 05-24-96: M. Kooper: ra and rb, the Schur matrices of A and B, 
c              added, as well as the vectors sa and sb, which contain the
c              innerproducts of ra with Z and rb with Z. This is added to be
c              enable to compute the eigenvectors in EIVEC
c     Modification 08-27-96: Different computation of eigenvectors, MvG
c
c     .. Parameters ..
c
      implicit none
      integer gmres, cgstab
      parameter ( gmres = 1, cgstab = 2 )
      integer kmax, jmax, jmin, method, m, l, maxnmv, maxstep, order
      integer testspace
      real*8 eps, lock
      double complex shift
      integer n, lwork
      double complex work(n,*), eivec(n,*), alpha(*), beta(*)
      logical wanted
c
c     .. Local Parameters ..
c
      logical loop, try, found, ok
      integer ldvs, ldzwork, iseed(4)
      parameter (ldvs=50, ldzwork=4*ldvs)
      double complex ma(ldvs,ldvs), mb(ldvs,ldvs),
     $     zma(ldvs,ldvs), zmb(ldvs,ldvs),
     $     vsl(ldvs,ldvs), vsr(ldvs,ldvs),
     $     ra(ldvs,ldvs), rb(ldvs, ldvs),
     $     zwork(4*ldvs), aconv(ldvs), bconv(ldvs)

      integer ldqkz
      parameter (ldqkz=ldvs)
      integer ipivqkz(ldqkz)
      double complex mqkz(ldqkz,ldqkz), invqkz(ldqkz,ldqkz), f(ldqkz)

      integer i, j, k, info, mxmv, step
      integer d, u, v, w, aux, av, bv, q, z, kz, itmp, tp
      integer solvestep
      real*8 rnrm, rwork(3*ldvs), deps
      real*8 dtmp
      double complex zalpha, zbeta, targeta, targetb, evcond
      double complex shifta, shiftb
      double complex zero, one
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))
c
c     .. Functions ..
c
      real*8 dznrm2
c
c     .. Data ..
c
      data iseed /3,3,1966,29/

c     ... Executable Statements
c
c...     Are there errors in the parameters?

      if ( kmax .gt. jmax )
     $   call error('jdqz: kmax greater than jmax')
      if ( jmax .gt. 50 )
     $   call error('jdqz: jmax greater than 50')
      if ( method .ne. 1 .and. method .ne. 2 )
     $   call error('jdqz: illegal choice for solution method')
      if ( order .lt. -2 .or. order .gt. 2 )
     $   call error('illegal value for order, must be between -2 and 2')

c...     d   = rhs, these pointers refer to the columns of the workspace
      d   = 1
c...     Workspace for jdqzmv
      tp  = d+1
c...     u   = pointer to Krylov space GMRES(m) or Bi-CSTAB(l)
      u   = tp+1
c...     v   = pointer to search space JDQZ with max dimension jmax
      if ( method .eq. gmres ) then
	 v   = u+m+1
      else if ( method .eq. cgstab ) then
	 v   = u+2*l+6
      end if
c...     w   = pointer to test subspace JDQZ with max dimension jmax
      w   = v+jmax
c...     av  = pointer to subspace AV with max dimension jmax
      av  = w+jmax
c...     bv  = pointer to subspace BV with max dimension jmax
      bv  = av+jmax
c...     aux =
      aux = bv+jmax
c...     q   = pointer to search Schur basis in JDQZ with max dimension kmax
      q   = aux+jmax
c...     z   = pointer to test Schur basis in JDQZ with max dimension kmax
      z   = q+kmax
c...     kz  = pointer to matrix K^{-1}Z_k
      kz  = z+kmax
      if (kz+kmax-1.gt.lwork) call error ('qz: memory fault')
c
c     --- initialize loop
c
      ok = .true.

      evcond = dcmplx(sqrt(abs(shift)**2+abs(one)**2))
      shifta = shift/evcond
      shiftb = one/evcond

      targeta = shifta
      targetb = shiftb

      zalpha = shifta
      zbeta = shiftb

      step = 0
      deps = dble(one)
      mxmv = 0

      solvestep = 0

      j = 0
      k = 0

c
c     --- loop
c
 100  continue
      loop = (k.lt.kmax.and.step.lt.maxstep)
      if (loop) then
	 step = step+1
         write(6,'(''iteration : '',i4)') step

	 solvestep = solvestep+1
	 if (j.eq.0) then
	    call zlarnv(2, iseed, n, work(1,v+j))
	    call zlarnv(2, iseed, n, work(1,w+j))
	    do i=1,n
	       dtmp = dble(work(i,v+j))
	       work(i,v+j) = dcmplx(dtmp,0d0)
	       dtmp = dble(work(i,w+j))
	       work(i,w+j) = dcmplx(dtmp,0d0)
	    enddo
	 else
	    mxmv = maxnmv
	    deps = 2d0**(-solvestep)
	    if (j.lt.jmin) then
	       mxmv = 1
	       call zgmres (n, work(1,v+j), work(1,d), m, deps,
     $              mxmv, zalpha, zbeta, k+1,
     $              work(1,kz), work(1,q), invqkz, ldqkz,
     $              ipivqkz, f, work(1,u), work(1,tp) )
	    elseif (method.eq.gmres) then
	       mxmv = m
	       call zgmres (n, work(1,v+j), work(1,d), m, deps,
     $              mxmv, zalpha, zbeta, k+1,
     $              work(1,kz), work(1,q), invqkz, ldqkz,
     $              ipivqkz, f, work(1,u), work(1,tp) )
	    elseif (method.eq.cgstab) then
	       call zcgstabl (n, work(1,v+j), work(1,d), l,
     $              deps, mxmv, zalpha,
     $              zbeta, k+1, work(1,kz), work(1,q), invqkz,
     $              ldqkz, ipivqkz, f, work(1,u), 2*l+6)
	    endif
	 endif
	 j = j+1

	 call zmgs (n, j-1, work(1,v), work(1,v+j-1), 1 )
	 call zmgs (n, k, work(1,q), work(1,v+j-1), 1 )

	 if (testspace.eq.1) then
 	    call jdqzmv (n, work(1,v+j-1), work(1,w+j-1), work(1,tp), 
     $           -dconjg(shiftb), dconjg(shifta))
	 elseif (testspace.eq.2) then
 	    call jdqzmv (n, work(1,v+j-1), work(1,w+j-1), work(1,tp),
     $           -dconjg(zbeta), dconjg(zalpha))
	 elseif (testspace.eq.3) then
 	    call jdqzmv (n, work(1,v+j-1), work(1,w+j-1), work(1,tp),
     $           shifta, shiftb)
	 endif

         call zmgs (n, j-1, work(1,w), work(1,w+j-1), 1 )
 	 call zmgs (n, k, work(1,z), work(1,w+j-1), 1 )


         call amul(n, work(1,v+j-1), work(1,av+j-1))
         call bmul(n, work(1,v+j-1), work(1,bv+j-1))

         call makemm (n, j, work(1,w), work(1,av), ma, zma, ldvs)
         call makemm (n, j, work(1,w), work(1,bv), mb, zmb, ldvs)

         call zgegs ('v', 'v', j, zma, ldvs, zmb, ldvs,
     $        alpha, beta, vsl, ldvs, vsr, ldvs, zwork,
     $        ldzwork, rwork, info)

         try = .true.
 200     continue
         if (try) then
c
c           --- Sort the Petrov pairs ---
c
            call qzsort (targeta, targetb, j, zma, zmb, vsl, vsr,
     $           ldvs, alpha, beta, order)

            zalpha = zma(1,1)
            zbeta = zmb(1,1)

            evcond = dcmplx(sqrt(abs(zalpha)**2+abs(zbeta)**2))
c
c            --- compute new q ---
c
            call zgemv ('n', n, j, one, work(1,v), n, vsr(1,1),
     $           1, zero, work(1,q+k), 1)
            call zmgs (n, k, work(1,q), work(1,q+k), 1 )
c
c           --- compute new z ---
c
            call zgemv ('n', n, j, one, work(1,w), n, vsl(1,1),
     $           1, zero, work(1,z+k), 1)
            call zmgs (n, k, work(1,z), work(1,z+k), 1 )
c
c           --- Make new qkz ---
c
            call zcopy (n, work(1,z+k), 1, work(1,kz+k), 1)
            call precon (n, work(1,kz+k))
            call mkqkz (n, k+1, work(1,q), work(1,kz), mqkz, invqkz,
     $           ldqkz, ipivqkz)
c
c           --- compute new (right) residual= beta Aq - alpha Bq and
c               orthogonalize this vector on Z.
c
            call jdqzmv (n, work(1,q+k), work(1,d), work(1,tp),
     $           zalpha, zbeta)
            call zmgs (n, k, work(1,z), work(1,d), 0 )

            rnrm = dznrm2 (n, work(1,d), 1)/dble(evcond)

        write(6,'(''mxmv : '',i4)') mxmv
        write(6,'(''lambda('',i2,''): ('',1p,e11.4,'','',e11.4,'' )'')') k,zalpha/zbeta

            if (rnrm.lt.lock.and.ok) then
               targeta = zalpha
               targetb = zbeta
               ok = .false.
            endif

            found = (rnrm.lt.eps.and.
     $           (j.gt.1.or.k.eq.kmax-1))
            try =  found

            if (found) then

c              --- increase the number of found evs by 1 ---
               k = k+1

c              --- store the eigenvalue
               aconv(k) = zalpha
               bconv(k) = zbeta

               solvestep = 0
               if (k.eq.kmax) goto 100
               call zgemm ('n', 'n', n, j-1, j, one, work(1,v), n,
     $              vsr(1,2), ldvs, zero, work(1,aux), n)
               itmp = v
               v = aux
               aux = itmp
               call zgemm ('n', 'n', n, j-1, j, one, work(1,av), n,
     $              vsr(1,2), ldvs, zero, work(1,aux), n)
               itmp = av
               av = aux
               aux = itmp
               call zgemm ('n', 'n', n, j-1, j, one, work(1,bv), n,
     $              vsr(1,2), ldvs, zero, work(1,aux), n)
               itmp = bv
               bv = aux
               aux = itmp
               call zgemm ('n', 'n', n, j-1, j, one, work(1,w), n,
     $              vsl(1,2), ldvs, zero, work(1,aux), n)
               itmp = w
               w = aux
               aux = itmp
               j = j-1
               call zlacpy ('a', j, j, zma(2,2), ldvs, ma, ldvs)
               call zlacpy ('a', j, j, ma, ldvs, zma, ldvs)
               call zlacpy ('a', j, j, zmb(2,2), ldvs, mb, ldvs)
               call zlacpy ('a', j, j, mb, ldvs, zmb, ldvs)
               call zlaset ('a', j, j, zero, one, vsr, ldvs)
               call zlaset ('a', j, j, zero, one, vsl, ldvs)
               targeta = shifta
               targetb = shiftb
               ok = .true.
               mxmv = 0
               deps = dble(one)
            else if (j.eq.jmax) then
               call zgemm ('n', 'n', n, jmin, j, one, work(1,v), n,
     $              vsr, ldvs, zero, work(1,aux), n)
               itmp = v
               v = aux
               aux = itmp
               call zgemm ('n', 'n', n, jmin, j, one, work(1,av), n,
     $              vsr, ldvs, zero, work(1,aux), n)
               itmp = av
               av = aux
               aux = itmp
               call zgemm ('n', 'n', n, jmin, j, one, work(1,bv), n,
     $              vsr, ldvs, zero, work(1,aux), n)
               itmp = bv
               bv = aux
               aux = itmp
               call zgemm ('n', 'n', n, jmin, j, one, work(1,w), n,
     $              vsl, ldvs, zero, work(1,aux), n)
               itmp = w
               w = aux
               aux = itmp
               j = jmin
               call zlacpy ('a', j, j, zma, ldvs, ma, ldvs)
               call zlacpy ('a', j, j, zmb, ldvs, mb, ldvs)
               call zlaset ('a', j, j, zero, one, vsr, ldvs)
               call zlaset ('a', j, j, zero, one, vsl, ldvs)
            endif
            goto 200
         endif
         goto 100
      endif

      write(6,'(''number of steps : '',i4)') step
c
c...     Did enough eigenpairs converge?
      kmax = k

      if ( wanted ) then
c...        Compute the Schur matrices if the eigenvectors are
c...        wanted, work(1,tp) is used for temporary storage
c...        Compute RA:
         call zlaset ('l', k, k, zero, zero, ra, ldvs)    
         do i = 1, k
            call amul( n, work(1,q+i-1), work(1,tp) )
            call zgemv( 'c', n, i, one, work(1,z), n, work(1,tp), 1, 
     $                  zero, ra(1,i), 1 )
         end do
c...        Compute RB:
         call zlaset ('l', k, k, zero, zero, rb, ldvs)    
         do i = 1, k
            call bmul( n, work(1,q+i-1), work(1,tp) )
            call zgemv( 'c', n, i, one, work(1,z), n, work(1,tp), 1, 
     $                  zero, rb(1,i), 1 )
         end do

c        --- The eigenvectors RA and RB  belonging to the found eigenvalues
c            are computed. The Schur vectors in VR and VS are replaced by the
c            eigenvectors of RA and RB
         call zgegv('N','V',k,ra, ldvs, rb, ldvs, alpha, beta,vsl, ldvs,
     $              vsr,ldvs, zwork,ldzwork,rwork,info)
c        --- Compute the eigenvectors belonging to the found eigenvalues
c            of A and put them in EIVEC
         call zgemm('n', 'n', n, k, k, one, work(1,q), n,
     $              vsr, ldvs, zero, eivec, n)
      else
c
c...        Store the Schurvectors in eivec:
         call zcopy( k*n, work(1,q), 1, eivec, 1 )
         call zcopy( k, aconv, 1, alpha, 1 )
         call zcopy( k, bconv, 1, beta, 1 )
      end if

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine jdqzmv (n, x, y, work, alpha, beta)
c
c     Coded by Diederik Fokkema
c
c     .. Parameters ..
c
      implicit none
      integer n
      double complex alpha, beta, x(*), y(*), work(*) 
c
c     .. Local ..
c
      integer i
c     .. Executable statements ..
c
      call amul( n, x, work )
      call bmul( n, x, y )
      do i = 1, n
         y(i) = beta*work(i) - alpha*y(i)
      end do
c
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zcgstabl (n, x, r, l, eps, mxmv,
     $     zalpha, zbeta, nk, kz, qq, invqkz, ldqkz, jpiv, f,
     $     work, lwork)
c
c     Programmer: Diederik R. Fokkema
c
c
      implicit none
c
c     .. Parameters ..
c
      integer l, n, mxmv, nk, lwork, ldqkz, jpiv(*)
      real*8  eps
      double complex x(*), r(*), kz(n,*), qq(n,*), work(n,*),
     $     zalpha, zbeta, invqkz(ldqkz,*), f(*)

c
c     .. Local ..
c
      integer mxl
      parameter (mxl = 32)

      integer i, j, k, m, ipiv(mxl), nmv, info

      integer u, q, w, rr, xp, bp
      logical rcomp, xadd

      real*8 maxnorm, delta, bpnorm, tr0norm, trlnorm
      double complex varrho, hatgamma

      real*8 rnrm, rnrm0, eps1, dznrm2
      double complex alpha, beta, omega, gamma, rho0, rho1,
     $     yr(mxl), ys(mxl), z(mxl,mxl), ztmp

      double complex zero, one
      parameter (zero = (0.0d0,0.0d0), one = (1.0d0,0.0d0))

      double complex zdotc
c
c     .. Executable statements ..
c
      if (mxmv.eq.0) return

      if (l.gt.mxl)
     $     call error ('l exceeds mxl (zcgstabl)')

      u  = 1
      q  = u + (l+1)
      rr = q + (l+1)
      xp = rr + 1
      bp = xp + 1
      w = bp + 1

      if (w.gt.lwork)
     $     call error ('workspace too small (zcgstabl)')

c
c     --- set x to zero and compute first residue
c
      call zzeros (n, x)
      call zscal (n, -one, r, 1)
      call zcopy (n, r, 1, work(1,rr), 1)
      call psolve (n, r, nk, qq, kz, invqkz, ldqkz,
     $     jpiv, f)

c
c     --- initialize loop
c
      nmv = 0

      rnrm0 = dznrm2 (n, r, 1)
      rnrm = rnrm0
      eps1 = eps * rnrm0

      call zcopy (n, x, 1, work(1,xp), 1)
      call zcopy (n, r, 1, work(1,bp), 1)
      maxnorm = 0d0
      bpnorm = rnrm
      rcomp = .false.
      xadd = .false.
      delta = 1d-2

      m = 0
      rho0 = one
      alpha = one
      beta = zero
      omega = one

      call zzeros (n*(l+1), work(1,u))
      call zzeros (n*(l+1), work(1,q))
      call zcopy (n, r, 1, work(1,q), 1)
c
c     --- loop
c
 1000 continue
      m = m + l
c
c     --- BiCG part
c
      rho0 = -omega * rho0
      do k=1,l
         rho1 = zdotc (n, work(1,rr), 1, work(1,q+k-1), 1)
         beta = alpha * (rho1 / rho0)
         rho0 = rho1
         beta = beta
         do j=0,k-1
            call zxpay (n, work(1,q+j), 1, (-beta), work(1,u+j), 1)
         enddo
         call jdqzmv (n, work(1,u+k-1), work(1,u+k), work(1,w),
     $        zalpha, zbeta)
         call psolve (n, work(1,u+k), nk, qq, kz,
     $        invqkz, ldqkz, jpiv, f)
         gamma = zdotc (n, work(1,rr), 1, work(1,u+k), 1)
         alpha = rho0 / gamma
         do j=0,k-1
            call zaxpy (n, (-alpha), work(1,u+j+1), 1, work(1,q+j), 1)
         enddo
         call jdqzmv (n, work(1,q+k-1), work(1,q+k), work(1,w),
     $        zalpha, zbeta)
         call psolve (n, work(1,q+k), nk, qq, kz,
     $        invqkz, ldqkz, jpiv, f)
         call zaxpy (n, alpha, work(1,u), 1, x, 1)
         rnrm = dznrm2 (n, work(1,q), 1)
         maxnorm = max (maxnorm, rnrm)
         nmv = nmv+2
      enddo
c
c     --- MR part + Maintaining the convergence
c
      do i=1,l-1
         do j=1,i
            ztmp = zdotc (n, work(1,q+i), 1, work(1,q+j), 1)
            z(i,j) = ztmp
            z(j,i) = dconjg(ztmp)
         enddo
         yr(i) = zdotc (n, work(1,q+i), 1, work(1,q), 1)
         ys(i) = zdotc (n, work(1,q+i), 1, work(1,q+l), 1)
      enddo
      call zgetrf (l-1, l-1, z, mxl, ipiv, info)
      call zgetrs ('n', l-1, 1, z, mxl, ipiv, yr, mxl, info)
      call zgetrs ('n', l-1, 1, z, mxl, ipiv, ys, mxl, info)
      call zcopy (n, work(1,q), 1, r, 1)
      call zgemv ('n', n, l-1, (-one), work(1,q+1), n, yr, 1, one,
     $     r, 1)
      call zgemv ('n', n, l-1, (-one), work(1,q+1), n, ys, 1, one,
     $     work(1,q+l), 1)
      tr0norm = dznrm2 (n, r, 1)
      trlnorm = dznrm2 (n, work(1,q+l), 1)
      varrho = zdotc (n, work(1,q+l), 1, r, 1) / (tr0norm*trlnorm)
      hatgamma = varrho/abs(varrho) * max(abs(varrho),7d-1)
      hatgamma = (tr0norm/trlnorm)*hatgamma
      yr(l) = zero
      ys(l) = -one
      call zaxpy (l, (-hatgamma), ys, 1, yr, 1)

      omega = yr(l)
      call zgemv ('n', n, l, one, work(1,q), n, yr, 1, one, x, 1)
      call zgemv ('n', n, l, (-one), work(1,u+1), n, yr, 1, one,
     $     work(1,u), 1)
      call zaxpy (n, (-hatgamma), work(1,q+l), 1, r, 1)
      call zcopy (n, r, 1, work(1,q), 1)
c
c     --- reliable update
c
      rnrm = dznrm2 (n, work(1,q), 1)
      maxnorm = max (maxnorm, rnrm)
      xadd = (rnrm.lt.delta*rnrm0.and.rnrm0.lt.maxnorm)
      rcomp = ((rnrm.lt.delta*maxnorm.and.rnrm0.lt.maxnorm).or.xadd)

      if (rcomp) then
         call jdqzmv (n, x, work(1,q), work(1,w),
     $        zalpha, zbeta)
         call psolve (n, work(1,q), nk, qq, kz,
     $        invqkz, ldqkz, jpiv, f)
         call zxpay (n, work(1,bp), 1, -one, work(1,q), 1)
         maxnorm = rnrm
         if (xadd) then
            call zaxpy (n, one, x, 1, work(1,xp), 1)
            call zzeros (n, x)
            call zcopy (n, work(1,q), 1, work(1,bp), 1)
            bpnorm = rnrm
         endif
      endif
         
      if (nmv.lt.mxmv .and. rnrm.gt.eps1) goto 1000

      call zaxpy (n, one, work(1,xp), 1, x, 1)
c
c     --- return
c
      mxmv = nmv
      eps = rnrm/rnrm0

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zgmres (n, x, r, mxm, eps, mxmv, 
     $     alpha, beta, k, kz, q, invqkz, ldqkz, ipiv, f, v, tp)
c
c     Programmer: Diederik R. Fokkema
c
c
      implicit none
c
c     .. Parameters ..
c
      integer mxm, mxmv, n, k, ldqkz, ipiv(*)
      real*8 eps
      double complex x(*), r(*), kz(n,*), q(n,*), v(n,*),
     $     alpha, beta, invqkz(ldqkz,*), f(*), tp(*)

ctex@ \begin{manpage}{ZGMRES} 
ctex@
ctex@ \subtitle{ZGMRES} 
ctex@    ZGMRES -- Generalized Minimal Residual
ctex@    iterative method for solving linear systems $\beta A-\alpha B = -r$.
ctex@    This subroutine in specilized for use in JDQZ.
ctex@ 
ctex@ \subtitle{Declaration}
ctex@ \function{subroutine zgmres (n, x, r, mxm, eps, mxmv, a, ka, b, kb,
ctex@   alpha, beta, k, kz, mqkz, zmqkz, ldvs, q,
ctex@   lu, klu, dlu, v)}
ctex@
ctex@ \subtitle{Parameters}
ctex@    \variable{integer n} 
ctex@       On entry, n specifies the dimension of the matrix A.
ctex@       
ctex@    \variable{x} 
ctex@       double complex array of size n. 
ctex@       On exit, x is overwritten by the approximate solution.
ctex@
ctex@    \variable{r}
ctex@       double complex array of size n. Before entry, the array r 
ctex@       must contain the righthandside of the linear problem Ax=r. 
ctex@       Changed on exit.
ctex@ 
ctex@    \variable{integer mxm} 
ctex@       On entry, mxm specifies the degree of the Minimal Residual
ctex@       polynomial.
ctex@
ctex@    \variable{{real*8} eps}
ctex@       On entry, eps specifies the stop tolerance. On exit, eps contains
ctex@       the relative norm of the last residual.
ctex@
ctex@    \variable{integer mxmv}
ctex@       On Entry, mxmv specifies the maximum number of matrix 
ctex@       multiplications. On exit, mxmv contains the number of matrix
ctex@       multiplications performed by the method.
ctex@
ctex@    \variable{{double complex} zalpha}
ctex@       On entry, zalpha specifies $\alpha$. Unchanged on exit.
ctex@ 
ctex@    \variable{{double complex} zbeta}
ctex@       On entry, zbeta specifies $\beta$. Unchanged on exit.
ctex@ 
ctex@    \variable{integer k}
ctex@       On entry, k specifies the number of columns of the arrays
ctex@       kz and q.
ctex@
ctex@    \variable{z}
ctex@       double complex array z, of size (n,k). On entry the array z
ctex@       must contain the preconditioned matrix Z.
ctex@
ctex@    \variable{mqkz}
ctex@       double complex array mqkz, of size (ldvs,k). On entry the array 
ctex@       mqkz must contain the matrix Q'*KZ.
ctex@
ctex@    \variable{zmqkz}
ctex@       double complex array zmqkz, of size (ldvs,k). Workspace. Used to
ctex@       copy mqkz.
ctex@
ctex@    \variable{q}
ctex@       double complex array q, of size (n,k). On entry the array q
ctex@       must contain the preconditioned matrix Q.
ctex@
ctex@    \variable{v}
ctex@       double complex array of size (n,mxm+1). Workspace.
ctex@
ctex@ \subtitle{Description}
ctex@    ***
ctex@
ctex@ \subtitle{See Also}
ctex@    ***
ctex@
ctex@ \subtitle{References}
ctex@    ***
ctex@ 
ctex@ \subtitle{Bugs}
ctex@    ***
ctex@
ctex@ \subtitle{Author}
ctex@     Diederik R.\ Fokkema
ctex@
ctex@ \end{manpage}
ctex@ \begin{verbatim}
ctex@    % actual code
ctex@ \end{verbatim}
ctex@
c
c     .. Local ..
c
      logical restrt, loop
      integer maxm
      parameter (maxm = 100)

      integer i, m, m1, nmv
      real*8 rnrm0, rnrm, eps1, c(maxm)
      double complex hh(maxm,maxm-1), rs(maxm), s(maxm), y(maxm), rcs
      double complex zero, one
      parameter (zero = (0.0d0,0.0d0), one = (1.0d0,0.0d0))

      real*8  dznrm2
      double complex ztmp, zdotc
c
c     .. Executable Statements ..
c
      if ( mxm .gt. maxm-1 )
     $     call error ('mxm larger than maxm (zgmres)')
c
c     --- Initialize first residue
c
      call zzeros (n, x)
      call zscal (n, -one, r, 1)
      call psolve (n, r, k, q, kz, invqkz, ldqkz,
     $     ipiv, f)
c
c     --- initialize loop
c
      rnrm0  = dznrm2 (n,r,1)
      rnrm = rnrm0
      eps1  = eps * rnrm

      nmv = 0

      call zcopy (n, r, 1, v(1,1), 1)
c         
c     --- outer loop
c
 1000 restrt = ( nmv .lt. mxmv .and. rnrm .gt. eps1 )
      if ( restrt ) then
         ztmp = one / rnrm
         call zscal (n, ztmp, v(1,1), 1)
         rs(1) = rnrm
c
c     --- inner loop
c
         m = 0
 2000    loop = (nmv.lt.mxmv .and. m.lt.mxm .and. rnrm.gt.eps1)
         if (loop) then
            m  = m + 1
            m1 = m + 1
            call jdqzmv (n, v(1,m), v(1,m1), tp, alpha, beta)
            call psolve (n, v(1,m1), k, q, kz, invqkz,
     $           ldqkz, ipiv, f)
            nmv = nmv+1 
            do i = 1,m
               ztmp = zdotc (n, v(1,i), 1, v(1,m1), 1)
               hh(i,m) = ztmp
               call zaxpy (n, (-ztmp), v(1,i), 1, v(1,m1), 1)
            enddo
            ztmp = dznrm2( n, v(1,m1), 1 )
            hh(m1,m) = ztmp
            call zscal ( n, (one/ztmp), v(1,m1), 1 )
            do i = 1,m-1
               call zrot (1, hh(i,m), 1, hh(i+1,m), 1, c(i), s(i))
            enddo
            call zlartg (hh(m,m), hh(m1,m), c(m), s(m), rcs )
            hh(m,m) = rcs
            hh(m1,m) = zero
            rs(m1) = zero
            call zrot (1, rs(m), 1, rs(m1), 1, c(m), s(m))
            rnrm = abs (rs(m1))
            goto 2000
         endif
c
c     --- compute approximate solution x
c
         call zcopy ( m, rs, 1, y, 1 )
         call ztrsv ( 'u', 'n', 'n', m, hh, maxm, y, 1 )
         call zgemv ( 'n', n, m, one, v, n, y, 1, one, x, 1 )
c
c     --- compute residual for restart
c
         call jdqzmv (n, x, v(1,2), tp, alpha, beta)
         call psolve (n, v(1,2), k, q, kz, invqkz,
     $        ldqkz, ipiv, f)
         call zcopy (n, r, 1, v(1,1), 1)
         call zaxpy (n, -one, v(1,2), 1, v(1,1), 1)
         rnrm = dznrm2 (n, v(1,1), 1)

         goto 1000
      endif
c
c     --- return
c
      eps = rnrm/rnrm0
      mxmv = nmv

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
