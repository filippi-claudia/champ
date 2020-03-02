!
! Copyright (C) 2003-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE regterg( nparm, nparmx, nvec, nvecx, evc, ethr, &
                    e, btype, notcnv, dav_iter, ipr, idtask )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an symmetric matrix, evc is a real vector
  ! ... (real wavefunctions with only half plane waves stored)
  !
  IMPLICIT NONE
  !
  include 'mpif.h'
  !
  INTEGER, INTENT(IN) :: nparm, nparmx, nvec, nvecx, ipr
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! printing level
  REAL*8, INTENT(INOUT) :: evc(nparmx,nvec)
    !  evc   contains the  refined estimates of the eigenvectors
  REAL*8, INTENT(IN) :: ethr
    ! energy threshold for convergence: root improvement is stopped,
    ! when two consecutive estimates of the root differ by less than ethr.
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  REAL*8, INTENT(OUT) :: e(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer  number of iterations performed
    ! number of unconverged roots
  !
  INTEGER, INTENT(IN) :: idtask
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 200
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, n, m, nb1, ibnd
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
    ! counter on the bands
  INTEGER :: ierr
  INTEGER :: i,j
  REAL*8, ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:), ew(:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! eigenvectors of the Hamiltonian
    ! eigenvalues of the reduced hamiltonian
  REAL*8, ALLOCATABLE :: hhr(:,:), ssr(:,:)
  REAL*8, ALLOCATABLE :: psi(:,:), hpsi(:,:), spsi(:,:)
  REAL*8, ALLOCATABLE :: res(:,:)
  REAL*8, ALLOCATABLE :: res_norm(:)

    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  LOGICAL, ALLOCATABLE :: conv(:)
    ! true if the root is converged
  REAL*8 :: empty_ethr 
    ! threshold for empty bands
  REAL*8, EXTERNAL :: ddot
  !
  ! EXTERNAL  h_psi_lin_d, s_psi_lin_d, g_psi_lin_d
    ! h_psi_lin_d(nparm,nvec,psi,hpsi)
    !     calculates H|psi> 
    ! s_psi_lin_d(nparm,nvec,psi,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (nparmx,nvec)
    ! g_psi_lin_d(nparm,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
  IF ( nvec > nvecx / 2 ) CALL fatal_error( 'regter: nvecx is too small')
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  ALLOCATE( psi(  nparmx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL fatal_error( 'regterg: cannot allocate psi ')
  ALLOCATE( hpsi( nparmx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL fatal_error( 'regterg: cannot allocate hpsi ')
  !
  ALLOCATE( spsi( nparmx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL fatal_error( ' regterg: cannot allocate spsi ')
  !
  ALLOCATE( sr( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL fatal_error( 'regterg: cannot allocate sr ')
  ALLOCATE( hr( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL fatal_error( 'regterg: cannot allocate hr ')
  ALLOCATE( vr( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL fatal_error( 'regterg: cannot allocate vr ')
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL fatal_error( 'regterg: cannot allocate ew ')
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL fatal_error( 'regterg: cannot allocate conv ')
  !
  ALLOCATE( ssr( nvecx, nvecx ), STAT=ierr )
  ALLOCATE( hhr( nvecx, nvecx ), STAT=ierr )
  ALLOCATE( res(nparmx, nvecx), STAT=ierr )
  ALLOCATE( res_norm(nvec), STAT=ierr)
  !
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !
  spsi = 0.D0
  !
  hpsi = 0.D0
  psi  = 0.D0
  psi(:,1:nvec) = evc(:,1:nvec)
  !
  ! ... hpsi contains h times the basis vectors
  !
  CALL h_psi_lin_d( nparm, nvec, psi, hpsi )
  !
  CALL s_psi_lin_d( nparm, nvec, psi, spsi )
  !
  ! ... hr contains the projection of the hamiltonian onto the reduced
  ! ... space vr contains the eigenvectors of hr
  !
  IF(idtask.eq.0) then
  !
  hr(:,:) = 0.D0
  sr(:,:) = 0.D0
  vr(:,:) = 0.D0
  !
  CALL DGEMM( 'T', 'N', nbase, nbase, nparm, 1.D0 , &
              psi, nparmx, hpsi, nparmx, 0.D0, hr, nvecx )
  !
  CALL DGEMM( 'T', 'N', nbase, nbase, nparm, 1.D0, &
              psi, nparmx, spsi, nparmx, 0.D0, sr, nvecx )
  !
  ! ... diagonalize the reduced hamiltonian
  !
  IF(ipr.gt.1) then
    do i=1,nbase
      write(6,'(''REG HR '',100e10.3)') (hr(i,j),j=1,nbase)
    enddo
    do i=1,nbase
      write(6,'(''REG SR '',100e10.3)') (sr(i,j),j=1,nbase)
    enddo
  ENDIF
  !
  CALL rdiaghg( nbase, nvec, hr, sr, nvecx, ew, vr )
  !
  write(6,'(''LIN_D: EIG '',100e15.6)') (ew(j),j=1,nvec)
  IF(ipr.gt.1) then
    do i=1,nbase
      write(6,'(''REG VEC'',100e10.3)') (vr(i,j),j=1,nvec)
    enddo
  ENDIF
  !
  ! Claudia
  !
  e(1:nvec) = ew(1:nvec)
  !
  END IF ! idtask.eq.0
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     !
     IF(idtask.eq.0) then
     !
     dav_iter = kter
     write( 6,'(''REG: -----------------------------'')')
     write(6,'(''REG: Iteration: '', I10)') kter
     !
     np = 0
     !
     DO n = 1, nvec
        !
        IF ( .NOT. conv(n) ) THEN
           !
           ! ... this root not yet converged ... 
           !
           np = np + 1
           !
           ! ... reorder eigenvectors so that coefficients for unconverged
           ! ... roots come first. This allows to use quick matrix-matrix 
           ! ... multiplications to set a new basis vector (see below)
           !
           IF ( np /= n ) vr(:,np) = vr(:,n)
           !
           ! ... for use in g_psi
           !
           ew(nbase+np) = e(n)
           !   
        END IF
        !
     END DO
     !
     IF(ipr.gt.1) then
       write(6,*) 'Not converged',(conv(n),n=1,nvec)
       !
       write(6,'(''Expand with basis vectors '',i5)') notcnv
     ENDIF
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL DGEMM( 'N', 'N', nparm, notcnv, nbase, 1.D0, &
                    spsi, nparmx, vr, nvecx, 0.D0, psi(1,nb1), nparmx )
     !
     DO np = 1, notcnv
        !
        psi(:,nbase+np) = - ew(nbase+np) * psi(:,nbase+np)
        !
     END DO
     !
     CALL DGEMM( 'N', 'N', nparm, notcnv, nbase, 1.D0, &
                 hpsi, nparmx, vr, nvecx, 1.D0, psi(1,nb1), nparmx )
     !
     ! ... approximate inverse iteration
     !
!    CALL g_psi_lin_d( nparm, notcnv, nb1, psi(1,nb1), ew(nb1) )
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in 
     ! ... order to improve numerical stability of subspace diagonalization 
     ! ... (rdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
     DO n = 1, notcnv
        !
        ew(n) = ddot( nparm, psi(1,nbase+n), 1, psi(1,nbase+n), 1 )
        !
     END DO
     !
     DO n = 1, notcnv
        !
        psi(:,nbase+n) = psi(:,nbase+n) / SQRT( ew(n) )
        !
     END DO
     !
     END IF ! idtask.eq.0
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     nb1=nbase+1
     !
     CALL h_psi_lin_d( nparm, notcnv, psi(1,nb1), hpsi(1,nb1) )
     !
     CALL s_psi_lin_d( nparm, notcnv, psi(1,nb1), spsi(1,nb1) )
     !
     IF(idtask.eq.0) then
     !
     ! ... update the reduced hamiltonian
     !
     CALL DGEMM( 'T', 'N', nbase+notcnv, notcnv, nparm, 1.D0, psi, &
                 nparmx, hpsi(1,nb1), nparmx, 0.D0, hr(1,nb1), nvecx )
     !
     CALL DGEMM( 'T', 'N', nbase+notcnv, notcnv, nparm, 1.D0, psi, &
                 nparmx, spsi(1,nb1), nparmx, 0.D0, sr(1,nb1), nvecx )
     !
     nbase = nbase + notcnv
     !
     DO n = 1, nbase
        !
        sr(n,n)=sr(n,n)+1.d-06
        !
        DO m = n + 1, nbase
           !
           hr(m,n) = hr(n,m)
           sr(m,n) = sr(n,m)
           !
        END DO
        !
     END DO
     !
     ! ... diagonalize the reduced hamiltonian
     !
     IF(ipr.gt.1) then
       do i=1,nbase
         write(6,'(''REG HR '',100e10.3)') (hr(i,j),j=1,nbase)
       enddo
       do i=1,nbase
         write(6,'(''REG SR '',100e10.3)') (sr(i,j),j=1,nbase)
       enddo
     ENDIF
     !
     CALL rdiaghg( nbase, nvec, hr, sr, nvecx, ew, vr )
     !
     ! write(6,'(''LIN_D: EIG '',100e15.6)') (ew(j),j=1,nvec)
     IF(ipr.gt.1) then
       do i=1,nbase
         write(6,'(''REG VEC'',100e10.3)') (vr(i,j),j=1,nvec)
       enddo
     ENDIF
     !
     ! ... test for convergence
     !
    WHERE( btype(1:nvec) == 1 )
       !
       conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < ethr ) )
       !
    ELSEWHERE
       !
        conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < empty_ethr ) )
       !
     END WHERE
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
     e(1:nvec) = ew(1:nvec)
     
     CALL DGEMM( 'N', 'N', nparm, nbase, nbase, 1.D0, &
                    spsi, nparmx, vr, nvecx, 0.D0, res, nparmx )
     !
     DO np = 1, nbase
        !
        res(:,np) = - ew(np) * res(:,np)
        !
     END DO
     !
     CALL DGEMM( 'N', 'N', nparm, nbase, nbase, 1.D0, &
                 hpsi, nparmx, vr, nvecx, 1.D0, res, nparmx )     

     res_norm = norm2(res(:,:nvec),1)
     write(6,'(''REG: EIG '',100e15.6)') (e(j),j=1,nvec)
     write(6,'(''REG: RES '',100e15.6)') (res_norm(j),j=1,nvec)

!     WHERE( btype(1:nvec) == 1 )
        !
!        conv(1:nvec) = (  res_norm < ethr )
        !
!     ELSEWHERE
        !
!        conv(1:nvec) = ( res_norm < empty_ethr )
        !
!     END WHERE

!     notcnv = COUNT(.not. conv(:))

     !
     END IF ! idtask.eq.0
     !
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     !
     call MPI_BCAST(notcnv,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(nbase,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(dav_iter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     !
     IF ( notcnv == 0 .OR. &
          nbase+notcnv > nparm .OR. &
          nbase+notcnv > nvecx .OR. dav_iter == maxter ) THEN
        !
        IF(idtask.eq.0) &
        CALL DGEMM( 'N', 'N', nparm, nvec, nbase, 1.D0, &
                    psi, nparmx, vr, nvecx, 0.D0, evc, nparmx )
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           EXIT iterate
           !
        ELSE IF ( dav_iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           WRITE( 6, '(5X,"WARNING: ",I5, &
                &   " eigenvalues not converged in regterg")' ) notcnv
           !
           EXIT iterate
           !
        END IF
        !
        IF(idtask.eq.0) then
        !
!       if(ipr.gt.1) write(6,'(''Refresh, notcnv,nvec,nbase '',3i4)') notcnv,nvec,nbase
        write(6,'(''Refresh, notcnv,nvec,nbase '',3i4)') notcnv,nvec,nbase
        !
        ! ... refresh psi, H*psi and S*psi
        !
        psi(:,1:nvec) = evc(:,1:nvec)
        !
        CALL DGEMM( 'N', 'N', nparm, nvec, nbase, 1.D0, spsi, &
                    nparmx, vr, nvecx, 0.D0, psi(1,nvec+1), nparmx )
        !
        spsi(:,1:nvec) = psi(:,nvec+1:nvec+nvec)
        !
        CALL DGEMM( 'N', 'N', nparm, nvec, nbase, 1.D0, hpsi, &
                    nparmx, vr, nvecx, 0.D0, psi(1,nvec+1), nparmx )
        !
        hpsi(:,1:nvec) = psi(:,nvec+1:nvec+nvec)
        !
        ! ... refresh the reduced hamiltonian
        !
        nbase = nvec
        !
        ! Claudia
!       CALL DGEMM( 'T', 'N', nbase, nbase, nparm, 1.D0 , &
!                   psi, nparmx, hpsi, nparmx, 0.D0, hr, nvecx )
!       !
!       CALL DGEMM( 'T', 'N', nbase, nbase, nparm, 1.D0, &
!             psi, nparmx, spsi, nparmx, 0.D0, sr, nvecx )
!       CALL rdiaghg( nbase, nvec, hr, sr, nvecx, ew, vr )
!       !
!       e(1:nvec) = ew(1:nvec)
        ! Claudia
        !
        hr(:,1:nbase) = 0.D0
        sr(:,1:nbase) = 0.D0
        vr(:,1:nbase) = 0.D0
        !
        DO n = 1, nbase
           !
           hr(n,n) = e(n)
           sr(n,n) = 1.D0
           vr(n,n) = 1.D0
           !
        END DO
        !
        END IF ! idtask.eq.0
        ! 
        call MPI_BCAST(nbase,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        ! 
     END IF
     !
  END DO iterate
  !
  DEALLOCATE( conv )
  DEALLOCATE( ew )
  DEALLOCATE( vr )
  DEALLOCATE( hr )
  DEALLOCATE( sr )
  !
  DEALLOCATE( spsi )
  !
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )  
  !
  RETURN
  !
END SUBROUTINE regterg

!
!----------------------------------------------------------------------------
SUBROUTINE rdiaghg( n, m, h, s, ldh, e, v )
  !----------------------------------------------------------------------------
  ! ... Hv=eSv, with H symmetric matrix, S overlap matrix.
  ! ... On output both matrix are unchanged
  !
  ! ... LAPACK version - uses both DSYGV and DSYGVX
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: n, m, ldh
    ! dimension of the matrix to be diagonalized
    ! number of eigenstates to be calculated
    ! leading dimension of h, as declared in the calling pgm unit
  REAL*8, INTENT(INOUT) :: h(ldh,n), s(ldh,n)
    ! matrix to be diagonalized
    ! overlap matrix
  !
  REAL*8, INTENT(OUT) :: e(n)
    ! eigenvalues
  REAL*8, INTENT(OUT) :: v(ldh,m)
    ! eigenvectors (column-wise)
  !
  INTEGER               :: i, j, lwork, nb, mm, info
    ! mm = number of calculated eigenvectors
  REAL*8              :: abstol
  REAL*8, PARAMETER   :: one = 1.d0
  REAL*8, PARAMETER   :: zero = 0.d0
  INTEGER,  ALLOCATABLE :: iwork(:), ifail(:)
  REAL*8, ALLOCATABLE :: work(:), sdiag(:), hdiag(:)
  LOGICAL               :: all_eigenvalues
  INTEGER,  EXTERNAL    :: ILAENV
    ! ILAENV returns optimal block size "nb"
  !
  ! ... save the diagonal of input S (it will be overwritten)
  !
  ALLOCATE( sdiag( n ) )
  DO i = 1, n
     sdiag(i) = s(i,i)
  END DO
  !
  all_eigenvalues = ( m == n )
  !
  ! ... check for optimal block size
  !
  nb = ILAENV( 1, 'DSYTRD', 'U', n, -1, -1, -1 )
  !
  IF ( nb < 5 .OR. nb >= n ) THEN
     !
     lwork = 8*n
     !
  ELSE
     !
     lwork = ( nb + 3 )*n
     !
  END IF
  !
  ALLOCATE( work( lwork ) )
  !
  IF ( all_eigenvalues ) THEN
     !
     ! ... calculate all eigenvalues
     !
     v(:,:) = h(:,:)
     !
     CALL DSYGV( 1, 'V', 'U', n, v, ldh, s, ldh, e, work, lwork, info )
     !
  ELSE
     !
     ! ... calculate only m lowest eigenvalues
     !
     ALLOCATE( iwork( 5*n ) )
     ALLOCATE( ifail( n ) )
     !
     ! ... save the diagonal of input H (it will be overwritten)
     !
     ALLOCATE( hdiag( n ) )
     DO i = 1, n
        hdiag(i) = h(i,i)
     END DO
     !
     abstol = 0.D0
    ! abstol = 2.D0*DLAMCH( 'S' )
     !
     CALL DSYGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
                  0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
                  work, lwork, iwork, ifail, info )
     !
     DEALLOCATE( ifail )
     DEALLOCATE( iwork )
     !
     ! ... restore input H matrix from saved diagonal and lower triangle
     !
     DO i = 1, n
        h(i,i) = hdiag(i)
        DO j = i + 1, n
           h(i,j) = h(j,i)
        END DO
        DO j = n + 1, ldh
           h(j,i) = 0.0d0
        END DO
     END DO
     !
     DEALLOCATE( hdiag )
     !
  END IF
  !
  DEALLOCATE( work )
  !
  IF ( info > n ) THEN
!    CALL fatal_error( 'rdiaghg: S matrix not positive definite' )
     write(6,*) 'rdiaghg: S matrix not positive definite' 
  ELSE IF ( info > 0 ) THEN
     CALL fatal_error( 'rdiaghg: eigenvectors failed to converge' )
  ELSE IF ( info < 0 ) THEN
     CALL fatal_error( 'rdiaghg: incorrect call to DSYGV*' )
  END IF
  
  ! ... restore input S matrix from saved diagonal and lower triangle
  !
  DO i = 1, n
     s(i,i) = sdiag(i)
     DO j = i + 1, n
        s(i,j) = s(j,i)
     END DO
     DO j = n + 1, ldh
        s(j,i) = 0.0d0
     END DO
  END DO
  !
  DEALLOCATE( sdiag )
  !
  RETURN
  !
END SUBROUTINE rdiaghg

