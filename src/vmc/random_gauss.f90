SUBROUTINE gasdev_s(harvest)
! Numerical Recipes routine for generating a single normal random deviate,
! adapted to use the compiler's random number generator.

IMPLICIT NONE
integer, parameter :: dp = selected_real_kind(15, 307)
REAL(kind=dp), INTENT(OUT) :: harvest

! Local variables
REAL          :: rsq, v1, v2
REAL, SAVE    :: g
LOGICAL, SAVE :: gaus_stored = .false.

IF (gaus_stored) THEN
   harvest = g
   gaus_stored = .false.
ELSE
   DO
      CALL RANDOM_NUMBER(v1)
      CALL RANDOM_NUMBER(v2)
      v1 = 2.0*v1 - 1.0
      v2 = 2.0*v2 - 1.0
      rsq = v1**2 + v2**2
      if (rsq > 0.0 .and. rsq < 1.0) EXIT
   END DO
   rsq = SQRT(-2.0*LOG(rsq)/rsq)
   harvest = v1*rsq
   g = v2*rsq
   gaus_stored = .true.
END IF

RETURN
END SUBROUTINE gasdev_s
