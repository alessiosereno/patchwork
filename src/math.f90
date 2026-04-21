! Spline (cubic, linear-extrapolated out of range), linear interp1, linspace.

 MODULE linspace_mod
   USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
   IMPLICIT NONE

 CONTAINS
  !   Return evenly spaced numbers over a specified interval.
  !   Returns `num` evenly spaced samples, calculated over the interval `[start, stop]`.
  !   Ported from the numpy routine.
  !   Author: Ivan Pribec
  !
  FUNCTION LINSPACE(START,END,NUM,ENDPOINT,STEP) RESULT(SAMPLES)

  ! PARAMETERS
    real(kind=8), INTENT(IN) :: START
      !! The starting value of the sequence.
    real(kind=8), INTENT(IN) :: END
      !! The end value of the sequence, unless `endpoint` is set to `.false.`.
      !! In that case, the sequence consists of all but the last of `num + 1`
      !! evenly spaced samples, so that `end` is excluded. Note that the
      !! step size changes when `endpoint` is `.false.`.
    INTEGER, INTENT(IN) :: NUM
      !! Number of samples to generate. Default value is 50.
    LOGICAL, INTENT(IN), OPTIONAL :: ENDPOINT
      !! If `.true.`, `end` is the last sample. Otherwise, it is not included. Default is `.true.`.
    real(kind=8), INTENT(OUT), OPTIONAL :: STEP
      !! If present, `step` is the size of spacing between samples.

  ! RETURNS
    real(kind=8) :: SAMPLES(NUM)
      !! There are `num` equally spaced samples in the closed interval `[start, stop]` or
      !! the half-open interval `[start, stop)` (depending on whether `endpoint` is `.true.` or `.false.`).

    INTEGER :: NUM_, I
    LOGICAL :: ENDPOINT_
    real(kind=8) :: STEP_

        NUM_ = NUM

        ENDPOINT_ = .TRUE.
        IF (PRESENT(ENDPOINT)) ENDPOINT_ = ENDPOINT

        ! FIND STEP SIZE
        IF (ENDPOINT_) THEN
            STEP_ = (END - START)/REAL(NUM_-1,DP)
        ELSE
            STEP_ = (END - START)/REAL(NUM_,DP)
        END IF

        IF (PRESENT(STEP)) STEP = STEP_

        DO I = 1, NUM_
            SAMPLES(I) = START + (I-1)*STEP_
        END DO
    END FUNCTION LINSPACE
 END MODULE linspace_mod

!********************************************************************!

 MODULE spline_mod
!  Usage: USE spline_mod and call YS=SPLINE(X,Y,XS)
!  XS and YS can be of any size. Linearly extrapolated outside range.
 IMPLICIT NONE
 CONTAINS
   SUBROUTINE SEVAL(N,U,V,X,Y,B,C,D)
!------------------------------------------------------------------------
!  EVALUATE A CUBIC SPLINE INTERPOLATION OF A DISCRETE FUNCTION F(X),
!  GIVEN IN N POINTS X(I), Y(I). THE B, C AND D COEFFICIENTS DEFINING
!  THE BEST CUBIC SPLINE FOR THE GIVEN POINTS, ARE CALCULATED BEFORE
!  BY THE SPLINE SUBROUTINE.
!
!  INPUTS:
!    N       NUMBER OF POINTS OF CURVE Y = F(X)
!    U       ABSCISSA OF POINT TO BE INTERPOLATED
!    X,Y     TABLES OF DIMENSION N, STORING THE COORDINATES
!            OF CURVE F(X)
!    B,C,D   TABLES STORING THE COEFFICIENTS DEFINING THE
!            CUBIC SPLINE
!
!   OUTPUTS:
!    SEVAL   INTERPOLATED VALUE
!            = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
!            WITH DX = U-X(I), U BETWEEN X(I) AND X(I+1)
!
!   REFERENCE :
!    FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!    COMPUTATIONS. PRENTICE-HALL,INC.
!------------------------------------------------------------------------
     IMPLICIT NONE
     REAL*8    :: B(N),C(N),D(N),X(N),Y(N),U,V,DX
     INTEGER*4 :: I,J,K,N

     I = 1

!    OUT OF RANGE
     IF (U.LT.X(I)) THEN
       V = Y(1)+(U-X(1))*(Y(2)-Y(1))/(X(2)-X(1))
       GOTO 99
     ELSEIF (U.GT.X(N)) THEN
       V = Y(N)+(U-X(N))*(Y(N)-Y(N-1))/(X(N)-X(N-1))
       GOTO 99
     ENDIF

!    BINARY SEARCH
     IF (I.GE.N) I = 1
     IF (U.LT.X(I)) GO TO 10
     IF (U.LE.X(I+1)) GO TO 30
  10 I = 1
     J = N+1
  20 K = (I+J)/2
     IF (U.LT.X(K)) J = K
     IF (U.GE.X(K)) I = K
     IF (J.GT.I+1) GO TO 20

!    SPLINE EVALUATION

  30 DX = U-X(I)
     V = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))

  99 CONTINUE
     RETURN
     END

   SUBROUTINE COEFF(N,X,Y,B,C,D)
!---------------------------------------------------------------------
!  THIS SUBROUTINE CALCULATES THE COEFFICIENTS B,C,D OF A CUBIC
!  SPLINE TO BEST APPROXIMATE A DISCRETE FUNCTION GIVEN BY N POINTS
!
!  INPUTS:
!   N       NUMBER OF GIVEN POINTS
!   X,Y     VECTORS OF DIMENSION N, STORING THE COORDINATES
!             OF FUNCTION F(X)
!
!  OUTPUTS:
!   A,B,C   VECTORS OF DIMENSION N, STORING THE COEFFICIENTS
!             OF THE CUBIC SPLINE
!
!  REFERENCE:
!   FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!   COMPUTATIONS. PRENTICE-HALL,INC.
!---------------------------------------------------------------------
   IMPLICIT NONE
   REAL*8    :: B(N),C(N),D(N),X(N),Y(N)
   REAL*8    :: T
   INTEGER*4 :: N,NM1,I,L

      NM1 = N-1
      IF (N.LT.2) RETURN
      IF (N.LT.3) GO TO 50

!     BUILD THE TRIDIAGONAL SYSTEM
!     B (DIAGONAL), D (UPPERDIAGONAL) , C (SECOND MEMBER)

      D(1) = X(2)-X(1)
      C(2) = (Y(2)-Y(1))/D(1)
      DO 10 I = 2,NM1
      D(I) = X(I+1)-X(I)
      B(I) = 2.D0*(D(I-1)+D(I))
      C(I+1) = (Y(I+1)-Y(I))/D(I)
      C(I) = C(I+1)-C(I)
   10 CONTINUE

!     CONDITIONS AT LIMITS
!     THIRD DERIVATIVES OBTAINED BY DIVIDED DIFFERENCES

      B(1) = -D(1)
      B(N) = -D(N-1)
      C(1) = 0.D0
      C(N) = 0.D0
      IF (N.EQ.3) GO TO 15
      C(1) = C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
      C(N) = C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
      C(1) = C(1)*D(1)*D(1)/(X(4)-X(1))
      C(N) = -C(N)*D(N-1)**2/(X(N)-X(N-3))

!     FORWARD ELIMINATION

   15 DO 20 I = 2,N
      T = D(I-1)/B(I-1)
      B(I) = B(I)-T*D(I-1)
      C(I) = C(I)-T*C(I-1)
   20 CONTINUE

!     BACK SUBSTITUTION

      C(N) = C(N)/B(N)
      DO 30 L = 1,NM1
      I = N-L
      C(I) = (C(I)-D(I)*C(I+1))/B(I)
   30 CONTINUE

!     COEFFICIENTS OF 3RD DEGREE POLYNOMIAL

      B(N) = (Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2.D0*C(N))
      DO 40 I = 1,NM1
      B(I) = (Y(I+1)-Y(I))/D(I)-D(I)*(C(I+1)+2.D0*C(I))
      D(I) = (C(I+1)-C(I))/D(I)
      C(I) = 3.D0*C(I)
   40 CONTINUE
      C(N) = 3.D0*C(N)
      D(N) = D(NM1)
      RETURN

!     CAS N = 2

   50 B(1) = (Y(2)-Y(1))/(X(2)-X(1))
      C(1) = 0.D0
      D(1) = 0.D0
      B(2) = B(1)
      C(2) = 0.D0
      D(2) = 0.D0
      RETURN
      END

!=======================================================!
    FUNCTION SPLINE(X,Y,xs) RESULT(ys)
    IMPLICIT NONE
!   WRAPPER: COEFF +  SEVAL
!   X, Y : input data arrays size N
!   xs   : points to evaluate size M
!   ys   : interpolated array size M
    REAL*8,DIMENSION(:),INTENT(IN)   :: X,Y
    REAL*8,DIMENSION(:),INTENT(IN)   :: xs
    INTEGER*4                        :: N,M
    REAL*8,DIMENSION(:),ALLOCATABLE  :: ys
    REAL*8,DIMENSION(:),ALLOCATABLE  :: B,C,D
    INTEGER*4                        :: I

    N = SIZE(X)
    M = SIZE(xs)

    ALLOCATE(B(N),C(N),D(N))
    ALLOCATE(ys(M))

!   SPLINE COEFFICIENTS
    CALL COEFF(N,X,Y,B,C,D)
!   INTERPOLATION
    DO I=1,M
      CALL SEVAL(N,XS(I),YS(I),X,Y,B,C,D)
    END DO

    END FUNCTION
 END MODULE spline_mod

!********************************************************!

 MODULE interp1_mod
!  Usage: USE interp1_mod and call interp1(X,Y,XS,YS).
!  XS and YS can be of any size. Constant-extrapolated outside range.
 IMPLICIT NONE
 CONTAINS
   SUBROUTINE SEVAL1(U,V,X,Y)
     IMPLICIT NONE
     REAL*8, DIMENSION(:), INTENT(IN) :: X, Y
     REAL*8    :: U,V,DX
     INTEGER*4 :: I,J,K,N

     I = 1
     N = SIZE(Y)

!    OUT OF RANGE
     IF (U.LT.X(I)) THEN
       V = Y(1)
       GOTO 99
     ELSEIF (U.GT.X(N)) THEN
       V = Y(N)
       GOTO 99
     ENDIF

!    BINARY SEARCH
     IF (I.GE.N) I = 1
     IF (U.LT.X(I)) GO TO 10
     IF (U.LE.X(I+1)) GO TO 30
  10 I = 1
     J = N+1
  20 K = (I+J)/2
     IF (U.LT.X(K)) J = K
     IF (U.GE.X(K)) I = K
     IF (J.GT.I+1) GO TO 20

!    LINEAR INTERPOLATION

  30 DX = U-X(I)
     V = Y(I)+DX*(Y(I+1)-Y(I))/(X(I+1)-X(I))

  99 CONTINUE
     RETURN
     END

!=======================================================!
    subroutine INTERP1(X,Y,xs,ys)
    IMPLICIT NONE
!   WRAPPER over SEVAL1
!   X, Y : input data arrays size N
!   xs   : points to evaluate size M
!   ys   : interpolated array size M
    REAL*8,DIMENSION(:),INTENT(IN)    :: X,Y
    REAL*8,DIMENSION(:),INTENT(IN)    :: xs
    INTEGER*4                         :: N,M
    REAL*8,DIMENSION(:),INTENT(INOUT) :: ys
    INTEGER*4                         :: I

    N = SIZE(X)
    M = SIZE(xs)

    DO I=1,M
      CALL SEVAL1(XS(I),YS(I),X,Y)
    END DO

    END subroutine
 END MODULE INTERP1_MOD
