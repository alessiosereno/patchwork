! Contiene (per il momento 12-04-2023):
!  * Interpolazione spline (con estrapolazione lineare out of range)
!  * Interpolazione lineare (con valori costanti out of range)
!  * Costanti PI, e, g0
!  * Linspace
!  * Matrice di rotazione 2D 
!  * Differenze finite 

  MODULE constants_mod
  IMPLICIT NONE

  ! Constants contain more digits than double precision, so that
  ! they are rounded correctly. Single letter constants contain underscore so
  ! that they do not clash with user variables ("e" and "i" are frequently used as
  ! loop variables)
  REAL*8,PARAMETER  :: PI = 3.1415926535897932384626433832795d0
  REAL*8,PARAMETER  :: e_ = 2.7182818284590452353602874713527d0
  REAL*8,PARAMETER  :: g0 = 9.80665d0

  END MODULE constants_mod

!********************************************************************!   

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
!  Come si usa:
!  Aggiungere "USE spline_mod" al main
!  Richiamare la funzione con YS=SPLINE(X,Y,XS)
!  La dimensione di YS e XS puo' essere qualsiasi
!  Fuori dall'intervallo estrapola linearmente
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

!=======================================================!
    FUNCTION SCALARSPLINE(X,Y,xs) RESULT(ys)
    IMPLICIT NONE
!   WRAPPER: COEFF +  SEVAL
!   X, Y : input data arrays size N
!   xs   : point to evaluate
!   ys   : interpolated point
    REAL*8,DIMENSION(:),INTENT(IN)   :: X,Y
    REAL*8,INTENT(IN)                :: xs
    INTEGER*4                        :: N
    REAL*8                           :: ys
    REAL*8,DIMENSION(:),ALLOCATABLE  :: B,C,D
    INTEGER*4                        :: I

    N = SIZE(X)

    ALLOCATE(B(N),C(N),D(N))

!   SPLINE COEFFICIENTS
    CALL COEFF(N,X,Y,B,C,D)
!   INTERPOLATION
    CALL SEVAL(N,XS,YS,X,Y,B,C,D)

    END FUNCTION
 END MODULE spline_mod

!********************************************************!

 MODULE interp1_mod
!  Come si usa:
!  Aggiungere "USE interp1_mod" al main
!  Richiamare la funzione con YS=interp1(X,Y,XS)
!  La dimensione di YS e XS puo' essere qualsiasi
!  Fuori dall'intervallo restituisce valori costanti
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
!   WRAPPER: COEFF +  SEVAL
!   X, Y : inpurprime = spline(x,r,xprime)
!   data arrays size N
!   xs   : points to evaluate size M
!   ys   : interpolated array size M
    REAL*8,DIMENSION(:),INTENT(IN)    :: X,Y
    REAL*8,DIMENSION(:),INTENT(IN)    :: xs
    INTEGER*4                         :: N,M
    REAL*8,DIMENSION(:),INTENT(INOUT) :: ys
    INTEGER*4                         :: I

    N = SIZE(X)
    M = SIZE(xs)

    !if (.not.allocated(ys)) ALLOCATE(ys(lbound(xs,1):ubound(xs,1)))

!   INTERPOLATION
    DO I=1,M
      CALL SEVAL1(XS(I),YS(I),X,Y)
    END DO

    END subroutine
 END MODULE INTERP1_MOD

!********************************************************************!   

 MODULE INTERP2_MOD
!  Come si usa:
!  Aggiungere "USE interp2_mod" al main
!  Richiamare la funzione con ZS=interp2(X,Y,XS,YS)
!  La dimensione di YS e XS puo' essere qualsiasi
!  Fuori dall'intervallo restituisce valori costanti
 IMPLICIT NONE
 CONTAINS

   SUBROUTINE BILIN(F1,F2,F3,F4,X1,X2,X3,X4,Y1,Y2,Y3,Y4,A,B,C,D)
   IMPLICIT NONE
   REAL*8   ::  F1,F2,F3,F4,X1,X2,X3,X4,Y1,Y2,Y3,Y4,A,B,C,D
   REAL*8   ::  Z1,Z2,Z3,Z4,DET31,DET32,DET33,DET34,DET4 
   REAL*8   ::  DET4A,DET4B,DET4C,DET4D 
!  LINEAR COEFFICIENT TO MAP F BETWEEN POINTS 1,2,3,4
!  F = A + B*X + C*Y + D*X*Y   
     Z1=X1*Y1
     Z2=X2*Y2
     Z3=X3*Y3
     Z4=X4*Y4

     DET31=X2*Y3*Z4+Y2*Z3*X4+Z2*X3*Y4-Z2*Y3*X4-Z3*Y4*X2-Z4*Y2*X3
     DET32=Y3*Z4+Y2*Z3+Z2*Y4-Z2*Y3-Z3*Y4-Z4*Y2
     DET33=X2*Z4+Z3*X4+Z2*X3-Z2*X4-Z3*X2-Z4*X3
     DET34=X2*Y3+Y2*X4+X3*Y4-Y3*X4-Y4*X2-Y2*X3
     DET4=DET31-X1*DET32-Y1*DET33-Z1*DET34

     DET31=X2*Y3*Z4+Y2*Z3*X4+Z2*X3*Y4-Z2*Y3*X4-Z3*Y4*X2-Z4*Y2*X3
     DET32=F2*Y3*Z4+Y2*Z3*F4+Z2*F3*Y4-Z2*Y3*F4-Z3*Y4*F2-Z4*Y2*F3
     DET33=X2*F3*Z4+F2*Z3*X4+Z2*X3*F4-Z2*F3*X4-Z3*F4*X2-Z4*F2*X3
     DET34=X2*Y3*F4+Y2*F3*X4+F2*X3*Y4-F2*Y3*X4-F3*Y4*X2-F4*Y2*X3
     DET4A=F1*DET31-X1*DET32-Y1*DET33-Z1*DET34

     DET31=F2*Y3*Z4+Y2*Z3*F4+Z2*F3*Y4-Z2*Y3*F4-Z3*Y4*F2-Z4*Y2*F3
     DET32=Y3*Z4+Y2*Z3+Z2*Y4-Z2*Y3-Z3*Y4-Z4*Y2
     DET33=F2*Z4+Z3*F4+Z2*F3-Z2*F4-Z3*F2-Z4*F3
     DET34=F2*Y3+Y2*F4+F3*Y4-Y3*F4-Y4*F2-Y2*F3
     DET4B=DET31-F1*DET32-Y1*DET33-Z1*DET34

     DET31=X2*F3*Z4+F2*Z3*X4+Z2*X3*F4-Z2*F3*X4-Z3*F4*X2-Z4*F2*X3
     DET32=F3*Z4+F2*Z3+Z2*F4-Z2*F3-Z3*F4-Z4*F2
     DET33=X2*Z4+Z3*X4+Z2*X3-Z2*X4-Z3*X2-Z4*X3
     DET34=X2*F3+F2*X4+X3*F4-F3*X4-F4*X2-F2*X3
     DET4C=DET31-X1*DET32-F1*DET33-Z1*DET34

     DET31=X2*Y3*F4+Y2*F3*X4+F2*X3*Y4-F2*Y3*X4-F3*Y4*X2-F4*Y2*X3
     DET32=Y3*F4+Y2*F3+F2*Y4-F2*Y3-F3*Y4-F4*Y2
     DET33=X2*F4+F3*X4+F2*X3-F2*X4-F3*X2-F4*X3
     DET34=X2*Y3+Y2*X4+X3*Y4-Y3*X4-Y4*X2-Y2*X3
     DET4D=DET31-X1*DET32-Y1*DET33-F1*DET34

     A=DET4A/DET4
     B=DET4B/DET4
     C=DET4C/DET4
     D=DET4D/DET4

     RETURN
   END

!=======================================================!

    FUNCTION INTERP2(X,Y,Z,XS,YS) RESULT(ZS)
    IMPLICIT NONE
!   WRAPPER
!   X, Y  : input coordinates, data arrays size NI,NJ
!   XS,YS : points to evaluate size MI,MJ
!   ZS    : interpolated array size MI,MJ
    REAL*8,DIMENSION(:,:),INTENT(IN)   :: X,Y,Z
    REAL*8,DIMENSION(:,:),ALLOCATABLE  :: DIST
    REAL*8,DIMENSION(:,:),INTENT(IN)   :: XS,YS
    INTEGER*4                          :: NI,NJ,MI,MJ
    REAL*8,DIMENSION(:,:),ALLOCATABLE  :: ZS
    INTEGER*4                          :: I,J,K,H,I1,I2,I3,I4,J1,J2,J3,J4
    REAL*8                             :: XP,YP,NUM,DEN
    REAL*8,DIMENSION(2)                :: DUM

    NI = SIZE(X,1)
    NJ = SIZE(X,2)
    MI = SIZE(XS,1)
    MJ = SIZE(XS,2)

    ALLOCATE(ZS(MI,MJ))

    DO I = 1,MI
      DO J = 1,MJ
!   XP,YP INTERPOLATION COORDINATES
        XP = XS(I,J)
        YP = YS(I,J)

!   DISTANCE OF EACH POINT OF DOMAIN FROM P
        DIST = SQRT( (X-XP)**2 + (Y-YP)**2 )

!   INVERSE DISTANCE WEIGHTED INTERPOLATION
!        NUM = 0.d0
!        DEN = 0.d0
!        DO K = 1,NI
!          DO H = 1,MI
!            NUM = NUM + 1.d0/DIST(K,H)**4 * Z(K,H)
!            DEN = DEN + 1.d0/DIST(K,H)**4
!          END DO
!        END DO
!        ZS(I,J) = NUM/DEN
!      END DO
!    END DO

!   FIND POINT OF MINIMUM DISTANCE
        dum = MINLOC(DIST)
        I1 = dum(1)
        J1 = dum(2)

!   FIND 2nd POINT
    J2 = J1
    IF ( ( X(I1+1,J1).GT.XP ) .AND. ( X(I1,J1).LT.XP ) ) THEN
      I2 = I1+1
    ELSE
      I2 = I1-1
    END IF

!   FIND 3rd POINT
    I3 = I1
    IF ( ( Y(I1+1,J1).GT.YP ) .AND. ( Y(I1,J1).LT.YP ) ) THEN
      J3 = J1+1
    ELSE
      J3 = J1-1
    END IF

!   FIND 4th POINT
    J4 = J3
    I4 = I2

!   RE-ARRANGE THE POINTS. 1: (I,J) ; 2: (I+1,J) ; 3: (I,J+1) ; 4: (I+1,J+1)
!    II1 = MINVAL(I1,I2) 

!   FUNCTION VALUES AT THE VERTICES
!    Z1 = Z(I1,J1)
!    Z2 = Z(I2,J2)
!    Z3 = Z(I3,J3)
!    Z4 = Z(I4,J4)
    END DO
    END DO

    END FUNCTION
 END MODULE INTERP2_MOD

!********************************************************************!   

 MODULE math_mod
 USE SPLINE_MOD
 USE LINSPACE_MOD
 USE CONSTANTS_MOD
 USE INTERP1_MOD
 USE INTERP2_MOD

 IMPLICIT NONE

 CONTAINS


 FUNCTION DERIVATIVE(X,Y) RESULT(DXDY)
 ! Finite differences DX/DY. 1st order at the boundary, 2nd internal points
   IMPLICIT NONE
 
   REAL*8,dimension(:),intent(IN) :: x, y
   REAL*8,dimension(size(x))      :: dxdy 
   INTEGER                        :: i, j, n 

   N = size(x) 

   do i = 1, n 
     if ( i .eq. 1) then
       dxdy(i) = (y(2)-y(1))/(x(2)-x(1))
     else if (i .eq. n) then
       dxdy(i) = (y(n)-y(n-1))/(x(n)-x(n-1))
     else
       dxdy(i) = (y(n+1)-y(n-1))/(x(n+1)-x(n-1))
     end if
   end do

   return 

 END FUNCTION DERIVATIVE

!********************************************************************!   

 FUNCTION RotMatrix(phi) result(RM)
  IMPLICIT NONE
  REAL*8,intent(in)      :: phi
  REAL*8,dimension(2,2)  :: RM

    RM(1,1) = COS(phi)
    RM(1,2) =-SIN(phi)
    RM(2,1) = SIN(phi)
    RM(2,2) = COS(phi)

 END FUNCTION RotMatrix

!********************************************************************!   

 FUNCTION CO_ATAN(Y,X) RESULT(A)
 ! ATAN2(Y,Z) but from 0 to 2*PI
   IMPLICIT NONE
 
   REAL*8,DIMENSION(:),INTENT(IN) :: X, Y
   REAL*8,DIMENSION(SIZE(X))      :: A
   INTEGER                        :: i, j, n 

   N = size(x) 
   A = ATAN2(Y,X)

   DO i = 1, N 
     IF ( A(i) .LE. 0.D0) THEN
       A(i) = A(i) + 2*PI
     END IF
   END DO

   RETURN 

 END FUNCTION CO_ATAN

!********************************************************************!   
 END MODULE math_mod

