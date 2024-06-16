 MODULE io_mod
   USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
   IMPLICIT NONE

 CONTAINS
 !------------------------------------------------------------------- 
   SUBROUTINE IMPORTDATA(FILENAME,X,Y)
     CHARACTER(100),INTENT(IN)             :: FILENAME
     INTEGER(KIND=4)                       :: FID    
     REAL(DP),ALLOCATABLE,INTENT(INOUT)    :: X(:),Y(:)
     REAL(DP),DIMENSION(100000)            :: xtemp,ytemp 
     INTEGER(KIND=4)                       :: I

     OPEN(NEWUNIT=FID,FILE=FILENAME,STATUS='OLD')
     DO I = 1, 100000
       READ(FID,*,END=100) xtemp(I),ytemp(I)
     END DO

100  CONTINUE

     ALLOCATE(X(I-1),Y(I-1))

     X = xtemp(1:I-1)
     Y = ytemp(1:I-1)

     CLOSE(FID)
   END SUBROUTINE IMPORTDATA
 !------------------------------------------------------------------- 
 !------------------------------------------------------------------- 
   SUBROUTINE READARRAY(FILENAME,X)
     CHARACTER(100),INTENT(IN)             :: FILENAME
     INTEGER(KIND=4)                       :: FID    
     REAL(DP),ALLOCATABLE,INTENT(INOUT)    :: X(:)
     REAL(DP),DIMENSION(100000)            :: xtemp
     INTEGER(KIND=4)                       :: I

     OPEN(NEWUNIT=FID,FILE=FILENAME,STATUS='OLD')
     DO I = 1, 100000
       READ(FID,*,END=100) xtemp(I)
     END DO

100  CONTINUE

     ALLOCATE(X(I-1))

     X = xtemp(1:I-1)

     CLOSE(FID)
   END SUBROUTINE READARRAY
 !------------------------------------------------------------------- 
 !------------------------------------------------------------------- 
   SUBROUTINE WRITEARRAY(FILENAME,X)
     CHARACTER(100),INTENT(IN)             :: FILENAME
     INTEGER(KIND=4)                       :: FID    
     REAL(DP),DIMENSION(:),INTENT(IN)      :: X(:)
     INTEGER(KIND=4)                       :: I, N

     N = SIZE(X)

     OPEN(NEWUNIT=FID,FILE=FILENAME,STATUS='UNKNOWN')

     DO I = 1, N
       WRITE(FID,*) X(I)
     END DO

     CLOSE(FID)
   END SUBROUTINE WRITEARRAY
 !------------------------------------------------------------------- 
 !------------------------------------------------------------------- 
   SUBROUTINE FLIPARRAY(X)
     REAL(DP),DIMENSION(:),INTENT(INOUT)   :: X(:)
     INTEGER(KIND=4)                       :: I, N
     REAL(DP),ALLOCATABLE                  :: XFLIP(:)

     N = SIZE(X)
     ALLOCATE(XFLIP(N))

     DO I = 1, N
       XFLIP(I) = X(N-I+1)
     END DO
 
     X = XFLIP
   
   END SUBROUTINE FLIPARRAY
 !------------------------------------------------------------------- 

 END MODULE io_mod

