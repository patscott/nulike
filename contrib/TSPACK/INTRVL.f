      INTEGER FUNCTION INTRVL (T,N,X)
      INTEGER N
      DOUBLE PRECISION T, X(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/10/92
C
C   This function returns the index of the left end of an
C interval (defined by an increasing sequence X) which
C contains the value T.  The method consists of first test-
C ing the interval returned by a previous call, if any, and
C then using a binary search if necessary.
C
C On input:
C
C       T = Point to be located.
C
C       N = Length of X.  N .GE. 2.
C
C       X = Array of length N assumed (without a test) to
C           contain a strictly increasing sequence of
C           values.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       INTRVL = Index I defined as follows:
C
C                  I = 1    if  T .LT. X(2) or N .LE. 2,
C                  I = N-1  if  T .GE. X(N-1), and
C                  X(I) .LE. T .LT. X(I+1) otherwise.
C
C Modules required by INTRVL:  None
C
C***********************************************************
C
      INTEGER IH, IL, K
      DOUBLE PRECISION TT
C
C     Commented out for thread safety by Pat Scott Jan 6 2016
C      SAVE IL
C      DATA IL/1/
      TT = T
C     Commented out for thread safety by Pat Scott Jan 6 2016
C      IF (IL .GE. 1  .AND.  IL .LT. N) THEN
C        IF (X(IL) .LE. TT  .AND.  TT .LT. X(IL+1)) GO TO 2
C      ENDIF
C
C Initialize low and high indexes.
C
      IL = 1
      IH = N
C
C Binary search:
C
    1 IF (IH .LE. IL+1) GO TO 2
        K = (IL+IH)/2
        IF (TT .LT. X(K)) THEN
          IH = K
        ELSE
          IL = K
        ENDIF
        GO TO 1
C
C X(IL) .LE. T .LT. X(IL+1)  or  (T .LT. X(1) and IL=1)
C                            or  (T .GE. X(N) and IL=N-1)
C
    2 INTRVL = IL
      RETURN
      END
