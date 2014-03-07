      SUBROUTINE ARCL2D (N,X,Y, T,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), T(N)
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
C   Given an ordered sequence of points (X,Y) defining a
C polygonal curve in the plane, this subroutine computes the
C sequence T of cumulative arc lengths along the curve:
C T(1) = 0 and, for 2 .LE. K .LE. N, T(K) is the sum of
C Euclidean distances between (X(I-1),Y(I-1)) and (X(I),Y(I))
C for I = 2,...,K.  A closed curve corresponds to X(1) =
C X(N) and Y(1) = Y(N), and more generally, duplicate points
C are permitted but must not be adjacent.  Thus, T contains
C a strictly increasing sequence of values which may be used
C as parameters for fitting a smooth curve to the sequence
C of points.
C
C On input:
C
C       N = Number of points defining the curve.  N .GE. 2.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the points.
C
C The above parameters are not altered by this routine.
C
C       T = Array of length at least N.
C
C On output:
C
C       T = Array containing cumulative arc lengths defined
C           above unless IER > 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N < 2.
C             IER = I if X(I-1) = X(I) and Y(I-1) = Y(I) for
C                     some I in the range 2,...,N.
C
C Modules required by ARCL2D:  None
C
C Intrinsic function called by ARCL2D:  SQRT
C
C***********************************************************
C
      INTEGER I, NN
      DOUBLE PRECISION DS
C
      NN = N
      IF (NN .LT. 2) GO TO 2
      T(1) = 0.D0
      DO 1 I = 2,NN
        DS = (X(I)-X(I-1))**2 + (Y(I)-Y(I-1))**2
        IF (DS .EQ. 0.D0) GO TO 3
        T(I) = T(I-1) + SQRT(DS)
    1   CONTINUE
      IER = 0
      RETURN
C
C N is outside its valid range.
C
    2 IER = 1
      RETURN
C
C Points I-1 and I coincide.
C
    3 IER = I
      RETURN
      END
