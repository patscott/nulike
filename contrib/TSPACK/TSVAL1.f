      SUBROUTINE TSVAL1 (N,X,Y,YP,SIGMA,IFLAG,NE,TE, V,
     .                   IER)
      INTEGER N, IFLAG, NE, IER
      DOUBLE PRECISION X(N), Y(N), YP(N), SIGMA(N), TE(NE),
     .                 V(NE)
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
C   This subroutine evaluates a Hermite interpolatory ten-
C sion spline H or its first or second derivative at a set
C of points TE.
C
C   Note that a large tension factor in SIGMA may cause
C underflow.  The result is assumed to be zero.  If not the
C default, this may be specified by either a compiler option
C or operating system option.
C
C On input:
C
C       N = Number of data points.  N .GE. 2.
C
C       X = Array of length N containing the abscissae.
C           These must be in strictly increasing order:
C           X(I) < X(I+1) for I = 1,...,N-1.
C
C       Y = Array of length N containing data values or
C           function values returned by Subroutine SMCRV.
C           Y(I) = H(X(I)) for I = 1,...,N.
C
C       YP = Array of length N containing first deriva-
C            tives.  YP(I) = HP(X(I)) for I = 1,...,N, where
C            HP denotes the derivative of H.
C
C       SIGMA = Array of length N-1 containing tension fac-
C               tors whose absolute values determine the
C               balance between cubic and linear in each
C               interval.  SIGMA(I) is associated with int-
C               erval (I,I+1) for I = 1,...,N-1.
C
C       IFLAG = Output option indicator:
C               IFLAG = 0 if values of H are to be computed.
C               IFLAG = 1 if first derivative values are to
C                         be computed.
C               IFLAG = 2 if second derivative values are to
C                         be computed.
C
C       NE = Number of evaluation points.  NE > 0.
C
C       TE = Array of length NE containing the evaluation
C            points.  The sequence should be strictly in-
C            creasing for maximum efficiency.  Extrapolation
C            is performed if a point is not in the interval
C            [X(1),X(N)].
C
C The above parameters are not altered by this routine.
C
C       V = Array of length at least NE.
C
C On output:
C
C       V = Array of function, first derivative, or second
C           derivative values at the evaluation points un-
C           less IER < 0.  If IER = -1, V is not altered.
C           If IER = -2, V may be only partially defined.
C
C       IER = Error indicator:
C             IER = 0  if no errors were encountered and
C                      no extrapolation occurred.
C             IER > 0  if no errors were encountered but
C                      extrapolation was required at IER
C                      points.
C             IER = -1 if N < 2, IFLAG < 0, IFLAG > 2, or
C                      NE < 1.
C             IER = -2 if the abscissae are not in strictly
C                      increasing order.  (This error will
C                      not necessarily be detected.)
C
C Modules required by TSVAL1:  HPPVAL, HPVAL, HVAL, INTRVL,
C                                SNHCSH
C
C***********************************************************
C
      INTEGER I, IERR, IFLG, NVAL, NX
      DOUBLE PRECISION HPPVAL, HPVAL, HVAL
C
      IFLG = IFLAG
      NVAL = NE
C
C Test for invalid input.
C
      IF (N .LT. 2  .OR.  IFLG .LT. 0  .OR.  IFLG .GT. 2
     .    .OR.  NVAL .LT. 1) GO TO 2
C
C Initialize the number of extrapolation points NX and
C   loop on evaluation points.
C
      NX = 0
      DO 1 I = 1,NVAL
        IF (IFLG .EQ. 0) THEN
          V(I) = HVAL (TE(I),N,X,Y,YP,SIGMA, IERR)
        ELSEIF (IFLG .EQ. 1) THEN
          V(I) = HPVAL (TE(I),N,X,Y,YP,SIGMA, IERR)
        ELSE
          V(I) = HPPVAL (TE(I),N,X,Y,YP,SIGMA, IERR)
        ENDIF
        IF (IERR .LT. 0) GO TO 3
        NX = NX + IERR
    1   CONTINUE
C
C No errors encountered.
C
      IER = NX
      RETURN
C
C N, IFLAG, or NE is outside its valid range.
C
    2 IER = -1
      RETURN
C
C X is not strictly increasing.
C
    3 IER = -2
      RETURN
      END
