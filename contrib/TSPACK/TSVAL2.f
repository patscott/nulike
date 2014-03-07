      SUBROUTINE TSVAL2 (N,T,X,Y,XP,YP,SIGMA,IFLAG,NE,
     .                   TE, VX,VY,IER)
      INTEGER N, IFLAG, NE, IER
      DOUBLE PRECISION T(N), X(N), Y(N), XP(N), YP(N),
     .                 SIGMA(N), TE(NE), VX(NE), VY(NE)
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
C   This subroutine returns values or derivatives of a pair
C of Hermite interpolatory tension splines H1 and H2 which
C form the components of a parametric planar curve C(t) =
C (H1(t),H2(t)).  Refer to Subroutines TSPBP and TSPSP.
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
C       T = Array of length N containing a strictly in-
C           creasing sequence of abscissae (parameter
C           values).  Refer to Subroutine ARCL2D.
C
C       X = Array of length N containing data values or
C           function values returned by Subroutine SMCRV.
C           X(I) = H1(T(I)) for I = 1,...,N.
C
C       Y = Array of length N containing data values or
C           function values returned by Subroutine SMCRV.
C           Y(I) = H2(T(I)) for I = 1,...,N.
C
C       XP = Array of length N containing first deriva-
C            tives.  XP(I) = H1P(T(I)) for I = 1,...,N,
C            where H1P denotes the derivative of H1.
C
C       YP = Array of length N containing first deriva-
C            tives.  YP(I) = H2P(T(I)) for I = 1,...,N,
C            where H2P denotes the derivative of H2.
C
C   Note that C(T(I)) = (X(I),Y(I)) and CP(T(I)) = (XP(I),
C YP(I)), I = 1,...,N, are data (control) points and deriva-
C tive (velocity) vectors, respectively.
C
C       SIGMA = Array of length N-1 containing tension fac-
C               tors whose absolute values determine the
C               balance between cubic and linear in each
C               interval.  SIGMA(I) is associated with int-
C               erval (I,I+1) for I = 1,...,N-1.
C
C       IFLAG = Output option indicator:
C               IFLAG = 0 if values of H1 and H2 (points on
C                         the curve) are to be computed.
C               IFLAG = 1 if first derivative vectors are to
C                         be computed.  Unit tangent vectors
C                         can be obtained by normalizing
C                         these to unit vectors.
C               IFLAG = 2 if second derivative (accelera-
C                         tion) vectors are to be computed.
C                         Given a unit tangent vector U and
C                         a second derivative vector V, the
C                         corresponding curvature vector
C                         can be computed as the cross
C                         product U X V X U.
C
C       NE = Number of evaluation points.  NE > 0.
C
C       TE = Array of length NE containing the evaluation
C            points.  The sequence should be strictly in-
C            creasing for maximum efficiency.  Extrapolation
C            is performed if a point is not in the interval
C            [T(1),T(N)].
C
C The above parameters are not altered by this routine.
C
C       VX,VY = Arrays of length at least NE.
C
C On output:
C
C       VX,VY = Arrays containing values, first derivatives,
C               or second derivatives of H1 and H2, respec-
C               tively, at the evaluation points (unless
C               IER < 0).  If IER = -1, VX and VY are not
C               altered.  If IER = -2, VX and VY may be only
C               partially defined.
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
C Modules required by TSVAL2:  HPPVAL, HPVAL, HVAL, INTRVL,
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
          VX(I) = HVAL (TE(I),N,T,X,XP,SIGMA, IERR)
          VY(I) = HVAL (TE(I),N,T,Y,YP,SIGMA, IERR)
        ELSEIF (IFLG .EQ. 1) THEN
          VX(I) = HPVAL (TE(I),N,T,X,XP,SIGMA, IERR)
          VY(I) = HPVAL (TE(I),N,T,Y,YP,SIGMA, IERR)
        ELSE
          VX(I) = HPPVAL (TE(I),N,T,X,XP,SIGMA, IERR)
          VY(I) = HPPVAL (TE(I),N,T,Y,YP,SIGMA, IERR)
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
C T is not strictly increasing.
C
    3 IER = -2
      RETURN
      END
