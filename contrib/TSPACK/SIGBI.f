      SUBROUTINE SIGBI (N,X,Y,YP,TOL,B,BMAX, SIGMA, ICFLG,
     .                  DSMAX,IER)
      INTEGER N, ICFLG(N), IER
      DOUBLE PRECISION X(N), Y(N), YP(N), TOL, B(5,N), BMAX,
     .                 SIGMA(N), DSMAX
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
C   Given a set of abscissae X with associated data values Y
C and derivatives YP, this subroutine determines the small-
C est (nonnegative) tension factors SIGMA such that the Her-
C mite interpolatory tension spline H(x) satisfies a set of
C user-specified constraints.
C
C   SIGBI may be used in conjunction with Subroutine YPC2
C (or YPC2P) in order to produce a C-2 interpolant which
C satisfies the constraints.  This is achieved by calling
C YPC2 with SIGMA initialized to the zero vector, and then
C alternating calls to SIGBI with calls to YPC2 until the
C change in SIGMA is small (refer to the parameter descrip-
C tions for SIGMA, DSMAX and IER), or the maximum relative
C change in YP is bounded by a tolerance (a reasonable value
C is .01).  A similar procedure may be used to produce a C-2
C shape-preserving smoothing curve (Subroutine SMCRV).
C
C   Refer to Subroutine SIGS for a means of selecting mini-
C mum tension factors to preserve shape properties of the
C data.
C
C On input:
C
C       N = Number of data points.  N .GE. 2.
C
C       X = Array of length N containing a strictly in-
C           creasing sequence of abscissae:  X(I) < X(I+1)
C           for I = 1,...,N-1.
C
C       Y = Array of length N containing data values (or
C           function values computed by SMCRV) associated
C           with the abscissae.  H(X(I)) = Y(I) for I =
C           1,...,N.
C
C       YP = Array of length N containing first derivatives
C            of H at the abscissae.  Refer to Subroutines
C            YPC1, YPC1P, YPC2, YPC2P, and SMCRV.
C
C       TOL = Tolerance whose magnitude determines how close
C             each tension factor is to its optimal value
C             when nonzero finite tension is necessary and
C             sufficient to satisfy a constraint.  Refer to
C             functions SIG0, SIG1, and SIG2.  TOL should be
C             set to 0 for optimal tension.
C
C       B = Array dimensioned 5 by N-1 containing bounds or
C           flags which define the constraints.  For I = 1
C           to N-1, column I defines the constraints associ-
C           ated with interval I (X(I),X(I+1)) as follows:
C
C             B(1,I) is an upper bound on H
C             B(2,I) is a lower bound on H
C             B(3,I) is an upper bound on HP
C             B(4,I) is a lower bound on HP
C             B(5,I) specifies the required sign of HPP
C
C           where HP and HPP denote the first and second
C           derivatives of H, respectively.  A null con-
C           straint is specified by abs(B(K,I)) .GE. BMAX
C           for K < 5, or B(5,I) = 0:   B(1,I) .GE. BMAX,
C           B(2,I) .LE. -BMAX, B(3,I) .GE. BMAX, B(4,I) .LE.
C           -BMAX, or B(5,I) = 0.  Any positive value of
C           B(5,I) specifies that H should be convex, a
C           negative values specifies that H should be con-
C           cave, and 0 specifies that no restriction be
C           placed on HPP.  Refer to Functions SIG0, SIG1,
C           and SIG2 for definitions of valid constraints.
C
C       BMAX = User-defined value of infinity which, when
C              used as an upper bound in B (or when its
C              negative is used as a lower bound), specifies
C              that no constraint is to be enforced.
C
C The above parameters are not altered by this routine.
C
C       SIGMA = Array of length N-1 containing minimum val-
C               ues of the tension factors.  SIGMA(I) is as-
C               sociated with interval (I,I+1) and SIGMA(I)
C               .GE. 0 for I = 1,...,N-1.  SIGMA should be
C               set to the zero vector if minimal tension
C               is desired, and should be unchanged from a
C               previous call in order to ensure convergence
C               of the C-2 iterative procedure.
C
C       ICFLG = Array of length .GE. N-1.
C
C On output:
C
C       SIGMA = Array containing tension factors for which
C               H(x) satisfies the constraints defined by B,
C               with the restriction that SIGMA(I) .LE. 85
C               for all I (unless the input value is larger).
C               The factors are as small as possible (within
C               the tolerance), but not less than their
C               input values.  If infinite tension is re-
C               quired in interval (X(I),X(I+1)), then
C               SIGMA(I) = 85 (and H is an approximation to
C               the linear interpolant on the interval),
C               and if no constraint is specified in the
C               interval, then SIGMA(I) = 0 (unless the
C               input value is positive), and thus H is
C               cubic.  Invalid constraints are treated as
C               null constraints.
C
C       ICFLG = Array of invalid constraint flags associated
C               with intervals.  For I = 1 to N-1, ICFLG(I)
C               is a 5-bit value b5b4b3b2b1, where bK = 1 if
C               and only if constraint K cannot be satis-
C               fied.  Thus, all constraints in interval I
C               are satisfied if and only if ICFLG(I) = 0
C               (and IER .GE. 0).
C
C       DSMAX = Maximum increase in a component of SIGMA
C               from its input value.  The increase is a
C               relative change if the input value is
C               positive, and an absolute change otherwise.
C
C       IER = Error indicator and information flag:
C             IER = I if no errors (other than invalid con-
C                     straints) were encountered and I
C                     components of SIGMA were altered from
C                     their input values for 0 .LE. I .LE.
C                     N-1.
C             IER = -1 if N < 2.  SIGMA and ICFLG are not
C                      altered in this case.
C             IER = -I if X(I) .LE. X(I-1) for some I in the
C                      range 2,...,N.  SIGMA(J) and ICFLG(J)
C                      are unaltered for J .GE. I-1 in this
C                      case.
C
C Modules required by SIGBI:  SIG0, SIG1, SIG2, SNHCSH,
C                               STORE
C
C Intrinsic functions called by SIGBI:  ABS, MAX, MIN
C
C***********************************************************
C
      INTEGER I, ICFK, ICNT, IERR, IFL, K, NM1
      DOUBLE PRECISION BMX, BND, DSIG, DSM, S, SBIG, SIG,
     .                 SIGIN
      DOUBLE PRECISION SIG0, SIG1, SIG2
C
      DATA SBIG/85.D0/
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 4
      BMX = BMAX
C
C Initialize change counter ICNT and maximum change DSM for
C   loop on intervals.
C
      ICNT = 0
      DSM = 0.D0
      DO 3 I = 1,NM1
        IF (X(I) .GE. X(I+1)) GO TO 5
        ICFLG(I) = 0
C
C Loop on constraints for interval I.  SIG is set to the
C   largest tension factor required to satisfy all five
C   constraints.  ICFK = 2**(K-1) is the increment for
C   ICFLG(I) when constraint K is invalid.
C
        SIG = 0.D0
        ICFK = 1
        DO 2 K = 1,5
          BND = B(K,I)
          IF (K .LT. 5  .AND.  ABS(BND) .GE. BMX) GO TO 1
          IF (K .LE. 2) THEN
            IFL = 3 - 2*K
            S = SIG0 (X(I),X(I+1),Y(I),Y(I+1),YP(I),YP(I+1),
     .                IFL,BND,TOL, IERR)
          ELSEIF (K .LE. 4) THEN
            IFL = 7 - 2*K
            S = SIG1 (X(I),X(I+1),Y(I),Y(I+1),YP(I),YP(I+1),
     .                IFL,BND,TOL, IERR)
          ELSE
            IF (BND .EQ. 0.D0) GO TO 1
            IFL = -1
            IF (BND .GT. 0.D0) IFL = 1
            S = SIG2 (X(I),X(I+1),Y(I),Y(I+1),YP(I),YP(I+1),
     .                IFL,TOL, IERR)
          ENDIF
          IF (IERR .EQ. -2) THEN
C
C   An invalid constraint was encountered.  Increment
C     ICFLG(I).
C
            ICFLG(I) = ICFLG(I) + ICFK
          ELSE
C
C   Update SIG.
C
            SIG = MAX(SIG,S)
          ENDIF
C
C   Bottom of loop on constraints K:  update ICFK.
C
    1     ICFK = 2*ICFK
    2     CONTINUE
C
C Bottom of loop on intervals:  update SIGMA(I), ICNT, and
C   DSM if necessary.
C
        SIG = MIN(SIG,SBIG)
        SIGIN = SIGMA(I)
        IF (SIG .GT. SIGIN) THEN
          SIGMA(I) = SIG
          ICNT = ICNT + 1
          DSIG = SIG-SIGIN
          IF (SIGIN .GT. 0.D0) DSIG = DSIG/SIGIN
          DSM = MAX(DSM,DSIG)
        ENDIF
    3   CONTINUE
C
C No errors (other than invalid constraints) encountered.
C
      DSMAX = DSM
      IER = ICNT
      RETURN
C
C N < 2.
C
    4 DSMAX = 0.D0
      IER = -1
      RETURN
C
C X(I) .GE. X(I+1).
C
    5 DSMAX = DSM
      IER = -(I+1)
      RETURN
      END
