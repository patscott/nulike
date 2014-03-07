      DOUBLE PRECISION FUNCTION SIG2 (X1,X2,Y1,Y2,Y1P,Y2P,
     .                                IFL,TOL, IER)
      INTEGER IFL, IER
      DOUBLE PRECISION X1, X2, Y1, Y2, Y1P, Y2P, TOL
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/08/92
C
C   Given a pair of abscissae with associated ordinates and
C slopes, this function determines the smallest (nonnega-
C tive) tension factor SIGMA such that the Hermite interpo-
C latory tension spline H(x) preserves convexity (or con-
C cavity) of the data;  i.e.,
C
C   Y1P .LE. S .LE. Y2P implies HPP(x) .GE. 0  or
C   Y1P .GE. S .GE. Y2P implies HPP(x) .LE. 0
C
C for all x in the open interval (X1,X2), where S = (Y2-Y1)/
C (X2-X1) and HPP denotes the second derivative of H.  Note,
C however, that infinite tension is required if Y1P = S or
C Y2P = S (unless Y1P = Y2P = S).
C
C On input:
C
C       X1,X2 = Abscissae.  X1 < X2.
C
C       Y1,Y2 = Values of H at X1 and X2.
C
C       Y1P,Y2P = Derivative values of H at X1 and X2.
C
C       IFL = Option indicator (sign of HPP):
C             IFL = -1 if HPP is to be bounded above by 0.
C             IFL = 1 if HPP is to be bounded below by 0
C                     (preserve convexity of the data).
C
C       TOL = Tolerance whose magnitude determines how close
C             SIGMA is to its optimal value when nonzero
C             finite tension is necessary and sufficient to
C             satisfy convexity or concavity.  In the case
C             of convexity, SIGMA is chosen so that 0 .LE.
C             HPPMIN .LE. abs(TOL), where HPPMIN is the min-
C             imum value of HPP in the interval.  In the
C             case of concavity, the maximum value of HPP
C             satisfies -abs(TOL) .LE. HPPMAX .LE. 0.  Thus,
C             the constraint is satisfied but possibly with
C             more tension than necessary.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and fin-
C                     ite tension is sufficient to satisfy
C                     the constraint.
C             IER = 1 if no errors were encountered but in-
C                     finite tension is required to satisfy
C                     the constraint.
C             IER = -1 if X2 .LE. X1 or abs(IFL) .NE. 1.
C             IER = -2 if the constraint cannot be satis-
C                      fied:  the sign of S-Y1P or Y2P-S
C                      does not agree with IFL.
C
C       SIG2 = Tension factor defined above unless IER < 0,
C              in which case SIG2 = -1.  If IER = 1, SIG2
C              is set to 85, resulting in an approximation
C              to the linear interpolant of the endpoint
C              values.  Note, however, that SIG2 may be
C              larger than 85 if IER = 0.
C
C Modules required by SIG2:  SNHCSH, STORE
C
C Intrinsic functions called by SIG2:  ABS, EXP, MAX, MIN,
C                                        SQRT
C
C***********************************************************
C
      INTEGER LUN, NIT
      DOUBLE PRECISION COSHM, D1, D2, DSIG, DUMMY, DX, EMS,
     .                 F, FP, FTOL, RTOL, S, SBIG, SIG,
     .                 SINHM, SSM, T, T1, TP1
      DOUBLE PRECISION STORE
C
      DATA SBIG/85.D0/, LUN/-1/
C
C Test for an errors in the input parameters.
C
      DX = X2 - X1
      IF (ABS(IFL) .NE. 1.D0  .OR.  DX .LE. 0.D0) GO TO 5
C
C Compute the slope and second differences, and test for
C   an invalid constraint.
C
      S = (Y2-Y1)/DX
      D1 = S - Y1P
      D2 = Y2P - S
      IF ((IFL .GT. 0.D0  .AND.  MIN(D1,D2) .LT. 0.D0)
     .    .OR.  (IFL .LT. 0.D0  .AND.
     .           MAX(D1,D2) .GT. 0.D0)) GO TO 6
C
C Test for infinite tension required.
C
      IF (D1*D2 .EQ. 0.D0  .AND.  D1 .NE. D2) GO TO 4
C
C Test for SIG = 0 sufficient.
C
      SIG = 0.D0
      IF (D1*D2 .EQ. 0.D0) GO TO 3
      T = MAX(D1/D2,D2/D1)
      IF (T .LE. 2.D0) GO TO 3
C
C Find a zero of F(SIG) = SIG*COSHM(SIG)/SINHM(SIG) - (T+1).
C   Since the derivative of F vanishes at the origin, a
C   quadratic approximation is used to obtain an initial
C   estimate for the Newton method.
C
      TP1 = T + 1.D0
      SIG = SQRT(10.D0*T-20.D0)
      NIT = 0
C
C   Compute an absolute tolerance FTOL = abs(TOL) and a
C     relative tolerance RTOL = 100*MACHEPS.
C
      FTOL = ABS(TOL)
      RTOL = 1.D0
    1 RTOL = RTOL/2.D0
        IF (STORE(RTOL+1.D0) .GT. 1.D0) GO TO 1
      RTOL = RTOL*200.D0
C
C Evaluate F and its derivative FP.
C
    2 IF (SIG .LE. .5D0) THEN
C
C   Use approximations designed to avoid cancellation error
C     in the hyperbolic functions.
C
        CALL SNHCSH (SIG, SINHM,COSHM,DUMMY)
        T1 = COSHM/SINHM
        FP = T1 + SIG*(SIG/SINHM - T1*T1 + 1.D0)
      ELSE
C
C   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid
C     overflow.
C
        EMS = EXP(-SIG)
        SSM = 1.D0 - EMS*(EMS+SIG+SIG)
        T1 = (1.D0-EMS)*(1.D0-EMS)/SSM
        FP = T1 + SIG*(2.D0*SIG*EMS/SSM - T1*T1 + 1.D0)
      ENDIF
C
      F = SIG*T1 - TP1
      IF (LUN .GE. 0) WRITE (LUN,100) SIG, F, FP
  100 FORMAT (1X,'SIG2 -- SIG = ',D15.8,', F(SIG) = ',
     .        D15.8/1X,29X,'FP(SIG) = ',D15.8)
      NIT = NIT + 1
C
C   Test for convergence.
C
      IF (FP .LE. 0.D0) GO TO 3
      DSIG = -F/FP
      IF (ABS(DSIG) .LE. RTOL*SIG  .OR.  (F .GE. 0.D0 .AND.
     .    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 3
C
C   Update SIG.
C
      SIG = SIG + DSIG
      GO TO 2
C
C No errors encountered, and SIGMA is finite.
C
    3 IER = 0
      SIG2 = SIG
      RETURN
C
C Infinite tension required.
C
    4 IER = 1
      SIG2 = SBIG
      RETURN
C
C X2 .LE. X1 or abs(IFL) .NE. 1.
C
    5 IER = -1
      SIG2 = -1.D0
      RETURN
C
C The constraint cannot be satisfied.
C
    6 IER = -2
      SIG2 = -1.D0
      RETURN
      END
