      DOUBLE PRECISION FUNCTION SIG0 (X1,X2,Y1,Y2,Y1P,Y2P,
     .                                IFL,HBND,TOL, IER)
      INTEGER IFL, IER
      DOUBLE PRECISION X1, X2, Y1, Y2, Y1P, Y2P, HBND, TOL
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/18/96
C
C   Given a pair of abscissae with associated ordinates and
C slopes, this function determines the smallest (nonnega-
C tive) tension factor SIGMA such that the Hermite interpo-
C latory tension spline H(x), defined by SIGMA and the data,
C is bounded (either above or below) by HBND for all x in
C (X1,X2).
C
C On input:
C
C       X1,X2 = Abscissae.  X1 < X2.
C
C       Y1,Y2 = Values of H at X1 and X2.
C
C       Y1P,Y2P = Derivative values of H at X1 and X2.
C
C       IFL = Option indicator:
C             IFL = -1 if HBND is a lower bound on H.
C             IFL = 1 if HBND is an upper bound on H.
C
C       HBND = Bound on H.  If IFL = -1, HBND .LE. min(Y1,
C              Y2).  If IFL = 1, HBND .GE. max(Y1,Y2).
C
C       TOL = Tolerance whose magnitude determines how close
C             SIGMA is to its optimal value when nonzero
C             finite tension is necessary and sufficient to
C             satisfy the constraint.  For a lower bound,
C             SIGMA is chosen so that HBND .LE. HMIN .LE.
C             HBND + abs(TOL), where HMIN is the minimum
C             value of H in the interval, and for an upper
C             bound, the maximum of H satisfies HBND -
C             abs(TOL) .LE. HMAX .LE. HBND.  Thus, the con-
C             straint is satisfied but possibly with more
C             tension than necessary.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and the
C                     constraint can be satisfied with fin-
C                     ite tension.
C             IER = 1 if no errors were encountered but in-
C                     finite tension is required to satisfy
C                     the constraint (e.g., IFL = -1, HBND
C                     = Y1, and Y1P < 0.).
C             IER = -1 if X2 .LE. X1 or abs(IFL) .NE. 1.
C             IER = -2 if HBND is outside its valid range
C                      on input.
C
C       SIG0 = Minimum tension factor defined above unless
C              IER < 0, in which case SIG0 = -1.  If IER =
C              1, SIG0 = 85, resulting in an approximation
C              to the linear interpolant of the endpoint
C              values.  Note, however, that SIG0 may be
C              larger than 85 if IER = 0.
C
C Modules required by SIG0:  SNHCSH, STORE
C
C Intrinsic functions called by SIG0:  ABS, DBLE, EXP, LOG,
C                                        MAX, MIN, SIGN,
C                                        SQRT
C
C***********************************************************
C
      INTEGER LUN, NIT
      DOUBLE PRECISION A, A0, AA, B, B0, BND, C, C1, C2,
     .                 COSHM, COSHMM, D, D0, D1PD2, D2,
     .                 DMAX, DSIG, DX, E, EMS, F, F0, FMAX,
     .                 FNEG, FTOL, R, RF, RSIG, RTOL, S, S1,
     .                 S2, SBIG, SCM, SIG, SINHM, SNEG,
     .                 SSINH, SSM, STOL, T, T0, T1, T2, TM,
     .                 Y1L, Y2L
      DOUBLE PRECISION STORE
C
      DATA SBIG/85.D0/,  LUN/-1/
C
C Store local parameters and test for errors.
C
      RF = DBLE(IFL)
      DX = X2 - X1
      IF (ABS(RF) .NE. 1.D0  .OR.  DX .LE. 0.D0) GO TO 8
      Y1L = Y1
      Y2L = Y2
      BND = HBND
C
C Test for a valid constraint.
C
      IF ((RF .LT. 0.D0  .AND.  MIN(Y1L,Y2L) .LT. BND)
     .    .OR.  (RF .GT. 0.D0  .AND.
     .           BND .LT. MAX(Y1L,Y2L))) GO TO 9
C
C Test for infinite tension required.
C
      S1 = Y1P
      S2 = Y2P
      IF ((Y1L .EQ. BND  .AND.  RF*S1 .GT. 0.D0)  .OR.
     .    (Y2L .EQ. BND  .AND.  RF*S2 .LT. 0.D0)) GO TO 7
C
C Test for SIG = 0 sufficient.
C
      SIG = 0.D0
      IF (RF*S1 .LE. 0.D0  .AND.  RF*S2 .GE. 0.D0) GO TO 6
C
C   Compute coefficients A0 and B0 of the Hermite cubic in-
C     terpolant H0(x) = Y2 - DX*(S2*R + B0*R**2 + A0*R**3/3)
C     where R = (X2-x)/DX.
C
      S = (Y2L-Y1L)/DX
      T0 = 3.D0*S - S1 - S2
      A0 = 3.D0*(S-T0)
      B0 = T0 - S2
      D0 = T0*T0 - S1*S2
C
C   H0 has local extrema in (X1,X2) iff S1*S2 < 0 or
C     (T0*(S1+S2) < 0 and D0 .GE. 0).
C
      IF (S1*S2 .GE. 0.D0  .AND.  (T0*(S1+S2) .GE. 0.D0
     .    .OR.  D0 .LT. 0.D0)) GO TO 6
      IF (A0 .EQ. 0.D0) THEN
C
C   H0 is quadratic and has an extremum at R = -S2/(2*B0).
C     H0(R) = Y2 + DX*S2**2/(4*B0).  Note that A0 = 0 im-
C     plies 2*B0 = S1-S2, and S1*S2 < 0 implies B0 .NE. 0.
C     Also, the extremum is a min iff HBND is a lower bound.
C
        F0 = (BND - Y2L - DX*S2*S2/(4.D0*B0))*RF
      ELSE
C
C   A0 .NE. 0 and H0 has extrema at R = (-B0 +/- SQRT(D0))/
C     A0 = S2/(-B0 -/+ SQRT(D0)), where the negative root
C     corresponds to a min.  The expression for R is chosen
C     to avoid cancellation error.  H0(R) = Y2 + DX*(S2*B0 +
C     2*D0*R)/(3*A0).
C
        T = -B0 - SIGN(SQRT(D0),B0)
        R = T/A0
        IF (RF*B0 .GT. 0.D0) R = S2/T
        F0 = (BND - Y2L - DX*(S2*B0+2.D0*D0*R)/(3.D0*A0))*RF
      ENDIF
C
C   F0 .GE. 0 iff SIG = 0 is sufficient to satisfy the
C     constraint.
C
      IF (F0 .GE. 0.D0) GO TO 6
C
C Find a zero of F(SIG) = (BND-H(R))*RF where the derivative
C   of H, HP, vanishes at R.  F is a nondecreasing function,
C   F(0) < 0, and F = FMAX for SIG sufficiently large.
C
C Initialize parameters for the secant method.  The method
C   uses three points:  (SG0,F0), (SIG,F), and (SNEG,FNEG),
C   where SG0 and SNEG are defined implicitly by DSIG = SIG
C   - SG0 and DMAX = SIG - SNEG.  SG0 is initially zero and
C   SNEG is initialized to a sufficiently large value that
C   FNEG > 0.  This value is used only if the initial value
C   of F is negative.
C
      FMAX = MAX(1.D-3,MIN(ABS(Y1L-BND),ABS(Y2L-BND)))
      T = MAX(ABS(Y1L-BND),ABS(Y2L-BND))
      SIG = DX*MAX(ABS(S1),ABS(S2))/T
      DMAX = SIG*(1.D0-T/FMAX)
      SNEG = SIG - DMAX
      IF (LUN .GE. 0  .AND.  RF .LT. 0.D0)
     .   WRITE (LUN,100) F0, FMAX, SNEG
      IF (LUN .GE. 0  .AND.  RF .GT. 0.D0)
     .   WRITE (LUN,110) F0, FMAX, SNEG
  100 FORMAT (//1X,'SIG0 (LOWER BOUND) -- F(0) = ',D15.8,
     .        ', FMAX = ',D15.8/1X,46X,'SNEG = ',D15.8/)
  110 FORMAT (//1X,'SIG0 (UPPER BOUND) -- F(0) = ',D15.8,
     .        ', FMAX = ',D15.8/1X,46X,'SNEG = ',D15.8/)
      DSIG = SIG
      FNEG = FMAX
      D2 = S2 - S
      D1PD2 = S2 - S1
      NIT = 0
C
C Compute an absolute tolerance FTOL = abs(TOL), and a
C   relative tolerance RTOL = 100*MACHEPS.
C
      FTOL = ABS(TOL)
      RTOL = 1.D0
    1 RTOL = RTOL/2.D0
        IF (STORE(RTOL+1.D0) .GT. 1.D0) GO TO 1
      RTOL = RTOL*200.D0
C
C Top of loop:  compute F.
C
    2 EMS = EXP(-SIG)
      IF (SIG .LE. .5D0) THEN
C
C   SIG .LE. .5:  use approximations designed to avoid can-
C                 cellation error (associated with small
C                 SIG) in the modified hyperbolic functions.
C
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        C1 = SIG*COSHM*D2 - SINHM*D1PD2
        C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
        A = C2 - C1
        AA = A/EMS
        E = SIG*SINHM - COSHMM - COSHMM
      ELSE
C
C   SIG > .5:  scale SINHM and COSHM by 2*EXP(-SIG) in order
C              to avoid overflow.
C
        TM = 1.D0 - EMS
        SSINH = TM*(1.D0+EMS)
        SSM = SSINH - 2.D0*SIG*EMS
        SCM = TM*TM
        C1 = SIG*SCM*D2 - SSM*D1PD2
        C2 = SIG*SSINH*D2 - SCM*D1PD2
        AA = 2.D0*(SIG*TM*D2 + (TM-SIG)*D1PD2)
        A = EMS*AA
        E = SIG*SSINH - SCM - SCM
      ENDIF
C
C   HP(R) = S2 - (C1*SINH(SIG*R) - C2*COSHM(SIG*R))/E = 0
C     for ESR = (-B +/- SQRT(D))/A = C/(-B -/+ SQRT(D)),
C     where ESR = EXP(SIG*R), A = C2-C1, D = B**2 - A*C, and
C     B and C are defined below.
C
      B = E*S2 - C2
      C = C2 + C1
      D = B*B - A*C
      F = 0.D0
      IF (AA*C .EQ. 0.D0  .AND.  B .EQ. 0.D0) GO TO 3
      F = FMAX
      IF (D .LT. 0.D0) GO TO 3
      T1 = SQRT(D)
      T = -B - SIGN(T1,B)
      RSIG = 0.D0
      IF (RF*B .LT. 0.D0  .AND.  AA .NE. 0.) THEN
        IF (T/AA .GT. 0.D0) RSIG = SIG + LOG(T/AA)
      ENDIF
      IF ((RF*B .GT. 0.D0  .OR.  AA .EQ. 0.D0)  .AND.
     .    C/T .GT. 0.D0) RSIG = LOG(C/T)
      IF ((RSIG .LE. 0.D0  .OR.  RSIG .GE. SIG)  .AND.
     .    B .NE. 0.D0) GO TO 3
C
C   H(R) = Y2 - DX*(B*SIG*R + C1 + RF*SQRT(D))/(SIG*E).
C
      F = (BND - Y2L + DX*(B*RSIG+C1+RF*T1)/(SIG*E))*RF
C
C   Update the number of iterations NIT.
C
    3 NIT = NIT + 1
      IF (LUN .GE. 0) WRITE (LUN,120) NIT, SIG, F
  120 FORMAT (1X,3X,I2,' -- SIG = ',D15.8,', F = ',
     .        D15.8)
      IF (F0*F .LT. 0.D0) THEN
C
C   F0*F < 0.  Update (SNEG,FNEG) to (SG0,F0) so that F
C     and FNEG always have opposite signs.  If SIG is
C     closer to SNEG than SG0, then swap (SNEG,FNEG) with
C     (SG0,F0).
C
        T1 = DMAX
        T2 = FNEG
        DMAX = DSIG
        FNEG = F0
        IF (ABS(DSIG) .GT. ABS(T1)) THEN
          DSIG = T1
          F0 = T2
        ENDIF
      ENDIF
C
C   Test for convergence.
C
      STOL = RTOL*SIG
      IF (ABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.D0  .AND.
     .    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 6
C
C   Test for F0 = F = FMAX or F < 0 on the first iteration.
C
      IF (F0 .NE. F  .AND.  (NIT .GT. 1  .OR.  F .GT. 0.D0))
     .   GO TO 5
C
C   F*F0 > 0 and either the new estimate would be outside of
C     the bracketing interval of length abs(DMAX) or F < 0
C     on the first iteration.  Reset (SG0,F0) to (SNEG,FNEG).
C
    4 DSIG = DMAX
      F0 = FNEG
C
C   Compute the change in SIG by linear interpolation be-
C     tween (SG0,F0) and (SIG,F).
C
    5 DSIG = -F*DSIG/(F-F0)
      IF (LUN .GE. 0) WRITE (LUN,130) DSIG
  130 FORMAT (1X,8X,'DSIG = ',D15.8)
      IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.
     .     DSIG*DMAX .GT. 0. ) GO TO 4
C
C   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
C     Note that DSIG and DMAX have opposite signs.
C
      IF (ABS(DSIG) .LT. STOL/2.D0)
     .  DSIG = -SIGN(STOL/2.D0,DMAX)
C
C   Bottom of loop:  update SIG, DMAX, and F0.
C
      SIG = SIG + DSIG
      DMAX = DMAX + DSIG
      F0 = F
      GO TO 2
C
C No errors encountered and SIGMA finite.
C
    6 IER = 0
      SIG0 = SIG
      RETURN
C
C Infinite tension required.
C
    7 IER = 1
      SIG0 = SBIG
      RETURN
C
C Error in an input parameter.
C
    8 IER = -1
      SIG0 = -1.D0
      RETURN
C
C Invalid constraint.
C
    9 IER = -2
      SIG0 = -1.D0
      RETURN
      END
