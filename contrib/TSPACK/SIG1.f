      DOUBLE PRECISION FUNCTION SIG1 (X1,X2,Y1,Y2,Y1P,Y2P,
     .                                IFL,HPBND,TOL, IER)
      INTEGER IFL, IER
      DOUBLE PRECISION X1, X2, Y1, Y2, Y1P, Y2P, HPBND, TOL
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/17/96
C
C   Given a pair of abscissae with associated ordinates and
C slopes, this function determines the smallest (nonnega-
C tive) tension factor SIGMA such that the derivative HP(x)
C of the Hermite interpolatory tension spline H(x), defined
C by SIGMA and the data, is bounded (either above or below)
C by HPBND for all x in (X1,X2).
C
C On input:
C
C       X1,X2 = Abscissae.  X1 < X2.
C
C       Y1,Y2 = Values of H at X1 and X2.
C
C       Y1P,Y2P = Values of HP at X1 and X2.
C
C       IFL = Option indicator:
C             IFL = -1 if HPBND is a lower bound on HP.
C             IFL = 1 if HPBND is an upper bound on HP.
C
C       HPBND = Bound on HP.  If IFL = -1, HPBND .LE.
C               min(Y1P,Y2P,S) for S = (Y2-Y1)/(X2-X1).  If
C               IFL = 1, HPBND .GE. max(Y1P,Y2P,S).
C
C       TOL = Tolerance whose magnitude determines how close
C             SIGMA is to its optimal value when nonzero
C             finite tension is necessary and sufficient to
C             satisfy the constraint.  For a lower bound,
C             SIGMA is chosen so that HPBND .LE. HPMIN .LE.
C             HPBND + abs(TOL), where HPMIN is the minimum
C             value of HP in the interval, and for an upper
C             bound, the maximum of HP satisfies HPBND -
C             abs(TOL) .LE. HPMAX .LE. HPBND.  Thus, the
C             constraint is satisfied but possibly with more
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
C                     the constraint (e.g., IFL = -1, HPBND
C                     = S, and Y1P > S).
C             IER = -1 if X2 .LE. X1 or abs(IFL) .NE. 1.
C             IER = -2 if HPBND is outside its valid range
C                      on input.
C
C       SIG1 = Minimum tension factor defined above unless
C              IER < 0, in which case SIG1 = -1.  If IER =
C              1, SIG1 = 85, resulting in an approximation
C              to the linear interpolant of the endpoint
C              values.  Note, however, that SIG1 may be
C              larger than 85 if IER = 0.
C
C Modules required by SIG1:  SNHCSH, STORE
C
C Intrinsic functions called by SIG1:  ABS, DBLE, EXP, MAX,
C                                        MIN, SIGN, SQRT
C
C***********************************************************
C
      INTEGER LUN, NIT
      DOUBLE PRECISION A, A0, B0, BND, C0, C1, C2, COSHM,
     .                 COSHMM, D0, D1, D1PD2, D2, DMAX,
     .                 DSIG, DX, E, EMS, EMS2, F, F0, FMAX,
     .                 FNEG, FTOL, RF, RTOL, S, S1, S2,
     .                 SBIG, SIG, SINH, SINHM, STOL, T0, T1,
     .                 T2, TM
      DOUBLE PRECISION STORE
C
      DATA SBIG/85.D0/,  LUN/-1/
C
C Store local parameters and test for errors.
C
      RF = DBLE(IFL)
      DX = X2 - X1
      IF (ABS(RF) .NE. 1.D0  .OR.  DX .LE. 0.D0) GO TO 7
      S1 = Y1P
      S2 = Y2P
      S = (Y2-Y1)/DX
      BND = HPBND
C
C Test for a valid constraint.
C
      IF ((RF .LT. 0.D0  .AND.  MIN(S1,S2,S) .LT. BND)
     .    .OR.  (RF .GT. 0.D0  .AND.
     .           BND .LT. MAX(S1,S2,S))) GO TO 8
C
C Test for infinite tension required.
C
      IF (S .EQ. BND  .AND.  (S1 .NE. S  .OR.  S2 .NE. S))
     .   GO TO 6
C
C Test for SIG = 0 sufficient.  The Hermite cubic interpo-
C   land H0 has derivative HP0(x) = S2 + 2*B0*R + A0*R**2,
C   where R = (X2-x)/DX.
C
      SIG = 0.D0
      T0 = 3.D0*S - S1 - S2
      B0 = T0 - S2
      C0 = T0 - S1
      A0 = -B0 - C0
C
C   HP0(R) has an extremum (at R = -B0/A0) in (0,1) iff
C     B0*C0 > 0 and the third derivative of H0 has the
C     sign of A0.
C
      IF (B0*C0 .LE. 0.D0  .OR.  A0*RF .GT. 0.D0) GO TO 5
C
C   A0*RF < 0 and HP0(R) = -D0/A0 at R = -B0/A0.
C
      D0 = T0*T0 - S1*S2
      F0 = (BND + D0/A0)*RF
      IF (F0 .GE. 0.D0) GO TO 5
C
C Find a zero of F(SIG) = (BND-HP(R))*RF, where HP has an
C   extremum at R.  F has a unique zero, F(0) = F0 < 0, and
C   F = (BND-S)*RF > 0 for SIG sufficiently large.
C
C Initialize parameters for the secant method.  The method
C   uses three points:  (SG0,F0), (SIG,F), and (SNEG,FNEG),
C   where SG0 and SNEG are defined implicitly by DSIG = SIG
C   - SG0 and DMAX = SIG - SNEG.  SG0 is initially zero and
C   SIG is initialized to the zero of (BND - (SIG*S-S1-S2)/
C   (SIG-2.))*RF -- a value for which F(SIG) .GE. 0 and
C   F(SIG) = 0 for SIG sufficiently large that 2*SIG is in-
C   significant relative to EXP(SIG).
C
      FMAX = (BND-S)*RF
      IF (LUN .GE. 0  .AND.  RF .LT. 0.D0)
     .  WRITE (LUN,100) F0, FMAX
      IF (LUN .GE. 0  .AND.  RF .GT. 0.D0)
     .  WRITE (LUN,110) F0, FMAX
  100 FORMAT (//1X,'SIG1 (LOWER BOUND) -- F(0) = ',D15.8,
     .        ', FMAX = ',D15.8/)
  110 FORMAT (//1X,'SIG1 (UPPER BOUND) -- F(0) = ',D15.8,
     .        ', FMAX = ',D15.8/)
      SIG = 2.D0 - A0/(3.D0*(BND-S))
      IF (STORE(SIG*EXP(-SIG)+.5D0) .EQ. .5D0) GO TO 5
      DSIG = SIG
      DMAX = -2.D0*SIG
      FNEG = FMAX
      D1 = S - S1
      D2 = S2 - S
      D1PD2 = D1 + D2
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
    2 IF (SIG .LE. .5D0) THEN
C
C   Use approximations designed to avoid cancellation error
C     (associated with small SIG) in the modified hyperbolic
C     functions.
C
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        C1 = SIG*COSHM*D2 - SINHM*D1PD2
        C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
        A = C2 - C1
        E = SIG*SINHM - COSHMM - COSHMM
      ELSE
C
C   Scale SINHM and COSHM by 2*EXP(-SIG) in order to avoid
C     overflow.
C
        EMS = EXP(-SIG)
        EMS2 = EMS + EMS
        TM = 1.D0 - EMS
        SINH = TM*(1.D0+EMS)
        SINHM = SINH - SIG*EMS2
        COSHM = TM*TM
        C1 = SIG*COSHM*D2 - SINHM*D1PD2
        C2 = SIG*SINH*D2 - COSHM*D1PD2
        A = EMS2*(SIG*TM*D2 + (TM-SIG)*D1PD2)
        E = SIG*SINH - COSHM - COSHM
      ENDIF
C
C   The second derivative of H(R) has a zero at EXP(SIG*R) =
C     SQRT((C2+C1)/A) and R is in (0,1) and well-defined
C     iff HPP(X1)*HPP(X2) < 0.
C
      F = FMAX
      T1 = A*(C2+C1)
      IF (T1 .GE. 0.D0) THEN
        IF (C1*(SIG*COSHM*D1 - SINHM*D1PD2) .LT. 0.D0) THEN
C
C   HP(R) = (B+SIGN(A)*SQRT(A*C))/E at the critical value
C     of R, where A = C2-C1, B = E*S2-C2, and C = C2+C1.
C     NOTE THAT RF*A < 0.
C
          F = (BND - (E*S2-C2 - RF*SQRT(T1))/E)*RF
        ENDIF
      ENDIF
C
C   Update the number of iterations NIT.
C
      NIT = NIT + 1
      IF (LUN .GE. 0) WRITE (LUN,120) NIT, SIG, F
  120 FORMAT (1X,3X,I2,' -- SIG = ',D15.8,', F = ',
     .        D15.8)
      IF (F0*F .LT. 0.D0) THEN
C
C   F0*F < 0.  Update (SNEG,FNEG) to (SG0,F0) so that F
C     and FNEG always have opposite signs.  If SIG is closer
C     to SNEG than SG0 and abs(F) < abs(FNEG), then swap
C     (SNEG,FNEG) with (SG0,F0).
C
        T1 = DMAX
        T2 = FNEG
        DMAX = DSIG
        FNEG = F0
        IF ( ABS(DSIG) .GT. ABS(T1)  .AND.
     .       ABS(F) .LT. ABS(T2)          ) THEN
C
          DSIG = T1
          F0 = T2
        ENDIF
      ENDIF
C
C   Test for convergence.
C
      STOL = RTOL*SIG
      IF (ABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.D0  .AND.
     .    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 5
      IF (F0*F .LT. 0  .OR.  ABS(F) .LT. ABS(F0)) GO TO 4
C
C   F*F0 > 0 and the new estimate would be outside of the
C     bracketing interval of length abs(DMAX).  Reset
C     (SG0,F0) to (SNEG,FNEG).
C
    3 DSIG = DMAX
      F0 = FNEG
C
C   Compute the change in SIG by linear interpolation be-
C     tween (SG0,F0) and (SIG,F).
C
    4 DSIG = -F*DSIG/(F-F0)
      IF (LUN .GE. 0) WRITE (LUN,130) DSIG
  130 FORMAT (1X,8X,'DSIG = ',D15.8)
      IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.
     .     DSIG*DMAX .GT. 0. ) GO TO 3
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
    5 IER = 0
      SIG1 = SIG
      RETURN
C
C Infinite tension required.
C
    6 IER = 1
      SIG1 = SBIG
      RETURN
C
C Error in an input parameter.
C
    7 IER = -1
      SIG1 = -1.D0
      RETURN
C
C Invalid constraint.
C
    8 IER = -2
      SIG1 = -1.D0
      RETURN
      END
