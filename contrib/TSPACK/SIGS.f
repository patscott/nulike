      SUBROUTINE SIGS (N,X,Y,YP,TOL, SIGMA, DSMAX,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), YP(N), TOL, SIGMA(N),
     .                 DSMAX
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
C   Given a set of abscissae X with associated data values Y
C and derivatives YP, this subroutine determines the small-
C est (nonnegative) tension factors SIGMA such that the Her-
C mite interpolatory tension spline H(x) preserves local
C shape properties of the data.  In an interval (X1,X2) with
C data values Y1,Y2 and derivatives YP1,YP2, the properties
C of the data are
C
C       Monotonicity:  S, YP1, and YP2 are nonnegative or
C                        nonpositive,
C  and
C       Convexity:     YP1 .LE. S .LE. YP2  or  YP1 .GE. S
C                        .GE. YP2,
C
C where S = (Y2-Y1)/(X2-X1).  The corresponding properties
C of H are constant sign of the first and second deriva-
C tives, respectively.  Note that, unless YP1 = S = YP2, in-
C finite tension is required (and H is linear on the inter-
C val) if S = 0 in the case of monotonicity, or if YP1 = S
C or YP2 = S in the case of convexity.
C
C   SIGS may be used in conjunction with Subroutine YPC2
C (or YPC2P) in order to produce a C-2 interpolant which
C preserves the shape properties of the data.  This is
C achieved by calling YPC2 with SIGMA initialized to the
C zero vector, and then alternating calls to SIGS with
C calls to YPC2 until the change in SIGMA is small (refer to
C the parameter descriptions for SIGMA, DSMAX and IER), or
C the maximum relative change in YP is bounded by a toler-
C ance (a reasonable value is .01).  A similar procedure may
C be used to produce a C-2 shape-preserving smoothing curve
C (Subroutine SMCRV).
C
C   Refer to Subroutine SIGBI for a means of selecting mini-
C mum tension factors to satisfy more general constraints.
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
C             sufficient to satisfy the constraint:
C             abs(TOL) is an upper bound on the magnitude
C             of the smallest (nonnegative) or largest (non-
C             positive) value of the first or second deriva-
C             tive of H in the interval.  Thus, the con-
C             straint is satisfied, but possibly with more
C             tension than necessary.  TOL should be set to
C             0 for optimal tension.
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
C On output:
C
C       SIGMA = Array containing tension factors for which
C               H(x) preserves the properties of the data,
C               with the restriction that SIGMA(I) .LE. 85
C               for all I (unless the input value is larger).
C               The factors are as small as possible (within
C               the tolerance), but not less than their
C               input values.  If infinite tension is re-
C               quired in interval (X(I),X(I+1)), then
C               SIGMA(I) = 85 (and H is an approximation to
C               the linear interpolant on the interval),
C               and if neither property is satisfied by the
C               data, then SIGMA(I) = 0 (unless the input
C               value is positive), and thus H is cubic in
C               the interval.
C
C       DSMAX = Maximum increase in a component of SIGMA
C               from its input value.  The increase is a
C               relative change if the input value is
C               nonzero, and an absolute change otherwise.
C
C       IER = Error indicator and information flag:
C             IER = I if no errors were encountered and I
C                     components of SIGMA were altered from
C                     their input values for 0 .LE. I .LE.
C                     N-1.
C             IER = -1 if N < 2.  SIGMA is not altered in
C                      this case.
C             IER = -I if X(I) .LE. X(I-1) for some I in the
C                      range 2,...,N.  SIGMA(J-1) is unal-
C                      tered for J = I,...,N in this case.
C
C Modules required by SIGS:  SNHCSH, STORE
C
C Intrinsic functions called by SIGS:  ABS, EXP, MAX, MIN,
C                                        SIGN, SQRT
C
C***********************************************************
C
      INTEGER I, ICNT, IP1, LUN, NIT, NM1, NITSTEP, NITMAX
      DOUBLE PRECISION A, C1, C2, COSHM, COSHMM, D0, D1,
     .                 D1D2, D1PD2, D2, DMAX, DSIG, DSM, DX,
     .                 E, EMS, EMS2, F, F0, FMAX, FNEG, FP,
     .                 FTOL, RTOL, S, S1, S2, SBIG, SCM,
     .                 SGN, SIG, SIGIN, SINHM, SSINH, SSM,
     .                 STOL, T, T0, T1, T2, TM, TP1
      DOUBLE PRECISION STORE
C
C The following line altered by Pat Scott Jan 31 2008
      DATA SBIG/85.D0/,  LUN/-1/, NITSTEP/200/, NITMAX/500/
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 9
C
C Compute an absolute tolerance FTOL = abs(TOL) and a
C   relative tolerance RTOL = 100*MACHEPS (now 200*MACHEPS).
C
      FTOL = ABS(TOL)
      RTOL = 1.D0
    1 RTOL = RTOL/2.D0
        IF (STORE(RTOL+1.D0) .GT. 1.D0) GO TO 1
C The following line altered by Pat Scott Jan 31 2008
      RTOL = RTOL*400.D0
C
C Initialize change counter ICNT and maximum change DSM for
C   loop on intervals.
C
      ICNT = 0
      DSM = 0.D0
      DO 8 I = 1,NM1
        IF (LUN .GE. 0) WRITE (LUN,100) I
  100   FORMAT (//1X,'SIGS -- INTERVAL',I4)
        IP1 = I + 1
        DX = X(IP1) - X(I)
        IF (DX .LE. 0.D0) GO TO 10
        SIGIN = SIGMA(I)
        IF (SIGIN .GE. SBIG) GO TO 8
C
C Compute first and second differences.
C
        S1 = YP(I)
        S2 = YP(IP1)
        S = (Y(IP1)-Y(I))/DX
        D1 = S - S1
        D2 = S2 - S
        D1D2 = D1*D2
C
C Test for infinite tension required to satisfy either
C   property.
C
        SIG = SBIG
        IF ((D1D2 .EQ. 0.D0  .AND.  S1 .NE. S2)  .OR.
     .      (S .EQ. 0.D0  .AND.  S1*S2 .GT. 0.D0)) GO TO 7
C
C Test for SIGMA = 0 sufficient.  The data satisfies convex-
C   ity iff D1D2 .GE. 0, and D1D2 = 0 implies S1 = S = S2.
C
        SIG = 0.D0
        IF (D1D2 .LT. 0.D0) GO TO 3
        IF (D1D2 .EQ. 0.D0) GO TO 7
        T = MAX(D1/D2,D2/D1)
        IF (T .LE. 2.D0) GO TO 7
        TP1 = T + 1.D0
C
C Convexity:  Find a zero of F(SIG) = SIG*COSHM(SIG)/
C   SINHM(SIG) - TP1.
C
C   F(0) = 2-T < 0, F(TP1) .GE. 0, the derivative of F
C     vanishes at SIG = 0, and the second derivative of F is
C     .2 at SIG = 0.  A quadratic approximation is used to
C     obtain a starting point for the Newton method.
C
        SIG = SQRT(10.D0*T-20.D0)
        NIT = 0
C
C   Top of loop:
C
    2   IF (SIG .LE. .5D0) THEN
          CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
          T1 = COSHM/SINHM
          FP = T1 + SIG*(SIG/SINHM - T1*T1 + 1.D0)
        ELSE
C
C   Scale SINHM and COSHM by 2*EXP(-SIG) in order to avoid
C     overflow with large SIG.
C
          EMS = EXP(-SIG)
          SSM = 1.D0 - EMS*(EMS+SIG+SIG)
          T1 = (1.D0-EMS)*(1.D0-EMS)/SSM
          FP = T1 + SIG*(2.D0*SIG*EMS/SSM - T1*T1 + 1.D0)
        ENDIF
C
        F = SIG*T1 - TP1
        IF (LUN .GE. 0) WRITE (LUN,110) SIG, F, FP
  110   FORMAT (5X,'CONVEXITY -- SIG = ',D15.8,
     .          ', F(SIG) = ',D15.8/1X,35X,'FP(SIG) = ',
     .          D15.8)
        NIT = NIT + 1
C
C   Test for convergence.
C
        IF (FP .LE. 0.D0) GO TO 7
        DSIG = -F/FP
        IF (ABS(DSIG) .LE. RTOL*SIG  .OR.  (F .GE. 0.D0
     .      .AND.  F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL)
     .    GO TO 7
C
C   Update SIG.
C
        SIG = SIG + DSIG
C The following lines added by Pat Scott Jan 31 2008
        IF (NIT .EQ. NITSTEP) THEN
          RTOL = RTOL * 2.D0
          WRITE(*,*) 'WARNING: RTOL in SIGS doubled to aid convergence'
        ENDIF
        IF (NIT .EQ. NITMAX) THEN
          WRITE(*,*) 'Error: SIGS did not converge'
          GO TO 9
        ENDIF
C end addition on Jan 31 2008
        GO TO 2
C
C Convexity cannot be satisfied.  Monotonicity can be satis-
C   fied iff S1*S .GE. 0 and S2*S .GE. 0 since S .NE. 0.
C
    3   IF (S1*S .LT. 0.D0  .OR.  S2*S .LT. 0.D0) GO TO 7
        T0 = 3.D0*S - S1 - S2
        D0 = T0*T0 - S1*S2
C
C SIGMA = 0 is sufficient for monotonicity iff S*T0 .GE. 0
C   or D0 .LE. 0.
C
        IF (D0 .LE. 0.D0  .OR.  S*T0 .GE. 0.D0) GO TO 7
C
C Monotonicity:  find a zero of F(SIG) = SIGN(S)*HP(R),
C   where HPP(R) = 0 and HP, HPP denote derivatives of H.
C   F has a unique zero, F(0) < 0, and F approaches abs(S)
C   as SIG increases.
C
C   Initialize parameters for the secant method.  The method
C     uses three points:  (SG0,F0), (SIG,F), and
C     (SNEG,FNEG), where SG0 and SNEG are defined implicitly
C     by DSIG = SIG - SG0 and DMAX = SIG - SNEG.
C
        SGN = SIGN(1.D0,S)
        SIG = SBIG
        FMAX = SGN*(SIG*S-S1-S2)/(SIG-2.D0)
        IF (FMAX .LE. 0.D0) GO TO 7
        STOL = RTOL*SIG
        F = FMAX
        F0 = SGN*D0/(3.D0*(D1-D2))
        FNEG = F0
        DSIG = SIG
        DMAX = SIG
        D1PD2 = D1 + D2
        NIT = 0
C
C   Top of loop:  compute the change in SIG by linear
C     interpolation.
C
    4   DSIG = -F*DSIG/(F-F0)
C   Modified by Pat Scott Oct 27 2017, to prevent 4-6 infinite loop due to roundoff error on Xenon Phi
        IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.
     .       DSIG*DMAX .GT. 0. ) THEN
          F0 = FNEG
          DSIG = -F*DMAX/(F-F0)
        ENDIF
        IF (LUN .GE. 0) WRITE (LUN,120) DSIG
  120   FORMAT (5X,'MONOTONICITY -- DSIG = ',D15.8)
C
C   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
C     Note that DSIG and DMAX have opposite signs.
C
        IF (ABS(DSIG) .LT. STOL/2.D0)
     .    DSIG = -SIGN(STOL/2.D0,DMAX)
C
C   Update SIG, F0, and F.
C
        SIG = SIG + DSIG
        F0 = F

C   Added by Pat Scott Oct 2017, to prevent F=NaN and SIG<0 due to roundoff error on Xenon Phi
        IF (SIG .LE. 0.D0) GO TO 7

        IF (SIG .LE. .5D0) THEN
C
C   Use approximations to the hyperbolic functions designed
C     to avoid cancellation error with small SIG.
C
          CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
          C1 = SIG*COSHM*D2 - SINHM*D1PD2
          C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
          A = C2 - C1
          E = SIG*SINHM - COSHMM - COSHMM
        ELSE
C
C   Scale SINHM and COSHM by 2*EXP(-SIG) in order to avoid
C     overflow with large SIG.
C
          EMS = EXP(-SIG)
          EMS2 = EMS + EMS
          TM = 1.D0 - EMS
          SSINH = TM*(1.D0+EMS)
          SSM = SSINH - SIG*EMS2
          SCM = TM*TM
          C1 = SIG*SCM*D2 - SSM*D1PD2
          C2 = SIG*SSINH*D2 - SCM*D1PD2
C
C   R is in (0,1) and well-defined iff HPP(X1)*HPP(X2) < 0.
C
          F = FMAX
          IF (C1*(SIG*SCM*D1 - SSM*D1PD2) .GE. 0.D0) GO TO 5
          A = EMS2*(SIG*TM*D2 + (TM-SIG)*D1PD2)
          IF (A*(C2+C1) .LT. 0.D0) GO TO 5
          E = SIG*SSINH - SCM - SCM
        ENDIF
C
        F = (SGN*(E*S2-C2) + SQRT(A*(C2+C1)))/E
C
C   Update number of iterations NIT.
C
    5   NIT = NIT + 1
        IF (LUN .GE. 0) WRITE (LUN,130) NIT, SIG, F
  130   FORMAT (1X,10X,I2,' -- SIG = ',D15.8,', F = ',
     .          D15.8)
C
C   Test for convergence.
C
        STOL = RTOL*SIG
        IF ( ABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.D0  .AND.
     .      F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL ) GO TO 7
        DMAX = DMAX + DSIG
        IF ( F0*F .GT. 0.D0  .AND.  ABS(F) .GE. ABS(F0) )
     .     GO TO 6
        IF (F0*F .LE. 0.D0) THEN
C
C   F and F0 have opposite signs.  Update (SNEG,FNEG) to
C     (SG0,F0) so that F and FNEG always have opposite
C     signs.  If SIG is closer to SNEG than SG0 and abs(F) <
C     abs(FNEG), then swap (SNEG,FNEG) with (SG0,F0).
C
          T1 = DMAX
          T2 = FNEG
          DMAX = DSIG
          FNEG = F0
          IF ( ABS(DSIG) .GT. ABS(T1)  .AND.
     .         ABS(F) .LT. ABS(T2)          ) THEN
C
            DSIG = T1
            F0 = T2
          ENDIF
        ENDIF
        GO TO 4
C
C   Bottom of loop:  F0*F > 0 and the new estimate would
C     be outside of the bracketing interval of length
C     abs(DMAX).  Reset (SG0,F0) to (SNEG,FNEG).
C
    6   DSIG = DMAX
        F0 = FNEG
        GO TO 4
C
C  Update SIGMA(I), ICNT, and DSM if necessary.
C
    7   SIG = MIN(SIG,SBIG)
        IF (SIG .GT. SIGIN) THEN
          SIGMA(I) = SIG
          ICNT = ICNT + 1
          DSIG = SIG-SIGIN
          IF (SIGIN .GT. 0.D0) DSIG = DSIG/SIGIN
          DSM = MAX(DSM,DSIG)
        ENDIF
    8   CONTINUE
C
C No errors encountered.
C
      DSMAX = DSM
      IER = ICNT
      RETURN
C
C N < 2.
C
    9 DSMAX = 0.D0
      IER = -1
      RETURN
C
C X(I+1) .LE. X(I).
C
   10 DSMAX = DSM
      IER = -IP1
      RETURN
      END
