      SUBROUTINE SIGBP (N,X,Y,XP,YP,TOL,BL,BU,
     .                  BMAX, SIGMA, DSMAX,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), XP(N), YP(N), TOL, BL(N),
     .                 BU(N), BMAX, SIGMA(N), DSMAX
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
C   Given an ordered sequence of points C(I) = (X(I),Y(I))
C with associated derivative vectors CP(I) = (XP(I),YP(I)),
C this subroutine determines the smallest (nonnegative) ten-
C sion factors SIGMA such that a parametric planar curve
C C(t) satisfies a set of user-specified constraints.  The
C components x(t) and y(t) of C(t) are the Hermite interpo-
C latory tension splines defined by the data and tension
C factors:  C(t(I)) = C(I) and C'(t(I)) = CP(I) for para-
C meter values t(1), t(2), ..., t(N).  In each subinterval
C [t1,t2], the signed perpendicular distance from the
C corresponding line segment C1-C2 to the curve C(t) is
C given by the vector cross product
C
C     d(t) = (C2-C1)/DC X (C(t)-C1)
C
C where DC = abs(C2-C1) is the length of the line segment.
C The associated tension factor SIGMA is chosen to satisfy
C an upper bound on the maximum of d(t) and a lower bound on
C the minimum of d(t) over t in [t1,t2].  Thus, the upper
C bound is associated with distance to the left of the line
C segment as viewed from C1 toward C2.  Note that the curve
C is assumed to be parameterized by arc length (Subroutine
C ARCL2D) so that t2-t1 = DC.  If this is not the case, the
C required bounds should be scaled by DC/(t2-t1) to obtain
C the input parameters BL and BU.
C
C   SIGBP may be used in conjunction with Subroutine YPC2
C (or YPC2P) in order to produce a C-2 interpolant which
C satisfies the constraints.  This is achieved by calling
C YPC2 with SIGMA initialized to the zero vector, and then
C alternating calls to SIGBP with calls to YPC2 until the
C change in SIGMA is small (refer to the parameter descrip-
C tions for SIGMA, DSMAX and IER), or the maximum relative
C change in YP is bounded by a tolerance (a reasonable value
C is .01).  A similar procedure may be used to produce a C-2
C shape-preserving smoothing curve (Subroutine SMCRV).
C
C On input:
C
C       N = Number of data points.  N .GE. 2.
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of the points C(I), I = 1 to N.
C
C       XP,YP = Arrays of length N containing the components
C               of the derivative (velocity) vectors CP(I).
C               Refer to Subroutines YPC1, YPC1P, YPC2,
C               YPC2P, and SMCRV.
C
C       TOL = Tolerance whose magnitude determines how close
C             each tension factor SIGMA is to its optimal
C             value when nonzero finite tension is necessary
C             and sufficient to satisfy a constraint.
C             SIGMA(I) is chosen so that BL(I) .LE. dmin
C             .LE. BL(I) + abs(TOL) and BU(I) - abs(TOL)
C             .LE. dmax .LE. BU(I), where dmin and dmax are
C             the minimum and maximum values of d(t) in the
C             interval [t(I),t(I+1)].  Thus, a large toler-
C             ance might increase execution efficiency but
C             may result in more tension than is necessary.
C             TOL may be set to 0 for optimal tension.
C
C       BL,BU = Arrays of length N-1 containing lower and
C               upper bounds, respectively, which define
C               the constraints as described above.  BL(I)
C               < 0 and BU(I) > 0 for I = 1 to N-1.  A null
C               straint is specified by BL(I) .LE. -BMAX or
C               BU(I) .GE. BMAX.
C
C       BMAX = User-defined value of infinity which, when
C              used as an upper bound in BU (or when its
C              negative is used as a lower bound in BL),
C              specifies that no constraint is to be en-
C              forced.
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
C               d(t) satisfies the constraints defined by
C               BL and BU, with the restriction that
C               SIGMA(I) .LE. 85 for all I (unless the input
C               value is larger).  The factors are as small
C               as possible (within the tolerance), but not
C               less than their input values.  If no con-
C               straint is specified in interval I, then
C               SIGMA(I) = 0 (unless the input value is
C               positive), and thus x(t) and y(t) are cubic
C               polynomials.
C
C       DSMAX = Maximum increase in a component of SIGMA
C               from its input value.  The increase is a
C               relative change if the input value is
C               positive, and an absolute change otherwise.
C
C       IER = Error indicator and information flag:
C             IER = I if no errors were encountered and I
C                     components of SIGMA were altered from
C                     their input values for 0 .LE. I .LE.
C                     N-1.
C             IER = -1 if N < 2.  SIGMA is not altered in
C                      this case.
C             IER = -I if BL(I-1) .GE. 0 or BU(I-1) .LE. 0
C                      for some I in the range 2 to N.
C                      SIGMA(J) is unaltered for J .GE. I-1
C                      in this case.
C
C Modules required by SIGBP:  SNHCSH, STORE
C
C Intrinsic functions called by SIGBP:  ABS, EXP, LOG, MAX,
C                                         MIN, SIGN, SQRT
C
C***********************************************************
C
      INTEGER I, ICNT, IP1, LUN, NIT, NM1
      DOUBLE PRECISION A, A1, A2, AA, B, B0, BHI, BLO, BMX,
     .                 C, COSHM, COSHMM, D, D0, DM, DMAX,
     .                 DP, DSIG, DSM, DX, DY, E, EB, EMS, F,
     .                 F0, FMAX, FNEG, FTOL, RM, RP, RSM,
     .                 RSP, RTOL, S, SBIG, SIG, SIGIN, SINH,
     .                 SINHM, SNEG, STOL, T, T1, T2, TM, V1,
     .                 V2, V2M1
      DOUBLE PRECISION STORE
C
      DATA SBIG/85.D0/,  LUN/-1/
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 8
      BMX = BMAX
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
C Initialize change counter ICNT and maximum change DSM for
C   loop on intervals.
C
      ICNT = 0
      DSM = 0.D0
      DO 7 I = 1,NM1
        IP1 = I + 1
        BLO = BL(I)
        BHI = BU(I)
        SIGIN = SIGMA(I)
        IF (LUN .GE. 0) WRITE (LUN,100) I, BLO, BHI, SIGIN
  100   FORMAT (//1X,'SIGBP -- INTERVAL',I4,', BL = ',D10.3,
     .          ', BU = ',D10.3,', SIGIN = ',D15.8)
        IF (BLO .GE. 0.D0  .OR.  BHI .LE. 0.D0) GO TO 9
        IF (SIGIN .GE. SBIG) GO TO 7
C
C Initialize SIG to 0 and test for a null constraint.
C
        SIG = 0.D0
        IF (BLO .LE. -BMX  .AND.  BHI .GE. BMX) GO TO 6
C
C Test for SIG = 0 sufficient.
C
C   The signed orthogonal distance is d(b) = b*(1-b)*
C     (b*V1 - (1-b)*V2), where b = (t2-t)/(t2-t1),
C     V1 = (C2-C1) X CP(1), and V2 = (C2-C1) X CP(2).
C
        DX = X(IP1) - X(I)
        DY = Y(IP1) - Y(I)
        V1 = DX*YP(I) - DY*XP(I)
        V2 = DX*YP(IP1) - DY*XP(IP1)
C
C   Set DP and DM to the maximum and minimum values of d(b)
C     for b in [0,1].  Note that DP .GE. 0 and DM .LE. 0.
C
        S = V1 + V2
        IF (S .EQ. 0.D0) THEN
C
C   The derivative d'(b) is zero at the midpoint b = .5.
C
          IF (V1 .GE. 0.D0) THEN
            DP = V1/4.D0
            DM = 0.D0
          ELSE
            DP = 0.D0
            DM = V1/4.D0
          ENDIF
        ELSE
C
C   Set RP/RM to the roots of the quadratic equation d'(b) =
C     (B0 +/- SQRT(D0))/(3*S) = V2/(B0 -/+ SQRT(D0)) = 0,
C     where B0 = V1 + 2*V2 and D0 = V1**2 + V1*V2 + V2**2.
C     The expression is chosen to avoid cancellation error.
C
          B0 = S + V2
          D0 = S*S - V1*V2
          T = B0 + SIGN(SQRT(D0),B0)
          IF (B0 .GE. 0.D0) THEN
            RP = T/(3.D0*S)
            RM = V2/T
          ELSE
            RP = V2/T
            RM = T/(3.D0*S)
          ENDIF
          IF (V1 .LE. 0.D0  .AND.  V2 .GE. 0.D0) THEN
C
C   The maximum is DP = 0 at the endpoints.
C
            DP = 0.D0
          ELSE
            DP = RP*(1.D0-RP)*(RP*S - V2)
          ENDIF
          IF (V1 .GE. 0.D0  .AND.  V2 .LE. 0.D0) THEN
C
C   The minimum is DM = 0 at the endpoints.
C
            DM = 0.D0
          ELSE
            DM = RM*(1.D0-RM)*(RM*S - V2)
          ENDIF
        ENDIF
C
C   SIG = 0 is sufficient to satisfy the constraints iff
C     DP .LE. BHI and DM .GE. BLO iff F0 .GE. 0.
C
        F0 = MIN(BHI-DP,DM-BLO)
        IF (F0 .GE. 0.D0) GO TO 6
C
C Find a zero of F(SIG) = min(BHI-DP,DM-BLO), where DP and
C   DM are the maximum and minimum values of d(b).  F is an
C   increasing function, F(0) = F0 < 0, and F = FMAX =
C   min(BHI,-BLO) for SIG sufficiently large.  Note that F
C   has a discontinuity in its first derivative if the
C   curves BHI-DP and DM-BLO (as functions of SIG) inter-
C   sect, and the rate of convergence of the zero finder is
C   reduced to linear if such an intersection occurs near
C   the zero of F.
C
C Initialize parameters for the secant method.  The method
C   uses three points:  (SG0,F0), (SIG,F), and (SNEG,FNEG),
C   where SG0 and SNEG are defined implicitly by DSIG = SIG
C   - SG0 and DMAX = SIG - SNEG.  SG0 is initially zero and
C   SNEG is initialized to a sufficiently large value that
C   FNEG > 0.  This value is used only if the initial value
C   of F is negative.
C
        T = MIN(BHI,-BLO)
        FMAX = MAX(1.D-3,T)
        SIG = MAX(ABS(V1),ABS(V2))/T
        DMAX = SIG*(1.D0-T/FMAX)
        SNEG = SIG - DMAX
        IF (LUN .GE. 0) WRITE (LUN,110) F0, FMAX, SNEG
  110   FORMAT (//1X,'F(0) = ',D15.8,', FMAX = ',D15.8,
     .          ', SNEG = ',D15.8/)
        DSIG = SIG
        FNEG = FMAX
        V2M1 = V2 - V1
        NIT = 0
C
C Top of loop:  compute F.
C
    2   EMS = EXP(-SIG)
        IF (SIG .LE. .5D0) THEN
C
C   SIG .LE. .5:  use approximations designed to avoid can-
C                 cellation error (associated with small
C                 SIG) in the modified hyperbolic functions.
C
          CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
          SINH = SINHM + SIG
          A1 = SIG*COSHM*V2 - SINHM*V2M1
          A2 = SIG*SINH*V2 - COSHM*V2M1
          A = A2 - A1
          AA = A/EMS
          E = SIG*SINHM - COSHMM - COSHMM
        ELSE
C
C   SIG > .5:  scale SINHM and COSHM by 2*EXP(-SIG) in order
C              to avoid overflow.
C
          TM = 1.D0 - EMS
          SINH = TM*(1.D0+EMS)
          SINHM = SINH - 2.D0*SIG*EMS
          COSHM = TM*TM
          A1 = SIG*COSHM*V2 - SINHM*V2M1
          A2 = SIG*SINH*V2 - COSHM*V2M1
          AA = 2.D0*(SIG*TM*V2 + (TM-SIG)*V2M1)
          A = EMS*AA
          E = SIG*SINH - COSHM - COSHM
        ENDIF
        IF (S .EQ. 0.D0) THEN
C
C   The derivative d'(b) is zero at the midpoint b = .5.
C
          EB = SIG*COSHM - SINHM - SINHM
          IF (V1 .GE. 0.D0) THEN
            DP = E*V1/(SIG*(SQRT(EB*EB-E*E)+EB))
            DM = 0.D0
          ELSE
            DP = 0.D0
            DM = E*V1/(SIG*(SQRT(EB*EB-E*E)+EB))
          ENDIF
        ELSE
C
C   d'(b)*DC = V2 - (A1*sinh(SIG*b) - A2*coshm(SIG*b))/E = 0
C     for ESB = (-B +/- sqrt(D))/A = C/(-B -/+ sqrt(D)),
C     where ESB = exp(SIG*b), A = A2-A1, D = B**2 - A*C, and
C     B and C are defined below.
C
          B = -COSHM*S
          C = A2 + A1
          D = B*B - A*C
          F = FMAX
          IF (D .LT. 0.D0) GO TO 3
          T1 = SQRT(D)
          T = -B - SIGN(T1,B)
C
          RSP = 0.D0
          IF (B .LT. 0.D0  .AND.  AA .NE. 0.D0) THEN
            IF (T/AA .GT. 0.D0) RSP = SIG + LOG(T/AA)
          ENDIF
          IF ((B .GT. 0.D0  .OR.  AA .EQ. 0.D0)  .AND.
     .        C/T .GT. 0.D0) RSP = LOG(C/T)
          IF ((RSP .LE. 0.D0  .OR.  RSP .GE. SIG)  .AND.
     .        B .NE. 0.D0) THEN
C
C   The maximum is DP = 0 at the endpoints.
C
            DP = 0.D0
          ELSE
            DP = -(B*RSP+A1+T1)/(SIG*E)
          ENDIF
C
          RSM = 0.D0
          IF (B .GT. 0.D0  .AND.  AA .NE. 0.D0) THEN
            IF (T/AA .GT. 0.D0) RSM = SIG + LOG(T/AA)
          ENDIF
          IF ((B .LT. 0.D0  .OR.  AA .EQ. 0.D0)  .AND.
     .        C/T .GT. 0.D0) RSM = LOG(C/T)
          IF ((RSM .LE. 0.D0  .OR.  RSM .GE. SIG)  .AND.
     .        B .NE. 0.D0) THEN
C
C   The minimum is DM = 0 at the endpoints.
C
            DM = 0.D0
          ELSE
            DM = -(B*RSM+A1-T1)/(SIG*E)
          ENDIF
        ENDIF
C
        F = MIN(BHI-DP,DM-BLO)
C
C   Update the number of iterations NIT.
C
    3   NIT = NIT + 1
        IF (LUN .GE. 0) WRITE (LUN,120) NIT, SIG, F
  120   FORMAT (1X,3X,I2,' -- SIG = ',D15.8,', F = ',
     .          D15.8)
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
     .      F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 6
C
C   Test for F0 = F = FMAX or F < 0 on the first iteration.
C
        IF (F0 .NE. F  .AND.  (NIT .GT. 1  .OR.
     .                         F .GT. 0.D0))    GO TO 5
C
C   F*F0 > 0 and either the new estimate would be outside of
C     the bracketing interval of length abs(DMAX) or F < 0
C     on the first iteration.  Reset (SG0,F0) to (SNEG,FNEG).
C
    4   DSIG = DMAX
        F0 = FNEG
C
C   Compute the change in SIG by linear interpolation be-
C     tween (SG0,F0) and (SIG,F).
C
    5   DSIG = -F*DSIG/(F-F0)
        IF (LUN .GE. 0) WRITE (LUN,130) DSIG
  130   FORMAT (1X,8X,'DSIG = ',D15.8)
        IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.
     .       DSIG*DMAX .GT. 0. ) GO TO 4
C
C   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
C     Note that DSIG and DMAX have opposite signs.
C
        IF (ABS(DSIG) .LT. STOL/2.D0)
     .    DSIG = -SIGN(STOL/2.D0,DMAX)
C
C   Bottom of loop:  update SIG, DMAX, and F0.
C
        SIG = SIG + DSIG
        DMAX = DMAX + DSIG
        F0 = F
        GO TO 2
C
C Bottom of loop on intervals:  update SIGMA(I), ICNT, and
C   DSM if necessary.
C
    6   SIG = MIN(SIG,SBIG)
        IF (SIG .GT. SIGIN) THEN
          SIGMA(I) = SIG
          ICNT = ICNT + 1
          DSIG = SIG-SIGIN
          IF (SIGIN .GT. 0.D0) DSIG = DSIG/SIGIN
          DSM = MAX(DSM,DSIG)
        ENDIF
    7   CONTINUE
C
C No errors encountered.
C
      DSMAX = DSM
      IER = ICNT
      RETURN
C
C N < 2.
C
    8 DSMAX = 0.D0
      IER = -1
      RETURN
C
C BL(I) .GE. 0 or BU(I) .LE. 0.
C
    9 DSMAX = DSM
      IER = -(IP1)
      RETURN
      END
