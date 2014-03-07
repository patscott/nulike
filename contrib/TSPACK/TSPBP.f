      SUBROUTINE TSPBP (N,X,Y,NCD,IENDC,PER,BL,BU,BMAX,
     .                  LWK, WK, T,XP,YP,SIGMA,IER)
      INTEGER N, NCD, IENDC, LWK, IER
      LOGICAL PER
      DOUBLE PRECISION X(N), Y(N), BL(N), BU(N), BMAX,
     .                 WK(LWK), T(N), XP(N), YP(N), SIGMA(N)
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
C   This subroutine computes a set of values which define a
C parametric planar curve C(t) = (H1(t),H2(t)) whose compo-
C nents are Hermite interpolatory tension splines.  The
C output values consist of parameters (knots) T computed by
C ARCL2D, knot derivative values XP and YP computed by Sub-
C routine YPC1, YPC1P, YPC2, or YPC2P, and tension factors
C SIGMA chosen (by Subroutine SIGBP) to satisfy user-
C specified bounds on the distance between C(t) and the
C polygonal curve associated with the data points (refer to
C BL and BU below).
C
C   Refer to Subroutine TSPSP for an alternative method of
C computing tension factors.
C
C   The tension splines may be evaluated by Subroutine
C TSVAL2 or Functions HVAL (values), HPVAL (first deriva-
C tives), HPPVAL (second derivatives), and TSINTL
C (integrals).
C
C On input:
C
C       N = Number of knots and data points.  N .GE. 2 and
C           N .GE. 3 if PER = TRUE.
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of an ordered sequence of data
C             points C(I), I = 1 to N, such that C(I) .NE.
C             C(I+1).  C(t) is constrained to pass through
C             these points.  In the case of a closed curve
C             (PER = TRUE), the first and last points should
C             coincide.  (X(N) and Y(N) are set to X(1) and
C             Y(1) if NCD = 1, but not altered if NCD = 2,
C             in this case.)
C
C       NCD = Number of continuous derivatives at the knots.
C             NCD = 1 or NCD = 2.  If NCD = 1, XP and YP are
C             computed by local monotonicity-constrained
C             quadratic fits.  Otherwise, a linear system is
C             solved for the derivative values which result
C             in second derivative continuity.  This re-
C             quires iterating on calls to YPC2 or YPC2P and
C             calls to SIGBP, and generally results in more
C             nonzero tension factors (hence more expensive
C             evaluation).
C
C       IENDC = End condition indicator for NCD = 2 and PER
C               = FALSE (or dummy parameter otherwise):
C               IENDC = 0 if XP(1), XP(N), YP(1), and YP(N)
C                         are to be computed by monotonicity-
C                         constrained parabolic fits (YPC1).
C               IENDC = 1 if the first derivatives of H1 at
C                         the left and right endpoints are
C                         user-specified in XP(1) and XP(N),
C                         respectively, and the first deriv-
C                         atives of H2 at the ends are
C                         specified in YP(1) and YP(N).
C               IENDC = 2 if the second derivatives of H1
C                         and H2 at the endpoints are user-
C                         specified in XP(1), XP(N), YP(1),
C                         and YP(N).
C               IENDC = 3 if the end conditions are to be
C                         computed by Subroutine ENDSLP and
C                         vary with SIGMA(1) and SIGMA(N-1).
C
C       PER = Logical variable with value TRUE if and only
C             a closed curve is to be constructed -- H1(t)
C             and H2(t) are to be periodic functions with
C             period T(N)-T(1), where T(1) and T(N) are the
C             parameter values associated with the first and
C             last data points.  It is assumed that X(N) =
C             X(1) and Y(N) = Y(1) in this case, and, on
C             output, XP(N) = XP(1) and YP(N) = YP(1).
C
C       BL,BU = Arrays of length N-1 containing (for each
C               knot subinterval) lower and upper bounds,
C               respectively, on the signed perpendicular
C               distance d(t) = (C2-C1)/DC X (C(t)-C1),
C               where C1 and C2 are the ordered data points
C               associated with the interval, and DC is the
C               interval length (and length of the line seg-
C               ment C1-C2).  Note that d(t) > 0 iff C(t)
C               lies strictly to the left of the line seg-
C               ment as viewed from C1 toward C2.  For I =
C               1 to N-1, SIGMA(I) is chosen to be as small
C               as possible within the constraint that
C               BL(I) .LE. d(t) .LE. BU(I) for all t in the
C               interval.  BL(I) < 0 and BU(I) > 0 for I = 1
C               to N-1.  A null constraint is specified by
C               BL(I) .LE. -BMAX or BU(I) .GE. BMAX.
C
C       BMAX = User-defined value of infinity which, when
C              used as an upper bound in BU (or when its
C              negative is used as a lower bound in BL),
C              specifies that no constraint is to be en-
C              forced.
C
C       LWK = Length of work space WK:
C             LWK GE 3N-3 if NCD = 2 and PER = FALSE
C             LWK GE 4N-4 if NCD = 2 and PER = TRUE
C
C   The above parameters, except possibly X(N) and Y(N), are
C not altered by this routine.
C
C       WK = Array of length .GE. LWK to be used as tempor-
C            ary work space.
C
C       T = Array of length .GE. N.
C
C       XP = Array of length .GE. N containing end condition
C            values in positions 1 and N if NCD = 2 and
C            IENDC = 1 or IENDC = 2.
C
C       YP = Array of length .GE. N containing end condition
C            values in positions 1 and N if NCD = 2 and
C            IENDC = 1 or IENDC = 2.
C
C       SIGMA = Array of length .GE. N-1.
C
C On output:
C
C       WK = Array containing convergence parameters in the
C            first two locations if IER > 0 (NCD = 2 and
C            no error was encountered):
C            WK(1) = Maximum relative change in a component
C                    of XP or YP on the last iteration.
C            WK(2) = Maximum relative change in a component
C                    of SIGMA on the last iteration.
C
C       T = Array containing parameter values computed by
C           Subroutine ARCL2D unless IER = -1 or IER = -2.
C           T is only partially defined if IER = -4.
C
C       XP = Array containing derivatives of H1 at the
C            knots.  XP is not altered if -5 < IER < 0,
C            and XP is only partially defined if IER = -6.
C
C       YP = Array containing derivatives of H2 at the
C            knots.  YP is not altered if -5 < IER < 0,
C            and YP is only partially defined if IER = -6.
C
C       SIGMA = Array containing tension factors for which
C               C(t) satisfies the constraints defined by
C               BL and BU.  SIGMA(I) is associated with
C               interval (T(I),T(I+1)) for I = 1,...,N-1.
C               SIGMA(I) is limited to 85 (in which case
C               C(t) is indistinguishable from the line
C               segment associated with the interval), and
C               if no constraint is specified in the
C               interval, then SIGMA(I) = 0, and thus H1 and
C               H2 are cubic functions of t.  SIGMA is not
C               altered if -5 < IER < 0 (unless IENDC in
C               invalid), and SIGMA is the zero vector if
C               IER = -6 or IENDC (if used) is invalid.
C
C       IER = Error indicator or iteration count:
C             IER = IC .GE. 0 if no errors were encountered
C                      and IC calls to SIGBP and IC+1 calls
C                      to YPC1, YPC1P, YPC2 or YPC2P were
C                      employed.  (IC = 0 if NCD = 1).
C             IER = -1 if N, NCD, or IENDC is outside its
C                      valid range.
C             IER = -2 if LWK is too small.
C             IER = -4 if a pair of adjacent data points
C                      coincide:  X(I) = X(I+1) and Y(I) =
C                      Y(I+1) for some I in the range 1 to
C                      N-1.
C             IER = -5 if BL(I) .GE. 0 or BU(I) .LE. 0 for
C                      some I in the range 1 to N-1.
C                      SIGMA(J) = 0 for J .GE. I in this
C                      case.
C             IER = -6 if invalid knots T were returned by
C                      ARCL2D.  This should not occur.
C
C Modules required by TSPBP:  ARCL2D, ENDSLP, SIGBP, SNHCSH,
C                               STORE, YPCOEF, YPC1, YPC1P,
C                               YPC2, YPC2P
C
C Intrinsic functions called by TSPBP:  ABS, MAX
C
C***********************************************************
C
      INTEGER I, ICNT, IERR, ITER, MAXIT, N2M1, NM1, NN
      LOGICAL LOOP2
      DOUBLE PRECISION DSMAX, DYP, DYPTOL, EX, EY, STOL,
     .                 XP1, XPN, YP1, YPN
C
      DATA STOL/0.D0/,  MAXIT/49/,  DYPTOL/.01D0/
C
C Convergence parameters:
C
C   STOL = Absolute tolerance for SIGBP.
C   MAXIT = Maximum number of YPC2/SIGBP iterations for each
C             loop in NCD = 2.
C   DYPTOL = Bound on the maximum relative change in a
C              component of XP or YP defining convergence
C              of the YPC2/SIGBP iteration when NCD = 2.
C
      NN = N
      NM1 = NN - 1
C
C Test for invalid input parameters N, NCD, or LWK.
C
      N2M1 = 2*NN - 1
      IF (NN .LT. 2  .OR.  (PER  .AND.  NN .LT. 3)  .OR.
     .    NCD .LT. 1  .OR.  NCD .GT. 2) GO TO 11
      IF ( NCD .EQ. 2  .AND.  (LWK .LT. 3*NM1  .OR.
     .     (PER  .AND.  LWK .LT. 4*NM1)) ) GO TO 12
C
C Compute the sequence of parameter values T.
C
      CALL ARCL2D (NN,X,Y, T,IERR)
      IF (IERR .GT. 0) GO TO 14
C
C Initialize iteration count ITER, and initialize SIGMA to
C   zeros.
C
      ITER = 0
      DO 1 I = 1,NM1
        SIGMA(I) = 0.D0
    1   CONTINUE
      IF (NCD .EQ. 1) THEN
C
C NCD = 1.
C
        IF (.NOT. PER) THEN
          CALL YPC1 (NN,T,X, XP,IERR)
          CALL YPC1 (NN,T,Y, YP,IERR)
        ELSE
          CALL YPC1P (NN,T,X, XP,IERR)
          CALL YPC1P (NN,T,Y, YP,IERR)
        ENDIF
        IF (IERR .NE. 0) GO TO 16
C
C   Compute tension factors.
C
        CALL SIGBP (NN,X,Y,XP,YP,STOL,BL,BU,
     .              BMAX, SIGMA, DSMAX,IERR)
        IF (IERR .LT. 0) GO TO 15
        GO TO 10
      ENDIF
C
C NCD = 2.
C
      IF (.NOT. PER) THEN
C
C   Nonperiodic case:  call YPC2 and test for IENDC invalid.
C
        XP1 = XP(1)
        XPN = XP(NN)
        CALL YPC2 (NN,T,X,SIGMA,IENDC,IENDC,XP1,XPN,
     .             WK, XP,IERR)
        YP1 = YP(1)
        YPN = YP(NN)
        CALL YPC2 (NN,T,Y,SIGMA,IENDC,IENDC,YP1,YPN,
     .             WK, YP,IERR)
        IF (IERR .EQ. 1) GO TO 11
        IF (IERR .GT. 1) GO TO 16
      ELSE
C
C   Periodic fit:  call YPC2P.
C
        CALL YPC2P (NN,T,X,SIGMA,WK, XP,IERR)
        CALL YPC2P (NN,T,Y,SIGMA,WK, YP,IERR)
        IF (IERR .NE. 0) GO TO 16
      ENDIF
      LOOP2 = .FALSE.
C
C   Iterate on calls to SIGBP and YPC2 (or YPC2P).  The
C     first 2N-2 WK locations are used to store the deriva-
C     tive estimates XP and YP from the previous iteration.
C
C   LOOP2 is TRUE iff tension factors are not allowed to
C         decrease between iterations (loop 1 failed to
C         converge with MAXIT iterations).
C   DYP is the maximum relative change in a component of XP
C       or YP.
C   ICNT is the number of tension factors which were altered
C        by SIGBP.
C   DSMAX is the maximum relative change in a component of
C         SIGMA.
C
    2 DO 6 ITER = 1,MAXIT
        DYP = 0.D0
        DO 3 I = 2,NM1
          WK(I) = XP(I)
          WK(NM1+I) = YP(I)
    3     CONTINUE
        CALL SIGBP (NN,X,Y,XP,YP,STOL,BL,BU,
     .              BMAX, SIGMA, DSMAX,ICNT)
        IF (ICNT .LT. 0) GO TO 15
        IF (.NOT. PER) THEN
          CALL YPC2 (NN,T,X,SIGMA,IENDC,IENDC,XP1,XPN,
     .               WK(N2M1), XP,IERR)
          CALL YPC2 (NN,T,Y,SIGMA,IENDC,IENDC,YP1,YPN,
     .               WK(N2M1), YP,IERR)
        ELSE
          CALL YPC2P (NN,T,X,SIGMA,WK(N2M1), XP,IERR)
          CALL YPC2P (NN,T,Y,SIGMA,WK(N2M1), YP,IERR)
        ENDIF
        DO 4 I = 2,NM1
          EX = ABS(XP(I)-WK(I))
          IF (WK(I) .NE. 0.D0) EX = EX/ABS(WK(I))
          EY = ABS(YP(I)-WK(NM1+I))
          IF (WK(NM1+I) .NE. 0.D0) EY = EY/ABS(WK(NM1+I))
          DYP = MAX(DYP,EX,EY)
    4     CONTINUE
        IF (ICNT .EQ. 0  .OR.  DYP .LE. DYPTOL) GO TO 7
        IF (.NOT. LOOP2) THEN
C
C   Loop 1:  reinitialize SIGMA to zeros.
C
          DO 5 I = 1,NM1
            SIGMA(I) = 0.D0
    5       CONTINUE
        ENDIF
    6   CONTINUE
C
C   The loop failed to converge within MAXIT iterations.
C
      ITER = MAXIT
      IF (.NOT. LOOP2) THEN
        LOOP2 = .FALSE.
        GO TO 2
      ENDIF
C
C Store convergence parameters.
C
    7 WK(1) = DYP
      WK(2) = DSMAX
      IF (LOOP2) ITER = ITER + MAXIT
C
C No error encountered.
C
   10 IER = ITER
      RETURN
C
C Invalid input parameter N, NCD, or IENDC.
C
   11 IER = -1
      RETURN
C
C LWK too small.
C
   12 IER = -2
      RETURN
C
C Adjacent duplicate data points encountered.
C
   14 IER = -4
      RETURN
C
C Invalid constraint encountered by SIGBP.
C
   15 IER = -5
      RETURN
C
C Error flag returned by YPC1, YPC1P, YPC2, or YPC2P:
C   T is not strictly increasing.
C
   16 IER = -6
      RETURN
      END
