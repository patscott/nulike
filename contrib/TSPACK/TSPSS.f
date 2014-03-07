      SUBROUTINE TSPSS (N,X,Y,PER,UNIFRM,W,SM,SMTOL,LWK, WK,
     .                  SIGMA,YS,YP, NIT,IER)
      INTEGER N, LWK, NIT, IER
      LOGICAL PER, UNIFRM
      DOUBLE PRECISION X(N), Y(N), W(N), SM, SMTOL, WK(LWK),
     .                 SIGMA(N), YS(N), YP(N)
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
C   This subroutine computes a set of parameter values which
C define a smoothing tension spline H(x).  The parameters
C consist of knot values YS and derivatives YP computed
C by Subroutine SMCRV, and tension factors SIGMA computed by
C Subroutine SIGS (unless UNIFRM = TRUE).  The Hermite
C interpolatory tension spline H(x) defined by the knot
C values and derivatives has two continuous derivatives and
C satisfies either natural or periodic end conditions.
C
C   The tension spline may be evaluated by Subroutine TSVAL1
C or Functions HVAL (values), HPVAL (first derivatives),
C HPPVAL (second derivatives), and TSINTL (integrals).
C
C On input:
C
C       N = Number of data points.  N .GE. 2 and N .GE. 3 if
C           PER = TRUE.
C
C       X = Array of length N containing a strictly in-
C           creasing sequence of abscissae:  X(I) < X(I+1)
C           for I = 1,...,N-1.
C
C       Y = Array of length N containing data values asso-
C           ciated with the abscissae.  If PER = TRUE, it is
C           assumed that Y(N) = Y(1).
C
C       PER = Logical variable with value TRUE if and only
C             H(x) is to be a periodic function with period
C             X(N)-X(1).  It is assumed without a test that
C             Y(N) = Y(1) in this case.  On output, YP(N) =
C             YP(1) and, more generally, the values and
C             first two derivatives of H at X(1) agree with
C             those at X(N).  If H(x) is one of the compo-
C             nents of a parametric curve, this option may
C             be used to obtained a closed curve.  If PER =
C             FALSE, H satisfies natural end conditions:
C             zero second derivatives at X(1) and X(N).
C
C       UNIFRM = Logical variable with value TRUE if and
C                only if constant (uniform) tension is to be
C                used.  The tension factor must be input in
C                SIGMA(1) in this case and must be in the
C                range 0 to 85.  If SIGMA(1) = 0, H(x) is
C                a cubic spline, and as SIGMA increases,
C                H(x) approaches piecewise linear.  If
C                UNIFRM = FALSE, tension factors are chosen
C                (by SIGS) to preserve local monotonicity
C                and convexity of the data.  This may re-
C                sult in a better fit than the case of
C                uniform tension, but requires an iteration
C                on calls to SMCRV and SIGS.
C
C       W = Array of length N containing positive weights
C           associated with the data values.  The recommend-
C           ed value of W(I) is 1/DY**2, where DY is the
C           standard deviation associated with Y(I).  If
C           nothing is known about the errors in Y, a con-
C           stant (estimated value) should be used for DY.
C           If PER = TRUE, it is assumed that W(N) = W(1).
C
C       SM = Positive parameter specifying an upper bound on
C            Q2(YS), where Q2(YS) is the weighted sum of
C            squares of deviations from the data (differ-
C            ences between YS and Y).  H(x) is linear (and
C            Q2 is minimized) if SM is sufficiently large
C            that the constraint is not active.  It is
C            recommended that SM satisfy N-SQRT(2N) .LE. SM
C            .LE. N+SQRT(2N) and SM = N is reasonable if
C            W(I) = 1/DY**2.
C
C       SMTOL = Parameter in the range (0,1) specifying the
C               relative error allowed in satisfying the
C               constraint:  the constraint is assumed to
C               be satisfied if SM*(1-SMTOL) .LE. Q2 .LE.
C               SM*(1+SMTOL).  A reasonable value for SMTOL
C               is SQRT(2/N).
C
C       LWK = Length of work space WK:
C             LWK .GE. 6N   if PER=FALSE  and  UNIFRM=TRUE
C             LWK .GE. 7N   if PER=FALSE  and  UNIFRM=FALSE
C             LWK .GE. 10N  if PER=TRUE   and  UNIFRM=TRUE
C             LWK .GE. 11N  if PER=TRUE   and  UNIFRM=FALSE
C
C The above parameters are not altered by this routine.
C
C       WK = Array of length at least LWK to be used as
C            temporary work space.
C
C       SIGMA = Array of length .GE. N-1 containing a ten-
C               sion factor (0 to 85) in the first position
C               if UNIFRM = TRUE.
C
C       YS = Array of length .GE. N.
C
C       YP = Array of length .GE. N.
C
C On output:
C
C       WK = Array containing convergence parameters in the
C            first two locations if NIT > 0:
C            WK(1) = Maximum relative change in a component
C                    of YS on the last iteration.
C            WK(2) = Maximum relative change in a component
C                    of SIGMA on the last iteration.
C
C       SIGMA = Array containing tension factors.  SIGMA(I)
C               is associated with interval (X(I),X(I+1))
C               for I = 1,...,N-1.  SIGMA is not altered if
C               N is invalid or -4 < IER < -1, and SIGMA is
C               constant if IER = -1 (and N is valid) or
C               IER = -4.
C
C       YS = Array of length N containing values of H at the
C            abscissae.  YS(N) = YS(1) if PER = TRUE.  YS is
C            not altered if IER < 0.
C
C       YP = Array of length N containing first derivative
C            values of H at the abscissae.  YP(N) = YP(1)
C            if PER = TRUE.  YP is not altered if IER < 0.
C
C       NIT = Number of iterations (calls to SIGS).  NIT = 0
C             if IER < 0 or UNIFRM = TRUE.  If NIT > 0,
C             NIT+1 calls to SMCRV were employed.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and the
C                     constraint is active:  Q2(YS) is ap-
C                     proximately equal to SM.
C             IER = 1 if no errors were encountered but the
C                     constraint is not active:  YS and YP
C                     are the values and derivatives of the
C                     linear function (constant function if
C                     PERIOD = TRUE) which minimizes Q2, and
C                     Q1 = 0 (refer to SMCRV).
C             IER = -1 if N, W, SM, or SMTOL is outside its
C                      valid range.
C             IER = -2 if LWK is too small.
C             IER = -3 if UNIFRM = TRUE and SIGMA(1) is out-
C                      side its valid range.
C             IER = -4 if the abscissae X are not strictly
C                      increasing.
C
C Modules required by TSPSS:  B2TRI or B2TRIP, SIGS, SMCRV,
C                               SNHCSH, STORE, YPCOEF
C
C Intrinsic functions called by TSPSS:  ABS, MAX
C
C***********************************************************
C
      INTEGER I, ICNT, IERR, ITER, MAXIT, NM1, NN
      DOUBLE PRECISION DSMAX, DYS, DYSTOL, E, SBIG, SIG,
     .                 STOL
C
      DATA SBIG/85.D0/
C
C Convergence parameters:
C
C   STOL = Absolute tolerance for SIGS
C   MAXIT = Maximum number of SMCRV/SIGS iterations
C   DYSTOL = Bound on the maximum relative change in a
C              component of YS defining convergence of
C              the SMCRV/SIGS iteration when UNIFRM = FALSE
C
      DATA STOL/0.D0/,  MAXIT/99/,  DYSTOL/.01D0/
C
C Initialize NIT, and test for invalid input parameters LWK
C   and SIGMA(1).
C
      NIT = 0
      NN = N
      NM1 = NN - 1
      IF (NN .LT. 2  .OR.  (PER  .AND.  NN .LT. 3)) GO TO 11
      IF (UNIFRM) THEN
        IF ( LWK .LT. 6*NN  .OR.
     .       (PER  .AND.  LWK .LT. 10*NN) ) GO TO 12
        SIG = SIGMA(1)
        IF (SIG .LT. 0.D0  .OR.  SIG .GT. SBIG) GO TO 13
      ELSE
        IF ( LWK .LT. 7*NN  .OR.
     .       (PER  .AND.  LWK .LT. 11*NN) ) GO TO 12
        SIG = 0.D0
      ENDIF
C
C Store uniform tension factors, or initialize SIGMA to
C   zeros.
C
      DO 1 I = 1,NM1
        SIGMA(I) = SIG
    1   CONTINUE
C
C Compute smoothing curve for uniform tension.
C
      CALL SMCRV (NN,X,Y,SIGMA,PER,W,SM,SMTOL,WK, YS,YP,IER)
      IF (IER .LE. -2) IER = -4
      IF (IER .LT. 0  .OR.  UNIFRM) RETURN
C
C   Iterate on calls to SIGS and SMCRV.  The first N-1 WK
C     locations are used to store the function values YS
C     from the previous iteration.
C
C   DYS is the maximum relative change in a component of YS.
C   ICNT is the number of tension factors which were
C        increased by SIGS.
C   DSMAX is the maximum relative change in a component of
C         SIGMA.
C
      DO 4 ITER = 1,MAXIT
        DYS = 0.D0
        DO 2 I = 2,NM1
          WK(I) = YS(I)
    2     CONTINUE
        CALL SIGS (NN,X,Y,YP,STOL, SIGMA, DSMAX,ICNT)
        CALL SMCRV (NN,X,Y,SIGMA,PER,W,SM,SMTOL,WK(NN), YS,
     .              YP,IERR)
        DO 3 I = 2,NM1
          E = ABS(YS(I)-WK(I))
          IF (WK(I) .NE. 0.D0) E = E/ABS(WK(I))
          DYS = MAX(DYS,E)
    3     CONTINUE
        IF (ICNT .EQ. 0  .OR.  DYS .LE. DYSTOL) GO TO 5
    4   CONTINUE
      ITER = MAXIT
C
C No error encountered.
C
    5 WK(1) = DYS
      WK(2) = DSMAX
      NIT = ITER
      IER = IERR
      RETURN
C
C Invalid input parameter N, W, SM, or SMTOL.
C
   11 IER = -1
      RETURN
C
C LWK too small.
C
   12 IER = -2
      RETURN
C
C UNIFRM = TRUE and SIGMA(1) outside its valid range.
C
   13 IER = -3
      RETURN
      END
