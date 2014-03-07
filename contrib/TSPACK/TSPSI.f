      SUBROUTINE TSPSI (N,X,Y,NCD,IENDC,PER,UNIFRM,LWK, WK,
     .                  YP,SIGMA, IER)
      INTEGER N, NCD, IENDC, LWK, IER
      LOGICAL PER, UNIFRM
      DOUBLE PRECISION X(N), Y(N), WK(LWK), YP(N), SIGMA(N)
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
C   This subroutine computes a set of parameter values which
C define a Hermite interpolatory tension spline H(x).  The
C parameters consist of knot derivative values YP computed
C by Subroutine YPC1, YPC1P, YPC2, or YPC2P, and tension
C factors SIGMA computed by Subroutine SIGS (unless UNIFRM =
C TRUE).  Alternative methods for computing SIGMA are pro-
C vided by Subroutine TSPBI and Functions SIG0, SIG1, and
C SIG2.
C
C   Refer to Subroutine TSPSS for a means of computing
C parameters which define a smoothing curve rather than an
C interpolatory curve.
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
C           ciated with the abscissae.  H(X(I)) = Y(I) for
C           I = 1,...,N.  If NCD = 1 and PER = TRUE, Y(N)
C           is set to Y(1).
C
C       NCD = Number of continuous derivatives at the knots.
C             NCD = 1 or NCD = 2.  If NCD = 1, the YP values
C             are computed by local monotonicity-constrained
C             quadratic fits.  Otherwise, a linear system is
C             solved for the derivative values which result
C             in second derivative continuity.  Unless
C             UNIFRM = TRUE, this requires iterating on
C             calls to YPC2 or YPC2P and calls to SIGS, and
C             generally results in more nonzero tension
C             factors (hence more expensive evaluation).
C
C       IENDC = End condition indicator for NCD = 2 and PER
C               = FALSE (or dummy parameter otherwise):
C               IENDC = 0 if YP(1) and YP(N) are to be com-
C                         puted by monotonicity-constrained
C                         parabolic fits to the first three
C                         and last three points, respective-
C                         ly.  This is identical to the
C                         values computed by YPC1.
C               IENDC = 1 if the first derivatives of H at
C                         X(1) and X(N) are user-specified
C                         in YP(1) and YP(N), respectively.
C               IENDC = 2 if the second derivatives of H at
C                         X(1) and X(N) are user-specified
C                         in YP(1) and YP(N), respectively.
C               IENDC = 3 if the end conditions are to be
C                         computed by Subroutine ENDSLP and
C                         vary with SIGMA(1) and SIGMA(N-1).
C
C       PER = Logical variable with value TRUE if and only
C             H(x) is to be a periodic function with period
C             X(N)-X(1).  It is assumed without a test that
C             Y(N) = Y(1) in this case.  On output, YP(N) =
C             YP(1).  If H(x) is one of the components of a
C             parametric curve, this option may be used to
C             obtained a closed curve.
C
C       UNIFRM = Logical variable with value TRUE if and
C                only if constant (uniform) tension is to be
C                used.  The tension factor must be input in
C                SIGMA(1) in this case and must be in the
C                range 0 to 85.  If SIGMA(1) = 0, H(x) is
C                piecewise cubic (a cubic spline if NCD =
C                2), and as SIGMA increases, H(x) approaches
C                the piecewise linear interpolant.  If
C                UNIFRM = FALSE, tension factors are chosen
C                (by SIGS) to preserve local monotonicity
C                and convexity of the data.  This often
C                improves the appearance of the curve over
C                the piecewise cubic fit.
C
C       LWK = Length of work space WK:  no work space is
C             needed if NCD = 1; at least N-1 locations
C             are required if NCD = 2; another N-1 locations
C             are required if PER = TRUE; and an additional
C             N-1 locations are required for the convergence
C             test if SIGS is called (UNIFRM = FALSE):
C
C             LWK GE 0    if NCD=1
C             LWK GE N-1  if NCD=2, PER=FALSE, UNIFRM=TRUE
C             LWK GE 2N-2 if NCD=2, PER=TRUE,  UNIFRM=TRUE
C             LWK GE 2N-2 if NCD=2, PER=FALSE, UNIFRM=FALSE
C             LWK GE 3N-3 if NCD=2, PER=TRUE,  UNIFRM=FALSE
C
C   The above parameters, except possibly Y(N), are not
C altered by this routine.
C
C       WK = Array of length at least LWK to be used as
C            temporary work space.
C
C       YP = Array of length .GE. N containing end condition
C            values in positions 1 and N if NCD = 2 and
C            IENDC = 1 or IENDC = 2.
C
C       SIGMA = Array of length .GE. N-1 containing a ten-
C               sion factor (0 to 85) in the first position
C               if UNIFRM = TRUE.
C
C On output:
C
C       WK = Array containing convergence parameters in the
C            first two locations if IER > 0 (NCD = 2 and
C            UNIFRM = FALSE):
C            WK(1) = Maximum relative change in a component
C                    of YP on the last iteration.
C            WK(2) = Maximum relative change in a component
C                    of SIGMA on the last iteration.
C
C       YP = Array containing derivatives of H at the
C            abscissae.  YP is not altered if -4 < IER < 0,
C            and YP is only partially defined if IER = -4.
C
C       SIGMA = Array containing tension factors.  SIGMA(I)
C               is associated with interval (X(I),X(I+1))
C               for I = 1,...,N-1.  SIGMA is not altered if
C               -4 < IER < 0 (unless IENDC is invalid), and
C               SIGMA is constant (not optimal) if IER = -4
C               or IENDC (if used) is invalid.
C
C       IER = Error indicator or iteration count:
C             IER = IC .GE. 0 if no errors were encountered
C                      and IC calls to SIGS and IC+1 calls
C                      to YPC1, YPC1P, YPC2 or YPC2P were
C                      employed.  (IC = 0 if NCD = 1).
C             IER = -1 if N, NCD, or IENDC is outside its
C                      valid range.
C             IER = -2 if LWK is too small.
C             IER = -3 if UNIFRM = TRUE and SIGMA(1) is out-
C                      side its valid range.
C             IER = -4 if the abscissae X are not strictly
C                      increasing.
C
C Modules required by TSPSI:  ENDSLP, SIGS, SNHCSH, STORE,
C                               YPCOEF, YPC1, YPC1P, YPC2,
C                               YPC2P
C
C Intrinsic functions called by TSPSI:  ABS, MAX
C
C***********************************************************
C
      INTEGER I, ICNT, IERR, ITER, MAXIT, NM1, NN
      DOUBLE PRECISION DSMAX, DYP, DYPTOL, E, SBIG, SIG,
     .                 STOL, YP1, YPN
C
      DATA SBIG/85.D0/
C
C Convergence parameters:
C
C   STOL = Absolute tolerance for SIGS
C   MAXIT = Maximum number of YPC2/SIGS iterations
C   DYPTOL = Bound on the maximum relative change in a
C              component of YP defining convergence of
C              the YPC2/SIGS iteration when NCD = 2 and
C              UNIFRM = FALSE
C
      DATA STOL/0.D0/,  MAXIT/99/,  DYPTOL/.01D0/
C
C Test for invalid input parameters (other than X and
C   IENDC).
C
      NN = N
      NM1 = NN - 1
      IF (NN .LT. 2  .OR.  (PER  .AND.  NN .LT. 3)  .OR.
     .    NCD .LT. 1  .OR.  NCD .GT. 2) GO TO 11
      IF (UNIFRM) THEN
        IF ( NCD .EQ. 2  .AND.  (LWK .LT. NM1  .OR.
     .       (PER  .AND.  LWK .LT. 2*NM1)) ) GO TO 12
        SIG = SIGMA(1)
        IF (SIG .LT. 0.D0  .OR.  SIG .GT. SBIG) GO TO 13
      ELSE
        IF ( NCD .EQ. 2  .AND.  (LWK .LT. 2*NM1  .OR.
     .       (PER  .AND.  LWK .LT. 3*NM1)) ) GO TO 12
        SIG = 0.D0
      ENDIF
C
C Initialize iteration count ITER, and store uniform
C   tension factors, or initialize SIGMA to zeros.
C
      ITER = 0
      DO 1 I = 1,NM1
        SIGMA(I) = SIG
    1   CONTINUE
      IF (NCD .EQ. 1) THEN
C
C NCD = 1.
C
        IF (.NOT. PER) THEN
          CALL YPC1 (NN,X,Y, YP,IERR)
        ELSE
          CALL YPC1P (NN,X,Y, YP,IERR)
        ENDIF
        IF (IERR .NE. 0) GO TO 14
        IF (.NOT. UNIFRM) THEN
C
C   Call SIGS for UNIFRM = FALSE.
C
          CALL SIGS (NN,X,Y,YP,STOL, SIGMA, DSMAX,IERR)
          IF (IERR .LT. 0) GO TO 15
        ENDIF
        GO TO 10
      ENDIF
C
C NCD = 2.
C
      IF (.NOT. PER) THEN
C
C   Nonperiodic case:  call YPC2 and test for IENDC or X
C     invalid.
C
        YP1 = YP(1)
        YPN = YP(NN)
        CALL YPC2 (NN,X,Y,SIGMA,IENDC,IENDC,YP1,YPN,
     .             WK, YP,IERR)
        IF (IERR .EQ. 1) GO TO 11
        IF (IERR .GT. 1) GO TO 14
      ELSE
C
C   Periodic fit:  call YPC2P.
C
        CALL YPC2P (NN,X,Y,SIGMA,WK, YP,IERR)
        IF (IERR .GT. 1) GO TO 14
      ENDIF
      IF (UNIFRM) GO TO 10
C
C   Iterate on calls to SIGS and YPC2 (or YPC2P).  The first
C     N-1 WK locations are used to store the derivative
C     estimates YP from the previous iteration.
C
C   DYP is the maximum relative change in a component of YP.
C   ICNT is the number of tension factors which were
C        increased by SIGS.
C   DSMAX is the maximum relative change in a component of
C         SIGMA.
C
      DO 4 ITER = 1,MAXIT
        DYP = 0.D0
        DO 2 I = 2,NM1
          WK(I) = YP(I)
    2     CONTINUE
        CALL SIGS (NN,X,Y,YP,STOL, SIGMA, DSMAX,ICNT)
        IF (ICNT .LT. 0) GO TO 15
        IF (.NOT. PER) THEN
          CALL YPC2 (NN,X,Y,SIGMA,IENDC,IENDC,YP1,YPN,
     .               WK(NN), YP,IERR)
        ELSE
          CALL YPC2P (NN,X,Y,SIGMA,WK(NN), YP,IERR)
        ENDIF
        DO 3 I = 2,NM1
          E = ABS(YP(I)-WK(I))
          IF (WK(I) .NE. 0.D0) E = E/ABS(WK(I))
          DYP = MAX(DYP,E)
    3     CONTINUE
        IF (ICNT .EQ. 0  .OR.  DYP .LE. DYPTOL) GO TO 5
    4   CONTINUE
      ITER = MAXIT
C
C Store convergence parameters in WK.
C
    5 WK(1) = DYP
      WK(2) = DSMAX
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
C UNIFRM = TRUE and SIGMA(1) outside its valid range.
C
   13 IER = -3
      RETURN
C
C Abscissae are not strictly increasing.
C
   14 IER = -4
      RETURN

C Added by Pat Scott Jan 31 2008
C SIGS didn't exit correctly (probably didn't converge).
C
   15 IER = -5
      RETURN
      END

