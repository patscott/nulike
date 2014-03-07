      SUBROUTINE TSPBI (N,X,Y,NCD,IENDC,PER,B,BMAX,LWK, WK,
     .                  YP, SIGMA,ICFLG,IER)
      INTEGER N, NCD, IENDC, LWK, ICFLG(N), IER
      LOGICAL PER
      DOUBLE PRECISION X(N), Y(N), B(5,N), BMAX, WK(LWK),
     .                 YP(N), SIGMA(N)
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
C factors SIGMA chosen to satisfy user-specified constraints
C (by Subroutine SIGBI).  Refer to Subroutine TSPSI for an
C alternative method of computing tension factors.
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
C             in second derivative continuity.  This re-
C             quires iterating on calls to YPC2 or YPC2P and
C             calls to SIGBI, and generally results in more
C             nonzero tension factors (hence more expensive
C             evaluation).
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
C       B = Array dimensioned 5 by N-1 containing bounds or
C           flags which define the constraints.  For I = 1
C           to N-1, column I defines the constraints associ-
C           ated with interval (X(I),X(I+1)) as follows:
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
C              used as an upper bound in B (or when when
C              its negative is used as a lower bound),
C              specifies that no constraint is to be en-
C              forced.
C
C       LWK = Length of work space WK:
C             LWK GE 2N-2 if NCD = 2 and PER = FALSE
C             LWK GE 3N-3 if NCD = 2 and PER = TRUE
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
C       SIGMA = Array of length .GE. N-1.
C
C       ICFLG = Array of length .GE. N-1.
C
C On output:
C
C       WK = Array containing convergence parameters in the
C            first two locations if IER > 0 (NCD = 2 and
C            no error other than invalid constraints was
C            encountered):
C            WK(1) = Maximum relative change in a component
C                    of YP on the last iteration.
C            WK(2) = Maximum relative change in a component
C                    of SIGMA on the last iteration.
C
C       YP = Array containing derivatives of H at the
C            abscissae.  YP is not altered if -3 < IER < 0,
C            and YP is only partially defined if IER = -4.
C
C       SIGMA = Array containing tension factors for which
C               H(x) satisfies the constraints defined by B.
C               SIGMA(I) is associated with interval (X(I),
C               X(I+1)) for I = 1,...,N-1.  If infinite ten-
C               sion is required in interval I, then
C               SIGMA(I) = 85 (and H is an approximation to
C               the linear interpolant on the interval),
C               and if no constraint is specified in the
C               interval, then SIGMA(I) = 0, and thus H is
C               cubic.  Invalid constraints are treated as
C               null constraints.  SIGMA is not altered if
C               -3 < IER < 0 (unless IENDC is invalid), and
C               SIGMA is the zero vector if IER = -4 or
C               IENDC (if used) is invalid.
C
C       ICFLG = Array of invalid constraint flags associated
C               with intervals.  For I = 1 to N-1, ICFLG(I)
C               is a 5-bit value b5b4b3b2b1, where bK = 1 if
C               and only if constraint K cannot be satis-
C               fied.  Thus, all constraints in interval I
C               are satisfied if and only if ICFLG(I) = 0
C               (and IER .GE. 0).  ICFLG is not altered if
C               IER < 0.
C
C       IER = Error indicator or iteration count:
C             IER = IC .GE. 0 if no errors were encountered
C                      (other than invalid constraints) and
C                      IC calls to SIGBI and IC+1 calls to
C                      YPC1, YPC1P, YPC2 or YPC2P were
C                      employed.  (IC = 0 if NCD = 1).
C             IER = -1 if N, NCD, or IENDC is outside its
C                      valid range.
C             IER = -2 if LWK is too small.
C             IER = -4 if the abscissae X are not strictly
C                      increasing.
C
C Modules required by TSPBI:  ENDSLP, SIG0, SIG1, SIG2,
C                               SIGBI, SNHCSH, STORE,
C                               YPCOEF, YPC1, YPC1P, YPC2,
C                               YPC2P
C
C Intrinsic functions called by TSPBI:  ABS, MAX
C
C***********************************************************
C
      INTEGER I, ICNT, IERR, ITER, MAXIT, NM1, NN
      LOGICAL LOOP2
      DOUBLE PRECISION DSMAX, DYP, DYPTOL, E, STOL, YP1, YPN
C
      DATA STOL/0.D0/,  MAXIT/49/,  DYPTOL/.01D0/
C
C Convergence parameters:
C
C   STOL = Absolute tolerance for SIGBI.
C   MAXIT = Maximum number of YPC2/SIGBI iterations for
C             each loop if NCD = 2.
C   DYPTOL = Bound on the maximum relative change in a
C              component of YP defining convergence of
C              the YPC2/SIGBI iteration when NCD = 2.
C
      NN = N
      NM1 = NN - 1
C
C Test for invalid input parameters N, NCD, or LWK.
C
      IF (NN .LT. 2  .OR.  (PER  .AND.  NN .LT. 3)  .OR.
     .    NCD .LT. 1  .OR.  NCD .GT. 2) GO TO 11
      IF ( NCD .EQ. 2  .AND.  (LWK .LT. 2*NM1  .OR.
     .     (PER  .AND.  LWK .LT. 3*NM1)) ) GO TO 12
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
          CALL YPC1 (NN,X,Y, YP,IERR)
        ELSE
          CALL YPC1P (NN,X,Y, YP,IERR)
        ENDIF
        IF (IERR .NE. 0) GO TO 14
C
C   Compute tension factors.
C
        CALL SIGBI (NN,X,Y,YP,STOL,B,BMAX, SIGMA, ICFLG,
     .              DSMAX,IERR)
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
      LOOP2 = .FALSE.
C
C   Iterate on calls to SIGBI and YPC2 (or YPC2P).  The
C     first N-1 WK locations are used to store the deriva-
C     tive estimates YP from the previous iteration.
C
C   LOOP2 is TRUE iff tension factors are not allowed to
C         decrease between iterations (loop 1 failed to
C         converge with MAXIT iterations).
C   DYP is the maximum relative change in a component of YP.
C   ICNT is the number of tension factors which were altered
C        by SIGBI.
C   DSMAX is the maximum relative change in a component of
C         SIGMA.
C
    2 DO 6 ITER = 1,MAXIT
        DYP = 0.D0
        DO 3 I = 2,NM1
          WK(I) = YP(I)
    3     CONTINUE
        CALL SIGBI (NN,X,Y,YP,STOL,B,BMAX, SIGMA, ICFLG,
     .              DSMAX,ICNT)
        IF (.NOT. PER) THEN
          CALL YPC2 (NN,X,Y,SIGMA,IENDC,IENDC,YP1,YPN,
     .               WK(NN), YP,IERR)
        ELSE
          CALL YPC2P (NN,X,Y,SIGMA,WK(NN), YP,IERR)
        ENDIF
        DO 4 I = 2,NM1
          E = ABS(YP(I)-WK(I))
          IF (WK(I) .NE. 0.D0) E = E/ABS(WK(I))
          DYP = MAX(DYP,E)
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
        LOOP2 = .TRUE.
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
C Abscissae are not strictly increasing.
C
   14 IER = -4
      RETURN
      END
