      SUBROUTINE TSPSP (N,ND,X,Y,Z,NCD,IENDC,PER,UNIFRM,
     .                  LWK, WK, T,XP,YP,ZP,SIGMA,IER)
      INTEGER N, ND, NCD, IENDC, LWK, IER
      LOGICAL PER, UNIFRM
      DOUBLE PRECISION X(N), Y(N), Z(N), WK(LWK), T(N),
     .                 XP(N), YP(N), ZP(N), SIGMA(N)
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
C parametric planar curve C(t) = (H1(t),H2(t)) or space
C curve C(t) = (H1(t),H2(t),H3(t)) whose components are Her-
C mite interpolatory tension splines.  The output values
C consist of parameters (knots) T computed by ARCL2D or
C ARCL3D, knot derivative values XP, YP, (and ZP) computed
C by Subroutine YPC1, YPC1P, YPC2, or YPC2P, and tension
C factors SIGMA computed by Subroutine SIGS (unless UNIFRM =
C TRUE).
C
C   Refer to Subroutine TSPSP for an alternative method of
C computing tension factors in the case of a planar curve.
C
C   The tension splines may be evaluated by Subroutine
C TSVAL2 (or TSVAL3) or Functions HVAL (values), HPVAL
C (first derivatives), HPPVAL (second derivatives), and
C TSINTL (integrals).
C
C On input:
C
C       N = Number of knots and data points.  N .GE. 2 and
C           N .GE. 3 if PER = TRUE.
C
C       ND = Number of dimensions:
C            ND = 2 if a planar curve is to be constructed.
C            ND = 3 if a space curve is to be constructed.
C
C       X,Y,Z = Arrays of length N containing the Cartesian
C               coordinates of an ordered sequence of data
C               points C(I), I = 1 to N, such that C(I) .NE.
C               C(I+1).  C(t) is constrained to pass through
C               these points.  Z is an unused dummy parame-
C               ter if ND = 2.  In the case of a closed curve
C               (PER = TRUE), the first and last points
C               should coincide.  In this case, X(N), Y(N),
C               (and Z(N)) are set to X(1), Y(1), (and Z(1))
C               if NCD = 1, but not altered if NCD = 2.
C
C       NCD = Number of continuous derivatives at the knots.
C             NCD = 1 or NCD = 2.  If NCD = 1, XP, YP, (and
C             ZP) are computed by local monotonicity-
C             constrained quadratic fits.  Otherwise, a
C             linear system is solved for the derivative
C             values which result in second derivative con-
C             tinuity.  Unless UNIFRM = FALSE, this requires
C             iterating on calls to YPC2 or YPC2P and calls
C             to SIGS, and generally results in more nonzero
C             tension factors (hence more expensive evalua-
C             tion).
C
C       IENDC = End condition indicator for NCD = 2 and PER
C               = FALSE (or dummy parameter otherwise):
C               IENDC = 0 if XP(1), XP(N), YP(1), YP(N) (and
C                         ZP(1) and ZP(N)) are to be com-
C                         puted by monotonicity-constrained
C                         parabolic fits (YPC1).
C               IENDC = 1 if the first derivatives of H1 at
C                         the left and right endpoints are
C                         user-specified in XP(1) and XP(N),
C                         respectively, the first deriva-
C                         tives of H2 at the ends are
C                         specified in YP(1) and YP(N), and,
C                         if ND = 3, the first derivatives
C                         of H3 are specified in ZP(1) and
C                         ZP(N).
C               IENDC = 2 if the second derivatives of H1,
C                         H2, (and H3) at the endpoints are
C                         user-specified in XP(1), XP(N),
C                         YP(1), YP(N), (ZP(1), and ZP(N)).
C               IENDC = 3 if the end conditions are to be
C                         computed by Subroutine ENDSLP and
C                         vary with SIGMA(1) and SIGMA(N-1).
C
C       PER = Logical variable with value TRUE if and only
C             a closed curve is to be constructed -- H1(t),
C             H2(t), (and H3(t)) are to be periodic func-
C             tions with period T(N)-T(1), where T(1) and
C             T(N) are the parameter values associated with
C             the first and last data points.  It is assumed
C             in this case that X(N) = X(1), Y(N) = Y(1)
C             and, if ND = 3, Z(N) = Z(1), and, on output,
C             XP(N) = XP(1), YP(N) = YP(1), (and ZP(N) =
C             ZP(1) if ND = 3).
C
C       UNIFRM = Logical variable with value TRUE if and
C                only if constant (uniform) tension is to be
C                used.  The tension factor must be input in
C                SIGMA(1) in this case and must be in the
C                range 0 to 85.  If SIGMA(1) = 0, H(t) is
C                piecewise cubic (a cubic spline if NCD =
C                2), and as SIGMA increases, H(t) approaches
C                the piecewise linear interpolant, where H
C                is H1, H2, or H3.  If UNIFRM = FALSE,
C                tension factors are chosen (by SIGS) to
C                preserve local monotonicity and convexity
C                of the data.  This often improves the
C                appearance of the curve over the piecewise
C                cubic fitting functions.
C
C       LWK = Length of work space WK:  no work space is
C             needed if NCD = 1; at least N-1 locations
C             are required if NCD = 2; another N-1 locations
C             are required if PER = TRUE; and an additional
C             ND*(N-1) locations are required for the con-
C             vergence test if SIGS is called (UNIFRM =
C             FALSE):
C               If NCD=1 then LWK = 0 (not tested).
C               If NCD=2 then
C
C             LWK GE N-1          if PER=FALSE, UNIFRM=TRUE
C             LWK GE 2N-2         if PER=TRUE,  UNIFRM=TRUE
C             LWK GE (ND+1)*(N-1) if PER=FALSE, UNIFRM=FALSE
C             LWK GE (ND+2)*(N-1) if PER=TRUE,  UNIFRM=FALSE
C
C   The above parameters, except possibly X(N), Y(N), and
C Z(N), are not altered by this routine.
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
C       ZP = Dummy argument if ND=2, or, if ND=3, array of
C            length .GE. N containing end condition values
C            in positions 1 and N if NCD = 2 and IENDC = 1
C            or IENDC = 2.
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
C                    of XP, YP, or ZP on the last iteration.
C            WK(2) = Maximum relative change in a component
C                    of SIGMA on the last iteration.
C
C       T = Array containing parameter values computed by
C           Subroutine ARCL2D or ARCL3D unless -4 < IER < 0.
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
C       ZP = Array containing derivatives of H3 at the knots
C            if ND=3.  ZP is not altered if -5 < IER < 0,
C            and ZP is only partially defined if IER = -6.
C
C       SIGMA = Array containing tension factors.  SIGMA(I)
C               is associated with interval (T(I),T(I+1))
C               for I = 1,...,N-1.  SIGMA is not altered if
C               -5 < IER < 0 (unless IENDC is invalid), and
C               SIGMA is constant (not optimal) if IER = -6
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
C             IER = -4 if a pair of adjacent data points
C                      coincide:  X(I) = X(I+1), Y(I) =
C                      Y(I+1), (and Z(I) = Z(I+1)) for some
C                      I in the range 1 to N-1.
C             IER = -6 if invalid knots T were returned by
C                      ARCL2D or ARCL3D.  This should not
C                      occur.
C
C Modules required by TSPSP:  ARCL2D, ARCL3D, ENDSLP, SIGS,
C                               SNHCSH, STORE, YPCOEF, YPC1,
C                               YPC1P, YPC2, YPC2P
C
C Intrinsic functions called by TSPSP:  ABS, MAX
C
C***********************************************************
C
      INTEGER I, ICNT, IERR, ITER, IW1, MAXIT, N2M2, NM1, NN
      LOGICAL SCURV
      DOUBLE PRECISION DSMAX, DYP, DYPTOL, EX, EY, EZ, SBIG,
     .                 SIG, STOL, XP1, XPN, YP1, YPN, ZP1,
     .                 ZPN
C
      DATA SBIG/85.D0/
C
C Convergence parameters:
C
C   STOL = Absolute tolerance for SIGS
C   MAXIT = Maximum number of YPC2/SIGS iterations
C   DYPTOL = Bound on the maximum relative change in a com-
C              ponent of XP, YP, or ZP defining convergence
C              of the YPC2/SIGS iteration when NCD = 2 and
C              UNIFRM = FALSE
C
      DATA STOL/0.D0/,  MAXIT/99/,  DYPTOL/.01D0/
C
C Test for invalid input parameters N, NCD, or LWK.
C
      NN = N
      NM1 = NN - 1
      N2M2 = 2*NM1
      IF (NN .LT. 2  .OR.  (PER  .AND.  NN .LT. 3)  .OR.
     .    NCD .LT. 1  .OR.  NCD .GT. 2) GO TO 11
      IF (UNIFRM) THEN
        IF ( NCD .EQ. 2  .AND.  (LWK .LT. NM1  .OR.
     .       (PER  .AND.  LWK .LT. N2M2)) ) GO TO 12
        SIG = SIGMA(1)
        IF (SIG .LT. 0.D0  .OR.  SIG .GT. SBIG) GO TO 13
      ELSE
        IF ( NCD .EQ. 2  .AND.  (LWK .LT. (ND+1)*NM1  .OR.
     .       (PER  .AND.  LWK .LT. (ND+2)*NM1)) ) GO TO 12
        SIG = 0.D0
      ENDIF
C
C Compute the sequence of parameter values T.
C
      SCURV = ND .EQ. 3
      IF (.NOT. SCURV) THEN
        CALL ARCL2D (NN,X,Y, T,IERR)
      ELSE
        CALL ARCL3D (NN,X,Y,Z, T,IERR)
      ENDIF
      IF (IERR .GT. 0) GO TO 14
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
          CALL YPC1 (NN,T,X, XP,IERR)
          CALL YPC1 (NN,T,Y, YP,IERR)
          IF (SCURV) CALL YPC1 (NN,T,Z, ZP,IERR)
        ELSE
          CALL YPC1P (NN,T,X, XP,IERR)
          CALL YPC1P (NN,T,Y, YP,IERR)
          IF (SCURV) CALL YPC1P (NN,T,Z, ZP,IERR)
        ENDIF
        IF (IERR .NE. 0) GO TO 16
        IF (.NOT. UNIFRM) THEN
C
C   Call SIGS for UNIFRM = FALSE.
C
          CALL SIGS (NN,T,X,XP,STOL, SIGMA, DSMAX,IERR)
          CALL SIGS (NN,T,Y,YP,STOL, SIGMA, DSMAX,IERR)
          IF (SCURV) CALL SIGS (NN,T,Z,ZP,STOL, SIGMA,
     .                           DSMAX,IERR)
        ENDIF
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
        IF (SCURV) THEN
          ZP1 = ZP(1)
          ZPN = ZP(NN)
          CALL YPC2 (NN,T,Z,SIGMA,IENDC,IENDC,ZP1,ZPN,
     .               WK, ZP,IERR)
        ENDIF
        IF (IERR .EQ. 1) GO TO 11
        IF (IERR .GT. 1) GO TO 16
      ELSE
C
C   Periodic fit:  call YPC2P.
C
        CALL YPC2P (NN,T,X,SIGMA,WK, XP,IERR)
        CALL YPC2P (NN,T,Y,SIGMA,WK, YP,IERR)
        IF (SCURV) CALL YPC2P (NN,T,Z,SIGMA,WK, ZP,IERR)
        IF (IERR .NE. 0) GO TO 16
      ENDIF
      IF (UNIFRM) GO TO 10
C
C   Iterate on calls to SIGS and YPC2 (or YPC2P).  The
C     first ND*(N-1) WK locations are used to store the
C     derivative estimates XP, YP, (and ZP) from the
C     previous iteration.  IW1 is the first free WK location
C     following the stored derivatives.
C
C   DYP is the maximum relative change in a component of XP,
C       YP, or ZP.
C   ICNT is the number of tension factors which were
C        increased by SIGS.
C   DSMAX is the maximum relative change in a component of
C         SIGMA.
C
      IW1 = ND*NM1 + 1
      DO 5 ITER = 1,MAXIT
        DYP = 0.D0
        DO 2 I = 2,NM1
          WK(I) = XP(I)
          WK(NM1+I) = YP(I)
    2     CONTINUE
        IF (SCURV) THEN
          DO 3 I = 2,NM1
            WK(N2M2+I) = ZP(I)
    3       CONTINUE
        ENDIF
        CALL SIGS (NN,T,X,XP,STOL, SIGMA, DSMAX,ICNT)
        CALL SIGS (NN,T,Y,YP,STOL, SIGMA, DSMAX,ICNT)
        IF (SCURV) CALL SIGS (NN,T,Z,ZP,STOL, SIGMA, DSMAX,
     .                        ICNT)
        IF (.NOT. PER) THEN
          CALL YPC2 (NN,T,X,SIGMA,IENDC,IENDC,XP1,XPN,
     .               WK(IW1), XP,IERR)
          CALL YPC2 (NN,T,Y,SIGMA,IENDC,IENDC,YP1,YPN,
     .               WK(IW1), YP,IERR)
          IF (SCURV) CALL YPC2 (NN,T,Z,SIGMA,IENDC,IENDC,
     .                          ZP1,ZPN,WK(IW1), ZP,IERR)
        ELSE
          CALL YPC2P (NN,T,X,SIGMA,WK(IW1), XP,IERR)
          CALL YPC2P (NN,T,Y,SIGMA,WK(IW1), YP,IERR)
          IF (SCURV) CALL YPC2P (NN,T,Z,SIGMA,WK(IW1), ZP,
     .                           IERR)
        ENDIF
        EZ = 0.D0
        DO 4 I = 2,NM1
          EX = ABS(XP(I)-WK(I))
          IF (WK(I) .NE. 0.D0) EX = EX/ABS(WK(I))
          EY = ABS(YP(I)-WK(NM1+I))
          IF (WK(NM1+I) .NE. 0.D0) EY = EY/ABS(WK(NM1+I))
          IF (SCURV) THEN
            EZ = ABS(ZP(I)-WK(N2M2+I))
            IF (WK(N2M2+I) .NE. 0.D0)
     .        EZ = EZ/ABS(WK(N2M2+I))
          ENDIF
          DYP = MAX(DYP,EX,EY,EZ)
    4     CONTINUE
        IF (ICNT .EQ. 0  .OR.  DYP .LE. DYPTOL) GO TO 6
    5   CONTINUE
      ITER = MAXIT
C
C Store convergence parameters in WK.
C
    6 WK(1) = DYP
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
C Adjacent duplicate data points encountered.
C
   14 IER = -4
      RETURN
C
C Error flag returned by YPC1, YPC1P, YPC2, or YPC2P:
C   T is not strictly increasing.
C
   16 IER = -6
      RETURN
      END
