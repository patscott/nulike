      SUBROUTINE SMCRV (N,X,Y,SIGMA,PERIOD,W,SM,SMTOL,
     .                  WK, YS,YP,IER)
      LOGICAL PERIOD
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), SIGMA(N), W(N), SM,
     .                 SMTOL, WK(N,10), YS(N), YP(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   10/05/98
C
C   Given a sequence of abscissae X with associated data
C values Y and tension factors SIGMA, this routine deter-
C mines a set of function values YS and first derivatives YP
C associated with a Hermite interpolatory tension spline
C H(x) which smoothes the data.  H(x) has two continuous
C derivatives for all x and satisfies either natural or per-
C iodic end conditions.  The values and derivatives are
C chosen to minimize a quadratic functional Q1(YS,YP) sub-
C ject to the constraint Q2(YS) .LE. SM for Q2(YS) =
C (Y-YS)**T*W*(Y-YS), where **T denotes transpose and W is a
C diagonal matrix of positive weights.
C
C   Functions HVAL, HPVAL, HPPVAL, and TSINTL may be called
C to compute values, derivatives, and integrals of H.  The
C function values YS must be used as data values in those
C subprograms.
C
C   The smoothing procedure is an extension of the method
C for cubic spline smoothing due to C. Reinsch:  Numer.
C Math., 10 (1967) and 16 (1971).  Q1 is defined as the sum
C of integrals over the intervals (X(I),X(I+1)) of HPP**2 +
C (SIGMA(I)/DX)**2*(HP-S)**2, where DX = X(I+1)-X(I), HP and
C HPP denote first and second derivatives of H, and S =
C (YS(I+1)-YS(I))/DX.  Introducing a smoothing parameter P,
C and assuming the constraint is active, the problem is
C equivalent to minimizing Q(P,YS,YP) = Q1(YS,YP) +
C P*(Q2(YS)-SM).  The secant method is used to find a zero
C of G(P) = 1/SQRT(Q2) - 1/SQRT(SM), where YS and YP satisfy
C the order 2N symmetric positive-definite linear system
C obtained by setting the gradient of Q (treated as a func-
C tion of YS and YP) to zero.
C
C   Note that the interpolation problem corresponding to
C YS = Y, SM = 0, and P infinite is solved by Subroutine
C YPC2 or YPC2P.
C
C On input:
C
C       N = Number of data points.  N .GE. 2 if PERIOD =
C           FALSE, and N .GE. 3 if PERIOD = TRUE.
C
C       X = Array of length N containing a strictly in-
C           creasing sequence of abscissae:  X(I) < X(I+1)
C           for I = 1,...,N-1.
C
C       Y = Array of length N containing data values assoc-
C           iated with the abscissae.  If PERIOD = TRUE, it
C           is assumed that Y(N) = Y(1).
C
C       SIGMA = Array of length N-1 containing tension
C               factors.  SIGMA(I) is associated with inter-
C               val (X(I),X(I+1)) for I = 1,...,N-1.  If
C               SIGMA(I) = 0, H is cubic, and as SIGMA in-
C               creases, H approaches linear in the inter-
C               val.
C
C       PERIOD = Periodic end condition flag:
C                PERIOD = .F. if H is to satisfy natural end
C                             conditions:  zero second der-
C                             ivatives at X(1) and X(N).
C                PERIOD = .T. if H is to satisfy periodic
C                             end conditions:  the values
C                             and first two derivatives at
C                             X(1) agree with those at X(N),
C                             and a period thus has length
C                             X(N)-X(1).
C
C       W = Array of length N containing positive weights
C           associated with the data values.  The recommend-
C           ed value of W(I) is 1/DY**2, where DY is the
C           standard deviation associated with Y(I).  If
C           nothing is known about the errors in Y, a con-
C           stant (estimated value) should be used for DY.
C           If PERIOD = TRUE, it is assumed that W(N) =
C           W(1).
C
C       SM = Positive parameter specifying an upper bound on
C            Q2(YS).  H(x) is linear (and Q2 is minimized)
C            if SM is sufficiently large that the constraint
C            is not active.  It is recommended that SM sat-
C            isfy N-SQRT(2N) .LE. SM .LE. N+SQRT(2N) and
C            SM = N is reasonable if W(I) = 1/DY**2.
C
C       SMTOL = Parameter in the range (0,1) specifying the
C               relative error allowed in satisfying the
C               constraint:  the constraint is assumed to
C               be satisfied if SM*(1-SMTOL) .LE. Q2 .LE.
C               SM*(1+SMTOL).  A reasonable value for SMTOL
C               is SQRT(2/N).
C
C The above parameters are not altered by this routine.
C
C       WK = Work space of length at least 6N if PERIOD =
C            FALSE, and 10N if PERIOD = TRUE.
C
C On output:
C
C       YS = Array of length N containing values of H at the
C            abscissae unless IER < 0.  YS(N) = YS(1) if
C            PERIOD = TRUE.
C
C       YP = Array of length N containing first derivative
C            values of H at the abscissae unless IER < 0.
C            YP(N) = YP(1) if PERIOD = TRUE.
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
C                     Q1 = 0.
C             IER = -1 if N, W, SM, or SMTOL is outside its
C                      valid range.  YS and YP are unaltered
C                      in this case.
C             IER = -I if X(I) .LE. X(I-1) for some I in the
C                      range 2,...,N.  YS and YP are unal-
C                      tered in this case.
C
C Modules required by SMCRV:  B2TRI or B2TRIP, SNHCSH,
C                               YPCOEF
C
C Intrinsic functions called by SMCRV:  ABS, SQRT
C
C***********************************************************
C
      INTEGER I, IERR, ITER, LUN, NM1, NN
      LOGICAL PER
      DOUBLE PRECISION C11, C12, C22, D, DMAX, DP, DX, G,
     .                 G0, GNEG, H0, HP, P, P0, Q2, Q2MAX,
     .                 Q2MIN, R1, R2, S, SD, SIG, WI, WIXI,
     .                 XI, YI
C
      DATA    LUN/-1/
      NN = N
      PER = PERIOD
C
C Test for errors, and compute the components of the system
C   (normal equations) for the weighted least squares linear
C   fit.
C
      IER = -1
      IF (NN .LT. 2  .OR.  (NN .LT. 3  .AND.  PER)  .OR.
     .    SM .LE. 0.D0  .OR.  SMTOL .LE. 0.D0  .OR.
     .    SMTOL .GE. 1.D0) RETURN
      C11 = 0.D0
      C12 = 0.D0
      C22 = 0.D0
      R1 = 0.D0
      R2 = 0.D0
      XI = X(1) - 1.D0
      DO 1 I = 1,NN
        WI = W(I)
        IF (WI .LE. 0.D0) RETURN
        IF (X(I) .LE. XI) THEN
          IER = -I
          RETURN
        ENDIF
        XI = X(I)
        YI = Y(I)
        C22 = C22 + WI
        R2 = R2 + WI*YI
        IF (.NOT. PER) THEN
          WIXI = WI*XI
          C11 = C11 + WIXI*XI
          C12 = C12 + WIXI
          R1 = R1 + WIXI*YI
        ENDIF
    1   CONTINUE
C
C Solve the system for (HP,H0), where HP is the derivative
C   (constant) and H0 = H(0).
C
      IF (PER) THEN
        H0 = R2/C22
        HP = 0.D0
      ELSE
        H0 = (C11*R2-C12*R1)/(C11*C22-C12*C12)
        HP = (R1 - C12*H0)/C11
      ENDIF
C
C Store function values and derivatives, and accumulate
C   Q2 = (Y-YS)**T*W*(Y-YS).
C
      Q2 = 0.D0
      DO 2 I = 1,NN
        YS(I) = HP*X(I) + H0
        YP(I) = HP
        Q2 = Q2 + W(I)*(Y(I)-YS(I))**2
    2   CONTINUE
C
C Compute bounds on Q2 defined by SMTOL, and test for the
C   constraint satisfied by the linear fit.
C
      Q2MIN = SM*(1.D0 - SMTOL)
      Q2MAX = SM*(1.D0 + SMTOL)
      IF (Q2 .LE. Q2MAX) THEN
C
C   The constraint is satisfied by a linear function.
C
        IER = 1
        IF (LUN .GE. 0) WRITE (LUN,100)
  100   FORMAT (///1X,'SMCRV -- THE CONSTRAINT IS NOT ',
     .          'ACTIVE AND THE FIT IS LINEAR.'/)
        RETURN
      ENDIF
C
C Compute the matrix components for the linear system.
C
      IER = 0
      NM1 = NN - 1
      DO 3 I = 1,NM1
        DX = X(I+1) - X(I)
        SIG = ABS(SIGMA(I))
        CALL YPCOEF (SIG,DX, D,SD)
        WK(I,1) = D
        WK(I,2) = SD
    3   CONTINUE
C
C Compute G0 = G(0), and print a heading.
C
      S = 1.D0/SQRT(SM)
      G0 = 1.D0/SQRT(Q2) - S
      IF (LUN .GE. 0) WRITE (LUN,110) SM, SMTOL, G0
  110 FORMAT (///1X,3X,'SMCRV -- SM = ',D10.4,', SMTOL = ',
     .        D14.8,', G(0) = ',D15.8///)
C
C G(P) is strictly increasing and concave, and G(0) < 0.
C
C Initialize parameters for the secant method.  The method
C   uses three points:  (P0,G0), (P,G), and (PNEG,GNEG),
C   where P0 and PNEG are defined implicitly by DP = P - P0
C   and DMAX = P - PNEG.
C
      P = 10.D0*SM
      DP = P
      DMAX = 0.D0
      ITER = 0
C
C Top of loop:  compute G and print a message.  For each
C               secant iteration, the following values are
C               printed:  P, G(P), and DP, where DP is the
C               change in P computed by linear interpolation
C               between the current point (P,G) and a previ-
C               ous point.
C
C
    4 IF (.NOT. PER) THEN
        CALL B2TRI (NN,X,Y,W,P,WK,WK(1,2),WK(1,3),WK(1,4),
     .              WK(1,5),WK(1,6), YS,YP,IERR)
      ELSE
        CALL B2TRIP (NN,X,Y,W,P,WK,WK(1,2),WK(1,3),WK(1,4),
     .               WK(1,5),WK(1,6),WK(1,7),WK(1,8),
     .               WK(1,9),WK(1,10), YS,YP,IERR)
      ENDIF
      Q2 = 0.D0
      DO 5 I = 1,NN
        Q2 = Q2 + W(I)*(Y(I)-YS(I))**2
    5   CONTINUE
      G = 1.D0/SQRT(Q2) - S
      ITER = ITER + 1
      IF (LUN .GE. 0) THEN
        P0 = P - DP
        IF (LUN .GE. 0) WRITE (LUN,120) ITER, P, G, P0, G0
  120   FORMAT (/1X,I2,' -- P = ',D15.8,',  G = ',D15.8/
     .          6X,'P0 = ',D15.8,', G0 = ',D15.8)
      ENDIF
C
C   Test for convergence.
C
      IF ( G .EQ. G0  .OR.  (Q2MIN .LE. Q2  .AND.
     .                       Q2 .LE. Q2MAX) )      RETURN
      IF (DMAX .NE. 0.D0  .OR.  G .GT. 0.D0) GO TO 6
C
C   Increase P until G(P) > 0.
C
      P = 10.D0*P
      DP = P
      GO TO 4
C
C   A bracketing interval [P0,P] has been found.
C
    6 IF (G0*G .LE. 0.D0) THEN
C
C   G0*G < 0.  Update (PNEG,GNEG) to (P0,G0) so that G
C     and GNEG always have opposite signs.
C
        DMAX = DP
        GNEG = G0
      ENDIF
C
C   Compute the change in P by linear interpolation between
C     (P0,G0) and (P,G).
C
    7 DP = -G*DP/(G-G0)
      IF (LUN .GE. 0) WRITE (LUN,130) DP
  130 FORMAT (1X ,5X,'DP = ',D15.8)
      IF (ABS(DP) .GT. ABS(DMAX)) THEN
C
C   G0*G > 0, and the new estimate would be outside of the
C     bracketing interval of length abs(DMAX).  Reset
C     (P0,G0) to (PNEG,GNEG).
C
        DP = DMAX
        G0 = GNEG
        GO TO 7
      ENDIF
C
C   Bottom of loop:  update P, DMAX, and G0.
C
      P = P + DP
      DMAX = DMAX + DP
      G0 = G
      GO TO 4
      END
