      DOUBLE PRECISION FUNCTION TSINTL (A,B,N,X,Y,YP,
     .                                  SIGMA, IER)
      INTEGER N, IER
      DOUBLE PRECISION A, B, X(N), Y(N), YP(N), SIGMA(N)
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
C   This function computes the integral from A to B of a
C Hermite interpolatory tension spline H.
C
C On input:
C
C       A,B = Lower and upper limits of integration, re-
C             spectively.  Note that -TSINTL(B,A,...) =
C             TSINTL(A,B,...).
C
C       N = Number of data points.  N .GE. 2.
C
C       X = Array of length N containing the abscissae.
C           These must be in strictly increasing order:
C           X(I) < X(I+1) for I = 1,...,N-1.
C
C       Y = Array of length N containing data values.
C           H(X(I)) = Y(I) for I = 1,...,N.
C
C       YP = Array of length N containing first deriva-
C            tives.  HP(X(I)) = YP(I) for I = 1,...,N, where
C            HP denotes the derivative of H.
C
C       SIGMA = Array of length N-1 containing tension fac-
C               tors whose absolute values determine the
C               balance between cubic and linear in each
C               interval.  SIGMA(I) is associated with int-
C               erval (I,I+1) for I = 1,...,N-1.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0  if no errors were encountered and
C                      X(1) .LE. T .LE. X(N) for T = A and
C                      T = B, or A = B.
C             IER = 1  if no errors were encountered but
C                      extrapolation was necessary:  A or B
C                      not in the interval (X(1),X(N)).
C             IER = -1 IF N < 2.
C             IER = -2 if the abscissae are not in strictly
C                      increasing order.  Only those in or
C                      adjacent to the interval of integra-
C                      tion are tested.
C
C       TSINTL = Integral of H from A to B, or zero if
C                IER < 0.
C
C Modules required by TSINTL:  INTRVL, SNHCSH
C
C Intrinsic functions called by TSINTL:  ABS, EXP, MAX, MIN
C
C***********************************************************
C
      INTEGER I, IL, ILP1, IMAX, IMIN, IP1, IU, IUP1
      DOUBLE PRECISION B1, B2, CM, CM1, CM2, CMM, CMM1,
     .                 CMM2, D1, D2, DX, E, E1, E2, EMS, S,
     .                 S1, S2, SB1, SB2, SBIG, SIG, SM, SM1,
     .                 SM2, SUM, T, TM, TP, U, XL, XU, Y1,
     .                 Y2
      INTEGER INTRVL
C
      DATA SBIG/85.D0/
      IF (N .LT. 2) GO TO 7
C
C Accumulate the integral from XL to XU in SUM.
C
      XL = MIN(A,B)
      XU = MAX(A,B)
      SUM = 0.D0
      IER = 0
      IF (XL .EQ. XU) GO TO 6
C
C Find left-end indexes of intervals containing XL and XU.
C   If XL < X(1) or XU > X(N), extrapolation is performed
C   using the leftmost or rightmost interval.
C
      IL = INTRVL (XL,N,X)
      IU = INTRVL (XU,N,X)
      IF (XL .LT. X(1)  .OR.  XU .GT. X(N)) IER = 1
      ILP1 = IL + 1
      IMIN = IL
      IF (XL .EQ. X(IL)) GO TO 2
C
C Compute the integral from XL to X(IL+1).
C
      DX = X(ILP1) - X(IL)
      IF (DX .LE. 0.D0) GO TO 8
      U = X(ILP1) - XL
      IF (U .EQ. 0.D0) GO TO 1
      B1 = U/DX
      Y2 = Y(ILP1)
      S = (Y2-Y(IL))/DX
      S2 = YP(ILP1)
      D1 = S - YP(IL)
      D2 = S2 - S
      SIG = ABS(SIGMA(IL))
      IF (SIG .LT. 1.D-9) THEN
C
C   SIG = 0.
C
        SUM = SUM + U*(Y2 - U*(6.D0*S2 - B1*(4.D0*D2 +
     .                 (3.D0*B1-4.D0)*(D1-D2)))/12.D0)
      ELSEIF (SIG .LE. .5D0) THEN
C
C   0 .LT. SIG .LE. .5.
C
        SB1 = SIG*B1
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB1, SM1,CM1,CMM1)
        E = SIG*SM - CMM - CMM
        SUM = SUM + U*(Y2 - S2*U/2.D0) + ((CM*CMM1-SM*SM1)*
     .             (D1+D2) + SIG*(CM*SM1-(SM+SIG)*CMM1)*D2)/
     .             ((SIG/DX)**2*E)
      ELSE
C
C   SIG > .5.
C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
          SUM = SUM + U*(Y2 - S*U/2.D0)
        ELSE
          E1 = EXP(-SB1)
          E2 = EXP(-SB2)
          EMS = E1*E2
          TM = 1.D0 - EMS
          TP = 1.D0 + EMS
          T = SB1*SB1/2.D0 + 1.D0
          E = TM*(SIG*TP - TM - TM)
          SUM = SUM +U*(Y2 - S2*U/2.D0)+(SIG*TM*(TP*T-E1-E2-
     .               TM*SB1)*D2 - (TM*(TM*T-E1+E2-TP*SB1) +
     .               SIG*(E1*EMS-E2+2.D0*SB1*EMS))*(D1+D2))/
     .               ((SIG/DX)**2*E)
        ENDIF
      ENDIF
C
C Test for XL and XU in the same interval.
C
    1 IMIN = ILP1
      IF (IL .EQ. IU) GO TO 5
C
C Add in the integral from X(IMIN) to X(J) for J =
C   Max(IL+1,IU).
C
    2 IMAX = MAX(IL,IU-1)
      DO 3 I = IMIN,IMAX
        IP1 = I + 1
        DX = X(IP1) - X(I)
        IF (DX .LE. 0.D0) GO TO 8
        SIG = ABS(SIGMA(I))
        IF (SIG .LT. 1.D-9) THEN
C
C   SIG = 0.
C
          SUM = SUM + DX*((Y(I)+Y(IP1))/2.D0 -
     .                    DX*(YP(IP1)-YP(I))/12.D0)
        ELSEIF (SIG .LE. .5D0) THEN
C
C   0 .LT. SIG .LE. .5.
C
          CALL SNHCSH (SIG, SM,CM,CMM)
          E = SIG*SM - CMM - CMM
          SUM = SUM + DX*(Y(I)+Y(IP1) - DX*E*(YP(IP1)-YP(I))
     .                /(SIG*SIG*CM))/2.D0
        ELSE
C
C   SIG > .5.
C
          EMS = EXP(-SIG)
          SUM = SUM + DX*(Y(I)+Y(IP1) - DX*(SIG*(1.D0+EMS)/
     .                (1.D0-EMS)-2.D0)*(YP(IP1)-YP(I))/
     .                (SIG*SIG))/2.D0
        ENDIF
    3   CONTINUE
C
C Add in the integral from X(IU) to XU if IU > IL.
C
      IF (IL .EQ. IU) GO TO 4
      IUP1 = IU + 1
      DX = X(IUP1) - X(IU)
      IF (DX .LE. 0.D0) GO TO 8
      U = XU - X(IU)
      IF (U .EQ. 0.D0) GO TO 6
      B2 = U/DX
      Y1 = Y(IU)
      S = (Y(IUP1)-Y1)/DX
      S1 = YP(IU)
      D1 = S - S1
      D2 = YP(IUP1) - S
      SIG = ABS(SIGMA(IU))
      IF (SIG .LT. 1.D-9) THEN
C
C   SIG = 0.
C
        SUM = SUM + U*(Y1 + U*(6.D0*S1 + B2*(4.D0*D1 +
     .                (4.D0-3.D0*B2)*(D1-D2)))/12.D0)
      ELSEIF (SIG .LE. .5D0) THEN
C
C   0 .LT. SIG .LE. .5.
C
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,CMM2)
        E = SIG*SM - CMM - CMM
        SUM = SUM + U*(Y1 + S1*U/2.D0) + ((CM*CMM2-SM*SM2)*
     .              (D1+D2) + SIG*(CM*SM2-(SM+SIG)*CMM2)*D1)
     .              /((SIG/DX)**2*E)
      ELSE
C
C   SIG > .5.
C
        SB2 = SIG*B2
        SB1 = SIG - SB2
        IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
          SUM = SUM + U*(Y1 + S*U/2.D0)
        ELSE
          E1 = EXP(-SB1)
          E2 = EXP(-SB2)
          EMS = E1*E2
          TM = 1.D0 - EMS
          TP = 1.D0 + EMS
          T = SB2*SB2/2.D0 + 1.D0
          E = TM*(SIG*TP - TM - TM)
          SUM = SUM +U*(Y1 + S1*U/2.D0)+(SIG*TM*(TP*T-E1-E2-
     .               TM*SB2)*D1 - (TM*(TM*T-E2+E1-TP*SB2) +
     .               SIG*(E2*EMS-E1+2.D0*SB2*EMS))*(D1+D2))/
     .               ((SIG/DX)**2*E)
        ENDIF
      ENDIF
      GO TO 6
C
C IL = IU and SUM contains the integral from XL to X(IL+1).
C   Subtract off the integral from XU to X(IL+1).  DX and
C   SIG were computed above.
C
    4 Y2 = Y(ILP1)
      S = (Y2-Y(IL))/DX
      S2 = YP(ILP1)
      D1 = S - YP(IL)
      D2 = S2 - S
C
    5 U = X(ILP1) - XU
      IF (U .EQ. 0.D0) GO TO 6
      B1 = U/DX
      IF (SIG .LT. 1.D-9) THEN
C
C   SIG = 0.
C
        SUM = SUM - U*(Y2 - U*(6.D0*S2 - B1*(4.D0*D2 +
     .              (3.D0*B1-4.D0)*(D1-D2)))/12.D0)
      ELSEIF (SIG .LE. .5D0) THEN
C
C   0 .LT. SIG .LE. .5.
C
        SB1 = SIG*B1
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB1, SM1,CM1,CMM1)
        E = SIG*SM - CMM - CMM
        SUM = SUM - U*(Y2 - S2*U/2.D0) - ((CM*CMM1-SM*SM1)*
     .              (D1+D2) + SIG*(CM*SM1-(SM+SIG)*CMM1)*D2)
     .              /((SIG/DX)**2*E)
      ELSE
C
C   SIG > .5.
C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
          SUM = SUM - U*(Y2 - S*U/2.D0)
        ELSE
          E1 = EXP(-SB1)
          E2 = EXP(-SB2)
          EMS = E1*E2
          TM = 1.D0 - EMS
          TP = 1.D0 + EMS
          T = SB1*SB1/2.D0 + 1.D0
          E = TM*(SIG*TP - TM - TM)
          SUM = SUM -U*(Y2 - S2*U/2.D0)-(SIG*TM*(TP*T-E1-E2-
     .               TM*SB1)*D2 - (TM*(TM*T-E1+E2-TP*SB1) +
     .               SIG*(E1*EMS-E2+2.D0*SB1*EMS))*(D1+D2))/
     .               ((SIG/DX)**2*E)
        ENDIF
      ENDIF
C
C No errors were encountered.  Adjust the sign of SUM.
C
    6 IF (XL .EQ. B) SUM = -SUM
      TSINTL = SUM
      RETURN
C
C N < 2.
C
    7 IER = -1
      TSINTL = 0.D0
      RETURN
C
C Abscissae not strictly increasing.
C
    8 IER = -2
      TSINTL = 0.D0
      RETURN
      END
