      DOUBLE PRECISION FUNCTION ENDSLP (X1,X2,X3,Y1,Y2,Y3,
     .                                  SIGMA)
      DOUBLE PRECISION X1, X2, X3, Y1, Y2, Y3, SIGMA
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
C   Given data values associated with a strictly increasing
C or decreasing sequence of three abscissae X1, X2, and X3,
C this function returns a derivative estimate at X1 based on
C the tension spline H(x) which interpolates the data points
C and has third derivative equal to zero at X1.  Letting S1
C denote the slope defined by the first two points, the est-
C mate is obtained by constraining the derivative of H at X1
C so that it has the same sign as S1 and its magnitude is
C at most 3*abs(S1).  If SIGMA = 0, H(x) is quadratic and
C the derivative estimate is identical to the value computed
C by Subroutine YPC1 at the first point (or the last point
C if the abscissae are decreasing).
C
C On input:
C
C       X1,X2,X3 = Abscissae satisfying either X1 < X2 < X3
C                  or X1 > X2 > X3.
C
C       Y1,Y2,Y3 = Data values associated with the abscis-
C                  sae.  H(X1) = Y1, H(X2) = Y2, and H(X3)
C                  = Y3.
C
C       SIGMA = Tension factor associated with H in inter-
C               val (X1,X2) or (X2,X1).
C
C Input parameters are not altered by this function.
C
C On output:
C
C       ENDSLP = (Constrained) derivative of H at X1, or
C                zero if the abscissae are not strictly
C                monotonic.
C
C Module required by ENDSLP:  SNHCSH
C
C Intrinsic functions called by ENDSLP:  ABS, EXP, MAX, MIN
C
C***********************************************************
C
      DOUBLE PRECISION COSHM1, COSHMS, DUMMY, DX1, DXS, S1,
     .                 SIG1, SIGS, T
C
      DX1 = X2 - X1
      DXS = X3 - X1
      IF (DX1*(DXS-DX1) .LE. 0.D0) GO TO 2
      SIG1 = ABS(SIGMA)
      IF (SIG1 .LT. 1.D-9) THEN
C
C SIGMA = 0:  H is the quadratic interpolant.
C
        T = (DX1/DXS)**2
        GO TO 1
      ENDIF
      SIGS = SIG1*DXS/DX1
      IF (SIGS .LE. .5D0) THEN
C
C 0 < SIG1 < SIGS .LE. .5:  compute approximations to
C   COSHM1 = COSH(SIG1)-1 and COSHMS = COSH(SIGS)-1.
C
        CALL SNHCSH (SIG1, DUMMY,COSHM1,DUMMY)
        CALL SNHCSH (SIGS, DUMMY,COSHMS,DUMMY)
        T = COSHM1/COSHMS
      ELSE
C
C SIGS > .5:  compute T = COSHM1/COSHMS.
C
        T = EXP(SIG1-SIGS)*((1.D0-EXP(-SIG1))/
     .                      (1.D0-EXP(-SIGS)))**2
      ENDIF
C
C The derivative of H at X1 is
C   T = ((Y3-Y1)*COSHM1-(Y2-Y1)*COSHMS)/
C       (DXS*COSHM1-DX1*COSHMS).
C
C ENDSLP = T unless T*S1 < 0 or abs(T) > 3*abs(S1).
C
    1 T = ((Y3-Y1)*T-Y2+Y1)/(DXS*T-DX1)
      S1 = (Y2-Y1)/DX1
      IF (S1 .GE. 0.D0) THEN
        ENDSLP = MIN(MAX(0.D0,T), 3.D0*S1)
      ELSE
        ENDSLP = MAX(MIN(0.D0,T), 3.D0*S1)
      ENDIF
      RETURN
C
C Error in the abscissae.
C
    2 ENDSLP = 0.D0
      RETURN
      END
