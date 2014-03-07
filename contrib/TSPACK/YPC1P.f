      SUBROUTINE YPC1P (N,X,Y, YP,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), YP(N)
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
C   This subroutine employs a three-point quadratic interpo-
C lation method to compute local derivative estimates YP
C associated with a set of N data points (X(I),Y(I)).  It
C is assumed that Y(N) = Y(1), and YP(N) = YP(1) on output.
C Thus, a Hermite interpolant H(x) defined by the data
C points and derivative estimates is periodic with period
C X(N)-X(1).  The derivative-estimation formula is the
C monotonicity-constrained parabolic fit described in the
C reference cited below:  H(x) is monotonic in intervals in
C which the data is monotonic.  The method is invariant
C under a linear scaling of the coordinates but is not
C additive.
C
C On input:
C
C       N = Number of data points.  N .GE. 3.
C
C       X = Array of length N containing a strictly in-
C           creasing sequence of abscissae:  X(I) < X(I+1)
C           for I = 1,...,N-1.
C
C       Y = Array of length N containing data values asso-
C           ciated with the abscissae.  Y(N) is set to Y(1)
C           on output unless IER = 1.
C
C   Input parameters, other than Y(N), are not altered by
C this routine.
C
C On output:
C
C       YP = Array of length N containing estimated deriv-
C            atives at the abscissae unless IER .NE. 0.
C            YP is not altered if IER = 1, and is only par-
C            tially defined if IER > 1.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N < 3.
C             IER = I if X(I) .LE. X(I-1) for some I in the
C                     range 2,...,N.
C
C Reference:  J. M. Hyman, "Accurate Monotonicity-preserving
C               Cubic Interpolation",  LA-8796-MS, Los
C               Alamos National Lab, Feb. 1982.
C
C Modules required by YPC1P:  None
C
C Intrinsic functions called by YPC1P:  ABS, MAX, MIN, SIGN
C
C***********************************************************
C
      INTEGER I, NM1
      DOUBLE PRECISION ASI, ASIM1, DXI, DXIM1, SGN, SI,
     .                 SIM1, T
C
      NM1 = N - 1
      IF (NM1 .LT. 2) GO TO 2
      Y(N) = Y(1)
C
C Initialize for loop on interior points.
C
      I = 1
      DXI = X(2) - X(1)
      IF (DXI .LE. 0.D0) GO TO 3
      SI = (Y(2)-Y(1))/DXI
C
C YP(I) = (DXIM1*SI+DXI*SIM1)/(DXIM1+DXI) subject to the
C   constraint that YP(I) has the sign of either SIM1 or
C   SI, whichever has larger magnitude, and abs(YP(I)) .LE.
C   3*min(abs(SIM1),abs(SI)).
C
      DO 1 I = 2,NM1
        DXIM1 = DXI
        DXI = X(I+1) - X(I)
        IF (DXI .LE. 0.D0) GO TO 3
        SIM1 = SI
        SI = (Y(I+1)-Y(I))/DXI
        T = (DXIM1*SI+DXI*SIM1)/(DXIM1+DXI)
        ASIM1 = ABS(SIM1)
        ASI = ABS(SI)
        SGN = SIGN(1.D0,SI)
        IF (ASIM1 .GT. ASI) SGN = SIGN(1.D0,SIM1)
        IF (SGN .GT. 0.D0) THEN
          YP(I) = MIN(MAX(0.D0,T),3.D0*MIN(ASIM1,ASI))
        ELSE
          YP(I) = MAX(MIN(0.D0,T),-3.D0*MIN(ASIM1,ASI))
        ENDIF
    1   CONTINUE
C
C YP(N) = YP(1), I = 1, and IM1 = N-1.
C
      DXIM1 = DXI
      DXI = X(2) - X(1)
      SIM1 = SI
      SI = (Y(2) - Y(1))/DXI
      T = (DXIM1*SI + DXI*SIM1)/(DXIM1+DXI)
      ASIM1 = ABS(SIM1)
      ASI = ABS(SI)
      SGN = SIGN(1.D0,SI)
      IF (ASIM1 .GT. ASI) SGN = SIGN(1.D0,SIM1)
      IF (SGN .GT. 0.D0) THEN
        YP(1) = MIN(MAX(0.D0,T), 3.D0*MIN(ASIM1,ASI))
      ELSE
        YP(1) = MAX(MIN(0.D0,T), -3.D0*MIN(ASIM1,ASI))
      ENDIF
      YP(N) = YP(1)
C
C No error encountered.
C
      IER = 0
      RETURN
C
C N is outside its valid range.
C
    2 IER = 1
      RETURN
C
C X(I+1) .LE. X(I).
C
    3 IER = I + 1
      RETURN
      END
