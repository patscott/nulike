      SUBROUTINE YPC1 (N,X,Y, YP,IER)
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
C associated with a set of data points.  The interpolation
C formula is the monotonicity-constrained parabolic method
C described in the reference cited below.  A Hermite int-
C erpolant of the data values and derivative estimates pre-
C serves monotonicity of the data.  Linear interpolation is
C used if N = 2.  The method is invariant under a linear
C scaling of the coordinates but is not additive.
C
C On input:
C
C       N = Number of data points.  N .GE. 2.
C
C       X = Array of length N containing a strictly in-
C           creasing sequence of abscissae:  X(I) < X(I+1)
C           for I = 1,...,N-1.
C
C       Y = Array of length N containing data values asso-
C           ciated with the abscissae.
C
C Input parameters are not altered by this routine.
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
C             IER = 1 if N < 2.
C             IER = I if X(I) .LE. X(I-1) for some I in the
C                     range 2,...,N.
C
C Reference:  J. M. Hyman, "Accurate Monotonicity-preserving
C               Cubic Interpolation",  LA-8796-MS, Los
C               Alamos National Lab, Feb. 1982.
C
C Modules required by YPC1:  None
C
C Intrinsic functions called by YPC1:  ABS, MAX, MIN, SIGN
C
C***********************************************************
C
      INTEGER I, NM1
      DOUBLE PRECISION ASI, ASIM1, DX2, DXI, DXIM1, S2, SGN,
     .                 SI, SIM1, T
C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 2
      I = 1
      DXI = X(2) - X(1)
      IF (DXI .LE. 0.D0) GO TO 3
      SI = (Y(2)-Y(1))/DXI
      IF (NM1 .EQ. 1) THEN
C
C Use linear interpolation for N = 2.
C
        YP(1) = SI
        YP(2) = SI
        IER = 0
        RETURN
      ENDIF
C
C N .GE. 3.  YP(1) = S1 + DX1*(S1-S2)/(DX1+DX2) unless this
C   results in YP(1)*S1 .LE. 0 or abs(YP(1)) > 3*abs(S1).
C
      I = 2
      DX2 = X(3) - X(2)
      IF (DX2 .LE. 0.D0) GO TO 3
      S2 = (Y(3)-Y(2))/DX2
      T = SI + DXI*(SI-S2)/(DXI+DX2)
      IF (SI .GE. 0.D0) THEN
        YP(1) = MIN(MAX(0.D0,T), 3.D0*SI)
      ELSE
        YP(1) = MAX(MIN(0.D0,T), 3.D0*SI)
      ENDIF
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
C YP(N) = SNM1 + DXNM1*(SNM1-SNM2)/(DXNM2+DXNM1) subject to
C   the constraint that YP(N) has the sign of SNM1 and
C   abs(YP(N)) .LE. 3*abs(SNM1).  Note that DXI = DXNM1 and
C   SI = SNM1.
C
      T = SI + DXI*(SI-SIM1)/(DXIM1+DXI)
      IF (SI .GE. 0.D0) THEN
        YP(N) = MIN(MAX(0.D0,T), 3.D0*SI)
      ELSE
        YP(N) = MAX(MIN(0.D0,T), 3.D0*SI)
      ENDIF
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
