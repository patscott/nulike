      SUBROUTINE YPC2P (N,X,Y,SIGMA,WK, YP,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), SIGMA(N), WK(*), YP(N)
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
C   This subroutine solves a linear system for a set of
C first derivatives YP associated with a Hermite interpola-
C tory tension spline H(x).  The derivatives are chosen so
C that H(x) has two continuous derivatives for all x, and H
C satisfies periodic end conditions:  first and second der-
C ivatives of H at X(1) agree with those at X(N), and thus
C the length of a period is X(N) - X(1).  It is assumed that
C Y(N) = Y(1), and Y(N) is not referenced.
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
C           ciated with the abscissae.  H(X(I)) = Y(I) for
C           I = 1,...,N.
C
C       SIGMA = Array of length N-1 containing tension
C               factors.  SIGMA(I) is associated with inter-
C               val (X(I),X(I+1)) for I = 1,...,N-1.  If
C               SIGMA(I) = 0, H is the Hermite cubic interp-
C               olant of the data values and computed deriv-
C               atives at X(I) and X(I+1), and if all
C               tension factors are zero, H is the C-2 cubic
C               spline interpolant which satisfies the end
C               conditions.
C
C The above parameters are not altered by this routine.
C
C       WK = Array of length at least 2N-2 to be used as
C            temporary work space.
C
C       YP = Array of length .GE. N.
C
C On output:
C
C       YP = Array containing derivatives of H at the
C            abscissae.  YP is not defined if IER .NE. 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N is outside its valid range.
C             IER = I if X(I) .LE. X(I-1) for some I in the
C                     range 2,...,N.
C
C Modules required by YPC2P:  SNHCSH, YPCOEF
C
C Intrinsic function called by YPC2P:  ABS
C
C***********************************************************
C
      INTEGER I, NM1, NM2, NM3, NN, NP1, NPI
      DOUBLE PRECISION D, D1, D2, DIN, DNM1, DX, R1, R2,
     .                 RNM1, S, SD1, SD2, SDNM1, SIG, YPNM1
C
      NN = N
      IF (NN .LT. 3) GO TO 4
      NM1 = NN - 1
      NM2 = NN - 2
      NM3 = NN - 3
      NP1 = NN + 1
C
C The system is order N-1, symmetric, positive-definite, and
C   tridiagonal except for nonzero elements in the upper
C   right and lower left corners.  The forward elimination
C   step zeros the subdiagonal and divides each row by its
C   diagonal entry for the first N-2 rows.  The superdiago-
C   nal is stored in WK(I), the negative of the last column
C   (fill-in) in WK(N+I), and the right hand side in YP(I)
C   for I = 1,...,N-2.
C
      I = NM1
      DX = X(NN) - X(NM1)
      IF (DX .LE. 0.D0) GO TO 5
      S = (Y(1)-Y(NM1))/DX
      SIG = ABS(SIGMA(NM1))
      CALL YPCOEF (SIG,DX, DNM1,SDNM1)
      RNM1 = (SDNM1+DNM1)*S
      I = 1
      DX = X(2) - X(1)
      IF (DX .LE. 0.D0) GO TO 5
      S = (Y(2)-Y(1))/DX
      SIG = ABS(SIGMA(1))
      CALL YPCOEF (SIG,DX, D1,SD1)
      R1 = (SD1+D1)*S
      D = DNM1 + D1
      WK(1) = SD1/D
      WK(NP1) = -SDNM1/D
      YP(1) = (RNM1+R1)/D
      DO 1 I = 2,NM2
        DX = X(I+1) - X(I)
        IF (DX .LE. 0.D0) GO TO 5
        S = (Y(I+1)-Y(I))/DX
        SIG = ABS(SIGMA(I))
        CALL YPCOEF (SIG,DX, D2,SD2)
        R2 = (SD2+D2)*S
        D = D1 + D2 - SD1*WK(I-1)
        DIN = 1.D0/D
        WK(I) = SD2*DIN
        NPI = NN + I
        WK(NPI) = -SD1*WK(NPI-1)*DIN
        YP(I) = (R1 + R2 - SD1*YP(I-1))*DIN
        SD1 = SD2
        D1 = D2
        R1 = R2
    1   CONTINUE
C
C The backward elimination step zeros the superdiagonal
C   (first N-3 elements).  WK(I) and YP(I) are overwritten
C   with the negative of the last column and the new right
C   hand side, respectively, for I = N-2, N-3, ..., 1.
C
      NPI = NN + NM2
      WK(NM2) = WK(NPI) - WK(NM2)
      DO 2 I = NM3,1,-1
        YP(I) = YP(I) - WK(I)*YP(I+1)
        NPI = NN + I
        WK(I) = WK(NPI) - WK(I)*WK(I+1)
    2   CONTINUE
C
C Solve the last equation for YP(N-1).
C
      YPNM1 = (R1 + RNM1 - SDNM1*YP(1) - SD1*YP(NM2))/
     .        (D1 + DNM1 + SDNM1*WK(1) + SD1*WK(NM2))
C
C Back substitute for the remainder of the solution
C   components.
C
      YP(NM1) = YPNM1
      DO 3 I = 1,NM2
        YP(I) = YP(I) + WK(I)*YPNM1
    3   CONTINUE
C
C YP(N) = YP(1).
C
      YP(N) = YP(1)
      IER = 0
      RETURN
C
C N is outside its valid range.
C
    4 IER = 1
      RETURN
C
C Abscissae out of order or duplicate points.
C
    5 IER = I + 1
      RETURN
      END
