      SUBROUTINE YPC2 (N,X,Y,SIGMA,ISL1,ISLN,BV1,BVN,
     .                 WK, YP,IER)
      INTEGER N, ISL1, ISLN, IER
      DOUBLE PRECISION X(N), Y(N), SIGMA(N), BV1, BVN,
     .                 WK(N), YP(N)
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
C that H(x) has two continuous derivatives for all x and H
C satisfies user-specified end conditions.
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
C       ISL1 = Option indicator for the condition at X(1):
C              ISL1 = 0 if YP(1) is to be estimated inter-
C                       nally by a constrained parabolic
C                       fit to the first three points.
C                       This is identical to the method used
C                       by Subroutine YPC1.  BV1 is not used
C                       in this case.
C              ISL1 = 1 if the first derivative of H at X(1)
C                       is specified by BV1.
C              ISL1 = 2 if the second derivative of H at
C                       X(1) is specified by BV1.
C              ISL1 = 3 if YP(1) is to be estimated inter-
C                       nally from the derivative of the
C                       tension spline (using SIGMA(1))
C                       which interpolates the first three
C                       data points and has third derivative
C                       equal to zero at X(1).  Refer to
C                       ENDSLP.  BV1 is not used in this
C                       case.
C
C       ISLN = Option indicator for the condition at X(N):
C              ISLN = 0 if YP(N) is to be estimated inter-
C                       nally by a constrained parabolic
C                       fit to the last three data points.
C                       This is identical to the method used
C                       by Subroutine YPC1.  BVN is not used
C                       in this case.
C              ISLN = 1 if the first derivative of H at X(N)
C                       is specified by BVN.
C              ISLN = 2 if the second derivative of H at
C                       X(N) is specified by BVN.
C              ISLN = 3 if YP(N) is to be estimated inter-
C                       nally from the derivative of the
C                       tension spline (using SIGMA(N-1))
C                       which interpolates the last three
C                       data points and has third derivative
C                       equal to zero at X(N).  Refer to
C                       ENDSLP.  BVN is not used in this
C                       case.
C
C       BV1,BVN = Boundary values or dummy parameters as
C                 defined by ISL1 and ISLN.
C
C The above parameters are not altered by this routine.
C
C       WK = Array of length at least N-1 to be used as
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
C             IER = 1 if N, ISL1, or ISLN is outside its
C                     valid range.
C             IER = I if X(I) .LE. X(I-1) for some I in the
C                     range 2,...,N.
C
C Modules required by YPC2:  ENDSLP, SNHCSH, YPCOEF
C
C Intrinsic function called by YPC2:  ABS
C
C***********************************************************
C
      INTEGER I, NM1, NN
      DOUBLE PRECISION D, D1, D2, DX, R1, R2, S, SD1, SD2,
     .                 SIG, YP1, YPN
      DOUBLE PRECISION ENDSLP
C
      NN = N
      IF (NN .LT. 2  .OR.  ISL1 .LT. 0  .OR.  ISL1 .GT. 3
     .    .OR.  ISLN .LT. 0  .OR.  ISLN .GT. 3) GO TO 3
      NM1 = NN - 1
C
C Set YP1 and YPN to the endpoint values.
C
      IF (ISL1 .EQ. 0) THEN
        IF (NN .GT. 2) YP1 = ENDSLP (X(1),X(2),X(3),Y(1),
     .                               Y(2),Y(3),0.D0)
      ELSEIF (ISL1 .NE. 3) THEN
        YP1 = BV1
      ELSE
        IF (NN .GT. 2) YP1 = ENDSLP (X(1),X(2),X(3),Y(1),
     .                               Y(2),Y(3),SIGMA(1))
      ENDIF
      IF (ISLN .EQ. 0) THEN
        IF (NN .GT. 2) YPN = ENDSLP (X(NN),X(NM1),X(NN-2),
     .                            Y(NN),Y(NM1),Y(NN-2),0.D0)
      ELSEIF (ISLN .NE. 3) THEN
        YPN = BVN
      ELSE
        IF (NN .GT. 2) YPN = ENDSLP (X(NN),X(NM1),X(NN-2),
     .                      Y(NN),Y(NM1),Y(NN-2),SIGMA(NM1))
      ENDIF
C
C Solve the symmetric positive-definite tridiagonal linear
C   system.  The forward elimination step consists of div-
C   iding each row by its diagonal entry, then introducing a
C   zero below the diagonal.  This requires saving only the
C   superdiagonal (in WK) and the right hand side (in YP).
C
      I = 1
      DX = X(2) - X(1)
      IF (DX .LE. 0.D0) GO TO 4
      S = (Y(2)-Y(1))/DX
      IF (NN .EQ. 2) THEN
        IF (ISL1 .EQ. 0  .OR.  ISL1 .EQ. 3) YP1 = S
        IF (ISLN .EQ. 0  .OR.  ISLN .EQ. 3) YPN = S
      ENDIF
C
C Begin forward elimination.
C
      SIG = ABS(SIGMA(1))
      CALL YPCOEF (SIG,DX, D1,SD1)
      R1 = (SD1+D1)*S
      WK(1) = 0.D0
      YP(1) = YP1
      IF (ISL1 .EQ. 2) THEN
        WK(1) = SD1/D1
        YP(1) = (R1-YP1)/D1
      ENDIF
      DO 1 I = 2,NM1
        DX = X(I+1) - X(I)
        IF (DX .LE. 0.D0) GO TO 4
        S = (Y(I+1)-Y(I))/DX
        SIG = ABS(SIGMA(I))
        CALL YPCOEF (SIG,DX, D2,SD2)
        R2 = (SD2+D2)*S
        D = D1 + D2 - SD1*WK(I-1)
        WK(I) = SD2/D
        YP(I) = (R1 + R2 - SD1*YP(I-1))/D
        D1 = D2
        SD1 = SD2
        R1 = R2
    1   CONTINUE
      D = D1 - SD1*WK(NM1)
      YP(NN) = YPN
      IF (ISLN .EQ. 2) YP(NN) = (R1 + YPN - SD1*YP(NM1))/D
C
C Back substitution:
C
      DO 2 I = NM1,1,-1
        YP(I) = YP(I) - WK(I)*YP(I+1)
    2   CONTINUE
      IER = 0
      RETURN
C
C Invalid integer input parameter.
C
    3 IER = 1
      RETURN
C
C Abscissae out of order or duplicate points.
C
    4 IER = I + 1
      RETURN
      END
