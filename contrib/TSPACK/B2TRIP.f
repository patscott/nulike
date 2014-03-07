      SUBROUTINE B2TRIP (N,X,Y,W,P,D,SD,T11,T12,T21,T22,U11,
     .                   U12,U21,U22, YS,YP,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), W(N), P, D(N), SD(N),
     .                 T11(N), T12(N), T21(N), T22(N),
     .                 U11(N), U12(N), U21(N), U22(N),
     .                 YS(N), YP(N)
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
C   This subroutine solves the order 2(N-1) symmetric posi-
C tive-definite linear system associated with minimizing the
C quadratic functional Q(YS,YP) (described in Subroutine
C SMCRV) with periodic end conditions.  The matrix is block
C tridiagonal except for nonzero blocks in the upper right
C and lower left corners.
C
C On input:
C
C       N = Number of data points.  N .GE. 3.
C
C       X = Array of length N containing a strictly in-
C           creasing sequence of abscissae.
C
C       Y,W = Arrays of length N-1 containing data values
C             and positive weights, respectively, associated
C             with the first N-1 abscissae.
C
C       P = Positive smoothing parameter defining Q.
C
C       D,SD = Arrays of length N-1 containing positive ma-
C              trix elements.  Letting DX and SIG denote the
C              width and tension factor associated with the
C              interval (X(I),X(I+1)), D(I) = SIG*(SIG*
C              COSHM(SIG) - SINHM(SIG))/(DX*E) and SD(I) =
C              SIG*SINHM(SIG)/(DX*E) where E = SIG*SINH(SIG)
C              - 2*COSHM(SIG).
C
C The above parameters are not altered by this routine.
C
C       T11,T12,T21,T22,U11,U12,U21,U22 = Arrays of length
C                                         N-2 used as temp-
C                                         orary work space.
C
C On output:
C
C       YS,YP = Arrays of length N containing solution com-
C               ponents:  function and derivative values,
C               respectively, at the abscissae.  YS(N) =
C               YS(1) and YP(N) = YP(1).
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N or P is outside its valid range
C                     on input.
C             Note that no test is made for a nonpositive
C             value of X(I+1)-X(I), W(I), D(I), or SD(I).
C
C Modules required by B2TRIP:  None
C
C***********************************************************
C
      INTEGER I, IM1, IP1, NM1, NM2, NM3, NN
      DOUBLE PRECISION D11I, D12I, D22I, DEN, DI, DIM1,
     .                 DNM1, DX, PP, R1, R2, S11I, S11IM1,
     .                 S11NM1, S12I, S12IM1, S12NM1, S22I,
     .                 S22IM1, S22NM1, SU11, SU12, SU21,
     .                 SU22, YPNM1, YSNM1
C
      NN = N
      NM1 = NN - 1
      NM2 = NN - 2
      NM3 = NN - 3
      PP = P
      IER = 1
      IF (NN .LT. 3  .OR.  PP .LE. 0.D0) RETURN
C
C The forward elimination step consists of scaling a row by
C   the inverse of its diagonal block and eliminating the
C   subdiagonal block for the first N-2 rows.  The super-
C   diagonal is stored in T, the negative of the last column
C   in U, and the right hand side in YS,YP.  For J = 11, 12,
C   and 22, SJI and SJIM1 denote the elements in position J
C   of the superdiagonal block in rows I and I-1, respect-
C   ively.  Similarly, DJI denotes an element in the diago-
C   nal block of row I.
C
C I = 1:
C
      DX = X(NN) - X(NM1)
      DNM1 = D(NM1)
      S22NM1 = SD(NM1)
      S12NM1 = -(DNM1 + S22NM1)/DX
      S11NM1 = 2.D0*S12NM1/DX
      DX = X(2) - X(1)
      DI = D(1)
      S22I = SD(1)
      S12I = (DI + S22I)/DX
      S11I = -2.D0*S12I/DX
      R1 = PP*W(1)
      D11I = R1 - S11NM1 - S11I
      D12I = S12I + S12NM1
      D22I = DNM1 + DI
      DEN = D11I*D22I - D12I*D12I
      T11(1) = (D22I*S11I + D12I*S12I)/DEN
      T12(1) = (D22I*S12I - D12I*S22I)/DEN
      T21(1) = -(D12I*S11I + D11I*S12I)/DEN
      T22(1) = (D11I*S22I - D12I*S12I)/DEN
      U11(1) = -(D22I*S11NM1 + D12I*S12NM1)/DEN
      U12(1) = (D12I*S22NM1 - D22I*S12NM1)/DEN
      U21(1) = (D12I*S11NM1 + D11I*S12NM1)/DEN
      U22(1) = (D12I*S12NM1 - D11I*S22NM1)/DEN
      R1 = R1*Y(1)/DEN
      YS(1) = D22I*R1
      YP(1) = -D12I*R1
      IF (NN .EQ. 3) GO TO 2
C
C I = 2,...,N-2:
C
      DO 1 I = 2,NM2
        IM1 = I - 1
        DIM1 = DI
        S22IM1 = S22I
        S12IM1 = S12I
        S11IM1 = S11I
        DX = X(I+1) - X(I)
        DI = D(I)
        S22I = SD(I)
        S12I = (DI + S22I)/DX
        S11I = -2.D0*S12I/DX
        R1 = PP*W(I)
        D11I = R1 - S11IM1 - S11I - (S11IM1*T11(IM1) -
     .         S12IM1*T21(IM1))
        D12I = S12I - S12IM1 - (S11IM1*T12(IM1) - S12IM1*
     .         T22(IM1))
        D22I = DIM1 + DI - (S12IM1*T12(IM1)+S22IM1*T22(IM1))
        DEN = D11I*D22I - D12I*D12I
        T11(I) = (D22I*S11I + D12I*S12I)/DEN
        T12(I) = (D22I*S12I - D12I*S22I)/DEN
        T21(I) = -(D12I*S11I + D11I*S12I)/DEN
        T22(I) = (D11I*S22I - D12I*S12I)/DEN
        SU11 = S11IM1*U11(IM1) - S12IM1*U21(IM1)
        SU12 = S11IM1*U12(IM1) - S12IM1*U22(IM1)
        SU21 = S12IM1*U11(IM1) + S22IM1*U21(IM1)
        SU22 = S12IM1*U12(IM1) + S22IM1*U22(IM1)
        U11(I) = (D12I*SU21 - D22I*SU11)/DEN
        U12(I) = (D12I*SU22 - D22I*SU12)/DEN
        U21(I) = (D12I*SU11 - D11I*SU21)/DEN
        U22(I) = (D12I*SU12 - D11I*SU22)/DEN
        R1 = R1*Y(I) - S11IM1*YS(IM1) + S12IM1*YP(IM1)
        R2 = -S12IM1*YS(IM1) - S22IM1*YP(IM1)
        YS(I) = (D22I*R1 - D12I*R2)/DEN
        YP(I) = (D11I*R2 - D12I*R1)/DEN
    1   CONTINUE
C
C The backward elimination step zeros the first N-3 blocks
C   of the superdiagonal.  For I = N-2,N-3,...,1, T(I) and
C   (YS(I),YP(I)) are overwritten with the negative of the
C   last column and the new right hand side, respectively.
C
    2 T11(NM2) = U11(NM2) - T11(NM2)
      T12(NM2) = U12(NM2) - T12(NM2)
      T21(NM2) = U21(NM2) - T21(NM2)
      T22(NM2) = U22(NM2) - T22(NM2)
      DO 3 I = NM3,1,-1
        IP1 = I + 1
        YS(I) = YS(I) - T11(I)*YS(IP1) - T12(I)*YP(IP1)
        YP(I) = YP(I) - T21(I)*YS(IP1) - T22(I)*YP(IP1)
        T11(I) = U11(I) - T11(I)*T11(IP1) - T12(I)*T21(IP1)
        T12(I) = U12(I) - T11(I)*T12(IP1) - T12(I)*T22(IP1)
        T21(I) = U21(I) - T21(I)*T11(IP1) - T22(I)*T21(IP1)
        T22(I) = U22(I) - T21(I)*T12(IP1) - T22(I)*T22(IP1)
    3   CONTINUE
C
C Solve the last equation for YS(N-1),YP(N-1).  SJI = SJNM2
C   and DJI = DJNM1.
C
      R1 = PP*W(NM1)
      D11I = R1 - S11I - S11NM1 + S11NM1*T11(1) -
     .       S12NM1*T21(1) + S11I*T11(NM2) - S12I*T21(NM2)
      D12I = -S12NM1 - S12I + S11NM1*T12(1) - S12NM1*T22(1)
     .       + S11I*T12(NM2) - S12I*T22(NM2)
      D22I = DI + DNM1 + S12NM1*T12(1) + S22NM1*T22(1) +
     .       S12I*T12(NM2) + S22I*T22(NM2)
      DEN = D11I*D22I - D12I*D12I
      R1 = R1*Y(NM1) - S11NM1*YS(1) + S12NM1*YP(1) -
     .     S11I*YS(NM2) + S12I*YP(NM2)
      R2 = -S12NM1*YS(1) - S22NM1*YP(1) - S12I*YS(NM2) -
     .     S22I*YP(NM2)
      YSNM1 = (D22I*R1 - D12I*R2)/DEN
      YPNM1 = (D11I*R2 - D12I*R1)/DEN
      YS(NM1) = YSNM1
      YP(NM1) = YPNM1
C
C Back substitute for the remainder of the solution
C   components.
C
      DO 4 I = 1,NM2
        YS(I) = YS(I) + T11(I)*YSNM1 + T12(I)*YPNM1
        YP(I) = YP(I) + T21(I)*YSNM1 + T22(I)*YPNM1
    4   CONTINUE
C
C YS(N) = YS(1) and YP(N) = YP(1).
C
      YS(NN) = YS(1)
      YP(NN) = YP(1)
      IER = 0
      RETURN
      END
