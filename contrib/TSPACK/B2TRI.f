      SUBROUTINE B2TRI (N,X,Y,W,P,D,SD,T11,T12,T21,T22, YS,
     .                  YP,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), W(N), P, D(N), SD(N),
     .                 T11(N), T12(N), T21(N), T22(N),
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
C   This subroutine solves the order 2N symmetric positive-
C definite block tridiagonal linear system associated with
C minimizing the quadratic functional Q(YS,YP) described in
C Subroutine SMCRV.
C
C On input:
C
C       N = Number of data points.  N .GE. 2.
C
C       X,Y,W = Arrays of length N containing abscissae,
C               data values, and positive weights, respect-
C               ively.  The abscissae must be strictly in-
C               creasing.
C
C       P = Positive smoothing parameter defining Q.
C
C       D,SD = Arrays of length N-1 containing positive ma-
C              trix entries.  Letting DX and SIG denote the
C              width and tension factor associated with the
C              interval (X(I),X(I+1)), D(I) = SIG*(SIG*
C              COSHM(SIG) - SINHM(SIG))/(DX*E) and SD(I) =
C              SIG*SINHM(SIG)/(DX*E) where E = SIG*SINH(SIG)
C              - 2*COSHM(SIG).
C
C The above parameters are not altered by this routine.
C
C       T11,T12,T21,T22 = Arrays of length N-1 used as
C                         temporary work space.
C
C On output:
C
C       YS,YP = Arrays of length N containing solution com-
C               ponents:  function and derivative values,
C               respectively, at the abscissae.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N or P is outside its valid range
C                     on input.
C             Note that no test is made for a nonpositive
C             value of X(I+1)-X(I), W(I), D(I), or SD(I).
C
C Modules required by B2TRI:  None
C
C***********************************************************
C
      INTEGER I, IM1, NM1, NN
      DOUBLE PRECISION D11I, D12I, D22I, DEN, DI, DIM1, DX,
     .                 PP, R1, R2, S11I, S11IM1, S12I,
     .                 S12IM1, S22I, S22IM1
C
      NN = N
      NM1 = NN - 1
      PP = P
      IER = 1
      IF (NN .LT. 2  .OR.  PP .LE. 0.D0) RETURN
C
C The forward elimination step consists of scaling a row by
C   the inverse of its diagonal block and eliminating the
C   subdiagonal block.  The superdiagonal is stored in T and
C   the right hand side in YS,YP.  For J = 11, 12, and 22,
C   SJI and SJIM1 denote the elements in position J of the
C   superdiagonal block in rows I and I-1, respectively.
C   Similarly, DJI denotes an element in the diagonal block
C   of row I.
C
C Initialize for I = 2.
C
      DX = X(2) - X(1)
      DIM1 = D(1)
      S22IM1 = SD(1)
      S12IM1 = (DIM1 + S22IM1)/DX
      S11IM1 = -2.D0*S12IM1/DX
      R1 = PP*W(1)
      D11I = R1 - S11IM1
      D12I = S12IM1
      D22I = DIM1
      DEN = D11I*D22I - D12I*D12I
      T11(1) = (D22I*S11IM1 + D12I*S12IM1)/DEN
      T12(1) = (D22I*S12IM1 - D12I*S22IM1)/DEN
      T21(1) = -(D12I*S11IM1 + D11I*S12IM1)/DEN
      T22(1) = (D11I*S22IM1 - D12I*S12IM1)/DEN
      R1 = R1*Y(1)/DEN
      YS(1) = D22I*R1
      YP(1) = -D12I*R1
C
C I = 2,...,N-1:
C
      DO 1 I = 2,NM1
        IM1 = I - 1
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
        R1 = R1*Y(I) - S11IM1*YS(IM1) + S12IM1*YP(IM1)
        R2 = -S12IM1*YS(IM1) - S22IM1*YP(IM1)
        YS(I) = (D22I*R1 - D12I*R2)/DEN
        YP(I) = (D11I*R2 - D12I*R1)/DEN
        DIM1 = DI
        S22IM1 = S22I
        S12IM1 = S12I
        S11IM1 = S11I
    1   CONTINUE
C
C I = N:
C
      R1 = PP*W(NN)
      D11I = R1 - S11IM1 - (S11IM1*T11(NM1)-S12IM1*T21(NM1))
      D12I = -S12IM1 - (S11IM1*T12(NM1) - S12IM1*T22(NM1))
      D22I = DIM1 - (S12IM1*T12(NM1) + S22IM1*T22(NM1))
      DEN = D11I*D22I - D12I*D12I
      R1 = R1*Y(NN) - S11IM1*YS(NM1) + S12IM1*YP(NM1)
      R2 = -S12IM1*YS(NM1) - S22IM1*YP(NM1)
      YS(NN) = (D22I*R1 - D12I*R2)/DEN
      YP(NN) = (D11I*R2 - D12I*R1)/DEN
C
C Back solve the system.
C
      DO 2 I = NM1,1,-1
        YS(I) = YS(I) - (T11(I)*YS(I+1) + T12(I)*YP(I+1))
        YP(I) = YP(I) - (T21(I)*YS(I+1) + T22(I)*YP(I+1))
    2   CONTINUE
      IER = 0
      RETURN
      END
