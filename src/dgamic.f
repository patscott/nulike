      DOUBLE PRECISION FUNCTION DGAMIC(A,X)
C***BEGIN PROLOGUE  DGAMIC
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  COMPLEMENTARY,COMPLEMENTARY INCOMPLETE GAMMA FUNCTION,
C             DOUBLE PRECISION,GAMMA,GAMMA FUNCTION,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Calculates the d.p. complementary incomplete Gamma function
C***DESCRIPTION
C
C Evaluate the complementary incomplete gamma function
C
C DGAMIC = integral from T = X to infinity of EXP(-T) * T**(A-1.)  .
C
C DGAMIC is evaluated for arbitrary real values of A and for non-
C negative values of X (even though DGAMIC is defined for X .LT. 0.0),
C except that for X = 0 and A .LE. 0.0, DGAMIC is undefined.  The
C function and both arguments are double precision.
C
C      A slight deterioration of 2 or 3 digits accuracy will occur when
C DGAMIC is very large or very small in absolute value, because log-
C arithmic variables are used.  Also, if the parameter A is very close
C to a negative INTEGER (but not a negative integer), there is a loss
C of accuracy, which is reported if the result is less than half
C machine precision.
C
C Ref. -- W. Gautschi, An Evaluation Procedure for Incomplete Gamma
C Functions, ACM Trans. Math. Software.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,D9GMIC,D9GMIT,D9LGIC,D9LGIT,DINT,DLGAMS,
C                    DLNGAM,XERCLR,XERROR
C***END PROLOGUE  DGAMIC
      DOUBLE PRECISION A, X, AEPS, AINTA, ALGAP1, ALNEPS, ALNGS, ALX,
     1  BOT, E, EPS, GSTAR, H, SGA, SGNG, SGNGAM, SGNGS, SQEPS, T,
     2  D1MACH, DLNGAM, D9GMIC, D9GMIT, D9LGIC, D9LGIT, DINT
      DATA EPS, SQEPS, ALNEPS, BOT / 4*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DGAMIC
      IF (EPS.NE.0.D0) GO TO 10
      EPS = 0.5D0*D1MACH(3)
      SQEPS = DSQRT (D1MACH(4))
      ALNEPS = -DLOG (D1MACH(3))
      BOT = DLOG (D1MACH(1))
C
 10   IF (X.LT.0.D0) CALL XERROR ( 'DGAMIC  X IS NEGATIVE', 21, 2, 2)
C
      IF (X.GT.0.D0) GO TO 20
      IF (A.LE.0.D0) CALL XERROR ( 'DGAMIC  X = 0 AND A LE 0 SO DGAMIC I
     1S UNDEFINED', 47, 3, 2)
C
      DGAMIC = DEXP (DLNGAM(A+1.D0) - DLOG(A))
      RETURN
C
 20   ALX = DLOG (X)
      SGA = 1.0D0
      IF (A.NE.0.D0) SGA = DSIGN (1.0D0, A)
      AINTA = DINT (A + 0.5D0*SGA)
      AEPS = A - AINTA
C
      IZERO = 0
      IF (X.GE.1.0D0) GO TO 40
C
      IF (A.GT.0.5D0 .OR. DABS(AEPS).GT.0.001D0) GO TO 30
      E = 2.0D0
      IF (-AINTA.GT.1.D0) E = 2.D0*(-AINTA+2.D0)/(AINTA*AINTA-1.0D0)
      E = E - ALX * X**(-0.001D0)
      IF (E*DABS(AEPS).GT.EPS) GO TO 30
C
      DGAMIC = D9GMIC (A, X, ALX)
      RETURN
C
 30   CALL DLGAMS (A+1.0D0, ALGAP1, SGNGAM)
      GSTAR = D9GMIT (A, X, ALGAP1, SGNGAM, ALX)
      IF (GSTAR.EQ.0.D0) IZERO = 1
      IF (GSTAR.NE.0.D0) ALNGS = DLOG (DABS(GSTAR))
      IF (GSTAR.NE.0.D0) SGNGS = DSIGN (1.0D0, GSTAR)
      GO TO 50
C
 40   IF (A.LT.X) DGAMIC = DEXP (D9LGIC(A, X, ALX))
      IF (A.LT.X) RETURN
C
      SGNGAM = 1.0D0
      ALGAP1 = DLNGAM (A+1.0D0)
      SGNGS = 1.0D0
      ALNGS = D9LGIT (A, X, ALGAP1)
C
C EVALUATION OF DGAMIC(A,X) IN TERMS OF TRICOMI-S INCOMPLETE GAMMA FN.
C
 50   H = 1.D0
      IF (IZERO.EQ.1) GO TO 60
C
      T = A*ALX + ALNGS
      IF (T.GT.ALNEPS) GO TO 70
      IF (T.GT.(-ALNEPS)) H = 1.0D0 - SGNGS*DEXP(T)
C
      IF (DABS(H).LT.SQEPS) CALL XERCLR
      IF (DABS(H).LT.SQEPS) CALL XERROR ( 'DGAMIC  RESULT LT HALF PRECIS
     1ION', 32, 1, 1)
C
 60   SGNG = DSIGN (1.0D0, H) * SGA * SGNGAM
      T = DLOG(DABS(H)) + ALGAP1 - DLOG(DABS(A))
      IF (T.LT.BOT) CALL XERCLR
      DGAMIC = SGNG * DEXP(T)
      RETURN
C
 70   SGNG = -SGNGS * SGA * SGNGAM
      T = T + ALGAP1 - DLOG(DABS(A))
      IF (T.LT.BOT) CALL XERCLR
      DGAMIC = SGNG * DEXP(T)
      RETURN
C
      END
      DOUBLE PRECISION FUNCTION D9LGIC(A,X,ALX)
C***BEGIN PROLOGUE  D9LGIC
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  DOUBLE PRECISION,GAMMA,INCOMPLETE GAMMA FUNCTION,
C             LOGARITHM INCOMPLETE GAMMA FUNCTION,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. log incomplete Gamma function for large X
C            and for A .LE. X.
C***DESCRIPTION
C
C Compute the log complementary incomplete gamma function for large X
C and for A .LE. X.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,XERROR
C***END PROLOGUE  D9LGIC
      DOUBLE PRECISION A, X, ALX, EPS, FK, P, R, S, T, XMA, XPA, D1MACH
      DATA EPS / 0.D0 /
C***FIRST EXECUTABLE STATEMENT  D9LGIC
      IF (EPS.EQ.0.D0) EPS = 0.5D0*D1MACH(3)
C
      XPA = X + 1.0D0 - A
      XMA = X - 1.D0 - A
C
      R = 0.D0
      P = 1.D0
      S = P
      DO 10 K=1,300
        FK = K
        T = FK*(A-FK)*(1.D0+R)
        R = -T/((XMA+2.D0*FK)*(XPA+2.D0*FK)+T)
        P = R*P
        S = S + P
        IF (DABS(P).LT.EPS*S) GO TO 20
 10   CONTINUE
      CALL XERROR ( 'D9LGIC  NO CONVERGENCE IN 300 TERMS OF CONTINUED FR
     1ACTION', 57, 1, 2)
C
 20   D9LGIC = A*ALX - X + DLOG(S/XPA)
C
      RETURN
      END
      SUBROUTINE DLGAMS(X,DLGAM,SGNGAM)
C***BEGIN PROLOGUE  DLGAMS
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7A
C***KEYWORDS  ABSOLUTE VALUE,DOUBLE PRECISION,GAMMA FUNCTION,LOGARITHM,
C             SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Calculates the log of the absolute value of the Gamma
C            function
C***DESCRIPTION
C
C DLGAMS(X,DLGAM,SGNGAM) calculates the double precision natural
C logarithm of the absolute value of the gamma function for
C double precision argument X and stores the result in double
C precision argument DLGAM.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DINT,DLNGAM
C***END PROLOGUE  DLGAMS
      DOUBLE PRECISION X, DLGAM, SGNGAM, DINT, DLNGAM
C***FIRST EXECUTABLE STATEMENT  DLGAMS
      DLGAM = DLNGAM(X)
      SGNGAM = 1.0D0
      IF (X.GT.0.D0) RETURN
C
      INT = DMOD (-DINT(X), 2.0D0) + 0.1D0
      IF (INT.EQ.0) SGNGAM = -1.0D0
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION D9LGIT(A,X,ALGAP1)
C***BEGIN PROLOGUE  D9LGIT
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  DOUBLE PRECISION,GAMMA,INCOMPLETE GAMMA FUNCTION,
C             LOGARITHM,SPECIAL FUNCTION,TRICOMI
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes  the log of Tricomi's incomplete Gamma function
C            with Perron's continued fraction for large X and A .GE. X.
C***DESCRIPTION
C
C Compute the log of Tricomi's incomplete gamma function with Perron's
C continued fraction for large X and for A .GE. X.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,XERROR
C***END PROLOGUE  D9LGIT
      DOUBLE PRECISION A, X, ALGAP1, AX, A1X, EPS, FK, HSTAR, P, R, S,
     1  SQEPS, T, D1MACH
      DATA EPS, SQEPS / 2*0.D0 /
C***FIRST EXECUTABLE STATEMENT  D9LGIT
      IF (EPS.NE.0.D0) GO TO 10
      EPS = 0.5D0*D1MACH(3)
      SQEPS = DSQRT (D1MACH(4))
C
 10   IF (X.LE.0.D0 .OR. A.LT.X) CALL XERROR ( 'D9LGIT  X SHOULD BE GT 0
     1.0 AND LE A', 35, 2, 2)
C
      AX = A + X
      A1X = AX + 1.0D0
      R = 0.D0
      P = 1.D0
      S = P
      DO 20 K=1,200
        FK = K
        T = (A+FK)*X*(1.D0+R)
        R = T/((AX+FK)*(A1X+FK)-T)
        P = R*P
        S = S + P
        IF (DABS(P).LT.EPS*S) GO TO 30
 20   CONTINUE
      CALL XERROR ( 'D9LGIT  NO CONVERGENCE IN 200 TERMS OF CONTINUED FR
     1ACTION', 57, 3, 2)
C
 30   HSTAR = 1.0D0 - X*S/A1X
      IF (HSTAR.LT.SQEPS) CALL XERROR ( 'D9LGIT  RESULT LESS THAN HALF P
     1RECISION', 39, 1, 1)
C
      D9LGIT = -X - ALGAP1 - DLOG(HSTAR)
      RETURN
C
      END
      DOUBLE PRECISION FUNCTION D9GMIT(A,X,ALGAP1,SGNGAM,ALX)
C***BEGIN PROLOGUE  D9GMIT
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  COMPLEMENTARY,COMPLEMENTARY INCOMPLETE GAMMA FUNCTION,
C             DOUBLE PRECISION,GAMMA,GAMMA FUNCTION,SPECIAL FUNCTON,
C             TRICOMI
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes d.p. Tricomi-s incomplete Gamma function for
C            small X.
C***DESCRIPTION
C
C Compute Tricomi's incomplete gamma function for small X.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DLNGAM,XERROR
C***END PROLOGUE  D9GMIT
      DOUBLE PRECISION A, X, ALGAP1, SGNGAM, ALX, AE, AEPS, ALGS, ALG2,
     1  BOT, EPS, FK, S, SGNG2, T, TE, D1MACH, DLNGAM
      DATA EPS, BOT / 2*0.D0 /
C***FIRST EXECUTABLE STATEMENT  D9GMIT
      IF (EPS.NE.0.D0) GO TO 10
      EPS = 0.5D0*D1MACH(3)
      BOT = DLOG (D1MACH(1))
C
 10   IF (X.LE.0.D0) CALL XERROR ( 'D9GMIT  X SHOULD BE GT 0', 24, 1, 2)
C
      MA = A + 0.5D0
      IF (A.LT.0.D0) MA = A - 0.5D0
      AEPS = A - DBLE(FLOAT(MA))
C
      AE = A
      IF (A.LT.(-0.5D0)) AE = AEPS
C
      T = 1.D0
      TE = AE
      S = T
      DO 20 K=1,200
        FK = K
        TE = -X*TE/FK
        T = TE/(AE+FK)
        S = S + T
        IF (DABS(T).LT.EPS*DABS(S)) GO TO 30
 20   CONTINUE
      CALL XERROR ( 'D9GMIT  NO CONVERGENCE IN 200 TERMS OF TAYLOR-S SER
     1IES', 54, 2, 2)
C
 30   IF (A.GE.(-0.5D0)) ALGS = -ALGAP1 + DLOG(S)
      IF (A.GE.(-0.5D0)) GO TO 60
C
      ALGS = -DLNGAM(1.D0+AEPS) + DLOG(S)
      S = 1.0D0
      M = -MA - 1
      IF (M.EQ.0) GO TO 50
      T = 1.0D0
      DO 40 K=1,M
        T = X*T/(AEPS-DBLE(FLOAT(M+1-K)))
        S = S + T
        IF (DABS(T).LT.EPS*DABS(S)) GO TO 50
 40   CONTINUE
C
 50   D9GMIT = 0.0D0
      ALGS = -DBLE(FLOAT(MA))*DLOG(X) + ALGS
      IF (S.EQ.0.D0 .OR. AEPS.EQ.0.D0) GO TO 60
C
      SGNG2 = SGNGAM * DSIGN (1.0D0, S)
      ALG2 = -X - ALGAP1 + DLOG(DABS(S))
C
      IF (ALG2.GT.BOT) D9GMIT = SGNG2 * DEXP(ALG2)
      IF (ALGS.GT.BOT) D9GMIT = D9GMIT + DEXP(ALGS)
      RETURN
C
 60   D9GMIT = DEXP (ALGS)
      RETURN
C
      END
      DOUBLE PRECISION FUNCTION D9GMIC(A,X,ALX)
C***BEGIN PROLOGUE  D9GMIC
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  COMPLEMENTARY,COMPLEMENTARY INCOMPLETE GAMMA FUNCTION,
C             DOUBLE PRECISION,GAMMA,GAMMA FUNCTION,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Calculates the d.p. incomplete Gamma function for A near a
C            negative integer and X small.
C***DESCRIPTION
C
C Compute the complementary incomplete gamma function for A near
C a negative integer and for small X.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DLNGAM,XERROR
C***END PROLOGUE  D9GMIC
      DOUBLE PRECISION A, X, ALX, ALNG, BOT, EPS, EULER, FK, FKP1, FM,
     1  S, SGNG, T, TE, D1MACH, DLNGAM
      DATA EULER / 0.5772156649 0153286060 6512090082 40 D0 /
      DATA EPS, BOT / 2*0.D0 /
C***FIRST EXECUTABLE STATEMENT  D9GMIC
      IF (EPS.NE.0.D0) GO TO 10
      EPS = 0.5D0*D1MACH(3)
      BOT = DLOG (D1MACH(1))
C
 10   IF (A.GT.0.D0) CALL XERROR ( 'D9GMIC  A MUST BE NEAR A NEGATIVE IN
     1TEGER', 41, 2, 2)
      IF (X.LE.0.D0) CALL XERROR ( 'D9GMIC  X MUST BE GT ZERO', 25, 3,2)
C
      M = -(A - 0.5D0)
      FM = M
C
      TE = 1.0D0
      T = 1.0D0
      S = T
      DO 20 K=1,200
        FKP1 = K + 1
        TE = -X*TE/(FM+FKP1)
        T = TE/FKP1
        S = S + T
        IF (DABS(T).LT.EPS*S) GO TO 30
 20   CONTINUE
      CALL XERROR ( 'D9GMIC  NO CONVERGENCE IN 200 TERMS OF CONTINUED FR
     1ACTION', 57, 4, 2)
C
 30   D9GMIC = -ALX - EULER + X*S/(FM+1.0D0)
      IF (M.EQ.0) RETURN
C
      IF (M.EQ.1) D9GMIC = -D9GMIC - 1.D0 + 1.D0/X
      IF (M.EQ.1) RETURN
C
      TE = FM
      T = 1.D0
      S = T
      MM1 = M - 1
      DO 40 K=1,MM1
        FK = K
        TE = -X*TE/FK
        T = TE/(FM-FK)
        S = S + T
        IF (DABS(T).LT.EPS*DABS(S)) GO TO 50
 40   CONTINUE
C
 50   DO 60 K=1,M
        D9GMIC = D9GMIC + 1.0D0/DBLE(FLOAT(K))
 60   CONTINUE
C
      SGNG = 1.0D0
      IF (MOD(M,2).EQ.1) SGNG = -1.0D0
      ALNG = DLOG(D9GMIC) - DLNGAM(FM+1.D0)
C
      D9GMIC = 0.D0
      IF (ALNG.GT.BOT) D9GMIC = SGNG * DEXP(ALNG)
      IF (S.NE.0.D0) D9GMIC = D9GMIC +
     1  DSIGN (DEXP(-FM*ALX+DLOG(DABS(S)/FM)), S)
C
      IF (D9GMIC.EQ.0.D0 .AND. S.EQ.0.D0) CALL XERROR (  'D9GMIC  RESULT
     1 UNDERFLOWS', 25, 1, 1)
      RETURN
C
      END
      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
C***BEGIN PROLOGUE  XERROR
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Processes an error (diagnostic) message.
C***DESCRIPTION
C     Abstract
C        XERROR processes a diagnostic message, in a manner
C        determined by the value of LEVEL and the current value
C        of the library error control flag, KONTRL.
C        (See subroutine XSETF for details.)
C
C     Description of Parameters
C      --Input--
C        MESSG - the Hollerith message to be processed, containing
C                no more than 72 characters.
C        NMESSG- the actual number of characters in MESSG.
C        NERR  - the error number associated with this message.
C                NERR must not be zero.
C        LEVEL - error category.
C                =2 means this is an unconditionally fatal error.
C                =1 means this is a recoverable error.  (I.e., it is
C                   non-fatal if XSETF has been appropriately called.)
C                =0 means this is a warning message only.
C                =-1 means this is a warning message which is to be
C                   printed at most once, regardless of how many
C                   times this call is executed.
C
C     Examples
C        CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2)
C        CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.',
C                    43,2,1)
C        CALL XERROR('ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL F
C    1ULLY COLLAPSED.',65,3,0)
C        CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1)
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  XERRWV
C***END PROLOGUE  XERROR
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERROR
      CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
      END
      SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
C***BEGIN PROLOGUE  XERRWV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Processes error message allowing 2 integer and two real
C            values to be included in the message.
C***DESCRIPTION
C     Abstract
C        XERRWV processes a diagnostic message, in a manner
C        determined by the value of LEVEL and the current value
C        of the library error control flag, KONTRL.
C        (See subroutine XSETF for details.)
C        In addition, up to two integer values and two real
C        values may be printed along with the message.
C
C     Description of Parameters
C      --Input--
C        MESSG - the Hollerith message to be processed.
C        NMESSG- the actual number of characters in MESSG.
C        NERR  - the error number associated with this message.
C                NERR must not be zero.
C        LEVEL - error category.
C                =2 means this is an unconditionally fatal error.
C                =1 means this is a recoverable error.  (I.e., it is
C                   non-fatal if XSETF has been appropriately called.)
C                =0 means this is a warning message only.
C                =-1 means this is a warning message which is to be
C                   printed at most once, regardless of how many
C                   times this call is executed.
C        NI    - number of integer values to be printed. (0 to 2)
C        I1    - first integer value.
C        I2    - second integer value.
C        NR    - number of real values to be printed. (0 to 2)
C        R1    - first real value.
C        R2    - second real value.
C
C     Examples
C        CALL XERRWV('SMOOTH -- NUM (=I1) WAS ZERO.',29,1,2,
C    1   1,NUM,0,0,0.,0.)
C        CALL XERRWV('QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM (
C    1R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN)
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  FDUMP,I1MACH,J4SAVE,XERABT,XERCTL,XERPRT,XERSAV,
C                    XGETUA
C***END PROLOGUE  XERRWV
      CHARACTER*(*) MESSG
      CHARACTER*20 LFIRST
      CHARACTER*37 FORM
      DIMENSION LUN(5)
C     GET FLAGS
C***FIRST EXECUTABLE STATEMENT  XERRWV
      LKNTRL = J4SAVE(2,0,.FALSE.)
      MAXMES = J4SAVE(4,0,.FALSE.)
C     CHECK FOR VALID INPUT
      IF ((NMESSG.GT.0).AND.(NERR.NE.0).AND.
     1    (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
         IF (LKNTRL.GT.0) CALL XERPRT('FATAL ERROR IN...',17)
         CALL XERPRT('XERROR -- INVALID INPUT',23)
         IF (LKNTRL.GT.0) CALL FDUMP
         IF (LKNTRL.GT.0) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.',
     1  29)
         IF (LKNTRL.GT.0) CALL XERSAV(' ',0,0,0,KDUMMY)
         CALL XERABT('XERROR -- INVALID INPUT',23)
         RETURN
   10 CONTINUE
C     RECORD MESSAGE
      JUNK = J4SAVE(1,NERR,.TRUE.)
      CALL XERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
C     LET USER OVERRIDE
      LFIRST = MESSG
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      CALL XERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
C     RESET TO ORIGINAL VALUES
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      LKNTRL = MAX0(-2,MIN0(2,LKNTRL))
      MKNTRL = IABS(LKNTRL)
C     DECIDE WHETHER TO PRINT MESSAGE
      IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
      IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN0(1,MAXMES)))
     1.OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES))
     2.OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1))
     3.OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX0(1,MAXMES)))) GO TO 100
         IF (LKNTRL.LE.0) GO TO 20
            CALL XERPRT(' ',1)
C           INTRODUCTION
            IF (LLEVEL.EQ.(-1)) CALL XERPRT
     1('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
            IF (LLEVEL.EQ.0) CALL XERPRT('WARNING IN...',13)
            IF (LLEVEL.EQ.1) CALL XERPRT
     1      ('RECOVERABLE ERROR IN...',23)
            IF (LLEVEL.EQ.2) CALL XERPRT('FATAL ERROR IN...',17)
   20    CONTINUE
C        MESSAGE
         CALL XERPRT(MESSG,LMESSG)
         CALL XGETUA(LUN,NUNIT)
         ISIZEI = LOG10(FLOAT(I1MACH(9))) + 1.0
         ISIZEF = LOG10(FLOAT(I1MACH(10))**I1MACH(11)) + 1.0
         DO 50 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
            DO 22 I=1,MIN(NI,2)
               WRITE (FORM,21) I,ISIZEI
   21          FORMAT ('(11X,21HIN ABOVE MESSAGE, I',I1,'=,I',I2,')   ')
               IF (I.EQ.1) WRITE (IUNIT,FORM) I1
               IF (I.EQ.2) WRITE (IUNIT,FORM) I2
   22       CONTINUE
            DO 24 I=1,MIN(NR,2)
               WRITE (FORM,23) I,ISIZEF+10,ISIZEF
   23          FORMAT ('(11X,21HIN ABOVE MESSAGE, R',I1,'=,E',
     1         I2,'.',I2,')')
               IF (I.EQ.1) WRITE (IUNIT,FORM) R1
               IF (I.EQ.2) WRITE (IUNIT,FORM) R2
   24       CONTINUE
            IF (LKNTRL.LE.0) GO TO 40
C              ERROR NUMBER
               WRITE (IUNIT,30) LERR
   30          FORMAT (15H ERROR NUMBER =,I10)
   40       CONTINUE
   50    CONTINUE
C        TRACE-BACK
         IF (LKNTRL.GT.0) CALL FDUMP
  100 CONTINUE
      IFATAL = 0
      IF ((LLEVEL.EQ.2).OR.((LLEVEL.EQ.1).AND.(MKNTRL.EQ.2)))
     1IFATAL = 1
C     QUIT HERE IF MESSAGE IS NOT FATAL
      IF (IFATAL.LE.0) RETURN
      IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX0(1,MAXMES))) GO TO 120
C        PRINT REASON FOR ABORT
         IF (LLEVEL.EQ.1) CALL XERPRT
     1   ('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
         IF (LLEVEL.EQ.2) CALL XERPRT
     1   ('JOB ABORT DUE TO FATAL ERROR.',29)
C        PRINT ERROR SUMMARY
         CALL XERSAV(' ',-1,0,0,KDUMMY)
  120 CONTINUE
C     ABORT
      IF ((LLEVEL.EQ.2).AND.(KOUNT.GT.MAX0(1,MAXMES))) LMESSG = 0
      CALL XERABT(MESSG,LMESSG)
      RETURN
      END
      SUBROUTINE XERCLR
C***BEGIN PROLOGUE  XERCLR
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Resets current error number to zero.
C***DESCRIPTION
C     Abstract
C        This routine simply resets the current error number to zero.
C        This may be necessary to do in order to determine that
C        a certain error has occurred again since the last time
C        NUMXER was referenced.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  7 June 1978
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  J4SAVE
C***END PROLOGUE  XERCLR
C***FIRST EXECUTABLE STATEMENT  XERCLR
      JUNK = J4SAVE(1,0,.TRUE.)
      RETURN
      END
      SUBROUTINE XGETUA(IUNITA,N)
C***BEGIN PROLOGUE  XGETUA
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Returns unit number(s) to which error messages are being
C            sent.
C***DESCRIPTION
C     Abstract
C        XGETUA may be called to determine the unit number or numbers
C        to which error messages are being sent.
C        These unit numbers may have been set by a call to XSETUN,
C        or a call to XSETUA, or may be a default value.
C
C     Description of Parameters
C      --Output--
C        IUNIT - an array of one to five unit numbers, depending
C                on the value of N.  A value of zero refers to the
C                default unit, as defined by the I1MACH machine
C                constant routine.  Only IUNIT(1),...,IUNIT(N) are
C                defined by XGETUA.  The values of IUNIT(N+1),...,
C                IUNIT(5) are not defined (for N .LT. 5) or altered
C                in any way by XGETUA.
C        N     - the number of units to which copies of the
C                error messages are being sent.  N will be in the
C                range from 1 to 5.
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  J4SAVE
C***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
C***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
      SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
C***BEGIN PROLOGUE  XERCTL
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Allows user control over handling of individual errors.
C***DESCRIPTION
C     Abstract
C        Allows user control over handling of individual errors.
C        Just after each message is recorded, but before it is
C        processed any further (i.e., before it is printed or
C        a decision to abort is made), a call is made to XERCTL.
C        If the user has provided his own version of XERCTL, he
C        can then override the value of KONTROL used in processing
C        this message by redefining its value.
C        KONTRL may be set to any value from -2 to 2.
C        The meanings for KONTRL are the same as in XSETF, except
C        that the value of KONTRL changes only for this message.
C        If KONTRL is set to a value outside the range from -2 to 2,
C        it will be moved back into that range.
C
C     Description of Parameters
C
C      --Input--
C        MESSG1 - the first word (only) of the error message.
C        NMESSG - same as in the call to XERROR or XERRWV.
C        NERR   - same as in the call to XERROR or XERRWV.
C        LEVEL  - same as in the call to XERROR or XERRWV.
C        KONTRL - the current value of the control flag as set
C                 by a call to XSETF.
C
C      --Output--
C        KONTRL - the new value of KONTRL.  If KONTRL is not
C                 defined, it will remain at its original value.
C                 This changed value of control affects only
C                 the current occurrence of the current message.
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERCTL
      CHARACTER*20 MESSG1
C***FIRST EXECUTABLE STATEMENT  XERCTL
      RETURN
      END
      SUBROUTINE FDUMP
C***BEGIN PROLOGUE  FDUMP
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Symbolic dump (should be locally written).
C***DESCRIPTION
C        ***Note*** Machine Dependent Routine
C        FDUMP is intended to be replaced by a locally written
C        version which produces a symbolic dump.  Failing this,
C        it should be replaced by a version which prints the
C        subprogram nesting list.  Note that this dump must be
C        printed on each of up to five files, as indicated by the
C        XGETUA routine.  See XSETUA and XGETUA for details.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  23 May 1979
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
      DOUBLE PRECISION FUNCTION DLNGAM(X)
C***BEGIN PROLOGUE  DLNGAM
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7A
C***KEYWORDS  ABSOLUTE VALUE,DOUBLE PRECISION,GAMMA FUNCTION,LOGARITHM,
C             SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. logarithm of the absolute value of the
C            Gamma function
C***DESCRIPTION
C
C DLNGAM(X) calculates the double precision logarithm of the
C absolute value of the gamma function for double precision
C argument X.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,D9LGMC,DGAMMA,DINT,XERROR
C***END PROLOGUE  DLNGAM
      EXTERNAL DGAMMA
      DOUBLE PRECISION X, DXREL, PI, SINPIY, SQPI2L, SQ2PIL, XMAX,
     1  Y, DINT, DGAMMA, D9LGMC, D1MACH,TEMP
      DATA SQ2PIL / 0.9189385332 0467274178 0329736405 62 D0 /
      DATA SQPI2L / +.2257913526 4472743236 3097614947 441 D+0    /
      DATA PI / 3.1415926535 8979323846 2643383279 50 D0 /
      DATA XMAX, DXREL / 2*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DLNGAM
      IF (XMAX.NE.0.D0) GO TO 10
      TEMP = 1.0D0/DLOG(D1MACH(2))
      XMAX = TEMP * D1MACH(2)
      DXREL = DSQRT (D1MACH(4))
C
 10   Y = DABS (X)
      IF (Y.GT.10.D0) GO TO 20
C
C DLOG (DABS (DGAMMA(X)) ) FOR DABS(X) .LE. 10.0
C
      DLNGAM = DLOG (DABS (DGAMMA(X)) )
      RETURN
C
C DLOG ( DABS (DGAMMA(X)) ) FOR DABS(X) .GT. 10.0
C
 20   IF (Y.GT.XMAX) CALL XERROR ( 'DLNGAM  DABS(X) SO BIG DLNGAM OVERFL
     1OWS', 39, 2, 2)
C
      IF (X.GT.0.D0) DLNGAM = SQ2PIL + (X-0.5D0)*DLOG(X) - X + D9LGMC(Y)
      IF (X.GT.0.D0) RETURN
C
      SINPIY = DABS (DSIN(PI*Y))
      IF (SINPIY.EQ.0.D0) CALL XERROR ( 'DLNGAM  X IS A NEGATIVE INTEGER
     1', 31, 3, 2)
C
      IF (DABS ((X-DINT(X-0.5D0))/X).LT.DXREL) CALL XERROR ( 'DLNGAM  AN
     1SWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER', 68, 1
     2, 1)
C
      DLNGAM = SQPI2L + (X-0.5D0)*DLOG(Y) - X - DLOG(SINPIY) - D9LGMC(Y)
      RETURN
C
      END
      DOUBLE PRECISION FUNCTION DGAMMA(X)
C***BEGIN PROLOGUE  DGAMMA
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7A
C***KEYWORDS  COMPLETE GAMMA FUNCTION,DOUBLE PRECISION,GAMMA FUNCTION,
C             SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. complete Gamma function.
C***DESCRIPTION
C
C DGAMMA(X) calculates the double precision complete gamma function
C for double precision argument X.
C
C Series for GAM        on the interval  0.          to  1.00000E+00
C                                        with weighted error   5.79E-32
C                                         log weighted error  31.24
C                               significant figures required  30.00
C                                    decimal places required  32.05
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,D9LGMC,DCSEVL,DGAMLM,DINT,INITDS,XERROR
C***END PROLOGUE  DGAMMA
      DOUBLE PRECISION X, GAMCS(42), DXREL, PI, SINPIY, SQ2PIL, XMAX,
     1  XMIN, Y,  DINT, D9LGMC, DCSEVL, D1MACH
C
      DATA GAM CS(  1) / +.8571195590 9893314219 2006239994 2 D-2      /
      DATA GAM CS(  2) / +.4415381324 8410067571 9131577165 2 D-2      /
      DATA GAM CS(  3) / +.5685043681 5993633786 3266458878 9 D-1      /
      DATA GAM CS(  4) / -.4219835396 4185605010 1250018662 4 D-2      /
      DATA GAM CS(  5) / +.1326808181 2124602205 8400679635 2 D-2      /
      DATA GAM CS(  6) / -.1893024529 7988804325 2394702388 6 D-3      /
      DATA GAM CS(  7) / +.3606925327 4412452565 7808221722 5 D-4      /
      DATA GAM CS(  8) / -.6056761904 4608642184 8554829036 5 D-5      /
      DATA GAM CS(  9) / +.1055829546 3022833447 3182350909 3 D-5      /
      DATA GAM CS( 10) / -.1811967365 5423840482 9185589116 6 D-6      /
      DATA GAM CS( 11) / +.3117724964 7153222777 9025459316 9 D-7      /
      DATA GAM CS( 12) / -.5354219639 0196871408 7408102434 7 D-8      /
      DATA GAM CS( 13) / +.9193275519 8595889468 8778682594 0 D-9      /
      DATA GAM CS( 14) / -.1577941280 2883397617 6742327395 3 D-9      /
      DATA GAM CS( 15) / +.2707980622 9349545432 6654043308 9 D-10     /
      DATA GAM CS( 16) / -.4646818653 8257301440 8166105893 3 D-11     /
      DATA GAM CS( 17) / +.7973350192 0074196564 6076717535 9 D-12     /
      DATA GAM CS( 18) / -.1368078209 8309160257 9949917230 9 D-12     /
      DATA GAM CS( 19) / +.2347319486 5638006572 3347177168 8 D-13     /
      DATA GAM CS( 20) / -.4027432614 9490669327 6657053469 9 D-14     /
      DATA GAM CS( 21) / +.6910051747 3721009121 3833697525 7 D-15     /
      DATA GAM CS( 22) / -.1185584500 2219929070 5238712619 2 D-15     /
      DATA GAM CS( 23) / +.2034148542 4963739552 0102605193 2 D-16     /
      DATA GAM CS( 24) / -.3490054341 7174058492 7401294910 8 D-17     /
      DATA GAM CS( 25) / +.5987993856 4853055671 3505106602 6 D-18     /
      DATA GAM CS( 26) / -.1027378057 8722280744 9006977843 1 D-18     /
      DATA GAM CS( 27) / +.1762702816 0605298249 4275966074 8 D-19     /
      DATA GAM CS( 28) / -.3024320653 7353062609 5877211204 2 D-20     /
      DATA GAM CS( 29) / +.5188914660 2183978397 1783355050 6 D-21     /
      DATA GAM CS( 30) / -.8902770842 4565766924 4925160106 6 D-22     /
      DATA GAM CS( 31) / +.1527474068 4933426022 7459689130 6 D-22     /
      DATA GAM CS( 32) / -.2620731256 1873629002 5732833279 9 D-23     /
      DATA GAM CS( 33) / +.4496464047 8305386703 3104657066 6 D-24     /
      DATA GAM CS( 34) / -.7714712731 3368779117 0390152533 3 D-25     /
      DATA GAM CS( 35) / +.1323635453 1260440364 8657271466 6 D-25     /
      DATA GAM CS( 36) / -.2270999412 9429288167 0231381333 3 D-26     /
      DATA GAM CS( 37) / +.3896418998 0039914493 2081663999 9 D-27     /
      DATA GAM CS( 38) / -.6685198115 1259533277 9212799999 9 D-28     /
      DATA GAM CS( 39) / +.1146998663 1400243843 4761386666 6 D-28     /
      DATA GAM CS( 40) / -.1967938586 3451346772 9510399999 9 D-29     /
      DATA GAM CS( 41) / +.3376448816 5853380903 3489066666 6 D-30     /
      DATA GAM CS( 42) / -.5793070335 7821357846 2549333333 3 D-31     /
      DATA PI / 3.1415926535 8979323846 2643383279 50 D0 /
      DATA SQ2PIL / 0.9189385332 0467274178 0329736405 62 D0 /
      DATA NGAM, XMIN, XMAX, DXREL / 0, 3*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DGAMMA
      IF (NGAM.NE.0) GO TO 10
      NGAM = INITDS (GAMCS, 42, 0.1*SNGL(D1MACH(3)) )
C
      CALL DGAMLM (XMIN, XMAX)
      DXREL = DSQRT (D1MACH(4))
C
 10   Y = DABS(X)
      IF (Y.GT.10.D0) GO TO 50
C
C COMPUTE GAMMA(X) FOR -XBND .LE. X .LE. XBND.  REDUCE INTERVAL AND FIND
C GAMMA(1+Y) FOR 0.0 .LE. Y .LT. 1.0 FIRST OF ALL.
C
      N = X
      IF (X.LT.0.D0) N = N - 1
      Y = X - DBLE(FLOAT(N))
      N = N - 1
      DGAMMA = 0.9375D0 + DCSEVL (2.D0*Y-1.D0, GAMCS, NGAM)
      IF (N.EQ.0) RETURN
C
      IF (N.GT.0) GO TO 30
C
C COMPUTE GAMMA(X) FOR X .LT. 1.0
C
      N = -N
      IF (X.EQ.0.D0) CALL XERROR ( 'DGAMMA  X IS 0', 14, 4, 2)
      IF (X.LT.0.0 .AND. X+DBLE(FLOAT(N-2)).EQ.0.D0) CALL XERROR ( 'DGAM
     1MA  X IS A NEGATIVE INTEGER', 31, 4, 2)
      IF (X.LT.(-0.5D0) .AND. DABS((X-DINT(X-0.5D0))/X).LT.DXREL) CALL
     1  XERROR ( 'DGAMMA  ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NE
     2GATIVE INTEGER', 68, 1, 1)
C
      DO 20 I=1,N
        DGAMMA = DGAMMA/(X+DBLE(FLOAT(I-1)) )
 20   CONTINUE
      RETURN
C
C GAMMA(X) FOR X .GE. 2.0 AND X .LE. 10.0
C
 30   DO 40 I=1,N
        DGAMMA = (Y+DBLE(FLOAT(I))) * DGAMMA
 40   CONTINUE
      RETURN
C
C GAMMA(X) FOR DABS(X) .GT. 10.0.  RECALL Y = DABS(X).
C
 50   IF (X.GT.XMAX) CALL XERROR ( 'DGAMMA  X SO BIG GAMMA OVERFLOWS',
     1  32, 3, 2)
C
      DGAMMA = 0.D0
      IF (X.LT.XMIN) CALL XERROR ( 'DGAMMA  X SO SMALL GAMMA UNDERFLOWS'
     1  , 35, 2, 1)
      IF (X.LT.XMIN) RETURN
C
      DGAMMA = DEXP ((Y-0.5D0)*DLOG(Y) - Y + SQ2PIL + D9LGMC(Y) )
      IF (X.GT.0.D0) RETURN
C
      IF (DABS((X-DINT(X-0.5D0))/X).LT.DXREL) CALL XERROR ( 'DGAMMA  ANS
     1WER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER'  , 61, 1, 1)
C
      SINPIY = DSIN (PI*Y)
      IF (SINPIY.EQ.0.D0) CALL XERROR ( 'DGAMMA  X IS A NEGATIVE INTEGER
     1', 31, 4, 2)
C
      DGAMMA = -PI/(Y*SINPIY*DGAMMA)
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION DCSEVL(X,A,N)
C***BEGIN PROLOGUE  DCSEVL
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C3A2
C***KEYWORDS  CHEBYSHEV,FNLIB,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Evaluate the double precision N-term Chebyshev series A
C            at X.
C***DESCRIPTION
C
C Evaluate the N-term Chebyshev series A at X.  Adapted from
C R. Broucke, Algorithm 446, C.A.C.M., 16, 254 (1973).
C W. Fullerton, C-3, Los Alamos Scientific Laboratory.
C
C       Input Arguments --
C X    double precision value at which the series is to be evaluated.
C A    double precision array of N terms of a Chebyshev series.  In
C      evaluating A, only half of the first coefficient is summed.
C N    number of terms in array A.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  DCSEVL
C
       DOUBLE PRECISION A(N),X,TWOX,B0,B1,B2
C***FIRST EXECUTABLE STATEMENT  DCSEVL
       IF(N.LT.1)CALL XERROR( 'DCSEVL  NUMBER OF TERMS LE 0', 28, 2,2)
       IF(N.GT.1000) CALL XERROR ( 'DCSEVL  NUMBER OF TERMS GT 1000',
     1   31, 3, 2)
       IF ((X.LT.-1.D0) .OR. (X.GT.1.D0)) CALL XERROR ( 'DCSEVL  X OUTSI
     1DE (-1,+1)', 25, 1, 1)
C
       TWOX = 2.0D0*X
       B1 = 0.D0
       B0=0.D0
       DO 10 I=1,N
         B2=B1
         B1=B0
         NI = N - I + 1
         B0 = TWOX*B1 - B2 + A(NI)
 10    CONTINUE
C
       DCSEVL = 0.5D0 * (B0-B2)
C
       RETURN
      END
      SUBROUTINE DGAMLM(XMIN,XMAX)
C***BEGIN PROLOGUE  DGAMLM
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7A,R2
C***KEYWORDS  COMPLETE GAMMA FUNCTION,DOUBLE PRECISION,GAMMA FUNCTION,
C             LIMITS,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. minimum and maximum bounds for X in
C            GAMMA(X).
C***DESCRIPTION
C
C Calculate the minimum and maximum legal bounds for X in gamma(X).
C XMIN and XMAX are not the only bounds, but they are the only non-
C trivial ones to calculate.
C
C             Output Arguments --
C XMIN   double precision minimum legal value of X in gamma(X).  Any
C        smaller value of X might result in underflow.
C XMAX   double precision maximum legal value of X in gamma(X).  Any
C        larger value of X might cause overflow.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,XERROR
C***END PROLOGUE  DGAMLM
      DOUBLE PRECISION XMIN, XMAX, ALNBIG, ALNSML, XLN, XOLD, D1MACH
C***FIRST EXECUTABLE STATEMENT  DGAMLM
      ALNSML = DLOG(D1MACH(1))
      XMIN = -ALNSML
      DO 10 I=1,10
        XOLD = XMIN
        XLN = DLOG(XMIN)
        XMIN = XMIN - XMIN*((XMIN+0.5D0)*XLN - XMIN - 0.2258D0 + ALNSML)
     1    / (XMIN*XLN+0.5D0)
        IF (DABS(XMIN-XOLD).LT.0.005D0) GO TO 20
 10   CONTINUE
      CALL XERROR ( 'DGAMLM  UNABLE TO FIND XMIN', 27, 1, 2)
C
 20   XMIN = -XMIN + 0.01D0
C
      ALNBIG = DLOG (D1MACH(2))
      XMAX = ALNBIG
      DO 30 I=1,10
        XOLD = XMAX
        XLN = DLOG(XMAX)
        XMAX = XMAX - XMAX*((XMAX-0.5D0)*XLN - XMAX + 0.9189D0 - ALNBIG)
     1    / (XMAX*XLN-0.5D0)
        IF (DABS(XMAX-XOLD).LT.0.005D0) GO TO 40
 30   CONTINUE
      CALL XERROR ( 'DGAMLM  UNABLE TO FIND XMAX', 27, 2, 2)
C
 40   XMAX = XMAX - 0.01D0
      XMIN = DMAX1 (XMIN, -XMAX+1.D0)
C
      RETURN
      END
      INTEGER FUNCTION I1MACH(I)
C***BEGIN PROLOGUE  I1MACH
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  840405   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns integer machine dependent constants
C***DESCRIPTION
C
C     This is the CMLIB version of I1MACH, the integer machine
C     constants subroutine originally developed for the PORT library.
C
C     I1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function
C     subroutine with one (input) argument, and can be called
C     as follows, for example
C
C          K = I1MACH(I)
C
C     where I=1,...,16.  The (output) value of K above is
C     determined by the (input) value of I.  The results for
C     various values of I are discussed below.
C
C  I/O unit numbers.
C    I1MACH( 1) = the standard input unit.
C    I1MACH( 2) = the standard output unit.
C    I1MACH( 3) = the standard punch unit.
C    I1MACH( 4) = the standard error message unit.
C
C  Words.
C    I1MACH( 5) = the number of bits per integer storage unit.
C    I1MACH( 6) = the number of characters per integer storage unit.
C
C  Integers.
C    assume integers are represented in the S-digit, base-A form
C
C               sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               where 0 .LE. X(I) .LT. A for I=0,...,S-1.
C    I1MACH( 7) = A, the base.
C    I1MACH( 8) = S, the number of base-A digits.
C    I1MACH( 9) = A**S - 1, the largest magnitude.
C
C  Floating-Point Numbers.
C    Assume floating-point numbers are represented in the T-digit,
C    base-B form
C               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               where 0 .LE. X(I) .LT. B for I=1,...,T,
C               0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
C    I1MACH(10) = B, the base.
C
C  Single-Precision
C    I1MACH(11) = T, the number of base-B digits.
C    I1MACH(12) = EMIN, the smallest exponent E.
C    I1MACH(13) = EMAX, the largest exponent E.
C
C  Double-Precision
C    I1MACH(14) = T, the number of base-B digits.
C    I1MACH(15) = EMIN, the smallest exponent E.
C    I1MACH(16) = EMAX, the largest exponent E.
C
C  To alter this function for a particular environment,
C  the desired set of DATA statements should be activated by
C  removing the C from column 1.  Also, the values of
C  I1MACH(1) - I1MACH(4) should be checked for consistency
C  with the local operating system.
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  I1MACH
C
      INTEGER IMACH(16),OUTPUT
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), SUN SPARCSTATIONS, SILICON GRAPHCS WORKSTATIONS, HP
C     9000 WORKSTATIONS, IBM RS/6000 WORKSTATIONS, AND 8087 BASED
C     MICROS (E.G. IBM PC AND AT&T 6300).
C
C === MACHINE = ATT.3B
C === MACHINE = ATT.6300
C === MACHINE = ATT.7300
C === MACHINE = HP.9000
C === MACHINE = IBM.PC
C === MACHINE = IBM.RS6000
C === MACHINE = IEEE.MOST-SIG-BYTE-FIRST
C === MACHINE = IEEE.LEAST-SIG-BYTE-FIRST
C === MACHINE = SGI
C === MACHINE = SUN
C === MACHINE = 68000
C === MACHINE = 8087
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -125 /
C      DATA IMACH(13) /  128 /
C      DATA IMACH(14) /   53 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /
C
C     MACHINE CONSTANTS FOR SUN WORKSTATIONS.  f77 WITH -r8 OPTION.
C
C === MACHINE = SUN.F77-WITH-R8-OPTION
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    6 /
C      DATA IMACH( 4) /    0 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   53 /
C      DATA IMACH(12) / -1021 /
C      DATA IMACH(13) /  1024 /
C      DATA IMACH(14) /   113 /
C      DATA IMACH(15) / -16382 /
C      DATA IMACH(16) /  16384 /
C
C     MACHINE CONSTANTS FOR SGI Origin 2000 with -r8 -d16 options.
C
C === MACHINE = SGI.ORIGIN.F77-WITH-R8-D16-OPTIONS
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    6 /
C      DATA IMACH( 4) /    0 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   53 /
C      DATA IMACH(12) / -1022 /
C      DATA IMACH(13) /  1024 /
C      DATA IMACH(14) /   107 /
C      DATA IMACH(15) /  -916 /
C      DATA IMACH(16) /  1025 /
C
C     MACHINE CONSTANTS FOR IBM RS/6000 WORKSTATIONS WITH -qautodbl=dblpad.
C
C === MACHINE = IBM.RS6000.XLF-WITH-AUTODBL-OPTION
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    6 /
C      DATA IMACH( 4) /    0 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   53 /
C      DATA IMACH(12) / -1021 /
C      DATA IMACH(13) /  1024 /
C      DATA IMACH(14) /   114 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C === MACHINE = AMDAHL
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  63 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  63 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C === MACHINE = BURROUGHS.1700
C      DATA IMACH( 1) /    7 /
C      DATA IMACH( 2) /    2 /
C      DATA IMACH( 3) /    2 /
C      DATA IMACH( 4) /    2 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   33 /
C      DATA IMACH( 9) / Z1FFFFFFFF /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -256 /
C      DATA IMACH(13) /  255 /
C      DATA IMACH(14) /   60 /
C      DATA IMACH(15) / -256 /
C      DATA IMACH(16) /  255 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C === MACHINE = BURROUGHS.5700
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  48 /
C      DATA IMACH( 6) /   6 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  39 /
C      DATA IMACH( 9) / O0007777777777777 /
C      DATA IMACH(10) /   8 /
C      DATA IMACH(11) /  13 /
C      DATA IMACH(12) / -50 /
C      DATA IMACH(13) /  76 /
C      DATA IMACH(14) /  26 /
C      DATA IMACH(15) / -50 /
C      DATA IMACH(16) /  76 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C === MACHINE = BURROUGHS.6700
C === MACHINE = BURROUGHS.7700
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  48 /
C      DATA IMACH( 6) /   6 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  39 /
C      DATA IMACH( 9) / O0007777777777777 /
C      DATA IMACH(10) /   8 /
C      DATA IMACH(11) /  13 /
C      DATA IMACH(12) / -50 /
C      DATA IMACH(13) /  76 /
C      DATA IMACH(14) /  26 /
C      DATA IMACH(15) / -32754 /
C      DATA IMACH(16) /  32780 /
C
C     MACHINE CONSTANTS FOR THE CONVEX C1, C2, C3 SERIES
C
C === MACHINE = CONVEX
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    0 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   53 /
C      DATA IMACH(15) / -1023 /
C      DATA IMACH(16) /  1023 /
C
C     MACHINE CONSTANTS FOR THE CONVEX C1, C2, C3 SERIES (NATIVE MODE)
C     WITH -P8 OPTION
C
C === MACHINE = CONVEX.P8
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /     0 /
C      DATA IMACH( 4) /     6 /
C      DATA IMACH( 5) /    64 /
C      DATA IMACH( 6) /     8 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    63 /
C      DATA IMACH( 9) / 9223372036854775807 /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    53 /
C      DATA IMACH(12) / -1023 /
C      DATA IMACH(13) /  1023 /
C      DATA IMACH(14) /   112 /
C      DATA IMACH(15) / -16383 /
C      DATA IMACH(16) /  16383 /
C
C     MACHINE CONSTANTS FOR THE CONVEX C1, C2, C3 SERIES (IEEE MODE)
C
C === MACHINE = CONVEX.IEEE
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    0 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -125 /
C      DATA IMACH(13) /  128 /
C      DATA IMACH(14) /   53 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /
C
C     MACHINE CONSTANTS FOR THE CONVEX C1, C2, C3 SERIES(IEEE MODE)
C     WITH -P8 OPTION
C
C === MACHINE = CONVEX.IEEE.P8
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /     0 /
C      DATA IMACH( 4) /     6 /
C      DATA IMACH( 5) /    32 /
C      DATA IMACH( 6) /     4 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    53 /
C      DATA IMACH(12) / -1021 /
C      DATA IMACH(13) /  1024 /
C      DATA IMACH(14) /    53 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /
C
C     MACHINE CONSTANTS FOR THE CYBER 170/180 SERIES USING NOS (FTN5).
C
C === MACHINE = CYBER.170.NOS
C === MACHINE = CYBER.180.NOS
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   60 /
C      DATA IMACH( 6) /   10 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   48 /
C      DATA IMACH( 9) / O"00007777777777777777" /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   48 /
C      DATA IMACH(12) / -974 /
C      DATA IMACH(13) / 1070 /
C      DATA IMACH(14) /   96 /
C      DATA IMACH(15) / -927 /
C      DATA IMACH(16) / 1070 /
C
C     MACHINE CONSTANTS FOR THE CDC 180 SERIES USING NOS/VE
C
C === MACHINE = CYBER.180.NOS/VE
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /     7 /
C      DATA IMACH( 4) /     6 /
C      DATA IMACH( 5) /    64 /
C      DATA IMACH( 6) /     8 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    63 /
C      DATA IMACH( 9) / 9223372036854775807 /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    47 /
C      DATA IMACH(12) / -4095 /
C      DATA IMACH(13) /  4094 /
C      DATA IMACH(14) /    94 /
C      DATA IMACH(15) / -4095 /
C      DATA IMACH(16) /  4094 /
C
C     MACHINE CONSTANTS FOR THE CYBER 205
C
C === MACHINE = CYBER.205
C      DATA IMACH( 1) /      5 /
C      DATA IMACH( 2) /      6 /
C      DATA IMACH( 3) /      7 /
C      DATA IMACH( 4) /      6 /
C      DATA IMACH( 5) /     64 /
C      DATA IMACH( 6) /      8 /
C      DATA IMACH( 7) /      2 /
C      DATA IMACH( 8) /     47 /
C      DATA IMACH( 9) / X'00007FFFFFFFFFFF' /
C      DATA IMACH(10) /      2 /
C      DATA IMACH(11) /     47 /
C      DATA IMACH(12) / -28625 /
C      DATA IMACH(13) /  28718 /
C      DATA IMACH(14) /     94 /
C      DATA IMACH(15) / -28625 /
C      DATA IMACH(16) /  28718 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C === MACHINE = CDC.6000
C === MACHINE = CDC.7000
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   60 /
C      DATA IMACH( 6) /   10 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   48 /
C      DATA IMACH( 9) / 00007777777777777777B /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   48 /
C      DATA IMACH(12) / -974 /
C      DATA IMACH(13) / 1070 /
C      DATA IMACH(14) /   96 /
C      DATA IMACH(15) / -927 /
C      DATA IMACH(16) / 1070 /
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C     USING THE 46 BIT INTEGER COMPILER OPTION
C
C === MACHINE = CRAY.46-BIT-INTEGER
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /   102 /
C      DATA IMACH( 4) /     6 /
C      DATA IMACH( 5) /    64 /
C      DATA IMACH( 6) /     8 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    46 /
C      DATA IMACH( 9) /  777777777777777777777B /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    47 /
C      DATA IMACH(12) / -8189 /
C      DATA IMACH(13) /  8190 /
C      DATA IMACH(14) /    94 /
C      DATA IMACH(15) / -8099 /
C      DATA IMACH(16) /  8190 /
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C     USING THE 64 BIT INTEGER COMPILER OPTION
C
C === MACHINE = CRAY.64-BIT-INTEGER
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /   102 /
C      DATA IMACH( 4) /     6 /
C      DATA IMACH( 5) /    64 /
C      DATA IMACH( 6) /     8 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    63 /
C      DATA IMACH( 9) /  777777777777777777777B /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    47 /
C      DATA IMACH(12) / -8189 /
C      DATA IMACH(13) /  8190 /
C      DATA IMACH(14) /    94 /
C      DATA IMACH(15) / -8099 /
C      DATA IMACH(16) /  8190 /C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C === MACHINE = DATA_GENERAL.ECLIPSE.S/200
C      DATA IMACH( 1) /   11 /
C      DATA IMACH( 2) /   12 /
C      DATA IMACH( 3) /    8 /
C      DATA IMACH( 4) /   10 /
C      DATA IMACH( 5) /   16 /
C      DATA IMACH( 6) /    2 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   15 /
C      DATA IMACH( 9) /32767 /
C      DATA IMACH(10) /   16 /
C      DATA IMACH(11) /    6 /
C      DATA IMACH(12) /  -64 /
C      DATA IMACH(13) /   63 /
C      DATA IMACH(14) /   14 /
C      DATA IMACH(15) /  -64 /
C      DATA IMACH(16) /   63 /
C
C     ELXSI  6400
C
C === MACHINE = ELSXI.6400
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /     6 /
C      DATA IMACH( 4) /     6 /
C      DATA IMACH( 5) /    32 /
C      DATA IMACH( 6) /     4 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    32 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    24 /
C      DATA IMACH(12) /  -126 /
C      DATA IMACH(13) /   127 /
C      DATA IMACH(14) /    53 /
C      DATA IMACH(15) / -1022 /
C      DATA IMACH(16) /  1023 /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7
C
C === MACHINE = HARRIS.220
C === MACHINE = HARRIS.SLASH6
C === MACHINE = HARRIS.SLASH7
C      DATA IMACH( 1) /       5 /
C      DATA IMACH( 2) /       6 /
C      DATA IMACH( 3) /       0 /
C      DATA IMACH( 4) /       6 /
C      DATA IMACH( 5) /      24 /
C      DATA IMACH( 6) /       3 /
C      DATA IMACH( 7) /       2 /
C      DATA IMACH( 8) /      23 /
C      DATA IMACH( 9) / 8388607 /
C      DATA IMACH(10) /       2 /
C      DATA IMACH(11) /      23 /
C      DATA IMACH(12) /    -127 /
C      DATA IMACH(13) /     127 /
C      DATA IMACH(14) /      38 /
C      DATA IMACH(15) /    -127 /
C      DATA IMACH(16) /     127 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C === MACHINE = HONEYWELL.600/6000
C === MACHINE = HONEYWELL.DPS.8/70
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /   43 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   63 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C === MACHINE = HP.2100.3_WORD_DP
C      DATA IMACH(1) /      5/
C      DATA IMACH(2) /      6 /
C      DATA IMACH(3) /      4 /
C      DATA IMACH(4) /      1 /
C      DATA IMACH(5) /     16 /
C      DATA IMACH(6) /      2 /
C      DATA IMACH(7) /      2 /
C      DATA IMACH(8) /     15 /
C      DATA IMACH(9) /  32767 /
C      DATA IMACH(10)/      2 /
C      DATA IMACH(11)/     23 /
C      DATA IMACH(12)/   -128 /
C      DATA IMACH(13)/    127 /
C      DATA IMACH(14)/     39 /
C      DATA IMACH(15)/   -128 /
C      DATA IMACH(16)/    127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C === MACHINE = HP.2100.4_WORD_DP
C      DATA IMACH(1) /      5 /
C      DATA IMACH(2) /      6 /
C      DATA IMACH(3) /      4 /
C      DATA IMACH(4) /      1 /
C      DATA IMACH(5) /     16 /
C      DATA IMACH(6) /      2 /
C      DATA IMACH(7) /      2 /
C      DATA IMACH(8) /     15 /
C      DATA IMACH(9) /  32767 /
C      DATA IMACH(10)/      2 /
C      DATA IMACH(11)/     23 /
C      DATA IMACH(12)/   -128 /
C      DATA IMACH(13)/    127 /
C      DATA IMACH(14)/     55 /
C      DATA IMACH(15)/   -128 /
C      DATA IMACH(16)/    127 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86 AND
C     THE INTERDATA 3230 AND INTERDATA 7/32.
C
C === MACHINE = IBM.360
C === MACHINE = IBM.370
C === MACHINE = XEROX.SIGMA.5
C === MACHINE = XEROX.SIGMA.7
C === MACHINE = XEROX.SIGMA.9
C === MACHINE = SEL.85
C === MACHINE = SEL.86
C === MACHINE = INTERDATA.3230
C === MACHINE = INTERDATA.7/32
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / Z7FFFFFFF /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  63 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  63 /
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C === MACHINE = INTERDATA.8/32.UNIX
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   6 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / Z'7FFFFFFF' /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  62 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  62 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C === MACHINE = PDP-10.KA
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    5 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / "377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   54 /
C      DATA IMACH(15) / -101 /
C      DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C === MACHINE = PDP-10.KI
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    5 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / "377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   62 /
C      DATA IMACH(15) / -128 /
C      DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C === MACHINE = PDP-11.32-BIT
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C === MACHINE = PDP-11.16-BIT
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   16 /
C      DATA IMACH( 6) /    2 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   15 /
C      DATA IMACH( 9) / 32767 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C === MACHINE = SEQUENT.BALANCE.8000
C      DATA IMACH( 1) /     0 /
C      DATA IMACH( 2) /     0 /
C      DATA IMACH( 3) /     7 /
C      DATA IMACH( 4) /     0 /
C      DATA IMACH( 5) /    32 /
C      DATA IMACH( 6) /     1 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    31 /
C      DATA IMACH( 9) /  2147483647 /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    24 /
C      DATA IMACH(12) /  -125 /
C      DATA IMACH(13) /   128 /
C      DATA IMACH(14) /    53 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. FTN COMPILER
C
C === MACHINE = UNIVAC.1100
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    1 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   60 /
C      DATA IMACH(15) /-1024 /
C      DATA IMACH(16) / 1023 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C
C === MACHINE = VAX.11/780
C      DATA IMACH(1) /    5 /
C      DATA IMACH(2) /    6 /
C      DATA IMACH(3) /    5 /
C      DATA IMACH(4) /    6 /
C      DATA IMACH(5) /   32 /
C      DATA IMACH(6) /    4 /
C      DATA IMACH(7) /    2 /
C      DATA IMACH(8) /   31 /
C      DATA IMACH(9) /2147483647 /
C      DATA IMACH(10)/    2 /
C      DATA IMACH(11)/   24 /
C      DATA IMACH(12)/ -127 /
C      DATA IMACH(13)/  127 /
C      DATA IMACH(14)/   56 /
C      DATA IMACH(15)/ -127 /
C      DATA IMACH(16)/  127 /
C
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
C
      I1MACH=IMACH(I)
      RETURN
C
      END
      SUBROUTINE XERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)
C***BEGIN PROLOGUE  XERSAV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Records that an error occurred.
C***DESCRIPTION
C     Abstract
C        Record that this error occurred.
C
C     Description of Parameters
C     --Input--
C       MESSG, NMESSG, NERR, LEVEL are as in XERROR,
C       except that when NMESSG=0 the tables will be
C       dumped and cleared, and when NMESSG is less than zero the
C       tables will be dumped and not cleared.
C     --Output--
C       ICOUNT will be the number of times this message has
C       been seen, or zero if the table has overflowed and
C       does not contain this message specifically.
C       When NMESSG=0, ICOUNT will not be altered.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  19 Mar 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
C***END PROLOGUE  XERSAV
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
      CHARACTER*20 MESTAB(10),MES
      DIMENSION NERTAB(10),LEVTAB(10),KOUNT(10)
      SAVE MESTAB,NERTAB,LEVTAB,KOUNT,KOUNTX
C     NEXT TWO DATA STATEMENTS ARE NECESSARY TO PROVIDE A BLANK
C     ERROR TABLE INITIALLY
      DATA KOUNT(1),KOUNT(2),KOUNT(3),KOUNT(4),KOUNT(5),
     1     KOUNT(6),KOUNT(7),KOUNT(8),KOUNT(9),KOUNT(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      DATA KOUNTX/0/
C***FIRST EXECUTABLE STATEMENT  XERSAV
      IF (NMESSG.GT.0) GO TO 80
C     DUMP THE TABLE
         IF (KOUNT(1).EQ.0) RETURN
C        PRINT TO EACH UNIT
         CALL XGETUA(LUN,NUNIT)
         DO 60 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C           PRINT TABLE HEADER
            WRITE (IUNIT,10)
   10       FORMAT (32H0          ERROR MESSAGE SUMMARY/
     1      51H MESSAGE START             NERR     LEVEL     COUNT)
C           PRINT BODY OF TABLE
            DO 20 I=1,10
               IF (KOUNT(I).EQ.0) GO TO 30
               WRITE (IUNIT,15) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
   15          FORMAT (1X,A20,3I10)
   20       CONTINUE
   30       CONTINUE
C           PRINT NUMBER OF OTHER ERRORS
            IF (KOUNTX.NE.0) WRITE (IUNIT,40) KOUNTX
   40       FORMAT (41H0OTHER ERRORS NOT INDIVIDUALLY TABULATED=,I10)
            WRITE (IUNIT,50)
   50       FORMAT (1X)
   60    CONTINUE
         IF (NMESSG.LT.0) RETURN
C        CLEAR THE ERROR TABLES
         DO 75 I=1,10
   70       KOUNT(I) = 0
   75    CONTINUE
         KOUNTX = 0
         RETURN
   80 CONTINUE
C     PROCESS A MESSAGE...
C     SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C     OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
      MES = MESSG
      DO 90 I=1,10
         II = I
         IF (KOUNT(I).EQ.0) GO TO 110
         IF (MES.NE.MESTAB(I)) GO TO 90
         IF (NERR.NE.NERTAB(I)) GO TO 90
         IF (LEVEL.NE.LEVTAB(I)) GO TO 90
         GO TO 100
   90 CONTINUE
C     THREE POSSIBLE CASES...
C     TABLE IS FULL
         KOUNTX = KOUNTX+1
         ICOUNT = 1
         RETURN
C     MESSAGE FOUND IN TABLE
  100    KOUNT(II) = KOUNT(II) + 1
         ICOUNT = KOUNT(II)
         RETURN
C     EMPTY SLOT FOUND FOR NEW MESSAGE
  110    MESTAB(II) = MES
         NERTAB(II) = NERR
         LEVTAB(II) = LEVEL
         KOUNT(II)  = 1
         ICOUNT = 1
         RETURN
      END
      SUBROUTINE XERABT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERABT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Aborts program execution and prints error message.
C***DESCRIPTION
C     Abstract
C        ***Note*** machine dependent routine
C        XERABT aborts the execution of the program.
C        The error message causing the abort is given in the calling
C        sequence, in case one needs it for printing on a dayfile,
C        for example.
C
C     Description of Parameters
C        MESSG and NMESSG are as in XERROR, except that NMESSG may
C        be zero, in which case no message is being supplied.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  19 MAR 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERABT
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERABT
      STOP
      END
      SUBROUTINE XERPRT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERPRT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Prints error messages.
C***DESCRIPTION
C     Abstract
C        Print the Hollerith message in MESSG, of length NMESSG,
C        on each file indicated by XGETUA.
C     Latest revision ---  19 MAR 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
C***END PROLOGUE  XERPRT
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
C     OBTAIN UNIT NUMBERS AND WRITE LINE TO EACH UNIT
C***FIRST EXECUTABLE STATEMENT  XERPRT
      CALL XGETUA(LUN,NUNIT)
      LENMES = LEN(MESSG)
      DO 20 KUNIT=1,NUNIT
         IUNIT = LUN(KUNIT)
         IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
         DO 10 ICHAR=1,LENMES,72
            LAST = MIN0(ICHAR+71 , LENMES)
            WRITE (IUNIT,'(1X,A)') MESSG(ICHAR:LAST)
   10    CONTINUE
   20 CONTINUE
      RETURN
      END
      FUNCTION J4SAVE(IWHICH,IVALUE,ISET)
C***BEGIN PROLOGUE  J4SAVE
C***REFER TO  XERROR
C     Abstract
C        J4SAVE saves and recalls several global variables needed
C        by the library error handling routines.
C
C     Description of Parameters
C      --Input--
C        IWHICH - Index of item desired.
C                = 1 Refers to current error number.
C                = 2 Refers to current error control flag.
C                 = 3 Refers to current unit number to which error
C                    messages are to be sent.  (0 means use standard.)
C                 = 4 Refers to the maximum number of times any
C                     message is to be printed (as set by XERMAX).
C                 = 5 Refers to the total number of units to which
C                     each error message is to be written.
C                 = 6 Refers to the 2nd unit for error messages
C                 = 7 Refers to the 3rd unit for error messages
C                 = 8 Refers to the 4th unit for error messages
C                 = 9 Refers to the 5th unit for error messages
C        IVALUE - The value to be set for the IWHICH-th parameter,
C                 if ISET is .TRUE. .
C        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
C                 given the value, IVALUE.  If ISET=.FALSE., the
C                 IWHICH-th parameter will be unchanged, and IVALUE
C                 is a dummy parameter.
C      --Output--
C        The (old) value of the IWHICH-th parameter will be returned
C        in the function value, J4SAVE.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C    Adapted from Bell Laboratories PORT Library Error Handler
C     Latest revision ---  23 MAY 1979
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
C***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION D9LGMC(X)
C***BEGIN PROLOGUE  D9LGMC
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  COMPLETE GAMMA FUNCTION,CORRECTION FACTOR,
C             DOUBLE PRECISION,GAMMA FUNCTION,LOGARITHM,
C             SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the  d.p. log Gamma correction factor for
C            X .GE. 10. so that DLOG(DGAMMA(X)) = DLOG(DSQRT(2*PI)) +
C            (X-5.)*DLOG(X) - X + D9LGMC(X)
C***DESCRIPTION
C
C Compute the log gamma correction factor for X .GE. 10. so that
C DLOG (DGAMMA(X)) = DLOG(DSQRT(2*PI)) + (X-.5)*DLOG(X) - X + D9lGMC(X)
C
C Series for ALGM       on the interval  0.          to  1.00000E-02
C                                        with weighted error   1.28E-31
C                                         log weighted error  30.89
C                               significant figures required  29.81
C                                    decimal places required  31.48
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DCSEVL,INITDS,XERROR
C***END PROLOGUE  D9LGMC
      DOUBLE PRECISION X, ALGMCS(15), XBIG, XMAX, DCSEVL, D1MACH
      DATA ALGMCS(  1) / +.1666389480 4518632472 0572965082 2 D+0      /
      DATA ALGMCS(  2) / -.1384948176 0675638407 3298605913 5 D-4      /
      DATA ALGMCS(  3) / +.9810825646 9247294261 5717154748 7 D-8      /
      DATA ALGMCS(  4) / -.1809129475 5724941942 6330626671 9 D-10     /
      DATA ALGMCS(  5) / +.6221098041 8926052271 2601554341 6 D-13     /
      DATA ALGMCS(  6) / -.3399615005 4177219443 0333059966 6 D-15     /
      DATA ALGMCS(  7) / +.2683181998 4826987489 5753884666 6 D-17     /
      DATA ALGMCS(  8) / -.2868042435 3346432841 4462239999 9 D-19     /
      DATA ALGMCS(  9) / +.3962837061 0464348036 7930666666 6 D-21     /
      DATA ALGMCS( 10) / -.6831888753 9857668701 1199999999 9 D-23     /
      DATA ALGMCS( 11) / +.1429227355 9424981475 7333333333 3 D-24     /
      DATA ALGMCS( 12) / -.3547598158 1010705471 9999999999 9 D-26     /
      DATA ALGMCS( 13) / +.1025680058 0104709120 0000000000 0 D-27     /
      DATA ALGMCS( 14) / -.3401102254 3167487999 9999999999 9 D-29     /
      DATA ALGMCS( 15) / +.1276642195 6300629333 3333333333 3 D-30     /
      DATA NALGM, XBIG, XMAX / 0, 2*0.D0 /
C***FIRST EXECUTABLE STATEMENT  D9LGMC
      IF (NALGM.NE.0) GO TO 10
      NALGM = INITDS (ALGMCS, 15, SNGL(D1MACH(3)) )
      XBIG = 1.0D0/DSQRT(D1MACH(3))
      XMAX = DEXP (DMIN1(DLOG(D1MACH(2)/12.D0), -DLOG(12.D0*D1MACH(1))))
C
 10   IF (X.LT.10.D0) CALL XERROR ( 'D9LGMC  X MUST BE GE 10', 23, 1, 2)
      IF (X.GE.XMAX) GO TO 20
C
      D9LGMC = 1.D0/(12.D0*X)
      IF (X.LT.XBIG) D9LGMC = DCSEVL (2.0D0*(10.D0/X)**2-1.D0, ALGMCS,
     1  NALGM) / X
      RETURN
C
 20   D9LGMC = 0.D0
      CALL XERROR ( 'D9LGMC  X SO BIG D9LGMC UNDERFLOWS', 34, 2, 1)
      RETURN
C
      END
      FUNCTION INITDS(DOS,NOS,ETA)
C***BEGIN PROLOGUE  INITDS
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C3A2
C***KEYWORDS  CHEBYSHEV,DOUBLE PRECISION,INITIALIZE,
C             ORTHOGONAL POLYNOMIAL,SERIES,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Initializes the d.p. properly normalized orthogonal
C            polynomial series to determine the number of terms needed
C            for specific accuracy.
C***DESCRIPTION
C
C Initialize the double precision orthogonal series DOS so that INITDS
C is the number of terms needed to insure the error is no larger than
C ETA.  Ordinarily ETA will be chosen to be one-tenth machine precision
C
C             Input Arguments --
C DOS    dble prec array of NOS coefficients in an orthogonal series.
C NOS    number of coefficients in DOS.
C ETA    requested accuracy of series.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  INITDS
C
      DOUBLE PRECISION DOS(NOS)
C***FIRST EXECUTABLE STATEMENT  INITDS
      IF (NOS.LT.1) CALL XERROR ( 'INITDS  NUMBER OF COEFFICIENTS LT 1',
     1 35, 2, 2)
C
      ERR = 0.
      DO 10 II=1,NOS
        I = NOS + 1 - II
        ERR = ERR + ABS(SNGL(DOS(I)))
        IF (ERR.GT.ETA) GO TO 20
 10   CONTINUE
C
 20   IF (I.EQ.NOS) CALL XERROR ( 'INITDS  ETA MAY BE TOO SMALL', 28,
     1  1, 2)
      INITDS = I
C
      RETURN
      END
