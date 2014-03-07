      SUBROUTINE YPCOEF (SIGMA,DX, D,SD)
      DOUBLE PRECISION SIGMA, DX, D, SD
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
C   This subroutine computes the coefficients of the deriva-
C tives in the symmetric diagonally dominant tridiagonal
C system associated with the C-2 derivative estimation pro-
C cedure for a Hermite interpolatory tension spline.
C
C On input:
C
C       SIGMA = Nonnegative tension factor associated with
C               an interval.
C
C       DX = Positive interval width.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       D = Component of the diagonal term associated with
C           the interval.  D = SIG*(SIG*COSHM(SIG) -
C           SINHM(SIG))/(DX*E), where SIG = SIGMA and E =
C           SIG*SINH(SIG) - 2*COSHM(SIG).
C
C       SD = Subdiagonal (superdiagonal) term.  SD = SIG*
C            SINHM(SIG)/E.
C
C Module required by YPCOEF:  SNHCSH
C
C Intrinsic function called by YPCOEF:  EXP
C
C***********************************************************
C
      DOUBLE PRECISION COSHM, COSHMM, E, EMS, SCM, SIG,
     .                 SINHM, SSINH, SSM
C
      SIG = SIGMA
      IF (SIG .LT. 1.D-9) THEN
C
C SIG = 0:  cubic interpolant.
C
        D = 4.D0/DX
        SD = 2.D0/DX
      ELSEIF (SIG .LE. .5D0) THEN
C
C 0 .LT. SIG .LE. .5:  use approximations designed to avoid
C                      cancellation error in the hyperbolic
C                      functions when SIGMA is small.
C
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        E = (SIG*SINHM - COSHMM - COSHMM)*DX
        D = SIG*(SIG*COSHM-SINHM)/E
        SD = SIG*SINHM/E
      ELSE
C
C SIG > .5:  scale SINHM and COSHM by 2*EXP(-SIG) in order
C            to avoid overflow when SIGMA is large.
C
        EMS = EXP(-SIG)
        SSINH = 1.D0 - EMS*EMS
        SSM = SSINH - 2.D0*SIG*EMS
        SCM = (1.D0-EMS)*(1.D0-EMS)
        E = (SIG*SSINH - SCM - SCM)*DX
        D = SIG*(SIG*SCM-SSM)/E
        SD = SIG*SSM/E
      ENDIF
      RETURN
      END

