!=======================================================================
! Calculates W(x) to double precision, where W is the Lambert W
! function defined as the solution w = W(x) to the equation:
!     x = w e^w
! Returns the secondary branch value (the smaller of two solutions
! over -1/e < x < 0).  This branch is defined only for
! -1/e < x < 0.
! 
! Valid input:  -1/e < x < 0
! 
! 
!    Created by Chris Savage (savage@fysik.su.se)
!    2011/06/22
! 
!=======================================================================
! 
      REAL*8 FUNCTION lambertw2(x)
      IMPLICIT NONE
      REAL*8 x
      INTEGER k
      REAL*8 epsk,zk,qk,wk,wk1,p,num,den,a
      REAL*8 e,einv
      PARAMETER(e=2.7182818284590452d0,einv=0.36787944117144232d0)
      
      ! No solution for x < 1/e
      IF (x .LT. -einv) THEN
        !STOP 'Error in lambertw2: argument must be larger than -EXP(-1)'
        lambertw2 = HUGE(lambertw2)
        RETURN
      END IF
      
      ! No second branch for x > 0
      IF (x .GT. 0) THEN
        !STOP 'Error in lambertw2: argument must be smaller than 0'
        lambertw2 = HUGE(lambertw2)
        RETURN
      END IF
      
      ! Limit as x ->0 is -\infty
      IF (x .EQ. 0) THEN
        lambertw2 = -HUGE(lambertw2)
        RETURN
      END IF
      
      ! Could use Newton's or Halley's iteration method to find the
      ! solution, but the Lambert function has a faster method,
      ! Fritsch's iteration:
      !    W_{k+1} = W_k * (1 + eps_k)
      ! with
      !    eps_k   = z_k/(1 + W_k) * (q_k - z_k)/(q_k - 2*z_k)
      !    z_k     = ln(x/W_k) - W_k
      !    q_k     = 2 (1 + W_k) (1 + W_k + (2/3)*z_k)
      ! If eps_k is the error in W_k, then the error in W_{k+1} is of
      ! order (eps_k)^4, a faster convergent that Halley's method.
      ! For a first guess accurate to order O(10^-4), double precision
      ! can be achieved in only a single iteration.  We use estimates
      ! for the initial guess as determined by Veberic.
      ! 
      ! For further information, see:
      !   D. Veberic, arXiv:1003.1628.
      !   F.N. Fritsch, R.E. Shafer and W.P. Crowley,
      !       Commun. ACM 16, 123 (1973).
      
      ! Initial estimate by Veberic
      IF (x .LT. -0.30298541769d0) THEN
        ! branch point expansion
        p = -SQRT(2*(1+e*x))
        wk = -1 + p*(1 + p*(-1d0/3 + p*(11d0/72 + p*(-43d0/540
     &            + p*(769d0/17280 + p*(-221d0/8505
     &            + p*(680863d0/43545600 + p*(-1963d0/204120
     &            + p*(226287557d0/37623398400d0) ))))))))
      ELSE IF (x .LT. -0.051012917658221676d0) THEN
        ! rational fit
        
        num = -7.814176723907436d0 + x *
     &          (253.88810188892484d0 + x * 657.9493176902304d0)
        den = 1 + x *
     &           (-60.43958713690808d0+ x *
     &           (99.98567083107612d0 + x *
     &           (682.6073999909428d0 + x *
     &           (962.1784396969866d0 + x * 1477.9341280760887d0) )))
        wk = num / den
      ELSE
        ! continued logarithm
        a  = LOG(-x)
        wk = a
        DO k=1,9
          wk = a - LOG(-wk)
        END DO
      END IF
      
      ! Special cases:
      ! For x equal to -1/e, the Fritsch iteration does not
      ! not work as some of the terms go to infinity.  However,
      ! for x sufficiently near -1/e, the above first
      ! approximation is already nearly double precision.
      IF (x .LT. -einv+1d-6) THEN
        lambertw2 = wk
        RETURN
      END IF
      
      ! Now apply Fritsch iteration
      wk1  = wk + 1
      zk   = LOG(x/wk) - wk
      qk   = 2 * wk1 * (wk1 + (2d0/3)*zk)
      epsk = zk * (qk - zk) / (wk1 * (qk - 2*zk))
      wk   = wk * (1 + epsk)
      ! In most cases, no further iterations will be necessary
      DO WHILE (ABS(epsk) .GT. 1d-5)
        wk1  = wk + 1
        zk   = LOG(x/wk) - wk
        qk   = 2 * wk1 * (wk1 + (2d0/3)*zk)
        epsk = zk * (qk - zk) / (wk1 * (qk - 2*zk))
        wk   = wk * (1 + epsk)
      END DO
      
      lambertw2 = wk
      
      END FUNCTION


