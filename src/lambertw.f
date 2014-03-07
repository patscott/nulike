!=======================================================================
! Calculates W(x) to double precision, where W is the Lambert W
! function defined as the solution w = W(x) to the equation:
!     x = w e^w
! Returns the principal branch value (the larger of two solutions
! over -1/e < x < 0; there is only one solution for x > 0).  The
! W function is undefined for x < -1/e.
! 
! Valid input:  -1/e < x < \infty
! 
! 
!    Created by Chris Savage (savage@fysik.su.se)
!    2011/06/22
! 
!=======================================================================
! 
      REAL*8 FUNCTION lambertw(x)
      IMPLICIT NONE
      REAL*8 x
      REAL*8 epsk,zk,qk,wk,wk1,p,num,den,a,b,ia
      REAL*8 e,einv
      PARAMETER(e=2.7182818284590452d0,einv=0.36787944117144232d0)
      
      ! No solution for x < 1/e
      IF (x .LT. -einv) THEN
        !STOP 'Error in lambertw: argument must be larger than -EXP(-1)'
        lambertw = -HUGE(lambertw)
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
      IF (x .LT. -0.32358170806015724d0) THEN
        ! branch point expansion
        p = SQRT(2*(1+e*x))
        wk = -1 + p*(1 + p*(-1d0/3 + p*(11d0/72 + p*(-43d0/540
     &            + p*(769d0/17280 + p*(-221d0/8505
     &            + p*(680863d0/43545600 + p*(-1963d0/204120
     &            + p*(226287557d0/37623398400d0) ))))))))
      ELSE IF (x .LT. 0.14546954290661823d0) THEN
        ! rational fit
        num = x * (1 + x *
     &            (5.931375839364438d0 + x *
     &            (11.392205505329132d0+ x *
     &            (7.338883399111118d0 + x * 0.6534490169919599d0) )))
        den = 1 + x *
     &           (6.931373689597704d0 + x *
     &           (16.82349461388016d0 + x *
     &           (16.43072324143226d0 + x * 5.115235195211697d0) ))
        wk = num / den
      ELSE IF (x .LT. 8.706658967856612d0) THEN
        ! rational fit
        num = x * (1 + x *
     &            (2.4450530707265568d0 + x *
     &            (1.3436642259582265d0 + x *
     &            (0.14844005539759195d0+ x * 8.047501729129999d-4) )))
        den = 1 + x *
     &           (3.4447089864860025d0 + x *
     &           (3.2924898573719523d0 + x *
     &           (0.9164600188031222d0 + x * 0.05306864044833221d0) ))
        wk = num / den
      ELSE
        ! asymptotic expansion
        a = LOG(x)
        b = LOG(a)
        ia = 1/a
        wk = a - b + b * ia *
     &              (1 + ia *
     &              (0.5d0*(-2 + b) + ia *
     &              (1/6d0*(6 + b*(-9 + b*2)) + ia *
     &              (1/12d0*(-12 + b*(36 + b*(-22 + b*3))) + ia *
     &               1/60d0*(60 + b*(-300 + b*(350 + b*(-125 + b*12))))
     &              ))))
      END IF
      
      ! Special cases:
      ! For x equal to 0 or -1/e, the Fritsch iteration does
      ! not work as some of the terms go to infinity.  However,
      ! for x sufficiently near 0 or -1/e, the above first
      ! approximation is already nearly double precision.
      IF ((ABS(x) .LT. 1d-7) .OR. (x .LT. -einv+1d-6) ) THEN
        lambertw = wk
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
      
      lambertw = wk
      
      END FUNCTION


