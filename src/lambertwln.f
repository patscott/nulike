!=======================================================================
! Calculates W(x=e^ln(x)) to double precision, where W is the Lambert W
! function defined as the solution w = W(x) to the equation:
!     x = w e^w
! Returns the principal branch value (the larger of two solutions
! over -1/e < x < 0; there is only one solution for x > 0).  The
! W function is undefined for x < -1/e.
! 
! This version takes ln(x) as the input to allow for cases where x
! is large.
! 
! Valid input:  -\infty < lnx < \infty
! 
! 
!    Created by Chris Savage (savage@fysik.su.se)
!    2011/06/22
! 
!=======================================================================
! 
      REAL*8 FUNCTION lambertwln(lnx)
      IMPLICIT NONE
      REAL*8 lnx
      REAL*8 epsk,zk,qk,wk,wk1,a,b,ia,lambertw
      REAL*8 e,einv
      PARAMETER(e=2.7182818284590452d0,einv=0.36787944117144232d0)
      
      ! Here, we only calculate W(x) for very large x.  If x is a
      ! not very large, we use the lambertw routine.
      IF (lnx .LT. 300d0) THEN
        lambertwln = lambertw(EXP(lnx))
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
      ! asymptotic expansion
      a = lnx
      b = LOG(a)
      ia = 1/a
      wk = a - b + b * ia *
     &            (1 + ia *
     &            (0.5d0*(-2 + b) + ia *
     &            (1/6d0*(6 + b*(-9 + b*2)) + ia *
     &            (1/12d0*(-12 + b*(36 + b*(-22 + b*3))) + ia *
     &             1/60d0*(60 + b*(-300 + b*(350 + b*(-125 + b*12))))
     &            ))))
      
      ! Now apply Fritsch iteration
      wk1  = wk + 1
      zk   = lnx - LOG(wk) - wk
      qk   = 2 * wk1 * (wk1 + (2d0/3)*zk)
      epsk = zk * (qk - zk) / (wk1 * (qk - 2*zk))
      wk   = wk * (1 + epsk)
      ! In most cases, no further iterations will be necessary
      DO WHILE (ABS(epsk) .GT. 1d-5)
        wk1  = wk + 1
        zk   = lnx - LOG(wk) - wk
        qk   = 2 * wk1 * (wk1 + (2d0/3)*zk)
        epsk = zk * (qk - zk) / (wk1 * (qk - 2*zk))
        wk   = wk * (1 + epsk)
      END DO
      
      lambertwln = wk
      
      END FUNCTION


