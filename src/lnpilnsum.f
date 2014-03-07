!=======================================================================
! LN(Poisson Integral) for Log-Normal systematic uncertainty (SUMs)
! 
! This routine determines the sums over the probability terms, where
! the probability is a Poissonian marginalized over a systematic
! uncertainty in the average expected number of signal events (e.g.
! from a systematic uncertainty in the detector sensitivity).  The
! systematic uncertainty is taken to be a log-normal distributed
! fractional variation in the signal contribution to the Poisson
! distribution average.  The probability (marginalized likelihood)
! is given by:
!     P(n) = I(n,thetab,thetas,sigma)
!          = 1 / sqrt{2 \pi \sigma^2}
!            \int_0^\infty d\epsilon / epsilon
!                          e^{-(ln(\epsilon))^2 / 2\sigma^2}
!                          (\theta_b + \epsilon\theta_s)^n
!                          e^{-(\theta_b + \epsilon\theta_s)} / n!
! 
! Inputs:
!     n         Number of observed events
!     thetab    Average number of background events
!     thetas    Average number of signal events
!     sigma     Approximate s.d. of epsilon, the fractional variation
!               in thetas, i.e. thetas -> epsilon thetas (epsilon has
!               an approximate average of 1).  This is actually the
!               sigma parameter of the log-normal distribution, which
!               is not technically the s.d. of the distribution,
!               although it approaches the s.d. as sigma becomes much
!               smaller than the average.
! Outputs:
!     lnpiln    Logarithm of P(n)
!     lnlesum   Logarithm of \Sum_{k \le n} P(k)
!               i.e. the logarithm of the CDF
!     lngesum   Logarithm of \Sum_{k \ge n} P(k)
! NOTE: Both sums include P(n), so they obey this relation:
!       EXP(lnlesum) + EXP(lngesum) = 1 + EXP(lnpiln)
! 
! For further details, see:
!   * Conrad, Botner, Hallgren & Perez de los Heros, Phys. Rev. D 67,
!     012002 (2003) [arXiv:hep-ex/0202013].
!   * Scott et al., JCAP 1001, 031 (2010) [arXiv:0909.3300
!     [astro-ph.CO]].
! 
! 
!    Created by Chris Savage (savage@fysik.su.se)
!    2011/11/02
! 
!=======================================================================
! 
      SUBROUTINE nulike_lnpilnsum(n,thetab,thetas,sigma,
     &                       lnpiln,lnlesum,lngesum)
      IMPLICIT NONE
      INTEGER n
      REAL*8 thetab,thetas,sigma,lnpiln,lnlesum,lngesum
      REAL*8 LN_PRECISION
      ! Set the precision scale to ~ 1d-15.
      ! The true precision of the sum is as low as ~ Nterm*1d-15,
      ! where Nterm is the number of terms in the sum (typically < n).
      PARAMETER(LN_PRECISION=-35d0)   ! Precision of ~ 1d-15
      
      ! Checks --------------------------------------------
      IF ((n .LT. 0) .OR. (thetab .LT. 0d0)
     &    .OR. (thetas .LT. 0d0) .OR. (sigma .LT. 0d0)) THEN
        WRITE(*,*) 'ERROR: nulike_lnpilnsum called with negative argument'
        STOP
      END IF
      
      ! Special cases -------------------------------------
      ! We first handle some special cases where the probability and
      ! sums simplify:
      !    thetas = 0
      !    sigma  = 0
      
      ! Special case: thetas = 0
      ! The integration drops out and we have a strictly Poisson
      ! probability with average thetab
      IF (thetas .EQ. 0d0) THEN
        CALL poissonsums(n,thetab,lnpiln,lnlesum,lngesum)
        RETURN
      END IF
      
      ! Special case: sigma = 0
      ! The integration drops out and we have a strictly Poisson
      ! probability with average thetab+thetas
      IF (sigma .EQ. 0d0) THEN
        CALL poissonsums(n,thetab+thetas,lnpiln,lnlesum,lngesum)
        RETURN
      END IF
      
      ! General case ---------------------------------------
      CALL numericalsums(n,thetab,thetas,sigma,lnpiln,lnlesum,lngesum)
      
      
      CONTAINS
      
      ! ----------------------------------------------------------------
      ! Subroutine to calculate the sum of Poisson probabilities:
      !    P(n|theta) = theta^n EXP(-theta) / n!
      ! 
      SUBROUTINE poissonsums(n,theta,lnpmf,lnlesum,lngesum)
      IMPLICIT NONE
      INTEGER n
      REAL*8 theta,lnpmf,lnlesum,lngesum
      INTEGER k,kmax
      REAL*8 lnp,lngamma
      
      ! Special case: theta = 0
      IF (theta .EQ. 0d0) THEN
        IF (n .EQ. 0) THEN
          lnpmf   = 0d0
          lnlesum = 0d0
          lngesum = 0d0
        ELSE
          lnpmf   = -HUGE(1d0)
          lnlesum = 0d0
          lngesum = -HUGE(1d0)
        END IF
        RETURN
      END IF
      
      ! Special case: n = 0
      IF (n .EQ. 0) THEN
        lnpmf   = -theta
        lnlesum = lnpmf
        lngesum = 0d0
        RETURN
      END IF
      
      ! General case
      lnpmf = n*LOG(theta) - theta - lngamma(n+1d0)
      
      ! The distribution peaks around n = theta; to avoid a loss of
      ! precision, we will do the lower sum for n < theta and the
      ! upper sum for n > theta.
      IF (n .LE. theta) THEN
        k   = n
        lnp = lnpmf
        lnlesum = lnpmf
        DO WHILE (k .GT. 0)
          k = k-1
          lnp = lnp + LOG((k+1)/theta)
          lnlesum = lnsum(lnlesum,lnp)
          IF (lnp - lnlesum .LE. LN_PRECISION) EXIT
        END DO
        lngesum = LOG(1d0 - (EXP(lnlesum) - EXP(lnpmf)))
      ELSE
        k   = n
        lnp = lnpmf
        lngesum = lnpmf
        ! Determine a conservative upper limit on sum terms.
        ! Will hopefully hit precision condition first, but include
        ! this upper limit in case we didn't think of some pathological
        ! cases.
        IF ((n-theta)**2 .GT. theta) THEN
          kmax = n + NINT(-LN_PRECISION * (n/(n-theta))) + 10
        ELSE
          kmax = n + NINT(-LN_PRECISION * (1d0+SQRT(theta))) + 10
        END IF
        DO WHILE (k .LT. kmax)
          k = k+1
          lnp = lnp + LOG(theta/k)
          lngesum = lnsum(lngesum,lnp)
          IF (lnp - lngesum .LE. LN_PRECISION) EXIT
        END DO
        lnlesum = LOG(1d0 - (EXP(lngesum) - EXP(lnpmf)))
      END IF
      
      END SUBROUTINE poissonsums
      
      
      ! ----------------------------------------------------------------
      ! Subroutine to calculate the sum of probabilities using the
      ! numerical integration routines to calculate those
      ! probabilities.  This routine should not be used for thetas=0
      ! or sigma=0.
      ! 
      SUBROUTINE numericalsums(n,thetab,thetas,sigma,
     &                         lnpmf,lnlesum,lngesum)
      IMPLICIT NONE
      INTEGER n
      REAL*8 thetab,thetas,sigma,lnpmf,lnlesum,lngesum
      INTEGER k
      REAL*8 lnp,nulike_lnpiln
      LOGICAL is_increasing
      
      ! Special case: n = 0
      IF (n .EQ. 0) THEN
        lnpmf   = nulike_lnpiln(n,thetab,thetas,sigma)
        lnlesum = lnpmf
        lngesum = 0d0
        RETURN
      END IF
      
      ! General case
      lnpmf = nulike_lnpiln(n,thetab,thetas,sigma)
      
      ! Here, we check if the probability is increasing w.r.t. n
      lnp = nulike_lnpiln(n+1,thetab,thetas,sigma)
      is_increasing = lnp .GT. lnpmf
      
      ! To avoid a loss of precision, we will do the lower sum if
      ! we are in the lower tail of lnlin and the upper sum if
      ! we are in the upper tail of lnlin.
      IF (is_increasing) THEN
        k   = n
        lnp = lnpmf
        lnlesum = lnpmf
        DO WHILE (k .GT. 0)
          k = k-1
          lnp = nulike_lnpiln(k,thetab,thetas,sigma)
          lnlesum = lnsum(lnlesum,lnp)
          IF (lnp - lnlesum .LE. LN_PRECISION) EXIT
        END DO
        lngesum = LOG(1d0 - (EXP(lnlesum) - EXP(lnpmf)))
      ELSE
        k   = n
        lnp = lnpmf
        lngesum = lnpmf
        DO
          k = k+1
          lnp = nulike_lnpiln(k,thetab,thetas,sigma)
          lngesum = lnsum(lngesum,lnp)
          IF (lnp - lngesum .LE. LN_PRECISION) EXIT
        END DO
        lnlesum = LOG(1d0 - (EXP(lngesum) - EXP(lnpmf)))
      END IF
      
      END SUBROUTINE numericalsums
      
      
      ! ----------------------------------------------------------------
      ! Function to calculate the quantity ln(a+b) given ln(a) and
      ! ln(b).
      ! 
      REAL*8 FUNCTION lnsum(lna,lnb)
      IMPLICIT NONE
      REAL*8 lna,lnb,lnx,r
      
      IF (lna .GE. lnb) THEN
        lnx = lna
        r   = EXP(lnb-lna)
      ELSE
        lnx = lnb
        r   = EXP(lna-lnb)
      END IF
      
      ! Below, we use an expansion of ln(1+r) if r is small.
      ! The sum is terminated at approximately double precision.
      IF (r .EQ. 0d0) THEN
        lnsum = lnx
      ELSE IF (r .GT. 0.01d0) THEN
        lnsum = lnx + LOG(1d0 + r)
      ELSE IF (r .GT. 0.001d0) THEN
        lnsum = lnx + r*(1d0
     &                   -r*((1d0/2d0)
     &                       -r*((1d0/3d0)
     &                           -r*((1d0/4d0)
     &                               -r*((1d0/5d0)
     &                                   -r*((1d0/6d0)
     &                                       -r*((1d0/7d0)
     &                                           -r*((1d0/8d0)
     &                                               -r*(1d0/9d0)
     &                  ))))))))
      ELSE IF (r .GT. 0.00001d0) THEN
        lnsum = lnx + r*(1d0
     &                   -r*((1d0/2d0)
     &                       -r*((1d0/3d0)
     &                           -r*((1d0/4d0)
     &                               -r*((1d0/5d0)
     &                                   -r*(1d0/6d0)
     &                  )))))
      ELSE
        lnsum = lnx + r*(1d0
     &                   -r*((1d0/2d0)
     &                       -r*((1d0/3d0)
     &                           -r*(1d0/4d0)
     &                  )))
      END IF
      
      END FUNCTION lnsum
      
      END SUBROUTINE


