!=======================================================================
! LN(Poisson Integral) for Normal systematic uncertainty
! 
! The logarithm of a Poissonian likelihood marginalized over a
! systematic uncertainty in the average expected number of signal
! events (e.g. from a systematic uncertainty in the detector
! sensitivity).  The systematic uncertainty is taken to be a
! normally distributed fractional variation in the signal
! contribution to the Poisson distribution average.  The
! marginalized likelihood is given by:
!     I(n,thetab,thetas,sigma)
!         = A / sqrt{2 \pi \sigma^2}
!           \int_0^\infty d\epsilon e^{-(\epsilon-1)^2 / 2\sigma^2}
!                         (\theta_b + \epsilon\theta_s)^n
!                         e^{-(\theta_b + \epsilon\theta_s)} / n!
! with
!     A^{-1} = 1/2 Erfc(-1/sqrt{2\sigma^2})
! a normalization factor accounting for the the fact that the
! Gaussian is only integrated down to zero instead of -\infty.
! 
! Inputs:
!     n         Number of observed events
!     thetab    Average number of background events
!     thetas    Average number of signal events
!     sigma     s.d. of epsilon, the fractional variation in thetas,
!               i.e. thetas -> epsilon thetas (epsilon has average 1)
! 
! For further details, see:
!   * Conrad, Botner, Hallgren & Perez de los Heros, Phys. Rev. D 67,
!     012002 (2003) [arXiv:hep-ex/0202013].
!   * Scott et al., JCAP 1001, 031 (2010) [arXiv:0909.3300
!     [astro-ph.CO]].
! 
! This routine uses extensive optimizations to perform the calculation
! quickly and is fairly robust in handling even extreme cases for the
! parameters.
! 
! 
!    Created by Chris Savage (savage@fysik.su.se)
!    2011/06/02
! 
!=======================================================================
! 
      REAL*8 FUNCTION nulike_lnpin(n,thetab,thetas,sigma)
      IMPLICIT NONE
      INTEGER n
      REAL*8 thetab,thetas,sigma
      ! Largest n for using recurrence relation (thetab=0 case).
      ! Above this, we will just use a numerical integration.
      INTEGER RECURRENCE_NMAX
      PARAMETER(RECURRENCE_NMAX=1024)
      
      ! Checks --------------------------------------------
      IF ((n .LT. 0) .OR. (thetab .LT. 0d0)
     &    .OR. (thetas .LT. 0d0) .OR. (sigma .LT. 0d0)) THEN
        WRITE(*,*) 'ERROR: nulike_lnpin called with negative argument'
        STOP
      END IF
      
      ! Special cases -------------------------------------
      ! Before attempting to do a numerical integration, we first
      ! handle some special cases:
      !    thetas = 0
      !    sigma  = 0
      
      ! Special case: thetas = 0
      ! The integration drops out and we have a strictly Poisson
      ! probability with average thetab
      IF (thetas .EQ. 0d0) THEN
        nulike_lnpin = lnpoisson(n,thetab)
        RETURN
      END IF
      
      ! Special case: sigma = 0
      ! The integration drops out and we have a strictly Poisson
      ! probability with average thetab+thetas
      IF (sigma .EQ. 0d0) THEN
        nulike_lnpin = lnpoisson(n,thetab+thetas)
        RETURN
      END IF
      
      ! General case ---------------------------------------
      ! In most cases, we will do a numerical integration.
      ! For the special case thetab=0, a recurrence relation can be
      ! used.  Use of the recurrence relation is fast as long as n is
      ! not extremely large and results can be cached to speed up
      ! subsequent calls with the same thetas and sigma and a larger n
      ! (the recurrence relation is simply continued using the previous
      ! result instead of starting from the beginning).  However, the
      ! recurrence is numerically unstable for theta*sigma^2 > 1.
      ! In that case or for larger n, we will perform a numerical
      ! integration.  The two methods are good to approximately
      ! double precision in the cases they are used.
      IF ((thetab .EQ. 0d0) .AND. (n .LE. RECURRENCE_NMAX)
     &                      .AND. (thetas*sigma**2 .LE. 1d0)) THEN
        nulike_lnpin = recurrence(n,thetas,sigma)
      ELSE
        nulike_lnpin = nintegrate(n,thetab,thetas,sigma)
      END IF
      
      
      CONTAINS
      
      ! ----------------------------------------------------------------
      ! Subroutine to calculate the logarithm of a Poisson probability:
      !    P(n|theta) = theta^n EXP(-theta) / n!
      ! 
      REAL*8 FUNCTION lnpoisson(n,theta)
      IMPLICIT NONE
      INTEGER n
      REAL*8 theta
      REAL*8 lngamma
      
      ! Special case: theta = 0
      !    P(0|0) = 1
      !    P(n|0) = 0   [n > 0]
      IF (theta .EQ. 0d0) THEN
        IF (n .EQ. 0) THEN
          lnpoisson = 0d0
        ELSE
          lnpoisson = -HUGE(1d0)
        END IF
        RETURN
      END IF
      
      ! Special case: n = 0
      !    P(0|theta) = EXP(-theta)
      IF (n .EQ. 0) THEN
        lnpoisson = -theta
        RETURN
      END IF
      
      ! General case
      lnpoisson = n*LOG(theta) - theta - lngamma(n+1d0)
      
      END FUNCTION lnpoisson
      
      
      ! ----------------------------------------------------------------
      ! Function to calculate the integral using a recurrence relation
      ! (applies only for thetab=0).
      ! 
      ! The integral I(n) = I(n,thetab=0,thetas=theta,sigma) satisfies
      ! a recurrence relation:
      !   I(n) = (1/n) [B1*I(n-1) + B2*I(n-2)]
      ! where:
      !   B1 = theta*epsilon0
      !   B2 = theta^2 * sigma^2
      !   epsilon0 = 1 - theta*sigma^2
      !   u0 = -1 / sqrt{2 sigma^2}
      !   x0 = -epsilon0 / sqrt{2 sigma^2}
      ! The recurrence relation can be started with:
      !   I(-1,theta,sigma)
      !       = A e^{-1/2sigma^2} / (\sqrt{2\pi sigma^2} theta)
      !       = [\sqrt{\pi/2} theta*sigma e^{u0^2} erfc(u0)]^{-1}
      !   I(0,theta,sigma)
      !       = A e^{-1/2sigma^2} e^{x0^2} erfc(x0) / 2
      ! Note I(-1,theta,sigma) is defined as above only to start the
      ! recurrence; it does not correspond to any actual integration.
      ! 
      ! NOTE: The recurrence is numerically unstable for epsilon0 < 0,
      !       so this routine should not be used in that case (use e.g.
      !       a numerical integration instead).
      ! 
      ! Below, we use ratios:
      !   r_n = I(n)/I(n-1)
      !     => r_n = (1/n) [B1 + B2/r_{n-1}]
      !        r_0 = \sqrt{\pi/2} theta*sigma e^{x0^2} erfc(x0)
      ! Then:
      !   I(n,theta) = I(-1,theta) \Pi_{k=0}^{n} r_k
      
      ! Because this routine may be called repeatedly with the same
      ! theta and sigma, we will cache the previous result so we can
      ! resume the iteration from the previous spot.
      ! 
      ! The result is typically accurate to near full double precision.
      ! 
      REAL*8 FUNCTION recurrence(n,theta,sigma)
      IMPLICIT NONE
      INTEGER n
      REAL*8 theta,sigma
      INTEGER k
      REAL*8 theta_prev,sigma_prev
      REAL*8 B1,B2,lnIk,rk
      DATA theta_prev / -1d0 /
      DATA sigma_prev / -1d0 /
      SAVE theta_prev,sigma_prev,B1,B2,k,lnIk,rk
      
      ! NOTE: this routine should not be called for any of these cases:
      !   *) theta = 0
      !   *) sigma = 0
      !   *) theta*sigma^2 > 1
      ! These cases should be checked for before calling this function.
      ! The last case is fine for n=0 as this is calculated explicitly.
      
      ! Cannot use cache: must initialize the recursive algorithm
      IF ((theta .NE. theta_prev) .OR. (sigma .NE. sigma_prev)
     &     .OR. (n .LT. k) ) THEN
        theta_prev = theta
        sigma_prev = sigma
        B1   = theta*(1-theta*sigma**2)
        B2   = (theta*sigma)**2
        k = 0
        ! Calculate I(0) and r_0
        CALL recurrence_init(n,theta,sigma,lnIk,rk)
      END IF
      
      ! Apply recurrence relation
      DO WHILE (k < n)
        k  = k+1
        rk = (B1 + B2/rk)/k
        lnIk = lnIk + LOG(rk)
        !WRITE(*,*) 'DEBUG:',k,rk,lnIk
      END DO
      
      recurrence = lnIk
      
      END FUNCTION recurrence
      
      
      ! ----------------------------------------------------------------
      ! Determines the quantities I(n=0) and r_0 = I(0)/I(-1), which
      ! are used to start the recurrence algorithm.  See the
      ! description for the function "recurrence" above for details.
      ! 
      ! The results are typically accurate to near full double
      ! precision.
      ! 
      SUBROUTINE recurrence_init(n,theta,sigma,lnI0,r0)
      IMPLICIT NONE
      INTEGER n
      REAL*8 theta,sigma,lnI0,r0
      REAL*8 eps0,u0,x0,y,z
      REAL*8 PI,SQRTPI,SQRT2,SQRTHALF
      PARAMETER(PI=3.1415926535897932d0)
      PARAMETER(SQRTPI=1.7724538509055160d0)
      PARAMETER(SQRT2=1.4142135623730950d0)
      PARAMETER(SQRTHALF=0.70710678118654752d0)

      ! Special case: theta = 0
      ! Note that r0 is irrelevant here as it is not needed in the
      ! recurrence relation
      IF (theta .EQ. 0d0) THEN
        lnI0 = 0d0
        r0   = 1d0
        RETURN
      END IF
      
      ! Special case: sigma = 0
      ! Note that r0 is irrelevant here as it is not needed in the
      ! recurrence relation
      IF (sigma .EQ. 0d0) THEN
        lnI0 = -theta
        r0   = 1d0
        RETURN
      END IF
      
      ! Some useful quantities
      eps0 = 1 - theta*sigma**2
      u0   = - 1 / (SQRT2*sigma)
      x0   = u0*eps0
      
      ! Calculate r_0
      ! Find z = sqrt(pi) e^{x0^2} erfc(x0), using asymptotic
      ! expansion of erfc for large x0:
      !   z -> (1/x0) [1 + \sum_{j=1} (-1)^j (2j-1)!!/(2x0^2)^j]
      ! For x0 > 25, double precision in six terms (j <= 6)
      IF (x0 .GT. 25d0) THEN
        y = 1 / x0**2
        z = (1-0.5d0*y
     &         *(1-1.5d0*y
     &             *(1-2.5d0*y
     &                 *(1-3.5d0*y
     &                     *(1-4.5d0*y
     &                         *(1-5.5d0*y))))))
     &      / x0
      ELSE
        z = SQRTPI*EXP(x0**2)*ERFC(x0)
      END IF
      r0 = SQRTHALF*theta*sigma*z
        
      ! Calculate I(0)
      IF (x0 .LT. -25d0) THEN
        ! Calculate I(0,theta) directly using:
        !   I(0,theta) = e^{x0^2-u0^2} erfc(x0) / erfc(u0)
        ! where x0^2 - u0^2 = -theta*(1-theta*sigma^2/2)
        ! Numerical issues if x0 > 25 (erfc(x0) -> 0)
        lnI0 =  LOG(ERFC(x0)/ERFC(u0))
     &          - theta*(1 - 0.5d0*theta*sigma**2)
      ELSE
        ! Calculate I(-1), then use I(0) = I(-1)*r_0
        ! Numerical issues if x0 < -25 (r_0 -> \infty)
        !z = SQRTPI*EXP(u0**2)*ERFC(u0)
        !lnI0 = -LOG(SQRTHALF*theta*sigma*z)
        ! I(-1,theta)
        lnI0 = -LOG(SQRTHALF*SQRTPI*theta*sigma*ERFC(u0)) - u0**2
        ! I(0,theta)
        lnI0 = lnI0 + LOG(r0)
      END IF
      
      n=n !No-warn operation that should get optimised away.

      END SUBROUTINE recurrence_init
      
      
      ! ----------------------------------------------------------------
      ! Function to calculate the integral numerically.
      ! 
      ! To obtain a better form for numerical integration, we choose
      ! a more suitable integration variable.  First, write the
      ! integral as:
      !   I(n,thetab,thetas,sigma)
      !       = A/\sqrt{2\pi sigma^2} \int_0^\infty d\epsilon e^{g(\epsilon)}
      !   g(eps) = n ln(thetab+eps*thetas) - ln(n!) - (thetab+eps*thetas)
      !            - (eps-1)^2/(2*sigma^2)
      ! Now, define eps0 and sigmau such that:
      !   g'(eps0)  = 0
      !   g''(eps0) = -1/sigmau^2
      ! Then we change the integration variable to:
      !   u = (eps - eps0)/sigmau
      ! Under this transformation:
      !   I  = A (sigmau/sigma) e^{g(eps0)} Iu
      !   Iu = 1/sqrt(2\pi) \int_u0^\infty du e^{f(u)}
      ! with:
      !   f(u) = g(eps) - g(eps0)
      !   u0   = -eps0/sigmau
      ! Note that:
      !   f(u) = -u^2/2 + O(u^3)
      ! so, to first order, the integration is now approximately a
      ! gaussian integration with unit variance.  The quantity Iu
      ! should be of order 1 unless u0 >> 1.
      ! 
      REAL*8 FUNCTION nintegrate(n,thetab,thetas,sigma)
      IMPLICIT NONE
      INTEGER n
      REAL*8 thetab,thetas,sigma
      REAL*8 sigma2,alpha,beta,eta,r,kappa,lambda,nkroot,w0,eps0,
     &       sigmau,sigmau2,u0,rho,Ainv
      REAL*8 lnC,g0,lnIu
      !REAL*8 lnIu1,lnIu2
      REAL*8 lngamma
      REAL*8 SQRT2
      PARAMETER(SQRT2=1.4142135623730950d0)
      
      ! NOTE: this routine should not be called for any of these cases:
      !   *) thetas = 0
      !   *) sigma  = 0
      ! These cases should be checked for before calling this function.
      
      ! The n=0 case can be written as:
      !   I(0,thetab,thetas,sigma) = EXP(-thetab) * I(0,0,thetas,sigma)
      ! An explicit calculation for second term is performed in the
      ! recurrence routine, with:
      !   I(0,0,thetas,sigma) = e^{x0^2-u0^2} erfc(x0) / erfc(u0)
      !   u0 = -1 / sqrt{2 sigma^2}
      !   x0 = (thetas*sigma^2 - 1) / sqrt{2 sigma^2}
      ! Evaluation of this formula is handled carefully there.
      IF (n .EQ. 0) THEN
        nintegrate = -thetab + recurrence(0,thetas,sigma)
        RETURN
      END IF
      
      ! Define several useful quantities
      sigma2 = sigma**2
      alpha  = sigma2*thetas
      beta   = sigma2*thetab
      eta    = sigma2*n
      r      = thetab/thetas
      kappa  = (1d0-alpha) + r
      lambda = (1d0-alpha) - r
      nkroot = SQRT(eta+0.25d0*kappa**2)
      
      ! We add an intermediate step here, making a change of variable:
      !   w = eps + thetab/thetas
      ! If we define W(w) = g(eps), then:
      !   W'(w)  = (eta/w + kappa - w) / sigma^2
      !   W''(w) = -(1 + eta/w^2) / sigma^2
      ! The solution to W'(w0) = 0 is:
      !   w0 = 1/2 kappa +/- sqrt{eta + 1/4 kappa^2}
      ! We want the positive root.  We are careful here to avoid
      ! cancelling between the terms:
      IF (kappa .GE. 0d0) THEN
        !w0 = SQRT(eta+0.25d0*kappa**2) + 0.5d0*kappa
        w0 = nkroot + 0.5d0*kappa
      ELSE
        !w0 = eta / (SQRT(eta+0.25d0*kappa**2) - 0.5d0*kappa)
        w0 = eta / (nkroot - 0.5d0*kappa)
      END IF
      ! Also careful with determining eps0 = w0 - r to avoid
      ! cancelling:
      IF (lambda .GE. 0d0) THEN
        !eps0 = SQRT(eta+0.25d0*kappa**2) + 0.5d0*lambda
        eps0 = nkroot + 0.5d0*lambda
      ELSE
        !eps0 = (eta + (1d0-alpha)*r) / (SQRT(eta+0.25d0*kappa**2) - 0.5d0*lambda)
        eps0 = (eta + (1d0-alpha)*r) / (nkroot - 0.5d0*lambda)
      END IF
      
      ! The scaling parameter and limit of integration
      IF (n .EQ. 0) THEN
        sigmau2 = sigma2
        sigmau  = sigma
      ELSE
        sigmau2 = sigma2 * w0**2 / (eta + w0**2)
        sigmau  = SQRT(sigmau2)
      END IF
      u0      = -eps0/sigmau
      
      ! Useful quantity
      IF (n .EQ. 0) THEN
        rho = 0d0
      ELSE
        rho = sigmau/w0
      END IF
      
      ! Now we determine:
      !   Iu = 1/sqrt{2\pi} \int_u0^\infty e^{f(u)}
      ! which is a well-behaved integral which can be evaluated
      ! numerically.
      ! 
      ! With the above definitions,
      !   f(u) = -u^2/2 + n [ln(1+rho*u) - rho*u + (rho*u)^2/2]
      ! Note that, due to canceling, the term in brackets is of
      ! order (rho*u)^3.
      ! 
      ! If n*rho^3 << 1, we can simply expand e^{n [...]} in u and
      ! perform the integral term by term.  These integrals have
      ! analytical solutions, allowing a numerical integration to
      ! be avoided (provided the series converges rapidly).
      IF ((ABS(n*rho**3) .LE. 1d-3) .AND. (u0 .LT. 25d0)) THEN
        lnIu = expand_lnIu(u0,n,rho)
      ELSE
        lnIu = nintegrate_lnIu(u0,n,rho)
      END IF
      
      ! Put everything together, divided as follows:
      !   I = [A sigmau/sigma] [e^{g(eps0)}] [Iu]
      !             C            e^g0         Iu
      ! Take log of above terms.
      Ainv = 0.5d0*ERFC(-1/(SQRT2*sigma))
      lnC  = LOG(sigmau/sigma) - LOG(Ainv)
      g0   = n*LOG(thetas*w0) - lngamma(n+1d0) - thetas*w0
     &       - (eps0-1)**2/(2*sigma2)
      nintegrate = lnC + g0 + lnIu
      
      ! TESTING
      !lnIu1 = nintegrate_lnIu(u0,n,rho)
      !lnIu2 = expand_lnIu(u0,n,rho)
      !WRITE(*,'(A,I5,3(1PG12.4),5(1PG))')
      !&    'DEBUG:',n,u0,rho,(n*rho**3)**1,
      !&    lnIu1,lnIu2,EXP(lnC+g0+lnIu1),EXP(lnC+g0+lnIu2),
      !&    EXP(lnIu2-lnIu1)-1d0
      
      END FUNCTION nintegrate
      
      
      ! ----------------------------------------------------------------
      ! Function to calculate the following integral numerically and
      ! return its logarithm:
      !   Iu   = 1/sqrt{2\pi} \int_u0^\infty e^{f(u)}
      !   f(u) = -u^2/2 + n [ln(1+rho*u) - rho*u + (rho*u)^2/2]
      ! Note that, due to canceling, the term in brackets is of
      ! order (rho*u)^3.
      ! 
      ! This integral, to an order of magnitude, scales as:
      !   Iu ~ (1/2)*erfc(u0/sqrt(2))
      ! so, for u0 not much larger than 1, it is within an an order of
      ! magnitude of 1 and should not have underflow/overflow issues
      ! (below, we treat the u0 > 0 case in a way that also avoids
      ! these issues).  For u0 < 0, will break this integral into two
      ! sections below: [u0,0] and [0,\infty) (we just use the single
      ! interval [u0,\infty) if u0 > 0).
      ! In each section, we further transform the integration variables
      ! to improve the convergence of the numerical integration.
      ! The methods here are partly based on the techniques found in
      ! Numerical recipes (in particular, see the "Quadrature by
      ! Variable Transformation" section).
      ! 
      ! The result is nearly always accurate to at least 10 significant
      ! digits, but is usually more accurate than that.  In some cases,
      ! it is accurate to nearly full double precision.
      ! 
      ! NOTE: Should avoid using this routine for n=0.
      ! 
      REAL*8 FUNCTION nintegrate_lnIu(u0,n,rho)
      IMPLICIT NONE
      INTEGER n
      REAL*8 u0,rho
      INTEGER K,KMAX,J,Nintervals1,Nintervals2
      PARAMETER(KMAX=30)
      REAL*8 h,t,tmax,c,q,u,dudt,f,fL,fR,sv,v,dvdt,hf,fu0,f1u0,f2u0,
     &       lndiffu0,sum,I1,I2,I1prev,I2prev
      REAL*8 PI,SQRT2PI
      PARAMETER(PI=3.1415926535897932d0)
      PARAMETER(SQRT2PI=2.5066282746310005d0)
      ! Relative error goal in integration.
      ! The integration routine below seems to converge very fast,
      ! so this can be set close to full double precision.
      REAL*8 reltol,stepreltol
      PARAMETER(reltol=1d-15)
      
      !WRITE(*,*) 'iIu arguments:',u0,n,rho
      
      ! The allowed difference between two loop iterations of the
      ! integration routine below at which the tolerance condition
      ! is considered met.  For well behaved functions, this type
      ! of integration as well as an appropriate variable
      ! transformation often leads to an approximate doubling
      ! of the number of significant digits with each loop
      ! (might require some burn in).  In this case, the difference
      ! between the current and previous iteration need only be
      ! less than the square root of the desired tolerance (the
      ! error is coming from previous iteration, while the current
      ! iteration should be at the desired tolerance).
      stepreltol = MAX(SQRT(reltol),1d-15)
      
      ! If the routine does not quite follow the "significant digits
      ! double every loop" behavior, we can use a more conservative
      ! difference between loops.
      !stepreltol = MAX(reltol**0.67,1d-15)
      
      ! Via comparison with Mathematica, accuracy does not appear to
      ! improve when setting stepreltol better than sqrt(reltol).
      
      ! Value of f(u) at lower bound of integration u0.  Used in
      ! a few places below.
      fu0 = -0.5d0*u0**2 + n*lndiff(rho*u0)
      
      ! Below, we divide the integration region into multiple parts:
      !   u0 > 0:       [u0,\infty)
      !   u0 < 0:       (u0,0] + [0,\infty)
      ! However, if f(u0) << -1, we can take the lower limit of
      ! integration to be -\infty, which allows for a slightly better
      ! numerical interation method:
      !   f(u0) << -1:  (-\infty,0] + [0,\infty)
      
      ! Numerical integral over [u0,\infty) ----------------
      ! This integration region does not contain the peak of f(u), so
      ! we make a further transformation:
      !   hf(v) = f(u) - f(u0)
      !   v     = (u-u0)/sv          => u = u0 + sv*v
      !   sv    = 1 / sqrt[|f'(u0)|^2 + |f''(u0)|]
      ! With the above transformation,
      !   hf(v=1) ~ -1
      ! and
      !   Iu = e^f(u0) sv \int_0^\infty dv e^hf(v)
      ! where the latter integral should be of order 1.
      ! Then, we use the transformation (see e.g. Numerical Recipes):
      !   v     = e^{t - e^{-t}}
      !   dv/dt = e^{t - e^{-t}} * (1 + e^{-t})
      !         = v * (1 + e^{-t})
      IF (u0 .GE. 0d0) THEN
        tmax = 5d0
        
        ! For small u0, the linear term dominates in the region of
        ! interest (where hf changes by ~1), while for large u0, the
        ! quadratic term dominates.  The form for sv here ensures
        ! the scaling is reasonable in either case.
        f1u0 = -u0*(1d0 - (n*rho**2) * (rho*u0)/(1d0+rho*u0))
        f2u0 = -1d0 + (n*rho**2) * (1d0 - 1d0/(1d0+rho*u0)**2)
        sv = 1/SQRT(f1u0**2 + ABS(f2u0))
        lndiffu0 = lndiff(rho*u0)
        
        ! Trapezoidal rule over 3 points with integrand zero
        ! at t = +/- tmax
        I1prev = 0d0
        Nintervals1 = 2
        h = 2*tmax / Nintervals1
        v = EXP(-1d0)
        u = u0 + sv*v
        dvdt = 2*v
        hf = -0.5d0*sv*v*(u+u0) + n*(lndiff(rho*u)-lndiffu0)
        sum = EXP(hf)*dvdt
        I1 = h * sum
        
        ! Double number of intervals until we reach desired tolerance
        DO K = 2, KMAX
          Nintervals1 = 2*Nintervals1
          h   = 2*tmax / Nintervals1
          sum = 0d0
          DO J = 1, Nintervals1
            t = (2*J-1)*h
            ! Points left of center (t -> -t)
            v = EXP(-t - EXP(t))
            dvdt = v * (1+EXP(t))
            u = u0 + sv*v
            hf = -0.5d0*sv*v*(u+u0) + n*(lndiff(rho*u)-lndiffu0)
            sum = sum + EXP(hf)*dvdt
            ! Points right of center
            v = EXP(t - EXP(-t))
            dvdt = v * (1+EXP(-t))
            u = u0 + sv*v
            hf = -0.5d0*sv*v*(u+u0) + n*(lndiff(rho*u)-lndiffu0)
            sum = sum + EXP(hf)*dvdt
          END DO
          I1prev = I1
          I1 = 0.5d0*I1prev + h*sum
          ! If we reach the desired tolerance, exit the loop
          IF (ABS(I1-I1prev) .LE. stepreltol*ABS(I1)) EXIT
        END DO
        !WRITE(*,*) 'I1fin: ',Nintervals1,EXP(fu0)*sv*I1
        !WRITE(*,*) 'Iu:    ',Nintervals1,EXP(fu0)*sv*I1,
       !&           EXP(fu0)*sv*I1/SQRT2PI
        nintegrate_lnIu = fu0 + LOG(sv*I1/SQRT2PI)
        RETURN
      END IF
      
      ! Numerical integral over [0,\infty) -----------------
      ! We use the transformation (see e.g. Numerical Recipes):
      !   u     = e^{t - e^{-t}}
      !   du/dt = e^{t - e^{-t}} * (1 + e^{-t})
      !         = u * (1 + e^{-t})
      tmax = 5d0
      
      ! Trapezoidal rule over 3 points with integrand zero
      ! at t = +/- tmax
      I1prev = 0d0
      Nintervals1 = 2
      h = 2*tmax / Nintervals1
      u = EXP(-1d0)
      dudt = 2*u
      f = -0.5d0*u**2 + n*lndiff(rho*u)
      sum = EXP(f)*dudt
      I1 = h * sum
      
      ! Double number of intervals until we reach desired tolerance
      DO K = 2, KMAX
        Nintervals1 = 2*Nintervals1
        h   = 2*tmax / Nintervals1
        sum = 0d0
        DO J = 1, Nintervals1
          t = (2*J-1)*h
          ! Points left of center (t -> -t)
          u = EXP(-t - EXP(t))
          dudt = u * (1+EXP(t))
          f = -0.5d0*u**2 + n*lndiff(rho*u)
          sum = sum + EXP(f)*dudt
          ! Points right of center
          u = EXP(t - EXP(-t))
          dudt = u * (1+EXP(-t))
          f = -0.5d0*u**2 + n*lndiff(rho*u)
          sum = sum + EXP(f)*dudt
        END DO
        I1prev = I1
        I1 = 0.5d0*I1prev + h*sum
        ! If we reach the desired tolerance, exit the loop
        IF (ABS(I1-I1prev) .LE. stepreltol*ABS(I1)) EXIT
      END DO
      !WRITE(*,*) 'I1inf: ',Nintervals1,I1
      
      ! Numerical integral over (-\infty,0] ----------------
      ! If EXP(f(u0)) is smaller than ~ double precision, then we
      ! treat the lower bound of integration as -\infty to improve
      ! convergence.
      ! We use the transformation (see e.g. Numerical Recipes):
      !   u     = - e^{t - e^{-t}}
      !   du/dt = - e^{t - e^{-t}} * (1 + e^{-t})
      !         = u * (1 + e^{-t})
      ! The following has the limits of integration reversed, so
      ! we need to change the sign afterwards.
      IF (fu0 .LE. -40d0) THEN
        tmax = 5d0
        
        ! Trapezoidal rule over 3 points with integrand zero
        ! at t = +/- tmax
        I2prev = 0d0
        Nintervals2 = 2
        h = 2*tmax / Nintervals2
        u = -EXP(-1d0)
        dudt = 2*u
        f = -0.5d0*u**2 + n*lndiff(rho*u)
        sum = EXP(f)*dudt
        I2 = h * sum
        
        ! Double number of intervals until we reach desired tolerance
        DO K = 2, KMAX
          Nintervals2 = 2*Nintervals2
          h   = 2*tmax / Nintervals2
          sum = 0d0
          DO J = 1, Nintervals2
            t = (2*J-1)*h
            ! Points left of center (t -> -t)
            u = -EXP(-t - EXP(t))
            dudt = u * (1+EXP(t))
            f = -0.5d0*u**2 + n*lndiff(rho*u)
            sum = sum + EXP(f)*dudt
            ! Points right of center
            u = -EXP(t - EXP(-t))
            dudt = u * (1+EXP(-t))
            f = -0.5d0*u**2 + n*lndiff(rho*u)
            sum = sum + EXP(f)*dudt
          END DO
          I2prev = I2
          I2 = 0.5d0*I2prev + h*sum
          ! If we reach the desired tolerance, exit the loop
          IF (ABS(I2-I2prev) .LE. stepreltol*ABS(I2)) EXIT
        END DO
        ! Flip sign due to reversed limits of integration in this
        ! region.
        I2 = -I2
        !WRITE(*,*) 'I2inf: ',Nintervals1,I2
        
      ! Numerical integral over [u0,0] ---------------------
      ! We use the transformation (see e.g. Numerical Recipes):
      !   u     = (u0/2) [1 - tanh(c sinh(t))]
      !         = u0 q / (1+q)
      !   du/dt = -(u0/2) sech^2(c sinh(t)) c cosh(t)
      !         = -2 u0 c cosh(t) q / (1+q)^2
      ! where:
      !   q = e^{-2 c sinh(t)}
      ! NOTE: The convergence here is likely slower than the
      !       above regions.  This transformation here may,
      !       in fact, not be close to optimal, but it seems
      !       to work reasonably well.
      ELSE
        ! Any c within (0,pi/2] is valid, but the largest value will
        ! be most useful in our case.
        c    = 0.5d0*PI
        tmax = 5d0
        
        ! Trapezoidal rule over 3 points with integrand zero
        ! at t = +/- tmax
        I2prev = 0d0
        Nintervals2 = 2
        h = 2*tmax / Nintervals2
        q = 1d0
        dudt = -0.5d0*u0*c
        u = u0 *  q/(1+q)
        f = -0.5d0*u**2 + n*lndiff(rho*u)
        sum = EXP(f)*dudt
        I2 = h * sum
        
        ! Double number of intervals until we reach desired tolerance
        DO K = 2, KMAX
          Nintervals2 = 2*Nintervals2
          h   = 2*tmax / Nintervals2
          sum = 0d0
          DO J = 1, Nintervals2
            t = (2*J-1)*h
            q = EXP(-2*c*SINH(t))
            dudt  = -2*u0*c*COSH(t) * q/(1+q)**2
            u = u0 * 1/(1+q)  ! Points left of center
            fL = -0.5d0*u**2 + n*lndiff(rho*u)
            u = u0 * q/(1+q)  ! Points right of center
            fR = -0.5d0*u**2 + n*lndiff(rho*u)
            sum = sum + (EXP(fL)+EXP(fR))*dudt
          END DO
          I2prev = I2
          I2 = 0.5d0*I2prev + h*sum
          ! If we reach the desired tolerance, exit the loop
          IF (ABS(I2-I2prev) .LE. stepreltol*ABS(I2)) EXIT
        END DO
        !WRITE(*,*) 'I2fin: ',Nintervals2,I2
      END IF
      
      !WRITE(*,*) 'Iu:    ',Nintervals1+Nintervals2,I1+I2,(I1+I2)/SQRT2PI
      
      ! Full integral --------------------------------------
      nintegrate_lnIu = LOG((I1+I2) / SQRT2PI)
      
      END FUNCTION nintegrate_lnIu
      
      
      ! ----------------------------------------------------------------
      ! Function to calculate the following integral and
      ! return its logarithm:
      !   Iu   = 1/sqrt{2\pi} \int_u0^\infty e^{f(u)}
      !   f(u) = -u^2/2 + n [ln(1+rho*u) - rho*u + (rho*u)^2/2]
      ! where the integral is determined by expanding the non- u^2
      ! terms in the exponential in a power series and evaluating
      ! the integrals analytically.
      ! NOTE: This routine is only for small values of n*rho^3 (for
      !       convergence) and u0 < 25 (for underflow issues).
      !       Approximately double precision should be achieved for
      !       n*rho^3 < 1e-3.
      ! 
      REAL*8 FUNCTION expand_lnIu(u0,n,rho)
      IMPLICIT NONE
      INTEGER n
      REAL*8 u0,rho
      INTEGER k
      REAL*8 C0,C1,Ckeven,Ckodd,Dk,term,sum
      INTEGER KMAX
      PARAMETER(KMAX=16)
      REAL*8 b(3:KMAX)
      REAL*8 SQRT2PI,SQRT2,SQRTHALF
      PARAMETER(SQRT2PI=2.5066282746310005d0)
      PARAMETER(SQRT2=1.4142135623730950d0)
      PARAMETER(SQRTHALF=0.70710678118654752d0)
      
      !WRITE(*,*) 'eIu arguments:',u0,n,rho
      
      ! The exponential argument can be written as:
      !   f(u) = -u^2 + \Sum_{k=3}^\infty b_k u^k
      !   b_k = n * (-1)^{k+1} rho^k / k
      ! The integral is then expanded in the form:
      !   1/sqrt{2\pi} \int_u0^\infty e^{-u^2/2} [1 + \Sum_k A_k u^k]
      !     = C_0 +  \Sum_k A_k C_k
      ! where:
      !   C_k = 1/sqrt{2\pi} \int_u0^\infty u^k e^{-u^2/2}
      ! Note that the sum above starts at k=3.  The factors A_k are
      ! functions of n & rho.  The C_k have analytic forms given by:
      !   C_0 = (1/2)*erfc(u0/sqrt{2})
      !   C_1 = 1/sqrt{2\pi} e^{-u0^2/2}
      !   C_{k+2} = (k+1)*C_k + u0^{k+1} C_1
      ! Below, we use D_k = u0^{k-1} C_1.
      
      ! Convergence:
      ! Terms in the expansion are decreasing in general, but not
      ! monotonically!  Care should thus be taken at which points to
      ! check for convergence.  Considerations:
      !   *) For u0 << 0, the odd terms are suppressed
      !   *) Terms scale as ~ k!! * n^Floor(k/3) * rho^k
      !   *) rho < 1/sqrt{n}  => n*rho^2 < 1
      
      ! k=0
      C0     = 0.5d0*ERFC(SQRTHALF*u0)
      Ckeven = C0
      sum    = C0
      
      ! Special case: no further expansion
      IF ((n .EQ. 0) .OR. (rho .EQ. 0d0)) THEN
        expand_lnIu = LOG(sum)
        RETURN
      END IF
      
      ! Calculate exponential argument expansion coefficients
      b(3) = rho**3
      DO K = 4,KMAX
        b(K) = (-rho)*b(K-1)
      END DO
      DO K = 3,KMAX
        b(K) = n*b(K)/K
      END DO
      
      ! k=1 (not present in expansion)
      Dk     = EXP(-0.5d0*u0**2) / SQRT2PI
      C1     = Dk
      Ckodd  = C1
      
      ! k=2 (not present in expansion)
      Dk     = u0*Dk
      Ckeven = (2-1)*Ckeven + Dk
      
      ! k=3
      Dk     = u0*Dk
      Ckodd  = (3-1)*Ckodd + Dk
      term   = b(3) * Ckodd
      sum    = sum + term
      
      ! k=4
      Dk     = u0*Dk
      Ckeven = (4-1)*Ckeven + Dk
      term   = b(4) * Ckeven
      sum    = sum + term
      
      ! k=5
      Dk     = u0*Dk
      Ckodd  = (5-1)*Ckodd + Dk
      term   = b(5) * Ckodd
      sum    = sum + term
      
      ! k=6
      Dk     = u0*Dk
      Ckeven = (6-1)*Ckeven + Dk
      term   = (b(6) + 0.5d0*b(3)**2) * Ckeven
      sum    = sum + term
      
      ! Check if converged here
      IF (ABS(term) .LT. 1d-16*sum) THEN
        expand_lnIu = LOG(sum)
        RETURN
      END IF
      
      ! k=7
      Dk     = u0*Dk
      Ckodd  = (7-1)*Ckodd + Dk
      term   = (b(7) + b(3)*b(4)) * Ckodd
      sum    = sum + term
      
      ! k=8
      Dk     = u0*Dk
      Ckeven = (8-1)*Ckeven + Dk
      term   = (b(8) + b(3)*b(5) + 0.5d0*b(4)**2) * Ckeven
      sum    = sum + term
      
      ! k=9
      Dk     = u0*Dk
      Ckodd  = (9-1)*Ckodd + Dk
      term   = (b(9) + b(3)*b(6) + b(4)*b(5)
     &          + (1d0/6d0)*b(3)**3) * Ckodd
      sum    = sum + term
      
      ! k=10
      Dk     = u0*Dk
      Ckeven = (10-1)*Ckeven + Dk
      term   = (b(10) + b(3)*b(7) + b(4)*b(6) + 0.5d0*b(5)**2
     &          + 0.5d0*b(3)**2*b(4)) * Ckeven
      sum    = sum + term
      
      ! k=11
      Dk     = u0*Dk
      Ckodd  = (11-1)*Ckodd + Dk
      term   = (b(11) + b(3)*b(8) + b(4)*b(7) + b(5)*b(6)
     &          + 0.5d0*b(3)**2*b(5) + 0.5d0*b(3)*b(4)**2) * Ckodd
      sum    = sum + term
      
      ! k = 12
      Dk     = u0*Dk
      Ckeven = (12-1)*Ckeven + Dk
      term   = (b(12) + b(3)*b(9) + b(4)*b(8) + b(5)*b(7)
     &                + 0.5d0*b(6)**2
     &          + 0.5d0*b(3)**2*b(6) + b(3)*b(4)*b(5)
     &                + (1d0/6d0)*b(4)**3
     &          + (1d0/24d0)*b(3)**4) * Ckeven
      sum  = sum + term
      
      ! Check if converged here
      IF (ABS(term) .LT. 1d-16*sum) THEN
        expand_lnIu = LOG(sum)
        RETURN
      END IF
      
      ! k = 13
      Dk     = u0*Dk
      Ckodd  = (13-1)*Ckodd + Dk
      term = (b(13) + b(3)*b(10) + b(4)*b(9) + b(5)*b(8) + b(6)*b(7)
     &        + 0.5d0*b(3)**2*b(7) + b(3)*b(4)*b(6) + 0.5d0*b(3)*b(5)**2
     &              + 0.5d0*b(4)**2*b(5)
     &        + (1d0/6d0)*b(3)**3*b(4)
     &       ) * Ckodd
      sum  = sum + term
      
      ! k = 14
      Dk     = u0*Dk
      Ckeven = (14-1)*Ckeven + Dk
      term   = (b(14) + b(3)*b(11) + b(4)*b(10) + b(5)*b(9) + b(6)*b(8)
     &                + 0.5d0*b(7)**2
     &          + 0.5d0*b(3)**2*b(8) + b(3)*b(4)*b(7) + b(3)*b(5)*b(6)
     &                + 0.5d0*b(4)**2*b(6) + 0.5d0*b(4)*b(5)**2
     &          + (1d0/6d0)*b(3)**3*b(5) + 0.25d0*b(3)**2*b(4)**2
     &         ) * Ckeven
      sum  = sum + term
      
      ! k = 15
      Dk     = u0*Dk
      Ckodd  = (15-1)*Ckodd + Dk
      term = (b(15) + b(3)*b(12) + b(4)*b(11) + b(5)*b(10) + b(6)*b(9)
     &              + b(7)*b(8)
     &        + 0.5d0*b(3)**2*b(9) + b(3)*b(4)*b(8) + b(3)*b(5)*b(7)
     &              + 0.5d0*b(3)*b(6)**2 + 0.5d0*b(4)**2*b(7)
     &              + b(4)*b(5)*b(6) + (1d0/6d0)*b(5)**3
     &        + (1d0/6d0)*b(3)**3*b(6) + 0.5d0*b(3)**2*b(4)*b(5)
     &              + (1d0/6d0)*b(3)*b(4)**3
     &        + (1d0/120d0)*b(3)**5
     &       ) * Ckodd
      sum  = sum + term
      
      ! k = 16
      Dk     = u0*Dk
      Ckeven = (16-1)*Ckeven + Dk
      term   = (b(16) + b(3)*b(13) + b(4)*b(12) + b(5)*b(11)
     &                + b(6)*b(10) + b(7)*b(9) + 0.5d0*b(8)**2
     &          + 0.5d0*b(3)**2*b(10) + b(3)*b(4)*b(9) + b(3)*b(5)*b(8)
     &                + b(3)*b(6)*b(7) + 0.5d0*b(4)**2*b(8)
     &                + b(4)*b(5)*b(7) + 0.5d0*b(4)*b(6)**2
     &                + 0.5d0*b(5)**2*b(6)
     &          + (1d0/6d0)*b(3)**3*b(7) + 0.5d0*b(3)**2*b(4)*b(6)
     &                + 0.25d0*b(3)**2*b(5)**2 + 0.5d0*b(3)*b(4)**2*b(5)
     &                + (1d0/24d0)*b(4)**4
     &          + (1d0/24d0)*b(3)**4*b(4)
     &         ) * Ckeven
      sum  = sum + term
      
      ! Remaining terms are of order n^5*rho^17 or n^6*rho^18 with
      ! Ck term of order k!!.  In the worst case (n*rho^2 = 1), the
      ! next term is of order ~ 100/n^3 relative to the sum.
      ! If n*rho^2 < 1e-3, then approximately double precision is
      ! achieved at this point.
      
      expand_lnIu = LOG(sum)
      
      END FUNCTION expand_lnIu
      
      
      ! ----------------------------------------------------------------
      ! Function to calculate the quantity:
      !   ln(1+z) - (z - (1/2)*z^2)
      ! Small z must be handled carefully as the two terms on the
      ! right correspond to the first two terms of the expansion
      ! of ln(1+z) about zero.
      ! 
      ! Precision is ~ 100*EPSILON.
      ! 
      REAL*8 FUNCTION lndiff(z)
      IMPLICIT NONE
      REAL*8 z,zabs,zk
      INTEGER k,klast
      
      ! Special case
      IF (z .LE. -1d0) THEN
        lndiff = -HUGE(lndiff)
        RETURN
      END IF
      
      zabs = ABS(z)
      
      ! If z is not too small, we can just evaluate explicitly
      IF (zabs .GT. 0.2d0) THEN
        lndiff = LOG(1+z) - z*(1 - 0.5d0*z)
        RETURN
      END IF
      
      ! We will use an expansion about zero, canceling the leading
      ! terms.  We keep terms in the expansion up to x^klast/klast.
      ! klast is chosen below to give almost double precision.
      ! Precision can be as poor as ~ 100*EPSILON, but is better for
      ! smaller z.  The expansion is:
      !     \sum_{k=3}^{\infty} (-1)^{k+1} z^k / k
      ! The precision is approximately:
      !     3*x^{klast-2} / (klast+1)
      ! where klast is the last term in the sum.
      
      IF (zabs .LT. 1d-4) THEN
        klast = 6
      ELSE IF (zabs .LT. 0.01d0) THEN
        klast = 10
      ELSE IF (zabs .LT. 0.05d0) THEN
        klast = 14
      ELSE IF (zabs .LT. 0.14d0) THEN
        klast = 18
      ELSE
        klast = 22
      END IF
      
      ! Use expansion about zero
      zk     = z**3
      lndiff = zk/3
      DO k = 4,klast
        zk     = -z*zk
        lndiff = lndiff + zk/k
      END DO
      
      END FUNCTION lndiff
      
      END FUNCTION


