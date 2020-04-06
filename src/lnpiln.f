!=======================================================================
! LN(Poisson Integral) for Log-Normal systematic uncertainty
!
! The logarithm of a Poissonian likelihood marginalized over a
! systematic uncertainty in the average expected number of signal
! events (e.g. from a systematic uncertainty in the detector
! sensitivity).  The systematic uncertainty is taken to be a
! log-normal distributed fractional variation in the signal
! contribution to the Poisson distribution average.  The
! marginalized likelihood is given by:
!     I(n,thetab,thetas,sigma)
!         = 1 / sqrt{2 \pi \sigma^2}
!           \int_0^\infty d\epsilon e^{-(ln(\epsilon))^2 / 2\sigma^2}
!                         (\theta_b + \epsilon\theta_s)^n
!                         e^{-(\theta_b + \epsilon\theta_s)}
!                         / (\epsilon n!)
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
!
! For further details, see:
!   * Conrad, Botner, Hallgren & Perez de los Heros, Phys. Rev. D 67,
!     012002 (2003) [arXiv:hep-ex/0202013].
!   * Scott et al., JCAP 1001, 031 (2010) [arXiv:0909.3300
!     [astro-ph.CO]].
!
! This routine uses a either a conditioned numerical integration or an
! analytically integrated Taylor series to perform the calculation
! quickly.  It has not been thoroughly tested for extreme cases in
! the parameters (lnpin is fairly robust in this respect); however,
! it works well for "reasonable" parameter values.
!
!
!    Created by Chris Savage (savage@fysik.su.se)
!    2011/06/02
!
!=======================================================================
!
      REAL*8 FUNCTION nulike_lnpiln(n,thetab,thetas,sigma)
      IMPLICIT NONE
      INTEGER n
      REAL*8 thetab,thetas,sigma

      ! Checks --------------------------------------------
      IF ((n .LT. 0) .OR. (thetab .LT. 0d0)
     &    .OR. (thetas .LT. 0d0) .OR. (sigma .LT. 0d0)) THEN
        WRITE(*,*) 'ERROR: nulike_lnpiln called with negative argument'
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
        nulike_lnpiln = lnpoisson(n,thetab)
        RETURN
      END IF

      ! Special case: sigma = 0
      ! The integration drops out and we have a strictly Poisson
      ! probability with average thetab+thetas
      IF (sigma .EQ. 0d0) THEN
        nulike_lnpiln = lnpoisson(n,thetab+thetas)
        RETURN
      END IF

      ! General case ---------------------------------------
      ! For the general case, we will do a numerical integration.
      nulike_lnpiln = nintegrate(n,thetab,thetas,sigma)


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
      ! Function to calculate the integral numerically.
      !
      ! To obtain a better form for numerical integration, we choose
      ! a more suitable integration variable.  First, write the
      ! integral as:
      !   I(n,thetab,thetas,sigma)
      !        = 1/\sqrt{2\pi sigma^2} \int_{-\infty}^\infty dy e^{g(y)}
      !   g(y) = n ln(thetab+thetas*e^y) - ln(n!) - (thetab+thetas*e^y)
      !            - y^2/(2*sigma^2)
      !   y    = ln(\epsilon)
      ! Now, define y0 and sigmau such that:
      !   g'(y0)  = 0
      !   g''(y0) = -1/sigmau^2
      ! Then we change the integration variable to:
      !   u = (y - y0)/sigmau
      ! Under this transformation:
      !   I  = (sigmau/sigma) e^{g(y0)} Iu
      !   Iu = 1/sqrt(2\pi) \int_{-\infty}^\infty du e^{f(u)}
      ! with:
      !   f(u) = g(y) - g(y0)
      ! Note that:
      !   f(u) = -u^2/2 + O(u^3)
      ! so, to first order, the integration is now approximately a
      ! gaussian integration with unit variance.  The quantity Iu
      ! should be of order 1.
      !
      ! NOTE: This routine has not been checked for proper handling
      !       under extreme cases (very large and/or very small
      !       parameters).
      !
      REAL*8 FUNCTION nintegrate(n,thetab,thetas,sigma)
      IMPLICIT NONE
      INTEGER n
      REAL*8 thetab,thetas,sigma
      REAL*8 y0,g0,sigmau,expy,kappa,rho,lnC,lnIu

      ! NOTE: this routine should not be called for any of these cases:
      !   *) thetas = 0
      !   *) sigma  = 0
      ! These cases should be checked for before calling this function.

      ! Attempt to find the maximum of g(y) and the corresponding y0, g(y0),
      ! and sigmau.
      IF (.NOT. find_maximum(n,thetab,thetas,sigma,y0,g0,sigmau)) THEN
        ! If attempting to find the maximum fails silently, this indicates
        ! that thetas is approximately zero and we should fall back to the
        ! standard Poisson case.
        nintegrate = lnpoisson(n,thetab)
        RETURN
      ENDIF

      ! Useful quantities
      expy  = EXP(y0)
      kappa = thetas*expy
      rho   = kappa / (thetab+kappa)    ! also expy / (r+expy)

      ! Now we determine:
      !   Iu = 1/sqrt{2\pi} \int_u0^\infty e^{f(u)}
      ! which is a well-behaved integral which can be evaluated
      ! numerically.
      !
      ! With the above definitions,
      !   f(u) = -u^2/2 - kappa [e^(dy) - 1 - dy - dy^2/2]
      !          + n [ln(1+rho*(e^dy-1)) - rho*dy - rho*(1-rho)/2 dy^2]
      ! where dy = sigmau*u.  Note that, due to canceling, the terms in
      ! brackets are each of order u^3.
      !
      ! If kappa, n and sigmau are sufficiently small, we can simply
      ! expand the non- u^2 part of the exponential in u and
      ! perform the integral term by term.  These integrals have
      ! analytical solutions, allowing a numerical integration to
      ! be avoided (provided the series converges rapidly).
      ! NOTE: The check is made against _one_ of the contributions to
      !       the next term in the sum.  This is only a guesstimate of
      !       where this expansion is appropriate!  It may be
      !       over/under-estimating the parameter space where the
      !       convergence is good.
      ! Note that, despite the comparison here, the expansion tends to
      ! give better than 1d-10 precision except in cases where the
      ! numerical integration also fails to give better.
      IF (((n*rho*(1-rho*(3-2*rho))-kappa)*sigmau**3)**6
     &    .LE. 1d-10) THEN
        lnIu = expand_lnIu(sigmau,n,rho,kappa)
        ! A very large value indicates the expansion failed;
        ! revert to the slower numerical integration in this case
        IF (lnIu .GT. 1d100) lnIu = nintegrate_lnIu(sigmau,n,rho,kappa)
      ELSE
        lnIu = nintegrate_lnIu(sigmau,n,rho,kappa)
      END IF

      ! Put everything together, divided as follows:
      !   I = [sigmau/sigma] [e^{g(eps0)}] [Iu]
      !           C            e^g0         Iu
      ! Take log of above terms.
      lnC  = LOG(sigmau/sigma)
      nintegrate = lnC + g0 + lnIu

      ! TESTING
      !lnIu1 = nintegrate_lnIu(sigmau,n,rho,kappa)
      !lnIu2 = expand_lnIu(sigmau,n,rho,kappa)
      !WRITE(*,'(A,I5,5(1PG12.4),5(1PG))')
      !&    'DEBUG:',n,y0,sigmau,rho,kappa,
      !&    ((n*rho*(1-rho*(3-2*rho))-kappa)*sigmau**3)**6,
      !&    lnIu1,lnIu2,EXP(lnC+g0+lnIu1),EXP(lnC+g0+lnIu2),
      !&    EXP(lnIu2-lnIu1)-1d0

      END FUNCTION nintegrate


      ! ----------------------------------------------------------------
      ! Function to calculate the following integral numerically and
      ! return its logarithm:
      !   Iu   = 1/sqrt{2\pi} \int_{-\infty}^\infty e^{f(u)}
      !   f(u) = -u^2/2 - kappa [e^(dy) - 1 - dy - dy^2/2]
      !          + n [ln(1+rho*(e^dy-1)) - rho*dy - rho*(1-rho)/2 dy^2]
      ! where dy = sigmau*u.  Note that, due to canceling, the terms in
      ! brackets are each of order u^3.
      !
      ! To a rough approximation, this is a Gaussian integration and
      ! should be of O(1); in many (most?) cases, this integral is
      ! actually very close to 1 (within 1% or better).  We break the
      ! integration region below into (-\infty,0] and [0,+\infty).
      ! In each section, we transform the integration variables
      ! to improve the convergence of the numerical integration.
      ! The methods here are partly based on the techniques found in
      ! Numerical recipes (in particular, see the "Quadrature by
      ! Variable Transformation" section).
      !
      ! The result is nearly always accurate to at least 10 significant
      ! digits, but is usually more accurate than that.  In some cases,
      ! it is accurate to nearly full double precision.
      !
      ! NOTE: This routine has not been checked for proper handling
      !       under extreme cases (very large and/or very small
      !       parameters).
      !
      ! NOTE2: The nearly unit value obtained for this integral
      !        suggests there might be a power expansion in the
      !        non-quadratic part of the exponential, i.e.:
      !            e^f(u) = e^{-u^2/2} [1 + \Sum_{k=3}^\infty C_k u^k]
      !        The terms in the expansion have closed forms for their
      !        integrals.  For a well-behaved power expansion, summing
      !        over the integrated terms in the power expansion can
      !        yield the value for the integral much more rapidly than
      !        the numerical integration done here.  This expansion
      !        has been implemented in the expand_lnIu routine below.
      !
      REAL*8 FUNCTION nintegrate_lnIu(sigmau,n,rho,kappa)
      IMPLICIT NONE
      INTEGER n
      REAL*8 sigmau,rho,kappa
      INTEGER K,KMAX,J,Nintervals1,Nintervals2
      !PARAMETER(KMAX=30)
      PARAMETER(KMAX=10)
      REAL*8 rho0,h,t,tmax,u,dudt,f,sum,I1,I2,I1prev,I2prev
      REAL*8 PI,SQRT2PI
      PARAMETER(PI=3.1415926535897932d0)
      PARAMETER(SQRT2PI=2.5066282746310005d0)
      ! Relative error goal in integration.
      ! The integration routine below seems to converge very fast,
      ! so this can be set close to full double precision.
      REAL*8 reltol,stepreltol
      PARAMETER(reltol=1d-15)

      !WRITE(*,'(A,G,I,2(G))') 'iIu arguments:',sigmau,n,rho,kappa

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

      ! For the case n=0, we use rho=0 below, which disables
      ! calculation of an irrelevant term
      IF (n .EQ. 0) THEN
        rho0 = 0d0
      ELSE
        rho0 = rho
      END IF

      ! Numerical integral over (-\infty,0] ----------------
      ! We use the transformation (see e.g. Numerical Recipes):
      !   u     = - e^{t - e^{-t}}
      !   du/dt = - e^{t - e^{-t}} * (1 + e^{-t})
      !         = u * (1 + e^{-t})
      ! The following has the limits of integration reversed, so
      ! we need to change the sign afterwards.
      tmax = 5d0

      ! Trapezoidal rule over 3 points with integrand zero
      ! at t = +/- tmax
      I1prev = 0d0
      Nintervals1 = 2
      h = 2*tmax / Nintervals1
      u = -EXP(-1d0)
      dudt = 2*u
      f = -0.5d0*u**2 - kappa*expdiff2(sigmau*u)
     &                + n*lnexpdiff2(rho0,sigmau*u)
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
          u = -EXP(-t - EXP(t))
          dudt = u * (1+EXP(t))
          f = -0.5d0*u**2 - kappa*expdiff2(sigmau*u)
     &                    + n*lnexpdiff2(rho0,sigmau*u)
          sum = sum + EXP(f)*dudt
          ! Points right of center
          u = -EXP(t - EXP(-t))
          dudt = u * (1+EXP(-t))
          f = -0.5d0*u**2 - kappa*expdiff2(sigmau*u)
     &                    + n*lnexpdiff2(rho0,sigmau*u)
          sum = sum + EXP(f)*dudt
        END DO
        I1prev = I1
        I1 = 0.5d0*I1prev + h*sum
        ! If we reach the desired tolerance, exit the loop
        IF (ABS(I1-I1prev) .LE. stepreltol*ABS(I1)) EXIT
      END DO
      ! Flip sign due to reversed limits of integration in this
      ! region.
      I1 = -I1
      !WRITE(*,*) 'I1inf: ',Nintervals1,I1

      ! Numerical integral over [0,\infty) -----------------
      ! We use the transformation (see e.g. Numerical Recipes):
      !   u     = e^{t - e^{-t}}
      !   du/dt = e^{t - e^{-t}} * (1 + e^{-t})
      !         = u * (1 + e^{-t})
      tmax = 5d0

      ! Trapezoidal rule over 3 points with integrand zero
      ! at t = +/- tmax
      I2prev = 0d0
      Nintervals2 = 2
      h = 2*tmax / Nintervals2
      u = EXP(-1d0)
      dudt = 2*u
      f = -0.5d0*u**2 - kappa*expdiff2(sigmau*u)
     &                + n*lnexpdiff2(rho0,sigmau*u)
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
          u = EXP(-t - EXP(t))
          dudt = u * (1+EXP(t))
          f = -0.5d0*u**2 - kappa*expdiff2(sigmau*u)
     &                    + n*lnexpdiff2(rho0,sigmau*u)
          sum = sum + EXP(f)*dudt
          ! Points right of center
          u = EXP(t - EXP(-t))
          dudt = u * (1+EXP(-t))
          f = -0.5d0*u**2 - kappa*expdiff2(sigmau*u)
     &                    + n*lnexpdiff2(rho0,sigmau*u)
          sum = sum + EXP(f)*dudt
        END DO
        I2prev = I2
        I2 = 0.5d0*I2prev + h*sum
        ! If we reach the desired tolerance, exit the loop
        IF (ABS(I2-I2prev) .LE. stepreltol*ABS(I2)) EXIT
      END DO
      !WRITE(*,*) 'I2inf: ',Nintervals2,I2

      ! Full integral --------------------------------------
      nintegrate_lnIu = LOG((I1+I2) / SQRT2PI)

      !WRITE(*,*) 'Iu:    ',Nintervals1+Nintervals2,I1+I2,(I1+I2)/SQRT2PI

      END FUNCTION nintegrate_lnIu


      ! ----------------------------------------------------------------
      ! Function to calculate the following integral and
      ! return its logarithm:
      !   Iu   = 1/sqrt{2\pi} \int_{-\infty}^\infty e^{f(u)}
      !   f(u) = -u^2/2 - kappa [e^(dy) - 1 - dy - dy^2/2]
      !          + n [ln(1+rho*(e^dy-1)) - rho*dy - rho*(1-rho)/2 dy^2]
      ! with dy = sigmau*u, where the integral is determined by
      ! expanding the non- u^2 terms in the exponential in a power
      ! series and evaluating the integrals analytically.  A very large
      ! number (> 1d100) indicates this series did not converge to a
      ! reasonable precision.
      !
      ! NOTE: This routine is only for small values of sigmau, rho,
      !       and/or kappa (for convergence).  The exact value of
      !       "small" is not well-defined, but:
      !         (((n*rho*(1-rho*(3-2*rho))-kappa)*sigmau**3)**6 < 1d-10
      !       is often reasonable, where the left side is _one_
      !       contribution to the next term in the series.  Other
      !       contributions may be signifcantly larger, so this is
      !       not a definitive estimate of convergence and there are
      !       instances where is grossly incorrect.  If the series
      !       does not converge to at least 1d-10 precision in the
      !       limited number of terms included in the calculation, a
      !       very large number (HUGE(1d0)) is returned instead as an
      !       indicator of an error.  In many cases, the precision
      !       is much better than 1d-10 but, in cases where it is not,
      !       the numerical integral also appears to have reduced
      !       precision.
      !
      REAL*8 FUNCTION expand_lnIu(sigmau,n,rho,kappa)
      IMPLICIT NONE
      INTEGER n
      REAL*8 sigmau,rho,kappa
      INTEGER K
      REAL*8 kfact,sigmauk,Ck,term,sum
      INTEGER KMAX
      PARAMETER(KMAX=16)
      REAL*8 P(0:KMAX-1),W(0:KMAX-1,2:KMAX-1),b(1:KMAX)
      ! W(j,i) = (-1)^j W_{i,j} where W_{i,j} are Worpitzky numbers
      ! Note recurrence relation:  W_{i,j} = (j+1)*W_{i-1,j} + j*W_{i-1,j-1}
      SAVE W
      ! Initialize W array
      !      j=    -0- -1- -2-   -3-  -4-   -5- -6--7-8-9-0-1-2-3-4-5-
      !DATA W(:,0) / 1,  0,  0,    0,   0,    0,  0, 0,0,0,0,0,0,0,0,0 /
      !DATA W(:,1) / 1, -1,  0,    0,   0,    0,  0, 0,0,0,0,0,0,0,0,0 /
      DATA W(:,2) / 1, -3,  2,    0,   0,    0,  0, 0,0,0,0,0,0,0,0,0 /
      DATA W(:,3) / 1, -7, 12,   -6,   0,    0,  0, 0,0,0,0,0,0,0,0,0 /
      DATA W(:,4) / 1,-15, 50,  -60,  24,    0,  0, 0,0,0,0,0,0,0,0,0 /
      DATA W(:,5) / 1,-31,180, -390, 360, -120,  0, 0,0,0,0,0,0,0,0,0 /
      DATA W(:,6) / 1,-63,602,-2100,3360,-2520,720, 0,0,0,0,0,0,0,0,0 /
      DATA W(:,7)  / 1, -127, 1932, -10206,  25200,   -31920,   20160,
     &                   -5040,       0,        0,      0,0,0,0,0,0 /
      DATA W(:,8)  / 1, -255, 6050, -46620, 166824,  -317520,  332640,
     &                 -181440,   40320,        0,      0,0,0,0,0,0 /
      DATA W(:,9)  / 1, -511,18660,-204630,1020600, -2739240, 4233600,
     &                -3780000, 1814400,  -362880,      0,0,0,0,0,0 /
      DATA W(:,10) / 1,-1023,57002,-874500,5921520,-21538440,46070640,
     &               -59875200,46569600,-19958400,3628800,0,0,0,0,0 /
      DATA W(:,11) / 1, -2047,  173052d0,   -3669006d0,   33105600d0,
     &              -158838240d0,      451725120d0,     -801496080d0,
     &               898128000d0,     -618710400d0,      239500800d0,
     &               -39916800d0,              0d0,              0d0,
     &                       0d0,              0d0 /
      DATA W(:,12) / 1, -4095,  523250d0,  -15195180d0,  180204024d0,
     &             -1118557440d0,     4115105280d0,    -9574044480d0,
     &             14495120640d0,   -14270256000d0,     8821612800d0,
     &             -3113510400d0,      479001600d0,              0d0,
     &                       0d0,              0d0 /
      DATA W(:,13) / 1, -8191, 1577940d0,  -62350470d0,  961800840d0,
     &             -7612364760d0,    35517081600d0,  -105398092800d0,
     &            207048441600d0,  -273158645760d0,   239740300800d0,
     &           -134399865600d0,    43589145600d0,    -6227020800d0,
     &                       0d0,              0d0 /
      DATA W(:,14) / 1,-16383, 4750202d0, -254135700d0, 5058406080d0,
     &            -50483192760d0,   294293759760d0, -1091804313600d0,
     &           2706620716800d0, -4595022432000d0,  5368729766400d0,
     &          -4249941696000d0,  2179457280000d0,  -653837184000d0,
     &             87178291200d0,              0d0 /
      DATA W(:,15) / 1,-32767,14283372d0,-1030793406d0,26308573200d0,
     &           -328191186960d0,  2362955474880d0,-10794490827120d0,
     &          33094020960000d0,-70309810771200d0,105006251750400d0,
     &        -110055327782400d0, 79332244992000d0,-37486665216000d0,
     &          10461394944000d0, -1307674368000d0 /

      !WRITE(*,'(A,G,I,2(G))') 'eIu arguments:',sigmau,n,rho,kappa

      ! The exponential argument can be written as:
      !   f(u) = -u^2 + \Sum_{k=3}^\infty b_k u^k
      !   b_k = (n*rho P_{k-1}.W{k-l} - kappa) sigmau^k / k!
      ! where:
      !   P_{k-1} = [1,rho,rho^2,...,rho^{k-1}]
      !   W_{k-1} = [W_{k-1,0}, (-1)^1 W_{k-1,1}, ..., (-1)^{k-1} W_{k-1,k-1}]
      ! and W_{i,j} are the Worpitzky numbers:
      !   W_{i,j} = \Sum_{v=0}^j (-1)^{v+j} (1+v)^i j! / (v!(j-v)!)
      ! The integral is then expanded in the form:
      !   1/sqrt{2\pi} \int_{-\infty}^\infty e^{-u^2/2} [1 + \Sum_k A_k u^k]
      !     = 1 +  \Sum_{k=3} A_k C_k
      ! where:
      !   C_k = 1/sqrt{2\pi} \int_{-\infty}^\infty u^k e^{-u^2/2}
      ! Note that the sum above starts at k=3.  The factors A_k are
      ! functions of sigma0, n, rho & kappa.  The C_k have analytic
      ! forms given by:
      !   C_0 = 1
      !   C_1 = 0
      !   C_{k+2} = (k+1)*C_k
      ! From the above, only the even terms contribute to the integral.

      ! Convergence:
      ! Terms in the expansion are decreasing in general, but not
      ! monotonically!  Care should thus be taken at which points to
      ! check for convergence.

      ! Special case: no expansion
      IF ((n .EQ. 0) .AND. (kappa .EQ. 0d0)) THEN
        expand_lnIu = 0d0
        RETURN
      END IF

      ! Calculate exponential argument expansion coefficients
      P(0) = 1d0
      P(1) = rho
      sigmauk = sigmau**2
      kfact = 2d0
      DO K = 3,KMAX
        P(K-1) = rho*P(K-2)
        sigmauk = sigmau*sigmauk
        kfact = K*kfact
        b(K) = (n*rho*DOT_PRODUCT(P(0:K-1),W(0:K-1,K-1)) - kappa )
     &         * sigmauk / kfact
      END DO

      ! Now sum over terms in expansion
      ! k = 0
      Ck  = 1d0
      sum = 1d0

      ! k = 2 (not present in expansion)
      !Ck  = (2-1)*Ck

      ! k = 4
      Ck   = (4-1)*Ck
      term = b(4) * Ck
      sum  = sum + term

      ! k = 6
      Ck   = (6-1)*Ck
      term = (b(6) + 0.5d0*b(3)**2) * Ck
      sum  = sum + term

      ! Check if converged here
      IF (ABS(term) .LT. 1d-16*sum) THEN
        expand_lnIu = LOG(sum)
        RETURN
      END IF

      ! k = 8
      Ck   = (8-1)*Ck
      term = (b(8) + b(3)*b(5) + 0.5d0*b(4)**2) * Ck
      sum  = sum + term

      ! k = 10
      Ck   = (10-1)*Ck
      term = (b(10) + b(3)*b(7) + b(4)*b(6) + 0.5d0*b(5)**2
     &        + 0.5d0*b(3)**2*b(4)) * Ck
      sum  = sum + term

      ! k = 12
      Ck   = (12-1)*Ck
      term = (b(12) + b(3)*b(9) + b(4)*b(8) + b(5)*b(7) + 0.5d0*b(6)**2
     &        + 0.5d0*b(3)**2*b(6) + b(3)*b(4)*b(5) + (1d0/6d0)*b(4)**3
     &        + (1d0/24d0)*b(3)**4) * Ck
      sum  = sum + term

      ! Check if converged here
      IF (ABS(term) .LT. 1d-16*sum) THEN
        expand_lnIu = LOG(sum)
        RETURN
      END IF

      ! k = 14
      Ck   = (14-1)*Ck
      term = (b(14) + b(3)*b(11) + b(4)*b(10) + b(5)*b(9) + b(6)*b(8)
     &              + 0.5d0*b(7)**2
     &        + 0.5d0*b(3)**2*b(8) + b(3)*b(4)*b(7) + b(3)*b(5)*b(6)
     &              + 0.5d0*b(4)**2*b(6) + 0.5d0*b(4)*b(5)**2
     &        + (1d0/6d0)*b(3)**3*b(5) + 0.25d0*b(3)**2*b(4)**2
     &       ) * Ck
      sum  = sum + term

      ! k = 16
      Ck   = (16-1)*Ck
      term = (b(16) + b(3)*b(13) + b(4)*b(12) + b(5)*b(11) + b(6)*b(10)
     &              + b(7)*b(9) + 0.5d0*b(8)**2
     &        + 0.5d0*b(3)**2*b(10) + b(3)*b(4)*b(9) + b(3)*b(5)*b(8)
     &              + b(3)*b(6)*b(7) + 0.5d0*b(4)**2*b(8)
     &              + b(4)*b(5)*b(7) + 0.5d0*b(4)*b(6)**2
     &              + 0.5d0*b(5)**2*b(6)
     &        + (1d0/6d0)*b(3)**3*b(7) + 0.5d0*b(3)**2*b(4)*b(6)
     &              + 0.25d0*b(3)**2*b(5)**2 + 0.5d0*b(3)*b(4)**2*b(5)
     &              + (1d0/24d0)*b(4)**4
     &        + (1d0/24d0)*b(3)**4*b(4)
     &       ) * Ck
      sum  = sum + term

      ! Remaining terms are of 18th order.  Hopefully, this is
      ! sufficient.  A guesstimate for the error is:
      !   (MAX[kappa,n*rho]*sigmau^3)^6
      ! If the last term is not sufficiently small, we will return
      ! a very large number as a flag to indicate convergence
      ! failure.
      IF ((ABS(term) .GT. 1d-10*sum) .OR. (sum .LE. 0d0)) THEN
        expand_lnIu = HUGE(1d0)
        RETURN
      END IF

      expand_lnIu = LOG(sum)

      END FUNCTION expand_lnIu


      ! ----------------------------------------------------------------
      ! Routine to find the location y=y0 and value g0=g(y0) of the
      ! global maximum of the  integrand e^g(y) for the integral in the
      ! nintegrate routine as described above, with:
      !   g(y) = n ln(thetab+thetas*e^y) - ln(n!) - (thetab+thetas*e^y)
      !            - y^2/(2*sigma^2)
      !        = -1/sigma^2 [y^2/2 + alpha*(r+e^y) - eta*ln(r+e^y)
      !                      + sigma^2*ln(n!) - eta*ln(thetas)]
      ! where:
      !   sigma2 = sigma^2
      !   alpha  = sigma^2 thetas
      !   beta   = sigma^2 thetab
      !   eta    = sigma^2 n
      !   r      = beta/alpha
      ! The routine also returns:
      !   sigmau =  1 / [-g''(y0)]^{1/2}
      !
      ! The maximum occurs at one of the zeros of:
      !   f(y) = -sigma^2 g'(y)
      !        = y + alpha*e^y - eta*e^y/(r+e^y)
      ! While there is typically only one zero, there can be up to
      ! three for some choices of parameters, such as (1+r)*alpha <<
      ! eta; care is taken here to find the highest (global) maximum.
      !
      ! Function returns true if and only if it successfully found the
      ! maximum; return value of false indicates that thetas ~ 0 to
      ! within the usable precision of this function.
      LOGICAL FUNCTION find_maximum(n,thetab,thetas,sigma,y0,g0,sigmau)
      IMPLICIT NONE
      INTEGER n
      REAL*8 thetab,thetas,sigma,y0,g0,sigmau
      INTEGER I,J,Nintervals
      REAL*8 sigma2,alpha,beta,eta,r,g2,expy,lnA,z0,
     &       yP,yG,ymin,ymax,p,q,D,t(3),tc,ttmp,
     &       yarr(0:4),farr(0:4),yi,ytmp,gtmp
      LOGICAL found,foundtmp
      COMPLEX*16 u,C1,C2,Z(3)
      REAL*8 lngamma,lambertw,lambertwln
      PARAMETER(C1=(-0.5d0,0.86602540378443865d0),
     &          C2=(-0.5d0,-0.86602540378443865d0))

      ! This is only set true when the function terminates successfully
      find_maximum = .false.

      !WRITE(*,'(A,I,3(G))') 'findmax arguments: ',n,thetab,thetas,sigma

      ! Define several useful quantities
      sigma2 = sigma**2
      alpha  = sigma2*thetas
      beta   = sigma2*thetab
      eta    = sigma2*n
      r      = thetab/thetas

      ! ---------------------------------------
      ! Special case: thetab=0
      IF (thetab .EQ. 0d0) THEN
        ! Maximum occurs at A = z0 e^z0 with:
        !   A  = alpha e^eta
        !   z0 = eta - y0
        ! Note: A will easily overflow.  We work with ln(A) instead.
        !A  = alpha * EXP(eta)
        !z0 = lambertw(A)
        lnA = eta + LOG(alpha)
        z0  = lambertwln(lnA)
        y0     = eta - z0
        sigmau = sigma / SQRT(1+z0)
        g0     = y0*(1-0.5d0*y0)/sigma2 + n*(LOG(thetas)+y0-1)
     &            - lngamma(n+1d0)
        RETURN
      END IF

      ! ---------------------------------------
      ! Special case: n=0
      IF (n .EQ. 0) THEN
        ! Maximum occurs at alpha = -y0 e^{-y0}.
        y0     = -lambertw(alpha)
        sigmau = sigma / SQRT(1-y0)
        g0     = y0*(1-0.5d0*y0)/sigma2 - thetab  - lngamma(n+1d0)
        RETURN
      END IF

      ! ---------------------------------------
      ! Special case: thetab+thetas=n
      IF (thetab+thetas .EQ. n) THEN
        y0     = 0d0
        sigmau = sigma / SQRT(1+alpha/(1+r))
        g0     = n*LOG(thetab+thetas) - (thetab+thetas)
     &            - lngamma(n+1d0)
        RETURN
      END IF

      ! ---------------------------------------
      ! General case
      ! We use a numerical search here.

      ! There are contributions to g(y) from a Gaussian and Poisson
      ! component.  Each component has a peak; the sum may have a
      ! single peak or two distinguishable peaks.  In any case, the
      ! peak(s) in g(y) should lie between the centers of the
      ! individual Gaussian and Poisson peaks, so we can use their
      ! locations as boundaries for the numerical search.
      !
      ! Gaussian peaks at y=0.
      yG = 0d0
      ! Poisson peaks at thetab + thetas*e^y = n.
      ! No peak for n <= thetab; use another lower bound instead.
      IF (n .GT. thetab) THEN
        yP = LOG((n-thetab)/thetas)
      ELSE
        yP = - lambertw(alpha)
      END IF

      ! Global boundaries for search
      ymin = MIN(yG,yP)
      ymax = MAX(yG,yP)

      ! If there is more than one extremum in g(y), they must be
      ! separated by points where g''(y)=0 (or f'(y)=0).  We will
      ! then check for the number of real roots of f'(y), which is
      ! a cubic equation in x = e^y.
      ! We use the Tschirnhaus transformation t = x + tc, with tc
      ! as defined below, to reduce the problem to a monic trinomial:
      !   t^3 + p*t + q = 0
      tc = (1+2*beta) / (3*alpha)
      p  = - (3*beta*eta + (1-beta)**2) / (3*alpha**2)
      q  = (2*(1-beta)**3 + 9*beta*eta*(1+2*beta)) / (27*alpha**3)

      ! Bail out here if any of these quantities is infinite, indicating
      ! that thetas is too small for the calculation to proceed (i.e. ~0).
      IF (ANY((/ABS(tc),ABS(p),ABS(q),ABS(0.25d0*q**2),ABS((1d0/27d0)*p**3)/).gt.HUGE(tc))) RETURN

      ! Now that we know D is not going to be NaN, we can proceed.
      D  = 0.25d0*q**2 + (1d0/27d0)*p**3

      ! We set up intervals to search, ensuring at most one extremum
      ! in each interval.
      Nintervals = 1
      yarr(0) = ymin
      yarr(1) = ymax

      ! For D>0, there is one real root of the cubic equation in terms
      ! of x, while for D<0, there are three real roots.  Since x=e^y,
      ! only positive roots are relevant here.  Due to the nature of
      ! the problem, at least one root will be negative.  For a
      ! double peak to occur in g(y), there must be two positive roots
      ! in x.  In this latter case, we break up the region into
      ! subintervals separated by these roots as these subintervals
      ! are then ensured to have at most one maximum in g(y).
      IF (D .LT. 0d0) THEN
        u = (-0.5d0*q - SQRT((1d0,0d0)*D))**(1d0/3d0)
        Z(1) = u - p/(3*u)
        Z(2) = C1*u - p/(3*C1*u)
        Z(3) = C2*u - p/(3*C2*u)
        t(1) = DBLE(Z(1))
        t(2) = DBLE(Z(2))
        t(3) = DBLE(Z(3))
        ! Sort
        DO I=1,2
          DO J=1,3-I
            IF (t(J+1) .LT. t(J)) THEN
              ttmp   = t(J)
              t(J)   = t(J+1)
              t(J+1) = ttmp
            END IF
          END DO
        END DO
        ! Check if positive root and within desired range
        DO I=1,3
          IF (t(I) .LE. tc) CYCLE
          ytmp = LOG(t(I)-tc)
          IF ((ytmp .GE. ymin) .AND. (ytmp .LE. ymax)) THEN
            Nintervals = Nintervals+1
            yarr(Nintervals)   = yarr(Nintervals-1)
            yarr(Nintervals-1) = ytmp
          END IF
        END DO
      END IF

      DO I=0,Nintervals
        ytmp = yarr(I)
        expy = EXP(ytmp)
        farr(I) = ytmp - eta*expy/(r+expy) + alpha*expy
      END DO

      g0 = -HUGE(1d0)
      found = .FALSE.

      ! Find potential maxima within each interval
      DO I=1,Nintervals
        ! Interval contains minimum or no extremum
        IF ((farr(I-1) .GT. 0d0) .OR. (farr(I) .LT. 0d0)) CYCLE
        IF ((yG .GE. yarr(I-1)) .AND. (yG .LE. yarr(I))) THEN
          yi = yG
        ELSE IF ((yP .GE. yarr(I-1)) .AND. (yP .LE. yarr(I))) THEN
          yi = yP
        ELSE
          yi = 0.5d0*(yarr(I-1)+yarr(I))
        END IF
        CALL find_zero(alpha,r,eta,yi,yarr(I-1),yarr(I),ytmp,foundtmp)
        IF (.NOT. foundtmp) CYCLE
        found = .TRUE.
        expy  = EXP(ytmp)
        gtmp  = - ( 0.5d0*ytmp**2 - eta*LOG(thetas*(r+expy))
     &              + alpha*(r+expy) )
     &            / sigma2
     &          - lngamma(n+1d0)
        IF (gtmp .GT. g0) THEN
          y0 = ytmp
          g0 = gtmp
        END IF
      END DO

      ! If we haven't found a zero, try again over the whole interval
      ! with larger boundaries.
      IF (.NOT. found) THEN
        CALL find_zero(alpha,r,eta,yi,ymin-1d0,ymax+1d0,y0,found)
        expy = EXP(y0)
        g0   = - ( 0.5d0*y0**2 - eta*LOG(thetas*(r+expy))
     &             + alpha*(r+expy) )
     &           / sigma2
     &         - lngamma(n+1d0)
      END IF

      ! Failed to find any maximum.
      ! This should not happen with the above algorithms except for
      ! possible finite precision effects at various boundaries,
      ! which hopefully the last resort case above will catch....
      IF (.NOT. found) THEN
        WRITE(*,*) 'ERROR: failed to find integrand maximum'
     &             // ' (nulike_lnpiln)'
        WRITE(*,*) 'Called with arguments: ',n,thetab,thetas,sigma,y0,g0,
     &   sigmau
        STOP('Please report this to Chris Savage (chris@savage.name)')
      END IF

      ! Now calculate other quantites at y0.
      expy   = EXP(y0)
      g2     = - (1d0 - eta*r*expy/(r+expy)**2 + alpha*expy) / sigma2
      sigmau = 1d0 / SQRT(-g2)

      find_maximum = .true.

      END FUNCTION find_maximum


      ! ----------------------------------------------------------------
      ! Function to find the location y0 for the zero of:
      !   f(y) = y + alpha*e^y - eta*e^y/(r+e^y)
      ! witin the range ymin < y < ymax, using a combination of the
      ! Newton/Halley methods and a bisection method.  The function
      ! should only be called for ymin and ymax such that f(ymin) and
      ! f(ymax) are of opposite signs, so that at least one zero is
      ! guaranteed to exist in the region and can be found through
      ! bisection.  If multiple zeros exist in the region, there is
      ! no guarantee which zero will be found: it will depend on both
      ! the boundaries and the initial point yi used to begin the
      ! Newton/Halley method.  The flag "found" indicates if a zero was
      ! found to within a reasonable tolerance within a reasonable
      ! number of iterations.
      !
      ! Note that f(y) = -sigma^2 g'(y) so that this function finds
      ! extrema of g(y).
      !
      SUBROUTINE find_zero(alpha,r,eta,yi,ymin,ymax,y0,found)
      IMPLICIT NONE
      REAL*8 alpha,r,eta,yi,ymin,ymax,y0
      LOGICAL found
      REAL*8 expy,f0,f1,f2,deltay,ya,yb,fa,fb,y_prev,f_prev
      INTEGER K,MAX_ITERATIONS,MAX_BISECTIONS
      PARAMETER(MAX_ITERATIONS=8,MAX_BISECTIONS=40)

      !WRITE(*,'(A,6(G))') 'findzero arguments:',alpha,r,eta,yi,ymin,ymax

      found = .FALSE.
      deltay = HUGE(1d0)
      y0 = yi

      ! Bracketing -----------------------------
      ! We bracket the zero of f(y) between [ya,yb].
      ya   = ymin
      expy = EXP(ya)
      fa   = ya - eta*expy/(r+expy) + alpha*expy
      yb   = ymax
      expy = EXP(yb)
      fb   = yb - eta*expy/(r+expy) + alpha*expy

      ! Newton/Halley's method -----------------
      ! Convergence for this technique appears to be very rapid for
      ! the vast majority of cases, often reaching as close to the
      ! zero as it can in only 2-4 iterations.  Due to finite
      ! precision effects, however, this routine can get stuck near
      ! the zero without converging to full double precision in y0.
      ! Note that stability has not been checked for all parameters
      ! and thus even approximate convergence is not guaranteed!
      ! We will keep bracketing the zero of g'(y) and use a bisection
      ! method below if this Newton's/Halley's method does not appear
      ! to converge.  To aid in the bracketing, it is useful to allow
      ! the Newton/Halley method to go through a few extra iterations
      ! than would otherwise be necessary.
      ! NOTE: Newton is more stable, but Halley is quicker.
      !       We will use the faster Halley as we can fall back to the
      !       safe bisection method if it fails.  The Newton method
      !       can be used instead by uncommenting the appropriate
      !       lines below).
      K = 0
      !WRITE(*,'(A,I3,6(G))') ' DEBUG N/H:',K,y0,0d0,ya,yb,fa,fb
      DO WHILE ((.NOT. found) .AND. (K .LT. MAX_ITERATIONS))
        K = K+1
        expy = EXP(y0)
        f0 = y0 - eta*expy/(r+expy) + alpha*expy
        f1 = 1d0 - eta*r*expy/(r+expy)**2 + alpha*expy
        f2 = - eta*r*(r-expy)*expy/(r+expy)**3 + alpha*expy
        ! Update bracketing
        IF ((f0*fa .GE. 0d0) .AND. (y0 .GT. ya)) THEN
          ya = y0
          fa = f0
        END IF
        IF ((f0*fb .GE. 0d0) .AND. (y0 .LT. yb)) THEN
          yb = y0
          fb = f0
        END IF
        ! Newton ---------------
        !deltay = - f0/f1
        ! Halley ---------------
        deltay = - 2*f0*f1 / (2*f1**2 - f0*f2)
        ! ----------------------
        y_prev = y0
        f_prev = f0
        y0 = y0 + deltay
        IF (ABS(deltay) .LE. 10*EPSILON(y0)*ABS(y0)) found=.TRUE.
        !IF (ABS(deltay) .LE. 10*EPSILON(y0)*ABS(y0)) THEN
        !  found = .TRUE.
        !  RETURN
        !END IF
        ! If we would go beyond the current constraints, Newton/Halley
        ! may be unstable.  We will go to the bisection method instead.
        IF ((y0 .LT. ya) .OR. (y0 .GT. yb)) EXIT
        !WRITE(*,'(A,I3,6(G))') ' DEBUG N/H:',K,y0,deltay,ya,yb,fa,fb
      END DO
      !WRITE(*,'(A,I3,6(G))') ' DEBUG N/H:',K,y0,deltay,ya,yb,fa,fb

      ! Bisection ------------------------------
      ! If the maximum was not found, there is either a numerical
      ! issue or convergence failed.  One (not uncommon) numerical
      ! issue is a search that gets stuck cycling between the same
      ! points around the maximum.  For both these cases, we use a
      ! bisection method starting with the interval [ya,yb], which
      ! may already very closely bracket the zero of f(y).
      K = 0
      !WRITE(*,'(A,I3,6(G))') ' DEBUG BIS:',K,y0,deltay,ya,yb,fa,fb
      DO WHILE ((.NOT. found) .AND. (K .LT. MAX_BISECTIONS))
        K = K+1
        deltay = 0.5d0*(yb-ya)
        y0 = ya + deltay
        expy = EXP(y0)
        f0 = y0 - eta*expy/(r+expy) + alpha*expy
        IF ((f0*fa .GE. 0d0) .AND. (y0 .GT. ya)) THEN
          ya = y0
          fa = f0
        END IF
        IF ((f0*fb .GE. 0d0) .AND. (y0 .LT. yb)) THEN
          yb = y0
          fb = f0
        END IF
        IF (ABS(deltay) .LE. 10*EPSILON(y0)*ABS(y0)) found=.TRUE.
        !IF (ABS(deltay) .LE. 10*EPSILON(y0)*ABS(y0)) THEN
        !  found = .TRUE.
        !  RETURN
        !END IF
        !WRITE(*,'(A,I3,6(G))') ' DEBUG BIS:',K,y0,deltay,ya,yb,fa,fb
      END DO
      !WRITE(*,'(A,I3,6(G))') ' DEBUG BIS:',K,y0,deltay,ya,yb,fa,fb

      ! Checks ---------------------------------
      ! Loosen the tolerance at this point rather than generate an
      ! error if it looks like we are at least close.
      IF (ABS(deltay) .LT. 1d6*EPSILON(y0)*ABS(y0)) found=.TRUE.

      END SUBROUTINE find_zero


      ! ----------------------------------------------------------------
      ! Function to calculate the quantity:
      !   e^z - 1
      ! Small z must be handled carefully as 1 is the first term in
      ! the expansion of e^z about zero.
      !
      ! Precision is ~ 2*EPSILON.
      !
      REAL*8 FUNCTION expdiff0(z)
      IMPLICIT NONE
      REAL*8 z,zabs
      INTEGER k,klast

      ! If z is not too small, we can just evaluate explicitly.
      ! Precision in this range is ~ few*EPSILON.
      IF ((z .LT. -0.5d0) .OR. (z .GT. 0.5d0)) THEN
        expdiff0 = EXP(z) - 1d0
        RETURN
      END IF

      zabs = ABS(z)

      ! We will use an expansion about zero, canceling the leading
      ! term.  We keep terms in the expansion up to x^klast/klast!.
      ! klast is chosen below to ensure double precision (within
      ! ~ few*EPSILON).  The expansion is:
      !     \sum_{k=1}^{\infty} z^k / k!
      ! The precision is approximately:
      !     x^{klast} / (klast+1)!
      ! where klast is the last term in the sum.

      IF (zabs .LT. 2.5d-4) THEN
        klast = 4
      ELSE IF (zabs .LT. 0.043d0) THEN
        klast = 8
      ELSE IF (zabs .LT. 0.28d0) THEN
        klast = 12
      ELSE
        klast = 16
      END IF

      ! Use expansion about zero
      expdiff0 = z / klast
      DO k = klast-1,1,-1
        expdiff0 = z * (1 + expdiff0) / k
      END DO

      END FUNCTION expdiff0


      ! ----------------------------------------------------------------
      ! Function to calculate the quantity:
      !   e^z - (1 + z + (1/2)*z^2)
      ! Small z must be handled carefully as the three terms on the
      ! right correspond to the first three terms of the expansion
      ! of e^z about zero.
      !
      ! Precision is ~ 10*EPSILON.
      !
      REAL*8 FUNCTION expdiff2(z)
      IMPLICIT NONE
      REAL*8 z,zabs
      INTEGER k,klast

      ! If z is not too small, we can just evaluate explicitly.
      ! Precision in this range is within ~ 10*EPSILON.
      IF ((z .LT. -0.7d0) .OR. (z .GT. 1.0d0)) THEN
        expdiff2 = EXP(z) - (1 + z + 0.5d0*z**2)
        RETURN
      END IF

      zabs = ABS(z)

      ! We will use an expansion about zero, canceling the leading
      ! terms.  We keep terms in the expansion up to x^klast/klast!.
      ! klast is chosen below to ensure near double precision (within
      ! ~ 10*EPSILON).  The expansion is:
      !     \sum_{k=3}^{\infty} z^k / k!
      ! The precision is approximately:
      !     6*x^{klast-2} / (klast+1)!
      ! where klast is the last term in the sum.

      IF (zabs .LT. 7d-4) THEN
        klast = 6
      ELSE IF (zabs .LT. 0.08d0) THEN
        klast = 10
      ELSE IF (zabs .LT. 0.45d0) THEN
        klast = 14
      ELSE IF (zabs .LT. 1.11d0) THEN
        klast = 18
      ELSE
        klast = 22
      END IF

      ! Use expansion about zero
      expdiff2 = z / klast
      DO k = klast-1,3,-1
        expdiff2 = z * (1 + expdiff2) / k
      END DO
      expdiff2 = z**2 * expdiff2 / 2

      END FUNCTION expdiff2


      ! ----------------------------------------------------------------
      ! Function to calculate the quantity:
      !   (e^z - 1)^2 - z^2
      ! Small z must be handled carefully as the first three terms
      ! of the expansion about zero cancel.
      !
      ! Precision is ~ 10*EPSILON.
      !
      REAL*8 FUNCTION expdiffsq(z)
      IMPLICIT NONE
      REAL*8 z,zabs
      INTEGER k,klast,twok

      ! If z is not too small, we can just evaluate explicitly.
      ! Precision in this range is ~ few*EPSILON.
      IF ((z .LT. -0.5d0) .OR. (z .GT. 0.5d0)) THEN
        expdiffsq = (EXP(z)-1)**2 - z**2
        RETURN
      END IF

      zabs = ABS(z)

      ! We will use an expansion about zero, canceling the leading
      ! terms.  We keep terms in the expansion up to x^klast/klast!.
      ! klast is chosen below to ensure double precision (within
      ! ~ few*EPSILON).  The expansion is:
      !     \sum_{k=3}^{\infty} (2^k - 2) z^k / k!
      ! The precision is approximately:
      !     (2^{klast+1}-2)*x^{klast-2} / (klast+1)!
      ! where klast is the last term in the sum.

      IF (zabs .LT. 4d-4) THEN
        klast = 6
      ELSE IF (zabs .LT. 0.06d0) THEN
        klast = 10
      ELSE IF (zabs .LT. 0.37d0) THEN
        klast = 14
      ELSE
        klast = 18
      END IF

      ! Use expansion about zero
      twok = 2**klast
      expdiffsq = (twok-2) * z / klast
      DO k = klast-1,3,-1
        twok = twok / 2
        expdiffsq = z * ((twok-2) + expdiffsq) / k
      END DO
      expdiffsq = z**2 * expdiffsq / 2

      END FUNCTION expdiffsq


      ! ----------------------------------------------------------------
      ! Function to calculate the quantity:
      !   ln(1+z)
      ! Small z must be handled carefully.  Note the expansion has
      ! no 0th order term, so this is equivalent to the difference
      ! with the 0th order expansion of ln(1+z), hence the name.
      !
      ! Precision is ~ 100*EPSILON.
      !
      REAL*8 FUNCTION lndiff0(z)
      IMPLICIT NONE
      REAL*8 z,zabs,zk
      INTEGER k,klast

      zabs = ABS(z)

      ! If z is not too small, we can just evaluate explicitly
      IF (zabs .GT. 0.01d0) THEN
        lndiff0 = LOG(1+z)
        RETURN
      END IF

      ! We will use an expansion about zero, keeping terms in the
      ! expansion up to x^klast/klast.
      ! klast is chosen below to give almost double precision.
      ! Precision can be as poor as ~ 100*EPSILON, but is better
      ! for smaller z.  The expansion is:
      !     \sum_{k=1}^{\infty} (-1)^{k+1} z^k / k
      ! The precision is approximately:
      !     x^klast / (klast+1)
      ! where klast is the last term in the sum.

      IF (zabs .LT. 1d-4) THEN
        klast = 4
      ELSE
        klast = 8
      END IF
      ! Go to klast=12 for |z| < 0.05
      ! Go to klast=16 for |z| < 0.1

      ! Use expansion about zero
      zk      = z
      lndiff0 = zk
      DO k = 2,klast
        zk      = -z*zk
        lndiff0 = lndiff0 + zk/k
      END DO

      END FUNCTION lndiff0


      ! ----------------------------------------------------------------
      ! Function to calculate the quantity:
      !   ln(1+z) - (z - (1/2)*z^2)
      ! Small z must be handled carefully as the two terms on the
      ! right correspond to the first two terms of the expansion
      ! of ln(1+z) about zero.
      !
      ! Precision is ~ 100*EPSILON.
      !
      REAL*8 FUNCTION lndiff2(z)
      IMPLICIT NONE
      REAL*8 z,zabs,zk
      INTEGER k,klast

      ! Special case
      IF (z .LE. -1d0) THEN
        lndiff2 = -HUGE(lndiff2)
        RETURN
      END IF

      zabs = ABS(z)

      ! If z is not too small, we can just evaluate explicitly
      IF (zabs .GT. 0.2d0) THEN
        lndiff2 = LOG(1+z) - z*(1 - 0.5d0*z)
        RETURN
      END IF

      ! We will use an expansion about zero, canceling the leading
      ! terms.  We keep terms in the expansion up to x^klast/klast.
      ! klast is chosen below to give almost double precision.
      ! Precision can be as poor as ~ 100*EPSILON, but is better
      ! for smaller z.  The expansion is:
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
      zk      = z**3
      lndiff2 = zk/3
      DO k = 4,klast
        zk      = -z*zk
        lndiff2 = lndiff2 + zk/k
      END DO

      END FUNCTION lndiff2


      ! ----------------------------------------------------------------
      ! Function to calculate the quantity:
      !   ln(1+rho*(e^z-1)) - rho*z - rho*(1-rho) z^2 / 2
      ! Small z must be handled carefully as the two terms on the
      ! right correspond to the first two terms of the expansion
      ! of the ln term about zero.
      !
      ! Precision is ~ 100*EPSILON.
      !
      REAL*8 FUNCTION lnexpdiff2(rho,z)
      IMPLICIT NONE
      REAL*8 rho,z,v

      ! Special case
      IF (rho .EQ. 0d0) THEN
        lnexpdiff2 = 0d0
        RETURN
      END IF

      ! Large z case:
      !   ln(1+rho*(e^z-1)) -> ln(rho+e^-z) + z
      IF (z .GT. 100) THEN
        lnexpdiff2 = LOG(rho+EXP(-z)) + (1-rho)*z*(1-0.5d0*rho*z)
        RETURN
      END IF

      ! Useful quantity
      v = rho*expdiff0(z)

      ! If v not too small (implies z not small), we should evaluate
      ! explicitly
      IF (ABS(v) .GT. 1d0) THEN
        lnexpdiff2 = LOG(1+v) - rho*z*(1+0.5d0*(1-rho)*z)
        RETURN
      END IF

      ! With v=rho*(e^z-1), we can rewrite the desired quantity as:
      !   ln(1+v) - rho*z - rho*(1-rho) z^2 / 2
      !       = ln(1+v) - v + v^2/2
      !         + rho*[e^z - (1+z+z^2/2)]
      !         - rho^2/2 [(e^z-1)^2 - z^2]
      ! where each line is of order z^3.
      lnexpdiff2 = lndiff2(v) + rho*(expdiff2(z)-0.5d0*rho*expdiffsq(z))

      END FUNCTION lnexpdiff2

      END FUNCTION


