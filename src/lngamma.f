!=======================================================================
! Calculates ln(gamma(x)) to double precision.
!     gamma(x) = \int_0^\infty dt t^{x-1} e^{-t}
! The argument must be positive.
! 
! 
!    Created by Chris Savage (savage@fysik.su.se)
!    2011/05/23
! 
!=======================================================================
! 
! Calculated using Lanczos' approximation.  Valid only for x > 0.
! 
! Lanczos's approximation:
!   Gamma(z+1) = (z+r+1/2)^{z+1/2} * EXP[-(z+r+1/2)] * sqrt(2*pi)
!                * (B_0 + \Sum_{k=0}^{N} B_k/(z+k) + eps)
! where r is some constant s.t. z+r+1/2 > 0.
! Coefficients in the series expansion are dependent on r and the
! number of terms at which the series is truncated.  The determination
! of these coefficients is too complicated to show here.
! 
! Note the formula above is valid for complex cases as well, as long as
! Re(x) > 1.  The following identity also applies for the complex case:
!   Gamma(1-z) = pi*z / Gamma(1+z) / sin(pi*z)
! 
! The coefficients are calculated as described in Section 6.7 of
! Glendon Pugh's thesis (2004), which provides an extensive discussion
! of the Lanczos approximation.  The coefficients are given for the
! following modified form of the formula:
!   Gamma(z+1) = 2 * sqrt(e/pi) * ((z+r+1/2)/e)^(z+1/2)
!                * (D_0 + \Sum_{k=0}^{N} D_k/(z+k) + eps)
! The quantities N and r are chosen from Table C.1 to achieve the
! desired precision:
!   Single precision:  N=4,  r=4.340882
!   Double precision:  N=10, r=10.900511  (see Table 8.5)
!   Quad precision:    N=21, r=22.618910  (see Table 9.4)
! 
      REAL*8 FUNCTION lngamma(x)
      IMPLICIT NONE
      REAL*8 x,x0,z,series_sum
      INTEGER K
      REAL*8 PI
      PARAMETER(PI=3.1415926535897932d0)
      ! Single precision ---------------
      !REAL*8 R
      !PARAMETER(R=4.340882)
      !INTEGER ND
      !PARAMETER(ND=4)
      !REAL*8 Dk(0:ND)
      !DATA Dk / +0.017549572,
      !&    +0.63894863,  -0.574296,  +0.1014372,   -0.0014755151  /)
      ! Double precision ---------------
      REAL*8 R
      PARAMETER(R=10.900511d0)
      INTEGER ND
      PARAMETER(ND=10)
      REAL*8 Dk(0:ND)
      DATA Dk / +2.4857408913875357d-5,
     &    +1.0514237858172197d0,  -3.4568709722201624d0,
     &    +4.5122770946689482d0,  -2.9828522532357666d0,
     &    +1.0563971157712671d0,  -1.9542877319164587d-1,
     &    +1.7097054340444122d-2, -5.7192611740430578d-4,
     &    +4.6339947335990564d-6, -2.7199490848860770d-9  /
      
      IF (x .LE. 0d0) THEN
        !lngamma = -HUGE(lngamma)
        !RETURN
        WRITE(*,*) 'x = ',x
        STOP 'ERROR: lngamma cannot be called with argument <= 0'
      END IF
      
      ! Error for x>0 should theoretically be smaller than 1e-15,
      ! but loss of precision in sum below might make it slightly worse
      ! (maybe ~ 1e-14?)
      IF (x .GE. 1d0) THEN
        ! Write as Gamma(x) = Gamma(1+x0)
        x0 = x - 1d0
      ELSE
        ! Special case: for 0<x<1, find Gamma(1+(1-x)) = Gamma(1+x0)
        ! and use identity Gamma(1-y) = pi*y / Gamma(1+y) / sin(pi*y)
        x0 = 1d0 - x
      END IF
      series_sum = Dk(0)
      DO K = 1,ND
        series_sum = series_sum + Dk(K)/(x0+K)
      END DO
      z = (x0+0.5d0)*(LOG(x0+R+0.5d0) - 1d0)
     &    + 0.5d0*(1d0 + LOG(4d0/PI))
     &    + LOG(series_sum)
      IF (x .LT. 1d0) THEN
        ! Special case: for 0<x<1, we found LogGamma(1+(1-x)).
        ! Now use identity Gamma(1-y) = pi*y / Gamma(1+y) / sin(pi*y)
        x0 = 1d0 - x
        z = LOG(PI*x0 / SIN(PI*x0)) - z
      END IF
      
      lngamma = z
      
      END FUNCTION

