***********************************************************************
*** nulike_bglikeprecomp calls the calculation of the p-value for the
*** background in IceCube calculations based on Poissonian statistics.
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: May, July 2011
***********************************************************************

      subroutine nulike_bglikeprecomp

      implicit none
      include 'nulike.h'

      real*8 nulike_pval

      BGpvalPoissonian = nulike_pval(sum(nEvents_inEAErrBins),sum(theta_BG),0.d0)
      pvalBGPoisComputed = .true.

      end subroutine nulike_bglikeprecomp
