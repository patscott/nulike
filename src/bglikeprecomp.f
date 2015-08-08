***********************************************************************
*** nulike_bglikeprecomp calls the calculation of the p-value for the
*** background in IceCube calculations based on Poissonian statistics.
***        
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: May, Jul, 2011
*** Modified: Jun 3, 7 2014 
***********************************************************************

      subroutine nulike_bglikeprecomp

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      real*8 nulike_pval

      BGpvalPoissonian(analysis) = nulike_pval(nEvents(analysis),
     & theta_BG(analysis),0.d0,1.d0/dsqrt(dble(nEvents(analysis))))
      pvalBGPoisComputed(analysis) = .true.

      end subroutine nulike_bglikeprecomp
