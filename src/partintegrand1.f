***********************************************************************
*** nulike_partialintegrand1 does the innermost integral in the inner
*** double integral required to be computed to obtain partial likelihoods.
*** This routine is used only with the 2014 likelihood.
***
*** Input:              x              Bjorken x
***                     dsdxdy         cross-section function
***                                    (see partials.f for details)
***
*** Output:             integrand      m^2 degrees^-1 [energy estimator]^-1
***
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Jun 17, 2014
***********************************************************************

      real*8 function nulike_partintegrand1(x,dsdxdy)

      implicit none
      include 'nucommon.h'
      include 'nuprep.h'

      real*8 x, dsdxdy, nulike_partintegrand2, nulike_simpson2
      external dsdxdy, nulike_partintegrand2

      xshare = x 
      nulike_partintegrand1 = nulike_simpson2(nulike_partintegrand2,
     & dsdxdy, 0.d0, 1.d0, eps_partials)

      end function nulike_partintegrand1
