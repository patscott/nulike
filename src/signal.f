***********************************************************************
*** nulike_signal computes the predicted number of neutrino events due
*** to neutralino annihilation, saving global variables
*** for later access by likelihood codes.
***        
*** Input:      muonyield       external double function that returns
***                             the differential muon/neutrino flux
***                             at the detector in units of m^-2 GeV^-1
***
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 22, 2011
*** Modified: March 6 2014
*** Modified: Jun 3, 2014
***********************************************************************


      subroutine nulike_signal(muonyield)

      implicit none
      include 'nulike.h'

      real*8 integral, eps, nulike_simpson, nulike_sigintegrand
      real*8 muonyield, upperLimit
      parameter (eps = 1.d-3)
      external nulike_sigintegrand, muonyield
 
      if (log10mwimp .lt. EAlogE_inEAErrBins(1)) then

        theta_Snu = 0.d0
        theta_Snubar = 0.d0
 
      else

        if (log10mwimp .lt. EAlogE_inEAErrBins(2)) then
          upperLimit = log10mwimp
        else
          upperLimit = EAlogE_inEAErrBins(2)
        endif

        ptypeshare = 1
        integral = nulike_simpson(nulike_sigintegrand,muonyield,
     &   EAlogE_inEAErrBins(1),upperLimit,eps)
        theta_Snu = integral * dlog(10.d0) * exp_time * annrate

        ptypeshare = 2
        integral = nulike_simpson(nulike_sigintegrand,muonyield,
     &   EAlogE_inEAErrBins(1),upperLimit,eps)
        theta_Snubar = integral * dlog(10.d0) * exp_time *annrate

      endif

      theta_S = theta_Snu + theta_Snubar 

      end subroutine nulike_signal

