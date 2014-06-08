***********************************************************************
*** nulike_signal computes the predicted number of neutrino events due
*** to neutralino annihilation, saving global variables
*** for later access by likelihood codes.
***        
*** Input:      muonyield       External double function that returns
***                             the differential muon/neutrino flux
***                             at the detector in units of m^-2 GeV^-1
***             annrate         Annihilation rate (s^-1) 
***             logmw           log_10(m_WIMP / GeV)
***             like            Likelihood type (2012 or 2014)
***
*** Output:     theta_S         predicted number of signal events.
*** 
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 22, 2011
*** Modified: March 6 2014
*** Modified: Jun 3, 6 2014
***********************************************************************


      double precision function nulike_signal(muonyield, annrate, logmw, like)

      implicit none
      include 'nulike.h'

      real*8 integral, eps, nulike_simpson, nulike_sigintegrand, logmw
      real*8 muonyield, upperLimit, theta_Snu, theta_Snubar, annrate
      integer like
      parameter (eps = 1.d-3)
      external muonyield, nulike_sigintegrand
 
      if (like .eq. 2012) then

      if (logmw .lt. effArea_logE(1,1,analysis)) then

        theta_Snu = 0.d0
        theta_Snubar = 0.d0
 
      else

        if (logmw .lt. effArea_logE(2,nBinsEA(analysis),analysis)) then
          upperLimit = logmw
        else
          upperLimit = effArea_logE(2,nBinsEA(analysis),analysis)
        endif

        ptypeshare = 1
        integral = nulike_simpson(nulike_sigintegrand,muonyield,
     &   effArea_logE(1,1,analysis),upperLimit,eps)
        theta_Snu = integral * dlog(10.d0) * exp_time(analysis) * annrate

        ptypeshare = 2
        integral = nulike_simpson(nulike_sigintegrand,muonyield,
     &   effArea_logE(1,1,analysis),upperLimit,eps)
        theta_Snubar = integral * dlog(10.d0) * exp_time(analysis) * annrate

      endif

      nulike_signal = theta_Snu + theta_Snubar 

      else if (like .eq. 2014) then

        write(*,*) '2014 likelihood not implemented'
        stop

      else 

        write(*,*) 'Unrecognised likelihood type in nulike_signal: ',like
        write(*,*) 'Quitting...'
        stop

      endif

      end function nulike_signal

