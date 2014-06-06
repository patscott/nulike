***********************************************************************
*** nulike_signal computes the predicted number of neutrino events due
*** to neutralino annihilation, saving global variables
*** for later access by likelihood codes.
***        
*** Input:      muonyield       external double function that returns
***                             the differential muon/neutrino flux
***                             at the detector in units of m^-2 GeV^-1
***             like            likelihood type (2012 or 2014)
***
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 22, 2011
*** Modified: March 6 2014
*** Modified: Jun 3, 6 2014
***********************************************************************


      subroutine nulike_signal(muonyield, like)

      implicit none
      include 'nulike.h'

      real*8 integral, eps, nulike_simpson, nulike_sigintegrand
      real*8 muonyield, upperLimit
      integer like
      parameter (eps = 1.d-3)
      external muonyield, nulike_sigintegrand
 
      if (like .eq. 2012) then

      if (log10mwimp .lt. effArea_logE(analysis,1,1)) then

        theta_Snu = 0.d0
        theta_Snubar = 0.d0
 
      else

        if (log10mwimp .lt. effArea_logE(analysis,2,nBinsEA(analysis))) then
          upperLimit = log10mwimp
        else
          upperLimit = effArea_logE(analysis,2,nBinsEA(analysis))
        endif

        ptypeshare = 1
        integral = nulike_simpson(nulike_sigintegrand,muonyield,
     &   effArea_logE(analysis,1,1),upperLimit,eps)
        theta_Snu = integral * dlog(10.d0) * exp_time(analysis) * annrate

        ptypeshare = 2
        integral = nulike_simpson(nulike_sigintegrand,muonyield,
     &   effArea_logE(analysis,1,1),upperLimit,eps)
        theta_Snubar = integral * dlog(10.d0) * exp_time(analysis) * annrate

      endif

      theta_S = theta_Snu + theta_Snubar 

      else if (like .eq. 2014) then

        write(*,*) '2014 likelihood not implemented'
        stop

      else 

        write(*,*) 'Unrecognised likelihood type in nulike_signal: ',like
        write(*,*) 'Quitting...'
        stop

      endif

      end subroutine nulike_signal

