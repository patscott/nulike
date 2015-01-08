***********************************************************************
*** nulike_signal computes the predicted number of neutrino events due
*** to neutralino annihilation.
***        
*** Input:      nuyield         External double function that returns
***                             the differential neutrino flux
***                             at the detector in units of m^-2 GeV^-1
***                             annihilation^-1
***             context         A c_ptr passed in to nuyield when it is
***                             called
***             annrate         Annihilation rate (s^-1) 
***             logmw           log_10(m_WIMP / GeV)
***             like            Likelihood version (2012 or 2014)
***
*** Output:     theta_S         predicted number of signal events.
*** 
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 22, 2011
*** Modified: March 6 2014
*** Modified: Jun 3, 6, 8 2014
***********************************************************************


      double precision function nulike_signal(nuyield, context, annrate, logmw, like)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      real*8 integral, nulike_simpson, nulike_sigintegrand, logmw
      real*8 nuyield, upperLimit, theta_Snu, theta_Snubar, annrate
      real*8 nulike_specangintegrand, eps2012, eps2014
      integer like
      type(c_ptr) context
      parameter (eps2012 = 1.d-3, eps2014 = 1e-3)
      external nuyield, nulike_sigintegrand, nulike_specangintegrand, nulike_simpson
 
      ! Short-circuit if the mass is too low to produce any observable events.
      if (logmw .lt. sens_logE(1,1,analysis)) then
        nulike_signal = 0.d0
        return
      endif
      

      ! Switch according to likelihood version.
      select case (like)

      ! 2012 likelihood, as per arXiv:1207.0810
      case (2012)

        if (logmw .lt. sens_logE(2,nSensBins(analysis),analysis)) then
          upperLimit = logmw
        else
          upperLimit = sens_logE(2,nSensBins(analysis),analysis)
        endif

        ! Neutrinos
        ptypeshare = 1
        integral = nulike_simpson(nulike_sigintegrand,nuyield,context,
     &   sens_logE(1,1,analysis),upperLimit,eps2012)
        theta_Snu = integral * dlog(10.d0) * exp_time(analysis) * annrate

        ! Anti-neutrinos
        ptypeshare = 2
        integral = nulike_simpson(nulike_sigintegrand,nuyield,context,
     &   sens_logE(1,1,analysis),upperLimit,eps2012)
        theta_Snubar = integral * dlog(10.d0) * exp_time(analysis) * annrate

        ! Total
        nulike_signal = theta_Snu + theta_Snubar 

      !2014 likelihood, as per arXiv:141x.xxxx
      case (2014)

        eventnumshare = 0 ! Use effective area from previous tabulation.
        integral = nulike_simpson(nulike_specangintegrand,nuyield,context,
     &   sens_logE(1,1,analysis),logmw,eps2014)
        nulike_signal = integral * dlog(10.d0) * exp_time(analysis) * annrate

      case default
        write(*,*) "Unrecognised likelihood version in nulike_signal."
        write(*,*) "Quitting..."
        stop

      end select


      end function nulike_signal

