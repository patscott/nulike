***********************************************************************
*** nulike_speclike returns the contribution of a single event to the 
*** unbinned likelihood, based on the number of hit DOMs and the
*** theoretical energy spectrum of incoming neutrinos.
***
*** input:  nchan      observed number of hit DOMs for this event
***         theta_S_input total predicted number of signal events
***                     within analysis window (cut cone and superbin)
***         f_S        signal fraction; percentage of predicted counts
***                     expected to be due to signal rather than back-
***                     ground.
***         logEmin    log10(Emin/GeV), where Emin is the lower energy
***                     boundary of the analysis window (i.e. superbin)  
***         logEmax    log10(Emax/GeV), where Emax is the upper energy
***                     boundary of the analysis window (i.e. superbin)
***         muonyield  external double function that returns
***                     the differential muon/neutrino flux
***                     at the detector in units of m^-2 GeV^-1
*** output:            ln(Likelihood / chan^-1)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 22, 2011
*** Modified: March 6, 2014
***********************************************************************

      double precision function nulike_speclike(nchan,theta_S_input,
     & f_S,reset,logEmin,logEmax,muonyield)

      implicit none
      include 'nulike.h'

      integer nchan
      logical reset, savedSpecLikeFlags(nchan_maxallowed)
      real*8 theta_S_input, f_S, nulike_simpson, upperLimit
      real*8 signalpartiallike, bgpartiallike, integral, nulike_bgspec
      real*8 nulike_specintegrand, eps, logEmin, logEmax
      real*8 savedSpecLikes(nchan_maxallowed), muonyield
      parameter (eps = 1.d-2)
      external nulike_specintegrand, muonyield
      save savedSpecLikeFlags, savedSpecLikes
      
      nchanshare = nchan
      thetashare = theta_S_input

      if (nchan .gt. nchan_maxallowed) 
     & stop 'nchan>nchan_maxallowed in nulike_speclike'

      !Reset saved spectral likelihoods if requested
      if (reset) then
        savedSpecLikes = 0.d0
        savedSpecLikeFlags = .false.
        reset = .false.
      endif

      !If a cached result is available, use it - otherwise, calculate...
      if (savedSpecLikeFlags(nchan)) then

        nulike_speclike = savedSpecLikes(nchan)
        return

      endif
         
      if (theta_S_input .ne. 0.d0 .and. log10mwimp .gt. logEmin) then

        if (log10mwimp .lt. logEmax) then
          upperLimit = log10mwimp
        else
          upperLimit = logEmax
        endif

        !Find the part of the spectral likelihood associated with the signal
        integral = nulike_simpson(nulike_specintegrand,muonyield,
     &   logEmin,upperLimit,eps)

      else

        integral = 0.d0

      endif    

      signalpartiallike = f_S * integral * dlog(10.d0)

      !Find the part associated with the background spectrum
      bgpartiallike = (1.d0-f_S) * nulike_bgspec(nchan)

      nulike_speclike = dlog(signalpartiallike + bgpartiallike)
      savedSpecLikeFlags(nchan) = .true.
      savedSpecLikes(nchan) = nulike_speclike

      end function nulike_speclike

