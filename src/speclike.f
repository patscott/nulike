***********************************************************************
*** nulike_speclike returns the contribution of a single event to the 
*** unbinned likelihood, based on the number of hit DOMs and the
*** theoretical energy spectrum of incoming neutrinos.
*** This routine is used only with the 2012 likelihood.
***
*** input:  nchan      observed number of hit DOMs for this event
***         theta_S    total predicted number of signal events
***                     within analysis window (cut cone)
***         f_S        signal fraction; percentage of predicted counts
***                     expected to be due to signal rather than back-
***                     ground.
***         annrate    Annihilation rate (s^-1) 
***         logmw      log_10(m_WIMP / GeV)
***         reset      Reset cached spectral likelihoods (these allow reuse
***                     for multiple events with the same spectral data).
***         logEmin    log10(Emin/GeV), where Emin is the lower energy
***                     boundary of the analysis energy range  
***         logEmax    log10(Emax/GeV), where Emax is the upper energy
***                     boundary of the analysis energy range
***         nuyield    external double function that returns
***                     the differential neutrino flux
***                     at the detector in units of m^-2 GeV^-1 
***                     annihilation^-1
***         context    A c_ptr passed in to nuyield when it is called
***
*** output:            ln(Likelihood / chan^-1)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 22, 2011
*** Modified: March 6, 2014
*** Modified: Jun 3, 2014
***********************************************************************

      double precision function nulike_speclike(nchan,theta_S,
     & f_S,annrate,logmw,reset,logEmin,logEmax,nuyield,context)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'
     
      logical reset, savedSpecLikeFlags(nchan_maxallowed)
      real*8 nchan, theta_S, f_S, annrate, nulike_simpson, upperLimit
      real*8 signalpartiallike, bgpartiallike, integral, nulike_bgspec
      real*8 nulike_specintegrand, eps, logEmin, logEmax
      real*8 savedSpecLikes(nchan_maxallowed), nuyield, logmw
      integer nchan_int
      type(c_ptr) context
      parameter (eps = 1.d-2)
      external nuyield, nulike_specintegrand
      save savedSpecLikeFlags, savedSpecLikes
      
      nchan_int = nint(nchan)
      nchanshare = nchan
      thetashare = theta_S
      annrateshare = annrate

      if (nchan_int .gt. nchan_maxallowed) 
     & stop 'nchan > nchan_maxallowed in nulike_speclike'

      !Reset saved spectral likelihoods if requested
      if (reset) then
        savedSpecLikes = 0.d0
        savedSpecLikeFlags = .false.
        reset = .false.
      endif

      !If a cached result is available, use it - otherwise, calculate...
      if (savedSpecLikeFlags(nchan_int)) then

        nulike_speclike = savedSpecLikes(nchan_int)
        return

      endif
         
      if (theta_S .ne. 0.d0 .and. logmw .gt. logEmin) then

        if (logmw .lt. logEmax) then
          upperLimit = logmw
        else
          upperLimit = logEmax
        endif

        !Find the part of the spectral likelihood associated with the signal
        integral = nulike_simpson(nulike_specintegrand,nuyield,context,
     &   logEmin,upperLimit,eps)

      else

        integral = 0.d0

      endif    

      signalpartiallike = f_S * integral * dlog(10.d0)

      !Find the part associated with the background spectrum
      bgpartiallike = (1.d0-f_S) * nulike_bgspec(nchan,2012)

      nulike_speclike = dlog(signalpartiallike + bgpartiallike)
      savedSpecLikeFlags(nchan_int) = .true.
      savedSpecLikes(nchan_int) = nulike_speclike

      end function nulike_speclike

