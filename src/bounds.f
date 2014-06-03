***********************************************************************
*** nulike_bounds returns counts, likelihoods and p-values from neutrino
*** telescope searches for dark matter annihilation in the Sun. 
***
*** input: mwimp         WIMP mass (GeV)
***
***        ann_rate      Annihilation rate (s^-1) 
***
***        liketype      Sets combination of data to use in likelihood
***                         calculations
***                        1 => Number of events only
***                        2 => Number of events and event arrival angles
***                        3 => Number of events and energy estimator (nchan 
***                             = number of lit DOMs)
***                        4 => Number of events, event arrival angles and 
***                             energy estimator
***
***        pvalFromRef   T => calculate the p-value with reference to a user-
***                           specified likelihood, assuming that the log-
***                           likelihood follows a chi^2 distriubution in
***                           the vicinity of the reference point (i.e.
***                           essentially assume that the reference corresponds
***                           to the best-fit log likelihood, and that Wilks' 
***                           theorem is satisfied).
***                      F => calculate the p-value with reference to the
***                           observed background
***                                                
***        referenceLike reference value of the (natural) log likelihood to
***                      use for the calculation of the p-value if 
***                      pvalFromRef=T; ignored otherwise.
***
***        dof           degrees of freedom to assume in the calculation of the
***                      p-value in the case that pvalFromRef=T; ignored 
***                      otherwise.  This should normally be set to the
***                      number of parameters over which the likelihood has
***                      not been maximised when the solution is perturbed 
***                      from the minimum (i.e. the number of free parameters
***                      in a scan minus the number profiled over).
***
***                              
*** output: Nsignal_predicted  Predicted number of signal events within 
***                      phi_cut (includes solar coronal BG)
***
***         NBG_expected Expected number of (non-solar) background
***                      events insided phi_cut
***
***         Ntotal_observed  Observed number of events inside phi_cut
***
***         lnlike       natural log of chosen likelihood
***
***         pvalue       derived p-value for chosen model
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Mar 20, 2011
*** Updated: Jul 21, 2011
***          Mar 6, 2014
***          Jun 3, 2014
***********************************************************************


      subroutine nulike_bounds(mwimp, ann_rate, muonyield, 
     & Nsignal_predicted, NBG_expected, Ntotal_observed, 
     & lnlike, pvalue, liketype, pvalFromRef, referenceLike, dof)

      implicit none
      include 'nulike.h'

      integer Ntotal_observed, liketype, j
      integer counted1, counted2, countrate
      real*8 Nsignal_predicted, NBG_expected, nulike_pval
      real*8 lnlike, pvalue, referenceLike, dof, DGAMIC, DGAMMA, muonyield
      real*8 nLikelihood, angularLikelihood, spectralLikelihood
      real*8 theta_tot, f_S, nulike_anglike, nulike_speclike, nulike_nlike
      real*8 deltalnlike, mwimp, ann_rate
      logical pvalFromRef, nulike_speclike_reset, doProfiling
      character (len=*) pref,f1,f2,f3,f4
      external muonyield
      !Hidden option for doing speed profiling
      parameter (doProfiling = .false.)
      parameter (pref = 'share/DarkSUSY/IC_data/')
      parameter (f1 = pref//'events_10deg_IC79.dat')
      parameter (f2 = pref//'energy_histograms_IC79.dat')
      parameter (f3 = pref//'BG_distributions_IC79.dat')
      parameter (f4 = pref//'nuEffArea_IC79.dat')
        
      if (doProfiling) call system_clock(counted1,countrate)

      !Check validity of user-provided liketype
      if (liketype .lt. 1 .or. liketype .gt. 4) then
        write(*,*) 'Invalid liketype provided to nulike_bounds.'
        write(*,*) 'Quitting...'
        stop
      endif
     
      !if nulike_init has not yet been called, call it with the default IC-79 options
      if (.not. nulike_init_called) then
        call nulike_init(f1,f2,f3,f4,20.d0,0.05d0,.true.,.true.)
      endif

      !Set internal WMIP mass and annihilation rate
      annrate = ann_rate
      log10mwimp = dlog10(mwimp)

      if (doProfiling) then
        call system_clock(counted2,countrate)
        write(*,*) 'Elapsed time initialising (s): ', 
     &   real(counted2 - counted1)/real(countrate)
      endif  

      !Calculate signal counts and spectrum
      call nulike_signal(muonyield)

      if (doProfiling) then
        call system_clock(counted1,countrate)
        write(*,*) 'Elapsed time on signal calc (s): ', 
     &   real(counted1 - counted2)/real(countrate)
      endif  

      !Calculate likelihood
      lnlike = 0.d0
      angularLikelihood = 0.d0
      spectralLikelihood = 0.d0
      theta_tot = theta_BG + theta_S
      f_S = theta_S / theta_tot
      if (liketype .ne. 1) then
        !Reset the saved spectral likelihoods next time nulike_speclike is run
        nulike_speclike_reset = .true.
        !Step through the individual events
        do j = 1, nEvents
          !Add in angular likelihood for this event
          if (liketype .eq. 2 .or. liketype .eq. 4) 
     &     angularLikelihood = angularLikelihood + nulike_anglike(
     &      events_cosphi(j),events_cosphiErr(j),f_S)
          !Add in spectral likelihood for this event
          if (liketype .eq. 3 .or. liketype .eq. 4) 
     &     spectralLikelihood = spectralLikelihood + nulike_speclike(
     &      events_nchan(j),theta_S,f_S,nulike_speclike_reset,
     &      effArea_logE(1,1),effArea_logE(2,nBinsEA),
     &      muonyield)
        enddo
      endif

      if (doProfiling) then
        call system_clock(counted2,countrate)
        write(*,*) 'Elapsed time on ang/spec likelihood calc (s): ', 
     &   real(counted2 - counted1)/real(countrate)
      endif  

      !Calculate the number likelihood
      nLikelihood = nulike_nlike(nEvents,
     & theta_tot,theta_S,EAErr,theoryErr)

      if (doProfiling) then
        call system_clock(counted1,countrate)
        write(*,*) 'Elapsed time on number likelihood calc (s): ', 
     &   real(counted1 - counted2)/real(countrate)
      endif  

      !Put together the number, angular and spectral likelihoods
      lnlike = lnlike + nLikelihood + angularLikelihood + 
     & spectralLikelihood
      
      !Calculate pvalue
      if (pValFromRef) then !compute p-value with reference to input log likelihood

        deltalnlike = max(0.d0, referenceLike-lnlike)
        pvalue = DGAMIC(dof*0.5d0,deltalnlike)/DGAMMA(dof*0.5d0)

      else                  !compute p-value with reference to background

        theta_tot = theta_BG + theta_S
        !p-value from Poissonian statistics
        if (.not. pvalBGPoisComputed) call nulike_bglikeprecomp
        pvalue = nulike_pval(nEvents, theta_tot, theta_S)
        pvalue = pvalue / BGpvalPoissonian
        
      endif

      if (doProfiling) then
        call system_clock(counted2,countrate)
        write(*,*) 'Elapsed time on pval calc (s): ',
     &   real(counted2 - counted1)/real(countrate)
      endif  

      !Export various counts
      NBG_expected = theta_BG
      Nsignal_predicted = theta_S
      Ntotal_observed = nEvents 

      end subroutine nulike_bounds

